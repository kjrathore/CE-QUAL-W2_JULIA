# ==============================================================================
# Simulation.jl
#
# Time-stepping drivers, tying together Core/InitGeometry.jl,
# Hydrodynamics/FreeSurface.jl, Hydrodynamics/Transport.jl, IO/
# BoundaryReader.jl, and IO/OutputWriter.jl into runnable end-to-end
# simulations. Two entry points:
#
# - `run_zero_flow_sanity_check!` -- the original "zero-flow sanity check"
#   scope confirmed with the user (2026-08-12): no inflow/outflow, reduced
#   physics (see Hydrodynamics/FreeSurface.jl module docstring for exactly
#   what's real vs. stubbed), goal is a stable free-surface over many
#   timesteps with real TSR CSV output.
# - `run_forced_simulation!` (2026-08-23) -- real boundary-forced driver,
#   wiring in `IO/BoundaryReader.jl`'s QIN/TIN/QOT/QDTR/TDTR loading and
#   `Hydrodynamics/Transport.jl`'s `apply_temperature_sources!`/
#   `temperature_transport!`, both previously validated only via standalone
#   scripts, never an actual driver loop.
#
# Neither is a faithful port of `w2_4_win.f90`'s full driver.
#
# NOT PORTED (flagged, not guessed away):
# - The real adaptive timestep machinery (DLTF/DLTMIN/DLTD breakpoints,
#   stability-based DLT selection, `w2_4_win.f90`'s `AUTO_STEPPING`) -- both
#   drivers use one FIXED `dlt = tc.DLTMAX[1]` for every step. Correct for a
#   reduced-physics zero-flow run (nothing destabilizes the timestep); NOT
#   yet proven safe for `run_forced_simulation!` under strong real forcing
#   over a long run -- port real DLT selection before trusting that case.
# - Constituent transport (`constituent_transport!`, generic constituents,
#   WaterQuality/* kinetics) -- `run_forced_simulation!` calls
#   `temperature_transport!` but not `constituent_transport!` (no `CIN`
#   Tier-1 IO yet, see `IO/BoundaryReader.jl`).
# ==============================================================================

"""
    run_zero_flow_sanity_check!(g, geom, net, tc; nsteps, output_dir, output_segments, base_name="tsr")

Runs `nsteps` fixed-`dlt` reduced-physics hydrodynamic steps
(`hydrodynamic_step!`) starting from the current state of `g`/`geom`
(already passed through `init_geometry!`, `compute_dlxrho!`, and
`allocate_hydro_state!` by the caller), writing one TSR CSV row per step
(plus the initial condition) for each segment in `output_segments` via
`OutputWriter`.

`dlt` is fixed at `tc.DLTMAX[1]` for the whole run -- see module docstring
for why the real adaptive-timestep logic isn't ported yet. `output_segments`
is caller-supplied on purpose (not derived from Detroit-specific structure
here) -- keeps this function usable for any waterbody/branch topology.

Returns `(g, geom, jday_final)`.
"""
function run_zero_flow_sanity_check!(g, geom, net, tc; nsteps::Int, output_dir::AbstractString,
                                      output_segments::Vector{Int}, base_name::AbstractString="tsr")
    dlt = tc.DLTMAX[1]
    jday = tc.TMSTRT

    writer = OutputWriter.open_tsr_files(output_dir, base_name, output_segments)
    try
        OutputWriter.write_tsr_row!(writer, g, geom, jday, dlt)  # initial condition
        for _ in 1:nsteps
            hydrodynamic_step!(g, geom, net, dlt)
            jday += dlt / 86400.0
            OutputWriter.write_tsr_row!(writer, g, geom, jday, dlt)
        end
    finally
        OutputWriter.close_tsr_files!(writer)
    end
    return (g, geom, jday)
end

"""
    run_forced_simulation!(g, geom, net, tc, boundary; nsteps, output_dir, output_segments, base_name="tsr")

Real boundary-forced driver, closing the gap flagged since 2026-08-15 --
`temperature_transport!`/`apply_temperature_sources!` and `IO/
BoundaryReader.jl`'s per-timestep boundary update were validated standalone
but never wired into an actual driver loop. Each step: interpolate boundary
conditions at the current `jday` (`BoundaryReader.update_boundary_
conditions!`), run one reduced-physics hydrodynamic step (QIN/QOT/QDTR
coupling, see `Hydrodynamics/FreeSurface.jl`), populate `g.TSS` with the
reduced-physics QIN/QOT/QDT heat sources (`Hydrodynamics/Transport.jl`'s
`apply_temperature_sources!`), solve temperature transport, then the real
`temperature.F90`/`update.F90` old<-new swap for temperature
(`g.HYD[:,:,4] .= g.T1`, mirroring the `H2 .= H1` etc. swap already inside
`hydrodynamic_step!`). Writes TSR rows including the `T2(C)` column (see
`IO/OutputWriter.jl`).

`boundary` is `IO/BoundaryReader.jl`'s `load_boundary_conditions` return
value. Caller must have already run `init_geometry!`, `compute_dlxrho!`,
`allocate_hydro_state!`, `allocate_transport_state!`, `compute_sf1x!`, and
set `geom.THETA`/`geom.UPWIND`/`geom.ULTIMATE` (Tier 1, not read by
`InputReader.jl` yet -- see `Hydrodynamics/Transport.jl`'s
`temperature_transport!` docstring) -- this function does not allocate
transport state itself, matching `run_zero_flow_sanity_check!`'s existing
"caller allocates, this function only steps" convention.

Uses `Hydrodynamics/AdaptiveTimestep.jl`'s `step_hydrodynamics_adaptive!`
as of 2026-08-24 (previously a single fixed `dlt = tc.DLTMAX[1]`) --
motivated by a REAL, found-not-guessed instability: once `Hydrodynamics/
Turbulence.jl`'s real TKE closure went in, its correctly-computed (given
no wind/friction forcing yet, see that file's module docstring) near-floor
`DZ` exposed that the explicit UPWIND advection scheme is NOT
unconditionally stable at `DLTMAX=200s` -- `T2` blew up to `>99°C` around
`jday≈4.125` in a real DET-forced run before this was wired in. The
previous placeholder constant `DZ` had been accidentally providing enough
numerical damping to mask this. `AdaptiveTimestepState` is created and
owned internally here (unlike hydro/transport state, which the caller
allocates) since it's intrinsic to this function's own timestepping loop,
not shared with any other driver. `constituent_transport!` is NOT called
here -- temperature only, matching this port's current Tier-1 boundary IO
scope (QIN/TIN/QOT/QDTR/TDTR, no `CIN` constituent loading yet).

`met` is `IO/MetReader.jl`'s `load_met_conditions` return value
(`Dict{Int,MetSeries}`, keyed by waterbody), optional (default `nothing`,
2026-09-10) -- when given, `MetReader.update_met_conditions!` runs each
step BEFORE `hydrodynamic_step!`/`step_hydrodynamics_adaptive!` (which
call `Hydrodynamics/Turbulence.jl`'s `compute_wind_stress!`/
`calculate_tke!` internally, so `g.WIND` must already be current for this
step). Omitting `met` leaves `g.WIND` at its safe zero default (no wind
forcing) -- same "off unless a caller explicitly loads real data"
discipline as `boundary`'s QIN/QOT/QDTR.

Returns `(g, geom, jday_final)`.
"""
function run_forced_simulation!(g, geom, net, tc, boundary; nsteps::Int, output_dir::AbstractString,
                                 output_segments::Vector{Int}, base_name::AbstractString="tsr",
                                 met=nothing)
    state = init_adaptive_timestep(tc)
    dlt = state.dlt
    jday = tc.TMSTRT

    writer = OutputWriter.open_tsr_files(output_dir, base_name, output_segments; include_temp=true)
    try
        OutputWriter.write_tsr_row!(writer, g, geom, jday, dlt)  # initial condition
        for _ in 1:nsteps
            BoundaryReader.update_boundary_conditions!(g, boundary, jday)
            met !== nothing && MetReader.update_met_conditions!(g, met, jday)
            accepted_dlt, _ = step_hydrodynamics_adaptive!(g, geom, net, tc, state, jday)
            apply_temperature_sources!(g, geom)
            temperature_transport!(g, geom, accepted_dlt)
            g.HYD[:, :, 4] .= g.T1
            jday += accepted_dlt / 86400.0
            OutputWriter.write_tsr_row!(writer, g, geom, jday, accepted_dlt)
        end
    finally
        OutputWriter.close_tsr_files!(writer)
    end
    return (g, geom, jday)
end
