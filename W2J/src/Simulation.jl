# ==============================================================================
# Simulation.jl
#
# First-cut time-stepping driver, tying together Core/InitGeometry.jl,
# Hydrodynamics/FreeSurface.jl, and IO/OutputWriter.jl into a runnable
# end-to-end simulation. See CLAUDE.md "MVP hydrodynamic run" steps 5-7 --
# this is the "zero-flow sanity check" scope confirmed with the user
# (2026-08-12): no inflow/outflow, reduced physics (see
# Hydrodynamics/FreeSurface.jl module docstring for exactly what's real vs.
# stubbed), goal is a stable free-surface over many timesteps with real TSR
# CSV output, not a faithful port of `w2_4_win.f90`'s full driver.
#
# NOT PORTED (flagged, not guessed away):
# - The real adaptive timestep machinery (DLTF/DLTMIN/DLTD breakpoints,
#   stability-based DLT selection, `w2_4_win.f90`'s `AUTO_STEPPING`) -- this
#   driver uses one FIXED `dlt = tc.DLTMAX[1]` for every step. Correct for a
#   reduced-physics zero-flow run (nothing destabilizes the timestep), wrong
#   once real forcing (structures, meteorology, non-uniform density) is
#   added -- port real DLT selection before trusting a non-sanity-check run.
# - Any boundary condition IO (inflow/outflow/withdrawal/tributary time
#   series) -- Tier 1, not built (CLAUDE.md "Open questions" / MVP step 8).
#   `g.QSS`/`g.UXBR`/`g.UYBR` stay at their allocated zero forever in this
#   driver.
# - Kinetics / constituent transport (WaterQuality/*, Hydrodynamics/
#   Transport.jl) -- both still stubs, not called here.
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
