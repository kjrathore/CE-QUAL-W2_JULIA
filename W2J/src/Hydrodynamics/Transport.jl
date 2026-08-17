# ==============================================================================
# Hydrodynamics/Transport.jl
#
# First-cut constituent transport (temperature + one generic constituent),
# ported from transport.f90 (explicit advective flux terms) + temperature.F90
# (implicit vertical-diffusion assembly and solve) + wqconstituents.F90 (the
# same assembly, reused for other constituents). See CLAUDE.md for the
# tracing session that found transport.f90 does NOT itself update
# concentrations -- its four ENTRY points (INTERPOLATION_MULTIPLIERS,
# HORIZONTAL_MULTIPLIERS1/HORIZONTAL_MULTIPLIERS, VERTICAL_MULTIPLIERS1/
# VERTICAL_MULTIPLIERS) only compute the explicit advective/dispersive flux
# arrays ADX(K,I)/ADZ(K,I); temperature.F90 and wqconstituents.F90 are the
# real callers that build the implicit AT/CT/VT/DT tridiagonal system and
# solve it (temperature.F90:540-590, wqconstituents.F90:539-584).
#
# SCOPE -- REDUCED PHYSICS, BY DELIBERATE CHOICE (confirmed with user,
# 2026-08-15), same discipline as Hydrodynamics/FreeSurface.jl:
#
#   REAL (computed from actual formulas, traced from source, not guessed):
#   - SF1X(I): transport.f90's INTERPOLATION_MULTIPLIERS, the one piece of
#     it actually needed under UPWIND -- see State.jl's SF1X docstring for
#     why it collapsed from KMX x IMX to a per-segment vector.
#   - ADX(K,I)/ADZ(K,I): transport.f90's HORIZONTAL_MULTIPLIERS/VERTICAL_
#     MULTIPLIERS, UPWIND branch only (`IF (UPWIND(JW))` in the real
#     source) -- the ULTIMATE/QUICKEST 3rd-order branch (SF2X..SF13X,
#     RATD/CURX1-3, RATZ/CURZ1-3, ~400 more lines) is NOT ported. Callers
#     MUST set `geom.ULTIMATE[jw] = false` (see W2Geometry's THETA/UPWIND/
#     ULTIMATE docstring) -- this file does not silently fall back, it
#     errors if asked to run with ULTIMATE(jw) true.
#   - AT/CT/VT/DT: temperature.F90:560-573's real implicit vertical-
#     diffusion tridiagonal assembly (THETA-blended Crank-Nicolson-style),
#     solved with `Solvers/Tridiagonal.jl`'s `thomas_solve!` instead of the
#     inlined BTA1/GMA1 Thomas algorithm the Fortran source repeats at
#     every call site (transport.f90's own TRIDIAG subroutine is a third,
#     unused-here copy of the identical algorithm) -- this project's
#     existing "shared numerical primitives ported once" rule
#     (Solvers/Tridiagonal.jl's module docstring), not a new decision.
#   - REAL SEQUENCING DEPENDENCY, traced not assumed: AT/CT/VT are NOT
#     recomputed per constituent -- wqconstituents.F90:555-556 calls
#     HORIZONTAL_MULTIPLIERS/VERTICAL_MULTIPLIERS again (constituent-
#     specific ADX/ADZ) but never recomputes AT/CT/VT, confirmed by reading
#     wqconstituents.F90:562-580 directly: it reuses whatever temperature.F90
#     last wrote. `constituent_transport!` in this file does the same --
#     it must be called in the same timestep AFTER `temperature_transport!`,
#     which is the only function that calls `compute_vertical_diffusion_
#     coefficients!`. Not enforced by a runtime assertion (no cheap way to
#     check "was called after" without extra state) -- documented instead,
#     matching how `branch_processing_order`'s dependency is a traced
#     requirement, not a guess.
#
#   STUBBED, EXPLICITLY ZERO/CONSTANT, FLAGGED (not silently omitted):
#   - W (vertical velocity): NOT YET COMPUTED anywhere in this port --
#     w2_4_win.f90's continuity-derived vertical-velocity update hasn't
#     been traced or ported. Always 0 here (matches the zero-flow scope
#     this whole port has been validated against so far).
#   - DZ (vertical eddy diffusivity): needs Hydrodynamics/Turbulence.jl
#     (still a stub). Reduced-physics: caller sets a uniform constant via
#     `allocate_transport_state!`'s `dz_const` keyword -- NOT zero by
#     default, unlike W/SB/ST, because a zero DZ would make the vertical-
#     diffusion solve trivially degenerate (AT=CT=0, VT=1, i.e. a no-op) --
#     see VALIDATED below for why a nonzero DZ was chosen deliberately.
#   - DX (horizontal dispersion coefficient): needs DXI (Tier 1, not read
#     by InputReader.jl yet). Reduced-physics: caller-set uniform constant,
#     same mechanism as DZ. Zero is a valid (if physically inert) choice
#     here, unlike DZ.
#   - ULTIMATE/QUICKEST advection (see above) -- UPWIND only.
#
# VALIDATED (test/runtests.jl, "Hydrodynamics/Transport"): under the existing
# zero-flow scenario (U=W=0 everywhere), horizontal and vertical ADVECTION
# are exactly zero regardless of DZ/DX (U=W=0 makes ADX/ADZ's formulas
# collapse to 0*anything), so this alone would be as degenerate a test as
# FreeSurface.jl's zero-flow check. The IMPLICIT vertical-diffusion solve
# is NOT degenerate under zero flow, though -- AT/CT/VT depend on DZ and
# geometry, not on U/W (the THETA*W terms vanish, but the DZ terms do not).
# So a nonzero DZ plus a non-uniform initial temperature profile (warm
# surface, cold bottom -- physically realistic, not an arbitrary choice)
# gives a genuine test: temperature should diffuse toward uniform over
# time, conserve total heat (mass-like conservation check), and never go
# unstable/NaN. This is a real test of the AT/CT/VT assembly and the
# thomas_solve! integration, not a tautology.
#
# TWO REAL BUGS FOUND AND FIXED VIA THIS CONSERVATION CHECK (2026-08-15,
# not hypothetical edge cases -- see `allocate_transport_state!`'s
# docstring for the full derivation): `geom.BB` is genuinely nonzero both
# above the current top active layer AND at the lakebed-interface row of
# the bottom active layer (raw bathymetry width, independent of where the
# model currently truncates the water column) -- a naively uniform `DZ`
# fill coupled the top/bottom active layers to phantom cells above the
# surface and below the bed, leaking heat. First version (whole-array
# fill): total collapse (a 200-step run drove a spike profile to near-
# uniform ~0.05 in ~200 steps, an order of magnitude too fast). Second
# version (bounded to KTWB:KB[i], but including KB[i] itself): much
# better but still a real, compounding ~24% loss over 3000 steps for a
# sharp step profile in a single segment. Final version (bounded to
# KTWB:KB[i]-1, excluding the lakebed-interface row): full-domain
# single-step heat conservation to machine precision (`~1e-16` relative),
# 3000-step conservation to `~1e-13` relative, with the step profile
# smoothing exactly as physical diffusion should (std `7.6 -> 0.0002`,
# converging to the volume-weighted mean temperature). See
# `test/runtests.jl`, "Hydrodynamics/Transport", for the codified version
# of this check.
# ==============================================================================

"""
    allocate_transport_state!(g; dz_const=1e-5, dx_const=0.0)

Sizes the transport solve arrays added to `W2Global` for this file (W, DZ,
DX, SF1X, ADX, ADZ, AT, CT, VT, DT), plus `T1`/`TSS`/`HYD` -- declared in
`Core/State.jl` since the original session but never allocated by any
caller until this file needed them (`T1` is the new-timestep temperature,
`HYD[:,:,4]` the old-timestep temperature per `temperature.F90:541`'s
`COLD => HYD(:,:,4)`, `TSS` the source term). Call after `IO/InputReader.
allocate_geometry!`, before any `*_transport!` call.

`dz_const`/`dx_const` set the REDUCED-PHYSICS uniform stand-ins for DZ
(vertical eddy diffusivity) and DX (horizontal dispersion) -- see module
docstring for why DZ defaults nonzero (`1e-5 m^2/s`, a physically modest
molecular-diffusion-scale value, deliberately small so the implicit solve
stays well-behaved without needing a real Turbulence.jl closure) while DX
defaults to zero (physically inert under this port's zero horizontal-flow
scope, since ADX's only other term is U-weighted and U=0 here too).

`dz_const` is filled ONLY for `KTWB[jw] <= k <= KB[i]-1`, one SEGMENT at a
time -- deliberately NOT the whole KMX x IMX array, and NOT even the full
`KTWB[jw]:KB[i]` active range. `DZ[k, i]` represents the diffusivity at the
BOTTOM interface of layer `k` (the coupling between layer `k` and `k+1`) --
so it is only meaningful for interior interfaces, `KTWB[jw] <= k <=
KB[i]-1`; `DZ[KB[i], i]` (the interface below the LAST active layer, i.e.
"through the lakebed") must be 0.

This was found empirically, not assumed, across two rounds of debugging
during validation:

  1. First bug: blanket-filled `DZ` across the whole `KMX x IMX` array.
     `geom.BB[k, i]` for `k < KTWB[jw]` is genuinely nonzero (raw bathymetry
     width at that depth, independent of current water level -- confirmed
     against real Detroit data, e.g. `BB[KTWB-1, i]` = `17.5`/`479.5`/
     `819.8` across different segments, not a rounding artifact), so a
     nonzero `DZ` there spuriously coupled the top active layer to a
     phantom cell above the free surface. Fixed by bounding the fill to
     `KTWB[jw]:KMX`.
  2. Second bug, found by a heat-conservation check that still showed a
     real ~24% loss over 3000 steps for a sharp step-function profile even
     after fix 1: `geom.BB[KB[i], i]` (the lakebed-interface width) is ALSO
     genuinely nonzero (confirmed: `BB[KB[i], i] = 5` for a real Detroit
     segment), so filling `DZ` all the way down to (and including) each
     segment's own `KB[i]` spuriously coupled the bottom active layer to a
     phantom cell below the lakebed. Fixed by stopping the fill at
     `KB[i]-1`.

Post-fix, a full-domain single-step heat-conservation check (weighted by
`BH2`, summed over every wet segment) matched to machine precision
(`~1e-16` relative), and a 3000-step run of a sharp step-function initial
profile conserved total heat to `~1e-13` relative while smoothing exactly
as pure vertical diffusion should (std of the profile decaying from `7.6`
to `0.0002`, converging to the volume-weighted mean) -- see
`test/runtests.jl`, "Hydrodynamics/Transport".
"""
function allocate_transport_state!(g; dz_const::Float64=1e-5, dx_const::Float64=0.0)
    kmx, imx = g.KMX, g.IMX
    g.T1 = zeros(Float64, kmx, imx)
    g.TSS = zeros(Float64, kmx, imx)
    g.HYD = zeros(Float64, kmx, imx, g.NHY)
    g.W = zeros(Float64, kmx, imx)
    g.DZ = zeros(Float64, kmx, imx)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            for i in (g.US[jb]-1):(g.DS[jb]+1)
                kb = g.KB[i]
                kb > kt && (g.DZ[kt:(kb-1), i] .= dz_const)
            end
        end
    end
    g.DX = fill(dx_const, kmx, imx)
    g.SF1X = zeros(Float64, imx)
    g.ADX = zeros(Float64, kmx, imx)
    g.ADZ = zeros(Float64, kmx, imx)
    g.AT = zeros(Float64, kmx, imx)
    g.CT = zeros(Float64, kmx, imx)
    g.VT = zeros(Float64, kmx, imx)
    g.DT = zeros(Float64, kmx, imx)
    return g
end

"""
    compute_sf1x!(g, geom)

transport.f90's INTERPOLATION_MULTIPLIERS -- the K-independent SF1X piece
only (see State.jl's SF1X docstring). Static per-segment geometry, computed
once (like `Hydrodynamics/FreeSurface.jl`'s `compute_dlxrho!`), not every
timestep. Range matches `compute_horizontal_advection!`'s own loop
(`I=IU,ID-1`, transport.f90:247) -- SF1X is never read outside that range.
"""
function compute_sf1x!(g, geom)
    for jb in 1:g.NBR
        iu, id = g.CUS[jb], g.DS[jb]
        for i in iu:(id-1)
            g.SF1X[i] = (geom.DLX[i+1] + geom.DLX[i]) * 0.5
        end
    end
    return g
end

"""
    compute_horizontal_advection!(g, geom, jb, kt, cold)

transport.f90's HORIZONTAL_MULTIPLIERS1 + HORIZONTAL_MULTIPLIERS, UPWIND
branch only (lines 233-244, 286-299), algebraically combined into one pass
-- the DX1/DX2/DX3 intermediate multiplier arrays are pure algebra;
substituting them directly into ADX's formula instead of storing them as
separate fields is the same simplification style as
`Hydrodynamics/FreeSurface.jl`'s merged HDG/HPG. `cold` is the constituent's
OLD-timestep concentration field (KMX x IMX) -- `g.HYD[:,:,4]` for
temperature, `g.C1S[:,:,jc]` for a generic constituent.

Writes `g.ADX[k, i]` for `i in iu:(id-1)` only (transport.f90's own loop
range) -- `ADX[k, iu-1]` is deliberately left at whatever `fill!` set it to
(the correct upstream-boundary flux for this port's zero-inflow scope, not
an oversight -- see module docstring).
"""
function compute_horizontal_advection!(g, geom, jb, kt, cold)
    iu, id = g.CUS[jb], g.DS[jb]
    for i in iu:(id-1)
        for k in kt:g.KB[i]
            dxs = g.DX[k, i] / g.SF1X[i]
            if g.U[k, i] >= 0.0
                c2, c3 = cold[k, i], cold[k, i+1]
                g.ADX[k, i] = (-dxs - g.U[k, i]) * c2 + dxs * c3
            else
                c1, c2 = cold[k, i], cold[k, i+1]
                g.ADX[k, i] = -dxs * c1 + (dxs - g.U[k, i]) * c2
            end
        end
    end
    return g
end

"""
    compute_vertical_advection!(g, geom, jb, kt, cold)

transport.f90's VERTICAL_MULTIPLIERS, UPWIND branch only (lines 465-472).
Writes `g.ADZ[k, i]` for `i in iu:id`, `k in kt:(KB(i)-1)`.
"""
function compute_vertical_advection!(g, geom, jb, kt, cold)
    iu, id = g.CUS[jb], g.DS[jb]
    for i in iu:id
        for k in kt:(g.KB[i]-1)
            c2 = g.W[k, i] >= 0.0 ? cold[k, i] : cold[k+1, i]
            g.ADZ[k, i] = -g.W[k, i] * c2
        end
    end
    return g
end

"""
    compute_vertical_diffusion_coefficients!(g, geom, jw, jb, kt, dlt)

temperature.F90:566-576 -- the implicit vertical-diffusion tridiagonal
coefficients. NOT constituent-specific (see module docstring's "REAL
SEQUENCING DEPENDENCY" note) -- computed once per branch per timestep by
`temperature_transport!`, reused as-is by `constituent_transport!`.
"""
function compute_vertical_diffusion_coefficients!(g, geom, jw, jb, kt, dlt)
    iu, id = g.CUS[jb], g.DS[jb]
    theta = geom.THETA[jw]
    for i in iu:id
        for k in kt:g.KB[i]
            g.AT[k, i] = -dlt / geom.BH1[k, i] * (geom.BB[k-1, i] * (g.DZ[k-1, i] / geom.AVH1[k-1, i] + theta * 0.5 * g.W[k-1, i]))
            g.CT[k, i] =  dlt / geom.BH1[k, i] * (geom.BB[k, i] * (theta * 0.5 * g.W[k, i] - g.DZ[k, i] / geom.AVH1[k, i]))
            g.VT[k, i] =  1.0 + dlt / geom.BH1[k, i] * (geom.BB[k, i] * (g.DZ[k, i] / geom.AVH1[k, i] + theta * 0.5 * g.W[k, i]) +
                                                          geom.BB[k-1, i] * (g.DZ[k-1, i] / geom.AVH1[k-1, i] - theta * 0.5 * g.W[k-1, i]))
        end
    end
    return g
end

"""
    assemble_transport_rhs!(g, geom, jw, jb, kt, dlt, cold, ss_b; ss_k=nothing)

temperature.F90:558-561 (`ss_k` omitted, matching the real source -- it has
no volumetric source term) / wqconstituents.F90:562-565 (`ss_k` given,
`CSSK`) -- the explicit RHS of the transport equation. `ss_b`/`ss_k` are
KMX x IMX slices of `g.TSS` (temperature) or `g.CSSB[:,:,jc]`/
`g.CSSK[:,:,jc]` (generic constituent).
"""
function assemble_transport_rhs!(g, geom, jw, jb, kt, dlt, cold, ss_b; ss_k=nothing)
    iu, id = g.CUS[jb], g.DS[jb]
    theta = geom.THETA[jw]
    for i in iu:id
        for k in kt:g.KB[i]
            g.DT[k, i] = (cold[k, i] * geom.BH2[k, i] / dlt +
                          (g.ADX[k, i] * geom.BHR1[k, i] - g.ADX[k, i-1] * geom.BHR1[k, i-1]) / geom.DLX[i] +
                          (1.0 - theta) * (g.ADZ[k, i] * geom.BB[k, i] - g.ADZ[k-1, i] * geom.BB[k-1, i]) +
                          ss_b[k, i] / geom.DLX[i]) * dlt / geom.BH1[k, i]
            ss_k !== nothing && (g.DT[k, i] += ss_k[k, i] * dlt)
        end
    end
    return g
end

"""
    transport_solve!(new_field, g, jb, kt)

Per-column implicit solve using the AT/CT/VT/DT already assembled by
`compute_vertical_diffusion_coefficients!`/`assemble_transport_rhs!`, via
`Solvers/Tridiagonal.jl`'s `thomas_solve!` (this project's shared primitive,
not the inlined BTA1/GMA1 Thomas algorithm temperature.F90/wqconstituents.F90
each repeat -- see module docstring). Writes into `new_field[kt:KB(i), i]`
for every `i` in the branch -- `T1` for temperature, `C1[:,:,jc]` for a
generic constituent.
"""
function transport_solve!(new_field, g, jb, kt)
    iu, id = g.CUS[jb], g.DS[jb]
    for i in iu:id
        kb = g.KB[i]
        thomas_solve!((@view g.AT[kt:kb, i]), (@view g.VT[kt:kb, i]), (@view g.CT[kt:kb, i]),
                      (@view g.DT[kt:kb, i]), (@view new_field[kt:kb, i]))
    end
    return new_field
end

"""
    temperature_transport!(g, geom, dlt)

Orchestrates one timestep of temperature transport for every waterbody/
branch: horizontal + vertical advection, the implicit vertical-diffusion
coefficients (AT/CT/VT, shared with `constituent_transport!` -- see module
docstring), the RHS, and the solve into `g.T1`. Old-timestep temperature is
read from `g.HYD[:,:,4]` (`COLD => HYD(:,:,4)`, temperature.F90:541),
source term from `g.TSS`.

Requires `geom.THETA`/`geom.UPWIND`/`geom.ULTIMATE` set explicitly (Tier 1,
not read by InputReader.jl yet) -- errors loudly if any waterbody has
`ULTIMATE[jw]` true (not ported, see module docstring) or `UPWIND[jw]`
false (nothing else is implemented).
"""
function temperature_transport!(g, geom, dlt)
    fill!(g.ADX, 0.0)
    fill!(g.ADZ, 0.0)
    for jw in 1:g.NWB
        geom.ULTIMATE[jw] && error("temperature_transport!: ULTIMATE(jw=$jw) is true, but only the UPWIND advection scheme is ported (Hydrodynamics/Transport.jl module docstring) -- set geom.ULTIMATE[jw] = false")
        geom.UPWIND[jw] || error("temperature_transport!: UPWIND(jw=$jw) is false, but no other advection scheme is implemented")
        kt = g.KTWB[jw]
        cold = @view g.HYD[:, :, 4]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            compute_horizontal_advection!(g, geom, jb, kt, cold)
            compute_vertical_advection!(g, geom, jb, kt, cold)
            compute_vertical_diffusion_coefficients!(g, geom, jw, jb, kt, dlt)
            assemble_transport_rhs!(g, geom, jw, jb, kt, dlt, cold, g.TSS)
            transport_solve!(g.T1, g, jb, kt)
        end
    end
    return g
end

"""
    constituent_transport!(g, geom, dlt, jc)

Orchestrates one timestep of transport for generic constituent index `jc`
(`g.C1S[:,:,jc]` old concentration, `g.CSSB[:,:,jc]`/`g.CSSK[:,:,jc]` source
terms, solved into `g.C1[:,:,jc]`). MUST be called after
`temperature_transport!` in the same timestep -- reuses `g.AT`/`g.CT`/`g.VT`
as-is rather than recomputing them, matching wqconstituents.F90:555-580's
real behavior (traced, not assumed -- see module docstring). Only advection
(`g.ADX`/`g.ADZ`) is recomputed per constituent, since it depends on `cold`.
"""
function constituent_transport!(g, geom, dlt, jc)
    fill!(g.ADX, 0.0)
    fill!(g.ADZ, 0.0)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        cold = @view g.C1S[:, :, jc]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            compute_horizontal_advection!(g, geom, jb, kt, cold)
            compute_vertical_advection!(g, geom, jb, kt, cold)
            assemble_transport_rhs!(g, geom, jw, jb, kt, dlt, cold, (@view g.CSSB[:, :, jc]); ss_k=(@view g.CSSK[:, :, jc]))
            transport_solve!((@view g.C1[:, :, jc]), g, jb, kt)
        end
    end
    return g
end