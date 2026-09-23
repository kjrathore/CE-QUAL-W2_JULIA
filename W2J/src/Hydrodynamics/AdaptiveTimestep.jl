# ==============================================================================
# Hydrodynamics/AdaptiveTimestep.jl
#
# Full-fidelity port of the real adaptive-timestep machinery (user's explicit
# choice via AskUserQuestion, 2026-08-23, over a reduced check-then-step
# alternative): two real pieces, both traced directly from source, not
# guessed --
#
# 1. `update.F90:138-163` -- schedule-based DLTMAX/DLTF breakpoint lookup
#    (grows DLT after a successful step, capped by a per-JDAY schedule read
#    from w2_con.csv, `tc.DLTD`/`tc.DLTMAX`/`tc.DLTF`/`tc.DLTINTER` --
#    already read by IO/InputReader.jl). `advance_dlt_schedule!`.
# 2. `w2_4_win.f90`'s "Task 2.2.6: Autostepping" -- a CFL-based stability
#    check computed AFTER taking a step; if the step was too big, the real
#    Fortran GOTOs back to a saved state snapshot and retries with a
#    smaller DLT. `snapshot_hydro_state`/`restore_hydro_state!`/
#    `compute_curmax`/`step_hydrodynamics_adaptive!` port this exactly,
#    including the retry, not a "check the CFL limit before stepping and
#    cap DLT proactively" shortcut that would avoid needing state
#    snapshotting -- the user explicitly chose full fidelity here.
#
# REDUCED PHYSICS in `compute_curmax` (flagged, not silently narrowed):
# - `TAU1`/`TAU2` (viscosity-limit terms, real formula needs `AX(JW)`
#   Tier-1 IO and `AZ` from `Turbulence.jl`, neither ported) are always 0
#   here -- a narrower (less conservative) stability bound than real
#   Fortran's, since dropping a positive term from the denominator makes
#   `DLTCAL` LARGER (less restrictive), not smaller.
# - `CELRTY` uses `DEPTHB(KB(I),I)` (the real bottom active layer), not the
#   real formula's `DEPTHB(KBI(I),I)` -- `KBI` (the one-layer-adjustment
#   snapshot) is allocated by `Core/InitGeometry.jl` but never populated,
#   an existing gap this port inherited, not introduced here.
# - `CELERITY_LIMIT(JW)`/`VISCOSITY_LIMIT(JW)` (Tier-1 per-waterbody flags
#   gating whether these terms apply at all) aren't read -- celerity is
#   ALWAYS included here (the more conservative choice: a stricter bound
#   is safe, just not necessarily matching every real control file's exact
#   intent).
# - The real formula's `W`-based term (`(|W(K,I)|*BB(K,I)+...)*DLX(I)`) is
#   omitted since `W` (vertical velocity) is always 0 in this port
#   (continuity-derived, not ported) -- dropping an always-zero term is not
#   a simplification, just skipped algebra.
#
# NOT PORTED: the real `220 CONTINUE`/`230 CONTINUE` GOTO targets also
# handle a genuine RUNTIME ERROR case (negative surface-layer thickness at
# `DLT <= DLTMIN`, i.e. even the minimum timestep is unstable) by writing
# diagnostics and halting the whole simulation -- this port's
# `step_hydrodynamics_adaptive!` instead just accepts `DLTMIN` and moves on
# (a real difference, flagged: this could silently produce a bad step under
# genuinely severe forcing that would have been a hard error in the real
# model, rather than erroring loudly itself -- port the hard-error path
# before trusting a run that hits `DLTMIN` repeatedly).
# ==============================================================================

const NONZERO_EPS = 1.0e-20   # w2modules.F90:100 "NONZERO=1.0D-20"

"""
    HydroSnapshot

Pre-step copies of every array `hydrodynamic_step!` mutates, needed to
retry a timestep from the correct baseline after a CFL-instability
detection. CRITICAL: includes `H2`/`BH2`/`BHR2`/`AVH2` (the "old" geometry
arrays), not just `H1`/`BH1`/... (the "new" ones) -- `hydrodynamic_step!`'s
own `H2 .= H1` swap at the top of a retried call would otherwise re-swap
from the FAILED step's `H1`, corrupting the "old" baseline for the retry.

Also includes `TKE`/`AZ`/`DZ`/`AZT` (`Hydrodynamics/Turbulence.jl`'s
`calculate_tke!`, added 2026-08-24, AFTER this struct was first written --
a real gap found while wiring adaptive stepping into `Simulation.jl`'s
`run_forced_simulation!` the same day, not by a failing test): `calculate_
tke!` runs inside `hydrodynamic_step!` and its TKE update is an explicit
`TKE[k,i,1] += dlt*(...)` INCREMENT, not a fresh recomputation -- without
snapshotting/restoring it, a failed attempt's partial TKE advance would
silently carry into the retry attempt and get double-counted, corrupting
`AZ`/`DZ` (and therefore the CFL check itself) on every retry.

Also includes `KTI` (`g.KTI`, `Vector{Int}`) -- added 2026-09-22 alongside
`Hydrodynamics/FreeSurface.jl`'s `recompute_top_layer_geometry!` KTI
crossing-adjustment port, same class of gap as the TKE one above: `KTI`
was static (never mutated) when this struct was first written, so it
wasn't snapshotted; now `recompute_top_layer_geometry!` can genuinely
advance it mid-timestep, and a failed retry attempt's advance must be
undone before the next attempt, or the retry would start from the WRONG
tracked sub-layer (and, via the `Z`-rescaling that accompanies a `KTI`
change, a wrong `Z` too -- though `Z` itself is already snapshotted above).
"""
struct HydroSnapshot
    Z::Vector{Float64}; ELWS::Vector{Float64}
    U::Matrix{Float64}
    H1::Matrix{Float64}; H2::Matrix{Float64}
    BH1::Matrix{Float64}; BH2::Matrix{Float64}
    BHR1::Matrix{Float64}; BHR2::Matrix{Float64}
    AVH1::Matrix{Float64}; AVH2::Matrix{Float64}
    VOL::Matrix{Float64}
    BI::Matrix{Float64}; BKT::Vector{Float64}
    TKE::Array{Float64,3}; AZ::Matrix{Float64}; DZ::Matrix{Float64}; AZT::Matrix{Float64}
    KTI::Vector{Int}
end

"""
    snapshot_hydro_state(g, geom) -> HydroSnapshot

Copies every array `hydrodynamic_step!` mutates. Call once BEFORE the
first attempt at a new timestep (not before every retry -- the snapshot
represents the true pre-timestep baseline, reused across every retry of
the SAME timestep).
"""
function snapshot_hydro_state(g, geom)
    HydroSnapshot(copy(geom.Z), copy(geom.ELWS), copy(g.U),
                  copy(geom.H1), copy(geom.H2), copy(geom.BH1), copy(geom.BH2),
                  copy(geom.BHR1), copy(geom.BHR2), copy(geom.AVH1), copy(geom.AVH2),
                  copy(g.VOL), copy(geom.BI), copy(geom.BKT),
                  copy(g.TKE), copy(g.AZ), copy(g.DZ), copy(g.AZT),
                  copy(g.KTI))
end

"""
    restore_hydro_state!(g, geom, s::HydroSnapshot)

Restores `g`/`geom` to the state captured by `snapshot_hydro_state`. Call
before EVERY retry attempt (including the first), matching the real
Fortran's `220 CONTINUE` restore block (`Z=SZ`, `ELWS=SELWS`, `U=SU`, ...).
"""
function restore_hydro_state!(g, geom, s::HydroSnapshot)
    geom.Z .= s.Z; geom.ELWS .= s.ELWS
    g.U .= s.U
    geom.H1 .= s.H1; geom.H2 .= s.H2
    geom.BH1 .= s.BH1; geom.BH2 .= s.BH2
    geom.BHR1 .= s.BHR1; geom.BHR2 .= s.BHR2
    geom.AVH1 .= s.AVH1; geom.AVH2 .= s.AVH2
    g.TKE .= s.TKE; g.AZ .= s.AZ; g.DZ .= s.DZ; g.AZT .= s.AZT
    g.VOL .= s.VOL
    geom.BI .= s.BI; geom.BKT .= s.BKT
    g.KTI .= s.KTI
    return (g, geom)
end

"""
    compute_curmax(g, geom, dlt) -> Float64

w2_4_win.f90's "Task 2.2.6: Autostepping" CFL-based stability limit -- the
smallest `DLTCAL` (locally stable timestep) across every active cell:

    QTOT(K,I) = (|U(K,I)|*BHR1(K,I) + |U(K,I-1)|*BHR1(K,I-1)
                + DLX(I)*|BH2(K,I)-BH1(K,I)|/DLT + |QSS(K,I)|) * 0.5
    CELRTY    = sqrt(|RHO(KB(I),I)-RHO(KT,I)| / 1000 * G * DEPTHB(KB(I),I) * 0.5)
    DLTCAL    = 1 / ((QTOT(K,I)/BH1(K,I) + CELRTY)/DLX(I) + NONZERO)

See module docstring for the reduced-physics simplifications (no TAU1/TAU2
viscosity terms, no W term, CELRTY uses KB not KBI, celerity always
included). Returns `Inf` if no active cell exists (avoids `minimum` over
an empty collection erroring -- shouldn't happen for a real grid).

**`QOT` addition (2026-09-23, this port's OWN scoping call, NOT in real
Fortran)** -- confirmed by reading `w2_4_win.f90:1454` directly: the real
`QTOT` formula never includes `QOUT`/`QOT` at all, relying instead on the
real `DOWNSTREAM_WITHDRAWAL` selective-withdrawal algorithm spreading a
branch's total outflow across multiple layers (by outlet elevation +
density), which never concentrates enough flow into any ONE cell to need
a separate CFL term. THIS port's reduced-physics `QOT` (a per-branch
scalar summed into a SINGLE bottom-active-layer cell, see
`Hydrodynamics/FreeSurface.jl`'s `solve_branch_free_surface!` module
docstring) creates exactly that concentration this real formula was never
designed to bound -- found via a real full-2017 DET run where the
withdrawal cell's explicit temperature-source term
(`Hydrodynamics/Transport.jl`'s `apply_temperature_sources!`,
`TSS[kb,id] -= QOT[jb]*cold[kb,id]`) grew unboundedly and overflowed to
`-Inf` after ~50 days, BECAUSE `dlt` was never shrunk to respect that
cell's actual volume-turnover rate (`QOT[jb]*dlt` exceeding the cell's own
volume in a single explicit step is an unconditional explicit-Euler
instability, independent of any hydrodynamic CFL condition). Added as an
extra term in `qtot` at exactly the branch's `DN_FLOW` outlet cell
(`i==id`, `k==kb`), analogous to the existing `|QSS(K,I)|` term -- same
units, same role (a concentrated volumetric flow the cell must be able to
"turn over" within `dlt`), just accounting for a flow this port routes
through `QOT` instead of `QSS`.
"""
function compute_curmax(g, geom, dlt)
    curmax = Inf
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            iu, id = g.CUS[jb], g.DS[jb]
            for i in iu:id
                kb = g.KB[i]
                celrty = sqrt(abs(g.RHO[kb, i] - g.RHO[kt, i]) / 1000.0 * G_GRAVITY * geom.DEPTHB[kb, i] * 0.5)
                for k in kt:kb
                    qtot = (abs(g.U[k, i]) * geom.BHR1[k, i] + abs(g.U[k, i-1]) * geom.BHR1[k, i-1] +
                            geom.DLX[i] * abs(geom.BH2[k, i] - geom.BH1[k, i]) / dlt +
                            abs(g.QSS[k, i])) * 0.5
                    (i == id && k == kb && g.DN_FLOW[jb]) && (qtot += abs(g.QOT[jb]) * 0.5)
                    dltcal = 1.0 / ((qtot / geom.BH1[k, i] + celrty) / geom.DLX[i] + NONZERO_EPS)
                    dltcal < curmax && (curmax = dltcal)
                end
            end
        end
    end
    return curmax
end

"""
    AdaptiveTimestepState

Persistent adaptive-timestep state carried across timesteps by the caller
(real Fortran's module-level `SAVE`d `DLTDP`/`DLTMAXX`/`DLTFF`). `dltdp` is
the current 1-based index into `tc.DLTD`/`DLTMAX`/`DLTF`; `dltmaxx`/`dltff`
are the current schedule ceiling/safety-factor (possibly interpolated, see
`advance_dlt_schedule!`); `dlt` is the timestep to use for the NEXT step.
"""
mutable struct AdaptiveTimestepState
    dltdp::Int
    dltmaxx::Float64
    dltff::Float64
    dlt::Float64
end

"""
    init_adaptive_timestep(tc) -> AdaptiveTimestepState

Starts at the first schedule breakpoint (`tc.DLTD[1]`/`DLTMAX[1]`/`DLTF[1]`)
with `dlt = tc.DLTMAX[1]` -- matches real Fortran's initial condition.
"""
init_adaptive_timestep(tc) = AdaptiveTimestepState(1, tc.DLTMAX[1], tc.DLTF[1], tc.DLTMAX[1])

"""
    advance_dlt_schedule!(state, tc, jday)

`update.F90:152-160` -- advances `state.dltdp` to the next breakpoint once
`jday` reaches it, then (if `tc.DLTINTER` is on and a next breakpoint still
exists) linearly interpolates `dltmaxx`/`dltff` between the current and
next breakpoint. REDUCED: when `dltdp` is already at the LAST breakpoint
(`tc.NDLT`), real Fortran's `DLTD(DLTDP+1)` would read past the end of a
`NDLT`-length array (`IO/InputReader.jl` sizes `tc.DLTD`/`DLTMAX`/`DLTF`
to exactly `NDLT`, no sentinel padding) -- this port instead holds
`dltmaxx`/`dltff` at the last breakpoint's values forever once reached, a
safe plateau rather than an out-of-bounds guess. Correct for every real
control file validated so far (Detroit/DET both have `NDLT=1`, so this
plateau is the ONLY behavior ever exercised) -- flagged for a future
multi-breakpoint file, not assumed correct there.
"""
function advance_dlt_schedule!(state::AdaptiveTimestepState, tc, jday::Real)
    n = tc.NDLT
    if state.dltdp < n && jday >= tc.DLTD[state.dltdp+1]
        state.dltdp += 1
        state.dltmaxx = tc.DLTMAX[state.dltdp]
        state.dltff = tc.DLTF[state.dltdp]
    end
    if state.dltdp < n && tc.DLTINTER
        j1, j2 = tc.DLTD[state.dltdp], tc.DLTD[state.dltdp+1]
        state.dltmaxx = tc.DLTMAX[state.dltdp] + (tc.DLTMAX[state.dltdp+1] - tc.DLTMAX[state.dltdp]) / (j2 - j1) * (jday - j1)
        state.dltff = tc.DLTF[state.dltdp] + (tc.DLTF[state.dltdp+1] - tc.DLTF[state.dltdp]) / (j2 - j1) * (jday - j1)
    end
    return state
end

"""
    step_hydrodynamics_adaptive!(g, geom, net, tc, state, jday) -> (accepted_dlt, next_dlt)

Full step-detect-rollback-retry loop, matching the real GOTO 220 pattern:
snapshots state once, then repeatedly restores + calls `hydrodynamic_step!`
+ checks `compute_curmax` against the just-used `dlt`, halving toward the
CFL-safe value (`state.dltff * curmax`, matching `update.F90:145`'s
`DLT = DMAX1(DLTMIN,DLTFF*CURMAX)`) and retrying until stable or `dlt`
bottoms out at `tc.DLTMIN`.

Once accepted, computes `next_dlt` for the FOLLOWING step: grown to the
CFL ceiling implied by this step's `curmax`, capped at `1.1x` the accepted
`dlt` (matching `update.F90:146`'s `DLT = DMIN1(DLT,1.1*DLTS)` -- limits
how fast `DLT` can grow, not just how far it can fall), advances the
schedule (`advance_dlt_schedule!`), then caps at the (possibly just-updated)
schedule ceiling `dltmaxx`.

Mutates `state` in place (both `dlt` and the schedule fields). Returns
`(accepted_dlt, next_dlt)` -- `accepted_dlt` is what the caller should pass
to any output-writing / boundary-condition-interpolation step for the
timestep just taken; `next_dlt` equals `state.dlt` after the call (returned
separately for a caller that wants it without re-reading `state`).
"""
function step_hydrodynamics_adaptive!(g, geom, net, tc, state::AdaptiveTimestepState, jday::Real)
    snap = snapshot_hydro_state(g, geom)
    dlt = state.dlt
    curmax = dlt
    while true
        restore_hydro_state!(g, geom, snap)
        hydrodynamic_step!(g, geom, net, dlt)
        curmax = compute_curmax(g, geom, dlt)
        (curmax >= dlt || dlt <= tc.DLTMIN) && break
        dlt = max(tc.DLTMIN, state.dltff * curmax)
    end
    accepted_dlt = dlt

    next_dlt = max(tc.DLTMIN, state.dltff * curmax)
    next_dlt = min(next_dlt, 1.1 * accepted_dlt)
    advance_dlt_schedule!(state, tc, jday)
    next_dlt = min(next_dlt, state.dltmaxx)
    state.dlt = next_dlt
    return (accepted_dlt, next_dlt)
end
