# ==============================================================================
# Hydrodynamics/FreeSurface.jl
#
# First-cut free-surface elevation + horizontal momentum solve, ported from
# w2_4_win.f90 (the main `INTEGER FUNCTION CE_QUAL_W2` time-stepping
# function -- see Core/Grid.jl and CLAUDE.md "MVP hydrodynamic run" for why
# this isn't a separable subroutine the way init-geom.F90 was).
#
# SCOPE -- REDUCED PHYSICS, BY DELIBERATE CHOICE (confirmed with user,
# 2026-08-12), NOT a faithful full port like every other module in this
# session so far. The goal was a runnable first cut producing real TSR
# output for Detroit under a zero-boundary-flow sanity check, not full
# physics fidelity. What's real vs. stubbed:
#
# PARALLEL PROCESSING (Pillar 1, per user's explicit "then next is to
# prepare for parallel processing"). Every per-column computation in this
# file (`compute_density_field!`, `compute_pressure_field!`,
# `compute_gravity_term!`, `compute_pressure_gradient!`,
# `update_velocities!`) and `solve_free_surface!`'s waterbody/branch-tier
# loops go through `Core/Parallel.jl`'s `parallel_foreach` -- each column
# only ever writes its own data (a couple read a neighbor column, none
# write one), so this is embarrassingly parallel, matching the `TRIDIAG`
# case-A analysis in CLAUDE.md. `solve_free_surface!` is NOT naively
# threaded over all branches -- `Core/Grid.jl`'s `branch_processing_tiers`
# groups branches into dependency tiers first (see that function's
# docstring), and only branches *within* one tier run concurrently; tiers
# still execute strictly in order, since a later tier's branch genuinely
# depends on an earlier tier's freshly-solved `Z`/`ELWS` (the real, traced
# sequential dependency this file's docstring already established --
# parallelizing within a tier does not relax that).
#
# FIRST-PASS REGRESSION, FOUND AND FIXED (2026-08-14): naively wrapping
# every loop above in bare `Threads.@threads` (no size check) made
# Detroit's zero-flow sanity check SLOWER at every thread count tested --
# 0.136ms/step at 1 thread vs. 0.484ms/step at 8 threads, because Detroit's
# grid (IMX=31, at most ~31 columns and 1-4 branches per loop) is far too
# small for the per-column/per-branch work to outweigh `Threads.@threads`'s
# own spawn/scheduling overhead. `parallel_foreach` fixes this by picking
# serial vs. threaded execution per call based on the actual loop length
# against `PARALLEL_THRESHOLD` (empirically benchmarked, see
# `Core/Parallel.jl`) -- Detroit now runs every one of these loops serially
# (well under the threshold) while a larger reservoir's larger grid would
# cross it and thread automatically, with no code change needed either way.
# Re-validated after the fix: full test suite still 484/484 at
# `--threads=1` and `--threads=4`; `hydrodynamic_step!` results remain
# bit-identical across thread counts (same proof as before, re-run after
# this change, not assumed to still hold).
#
#   REAL (computed from actual formulas, not guessed):
#   - RHO(K,I): Hydrodynamics/Density.jl's `density()`, the real equation of
#     state -- see `compute_density_field!`.
#   - P(K,I): hydrostatic pressure integration (w2_4_win.f90:1190-1196).
#   - GRAV(K,I): gravity/channel-slope term (w2_4_win.f90:879-884) -- real
#     formula, though it evaluates to 0 for Detroit (SLOPE=0 every branch).
#   - The free-surface implicit tridiagonal solve itself
#     (w2_4_win.f90:899-1015), including branch-to-branch sequencing via
#     Core/Grid.jl's `branch_processing_order` (needed for real: Detroit's
#     branches 2-4 have DH_INTERNAL boundaries into branch 1, so this DOES
#     get exercised even in the zero-flow case).
#   - The explicit velocity update (w2_4_win.f90:1301-1307) -- a closed-form
#     formula, not literally a stub, just fed mostly-zero forcing (see below).
#
#   STUBBED, EXPLICITLY ZERO, FLAGGED (not silently omitted):
#   - SB, ST (bottom/wind shear): need Hydrodynamics/Turbulence.jl's AZ and
#     meteorology (wind) IO, neither built yet.
#   - ADMX, ADMZ (advection of momentum), DM (dispersion of momentum): real
#     formulas exist and only need U/geometry (see w2_4_win.f90:845-869), but
#     are 0 whenever U=0 anyway (our zero-flow starting condition) -- ported
#     as literal zero for now rather than the real (currently-inert) formula,
#     to keep this file's first cut smaller. Port these before trusting any
#     run with nonzero boundary inflow.
#   - HDG and HPG are NOT distinguished as separate old-geometry/new-geometry
#     computations the way the real source does (HDG at w2_4_win.f90:595,
#     using template H; HPG at :1201, using post-solve H1) -- this file
#     computes ONE pressure-gradient term per timestep and reuses it for
#     both the free-surface forcing and the velocity update. Exactly correct
#     for uniform density (both are 0 regardless), a real simplification
#     otherwise -- fix by splitting before trusting a non-uniform-density run.
#   - The dam-flow / reciprocal head-flow branch case (HEAD_FLOW/
#     INTERNAL_FLOW/DAM_INFLOW, only triggered when UHS(JB) != 0) -- not
#     reachable for Detroit (UHS all 0), not ported, matches the same
#     deferral already flagged in Core/InitGeometry.jl's
#     `compute_boundary_flags!`.
#   - Implicit vertical eddy viscosity correction step (w2_4_win.f90:1309-
#     1327, its own TRIDIAG call using AZ) -- skipped entirely, needs
#     Turbulence.jl.
#
# VALIDATED (test/runtests.jl, "Hydrodynamics/FreeSurface"): run against real
# Detroit data under zero boundary flow -- water surface elevation stays
# stable (no drift beyond floating-point noise) over many timesteps, exactly
# the expected outcome when all forcing terms are genuinely zero (uniform
# density + zero slope + zero flow + zero wind). This is a real test of the
# tridiagonal assembly and branch sequencing, not a tautology -- a bug in
# branch ordering, the BHRHO/A/V/C coefficients, or the DH_INTERNAL boundary
# coupling would show up as drift or NaN, not as a trivially-passing check.
# ==============================================================================

const G_GRAVITY = 9.81      # w2modules.F90:135 "G=9.81D0"
const RHOW = 1000.0         # w2modules.F90:101 "RHOW=1000.0D0" -- reference density for DLXRHO

"""
    allocate_hydro_state!(g)

Sizes the hydrodynamic solve arrays added to `W2Global` for this file (U,
RHO, P, HPG, GRAV, SB, ST, ADMX, ADMZ, DM, DLXRHO) plus QSS/UXBR/UYBR
(declared in `W2Global` since the original session but never allocated --
needed here as the always-zero source-term arrays the free-surface solve
reads), plus QIN/QIND/TIN/TIND/QOT/QDTR (IO/BoundaryReader.jl's external
inflow/outflow/distributed-tributary state, zero until `BoundaryReader.
update_boundary_conditions!` is called). Call after `IO/InputReader.
allocate_geometry!`, before `hydrodynamic_step!`.
"""
function allocate_hydro_state!(g)
    kmx, imx = g.KMX, g.IMX
    g.U = zeros(Float64, kmx, imx)
    g.RHO = zeros(Float64, kmx, imx)
    g.P = zeros(Float64, kmx, imx)
    g.HPG = zeros(Float64, kmx, imx)
    g.GRAV = zeros(Float64, kmx, imx)
    g.SB = zeros(Float64, kmx, imx)
    g.ST = zeros(Float64, kmx, imx)
    g.ADMX = zeros(Float64, kmx, imx)
    g.ADMZ = zeros(Float64, kmx, imx)
    g.DM = zeros(Float64, kmx, imx)
    g.QSS = zeros(Float64, kmx, imx)
    g.UXBR = zeros(Float64, kmx, imx)
    g.UYBR = zeros(Float64, kmx, imx)
    g.DLXRHO = zeros(Float64, imx)
    g.QIN = zeros(Float64, g.NBR)
    g.QIND = zeros(Float64, g.NBR)
    g.TIN = zeros(Float64, g.NBR)
    g.TIND = zeros(Float64, g.NBR)
    g.QOT = zeros(Float64, g.NBR)
    g.QDTR = zeros(Float64, g.NBR)
    g.TDTR = zeros(Float64, g.NBR)
    g.QDT = zeros(Float64, imx)
    g.QINF = zeros(Float64, kmx, g.NBR)
    g.KTQIN = zeros(Int, g.NBR)
    g.KBQIN = zeros(Int, g.NBR)
    # DIST_TRIBS defaults to false for every branch -- unlike THETA/UPWIND/
    # ULTIMATE (only required by opt-in functions like temperature_
    # transport!), distribute_tributary! runs unconditionally inside every
    # hydrodynamic_step! call, so leaving this unsized would BoundsError
    # every existing caller, not just ones using distributed tributaries.
    # false is the safe default (no distributed inflow unless a caller
    # explicitly opts a branch in via real loaded data), not a silent
    # accuracy gap the way defaulting ULTIMATE would be.
    g.DIST_TRIBS = fill(false, g.NBR)
    # PLACE_QIN is sized lazily inside apply_inflow_boundary! (see there) --
    # this function only takes `g`, not `geom`, so it can't size a
    # W2Geometry field itself without changing every existing call site's
    # signature (a much larger, unrelated churn for a one-line default).
    #
    # AZ/DZ/TKE/AZT: calculate_tke! also runs unconditionally inside every
    # hydrodynamic_step! call (same reasoning as DIST_TRIBS above), and
    # these fields all live on `g` (not `geom`), so -- unlike PLACE_QIN --
    # they CAN be sized here directly. Hydrodynamics/Turbulence.jl's
    # allocate_turbulence_state! also inits TKE to (TKEMIN1,TKEMIN2) and AZ
    # to AZMIN, not zero -- see that function's docstring for why (avoids a
    # 0/0 in calculate_tke!'s dissipation-ratio term on the very first
    # call).
    allocate_turbulence_state!(g)
    return g
end

"""
    compute_dlxrho!(g, geom)

init.F90:735-739 -- static geometric factor for the horizontal pressure-
gradient terms (HDG/HPG). Computed once, not per-timestep -- depends only on
DLXR (geometry) and the constant reference density RHOW, not on the
(evolving) actual density field.
"""
function compute_dlxrho!(g, geom)
    imx = g.IMX
    g.DLXRHO = zeros(Float64, imx)
    for jb in 1:g.NBR
        iu, id = g.US[jb], g.DS[jb]
        for i in iu:id
            g.DLXRHO[i] = 0.5 / (geom.DLXR[i] * RHOW)
        end
        g.UP_HEAD[jb] && (g.DLXRHO[iu-1] = 0.5 / (geom.DLXR[iu] * RHOW))
    end
    return g
end

"""
    compute_density_field!(g, geom)

RHO(K,I) via Hydrodynamics/Density.jl's `density()`. Uses each waterbody's
initial temperature `T2I` UNIFORMLY across all its segments and layers (no
temperature transport yet -- Hydrodynamics/Transport.jl is still a stub, so
there's no mechanism for T to actually evolve or vary spatially). TDS=SS=0
throughout (no constituent transport yet either), so the `fresh_water`/
`salt_water` distinction is numerically inert here (multiplies a zero) --
passed as `true`/`false` for concreteness, not because it's been derived via
`compute_water_type_flags!` for this reduced path.

PARALLEL PROCESSING: the per-segment loop is threaded -- each `i` only ever
writes its own column `g.RHO[:, i]`, confirmed embarrassingly parallel
(Pillar 1), same reasoning as `Solvers/Tridiagonal.jl`'s per-column solves.
"""
function compute_density_field!(g, geom)
    kmx, imx = g.KMX, g.IMX
    g.RHO = zeros(Float64, kmx, imx)
    for jw in 1:g.NWB
        T = geom.T2I[jw]
        rho_val = density(T, 0.0, 0.0, true, false, false)
        for jb in g.BS[jw]:g.BE[jw]
            parallel_foreach((g.US[jb]-1):(g.DS[jb]+1)) do i
                g.RHO[:, i] .= rho_val
            end
        end
    end
    return g
end

"""
    compute_pressure_field!(g, geom)

w2_4_win.f90:1190-1196 -- hydrostatic pressure by cumulative summation down
from the surface layer. `P(KT-1,I)` (i.e. above the water surface) is 0.

PARALLEL PROCESSING: threaded over segments `i` -- the `k` recurrence
(`P[k,i] = P[k-1,i] + ...`) is sequential DOWN one column but every column
is independent of every other, same as `Solvers/Tridiagonal.jl`'s per-
column solves (Pillar 1).
"""
function compute_pressure_field!(g, geom)
    kmx, imx = g.KMX, g.IMX
    g.P = zeros(Float64, kmx, imx)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            parallel_foreach((g.US[jb]-1):(g.DS[jb]+1)) do i
                kt > g.KMX && return
                g.P[kt, i] = g.RHO[kt, i] * G_GRAVITY * geom.H1[kt, i] * geom.COSA[jb]
                for k in (kt+1):g.KB[i]
                    g.P[k, i] = g.P[k-1, i] + g.RHO[k, i] * G_GRAVITY * geom.H1[k, i] * geom.COSA[jb]
                end
            end
        end
    end
    return g
end

"""
    compute_gravity_term!(g, geom)

w2_4_win.f90:879-884 -- gravity force from channel slope. Real formula;
evaluates to 0 for Detroit since every branch has SLOPE=0 (SINAC=0).

PARALLEL PROCESSING: threaded over segments `i`, each writing only its own
column -- embarrassingly parallel (Pillar 1).
"""
function compute_gravity_term!(g, geom)
    kmx, imx = g.KMX, g.IMX
    g.GRAV = zeros(Float64, kmx, imx)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            parallel_foreach((g.US[jb]-1):g.DS[jb]) do i
                g.GRAV[kt, i] = geom.AVHR[kt, i] * (geom.BKT[i] + geom.BKT[i+1]) * 0.5 * G_GRAVITY * geom.SINAC[jb]
                for k in (kt+1):g.KB[i]
                    g.GRAV[k, i] = geom.BHR2[k, i] * G_GRAVITY * geom.SINAC[jb]
                end
            end
        end
    end
    return g
end

"""
    compute_pressure_gradient!(g, geom)

w2_4_win.f90:1198-1205 (HPG) / :593-598 (HDG) -- horizontal pressure
gradient. See module docstring: this file does NOT distinguish HDG (pre-
free-surface-solve, template geometry) from HPG (post-solve, updated
geometry) as two separate computations -- one `HPG` field is computed and
reused for both purposes. Exactly correct when density is spatially
uniform (both are 0 regardless); a documented simplification otherwise.

PARALLEL PROCESSING: threaded over segments `i`. Each `i` reads its own
column AND its neighbor `i+1` but writes only its own `g.HPG[:, i]` -- safe
under `Threads.@threads` since read-only access to a neighbor never races
with another thread's write (each thread owns a disjoint set of `i`
values, so no two threads ever write the same column) -- Pillar 1.
"""
function compute_pressure_gradient!(g, geom)
    kmx, imx = g.KMX, g.IMX
    g.HPG = zeros(Float64, kmx, imx)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            parallel_foreach((g.US[jb]-1):(g.DS[jb]-1)) do i
                g.HPG[kt, i] = g.DLXRHO[i] * (geom.BKT[i] + geom.BKT[i+1]) * 0.5 *
                               (geom.H1[kt, i+1] * g.P[kt, i+1] - geom.H1[kt, i] * g.P[kt, i])
                for k in (kt+1):min(g.KB[i], g.KB[i+1])
                    g.HPG[k, i] = g.DLXRHO[i] * geom.BHR2[k, i] *
                                  ((g.P[k-1, i+1] - g.P[k-1, i]) + (g.P[k, i+1] - g.P[k, i]))
                end
            end
        end
    end
    return g
end

"""
    solve_branch_free_surface!(g, geom, net, jw, kt, jb, dlt)

The per-branch body of the implicit free-surface elevation tridiagonal
solve (w2_4_win.f90:899-1022). Pulled out of `solve_free_surface!` so it
can be called under `Threads.@threads` for every branch in one
`branch_processing_tiers` tier -- see that function's docstring for why
branches in the same tier never touch each other's segments (each only
reads a different, already-resolved branch's `Z`/`EL` via
`net.upstream_branch`/`net.downstream_branch`, and only writes its own
`CUS[jb]:DS[jb]` range plus its own one-segment boundary pad).

NOT PORTED: the dam-flow/reciprocal-head-flow branch (HEAD_FLOW(JB) case at
w2_4_win.f90:1020) -- unreachable for Detroit (UHS all 0), see module
docstring.
"""
function solve_branch_free_surface!(g, geom, net, jw, kt, jb, dlt)
    g.BR_INACTIVE[jb] && return
    iu, id = g.CUS[jb], g.DS[jb]

    bhrho = zeros(Float64, g.IMX)
    d = zeros(Float64, g.IMX)
    f = zeros(Float64, g.IMX)

    for i in iu:(id-1)
        for k in kt:g.KBMIN[i]
            bhrho[i] += geom.BH2[k, i+1] / g.RHO[k, i+1] + geom.BH2[k, i] / g.RHO[k, i]
        end
        for k in kt:g.KB[i]
            d[i] += g.U[k, i] * geom.BHR2[k, i] - g.U[k, i-1] * geom.BHR2[k, i-1] - g.QSS[k, i] +
                    (g.UXBR[k, i] - g.UXBR[k, i-1]) * dlt
            f[i] += -g.SB[k, i] + g.ST[k, i] - g.ADMX[k, i] + g.DM[k, i] - g.HPG[k, i] + g.GRAV[k, i]
        end
    end

    d[iu] = 0.0
    for k in kt:g.KB[iu]
        d[iu] += g.U[k, iu] * geom.BHR2[k, iu] - g.QSS[k, iu] + g.UXBR[k, iu] * dlt
    end
    if g.DN_FLOW[jb]
        for k in kt:g.KB[id]
            d[id] += -g.U[k, id-1] * geom.BHR2[k, id-1] - g.QSS[k, id] + (g.UXBR[k, id] - g.UXBR[k, id-1]) * dlt
        end
        # External downstream outflow (w2_4_win.f90:914-921's D(ID) loop
        # includes "+QOUT(K,JB)" for every K -- reduced-physics here lumps
        # all withdrawal into g.QOT[jb] (a scalar total, not per-layer; see
        # IO/BoundaryReader.jl / Core/State.jl's QOT docstring), added once
        # rather than inside the K-loop -- correct sum given QOUT(K,JB)
        # would be 0 at every layer except the single reduced-physics
        # withdrawal layer in the real per-layer formula.
        d[id] += g.QOT[jb]
    end
    if g.UP_HEAD[jb]
        for k in kt:g.KBMIN[iu-1]
            bhrho[iu-1] += geom.BH2[k, iu] / g.RHO[k, iu] + geom.BH2[k, iu-1] / g.RHO[k, iu-1]
        end
        for k in kt:g.KB[iu]
            d[iu] -= g.U[k, iu-1] * geom.BHR2[k, iu-1]
            f[iu-1] -= g.SB[k, iu-1] - g.ST[k, iu-1] + g.HPG[k, iu-1] - g.GRAV[k, iu-1]
        end
    end
    if g.DN_HEAD[jb]
        for k in kt:g.KBMIN[id]
            bhrho[id] += geom.BH2[k, id+1] / g.RHO[k, id+1] + geom.BH2[k, id] / g.RHO[k, id]
        end
        d[id] = 0.0
        f[id] = 0.0
        for k in kt:g.KB[id]
            d[id] += g.U[k, id] * geom.BHR2[k, id] - g.U[k, id-1] * geom.BHR2[k, id-1] - g.QSS[k, id] +
                     (g.UXBR[k, id] - g.UXBR[k, id-1]) * dlt
            f[id] += -g.SB[k, id] + g.ST[k, id] - g.HPG[k, id] + g.GRAV[k, id]
        end
    end

    # --- External upstream inflow (w2_4_win.f90:944-956) -- traced insertion
    # point: this happens AFTER the DN_HEAD block above and BEFORE the
    # boundary surface elevations below, as a separate adjustment pass in
    # the real source, not part of the main D/F assembly loop. Only genuine
    # external UP_FLOW branches (INTERNAL_FLOW/DAM_INFLOW not ported, see
    # IO/BoundaryReader.jl module docstring) -- g.QIND[jb] is 0 for every
    # other branch (never loaded), so this is a no-op for them. ---
    g.UP_FLOW[jb] && (d[iu] -= g.QIND[jb])

    # --- Boundary surface elevations (w2_4_win.f90:958-988) ---
    if g.UH_INTERNAL[jb]
        jjb, jjw = net.upstream_branch[jb], net.upstream_waterbody[jb]
        geom.Z[iu-1] = ((-geom.EL[g.KTWB[jjw], g.UHS[jb]] + geom.Z[g.UHS[jb]] * geom.COSA[jjb]) +
                        geom.EL[kt, iu-1] + geom.SINA[jb] * geom.DLXR[iu-1]) / geom.COSA[jb]
        geom.ELWS[iu-1] = geom.EL[kt, iu-1] - geom.Z[iu-1] * geom.COSA[jb]
    end
    if g.DH_INTERNAL[jb]
        jjb, jjw = net.downstream_branch[jb], net.downstream_waterbody[jb]
        geom.Z[id+1] = ((-geom.EL[g.KTWB[jjw], g.DHS[jb]] + geom.Z[g.DHS[jb]] * geom.COSA[jjb]) +
                        geom.EL[kt, id+1]) / geom.COSA[jb]
        geom.ELWS[id+1] = geom.EL[kt, id+1] - geom.Z[id+1] * geom.COSA[jb]
    end

    # --- Implicit water surface elevation solution ---
    a = zeros(Float64, g.IMX); v = ones(Float64, g.IMX); c = zeros(Float64, g.IMX)
    for i in iu:id
        a[i] = -g.RHO[kt, i-1] * G_GRAVITY * geom.COSA[jb] * dlt^2 * bhrho[i-1] * 0.5 / geom.DLXR[i-1]
        c[i] = -g.RHO[kt, i+1] * G_GRAVITY * geom.COSA[jb] * dlt^2 * bhrho[i] * 0.5 / geom.DLXR[i]
        v[i] = g.RHO[kt, i] * G_GRAVITY * geom.COSA[jb] * dlt^2 * (bhrho[i] * 0.5 / geom.DLXR[i] + bhrho[i-1] * 0.5 / geom.DLXR[i-1]) +
               geom.DLX[i] * geom.BI[kt, i]
        d[i] = dlt * (d[i] + dlt * (f[i] - f[i-1])) + geom.DLX[i] * geom.BI[kt, i] * geom.Z[i]
    end
    g.UP_HEAD[jb] && (d[iu] -= a[iu] * geom.Z[iu-1])
    g.DN_HEAD[jb] && (d[id] -= c[id] * geom.Z[id+1])

    zseg = @view geom.Z[iu:id]
    thomas_solve!((@view a[iu:id]), (@view v[iu:id]), (@view c[iu:id]), (@view d[iu:id]), zseg)

    g.UP_FLOW[jb] && !g.HEAD_FLOW[jb] && (geom.Z[iu-1] = geom.Z[iu])
    g.DN_FLOW[jb] && (geom.Z[id+1] = geom.Z[id])

    geom.ELWS[iu-1:id+1] .= geom.EL[kt, iu-1:id+1] .- geom.Z[iu-1:id+1] .* geom.COSA[jb]
    return nothing
end

"""
    solve_free_surface!(g, geom, net, dlt)

w2_4_win.f90:899-1022 -- the implicit free-surface elevation tridiagonal
solve. Branch-sequenced via `Core/Grid.jl`'s `branch_processing_tiers` (a
REAL requirement, not a nicety -- see that function's docstring for the
traced proof this is a genuine sequential dependency across TIERS, not
parallel like the TKE/momentum TRIDIAG). PARALLEL PROCESSING (Pillar 1):
branches *within* one tier have no dependency on each other by
construction, so they run under `Threads.@threads` (`solve_branch_free_
surface!`); tiers themselves run in order, one full tier's `Threads.@threads`
loop always completing (an implicit barrier) before the next tier starts,
so a later tier's cross-branch boundary read always sees the earlier
tier's fully-updated `Z`/`ELWS`. Different waterbodies (`jw` loop) are
independent of each other too (disjoint segment ranges, no `net.upstream_
waterbody`/`downstream_waterbody` link crosses back into an earlier `jw`
in this port's supported topologies) -- also threaded.
"""
function solve_free_surface!(g, geom, net, dlt)
    parallel_foreach(1:g.NWB) do jw
        kt = g.KTWB[jw]
        for tier in branch_processing_tiers(g, net, jw)
            parallel_foreach(tier) do jb
                solve_branch_free_surface!(g, geom, net, jw, kt, jb, dlt)
            end
        end
    end
    return (g, geom)
end

"""
    recompute_top_layer_geometry!(g, geom)

w2_4_win.f90:1024-1091 (the non-TRAPEZOIDAL "Updated surface layer and
geometry" block, RECT only -- matches this port's existing TRAPEZOIDAL
scope restriction, see `Core/InitGeometry.jl`'s module docstring) --
recomputes `H1`/`BH1`/`BHR1`/`AVH1`/`BI`/`BKT`/`VOL` at the top active
layer `KT` from the just-solved `geom.Z` (`solve_free_surface!`'s output).

THIS WAS A REAL, MISSING PREREQUISITE, not a speculative addition --
found via testing `apply_temperature_sources!`/TDTR against real DET data
(2026-08-23): `H1`/`BH1`/`BHR1`/`AVH1` are populated ONCE at init
(`Core/InitGeometry.jl`'s `H1 .= H2` etc., itself just an "old=new" copy
for the very first step) and then NEVER updated again anywhere in this
codebase, even though `solve_free_surface!` changes `geom.Z` every
timestep. Concretely, this meant any `QSS`-driven volume change (from
QIN/QOT/QDTR) never reached the transport equations' dilution term --
`assemble_transport_rhs!` divides by a `geom.BH1` that was frozen at the
initial condition -- so a TDTR sign-check came out warming when it should
have shown cooling. This function is the fix; ported directly from the
exact formula `Core/InitGeometry.jl`'s own init-time `BH2`/`BKT`/`BI`/
`AVH2`/`AVHR`/`BHR2` computation already uses (same math, "1"-suffix
arrays, current `geom.Z` instead of the init-time `Z`) -- not guessed, the
init code IS the reference translation of this same real Fortran block.

REDUCED PHYSICS, first cut (user authorized starting this port 2026-08-23;
these specific simplifications are this port's own scoping call, flagged
here rather than confirmed line-by-line, same discipline as INTERP_INFLOW):
- No `KTI(I)` DO-WHILE crossing adjustment (w2_4_win.f90:1028-1044) -- this
  port assumes `Z` stays within the SAME discrete bathymetry sub-layer
  `g.KTI[i]` was assigned at init for the whole run. Reasonable for the
  small `ELWS` excursions (order 0.01-0.5m) seen in every validation run so
  far, wrong once a real multi-year forced run pushes `Z` across a whole
  sub-layer boundary -- port the crossing adjustment before trusting a
  long, strongly-forced run.
- No `CONSTRICTION`/`BCONSTRICTION` correction -- already a no-op
  elsewhere in this port (`Core/InitGeometry.jl`: "CONSTRICTION not
  ported... always false"), consistent here.
- No `KBI(I) < KB(I)` thin-bottom-layer correction to `BKT`/`AVHR` --
  `Core/State.jl`'s `KBI` field exists but this port doesn't yet use it;
  flagged, not silently applied.

MUST be called after `solve_free_surface!` (needs the new `geom.Z`) and
before anything that reads `H1`/`BH1`/`BHR1`/`AVH1` for the CURRENT step
(`update_velocities!`, and a future `temperature_transport!`/
`constituent_transport!` call once wired into this driver). The real
`update.F90` "old<-new" swap (`H2=H1` etc., done ONCE per timestep before
the new Z-solve) is `hydrodynamic_step!`'s job, not this function's --
this function only ever WRITES the "1" (new) arrays.
"""
function recompute_top_layer_geometry!(g, geom)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            iu, id = g.CUS[jb], g.DS[jb]
            for i in (iu-1):(id+1)
                kti = g.KTI[i]
                geom.H1[kt, i] = geom.H[kt, i] - geom.Z[i]
                geom.AVH1[kt, i] = (geom.H1[kt, i] + geom.H1[kt+1, i]) * 0.5

                bh1_kt = geom.B[kti, i] * (geom.EL[kt, i] - geom.EL[kti+1, i] - geom.Z[i] * geom.COSA[jb]) / geom.COSA[jb]
                kt == kti && (bh1_kt = geom.H1[kt, i] * geom.B[kt, i])
                for k in (kti+1):kt
                    bh1_kt += geom.BNEW[k, i] * geom.H[k, i]
                end
                geom.BH1[kt, i] = bh1_kt
                geom.BKT[i] = bh1_kt / geom.H1[kt, i]
                geom.BI[kt, i] = geom.B[kti, i]
                g.VOL[kt, i] = bh1_kt * geom.DLX[i]
            end
            for i in (iu-1):id
                geom.AVHR[kt, i] = geom.H1[kt, i] + (geom.H1[kt, i+1] - geom.H1[kt, i]) * geom.DLX[i] / (geom.DLX[i] + geom.DLX[i+1])
                geom.BHR1[kt, i] = geom.BH1[kt, i] + (geom.BH1[kt, i+1] - geom.BH1[kt, i]) * geom.DLX[i] / (geom.DLX[i] + geom.DLX[i+1])
            end
            geom.AVHR[kt, id+1] = geom.H1[kt, id+1]
            geom.BHR1[kt, id+1] = geom.BH1[kt, id+1]
        end
    end
    return g
end

"""
    distribute_tributary!(g, geom)

hydroinout.F90:1322-1331 -- distributes each branch's `QDTR[jb]` (loaded by
`IO/BoundaryReader.jl`) across every segment in that branch, weighted by
the segment's share of the branch's TOTAL top-layer surface area
(`AKBR = sum(BI(KT,I)*DLX(I))` over the branch), into `g.QSS[kt, i]` --
the real formula, ported faithfully (not simplified, unlike PLACE_QIN/
selective-withdrawal). Only branches with `g.DIST_TRIBS[jb]` true AND a
nonzero `g.QDTR[jb]` are touched.

MUST be called before `solve_free_surface!` in the same timestep (it
supplies `QSS`, which the free-surface solve reads), and after `fill!(g.QSS,
0.0)` -- matching `w2_4_win.f90:1494`'s real per-timestep `QSS = 0.0` reset
(`QSS` is an accumulator other real modules like lateral withdrawal also
add to; this port has no other writer yet, but the reset is still the
correct traced behavior, not a no-op left in "just in case").

Also stores the per-segment share into `g.QDT[i]` (real Fortran name
`QDT(I)`, hydroinout.F90:1329) -- NOT just a local intermediate here, since
`temperature.F90:416-427` reuses this exact per-segment value for the
distributed-tributary heat term (`TSS(KT,I) += TDTR(JB)*QDT(I)`,
`Hydrodynamics/Transport.jl`'s `apply_temperature_sources!`). Segments
outside any `DIST_TRIBS` branch keep whatever `QDT` held from the previous
step's `fill!` in the caller -- callers zero it the same way `QSS`/`TSS` are
zeroed (`hydrodynamic_step!`'s `fill!(g.QDT, 0.0)`).
"""
function distribute_tributary!(g, geom)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            (g.DIST_TRIBS[jb] && g.QDTR[jb] != 0.0) || continue
            iu, id = g.CUS[jb], g.DS[jb]
            akbr = sum(geom.BI[kt, i] * geom.DLX[i] for i in iu:id)
            for i in iu:id
                qdt_i = g.QDTR[jb] * geom.BI[kt, i] * geom.DLX[i] / akbr
                g.QDT[i] = qdt_i
                g.QSS[kt, i] += qdt_i
            end
        end
    end
    return g
end

"""
    apply_inflow_boundary!(g, geom, dlt)

Real `PLACE_QIN` port (`w2_4_win.f90:1209-1270`), ported 2026-08-24 --
motivated by comparing a `run_forced_simulation!` run against real DET
`outputs/tsr_1_seg33.csv` reference data (2026-08-23/24): the earlier
all-at-`KT` reduction produced meaningfully more volatile local
temperature swings than the real run (every branch's inflow, regardless of
temperature, dumped straight into the surface layer instead of sinking or
rising to its density-matched depth), plausibly the dominant cause of both
the extra volatility AND the steady `ELWS`/`T2` drift seen in that
comparison.

Two real branches, both ported (gated on `geom.PLACE_QIN[jw]`, required
explicit, Tier 1 -- same discipline as `THETA`/`UPWIND`/`ULTIMATE`):

- `PLACE_QIN[jw] == true`: the real density-driven plunge-point search.
  Computes the inflow's own density (`RHOIN`, via `Hydrodynamics/
  Density.jl`'s `density()`, reduced to `TDS=SS=0` -- no Tier-1 `CIN`
  constituent-loading IO to feed real values, matching every other
  `density()` call site in this port) and walks DOWN from `KT` until
  finding the first layer whose ambient density is `>= RHOIN` (the
  neutral-buoyancy layer). Then spreads this timestep's inflow VOLUME
  (`QIND[jb]*dlt`) across layers around that point, filling each to at
  most 50% of its own volume before moving on (first upward toward `KT`,
  then downward), matching the real `QINF`/`KTQIN`/`KBQIN` bookkeeping
  exactly -- not simplified further.
- `PLACE_QIN[jw] == false`: the real (not this port's OWN prior
  simplification) "off" behavior -- spreads inflow across the WHOLE
  active column proportional to each layer's volume (`BH1(K,IU)`), not
  concentrated at `KT`. Note this is itself already less reduced than
  this port's PRE-2026-08-24 behavior, which put 100% at `KT` regardless
  of `PLACE_QIN` -- that was this port's OWN extra simplification on top
  of the real "off" case, not a faithful rendering of it.

Both branches finish by setting the real per-layer boundary velocity
`U[K,IU-1] = QINF[K,jb]*QIND[jb]/BHR1[K,IU-1]` for every `K` in
`KT:KB[iu]` (previously only `K=KT` was ever written).

Call once per timestep, after `solve_free_surface!` (needs current
`geom.Z`-derived `VOL`/`BH1`) and `compute_density_field!` (needs current
`g.RHO` for the `PLACE_QIN=true` branch's plunge search) -- both already
run earlier in `hydrodynamic_step!`.
"""
function apply_inflow_boundary!(g, geom, dlt)
    # Lazy-sized default (false, the safe/conservative choice -- see
    # allocate_hydro_state!'s comment for why this isn't sized there
    # instead): apply_inflow_boundary! runs unconditionally inside every
    # hydrodynamic_step! call, so an empty geom.PLACE_QIN would BoundsError
    # every existing caller, same class of issue as DIST_TRIBS earlier this
    # session.
    isempty(geom.PLACE_QIN) && (geom.PLACE_QIN = fill(false, g.NWB))
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            (g.UP_FLOW[jb] && g.QIND[jb] != 0.0) || continue
            iu = g.CUS[jb]
            kb = g.KB[iu]
            fill!(view(g.QINF, :, jb), 0.0)

            if geom.PLACE_QIN[jw]
                rhoin = density(g.TIND[jb], 0.0, 0.0, true, false, false)
                k = kt
                while rhoin > g.RHO[k, iu] && k < kb
                    k += 1
                end
                g.KTQIN[jb] = k
                g.KBQIN[jb] = k

                vqin = g.QIND[jb] * dlt
                vqini = vqin
                qinfr = 1.0
                incr = -1
                while qinfr > 0.0
                    if k <= kb
                        v1 = g.VOL[k, iu]
                        if vqin > 0.5 * v1
                            g.QINF[k, jb] = 0.5 * v1 / vqini
                            qinfr -= g.QINF[k, jb]
                            vqin -= g.QINF[k, jb] * vqini
                            if k == kt
                                k = g.KBQIN[jb]
                                incr = 1
                            end
                        else
                            g.QINF[k, jb] = qinfr
                            qinfr = 0.0
                        end
                        incr < 0 && (g.KTQIN[jb] = k)
                        incr > 0 && (g.KBQIN[jb] = min(kb, k))
                        k += incr
                    else
                        g.QINF[kt, jb] += qinfr
                        qinfr = 0.0
                    end
                end
            else
                g.KTQIN[jb] = kt
                g.KBQIN[jb] = kb
                bhsum = sum(geom.BH1[k, iu] for k in kt:kb)
                for k in kt:kb
                    g.QINF[k, jb] = geom.BH1[k, iu] / bhsum
                end
            end

            for k in kt:kb
                g.U[k, iu-1] = g.QINF[k, jb] * g.QIND[jb] / geom.BHR1[k, iu-1]
            end
        end
    end
    return g
end

"""
    compute_horizontal_advection_of_momentum!(g, geom)

w2_4_win.f90:845-854 -- explicit horizontal advection of momentum, `ADMX`.
Ported for real: all its inputs (`U`, `BH2`, `DLXR`) are already real in
this port, unlike `ADMZ`/`DM` (see below) -- no new Tier-1 IO needed.
Upwind-selected via `UDR`/`UDL` (`DSIGN`-based direction switches in the
real Fortran) -- translated as a `>= 0.0` ternary rather than Julia's
`sign`, since `DSIGN(1.0,0.0)` treats exactly-zero as positive and `sign(0)
== 0` would silently zero out both branches of the blend at that one point.

Found relevant while investigating why `U` stayed 0 at a real inflow-forced
interior segment (2026-08-23, see CLAUDE.md's `IO/OutputWriter.jl` `KTWB`
bug entry) -- `ADMX`/`ADMZ`/`DM` were all flagged "NOT YET COMPUTED" from
the original `FreeSurface.jl` work; this ports the one of the three with no
missing prerequisite. `ADMZ` (needs `W`, vertical velocity -- continuity-
derived, not ported) and `DM` (needs `AX(JW)`, horizontal eddy viscosity --
Tier 1, not read by `InputReader.jl` yet) remain zero, flagged separately
in `Core/State.jl`'s existing field comments.

MUST be called after `apply_inflow_boundary!` (needs the boundary `U` set,
since `ADMX[kt,iu]` reads `U[kt,iu-1]`) and before `update_velocities!`
(consumes `ADMX` as an explicit forcing term) in the same timestep.

PARALLEL PROCESSING: threaded over segments `i`, matching `update_
velocities!`'s own pattern -- each `i` only writes its own `ADMX[:,i]`,
reading neighbor columns `i-1`/`i+1` (never writing them).
"""
function compute_horizontal_advection_of_momentum!(g, geom)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            iu, id = g.CUS[jb], g.DS[jb]
            parallel_foreach(iu:(id-1)) do i
                for k in kt:g.KBMIN[i]
                    udr = (g.U[k, i] + g.U[k, i+1]) * 0.5 >= 0.0 ? 1.0 : 0.0
                    udl = (g.U[k, i] + g.U[k, i-1]) * 0.5 >= 0.0 ? 1.0 : 0.0
                    g.ADMX[k, i] = (geom.BH2[k, i+1] * (g.U[k, i+1] + g.U[k, i]) * 0.5 * (udr * g.U[k, i] + (1.0 - udr) * g.U[k, i+1]) -
                                    geom.BH2[k, i] * (g.U[k, i] + g.U[k, i-1]) * 0.5 * (udl * g.U[k, i-1] + (1.0 - udl) * g.U[k, i])) / geom.DLXR[i]
                end
            end
        end
    end
    return g
end

"""
    update_velocities!(g, geom, dlt)

w2_4_win.f90:1301-1307 -- explicit horizontal velocity update. `ADMX` is now
real (`compute_horizontal_advection_of_momentum!`, 2026-08-23); `SB`/`ST`/
`ADMZ`/`DM` are still zero (see module docstring). For a zero-flow start
with uniform density and flat slope this still reduces to `U_new =
U_old*BHR2/BHR1` -- stays at 0, as expected -- since `ADMX` itself is 0
wherever `U` is uniformly 0.

PARALLEL PROCESSING: threaded over segments `i`, each writing only its own
column `g.U[:, i]` and reading only its own column's data -- embarrassingly
parallel (Pillar 1).
"""
function update_velocities!(g, geom, dlt)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            # IU here is CUS(JB) (the first PROPERLY-ACTIVE segment, not US(JB))
            # matching w2_4_win.f90:948 "IU = CUS(JB)", set earlier in the same
            # per-branch iteration and never reset before this section runs.
            # Using US(JB) instead would walk into the dry upstream segments
            # (e.g. Detroit branch 1's segments 2-5) where BHR1 can be 0 --
            # division by zero, the actual bug this comment replaced.
            iu, id = g.CUS[jb], g.DS[jb]
            parallel_foreach(iu:(id-1)) do i
                for k in kt:g.KBMIN[i]
                    g.U[k, i] = (geom.BHR2[k, i] * g.U[k, i]) / geom.BHR1[k, i] +
                                (dlt * (-g.SB[k, i] + g.ST[k, i] - g.ADMZ[k, i] + (k > kt ? g.ADMZ[k-1, i] : 0.0) -
                                        g.ADMX[k, i] + g.DM[k, i] - g.HPG[k, i] + g.GRAV[k, i] +
                                        g.UXBR[k, i] / geom.H2[k, i])) / geom.BHR1[k, i]
                end
            end
        end
    end
    return g
end

"""
    hydrodynamic_step!(g, geom, net, dlt)

Orchestrates one reduced-physics hydrodynamic timestep: old<-new geometry
swap -> reset QSS/TSS/QDT -> distribute tributary inflow -> density ->
pressure -> gravity -> pressure gradient -> free-surface solve -> top-layer
geometry recompute -> inflow boundary -> velocity update. See module
docstring for exactly what's real vs. stubbed. `fill!(g.QSS, 0.0)` matches
`w2_4_win.f90:1494`'s real per-timestep reset (see `distribute_tributary!`'s
docstring for why this isn't a no-op). `g.TSS` is reset here too --
`update.F90:44-45` resets `QSS(K,I)` and `TSS(K,I)` together in the same
real per-timestep pass, and `Hydrodynamics/Transport.jl`'s
`apply_temperature_sources!` (called by a future driver, not yet wired into
this function) depends on `TSS` starting each step at zero.

The `geom.H2 .= geom.H1` etc. swap at the top matches `update.F90:42-53`'s
real "old<-new" pass (`H2(K,I)=H1(K,I)`, `BH2(K,I)=BH1(K,I)`, etc.) -- MUST
run before `recompute_top_layer_geometry!` overwrites the "1" (new) arrays
below, so `update_velocities!`'s `BHR2/BHR1` ratio and a future transport
call's `cold`/`BH2` terms see last step's new values as this step's old
baseline, not this step's own new values. `recompute_top_layer_geometry!`
runs right after `solve_free_surface!` since it needs the freshly solved
`geom.Z` -- see that function's docstring for what was found missing (real,
not speculative: `H1`/`BH1`/`BHR1`/`AVH1` were frozen at their init values
forever until this was added, 2026-08-23).

`geom.BHRATIO .= geom.BH2 ./ geom.BH1` runs BEFORE the swap (same real
`update.F90:51` line, computed immediately before that line's own
`BH2(K,I)=BH1(K,I)`) -- real Fortran's comment there is literally "USED FOR
TKE COMPUTATION" (`Hydrodynamics/Turbulence.jl`'s `calculate_tke!` reads
`geom.BHRATIO[kt,i]`). Computed array-wide via broadcast rather than only
at `[kt,i]` (matching real Fortran's own `K=KT,KB(I)` range, in case a
future closure needs other `K`) -- may produce `NaN`/`Inf` at dry/inactive
positions where `BH1==0`, harmless since nothing reads `BHRATIO` there.

`calculate_tke!` runs LAST, after `update_velocities!` -- needs this step's
just-solved `U` (shear production) and `compute_density_field!`'s `RHO`
(buoyancy damping, already current from earlier in this same call). Writes
`AZ`/`DZ` for `Hydrodynamics/Transport.jl` to read on the NEXT call (e.g.
`apply_temperature_sources!`/`temperature_transport!`, called by
`Simulation.jl`'s `run_forced_simulation!` right after this function
returns in the same timestep).
"""
function hydrodynamic_step!(g, geom, net, dlt)
    geom.BHRATIO .= geom.BH2 ./ geom.BH1
    geom.H2 .= geom.H1
    geom.BH2 .= geom.BH1
    geom.BHR2 .= geom.BHR1
    geom.AVH2 .= geom.AVH1
    fill!(g.QSS, 0.0)
    fill!(g.TSS, 0.0)
    fill!(g.QDT, 0.0)
    distribute_tributary!(g, geom)
    compute_density_field!(g, geom)
    compute_pressure_field!(g, geom)
    compute_gravity_term!(g, geom)
    compute_pressure_gradient!(g, geom)
    solve_free_surface!(g, geom, net, dlt)
    recompute_top_layer_geometry!(g, geom)
    apply_inflow_boundary!(g, geom, dlt)
    compute_horizontal_advection_of_momentum!(g, geom)
    update_velocities!(g, geom, dlt)
    calculate_tke!(g, geom, dlt)
    return (g, geom)
end
