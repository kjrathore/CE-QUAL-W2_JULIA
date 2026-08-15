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
reads). Call after `IO/InputReader.allocate_geometry!`, before
`hydrodynamic_step!`.
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
    update_velocities!(g, geom, dlt)

w2_4_win.f90:1301-1307 -- explicit horizontal velocity update. With SB/ST/
ADMX/ADMZ/DM all zero (see module docstring) and HPG/GRAV zero for Detroit's
uniform density and flat slope, this reduces to `U_new = U_old*BHR2/BHR1`
-- stays at 0 given a zero-flow start, as expected.

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

Orchestrates one reduced-physics hydrodynamic timestep: density -> pressure
-> gravity -> pressure gradient -> free-surface solve -> velocity update.
See module docstring for exactly what's real vs. stubbed.
"""
function hydrodynamic_step!(g, geom, net, dlt)
    compute_density_field!(g, geom)
    compute_pressure_field!(g, geom)
    compute_gravity_term!(g, geom)
    compute_pressure_gradient!(g, geom)
    solve_free_surface!(g, geom, net, dlt)
    update_velocities!(g, geom, dlt)
    return (g, geom)
end
