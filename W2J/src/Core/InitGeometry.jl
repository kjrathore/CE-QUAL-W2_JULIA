# ==============================================================================
# Core/InitGeometry.jl
#
# Ports INITGEOM (w2_source_5-21-2026/init-geom.F90, 641 lines) plus the small
# upstream dependency it assumes is already computed, ALPHA/SINA/COSA/SINAC
# and the branch boundary-condition flags (both from init.F90, since
# init-geom.F90 only USEs GEOMC/GLOBAL and reads them as already populated).
#
# WHY THIS FILE EXISTS: this is the single most foundational missing piece
# for any hydrodynamics work -- EL(K,I) (layer elevation) and KB(I) (bottom
# active layer index) are referenced by nearly every downstream routine
# (withdrawal/structure/tributary placement, the free-surface solve, volume
# and area calcs for kinetics). See CLAUDE.md "Open questions" -- the
# `H(K,JW)` shape mismatch flagged there is resolved here (see
# `compute_layer_elevations!`).
#
# SCOPE -- DELIBERATELY NOT A FULL PORT. `init-geom.F90` is 641 lines
# covering geometry finalization, boundary-condition bookkeeping, AND
# several concerns this pass explicitly defers (flagged at each site, not
# silently skipped):
#   - RESTART_IN is assumed false throughout (no restart-file reading exists
#     in W2J yet) -- this only skips a branch that recomputes KTI/water-
#     surface state from a restart file instead of from ELWS.
#   - TRAPEZOIDAL(JW) is now read for real (GRIDC, via InputReader.jl's
#     initial-condition block) and computed by `compute_trapezoidal_flag!` --
#     NOT hardcoded to RECT. Only the non-trapezoidal area/volume branch
#     (init-geom.F90:489-517) is ported, though; `compute_areas_volumes!`
#     ERRORS loudly if a waterbody's GRIDC=="TRAP" rather than silently
#     computing wrong (RECT-formula) areas for it.
#   - The w2_constriction.csv optional-file read (init-geom.F90:309-327) is
#     not ported -- CONSTRICTION stays all-false, BCONSTRICTION all-zero,
#     matching any control file (like Detroit's) that doesn't use it.
#   - The cross-branch boundary COLUMN WIDTH interpolation for internal
#     branch junctions (init-geom.F90:330-450, the ELR/ELL/B11 elevation-
#     matching logic) is APPROXIMATED, not faithfully ported: this pass
#     copies the adjacent internal branch's column width at the same layer
#     index instead of the elevation-interpolated value. Flagged loudly in
#     `extend_boundary_widths!` -- this is a real accuracy gap for branched
#     waterbodies (Detroit has 4 branches within 1 waterbody, so it IS
#     exercised), not a cosmetic one. Fix by porting the real interpolation
#     before trusting boundary-adjacent cell widths in a branched case.
#   - JBTR/JBWD (which branch a tributary/withdrawal belongs to,
#     init-geom.F90:211-220) are NOT computed -- they need ITR/IWD (segment
#     locations), which InputReader.jl doesn't read yet (Tier 1: the
#     structures/withdrawals block of w2_con.csv). Will be added when that
#     lands.
#   - NTAC/NTACMX/NTACMN (active cell counts) and KBR (snapshot-output
#     bookkeeping) are output/diagnostic-only in the original and are
#     skipped entirely.
#   - The DAM_FLOW / reciprocal head-flow branch-linkage detection in
#     init.F90:212-238 (triggered only when UHS(JB) is nonzero) is not
#     ported -- Detroit has UHS=0 for every branch (external upstream inflow
#     everywhere), so this never exercises against the one real dataset
#     available. Flagged in `compute_boundary_flags!` rather than guessed at.
#
# VALIDATION: see test/runtests.jl, "Core/InitGeometry: init_geometry!" --
# checked against real Detroit data (EL(KMX,:) == ELBOT, KT/KB sanity,
# positive volumes). The non-zero-slope branch-walking path in
# `compute_layer_elevations!` is ported but UNVALIDATED -- Detroit's single
# waterbody has SLOPE=0 on every branch (confirmed, see CLAUDE.md validated
# values), so the ZERO_SLOPE branch is the only one exercised by the
# available test data. Do not trust the sloped-reservoir path without a
# second test case.
# ==============================================================================

"""
    compute_angles!(geom)

init.F90:415-418 -- converts branch slope to angle/trig quantities used
throughout geometry and transport. Must run after SLOPE/SLOPEC are read
(InputReader.jl Phase A) and before `compute_layer_elevations!`.
"""
function compute_angles!(geom)
    geom.ALPHA = atan.(geom.SLOPE)
    geom.SINA = sin.(geom.ALPHA)
    geom.SINAC = sin.(atan.(geom.SLOPEC))
    geom.COSA = cos.(geom.ALPHA)
    return geom
end

"""
    compute_trapezoidal_flag!(geom)

init.F90:153 -- `TRAPEZOIDAL(JW) = (GRIDC(JW) == "TRAP")`. Depends only on
GRIDC (read by InputReader.jl's initial-condition block), unlike
`FRESH_WATER`/`SALT_WATER` below which also need Tier-1 fields not read yet
-- safe to compute unconditionally as part of `init_geometry!`.
"""
function compute_trapezoidal_flag!(geom)
    geom.TRAPEZOIDAL = [strip(uppercase(gc)) == "TRAP" for gc in geom.GRIDC]
    return geom
end

"""
    compute_water_type_flags!(g, geom; constituents::Bool, tds_active::Bool)

init.F90:172-173 -- `FRESH_WATER(JW)`/`SALT_WATER(JW)` = `CONSTITUENTS .AND.
WTYPEC(JW)==<match> .AND. CAC(NTDS)=="ON"`.

`constituents`/`tds_active` are REQUIRED explicit keyword arguments, not
defaulted to `true` and not silently inferred. The real Fortran dependency
is `CONSTITUENTS = (CCC=="ON")` (the water-quality-computations master
toggle, in the "CST COMP" control-file section) and `CAC(NTDS)=="ON"` (is
the TDS constituent itself active, in the CST block) -- neither is read
into `W2Global`/`W2Names` yet (Tier 1; both values, however, ARE already
visible in `w2_config.yaml` via `tools/xlsx_to_yaml`, under "CST COMP" and
"CST - ..." -> "TDS" -> "active", so a caller with that YAML loaded can
supply them today without guessing). Defaulting these to `true` would
silently diverge from the real model for any control file with water
quality off or TDS untracked -- not this project's style; NOT called from
`init_geometry!`'s automatic pipeline for that reason. This is separate
from geometry itself (`EL`/`KB`/`VOL`, which don't need it) -- it exists
for `Hydrodynamics/Density.jl`, which does.
"""
function compute_water_type_flags!(g, geom; constituents::Bool, tds_active::Bool)
    active = constituents && tds_active
    geom.FRESH_WATER = [active && strip(uppercase(wt)) == "FRESH" for wt in geom.WTYPEC]
    geom.SALT_WATER = [active && strip(uppercase(wt)) == "SALT" for wt in geom.WTYPEC]
    return (g, geom)
end

"""
    compute_boundary_flags!(g, geom)

init.F90:205-249 -- derives per-branch boundary-condition flags from
UHS/DHS/UQB/DQB, and the per-waterbody ZERO_SLOPE flag from SLOPE. Must run
after the branch-grid and waterbody-location blocks are read.

NOT PORTED: the DAM_FLOW / reciprocal head-flow detection (init.F90:212-238,
only entered when UP_HEAD(JB), i.e. UHS(JB) != 0) -- see module docstring.
UP_FLOW/UP_HEAD/UH_INTERNAL below are the *unconditional* values from
init.F90:208-211 only; they are not corrected for dam-flow/head-flow
reciprocal connections the way the real subroutine does when UHS(JB) != 0.
"""
function compute_boundary_flags!(g, geom)
    nbr = g.NBR
    g.UP_FLOW = g.UHS .== 0
    g.DN_FLOW = g.DHS .== 0
    g.UP_HEAD = g.UHS .!= 0
    g.UH_INTERNAL = g.UHS .> 0
    g.DH_INTERNAL = g.DHS .> 0
    g.DN_HEAD = g.DHS .!= 0
    g.UH_EXTERNAL = g.UHS .== -1
    g.DH_EXTERNAL = g.DHS .== -1
    g.UQ_EXTERNAL = g.UHS .== 0
    g.DQ_EXTERNAL = g.DHS .== 0
    g.DQ_INTERNAL = g.DQB .> 0
    g.UQ_INTERNAL = g.UQB .> 0   # NOT corrected for .AND. .NOT. DAM_INFLOW(JB) -- see docstring
    g.HEAD_FLOW = fill(false, nbr)       # only ever set true by the not-yet-ported dam-flow branch
    g.INTERNAL_FLOW = fill(false, nbr)   # ditto
    g.DAM_INFLOW = fill(false, nbr)      # ditto
    g.DAM_OUTFLOW = fill(false, nbr)     # ditto

    g.ZERO_SLOPE = [all(geom.SLOPE[jb] == 0.0 for jb in g.BS[jw]:g.BE[jw]) for jw in 1:g.NWB]
    return g
end

"""
    compute_layer_elevations!(g, geom)

init-geom.F90:17-126 -- computes EL(K,I), the elevation of the TOP of layer K
at segment I, for every waterbody. Two cases:

  - ZERO_SLOPE(JW): each segment's layers stack directly on ELBOT(JW),
    independent of branch connectivity (init-geom.F90:19-25). Simple and
    validated against Detroit.
  - Non-zero slope: walks the branch connectivity graph outward from the
    downstream-most branch (JBDN(JW)), propagating elevations across branch
    junctions (init-geom.F90:26-124). Ported faithfully but UNVALIDATED --
    see module docstring.

Resolves the `H(K,JW)` shape question flagged in CLAUDE.md: `geom.H` is
KMX x IMX in this port (matching how BathymetryReader.jl already populates
it, `bathy.H` per waterbody broadcast across that waterbody's segments), NOT
KMX x NWB like the raw Fortran declaration -- `H(K,JW)` in the formulas below
is read from `geom.H[K, i]` for the same reason (any segment `i` in
waterbody `JW` has identical `H[:, i]`, since BathymetryReader fills it that
way). This sidesteps the mismatch by never needing the per-waterbody-only
shape.
"""
function compute_layer_elevations!(g, geom)
    kmx = g.KMX
    npoint = fill(0, g.NBR)

    for jw in 1:g.NWB
        if g.ZERO_SLOPE[jw]
            for i in (g.US[g.BS[jw]]-1):(g.DS[g.BE[jw]]+1)
                geom.EL[kmx, i] = geom.ELBOT[jw]
                for k in (kmx-1):-1:1
                    geom.EL[k, i] = geom.EL[k+1, i] + geom.H[k, i]
                end
            end
            continue
        end

        # --- Non-zero-slope branch-walking case (init-geom.F90:26-124) ---
        jb = g.JBDN[jw]
        geom.EL[kmx, g.DS[jb]+1] = geom.ELBOT[jw]
        npoint .= 0
        npoint[jb] = 1
        nnbp = 1
        ncbp = 0
        ninternal = 0
        nup = 0
        jjb = jb

        while nnbp <= (g.BE[jw] - g.BS[jw] + 1)
            ncbp += 1
            if ninternal == 0
                if nup == 0
                    for i in g.DS[jb]:-1:g.US[jb]
                        if i != g.DS[jb]
                            geom.EL[kmx, i] = geom.EL[kmx, i+1] + geom.SINA[jb] * (geom.DLX[i] + geom.DLX[i+1]) * 0.5
                        else
                            geom.EL[kmx, i] = geom.EL[kmx, i+1]
                        end
                        for k in (kmx-1):-1:1
                            geom.EL[k, i] = geom.EL[k+1, i] + geom.H[k, i] * geom.COSA[jb]
                        end
                    end
                else
                    for i in g.US[jb]:g.DS[jb]
                        if i != g.US[jb]
                            geom.EL[kmx, i] = geom.EL[kmx, i-1] - geom.SINA[jb] * (geom.DLX[i] + geom.DLX[i-1]) * 0.5
                        else
                            geom.EL[kmx, i] = geom.EL[kmx, i-1]
                        end
                        for k in (kmx-1):-1:1
                            geom.EL[k, i] = geom.EL[k+1, i] + geom.H[k, i] * geom.COSA[jb]
                        end
                    end
                    nup = 0
                end
                for k in kmx:-1:1
                    geom.EL[k, g.US[jb]-1] = g.UP_HEAD[jb] ?
                        geom.EL[k, g.US[jb]] + geom.SINA[jb] * geom.DLX[g.US[jb]] :
                        geom.EL[k, g.US[jb]]
                    geom.EL[k, g.DS[jb]+1] = g.DN_HEAD[jb] ?
                        geom.EL[k, g.DS[jb]] - geom.SINA[jb] * geom.DLX[g.DS[jb]] :
                        geom.EL[k, g.DS[jb]]
                end
            else
                for k in (kmx-1):-1:1
                    geom.EL[k, g.UHS[jjb]] = geom.EL[k+1, g.UHS[jjb]] + geom.H[k, g.UHS[jjb]] * geom.COSA[jb]
                end
                for i in (g.UHS[jjb]+1):g.DS[jb]
                    geom.EL[kmx, i] = geom.EL[kmx, i-1] - geom.SINA[jb] * (geom.DLX[i] + geom.DLX[i-1]) * 0.5
                    for k in (kmx-1):-1:1
                        geom.EL[k, i] = geom.EL[k+1, i] + geom.H[k, i] * geom.COSA[jb]
                    end
                end
                for i in (g.UHS[jjb]-1):-1:g.US[jb]
                    geom.EL[kmx, i] = geom.EL[kmx, i+1] + geom.SINA[jb] * (geom.DLX[i] + geom.DLX[i+1]) * 0.5
                    for k in (kmx-1):-1:1
                        geom.EL[k, i] = geom.EL[k+1, i] + geom.H[k, i] * geom.COSA[jb]
                    end
                end
                ninternal = 0
            end

            nnbp == (g.BE[jw] - g.BS[jw] + 1) && break

            # Find next branch connected to the furthest-downstream branch so far.
            found = false
            for jb2 in g.BS[jw]:g.BE[jw]
                npoint[jb2] == 1 && continue
                for jjb2 in g.BS[jw]:g.BE[jw]
                    npoint[jjb2] != 1 && continue
                    if g.DHS[jb2] >= g.US[jjb2] && g.DHS[jb2] <= g.DS[jjb2]
                        npoint[jb2] = 1
                        geom.EL[kmx, g.DS[jb2]+1] = geom.EL[kmx, g.DHS[jb2]] + geom.SINA[jb2] * (geom.DLX[g.DS[jb2]] + geom.DLX[g.DHS[jb2]]) * 0.5
                        nnbp += 1
                        jb, found = jb2, true
                        break
                    elseif g.UHS[jjb2] == g.DS[jb2]
                        npoint[jb2] = 1
                        geom.EL[kmx, g.DS[jb2]+1] = geom.EL[kmx, g.US[jjb2]] + (geom.SINA[jjb2] * geom.DLX[g.US[jjb2]] + geom.SINA[jb2] * geom.DLX[g.DS[jb2]]) * 0.5
                        nnbp += 1
                        jb, found = jb2, true
                        break
                    elseif g.UHS[jjb2] >= g.US[jb2] && g.UHS[jjb2] <= g.DS[jb2]
                        npoint[jb2] = 1
                        geom.EL[kmx, g.UHS[jjb2]] = geom.EL[kmx, g.US[jjb2]] + geom.SINA[jjb2] * geom.DLX[g.US[jjb2]] * 0.5
                        nnbp += 1
                        ninternal = 1
                        jb, jjb, found = jb2, jjb2, true
                        break
                    elseif g.UHS[jb2] >= g.US[jjb2] && g.UHS[jb2] <= g.DS[jjb2]
                        npoint[jb2] = 1
                        geom.EL[kmx, g.US[jb2]-1] = geom.EL[kmx, g.UHS[jb2]] - geom.SINA[jb2] * geom.DLX[g.US[jb2]] * 0.5
                        nnbp += 1
                        nup = 1
                        jb, found = jb2, true
                        break
                    end
                end
                found && break
            end
        end
    end
    return g
end

"""
    allocate_init_geometry!(g, geom)

Sizes every array `init_geometry!` writes into that `IO/InputReader.jl`'s
`allocate_geometry!` doesn't already handle (that one only covers DLX/ELWS/B
for the bathymetry-read stage). Call after `allocate_geometry!` and
`BathymetryReader.read_bathymetry!`, before `init_geometry!`.
"""
function allocate_init_geometry!(g, geom)
    imx, kmx, nbr, nwb = g.IMX, g.KMX, g.NBR, g.NWB

    geom.EL = zeros(Float64, kmx, imx)
    geom.DLXR = zeros(Float64, imx)
    geom.Z = zeros(Float64, imx)
    # H2(:,I) = H(:,JW) baseline copy (input.F90:2265,2275) -- Fortran does this
    # per-segment during the bathymetry read itself, immediately after B(K,I) is
    # read. Here geom.H (already populated by BathymetryReader.read_bathymetry!,
    # which runs before this function) is copied wholesale instead, since
    # geom.H2 isn't allocated until this point in the port's call ordering.
    # Without this, every layer except KT (set later by
    # compute_water_surface_layer!) stays at zero forever, causing 0.0/0.0 NaNs
    # downstream (e.g. Hydrodynamics/FreeSurface.jl's UXBR/H2).
    geom.H1 = zeros(Float64, kmx, imx); geom.H2 = copy(geom.H)
    geom.BH1 = zeros(Float64, kmx, imx); geom.BH2 = zeros(Float64, kmx, imx)
    geom.BHR1 = zeros(Float64, kmx, imx); geom.BHR2 = zeros(Float64, kmx, imx)
    geom.AVHR = zeros(Float64, kmx, imx); geom.BHRATIO = ones(Float64, kmx, imx)
    geom.BI = zeros(Float64, kmx, imx); geom.BB = zeros(Float64, kmx, imx)
    geom.BH = zeros(Float64, kmx, imx); geom.BHR = zeros(Float64, kmx, imx)
    geom.BR = zeros(Float64, kmx, imx)
    geom.AVH1 = zeros(Float64, kmx, imx); geom.AVH2 = zeros(Float64, kmx, imx)
    geom.BNEW = zeros(Float64, kmx, imx)
    geom.DEPTHB = zeros(Float64, kmx, imx); geom.DEPTHM = zeros(Float64, kmx, imx)
    geom.FETCHU = zeros(Float64, imx, nbr); geom.FETCHD = zeros(Float64, imx, nbr)
    geom.HSEG = zeros(Float64, kmx, imx)
    geom.BKT = zeros(Float64, imx)
    geom.CONSTRICTION = fill(false, kmx, imx)
    geom.BCONSTRICTION = zeros(Float64, imx)
    geom.ONE_LAYER = fill(false, imx)
    geom.KBI = zeros(Int, imx)

    geom.JBUH = zeros(Int, nbr); geom.JBDH = zeros(Int, nbr)
    geom.JWUH = zeros(Int, nbr); geom.JWDH = zeros(Int, nbr)

    g.KB = zeros(Int, imx); g.KTI = fill(2, imx); g.KBMIN = zeros(Int, imx)
    g.KTWB = zeros(Int, nwb); g.KBMAX = zeros(Int, nwb)
    g.ELKT = zeros(Float64, nwb)
    g.CDHS = zeros(Int, nbr)
    g.CUS = zeros(Int, nbr)
    g.BR_INACTIVE = fill(false, nbr)
    g.VOL = zeros(Float64, kmx, imx)

    return (g, geom)
end

"""
    compute_water_surface_layer!(g, geom)

init-geom.F90:128-197 -- HMIN/HMAX/HMAX2, the top active layer KTI(I)/KTWB(JW)
(found by walking down from layer 2 until EL drops below ELWS), Z(I) (the
partial-layer thickness at the surface), the H2 correction for a water
surface spanning multiple layers, and ELKT(JW).

NOT PORTED: the `.NOT. RESTART_IN` guard is assumed always true (W2J has no
restart-file reader yet) -- see module docstring.
"""
function compute_water_surface_layer!(g, geom)
    kmx = g.KMX
    hmin, hmax = Inf, -Inf
    for jw in 1:g.NWB, k in (kmx-1):-1:1
        hmin = min(hmin, geom.H[k, g.US[g.BS[jw]]])
        hmax = max(hmax, geom.H[k, g.US[g.BS[jw]]])
    end
    g.HMIN, g.HMAX = hmin, hmax
    g.HMAX2 = hmax^2

    for jw in 1:g.NWB
        for jb in g.BS[jw]:g.BE[jw]
            for i in (g.US[jb]-1):(g.DS[jb]+1)
                kti = 2
                while geom.EL[kti, i] > geom.ELWS[i]
                    kti += 1
                    if kti == kmx
                        kti = kmx - 1
                        geom.ELWS[i] = geom.EL[kti, i]
                        break
                    end
                end
                geom.Z[i] = (geom.EL[kti, i] - geom.ELWS[i]) / geom.COSA[jb]
                ktmax = max(2, kti)
                g.KTWB[jw] = max(ktmax, g.KTWB[jw])
                geom.EL[kti, i] != geom.ELWS[i] && (kti -= 1)
                g.KTI[i] = max(kti, 2)
            end
        end

        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            for i in (g.US[jb]-1):(g.DS[jb]+1)
                geom.H2[kt, i] = geom.H[kt, i] - geom.Z[i]
                k = g.KTI[i] + 1
                while kt > k
                    geom.Z[i] -= geom.H[k, i]
                    geom.H2[kt, i] = geom.H[kt, i] - geom.Z[i]
                    k += 1
                end
            end
        end
        g.ELKT[jw] = geom.EL[g.KTWB[jw], g.DS[g.BE[jw]]] - geom.Z[g.DS[g.BE[jw]]] * geom.COSA[g.BE[jw]]
    end
    return (g, geom)
end

"""
    compute_bottom_layers!(g, geom, net)

init-geom.F90:164-283 -- KB(I) (first layer from the top with zero width,
minus one), boundary-cell KB copies, and the cross-branch bottom-layer
matching for internal head boundaries (UH_INTERNAL/DH_INTERNAL). Uses
`net::Core/Grid.BranchNetwork` (built by `build_branch_network`) for the
branch-linkage lookup instead of searching US/DS inline -- see Core/Grid.jl's
module docstring for why this was pulled out into its own file.

NOT PORTED: JBTR/JBWD (tributary/withdrawal branch mapping, init-geom.F90:
211-220) -- needs ITR/IWD, not read by InputReader.jl yet (Tier 1). See
module docstring.
"""
function compute_bottom_layers!(g, geom, net)
    for jw in 1:g.NWB
        for jb in g.BS[jw]:g.BE[jw]
            for i in (g.US[jb]-1):(g.DS[jb]+1)
                k = 2
                while geom.B[k, i] > 0.0
                    g.KB[i] = k
                    k += 1
                    k > g.KMX && break
                end
                g.KBMAX[jw] = max(g.KBMAX[jw], g.KB[i])
            end
            g.KB[g.US[jb]-1] = g.KB[g.US[jb]]
            g.KB[g.DS[jb]+1] = g.KB[g.DS[jb]]
        end
    end

    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            iu, id = g.US[jb], g.DS[jb]

            g.UH_EXTERNAL[jb] && (g.KB[iu-1] = g.KB[iu])
            g.DH_EXTERNAL[jb] && (g.KB[id+1] = g.KB[id])

            if g.UH_INTERNAL[jb]
                jbuh = net.upstream_branch[jb] == 0 ? nothing : net.upstream_branch[jb]
                jbuh !== nothing && (geom.JBUH[jb] = jbuh)
                if jbuh !== nothing && jbuh >= g.BS[jw] && jbuh <= g.BE[jw]
                    g.KB[iu-1] = min(g.KB[g.UHS[jb]], g.KB[iu])
                elseif jbuh !== nothing
                    if geom.EL[g.KB[iu], iu] >= geom.EL[g.KB[g.UHS[jb]], g.UHS[jb]]
                        g.KB[iu-1] = g.KB[iu]
                    else
                        for k in kt:g.KB[iu]
                            if geom.EL[g.KB[g.UHS[jb]], g.UHS[jb]] >= geom.EL[k, iu]
                                g.KB[iu-1] = k
                                break
                            end
                        end
                    end
                end
            end
            if g.DH_INTERNAL[jb]
                jbdh = net.downstream_branch[jb] == 0 ? nothing : net.downstream_branch[jb]
                jbdh !== nothing && (geom.JBDH[jb] = jbdh)
                if jbdh !== nothing && jbdh >= g.BS[jw] && jbdh <= g.BE[jw]
                    g.KB[id+1] = min(g.KB[g.DHS[jb]], g.KB[id])
                elseif jbdh !== nothing
                    if geom.EL[g.KB[id], id] >= geom.EL[g.KB[g.DHS[jb]], g.DHS[jb]]
                        g.KB[id+1] = g.KB[id]
                    else
                        for k in kt:g.KB[id]
                            if geom.EL[g.KB[g.DHS[jb]], g.DHS[jb]] >= geom.EL[k, id]
                                g.KB[id+1] = k
                                break
                            end
                        end
                    end
                end
            end

            geom.DLX[iu-1] = geom.DLX[iu]
            geom.DLX[id+1] = geom.DLX[id]

            for i in (iu-1):id
                g.KBMIN[i] = min(g.KB[i], g.KB[i+1])
                geom.DLXR[i] = (geom.DLX[i] + geom.DLX[i+1]) * 0.5
            end
            g.KBMIN[id+1] = g.KBMIN[id]
            geom.DLXR[id+1] = geom.DLX[id]
        end
    end
    return (g, geom)
end

"""
    extend_boundary_widths!(g, geom)

init-geom.F90:330-452 -- extends B(K,I) across the one-segment boundary pad
on each end of a branch (copying the adjacent real column, B(1,I)=B(2,I),
and holding B constant below KB(I)).

APPROXIMATED, not faithfully ported, for internal branch junctions
(UH_INTERNAL/DH_INTERNAL and not in the same waterbody): the real Fortran
elevation-matches columns across the junction (init-geom.F90's ELR/ELL/B11
interpolation, ~90 lines). This port instead copies the connected branch's
column directly at the same layer index. Flagged loudly here and in the
module docstring -- this is a real accuracy gap for branched waterbodies
(Detroit's 4-branches-in-1-waterbody topology DOES exercise this), not a
cosmetic shortcut. Fix before trusting boundary-adjacent cell widths.
"""
function extend_boundary_widths!(g, geom)
    kmx = g.KMX
    for jw in 1:g.NWB
        for jb in g.BS[jw]:g.BE[jw]
            iu, id = g.US[jb], g.DS[jb]
            for i in (iu-1):(id+1)
                geom.B[1, i] = geom.B[2, i]
                for k in (g.KB[i]+1):kmx
                    geom.B[k, i] = geom.B[g.KB[i], i]
                end
            end
        end
    end

    for jw in 1:g.NWB
        for jb in g.BS[jw]:g.BE[jw]
            iu, id = g.US[jb], g.DS[jb]
            for k in 1:(kmx-1)
                geom.B[k, iu-1] = geom.B[k, iu]
                if g.UH_INTERNAL[jb]
                    geom.B[k, iu-1] = geom.B[k, g.UHS[jb]]   # approximated -- see docstring
                end
            end
            for k in 1:(kmx-1)
                geom.B[k, id+1] = geom.B[k, id]
                if g.DH_INTERNAL[jb]
                    geom.B[k, id+1] = geom.B[k, g.DHS[jb]]   # approximated -- see docstring
                end
            end
        end
    end
    geom.BNEW .= geom.B
    geom.KBI .= g.KB
    return (g, geom)
end

"""
    compute_areas_volumes!(g, geom)

init-geom.F90:456-630 -- upstream active segment / single-layer adjustment
(CUS, ONE_LAYER, BR_INACTIVE), column widths and areas (BH2/BH/BB/BKT/BI),
average heights and reference widths (AVH2/AVHR/BR/BHR/BHR2), cell volumes
and depths (VOL/DEPTHB/DEPTHM), the initial old-timestep copies (H1/BH1/
BHR1/AVH1), the temporary downstream head segment (CDHS), wind fetch
lengths (FETCHU/FETCHD), and segment heights (HSEG).

NOT PORTED: the TRAPEZOIDAL(JW) branch (init-geom.F90:518-536, uses
GRID_AREA1). Rather than silently computing wrong (RECT-formula) areas for
a TRAP-gridded waterbody, this ERRORS if `geom.TRAPEZOIDAL[jw]` is true --
call `compute_trapezoidal_flag!` (already run by `init_geometry!`) before
this so the check is meaningful, not skipped. See module docstring.
"""
function compute_areas_volumes!(g, geom)
    kmx = g.KMX
    for jw in 1:g.NWB
        if !isempty(geom.TRAPEZOIDAL) && geom.TRAPEZOIDAL[jw]
            error("compute_areas_volumes!: waterbody $jw has GRIDC=\"TRAP\" (trapezoidal grid), " *
                  "but only the RECT-grid area/volume formulas (init-geom.F90:489-517) are ported -- " *
                  "see Core/InitGeometry.jl module docstring. Refusing to silently compute wrong areas.")
        end
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            iu, id = g.US[jb], g.DS[jb]

            if geom.SLOPE[jb] != 0.0
                for i in (iu-1):(id+1)
                    if geom.KBI[i] < kt
                        for k in (geom.KBI[i]+1):kt
                            geom.BNEW[k, i] = 1.0e-6
                        end
                        g.KB[i] = kt
                    end
                end
            end

            iut = iu
            for i in iu:id
                (g.KB[i] - kt < geom.NL[jb] - 1) && (iut = i + 1)
                geom.ONE_LAYER[i] = (kt == g.KB[i])
            end
            for i in (iu-1):id
                g.KBMIN[i] = min(g.KB[i], g.KB[i+1])
            end
            g.KBMIN[id+1] = g.KBMIN[id]

            g.CUS[jb] = iut
            iut >= id && (g.BR_INACTIVE[jb] = true)

            # --- areas and bottom widths (RECT grid only -- see docstring) ---
            for i in (iu-1):(id+1)
                for k in 1:(kmx-1)
                    geom.BH2[k, i] = geom.B[k, i] * geom.H[k, i]
                    geom.BH[k, i] = geom.B[k, i] * geom.H[k, i]
                    geom.BB[k, i] = geom.B[k, i] - (geom.B[k, i] - geom.B[k+1, i]) / (0.5 * (geom.H[k, i] + geom.H[k+1, i])) * geom.H[k, i] * 0.5
                end
                geom.BH[kmx, i] = geom.BH[kmx-1, i]
            end

            for i in (iu-1):(id+1)
                geom.BH2[kt, i] = geom.B[g.KTI[i], i] * (geom.EL[kt, i] - geom.EL[g.KTI[i]+1, i] - geom.Z[i] * geom.COSA[jb]) / geom.COSA[jb]
                kt == g.KTI[i] && (geom.BH2[kt, i] = geom.H2[kt, i] * geom.B[kt, i])
                for k in (g.KTI[i]+1):kt
                    geom.BH2[kt, i] += geom.BNEW[k, i] * geom.H[k, i]
                end
                geom.BKT[i] = geom.BH2[kt, i] / geom.H2[kt, i]
                geom.BI[kt, i] = geom.B[g.KTI[i], i]
            end

            for i in (iu-1):id
                for k in 1:(kmx-1)
                    geom.AVH2[k, i] = (geom.H2[k, i] + geom.H2[k+1, i]) * 0.5
                    geom.AVHR[k, i] = geom.H2[k, i] + (geom.H2[k, i+1] - geom.H2[k, i]) * geom.DLX[i] / (geom.DLX[i] + geom.DLX[i+1])
                end
                geom.AVH2[kmx, i] = geom.H2[kmx, i]
                for k in 1:kmx
                    geom.BR[k, i] = geom.B[k, i] + (geom.B[k, i+1] - geom.B[k, i]) / (0.5 * (geom.DLX[i] + geom.DLX[i+1])) * 0.5 * geom.DLX[i]
                    geom.BHR[k, i] = geom.BH[k, i] + (geom.BH[k, i+1] - geom.BH[k, i]) / (0.5 * (geom.DLX[i] + geom.DLX[i+1])) * 0.5 * geom.DLX[i]
                    geom.BHR2[k, i] = geom.BH2[k, i] + (geom.BH2[k, i+1] - geom.BH2[k, i]) / (0.5 * (geom.DLX[i] + geom.DLX[i+1])) * 0.5 * geom.DLX[i]
                    # CONSTRICTION not ported (no w2_constriction.csv reader) -- always false, no-op here.
                end
            end
            for k in 1:(kmx-1)
                geom.AVH2[k, id+1] = (geom.H2[k, id+1] + geom.H2[k+1, id+1]) * 0.5
                geom.BR[k, id+1] = geom.B[k, id+1]
                geom.BHR[k, id+1] = geom.BH[k, id+1]
            end
            geom.AVH2[kmx, id+1] = geom.H2[kmx, id+1]
            geom.AVHR[kt, id+1] = geom.H2[kt, id+1]
            geom.BHR2[kt, id+1] = geom.BH2[kt, id+1]

            iut = iu
            g.UP_HEAD[jb] && (iut = iu - 1)
            for i in iut:id
                for k in 1:(kmx-1)
                    g.VOL[k, i] = geom.B[k, i] * geom.H2[k, i] * geom.DLX[i]
                end
                g.VOL[kt, i] = geom.BH2[kt, i] * geom.DLX[i]
                geom.DEPTHB[kt, i] = geom.H2[kt, i]
                geom.DEPTHM[kt, i] = geom.H2[kt, i] * 0.5
                for k in (kt+1):kmx
                    geom.DEPTHB[k, i] = geom.DEPTHB[k-1, i] + geom.H2[k, i]
                    geom.DEPTHM[k, i] = geom.DEPTHM[k-1, i] + (geom.H2[k-1, i] + geom.H2[k, i]) * 0.5
                end
            end
        end
    end
    geom.H1 .= geom.H2
    geom.BH1 .= geom.BH2
    geom.BHR1 .= geom.BHR2
    geom.AVH1 .= geom.AVH2

    # Temporary downstream head segment (init-geom.F90:586-594).
    for jb in 1:g.NBR
        if g.DHS[jb] > 0
            jbdh = geom.JBDH[jb]
            if jbdh != 0
                g.CDHS[jb] = max(g.DHS[jb], g.CUS[jbdh])
            end
        end
    end

    # Wind fetch lengths (init-geom.F90:611-616) and segment heights (620-630).
    for jw in 1:g.NWB
        for jb in g.BS[jw]:g.BE[jw]
            for i in g.US[jb]:g.DS[jb]
                geom.FETCHD[i, jb] = geom.FETCHD[i-1 < 1 ? i : i-1, jb] + geom.DLX[i]
            end
            for i in g.DS[jb]:-1:g.US[jb]
                geom.FETCHU[i, jb] = (i == g.DS[jb] ? 0.0 : geom.FETCHU[i+1, jb]) + geom.DLX[i]
            end
            for i in (g.US[jb]-1):(g.DS[jb]+1)
                for k in min(kmx - 1, g.KB[i]):-1:2
                    geom.HSEG[k, i] = geom.HSEG[k+1, i] + geom.H2[k, i]
                end
            end
        end
    end
    return (g, geom)
end

"""
    init_geometry!(g, geom)

Top-level driver -- runs the full sequence above in the same order as the
real INITGEOM. Call after `IO/InputReader.allocate_geometry!` and
`IO/BathymetryReader.read_bathymetry!` have populated DLX/ELWS/B/H for every
waterbody, and after `compute_angles!`/`compute_boundary_flags!` have run
(both included here for convenience, but exposed separately since they only
need SLOPE/UHS/DHS -- not the bathymetry data -- and might be useful earlier
in a future input-validation pass).
"""
function init_geometry!(g, geom)
    net = build_branch_network(g)
    compute_angles!(geom)
    compute_boundary_flags!(g, geom)
    compute_trapezoidal_flag!(geom)
    allocate_init_geometry!(g, geom)
    compute_layer_elevations!(g, geom)
    compute_water_surface_layer!(g, geom)
    compute_bottom_layers!(g, geom, net)
    extend_boundary_widths!(g, geom)
    compute_areas_volumes!(g, geom)
    return (g, geom, net)
end
