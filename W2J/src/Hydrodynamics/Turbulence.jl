# ==============================================================================
# Hydrodynamics/Turbulence.jl
#
# Real port of az.f90's CALCULATE_TKE (2026-08-24) -- the turbulent-kinetic-
# energy closure real DET's own control file specifies (`AZC='TKE'`,
# confirmed by reading w2_con.csv directly, not assumed -- see CLAUDE.md).
# Motivated by the first real Fortran-output comparison (2026-08-23/24):
# this port's DZ (vertical eddy diffusivity, consumed by Hydrodynamics/
# Transport.jl) was a caller-set uniform constant, the leading suspect for
# the ELWS/T2 divergence found in that comparison.
#
# SCOPE: only `CALCULATE_TKE` is ported, not `CALCULATE_TKE1` (the
# Chapman-Cole/Gould variant, `AZC='TKE1'`) or the algebraic closures
# (NICK/RNG/PARAB/W2N, the non-TKE branch of `CALCULATE_AZ`) -- real DET
# uses `AZC='TKE'` specifically, confirmed by reading `w2_con.csv` row 125
# directly (context: MANN->FRICC, 0.001->Z0, then the real AZC/AZSLC/AZMAX/
# TKEBC/EROUGH/ARODI read order from input.F90:832-837 matches exactly).
#
# REDUCED PHYSICS:
# - The surface boundary-condition term (`USTAR`, wind-driven) is now REAL
#   as of 2026-09-10 -- `IO/MetReader.jl` reads real `MET.csv` wind speed
#   and `compute_wind_stress!` derives `WIND10(I)`/`CZ(I)` from it (see
#   that function's own docstring for its own, smaller reductions: no
#   log-law height rescale, no fetch correction, `WSC(I)=1`). Before this,
#   `USTAR` was hardcoded to 0 (confirmed with user 2026-08-24 as "port TKE
#   now, MET later").
# - The bottom boundary-condition term and the interior lateral-friction
#   production terms (`PRHE`/`PRHK`) all need `GC2 = G/FRIC(I)^2`, which
#   needs `FRIC(I)` (per-segment bottom friction, Tier 1, not read) or
#   `MANNINGS_N(JW)` (also not read). `GC2 = 0` here -- NOTE this is not
#   an invented shortcut: real Fortran's own `IF (FRIC(I) /= 0.0) GC2 = ...`
#   branch literally produces `GC2 = 0.0` whenever `FRIC(I) == 0.0`, so
#   this reduction lands on a real, Fortran-supported code path, not a new
#   one. STILL NOT PORTED (unlike wind, above).
# - Net effect: TKE production now comes from BOTH interior velocity shear
#   (`PRDK`) AND real wind-driven surface production, damped by buoyancy
#   (`BOUK`) and dissipation -- still no bottom-friction turbulence
#   generation (needs `FRIC`, above). Wind is usually the dominant real
#   surface-mixing mechanism in a reservoir, so this closes the single
#   biggest gap flagged when this file was first written.
#
# NOT PORTED: `CALCULATE_TKE1`, the non-TKE algebraic closures, `TKELATPRD`
# lateral-friction production (needs the same missing FRIC data),
# `STRICKON`/`CALCFRIC` (Strickler roughness-adjustment, needs FRIC too).
#
# NUMERICAL-STABILITY STOPGAP (added 2026-08-24, same day as the port
# above, NOT a physical value -- confirmed with user via AskUserQuestion):
# with no wind/friction forcing, TKE correctly decays to its absolute
# floor, driving `DZ` down to ~2e-7 -- far below the previous placeholder
# constant (`dz_const=1e-3`, a caller-set value with no physical basis
# either, but empirically large enough to numerically stabilize
# `Hydrodynamics/Transport.jl`'s explicit UPWIND advection scheme). Without
# that accidental damping, a real DET-forced `run_forced_simulation!` run
# blew `T2` up to >99°C around `jday≈4.125` -- confirmed NOT a CFL/momentum
# instability (traced: `compute_curmax` stayed comfortably high throughout,
# and wiring `step_hydrodynamics_adaptive!` into the driver the same day
# did NOT prevent this blow-up, since that check has no thermal-transport
# stability criterion at all -- a real, separate, still-OPEN gap, not
# fixed here). `DZ_STABILITY_FLOOR` restores the old constant's magnitude
# as an explicit floor UNDER the real TKE-derived value, so `calculate_
# tke!`'s spatially-varying, shear-responsive `DZ` is still used whenever
# it exceeds this floor -- only clamped when it would otherwise destabilize
# the (still explicit, still un-CFL-checked) advection scheme. Remove this
# once either MET/FRIC give TKE real forcing (so its floor value stops
# being pathologically small) or a real thermal-transport stability check
# is added to the adaptive timestep.
# ==============================================================================

const AZMIN = 1.4e-6      # w2modules.F90:101 "AZMIN=1.4D-6"
const TKEMIN1 = 1.25e-7   # az.f90:8 "TKEMIN1=1.25D-7"
const TKEMIN2 = 1.0e-9    # az.f90:8 "TKEMIN2=1.0D-9"
const FRAZDZ = 0.14       # w2modules.F90:100 "FRAZDZ=0.14D0"
const DZMIN = 1.4e-7      # w2modules.F90:100 "DZMIN=1.4D-7"
const DZ_STABILITY_FLOOR = 1.0e-3  # NOT a real Fortran constant -- numerical-stability stopgap, see module docstring
const TKE_SIG = (1.0, 1.3) # az.f90:20 "SIG(1)=1.0D0; SIG(2)=1.3D0"

"""
    allocate_turbulence_state!(g)

Sizes `g.AZ`/`g.DZ`/`g.TKE`/`g.AZT`. `TKE` is initialized to
`(TKEMIN1, TKEMIN2)` (not zero) -- `calculate_tke!`'s dissipation-ratio
term `TKE(K,I,2)/TKE(K,I,1)` would divide by zero on the very first call
otherwise, before the real per-step `MAX(...,TKEMIN)` floor has ever run.
`AZ` is initialized to `AZMIN` for the same reason (`compute_curmax`/
`update_velocities!` and this file's own shear-production term all read
`AZ` before `calculate_tke!` necessarily has written every entry it will
ever write, e.g. `AZ[kb,i]` before the first call).

Also lazily sizes `g.AT`/`g.CT`/`g.VT`/`g.DT` if empty -- `calculate_tke!`
reuses these as its TRIDIAG buffer (see that function's docstring), but
they're normally sized by `Hydrodynamics/Transport.jl`'s `allocate_
transport_state!`, which a caller that only wants hydrodynamics (not
temperature transport) never calls. Found as a real `BoundsError` (`0×0
Matrix` at a real index) across every test calling `hydrodynamic_step!`
without also calling `allocate_transport_state!` -- same failure SHAPE as
the `DIST_TRIBS`/`PLACE_QIN` unconditional-inner-call bugs earlier this
session, different root cause (a cross-module shared-buffer dependency,
not a missing default). Safe to size here even if `allocate_transport_
state!` runs later too -- both just zero/allocate fresh arrays, no data to
lose at this point.
"""
function allocate_turbulence_state!(g)
    kmx, imx = g.KMX, g.IMX
    g.AZ = fill(AZMIN, kmx, imx)
    g.DZ = fill(DZ_STABILITY_FLOOR, kmx, imx)
    g.TKE = zeros(Float64, kmx, imx, 2)
    g.TKE[:, :, 1] .= TKEMIN1
    g.TKE[:, :, 2] .= TKEMIN2
    g.AZT = fill(AZMIN, kmx, imx)
    isempty(g.AT) && (g.AT = zeros(Float64, kmx, imx))
    isempty(g.CT) && (g.CT = zeros(Float64, kmx, imx))
    isempty(g.VT) && (g.VT = zeros(Float64, kmx, imx))
    isempty(g.DT) && (g.DT = zeros(Float64, kmx, imx))
    # AZMAX(JW) (per-waterbody eddy-viscosity clamp, part of the same
    # AZC/AZSLC/AZMAX/TKEBC/EROUGH/ARODI block as AZC itself) is Tier 1,
    # not read by InputReader.jl -- defaults to a generous upper bound
    # (100.0, matching this port's existing DZMAX=100.0 Fortran-default
    # convention, see W2Global's own constructor comment) rather than a
    # narrow guess, since AZMAX is a pure safety CLAMP, not a physics
    # on/off switch like DIST_TRIBS/PLACE_QIN -- a generous default is
    # unlikely to ever bind, not a silent accuracy gap in the same way.
    isempty(g.AZMAX) && (g.AZMAX = fill(100.0, g.NWB))
    # Meteorology (IO/MetReader.jl) -- same "unconditional inner call"
    # reasoning as DIST_TRIBS/PLACE_QIN: compute_wind_stress! runs inside
    # every hydrodynamic_step! call (via calculate_tke!), so these must be
    # sized here regardless of whether a caller ever loads real MET data.
    # WIND=0 is the safe default (no wind forcing unless a caller explicitly
    # loads+updates real MET conditions) -- same "off by default" discipline
    # as QIN/QOT/QDTR. WSC=1.0 (no sheltering) is the real Fortran default
    # absent a per-segment wind-sheltering file (Tier 1, not read).
    isempty(g.WIND) && (g.WIND = zeros(Float64, g.NWB))
    isempty(g.PHI) && (g.PHI = zeros(Float64, g.NWB))
    isempty(g.TAIR) && (g.TAIR = zeros(Float64, g.NWB))
    isempty(g.TDEW) && (g.TDEW = zeros(Float64, g.NWB))
    isempty(g.CLOUD) && (g.CLOUD = zeros(Float64, g.NWB))
    isempty(g.SRO) && (g.SRO = zeros(Float64, g.NWB))
    isempty(g.WIND10) && (g.WIND10 = zeros(Float64, imx))
    isempty(g.CZ) && (g.CZ = zeros(Float64, imx))
    isempty(g.WSC) && (g.WSC = fill(1.0, imx))
    return g
end

"""
    compute_wind_stress!(g, geom)

`w2_4_win.f90:601-632`'s "Adjusted wind speed and surface wind shear drag
coefficient" block -- computes `WIND10(I)` (wind speed adjusted toward a
10m reference) and `CZ(I)` (the wind drag coefficient), the two values
`calculate_tke!`'s `USTAR` term consumes.

REDUCED PHYSICS: the real `WIND10(I) = WIND(JW)*WSC(I)*ln(10/Z0(JW))/
ln(WINDH(JW)/Z0(JW))` log-law height rescale needs `WINDH(JW)` (the
height at which wind was measured) and `Z0(JW)` (surface roughness) --
both Tier 1, not read by `InputReader.jl`. Here `WIND10(I) = WIND(JW)*
WSC(I)` (the rescale factor dropped, i.e. assumed `== 1`), a reasonable
assumption for real DET's `MET.csv` (a NASA POWER API export -- POWER's
WS10M product is already a 10m reference value, so a real WINDH=10/Z0
rescale would very nearly cancel anyway). The `FETCH_CALC(JW)` fetch-based
correction (needs per-segment fetch distances, Tier 1) is also NOT
applied. `CZ(I)`'s own piecewise formula has no Tier-1 dependency and is
ported exactly.

Call once per timestep, before `calculate_tke!` (needs current `WIND10`/
`CZ`) and after `IO/MetReader.jl`'s `update_met_conditions!` (needs
current `g.WIND`).
"""
function compute_wind_stress!(g, geom)
    for jw in 1:g.NWB
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            iu, id = g.CUS[jb], g.DS[jb]
            for i in max(1, iu-1):min(g.IMX, id+1)
                wind10 = g.WIND[jw] * g.WSC[i]
                g.WIND10[i] = wind10
                g.CZ[i] = if wind10 >= 15.0
                    0.0026
                elseif wind10 >= 4.0
                    0.0005 * sqrt(wind10)
                elseif wind10 >= 0.5
                    0.0044 * wind10^(-1.15)
                else
                    0.01
                end
            end
        end
    end
    return g
end

"""
    calculate_tke!(g, geom, dlt)

az.f90's `CALCULATE_TKE`, ported directly (see module docstring for the
real-vs-reduced breakdown). Per active segment `i`: sets the top/bottom
boundary-condition TKE values (real formula, but structurally 0 here --
see module docstring), advances the interior TKE/dissipation via explicit
production/dissipation terms, solves the implicit vertical-diffusion
system for both TKE and dissipation via `Solvers/Tridiagonal.jl`'s
`thomas_solve!` (reusing `g.AT`/`g.CT`/`g.VT`/`g.DT`, the same shared
buffer `Hydrodynamics/Transport.jl` uses -- safe since nothing runs
concurrently with this function and `temperature_transport!`/
`constituent_transport!` each recompute those arrays fresh before reading
them), then derives `AZ`/`DZ` from the solved TKE/dissipation
(`AZT = 0.09*TKE1^2/TKE2`, `AZ` = interface-averaged `AZT` clamped to
`[AZMIN, AZMAX(jw)]`, `DZ = max(DZMIN, FRAZDZ*AZ)`).

MUST be called after `compute_density_field!` (needs current `g.RHO` for
the buoyancy term `BOUK`) and after real `U` is current for this timestep
(needs `g.U` for shear production) -- i.e. after `update_velocities!` in
`hydrodynamic_step!`'s existing order, or equivalently at the START of a
driver's next call before `apply_temperature_sources!`/
`temperature_transport!` read `DZ`. NOT yet wired into any driver -- same
"port the mechanism, then wire it into the driver" two-step pattern as
`apply_temperature_sources!` earlier this session.
"""
function calculate_tke!(g, geom, dlt)
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            iu, id = g.CUS[jb], g.DS[jb]
            for i in iu:id
                kb = g.KB[i]

                gc2 = 0.0     # FRIC(I)/MANNINGS_N(JW) still not read -- see module docstring
                # USTAR: real formula, now using real WIND10(I)/CZ(I) (IO/MetReader.jl +
                # compute_wind_stress!, 2026-09-10) -- previously always 0 (no MET at all).
                ustar = sqrt(1.25 * g.CZ[i] * g.WIND10[i]^2 / g.RHO[kt, i])
                ustarbkt = sqrt(gc2) * abs(0.5 * (g.U[kt, i] + g.U[kt, i-1]))
                g.TKE[kt, i, 1] = (3.33 * (ustar^2 + ustarbkt^2)) * geom.BHRATIO[kt, i]
                g.TKE[kt, i, 2] = (ustar^3 + ustarbkt^3) * 5.0 / geom.H1[kt, i] * geom.BHRATIO[kt, i]

                for k in (kt+1):(kb-1)
                    bouk = max(g.AZ[k, i] * G_GRAVITY * (g.RHO[k+1, i] - g.RHO[k, i]) / (geom.H[k, i] * RHOW), 0.0)
                    prdk = g.AZ[k, i] * (0.5 * (g.U[k, i] + g.U[k, i-1] - g.U[k+1, i] - g.U[k+1, i-1]) /
                                          (geom.H[k, i] * 0.5 + geom.H[k+1, i] * 0.5))^2
                    unst = prdk - g.TKE[k, i, 2]
                    unse = 1.44 * g.TKE[k, i, 2] / g.TKE[k, i, 1] * prdk -
                           1.92 * (g.TKE[k, i, 2] / g.TKE[k, i, 1] * g.TKE[k, i, 2])
                    g.TKE[k, i, 1] += dlt * (unst - bouk)
                    g.TKE[k, i, 2] += dlt * unse
                end

                ustarb = sqrt(gc2) * abs(0.5 * (g.U[kb, i] + g.U[kb, i-1]))
                g.TKE[kb, i, 1] = 0.5 * (3.33 * ustarb^2 + g.TKE[kb, i, 1])
                g.TKE[kb, i, 2] = 0.5 * (ustarb^3 * 5.0 / geom.H[kb, i] + g.TKE[kb, i, 2])

                for j in 1:2
                    k = kt
                    g.AT[k, i] = 0.0; g.CT[k, i] = 0.0; g.VT[k, i] = 1.0; g.DT[k, i] = g.TKE[k, i, j]
                    for k in (kt+1):(kb-1)
                        g.AT[k, i] = -dlt / geom.BH1[k, i] * geom.BB[k-1, i] / TKE_SIG[j] * g.AZ[k-1, i] / geom.AVH1[k-1, i]
                        g.CT[k, i] = -dlt / geom.BH1[k, i] * geom.BB[k, i] / TKE_SIG[j] * g.AZ[k, i] / geom.AVH1[k, i]
                        g.VT[k, i] = 1.0 - g.AT[k, i] - g.CT[k, i]
                        g.DT[k, i] = g.TKE[k, i, j]
                    end
                    k = kb
                    g.AT[k, i] = 0.0; g.CT[k, i] = 0.0; g.VT[k, i] = 1.0; g.DT[k, i] = g.TKE[k, i, j]
                    thomas_solve!((@view g.AT[kt:kb, i]), (@view g.VT[kt:kb, i]), (@view g.CT[kt:kb, i]),
                                  (@view g.DT[kt:kb, i]), (@view g.TKE[kt:kb, i, j]))
                end

                for k in kt:kb
                    g.TKE[k, i, 1] = max(g.TKE[k, i, 1], TKEMIN1)
                    g.TKE[k, i, 2] = max(g.TKE[k, i, 2], TKEMIN2)
                    g.AZT[k, i] = 0.09 * g.TKE[k, i, 1]^2 / g.TKE[k, i, 2]
                end
                for k in kt:(kb-1)
                    az = 0.5 * (g.AZT[k, i] + g.AZT[k+1, i])
                    az = clamp(az, AZMIN, g.AZMAX[jw])
                    g.AZ[k, i] = az
                    g.DZ[k, i] = max(DZ_STABILITY_FLOOR, FRAZDZ * az)  # DZMIN alone is too small to be numerically stable here -- see module docstring
                end
                g.AZ[kb, i] = AZMIN
                g.AZT[kb, i] = AZMIN
            end
        end
    end
    return g
end
