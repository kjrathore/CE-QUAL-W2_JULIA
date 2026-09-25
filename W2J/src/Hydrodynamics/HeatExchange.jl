# ==============================================================================
# Hydrodynamics/HeatExchange.jl
#
# Real port of heat-exchange.f90's EQUILIBRIUM_TEMPERATURE-method surface
# heat budget + temperature.F90's solar-penetration/sediment-exchange
# coupling -- 2026-09-23, motivated directly by the ~8C summer cold bias
# found comparing a full 2017 DET run against real observed data (see
# CLAUDE.md "FOURTH fix..." section): this port had NO surface heat source
# at all before this file (SB/ST flagged "NOT YET COMPUTED" since
# FreeSurface.jl's very first reduced-physics scoping in 2026-08-12), so
# the model could only ever gain heat from boundary inflow.
#
# SCOPE, confirmed with user via AskUserQuestion (2026-09-23):
# - SHORT_WAVE_RADIATION: the SIMPLE cloud-based formula only (solar
#   altitude from latitude/longitude/declination, corrected by CLOUD^2),
#   NOT the full MEEUS/Bird astronomical model (~300 extra lines of solar-
#   position/atmospheric-turbidity code, heat-exchange.f90:183-496) --
#   real Fortran's I_SW flag selects between them; this port always uses
#   the simple branch.
# - EQUILIBRIUM_TEMPERATURE method (a single lumped ET/CSHE linearization,
#   temperature.F90:114-116), NOT the full TERM_BY_TERM budget
#   (SURFACE_TERMS' separate RE/RC/RB terms, heat-exchange.f90:116-167) --
#   fewer real formulas to port and validate, matching the user's "simple"
#   preference. This means RH_EVAP/SURFACE_TERMS are NOT ported here.
# - Solar penetration through the water column (Beer-Lambert attenuation,
#   temperature.F90:120-147) IS ported -- it's directly coupled to the
#   same code block (the equilibrium method implicitly assumes ALL
#   shortwave is absorbed at the surface, so this correction is needed to
#   avoid double-counting/mis-locating solar heating vertically) and isn't
#   materially extra scope.
# - Sediment/water heat exchange (temperature.F90:150-160) IS ported -- a
#   small, separate additive term already sitting in the same loop.
#
# TIER-1 IO GAP, same discipline as THETA/UPWIND/PLACE_QIN (confirmed with
# user): the real w2_con.csv row layout for SLHTC/SROC/RHEVC/METIC/FETCHC
# is AMBIGUOUS in DET's actual file (8 data rows found where the Fortran
# source's read order -- input.F90:793-801 -- expects 9; one flag appears
# dropped from this specific CSV export, and guessing which risks silently
# wrong physics). Per user's explicit choice, these are NOT parsed here;
# instead:
#   - AFW/BFW/CFW/WINDH ARE CSV-confirmed for DET: found immediately after
#     the ambiguous flags (9.2, 0.46, 2, 6) and independently corroborated
#     as recognizable standard real values for these coefficients (not
#     guessed) -- lazily defaulted to these in `compute_equilibrium_
#     temperature!`.
#   - RH_EVAP(JW) is inferred FALSE (not read) from AFW/BFW/CFW's presence
#     -- if RH_EVAP were true, those values would be unused; real Fortran
#     wouldn't provide specific calibrated-looking values for an unused
#     path (see SURFACE_TERMS' `IF (RH_EVAP(JW))` branch, NOT ported here
#     anyway since only EQUILIBRIUM_TEMPERATURE is used).
#   - EXH2O (light extinction), BETA (fraction of shortwave NOT absorbed
#     at the immediate surface), CC_SW (cloud-cover correction to
#     shortwave), CBHE/TSED/TSEDF (sediment exchange coefficient/
#     temperature/absorbed-fraction) have NO confirmed real DET value --
#     lazily defaulted to flagged PLACEHOLDER values (see each default
#     below) pending real Tier-1 IO or, fittingly, Enzyme.jl-based
#     calibration -- these are exactly the kind of coefficients a real W2
#     application calibrates by hand anyway, so approximate-then-calibrate
#     is a coherent plan, not a shortcut.
# - CC_SW's real read path (time-varying-data.f90:162, a SEPARATE small
#   file this port doesn't read) confirms it's genuinely not part of the
#   main w2_con.csv at all -- not a gap in THIS file's parsing, a whole
#   separate Tier-1 input this port doesn't have.
#
# WIND2(I): real Fortran rescales WIND(JW) from the met-station height
# WINDH(JW) to a 2m reference via a log law using Z0(JW) (roughness
# length, w2_4_win.f90:433) -- `Hydrodynamics/Turbulence.jl`'s
# `compute_wind_stress!` already skips this exact rescale for WIND10
# (assumes input wind is already near-surface-reference), so WIND2 here
# just reuses WIND10 rather than re-deriving a second, inconsistent
# reduction.
#
# NOT PORTED: ICE(I) (ice cover suppresses heat exchange entirely in real
# Fortran, `IF (.NOT. ICE(I))`) -- this port has no ice model, so heat
# exchange always runs, even below 0C. A real difference under genuinely
# cold conditions, not expected to matter for DET's 2017 comparison
# (surface rarely stays below freezing long in the observed data) but
# flagged, not silently assumed harmless.
# ==============================================================================

const RHOWCP = 1000.0 * 4186.0   # init.F90:732 "RHOWCP = RHOW*CP", CP=4186 J/kg/degC (w2modules.F90)
const W_M2_TO_BTU_FT2_DAY = 7.60796
const BTU_FT2_DAY_TO_W_M2 = 0.1314
const FLUX_BR_TO_FLUX_SI = 0.23659

deg_f(c) = c * 1.8 + 32.0
deg_c(f) = (f - 32.0) * 5.0 / 9.0

"""
    compute_short_wave_radiation!(g, geom, jday)

heat-exchange.f90:25-57's non-MEEUS branch -- solar altitude from
latitude/longitude/day-of-year, corrected by cloud cover, into
`g.SRON[jw]` (W/m^2). Real Fortran's `Met_Regions` branch not ported
(this port has no region concept, matching every other MetReader.jl
consumer). `geom.CC_SW` lazily defaulted to 0.0 (PLACEHOLDER -- see
module docstring; real Fortran reads this from a separate file this port
doesn't have) -- with CC_SW=0, cloud cover has NO effect on shortwave,
a real, flagged simplification, not a silent one.
"""
function compute_short_wave_radiation!(g, geom, jday)
    isempty(geom.CC_SW) && (geom.CC_SW = fill(0.0, g.NWB))
    isempty(g.SRON) && (g.SRON = zeros(Float64, g.NWB))
    for jw in 1:g.NWB
        local_lon = geom.LONGIT[jw]
        standard = 15.0 * floor(local_lon / 15.0)
        hour = (jday - floor(jday)) * 24.0
        iday = jday - floor(jday / 365.0) * 365.0
        taud = 2.0 * pi * (iday - 1.0) / 365.0
        eqtnew = 0.170 * sin(4.0 * pi * (iday - 80.0) / 373.0) - 0.129 * sin(2.0 * pi * (iday - 8.0) / 355.0)
        hh = 0.261799 * (hour - (local_lon - standard) * 0.0666667 + eqtnew - 12.0)
        decl = 0.006918 - 0.399912 * cos(taud) + 0.070257 * sin(taud) - 0.006758 * cos(2taud) +
               0.000907 * sin(2taud) - 0.002697 * cos(3taud) + 0.001480 * sin(3taud)
        lat_rad = geom.LAT[jw] * 0.0174533
        sinal = sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(hh)
        a0 = 57.2957795 * asin(clamp(sinal, -1.0, 1.0))
        cloud = g.CLOUD[jw]
        if a0 > 0.0
            g.SRON[jw] = (1.0 - geom.CC_SW[jw] * cloud * cloud) * 24.0 *
                         (2.044a0 + 0.1296a0^2 - 1.941e-3a0^3 + 7.591e-6a0^4) * BTU_FT2_DAY_TO_W_M2
        else
            g.SRON[jw] = 0.0
        end
    end
    return g
end

"""
    compute_equilibrium_temperature!(g, geom, jw)

heat-exchange.f90:66-110 -- the classic Edinger et al. equilibrium-
temperature linearization: iterates `ET`/`CSHE` (equilibrium temperature,
heat-exchange coefficient) to a fixed point (real Fortran's own
convergence tolerance/iteration cap: `|ETP-ET|<=0.05`, max 10 iterations).
Computed once per waterbody here (not per-segment like real Fortran's
`ET(I)`/`CSHE(I)`) since every input (TDEW/TAIR/SRON per-waterbody,
WIND2 uniform under this port's WSC=1.0 reduction) is currently uniform
across a waterbody's segments -- `apply_surface_heat_exchange!` broadcasts
the result to every segment, flagged as a reduction that would need
revisiting if WSC (wind-sheltering) is ever made genuinely per-segment.

`geom.AFW`/`BFW`/`CFW`/`WINDH` lazily defaulted to DET's real CSV-
confirmed values (9.2, 0.46, 2.0, 6.0 -- see module docstring). RH_EVAP
inferred false (not ported): always uses the `FW = ACONV*AFW+BCONV*BFW*
WIND2^CFW` wind function, matching real Fortran's `ELSE` branch.

`geom.CSHE_MULT` -- added 2026-09-24, a flagged EMPIRICAL CORRECTION, NOT
a real Fortran mechanism (lazily defaulted to 1.0, i.e. inert unless a
caller sets it). Found via a real, decisive diagnostic against observed
DET 2017 surface temperature: amplifying `CSHE` alone (leaving `ET`
untouched -- multiplying only the FINAL output, not the internal ET
fixed-point iteration, so the equilibrium TARGET is unaffected, only the
RATE of approach to it) by 8x collapsed the surface bias from -7.44degC
to -0.02degC and surface RMSE from 7.96degC to 2.35degC -- confirming the
heat-exchange RATE, not the solar input magnitude (checked and ruled out
separately -- the simple clear-sky formula was NOT found to systematically
underestimate real observed solar), was the dominant remaining error
source. Root mechanism NOT fully confirmed (checked the `FLUX_BR_TO_
FLUX_SI` unit conversion by hand -- numerically correct, not a units bug;
the more likely culprit is `AFW`/`BFW`/`CFW`, inferred from an ambiguous
w2_con.csv parse in the first heat-exchange session and never fully
confirmed) -- `CSHE_MULT` is deliberately a calibratable multiplier
instead of a guessed replacement value, to be set via real Bayesian
optimization (`tools/calibration/bayesopt_calibrate.jl`), not hardcoded.
"""
function compute_equilibrium_temperature!(g, geom, jw)
    isempty(geom.AFW) && (geom.AFW = fill(9.2, g.NWB))
    isempty(geom.BFW) && (geom.BFW = fill(0.46, g.NWB))
    isempty(geom.CFW) && (geom.CFW = fill(2.0, g.NWB))
    isempty(geom.WINDH) && (geom.WINDH = fill(6.0, g.NWB))
    isempty(geom.CSHE_MULT) && (geom.CSHE_MULT = fill(1.0, g.NWB))
    isempty(g.WIND2) && (g.WIND2 = zeros(Float64, g.IMX))
    isempty(g.ET) && (g.ET = zeros(Float64, g.IMX))
    isempty(g.CSHE) && (g.CSHE = zeros(Float64, g.IMX))

    iu, id = g.CUS[g.BS[jw]], g.DS[g.BE[jw]]
    tdew_f = deg_f(g.TDEW[jw])
    tair_f = deg_f(g.TAIR[jw])
    sro_br = g.SRON[jw] * W_M2_TO_BTU_FT2_DAY   # SHADE=1.0 (no shading, Tier 1 not ported elsewhere either)

    aconv = W_M2_TO_BTU_FT2_DAY
    cfw = geom.CFW[jw]
    bconv = cfw == 1.0 ? 3.401062 : (cfw == 2.0 ? 1.520411 : error("compute_equilibrium_temperature!: CFW(jw=$jw)=$cfw has no known BCONV (only 1.0/2.0 ported, matching heat-exchange.f90's own gap: \"CFW not determined for other values\")"))

    for i in iu:id
        g.WIND2[i] = g.WIND10[i]
        wind2 = g.WIND2[i]
        fw = aconv * geom.AFW[jw] + bconv * geom.BFW[jw] * wind2^cfw

        et = tdew_f
        tstar = (et + tdew_f) * 0.5
        beta_coef = 0.255 - 8.5e-3 * tstar + 2.04e-4 * tstar * tstar
        cshe = 15.7 + (0.26 + beta_coef) * fw
        ra = 3.1872e-8 * (tair_f + 459.67)^4
        etp = (sro_br + ra - 1801.0) / cshe + (cshe - 15.7) * (0.26 * tair_f + beta_coef * tdew_f) / (cshe * (0.26 + beta_coef))
        j = 0
        while abs(etp - et) > 0.05 && j < 10
            et = etp
            tstar = (et + tdew_f) * 0.5
            beta_coef = 0.255 - 8.5e-3 * tstar + 2.04e-4 * tstar * tstar
            cshe = 15.7 + (0.26 + beta_coef) * fw
            etp = (sro_br + ra - 1801.0) / cshe + (cshe - 15.7) * (0.26 * tair_f + beta_coef * tdew_f) / (cshe * (0.26 + beta_coef))
            j += 1
        end
        g.ET[i] = deg_c(etp)
        g.CSHE[i] = cshe * FLUX_BR_TO_FLUX_SI / RHOWCP * geom.CSHE_MULT[jw]
    end
    return g
end

"""
    apply_surface_heat_exchange!(g, geom)

temperature.F90:100-160's surface + solar-penetration + sediment coupling
into `g.TSS`, using the just-computed `ET`/`CSHE`/`SRON` and last step's
temperature (`g.HYD[:,:,4]`, matching `apply_temperature_sources!`'s own
`cold` convention). MUST be called after `compute_short_wave_radiation!`
and `compute_equilibrium_temperature!` for the same step, and (like
`apply_temperature_sources!`) before `temperature_transport!`.

`geom.EXH2O`/`BETA`/`CBHE`/`TSED`/`TSEDF` lazily defaulted to flagged
PLACEHOLDER values (no confirmed real DET value, see module docstring):
`EXH2O=0.35` (1/m, a commonly-cited mid-range CE-QUAL-W2 default for a
moderately turbid reservoir), `BETA=0.6` (a commonly-cited default
fraction NOT absorbed at the immediate surface), `CBHE=0.3` (W/m^2/degC,
a commonly-cited mild sediment heat-exchange default), `TSED=12.0`degC
(a plausible annual-mean-ish placeholder), `TSEDF=0.0` (no direct solar
absorption into sediment). `geom.GAMMA[k,i]` filled uniformly as
`EXH2O[jw]` per waterbody (temperature.F90:15's `READ_EXTINCTION`
branch, reduced -- no per-cell algae/organic-matter extinction
contribution, matching `EXOM`/`EXSS` not being ported).
"""
function apply_surface_heat_exchange!(g, geom)
    isempty(geom.EXH2O) && (geom.EXH2O = fill(0.35, g.NWB))
    isempty(geom.BETA) && (geom.BETA = fill(0.6, g.NWB))
    isempty(geom.CBHE) && (geom.CBHE = fill(0.3, g.NWB))
    isempty(geom.TSED) && (geom.TSED = fill(12.0, g.NWB))
    isempty(geom.TSEDF) && (geom.TSEDF = fill(0.0, g.NWB))
    if size(geom.GAMMA) != (g.KMX, g.IMX)
        geom.GAMMA = zeros(Float64, g.KMX, g.IMX)
    end

    cold = @view g.HYD[:, :, 4]
    for jw in 1:g.NWB
        kt = g.KTWB[jw]
        for jb in g.BS[jw]:g.BE[jw]
            g.BR_INACTIVE[jb] && continue
            iu, id = g.CUS[jb], g.DS[jb]
            for i in iu:id
                geom.GAMMA[kt:g.KB[i], i] .= geom.EXH2O[jw]

                heatex = (g.ET[i] - cold[kt, i]) * g.CSHE[i] * geom.BI[kt, i] * geom.DLX[i]
                g.TSS[kt, i] += heatex

                sroout = (1.0 - geom.BETA[jw]) * (g.SRON[jw] / RHOWCP) * geom.BI[kt, i] * geom.DLX[i] *
                         exp(-geom.GAMMA[kt, i] * geom.DEPTHB[kt, i])
                g.TSS[kt, i] -= sroout

                kb = g.KB[i]
                # NOTE: geom.BI is only ever populated at row `kt` (the
                # partial-width top layer) -- Core/InitGeometry.jl leaves
                # BI[k,i] for k>kt at its zeros(...) init value forever.
                # Real Fortran's BI(K,I) for K>KT genuinely equals B(K,I)
                # (only the top layer is ever partial), so B is used below
                # for every k>kt reference, not the all-zero BI.
                sroSed = kt == kb ? sroout * geom.TSEDF[jw] :
                         sroout * (1.0 - geom.B[kt+1, i] / geom.BI[kt, i]) * geom.TSEDF[jw]
                g.TSS[kt, i] += sroSed

                sroin = kt == kb ? 0.0 : sroout * geom.B[kt+1, i] / geom.BI[kt, i]
                for k in (kt+1):kb
                    sroout_k = sroin * exp(-geom.GAMMA[k, i] * geom.H1[k, i])
                    sronet = sroin - sroout_k
                    sroSed_k = k != kb ? sroout_k * (1.0 - geom.B[k+1, i] / geom.B[k, i]) * geom.TSEDF[jw] :
                               sroout_k * geom.TSEDF[jw]
                    g.TSS[k, i] += sronet + sroSed_k
                    sroin = k != kb ? sroout_k * geom.B[k+1, i] / geom.B[k, i] : 0.0
                end

                # Sediment/water exchange, every active layer -- BI only
                # meaningful at kt (see note above), B used for k>kt.
                tflux_kt = kt == kb ?
                    geom.CBHE[jw] / RHOWCP * (geom.TSED[jw] - cold[kt, i]) * geom.BI[kt, i] * geom.DLX[i] :
                    geom.CBHE[jw] / RHOWCP * (geom.TSED[jw] - cold[kt, i]) * (geom.BI[kt, i] - geom.B[kt+1, i]) * geom.DLX[i]
                g.TSS[kt, i] += tflux_kt
                for k in (kt+1):kb
                    tflux = k == kb ?
                        geom.CBHE[jw] / RHOWCP * (geom.TSED[jw] - cold[k, i]) * geom.B[k, i] * geom.DLX[i] :
                        geom.CBHE[jw] / RHOWCP * (geom.TSED[jw] - cold[k, i]) * (geom.B[k, i] - geom.B[k+1, i]) * geom.DLX[i]
                    g.TSS[k, i] += tflux
                end
            end
        end
    end
    return g
end
