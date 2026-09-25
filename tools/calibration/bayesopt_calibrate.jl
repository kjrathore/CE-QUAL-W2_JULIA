# ==============================================================================
# tools/calibration/bayesopt_calibrate.jl
#
# Gradient-free Bayesian calibration of the heat-exchange placeholder
# coefficients flagged in `Hydrodynamics/HeatExchange.jl` (EXH2O, BETA,
# CBHE, TSED, CSHE_MULT) against real observed DET 2017 temperature-at-
# depth data.
#
# CSHE_MULT added 2026-09-24: a real, decisive diagnostic (surface-only
# bias/RMSE broken out from the aggregate, then an amplification test)
# found the heat-exchange RATE coefficient (CSHE), not the solar input
# magnitude, was the dominant remaining error source -- amplifying CSHE
# alone by 8x (a blunt manual test) collapsed the surface bias from
# -7.44degC to -0.02degC. `CSHE_MULT` makes that a real, properly-
# calibrated parameter instead of a hand-picked "8". See `Hydrodynamics/
# HeatExchange.jl`'s `compute_equilibrium_temperature!` docstring for the
# full mechanism/diagnostic writeup.
#
# WHY GRADIENT-FREE, NOT ENZYME (2026-09-23 decision, see CLAUDE.md):
# differentiating through the real hydrodynamic+transport pipeline with
# Enzyme.jl proved unstable on this Windows/Julia 1.11.3/Enzyme 0.13.204
# setup -- repeated process crashes and one NaN-gradient bug (found and
# fixed: `Solvers/Tridiagonal.jl`'s `thomas_solve!` used `similar()` for
# scratch arrays, a known Enzyme reverse-mode pitfall, now `zeros()`) did
# not fully resolve it. With only 4 target parameters, a gradient-free
# Bayesian approach (Gaussian-process surrogate + Expected Improvement,
# via Surrogates.jl) sidesteps the whole stability question and is well-
# suited to this problem size (few dimensions, expensive-but-not-too-
# expensive objective evaluations -- each one a real ~365-day DET run).
# Revisit Enzyme later (a different version, or a more Enzyme-mature
# platform like Linux/WSL) per user's explicit plan.
#
# OBJECTIVE: for a given (EXH2O, BETA, CBHE, TSED), run the real 2017 DET
# forced simulation (same physics/forcing as `run_det_2017.jl`, but kept
# in-process -- no CSV round-trip per evaluation, since Bayesian
# optimization needs tens of evaluations and file I/O for a ~14,600-row
# profile + re-parsing it each time would add real overhead) and compute
# RMSE against `Experiments/DET/obs/DET_temp.csv` (same matching logic as
# `compare_det_2017.jl`: nearest-day, nearest-elevation-below-surface).
# Returns RMSE (to MINIMIZE).
#
# BOUNDS: physically plausible ranges bracketing this port's own flagged
# placeholder defaults (EXH2O=0.35, BETA=0.6, CBHE=0.3, TSED=12.0) --
# EXH2O in [0.1,1.0] 1/m, BETA in [0.3,0.9], CBHE in [0.05,1.0] W/m^2/degC,
# TSED in [5,20] degC. Not independently re-derived from source (no real
# DET value exists to check against, that's exactly the gap being
# calibrated) -- wide enough to let the optimizer move meaningfully, not
# so wide the surrogate has to cover an implausible region.
# ==============================================================================

using W2J
using Surrogates

const CON = "K:/Git_repos/CE-QUAL-W2_JULIA/Experiments/DET/w2_con.csv"
const BASE = "K:/Git_repos/CE-QUAL-W2_JULIA/Experiments/DET"
const SEG = 33
const JDAY_START = 367.0
const JDAY_END = 731.0
const NSTEPS_CAP = 2_000_000

"""
    run_year_profile(exh2o, beta, cbhe, tsed, cshe_mult) -> (jday, el, t)

Runs the real 2017 DET forced simulation with the given heat-exchange
coefficients, returning the daily full-column profile at `SEG` as three
parallel vectors (jday, elevation_m, temp_C) -- the same data
`run_det_2017.jl`'s `profile_seg33.csv` carries, kept in memory instead
of round-tripping through a file.
"""
function run_year_profile(exh2o, beta, cbhe, tsed, cshe_mult)
    g, geom, tc = W2J.InputReader.read_control_file(CON; debug=false)
    W2J.InputReader.allocate_geometry!(g, geom)
    bthfn = joinpath(BASE, "inputs", "DET_bathy_v2.csv")
    W2J.BathymetryReader.read_bathymetry!(geom, g, bthfn, 1; debug=false)
    g, geom, net = W2J.init_geometry!(g, geom)
    W2J.compute_dlxrho!(g, geom)
    W2J.allocate_hydro_state!(g)
    W2J.allocate_transport_state!(g; dz_const=1e-3)
    W2J.compute_sf1x!(g, geom)
    geom.THETA = fill(0.55, g.NWB)
    geom.UPWIND = fill(true, g.NWB)
    geom.ULTIMATE = fill(false, g.NWB)
    geom.PLACE_QIN = fill(true, g.NWB)
    g.DIST_TRIBS[1] = true

    geom.EXH2O = fill(exh2o, g.NWB)
    geom.BETA = fill(beta, g.NWB)
    geom.CC_SW = fill(0.0, g.NWB)
    geom.CBHE = fill(cbhe, g.NWB)
    geom.TSED = fill(tsed, g.NWB)
    geom.TSEDF = fill(0.0, g.NWB)
    geom.CSHE_MULT = fill(cshe_mult, g.NWB)

    kt = g.KTWB[1]
    for i in 1:g.IMX
        kbi = g.KB[i]
        kbi < kt && continue
        for k in kt:kbi
            g.HYD[k, i, 4] = geom.T2I[1]
        end
    end
    g.T1 .= g.HYD[:, :, 4]

    boundary = W2J.BoundaryReader.load_boundary_conditions(CON, BASE, g)
    met = W2J.MetReader.load_met_conditions(CON, BASE, g)

    state = W2J.init_adaptive_timestep(tc)
    jday = JDAY_START
    last_profile_day = floor(Int, jday) - 1

    prof_jday = Float64[]; prof_el = Float64[]; prof_t = Float64[]
    function record!(jd)
        for k in kt:g.KB[SEG]
            push!(prof_jday, jd); push!(prof_el, geom.EL[k, SEG]); push!(prof_t, g.T1[k, SEG])
        end
    end
    record!(jday)

    step = 0
    while jday < JDAY_END && step < NSTEPS_CAP
        step += 1
        W2J.BoundaryReader.update_boundary_conditions!(g, boundary, jday)
        W2J.MetReader.update_met_conditions!(g, met, jday)
        accepted_dlt, _ = W2J.step_hydrodynamics_adaptive!(g, geom, net, tc, state, jday)
        W2J.apply_temperature_sources!(g, geom)
        W2J.compute_short_wave_radiation!(g, geom, jday)
        for jw in 1:g.NWB
            W2J.compute_equilibrium_temperature!(g, geom, jw)
        end
        W2J.apply_surface_heat_exchange!(g, geom)
        W2J.temperature_transport!(g, geom, accepted_dlt)
        g.HYD[:, :, 4] .= g.T1
        jday += accepted_dlt / 86400.0

        today = floor(Int, jday)
        if today > last_profile_day
            record!(jday)
            last_profile_day = today
        end
    end
    return (jday=prof_jday, el=prof_el, t=prof_t)
end

"""
    model_temp_at(profile, jday, target_el) -> Union{Float64,Nothing}

Same matching logic as `compare_det_2017.jl`'s function of the same name,
operating on in-memory profile vectors instead of a loaded CSV.
"""
function model_temp_at(profile, jday::Real, target_el::Real, days_sorted)
    isempty(days_sorted) && return nothing
    (jday < days_sorted[1] - 1 || jday > days_sorted[end] + 1) && return nothing
    nearest_day = days_sorted[argmin(abs.(days_sorted .- jday))]
    mask = profile.jday .== nearest_day
    els = profile.el[mask]
    ts = profile.t[mask]
    isempty(els) && return nothing
    idx = argmin(abs.(els .- target_el))
    return ts[idx]
end

"""
    rmse_against_observed(profile, obs) -> Float64

Same matching/scoring logic as `compare_det_2017.jl`'s `main`, returning
just the scalar RMSE (the objective Bayesian optimization minimizes).
"""
function rmse_against_observed(profile, obs)
    days = sort(unique(profile.jday))
    elws_by_day = Dict{Float64,Float64}()
    for d in days
        mask = profile.jday .== d
        elws_by_day[d] = maximum(profile.el[mask])
    end
    days_avail = collect(keys(elws_by_day))

    diffs = Float64[]
    for i in 1:length(obs.jday)
        jd = obs.jday[i]
        isempty(days_avail) && continue
        nearest_day = days_avail[argmin(abs.(days_avail .- jd))]
        abs(nearest_day - jd) > 2.0 && continue
        elws = elws_by_day[nearest_day]
        target_el = elws - obs.depth_m[i]
        mt = model_temp_at(profile, jd, target_el, days)
        mt === nothing && continue
        push!(diffs, mt - obs.temp_c[i])
    end
    isempty(diffs) && return Inf
    return sqrt(sum(diffs .^ 2) / length(diffs))
end

function main(; n_init=8, maxiters=10, num_new_samples=4)
    obs = W2J.ObservedDataReader.read_det_temp_profile(
        joinpath("K:/Git_repos/CE-QUAL-W2_JULIA", "Experiments/DET/obs/DET_temp.csv"); base_year=2016)
    println("Loaded ", length(obs.jday), " real observed points.")

    lb = (0.1, 0.3, 0.05, 5.0, 1.0)
    ub = (1.0, 0.9, 1.0, 20.0, 12.0)

    eval_count = Ref(0)
    function objective(x)
        exh2o, beta, cbhe, tsed, cshe_mult = x
        eval_count[] += 1
        t0 = time()
        profile = run_year_profile(exh2o, beta, cbhe, tsed, cshe_mult)
        rmse = rmse_against_observed(profile, obs)
        println("[eval $(eval_count[])] EXH2O=$(round(exh2o,digits=3)) BETA=$(round(beta,digits=3)) ",
                "CBHE=$(round(cbhe,digits=3)) TSED=$(round(tsed,digits=2)) CSHE_MULT=$(round(cshe_mult,digits=2))  ",
                "RMSE=$(round(rmse,digits=4))  (", round(time()-t0, digits=1), "s)")
        flush(stdout)
        return rmse
    end

    println("Generating $n_init initial samples...")
    xs = Surrogates.sample(n_init, lb, ub, SobolSample())
    ys = objective.(xs)

    println("Building Kriging surrogate and running EI optimization ($maxiters iters)...")
    krig = Kriging(xs, ys, lb, ub)
    surrogate_optimize!(objective, EI(), lb, ub, krig, SobolSample();
                         maxiters=maxiters, num_new_samples=num_new_samples)

    best_idx = argmin(krig.y)
    best_x = krig.x[best_idx]
    best_y = krig.y[best_idx]
    println()
    println("=== BEST FOUND after $(eval_count[]) evaluations ===")
    println("EXH2O=", best_x[1], "  BETA=", best_x[2], "  CBHE=", best_x[3], "  TSED=", best_x[4],
            "  CSHE_MULT=", best_x[5])
    println("RMSE=", best_y)
    return (best_x=best_x, best_y=best_y, all_x=krig.x, all_y=krig.y)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
