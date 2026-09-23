# ==============================================================================
# tools/calibration/run_det_2017.jl
#
# First calibration-prep run: simulates real DET (Detroit Lake, OR) for
# calendar year 2017 with real boundary/MET forcing, and captures a daily
# full-water-column temperature profile at segment 33 (the dam segment,
# user-specified 2026-09-22 -- matches the real Fortran reference run's own
# choice, `Experiments/DET/outputs/tsr_1_seg33.csv`) for comparison against
# real observed thermistor-string data (`Experiments/DET/obs/DET_temp.csv`,
# downloaded from USACE's NWD "dataquery" tool).
#
# WHY 2017, NOT 2016: the observed data has a real, user-confirmed multi-
# month gap in 2016 (Jun-Sep); 2017 is clean aside from a Jan-~Feb 24 gap
# (confirmed via IO/ObservedDataReader.jl's own validation: first real 2017
# reading is at JDAY=421.375, not JDAY=367). The simulation still starts at
# JDAY=367 (2017-Jan-01) with a uniform-temperature initial condition
# (T2I=5.0, DET's own real IC value, reused here as a reasonable winter
# approximation since no real profile exists for that exact date) -- the
# Jan-Feb window before real data begins acts as an in-model spin-up, not a
# claim that the IC itself is validated for 2017-Jan-1 specifically.
#
# NOT a full year from JDAY=1 (2016-Jan-1): a continuous 2-year run was
# considered but rejected on cost grounds (roughly 2x the per-step cost of
# TKE's own TRIDIAG solve, real MET/boundary forcing, and adaptive-timestep
# retries, over an already-long single-year run) -- starting fresh at
# 2017-Jan-1 is the deliberate, reduced-cost choice for this FIRST
# calibration-prep run, not an oversight.
#
# OUTPUT: two files under `out_dir` --
#   - `tsr_seg<seg>.csv` (existing OutputWriter.jl format, JDAY/DLT/ELWS/T2/U
#     at the segment's own top active layer, one row per timestep)
#   - `profile_seg<seg>.csv` (NEW, this script only -- one row per (day, K)
#     pair: JDAY, K, EL_m (real layer elevation), T_C -- captured once per
#     calendar day, not every raw timestep, to keep file size manageable:
#     365 days x ~40 layers ~= 14,600 rows, vs. millions for a per-timestep
#     dump).
# ==============================================================================

using W2J

function main(; nsteps_limit::Union{Int,Nothing}=nothing, out_dir::AbstractString)
    con = "K:/Git_repos/CE-QUAL-W2_JULIA/Experiments/DET/w2_con.csv"
    base = "K:/Git_repos/CE-QUAL-W2_JULIA/Experiments/DET"
    seg = 33

    g, geom, tc = W2J.InputReader.read_control_file(con; debug=false)
    W2J.InputReader.allocate_geometry!(g, geom)
    bthfn = joinpath(base, "inputs", "DET_bathy_v2.csv")
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

    kt = g.KTWB[1]
    for i in 1:g.IMX
        kbi = g.KB[i]
        kbi < kt && continue
        for k in kt:kbi
            g.HYD[k, i, 4] = geom.T2I[1]
        end
    end
    g.T1 .= g.HYD[:, :, 4]

    boundary = W2J.BoundaryReader.load_boundary_conditions(con, base, g)
    met = W2J.MetReader.load_met_conditions(con, base, g)

    jday_start = 367.0   # 2017-Jan-01, relative to YEAR=2016 (TMSTRT=1 convention)
    jday_end = 731.0     # 2017-Dec-31
    dlt = tc.DLTMAX[1]
    # Loop by jday, not a fixed step count -- a fixed budget assuming average
    # dlt=DLTMAX undercounts once compute_curmax's QOT-aware CFL term
    # (2026-09-23, see AdaptiveTimestep.jl) legitimately shrinks dlt near the
    # withdrawal cell more often than DLTMAX would suggest. nsteps_cap is a
    # generous safety bound only, not the real stopping condition.
    nsteps_cap = nsteps_limit !== nothing ? nsteps_limit : 2_000_000

    isdir(out_dir) || mkpath(out_dir)
    tsr_writer = W2J.OutputWriter.open_tsr_files(out_dir, "tsr", [seg]; include_temp=true)
    profile_path = joinpath(out_dir, "profile_seg$(seg).csv")
    profile_io = open(profile_path, "w")
    println(profile_io, "JDAY,K,EL_m,T_C")

    state = W2J.init_adaptive_timestep(tc)
    jday = jday_start
    last_profile_day = floor(Int, jday) - 1

    function write_profile!(jd)
        for k in kt:g.KB[seg]
            println(profile_io, "$(jd),$(k),$(geom.EL[k, seg]),$(g.T1[k, seg])")
        end
    end

    W2J.OutputWriter.write_tsr_row!(tsr_writer, g, geom, jday, dlt)
    write_profile!(jday)

    step = 0
    while jday < jday_end && step < nsteps_cap
        step += 1
        W2J.BoundaryReader.update_boundary_conditions!(g, boundary, jday)
        W2J.MetReader.update_met_conditions!(g, met, jday)
        accepted_dlt, _ = W2J.step_hydrodynamics_adaptive!(g, geom, net, tc, state, jday)
        W2J.apply_temperature_sources!(g, geom)
        W2J.temperature_transport!(g, geom, accepted_dlt)
        g.HYD[:, :, 4] .= g.T1
        jday += accepted_dlt / 86400.0
        W2J.OutputWriter.write_tsr_row!(tsr_writer, g, geom, jday, accepted_dlt)

        today = floor(Int, jday)
        if today > last_profile_day
            write_profile!(jday)
            last_profile_day = today
        end

        if step % 5000 == 0
            println("step=$step / cap=$nsteps_cap  jday=$(round(jday, digits=2))  ELWS[seg]=$(round(geom.ELWS[seg], digits=3))  T1[kt,seg]=$(round(g.T1[kt, seg], digits=3))")
            flush(stdout)
        end
    end

    W2J.OutputWriter.close_tsr_files!(tsr_writer)
    close(profile_io)
    println("Done. jday_final=$jday  wrote $(out_dir)")
    return (g, geom, jday)
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Experiments/ is read-only (user instruction) -- model output goes to
    # tools/calibration/output/, NOT under Experiments/, even though the
    # real observed data it gets compared against lives there.
    out_dir = length(ARGS) >= 1 ? ARGS[1] : "K:/Git_repos/CE-QUAL-W2_JULIA/tools/calibration/output/det_2017"
    nsteps_limit = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
    main(; out_dir=out_dir, nsteps_limit=nsteps_limit)
end
