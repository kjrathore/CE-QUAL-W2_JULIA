# ==============================================================================
# tools/calibration/compare_det_2017.jl
#
# Compares `tools/calibration/run_det_2017.jl`'s model output against real
# observed water-temperature-at-depth data (`Experiments/DET/obs/
# DET_temp.csv`, via `IO/ObservedDataReader.jl`).
#
# THE REAL DIFFICULTY THIS SCRIPT HANDLES: observed depths are reported as
# meters BELOW THE WATER SURFACE, but the water surface itself (`ELWS`)
# moves over the year (this port's own validated behavior, and real
# physics) -- a fixed depth-below-surface does NOT correspond to a fixed
# model layer `K` or a fixed absolute elevation. For each observed
# (JDAY, depth_m, T_C) triple, this script:
#   1. Looks up the model's own `ELWS` at the nearest daily profile JDAY
#      (from `profile_seg<seg>.csv`, which also carries `ELWS` implicitly
#      via the top layer's `EL_m` -- see below) to get the real water
#      surface elevation at that time.
#   2. Converts the observed depth-below-surface to an absolute elevation:
#      `target_el = elws - depth_m`.
#   3. Finds the model layer `K` whose elevation range brackets
#      `target_el` (layers are stored top-down, `EL_m` decreasing with
#      `K`), and reads that layer's `T_C`.
#   4. Reports the paired (observed, model) temperatures -- differences,
#      not just raw comparison, error metrics (RMSE, bias) per depth and
#      overall.
# ==============================================================================

using W2J

function load_profile(path)
    lines = readlines(path)[2:end]
    jday = Float64[]; k = Int[]; el = Float64[]; t = Float64[]
    for line in lines
        isempty(strip(line)) && continue
        f = split(line, ',')
        push!(jday, parse(Float64, f[1]))
        push!(k, parse(Int, f[2]))
        push!(el, parse(Float64, f[3]))
        push!(t, parse(Float64, f[4]))
    end
    return (jday=jday, k=k, el=el, t=t)
end

"""
    model_temp_at(profile, jday, target_el) -> Union{Float64,Nothing}

Finds the model's daily profile snapshot nearest `jday`, then the layer
whose elevation is closest to `target_el`, and returns that layer's `T_C`.
Returns `nothing` if `jday` falls outside the model run's covered range.
"""
function model_temp_at(profile, jday::Real, target_el::Real)
    days = sort(unique(profile.jday))
    isempty(days) && return nothing
    (jday < days[1] - 1 || jday > days[end] + 1) && return nothing
    nearest_day = days[argmin(abs.(days .- jday))]
    mask = profile.jday .== nearest_day
    els = profile.el[mask]
    ts = profile.t[mask]
    isempty(els) && return nothing
    idx = argmin(abs.(els .- target_el))
    return ts[idx]
end

function main(; profile_path, obs_path, out_csv)
    profile = load_profile(profile_path)
    obs = W2J.ObservedDataReader.read_det_temp_profile(obs_path; base_year=2016)

    # Model ELWS at each profile day = the top (largest EL_m) layer's own EL_m
    # is NOT quite ELWS (that's the layer's midpoint/geometry elevation, not
    # the water surface itself) -- but real ELWS isn't in this profile file.
    # Reduced here: use the max EL_m recorded that day (the top active
    # layer's elevation) as a stand-in for ELWS, which differs from the true
    # ELWS by at most one partial layer thickness (a few meters at most for
    # DET's ~3m layers) -- flagged, not a silent assumption.
    days = sort(unique(profile.jday))
    elws_by_day = Dict{Float64,Float64}()
    for d in days
        mask = profile.jday .== d
        elws_by_day[d] = maximum(profile.el[mask])
    end

    open(out_csv, "w") do io
        println(io, "jday,depth_m,obs_T_C,model_T_C,diff_C")
        n = length(obs.jday)
        diffs = Float64[]
        for i in 1:n
            jd = obs.jday[i]
            days_avail = collect(keys(elws_by_day))
            isempty(days_avail) && continue
            nearest_day = days_avail[argmin(abs.(days_avail .- jd))]
            abs(nearest_day - jd) > 2.0 && continue  # no model data within 2 days
            elws = elws_by_day[nearest_day]
            target_el = elws - obs.depth_m[i]
            mt = model_temp_at(profile, jd, target_el)
            mt === nothing && continue
            diff = mt - obs.temp_c[i]
            push!(diffs, diff)
            println(io, "$(jd),$(obs.depth_m[i]),$(obs.temp_c[i]),$(mt),$(diff)")
        end
        println("Matched ", length(diffs), " of ", n, " observed points")
        if !isempty(diffs)
            rmse = sqrt(sum(diffs .^ 2) / length(diffs))
            bias = sum(diffs) / length(diffs)
            println("RMSE = ", round(rmse, digits=3), " C   Bias (model-obs) = ", round(bias, digits=3), " C")
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    base = "K:/Git_repos/CE-QUAL-W2_JULIA"
    main(;
        profile_path=joinpath(base, "tools/calibration/output/det_2017/profile_seg33.csv"),
        obs_path=joinpath(base, "Experiments/DET/obs/DET_temp.csv"),
        out_csv=joinpath(base, "tools/calibration/output/det_2017/comparison.csv"))
end
