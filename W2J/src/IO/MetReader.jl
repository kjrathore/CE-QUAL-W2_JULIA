# ==============================================================================
# IO/MetReader.jl
#
# Reads meteorology (MET.csv), the Tier-1 gap flagged repeatedly since the
# TKE turbulence closure went in (2026-08-24): `Hydrodynamics/Turbulence.jl`'s
# `calculate_tke!` needs real wind (`WIND10`/`CZ`) for its surface production
# term, previously always 0. Same `$`-format convention as `IO/
# BoundaryReader.jl` (`$`-comment, `#`-comment, header, then data rows) --
# confirmed against real DET `inputs/MET.csv` (NASA POWER API export),
# columns `JDAY,TAIR.C,TDEW.C,WS.MS,WD.RAD,CLOUD,GHI.WM2`, matching the real
# read order `time-varying-data.f90:278` (`TAIRNX,TDEWNX,WINDNX,PHINX,
# CLOUDNX,SRONX`) 1:1.
#
# SCOPE -- narrower than the full real MET pipeline:
#   REAL: the file format/column order, linear interpolation (matching
#   `time-varying-data.f90:2779`'s `WIND(JW) = (1-RATIO)*WINDNX+RATIO*WINDO`
#   exactly -- same formula already validated for QIN/TIN), and `Hydrodynamics/
#   Turbulence.jl`'s `compute_wind_stress!` (`WIND10`/`CZ`, the two values
#   `calculate_tke!` actually consumes).
#
#   NOT YET USED (loaded into `g.TAIR`/`TDEW`/`CLOUD`/`SRO` for a future
#   session, but nothing reads them yet): real surface heat exchange
#   (`SHORT_WAVE_RADIATION`/`EQUILIBRIUM_TEMPERATURE`/`SURFACE_TERMS`,
#   `heat-exchange.f90`, still a stub) would consume these -- this port has
#   no evaporation/heat-flux physics yet, so `TAIR`/`TDEW`/`CLOUD`/`SRO`
#   are stored but inert for now.
#
#   REDUCED in `compute_wind_stress!` (see that function's docstring):
#   `WSC(I)` (per-segment wind-sheltering coefficient, Tier 1) defaults to
#   1.0; the log-law height rescale (`WINDH(JW)`/`Z0(JW)`, both Tier 1) is
#   skipped entirely, assuming the input wind speed is already a 10m
#   reference value (a reasonable assumption for the real DET data, a NASA
#   POWER API export -- POWER's WS10M product is defined at 10m); wind
#   DIRECTION (`PHI`/`PHI0`) is loaded but not consumed -- `calculate_tke!`'s
#   `USTAR` term only needs wind SPEED (`WIND10`/`CZ`), not direction; the
#   fetch-based correction (`FETCH_CALC(JW)`, Tier 1) is not applied.
# ==============================================================================

module MetReader

export MetSeries, read_met_series, interpolate_met, find_met_filename,
       load_met_conditions, update_met_conditions!

"""
    MetSeries

A loaded, in-memory meteorology time series -- `JDAY` plus the six real
columns (`TAIR`, `TDEW`, `WIND`, `PHI`, `CLOUD`, `SRO`), matching real
Fortran's `TAIRNX`/`TDEWNX`/`WINDNX`/`PHINX`/`CLOUDNX`/`SRONX` 1:1.
"""
struct MetSeries
    jday::Vector{Float64}
    tair::Vector{Float64}
    tdew::Vector{Float64}
    wind::Vector{Float64}
    phi::Vector{Float64}
    cloud::Vector{Float64}
    sro::Vector{Float64}
end

"""
    read_met_series(path) -> MetSeries

Reads a `\$`-format MET file (same convention as `IO/BoundaryReader.jl`'s
`read_boundary_series`). Rows with any missing value are skipped rather
than erroring -- same discipline as `read_boundary_series`, since real
sensor/derived-product data can have gaps.
"""
function read_met_series(path::AbstractString)
    lines = readlines(path)
    jday = Float64[]; tair = Float64[]; tdew = Float64[]; wind = Float64[]
    phi = Float64[]; cloud = Float64[]; sro = Float64[]
    for line in lines[4:end]
        isempty(strip(line)) && continue
        fields = split(line, ',')
        length(fields) < 7 && continue
        vals = strip.(fields[1:7])
        any(isempty, vals) && continue
        push!(jday, parse(Float64, vals[1]))
        push!(tair, parse(Float64, vals[2]))
        push!(tdew, parse(Float64, vals[3]))
        push!(wind, parse(Float64, vals[4]))
        push!(phi, parse(Float64, vals[5]))
        push!(cloud, parse(Float64, vals[6]))
        push!(sro, parse(Float64, vals[7]))
    end
    return MetSeries(jday, tair, tdew, wind, phi, cloud, sro)
end

"""
    interpolate_met(series, jday) -> (tair, tdew, wind, phi, cloud, sro)

Linear interpolation per column, matching `time-varying-data.f90:2779`'s
`WIND(JW) = (1-RATIO)*WINDNX+RATIO*WINDO` exactly (the same formula already
validated for QIN/TIN in `IO/BoundaryReader.jl`'s `interpolate_series` --
duplicated here rather than shared since it operates on 6 columns at once
off one shared `jday` index, not one column at a time).
"""
function interpolate_met(series::MetSeries, jday::Real)
    n = length(series.jday)
    n == 0 && error("interpolate_met: empty series")
    if jday <= series.jday[1]
        i = 1
    elseif jday >= series.jday[n]
        i = n
    else
        i = searchsortedlast(series.jday, jday)
        j1, j2 = series.jday[i], series.jday[i+1]
        ratio = (jday - j1) / (j2 - j1)
        return ((1 - ratio) * series.tair[i] + ratio * series.tair[i+1],
                (1 - ratio) * series.tdew[i] + ratio * series.tdew[i+1],
                (1 - ratio) * series.wind[i] + ratio * series.wind[i+1],
                (1 - ratio) * series.phi[i] + ratio * series.phi[i+1],
                (1 - ratio) * series.cloud[i] + ratio * series.cloud[i+1],
                (1 - ratio) * series.sro[i] + ratio * series.sro[i+1])
    end
    return (series.tair[i], series.tdew[i], series.wind[i], series.phi[i], series.cloud[i], series.sro[i])
end

"""
    find_met_filename(con_path, nwb) -> Vector{String}

Searches `w2_con.csv` for the per-waterbody filename block (a
"WB1,WB2,...,WB<n>" header row, e.g. `w2_con.csv:869` in the DET
experiment) followed by real `.csv` rows, and reads `METFN` at
`wb_row+2` -- row order confirmed against the real DET control file
(`wb_row+1`=BTHFN, `wb_row+2`=METFN, `wb_row+3`=EXTFN). SEARCH-based, not
positional, and requiring a `.csv` filename on the row below (not just the
header text) for the same reason as `IO/BoundaryReader.jl`'s
`find_boundary_filenames`: "WB1,WB2,..." is a generic per-waterbody header
reused by multiple unrelated blocks in the real file.
"""
function find_met_filename(con_path::AbstractString, nwb::Int)
    lines = readlines(con_path)
    is_real_file(s) = endswith(lowercase(strip(s)), ".csv")
    wb_row = findfirst(eachindex(lines)) do i
        fields = strip.(split(lines[i], ','))
        (length(fields) >= 1 && fields[1] == "WB1") &&
            i + 1 <= length(lines) && any(is_real_file, split(lines[i+1], ','))
    end
    wb_row === nothing && error("find_met_filename: no \"WB1,WB2,...\" row followed by real .csv filenames found in $con_path")
    metfn = strip.(split(lines[wb_row+2], ','))[1:nwb]
    return [is_real_file(f) ? f : "" for f in metfn]
end

"""
    load_met_conditions(con_path, base_dir, g) -> Dict{Int,MetSeries}

Loads the MET file for every waterbody with a real (non-"not used")
filename. Returns `Dict{Int,MetSeries}` keyed by waterbody number --
deliberately NOT stored inside `W2Global`, same "separate object threaded
through calls" pattern as `IO/BoundaryReader.jl`'s `load_boundary_
conditions`. Branch filenames use Windows-style backslashes
(`inputs\\MET.csv`) -- normalized via `joinpath` on the split components.
"""
function load_met_conditions(con_path::AbstractString, base_dir::AbstractString, g)
    fn = find_met_filename(con_path, g.NWB)
    met = Dict{Int,MetSeries}()
    for jw in 1:g.NWB
        f = fn[jw]
        isempty(f) && continue
        path = joinpath(base_dir, split(f, '\\')...)
        met[jw] = read_met_series(path)
    end
    return met
end

"""
    update_met_conditions!(g, met, jday)

Sets `g.TAIR[jw]`/`g.TDEW[jw]`/`g.WIND[jw]`/`g.PHI[jw]`/`g.CLOUD[jw]`/
`g.SRO[jw]` for every waterbody in `met` (from `load_met_conditions`) at
the given `jday`. Call once per timestep, before `Hydrodynamics/
Turbulence.jl`'s `compute_wind_stress!` (which reads `g.WIND`).
"""
function update_met_conditions!(g, met::Dict{Int,MetSeries}, jday::Real)
    for (jw, series) in met
        tair, tdew, wind, phi, cloud, sro = interpolate_met(series, jday)
        g.TAIR[jw] = tair; g.TDEW[jw] = tdew; g.WIND[jw] = wind
        g.PHI[jw] = phi; g.CLOUD[jw] = cloud; g.SRO[jw] = sro
    end
    return g
end

end # module MetReader
