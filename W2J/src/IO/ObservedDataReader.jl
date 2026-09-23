# ==============================================================================
# IO/ObservedDataReader.jl
#
# Reads real observed reservoir data exported from USACE's NWD "dataquery"
# tool (www.nwd-wc.usace.army.mil/dd/common/dataquery/www/) -- the first
# calibration-target data this project has, per the user's explicit
# priority (Enzyme.jl-based parameter calibration, "Three Pillars" #2).
#
# FORMAT -- confirmed against the real downloaded file
# (`Experiments/DET/obs/DET_temp.csv`), NOT the `$`-format convention used
# elsewhere in this project (`IO/BoundaryReader.jl`/`IO/MetReader.jl`) --
# this is a raw USACE DSS export: `Date Time,<series1>,<series2>,...` header,
# then `DD-Mon-YYYY HH:MM,val1,val2,...` rows.
#
# REAL QUIRK, found via direct inspection, not guessed: the header's first
# depth-series label ("D0.5ft") is itself comma-broken by the export tool
# into TWO header fields ("DET_S1-D0" then "5ft.Temp-Water...") -- a
# decimal-comma artifact in the LABEL text only, not in the numeric data
# rows (which use a period decimal point throughout, confirmed: row 2's
# field count is 9 = 1 date + 8 values, one less than the header's 10
# comma-split fields). This reader does NOT trust the header for column
# order -- the real depth order (confirmed against the file directly) is
# hardcoded as `[0.5, 100, 10, 120, 20, 40, 60, 80]` feet, matching the data
# columns' actual left-to-right order.
#
# UNITS: real USACE export is Fahrenheit and feet-below-surface --
# converted to Celsius and meters-below-surface here, matching this port's
# existing SI-everywhere convention (`Hydrodynamics/Density.jl` etc. all
# take Celsius).
# ==============================================================================

module ObservedDataReader

export ObservedProfile, read_det_temp_profile, jday_from_datetime

const FT_TO_M = 0.3048
const DET_TEMP_DEPTHS_FT = [0.5, 100.0, 10.0, 120.0, 20.0, 40.0, 60.0, 80.0]

const MONTH_ABBR = Dict("Jan" => 1, "Feb" => 2, "Mar" => 3, "Apr" => 4, "May" => 5, "Jun" => 6,
                         "Jul" => 7, "Aug" => 8, "Sep" => 9, "Oct" => 10, "Nov" => 11, "Dec" => 12)
const DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

"""
    is_leap_year(y) -> Bool

Standard Gregorian leap-year rule (divisible by 4, not by 100 unless also
by 400) -- needed to convert real calendar dates into this port's `JDAY`
convention (`1.0` = `YEAR`-Jan-1 00:00, matching `IO/InputReader.jl`'s
`TMSTRT`/`YEAR`).
"""
is_leap_year(y::Int) = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0

"""
    jday_from_datetime(dt_str, base_year) -> Float64

Parses a USACE-export date-time string (`"DD-Mon-YYYY HH:MM"`, e.g.
`"09-Feb-2016 22:00"`) into this port's `JDAY` convention relative to
`base_year` (the real Fortran/`IO/InputReader.jl` convention: `JDAY=1.0`
is `base_year`-Jan-1 at 00:00, so a date one full year later than
`base_year` lands past `JDAY=365`/`366` depending on `base_year`'s own
leap-year status -- e.g. `base_year=2016` (a leap year, 366 days) puts
2017-Jan-1 at `JDAY=367.0`).
"""
function jday_from_datetime(dt_str::AbstractString, base_year::Int)
    parts = split(strip(dt_str))
    length(parts) != 2 && error("jday_from_datetime: expected \"DD-Mon-YYYY HH:MM\", got \"$dt_str\"")
    date_part, time_part = parts
    dfields = split(date_part, '-')
    length(dfields) != 3 && error("jday_from_datetime: bad date \"$date_part\"")
    day = parse(Int, dfields[1])
    month = MONTH_ABBR[dfields[2]]
    year = parse(Int, dfields[3])
    tfields = split(time_part, ':')
    hour = parse(Int, tfields[1])
    minute = length(tfields) > 1 ? parse(Int, tfields[2]) : 0

    days_before_year = 0
    for y in base_year:(year-1)
        days_before_year += is_leap_year(y) ? 366 : 365
    end
    days_before_month = sum(DAYS_IN_MONTH[1:(month-1)])
    (month > 2 && is_leap_year(year)) && (days_before_month += 1)
    day_of_year = days_before_month + day  # 1-based day within `year`

    return Float64(days_before_year + day_of_year) + hour / 24.0 + minute / 1440.0
end

"""
    ObservedProfile

Real observed water-temperature-at-depth data, one row per (timestamp,
depth) pair -- a "long" table, not one column per depth, so downstream
comparison code doesn't need to know the depth count/order up front.
`depth_m` is meters BELOW the water surface at the time of the reading
(not an absolute elevation -- the real water surface moves over a year,
see `IO/ObservedDataReader.jl`'s module docstring and the comparison
script that converts this to an absolute elevation using the model's own
`ELWS` at the matching `jday`).
"""
struct ObservedProfile
    jday::Vector{Float64}
    depth_m::Vector{Float64}
    temp_c::Vector{Float64}
end

"""
    read_det_temp_profile(path; base_year) -> ObservedProfile

Reads a USACE `DET_temp.csv`-format export (see module docstring for the
real format/quirks). Rows with any missing/blank value are skipped (same
discipline as `IO/BoundaryReader.jl`'s `read_boundary_series` -- real
sensor data has gaps, confirmed for this exact file: 2016 has a real
multi-month gap, user-reported and independently confirmed via direct
inspection). Returns one `ObservedProfile` entry per (timestamp, depth)
pair across all 8 depths -- `length(result.jday) == 8 * (number of
complete data rows)`.
"""
function read_det_temp_profile(path::AbstractString; base_year::Int)
    lines = readlines(path)
    jday = Float64[]; depth_m = Float64[]; temp_c = Float64[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, ',')
        length(fields) < 9 && continue
        vals = strip.(fields)
        any(isempty, vals[2:9]) && continue
        jd = jday_from_datetime(vals[1], base_year)
        for k in 1:8
            temp_f = parse(Float64, vals[k+1])
            push!(jday, jd)
            push!(depth_m, DET_TEMP_DEPTHS_FT[k] * FT_TO_M)
            push!(temp_c, (temp_f - 32.0) * 5.0 / 9.0)
        end
    end
    return ObservedProfile(jday, depth_m, temp_c)
end

end # module ObservedDataReader
