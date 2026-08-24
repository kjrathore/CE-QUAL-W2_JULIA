# ==============================================================================
# IO/BoundaryReader.jl
#
# First-cut Tier-1 boundary-condition IO: external upstream inflow quantity
# and temperature (QIN/TIN) and external downstream outflow (QOT), per
# branch. Traced from time-varying-data.f90 (READ_INPUT_DATA /
# INTERPOLATE_INPUTS entries), w2_4_win.f90 (where QIN(JB)/QOUT(K,JB)
# actually couple into the free-surface solve and the boundary velocity),
# and withdrawal.f90 (where QOUT(K,JB) really comes from -- DOWNSTREAM_
# WITHDRAWAL's selective-withdrawal algorithm, not ported, see below) --
# not guessed -- see Hydrodynamics/FreeSurface.jl's module docstring for
# the exact traced insertion points.
#
# SCOPE -- REDUCED, BY DELIBERATE CHOICE (confirmed with user, 2026-08-17),
# same discipline as every other module in this project:
#
#   REAL (traced from source, not guessed):
#   - The boundary-filename block in w2_con.csv: a "BR1,BR2,...,BR<n>" label
#     row followed by one row per boundary type (QINFN/TINFN/CINFN/QOTFN/
#     QDTFN/TDTFN/CDTFN/PRECIPFN), one filename per branch, "<name> - not
#     used" for branches without that boundary type. This file reads the
#     QINFN/TINFN/QOTFN rows -- the rest (QDT, TDT, CIN, PRECIP) are
#     located and skipped, not parsed, see `find_boundary_filenames`.
#   - The time-series file format ($-comment, #-comment, header, then
#     `JDAY,VALUE` rows) -- the same $-prefixed convention already used by
#     BathymetryReader.jl, confirmed against real DET/JPP/OOLO boundary
#     files (QIN_*.csv, TIN_*.csv), not assumed from Detroit's OWN boundary
#     files (DetroitReservoir/InputFiles/*.npt) -- those are the LEGACY
#     fixed-width format this project has always excluded (see
#     IO/InputReader.jl's module docstring, "CSV format only"). DET's
#     boundary files are a genuinely better fit for this project's existing
#     format support than Detroit's own.
#   - The interpolation formula (`interpolate_series`): matches
#     time-varying-data.f90:2957's `QRATIO = (NXQIN1-JDAY)/(NXQIN1-NXQIN2);
#     QIND = (1-QRATIO)*QINNX + QRATIO*QINO` exactly -- standard piecewise-
#     linear interpolation between the two bracketing (JDAY, value) points.
#
#   SIMPLIFIED, BY DELIBERATE DESIGN (not a missing feature, a genuine
#   redesign -- the real Fortran behavior doesn't apply to this port's use
#   case):
#   - The whole time series is read into memory upfront and interpolated by
#     direct lookup (`searchsortedlast`), rather than replicating Fortran's
#     incremental single-record-buffered read (NXQIN1/NXQIN2 double
#     buffering) with live RESSIM-coupling file polling (`SLEEPQQ`, see
#     CLAUDE.md's "Key Fortran analysis facts" -- 36 SLEEPQQ calls in this
#     file). W2J targets standalone historical-simulation runs, not live
#     coupling to a running reservoir-operations model, so there is no
#     "wait for more data to arrive" case to handle -- the whole file
#     already exists on disk. Mathematically identical result for the
#     historical-run case; genuinely different (and correctly so) for the
#     live-coupling case this project doesn't target.
#   - `INTERP_INFLOW(JB)` (a per-branch Tier-1 control-file flag gating
#     whether QIN is interpolated or held at a step function) is NOT read
#     -- this file always interpolates. Unlike ULTIMATE/UPWIND (a
#     structurally different numerical scheme) or THETA (changes the
#     implicit/explicit blend), this only affects the smoothness of a
#     forcing time series, not which physics regime is active -- a lower-
#     risk simplification, but still flagged here rather than silently
#     assumed.
#   - Out-of-range JDAY (before the first record or after the last) clamps
#     to the nearest endpoint value rather than replicating Fortran's
#     TMEND-proximity RESSIM-polling logic (time-varying-data.f90:2098-
#     2126) -- reasonable for a bounded historical run where the control
#     file's own TMEND should stay within the boundary file's date range.
#
#   - `QOUT(K,JB)` (real: KMX x NBR, one value per withdrawal-active layer,
#     selected by withdrawal.f90's DOWNSTREAM_WITHDRAWAL entry -- a
#     selective-withdrawal algorithm that finds the layer(s) matching the
#     outlet elevation AND the water column's density stratification, not
#     ported) is REDUCED here to `Core/State.jl`'s `g.QOT[jb]`, a per-branch
#     SCALAR total (confirmed with user, 2026-08-17). Multiple real outlet
#     structures at one branch (e.g. DET's `QOT.csv` has separate
#     `POWER.cms`/`SPILLWAY.cms` columns) are summed into that one total by
#     `read_boundary_series_summed` rather than tracked/placed separately.
#     Unlike QIN (which needs a real boundary VELOCITY set, since IU-1 is a
#     pad cell outside the domain with no other equation for it), QOT only
#     needs the continuity coupling (`w2_4_win.f90:914-921`'s `+QOUT(K,JB)`
#     term in the D(ID) assembly) -- confirmed by reading withdrawal.f90:230
#     (`WHERE (QOUT(:,JB)==0.0) U(:,ID)=0.0`): U(ID) for a real interior-
#     adjacent segment is computed by the normal momentum equation, not
#     prescribed, so no `apply_outflow_boundary!`-equivalent function exists
#     -- the `d[id] += g.QOT[jb]` term in `solve_branch_free_surface!` is
#     the whole of it.
#   - `QDTR(JB)` (distributed tributary inflow -- a single per-branch value
#     spread laterally along the WHOLE branch, e.g. runoff entering all
#     along a reservoir arm rather than at one point) is traced to
#     hydroinout.F90:1322-1331: gated on `DIST_TRIBS(JB)` (a Tier-1 flag,
#     required explicit here too -- `Core/State.jl`'s `DIST_TRIBS`),
#     distributed across every segment in the branch proportional to that
#     segment's share of the branch's TOTAL top-layer surface area
#     (`AKBR = sum(BI(KT,I)*DLX(I))`), and added directly into `QSS` --
#     the SAME source/sink array `solve_branch_free_surface!` already
#     reads (`-QSS` in its D(I) assembly). This needed NO new coupling
#     point in `FreeSurface.jl`, unlike QIN/QOT -- just correct population,
#     via `Hydrodynamics/FreeSurface.jl`'s `distribute_tributary!`, ported
#     faithfully (real area-weighting formula, not simplified).
#
#   - `TDTR(JB)` (distributed tributary temperature) is the temperature
#     counterpart of QDTR, gated on the same `DIST_TRIBS(JB)` flag (real:
#     temperature.F90:416-427, `Hydrodynamics/Transport.jl`'s
#     `apply_temperature_sources!` -- NOT coupled into the free-surface
#     solve like QIN/QOT/QDTR, since temperature only affects the transport
#     pass). `Core/State.jl`'s `QDT` (the per-segment share `distribute_
#     tributary!` already computes) is reused for this rather than
#     recomputed -- matching the real Fortran's own reuse.
#
#   NOT PORTED YET (Tier 1, future work):
#   - CIN/PRECIP/MET -- constituent loading, precipitation, meteorology.
#   - The density-driven plunge-point inflow layer placement (`PLACE_QIN`,
#     w2_4_win.f90:1210-1270) -- confirmed with user (2026-08-17) as a
#     REDUCED-PHYSICS first cut: inflow velocity is placed entirely at the
#     top active layer KT (`Hydrodynamics/FreeSurface.jl`'s
#     `apply_inflow_boundary!`), not distributed to the layer(s) matching
#     the inflow's density. Physically wrong for cold/dense inflows
#     (they'd actually plunge below the surface), but a real, flagged
#     simplification, not a silent one -- port PLACE_QIN before trusting
#     stratification results from a real inflow-forced run.
#   - The real selective-withdrawal layer selection for QOT (see above) --
#     `g.QOT[jb]` is withdrawn from the bottom active layer only
#     (`solve_branch_free_surface!`'s `d[id] += g.QOT[jb]` at K=KB(ID),
#     confirmed with user, 2026-08-17), not the real density/elevation-
#     matched zone.
#   - INTERNAL_FLOW/DAM_INFLOW/HEAD_FLOW branches -- only genuine external
#     UP_FLOW/DN_FLOW branches with a real (non-"not used") filename are
#     loaded.
# ==============================================================================

module BoundaryReader

export BoundarySeries, read_boundary_series, interpolate_series,
       find_boundary_filenames, load_boundary_conditions, update_boundary_conditions!

"""
    BoundarySeries

A loaded, in-memory `(JDAY, value)` time series, read once upfront (see
module docstring for why this differs from Fortran's incremental
buffered-read pattern).
"""
struct BoundarySeries
    jday::Vector{Float64}
    value::Vector{Float64}
end

"""
    read_boundary_series(path) -> BoundarySeries

Reads a `\$`-format boundary time-series file (`\$`-comment line, `#`-comment
line, header line, then `JDAY,VALUE` rows) -- confirmed against real
`QIN_*.csv`/`TIN_*.csv` files from the DET/JPP/OOLO experiments. Only the
first two columns are read (`JDAY` and one value column) -- multi-column
files (e.g. `QOT.csv`'s `POWER.cms,SPILLWAY.cms`) are not yet supported,
see module docstring.

Rows with a missing value (`JDAY,` with nothing after the comma) are
skipped rather than erroring -- found via real data, not a hypothetical:
`QIN_14178000.csv`/`QIN_14179000.csv` (real USGS gauge records) have
genuine sensor gaps (e.g. `82,` with no reading). Skipping means
`interpolate_series` bridges the gap linearly using the nearest valid
points on either side -- a reasonable, flagged simplification, since how
the real Fortran source handles a missing record in this file format
hasn't been traced.
"""
function read_boundary_series(path::AbstractString)
    lines = readlines(path)
    jday = Float64[]
    value = Float64[]
    for line in lines[4:end]  # skip $-comment, #-comment, header
        isempty(strip(line)) && continue
        fields = split(line, ',')
        length(fields) < 2 && continue
        isempty(strip(fields[1])) && continue
        isempty(strip(fields[2])) && continue
        push!(jday, parse(Float64, strip(fields[1])))
        push!(value, parse(Float64, strip(fields[2])))
    end
    return BoundarySeries(jday, value)
end

"""
    read_boundary_series_summed(path) -> BoundarySeries

Same `\$`-format as `read_boundary_series`, but for files with MULTIPLE value
columns after `JDAY` (e.g. `QOT.csv`'s `JDAY,POWER.cms,SPILLWAY.cms` -- two
separate real outlet structures at one branch). Sums every value column
into one total per row -- the reduced-physics "lump structures into one
total" scope confirmed with user (2026-08-17), see module docstring. A
blank individual column contributes 0 to that row's sum (not skipped
entirely, unlike a wholly-blank row) -- reasonable for "this gate had no
flow that day," but flagged since it hasn't been confirmed against a real
example of a partially-blank multi-column row. A row is skipped only if
EVERY value column is blank (matching `read_boundary_series`'s single-
column gap handling).
"""
function read_boundary_series_summed(path::AbstractString)
    lines = readlines(path)
    jday = Float64[]
    value = Float64[]
    for line in lines[4:end]
        isempty(strip(line)) && continue
        fields = split(line, ',')
        length(fields) < 2 && continue
        isempty(strip(fields[1])) && continue
        cols = [strip(f) for f in fields[2:end] if !isempty(strip(f))]
        isempty(cols) && continue
        push!(jday, parse(Float64, strip(fields[1])))
        push!(value, sum(parse(Float64, c) for c in cols))
    end
    return BoundarySeries(jday, value)
end

"""
    interpolate_series(series, jday) -> Float64

time-varying-data.f90:2957's `QRATIO`/`QIND` formula, evaluated by direct
lookup against the whole in-memory series instead of Fortran's incremental
double-buffer (see module docstring). Clamps to the nearest endpoint value
outside the series' own JDAY range.
"""
function interpolate_series(series::BoundarySeries, jday::Real)
    n = length(series.jday)
    n == 0 && error("interpolate_series: empty series")
    jday <= series.jday[1] && return series.value[1]
    jday >= series.jday[n] && return series.value[n]
    i = searchsortedlast(series.jday, jday)
    j1, j2 = series.jday[i], series.jday[i+1]
    v1, v2 = series.value[i], series.value[i+1]
    ratio = (jday - j1) / (j2 - j1)
    return (1.0 - ratio) * v1 + ratio * v2
end

"""
    find_boundary_filenames(con_path, nbr) -> (qinfn, tinfn, qotfn, qdtfn, tdtfn)

Searches `w2_con.csv` for the "BR1,BR2,...,BR<n>" label row (the header of
the per-branch boundary-filename block, e.g. `w2_con.csv:883` in the DET
experiment) and reads the rows below it: QINFN (`br_row+1`), TINFN
(`br_row+2`), CINFN (`br_row+3`, skipped -- constituent loading not ported),
QOTFN (`br_row+4`), QDTFN (`br_row+5`), TDTFN (`br_row+6`) -- row order
confirmed against the real DET control file (`w2_con.csv:882-892`), not
assumed. SEARCH-based, not positional, deliberately: this block sits ~800
rows past where `IO/InputReader.jl`'s Phase A parsing stops (structures,
withdrawals, the full constituent block, output control, and kinetics rates
all sit between), and there is no value in requiring all of that to be
understood first just to reach one filename block. A row is a "not used"
placeholder (not a real file) if it doesn't end in `.csv` after trimming.
"""
function find_boundary_filenames(con_path::AbstractString, nbr::Int)
    lines = readlines(con_path)
    is_real_file(s) = endswith(lowercase(strip(s)), ".csv")

    # "BR1,BR2,..." is a generic per-branch header reused by several
    # unrelated blocks in w2_con.csv (confirmed against real DET data: 5
    # separate "BR1,BR2,BR3,BR4" rows exist -- branch US values, ON/OFF
    # flags, etc., not just this filename block). A bare header match alone
    # picked the WRONG one (the first, at row 47, branch US values) --
    # found via this function actually being run against real data and
    # returning empty filenames, not a hypothetical. The real discriminator
    # is the row immediately below: only the genuine filename block's data
    # row contains ".csv" -- so require that too, not just the header text.
    br_row = findfirst(eachindex(lines)) do i
        fields = strip.(split(lines[i], ','))
        (length(fields) >= 1 && fields[1] == "BR1") &&
            i + 1 <= length(lines) && any(is_real_file, split(lines[i+1], ','))
    end
    br_row === nothing && error("find_boundary_filenames: no \"BR1,BR2,...\" row followed by real .csv filenames found in $con_path")

    qinfn = strip.(split(lines[br_row+1], ','))[1:nbr]
    tinfn = strip.(split(lines[br_row+2], ','))[1:nbr]
    qotfn = strip.(split(lines[br_row+4], ','))[1:nbr]
    qdtfn = strip.(split(lines[br_row+5], ','))[1:nbr]
    tdtfn = strip.(split(lines[br_row+6], ','))[1:nbr]
    return (qinfn = [is_real_file(f) ? f : "" for f in qinfn],
            tinfn = [is_real_file(f) ? f : "" for f in tinfn],
            qotfn = [is_real_file(f) ? f : "" for f in qotfn],
            qdtfn = [is_real_file(f) ? f : "" for f in qdtfn],
            tdtfn = [is_real_file(f) ? f : "" for f in tdtfn])
end

"""
    load_boundary_conditions(con_path, base_dir, g) -> (inflow=Dict, outflow=Dict, dist_trib=Dict, dist_trib_temp=Dict)

Loads external upstream inflow (QIN/TIN, gated on `g.UP_FLOW[jb]`), external
downstream outflow (QOT, gated on `g.DN_FLOW[jb]`, summed across however
many real outlet-structure columns that branch's `QOTFN` file has -- see
`read_boundary_series_summed`), distributed tributary inflow (QDTR, gated
on `g.DIST_TRIBS[jb]` -- required explicit, see `Core/State.jl`'s QDTR
docstring), and distributed tributary temperature (TDTR, gated the same
way -- real Fortran couples it via `Hydrodynamics/Transport.jl`'s
`apply_temperature_sources!`, not the free-surface solve) for every branch
with a real (non-"not used") filename. Returns four dicts,
`inflow::Dict{Int,@NamedTuple{qin,tin}}`, `outflow::Dict{Int,BoundarySeries}`,
`dist_trib::Dict{Int,BoundarySeries}`, `dist_trib_temp::Dict{Int,
BoundarySeries}`, keyed by branch number -- deliberately NOT stored inside
`W2Global` (no Fortran equivalent holds a loaded-series cache; this follows
the same "separate object threaded through calls" pattern as `Core/Grid.jl`'s
`BranchNetwork`, not a new global). Branch filenames use Windows-style
backslashes (`inputs\\QIN_....csv`) in the real control files --
normalized to the platform path separator via `joinpath` on the split
components.

Does NOT filter on `INTERNAL_FLOW`/`DAM_INFLOW` (Tier 1, not read by
`InputReader.jl` yet) -- caller must not apply these to a branch where
either is true, since `w2_4_win.f90:950-955` computes `QIN(JB)` from
upstream velocity in that case instead of reading a file at all.
`g.DIST_TRIBS` must be explicitly sized to `NBR` and set by the caller
before this runs (a `BoundsError` on an empty `Bool[]` is the intended
failure mode for "forgot to set it", not a silent skip).
"""
function load_boundary_conditions(con_path::AbstractString, base_dir::AbstractString, g)
    fn = find_boundary_filenames(con_path, g.NBR)
    inflow = Dict{Int,@NamedTuple{qin::BoundarySeries, tin::BoundarySeries}}()
    outflow = Dict{Int,BoundarySeries}()
    dist_trib = Dict{Int,BoundarySeries}()
    dist_trib_temp = Dict{Int,BoundarySeries}()
    for jb in 1:g.NBR
        if g.UP_FLOW[jb]
            qf, tf = fn.qinfn[jb], fn.tinfn[jb]
            if !isempty(qf) && !isempty(tf)
                qpath = joinpath(base_dir, split(qf, '\\')...)
                tpath = joinpath(base_dir, split(tf, '\\')...)
                inflow[jb] = (qin = read_boundary_series(qpath), tin = read_boundary_series(tpath))
            end
        end
        if g.DN_FLOW[jb]
            of = fn.qotfn[jb]
            if !isempty(of)
                opath = joinpath(base_dir, split(of, '\\')...)
                outflow[jb] = read_boundary_series_summed(opath)
            end
        end
        if g.DIST_TRIBS[jb]
            df = fn.qdtfn[jb]
            if !isempty(df)
                dpath = joinpath(base_dir, split(df, '\\')...)
                dist_trib[jb] = read_boundary_series(dpath)
            end
            tdf = fn.tdtfn[jb]
            if !isempty(tdf)
                tdpath = joinpath(base_dir, split(tdf, '\\')...)
                dist_trib_temp[jb] = read_boundary_series(tdpath)
            end
        end
    end
    return (inflow = inflow, outflow = outflow, dist_trib = dist_trib, dist_trib_temp = dist_trib_temp)
end

"""
    update_boundary_conditions!(g, boundary, jday)

Sets `g.QIN[jb]`/`g.QIND[jb]`/`g.TIN[jb]`/`g.TIND[jb]` for every branch in
`boundary.inflow`, `g.QOT[jb]` for every branch in `boundary.outflow`, and
`g.QDTR[jb]`/`g.TDTR[jb]` for every branch in `boundary.dist_trib`/
`boundary.dist_trib_temp` (from `load_boundary_conditions`) at the given
`jday`. Call once per timestep, before `Hydrodynamics/FreeSurface.jl`'s
`solve_free_surface!` (and before `distribute_tributary!`, which actually
applies `QDTR` into `QSS`, and before `Hydrodynamics/Transport.jl`'s
`apply_temperature_sources!`, which reads `TDTR`). `QIN==QIND`/`TIN==TIND`
here (see `Core/State.jl`'s QIN/QIND docstring for why -- the INTERNAL_FLOW/
DAM_INFLOW cases that would make them differ aren't ported).
"""
function update_boundary_conditions!(g, boundary, jday::Real)
    for (jb, series) in boundary.inflow
        q = interpolate_series(series.qin, jday)
        t = interpolate_series(series.tin, jday)
        g.QIN[jb] = q; g.QIND[jb] = q
        g.TIN[jb] = t; g.TIND[jb] = t
    end
    for (jb, series) in boundary.outflow
        g.QOT[jb] = interpolate_series(series, jday)
    end
    for (jb, series) in boundary.dist_trib
        g.QDTR[jb] = interpolate_series(series, jday)
    end
    for (jb, series) in boundary.dist_trib_temp
        g.TDTR[jb] = interpolate_series(series, jday)
    end
    return g
end

end # module BoundaryReader
