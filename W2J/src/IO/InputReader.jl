# ==============================================================================
# IO/InputReader.jl
#
# Reads the CE-QUAL-W2 control file (w2_con.csv, the .csv branch ONLY -- the
# legacy .npt fixed-format branch in input.F90 is not ported, since every
# control file we've been given is .csv) into the W2Core structs.
#
# SCOPE OF THIS FILE (Phase A -- "base module stack" for geometry):
#   1. Header dimensions (NWB, NBR, IMX, KMX, NPROC, CLOSEC, ...)
#   2. Derived constituent numbering (mirrors input.F90 lines ~68-148)
#   3. Time/stability control block (TMSTRT..DLTADD)
#   4. Grid definition block (US/DS/UHS/DHS/NL/SLOPE/SLOPEC per branch,
#      LAT/LONGIT/ELBOT/BS/BE/JBDN per waterbody)
#   5. Initial condition block, per waterbody (T2I/ICEI/WTYPEC/GRIDC,
#      input.F90:703) -- added when porting the density equation of state
#      (Hydrodynamics/Density.jl), which needs WTYPEC to derive FRESH_WATER/
#      SALT_WATER, and resolves the TRAPEZOIDAL(JW) assumption flagged in
#      Core/InitGeometry.jl.
#
# NOT YET IN SCOPE (deferred to Tier 1, see project notes):
#   Meteorology cards, ice, turbulence-closure control, inflow/outflow cards,
#   structures, withdrawals, the full constituent (CST) block, output control
#   cards, kinetics rate constants. These all sit POSITIONALLY between the
#   initial condition block and the bathymetry filename (BTHFN) in
#   w2_con.csv -- there is no way to "skip past" them without knowing each
#   section's exact row count, which is exactly the Tier 1 parsing work.
#
# CONSEQUENCE: BTHFN (the bathymetry filename) is NOT read from w2_con.csv
# yet -- read_control_file does not reach that row. Pass the bathymetry file
# path explicitly to read_bathymetry (see BathymetryReader.jl) until Tier 1
# closes this gap.
#
# VALIDATION NOTE: the row-offset logic below (every skip!() and read_row!()
# call) was independently verified against the real Detroit w2_con.csv using
# a Python harness BEFORE this file was written, since no Julia runtime is
# available in the environment this was drafted in. Confirmed values:
#   NWB=1 NBR=4 IMX=31 KMX=117 NPROC=1 CLOSEC=OFF
#   US=[2,14,22,28] DS=[11,19,25,30] BS=1 BE=4 LAT=45.7299 LONGIT=122.177
# This file's STRUCTURE (the cursor logic, field order) is verified; running
# it through an actual Julia interpreter has not been done -- review before
# trusting it on a control file this hasn't been checked against.
# ==============================================================================

module InputReader

using ..W2Core: W2Global, W2Geometry, W2Names, W2TimeControl

export read_control_file, allocate_geometry!

"""
    allocate_geometry!(g::W2Global, geom::W2Geometry)

Sizes the per-segment (`IMX`) and per-segment/per-layer (`KMX x IMX`) arrays
of `geom` that `BathymetryReader.read_bathymetry!` writes into. Mirrors the
Fortran ALLOCATE step in input.F90 that follows the dimension read -- must be
called after `read_control_file` (needs `g.IMX`/`g.KMX`) and before
`read_bathymetry!`.

Only allocates the fields BathymetryReader currently populates (DLX, ELWS, B,
H). Other W2Geometry fields (H1, H2, BH1, ...) stay zero-length until
Core/InitGeometry.jl's `allocate_init_geometry!` sizes them -- see that
file's module docstring for the resolved H(K,JW) shape decision.
"""
function allocate_geometry!(g::W2Global, geom::W2Geometry)
    geom.DLX  = zeros(Float64, g.IMX)
    geom.ELWS = zeros(Float64, g.IMX)
    geom.B    = zeros(Float64, g.KMX, g.IMX)
    geom.H    = zeros(Float64, g.KMX, g.IMX)
    return geom
end

# ------------------------------------------------------------------------------
# A tiny stateful cursor over the CSV rows, mirroring Fortran's sequential
# READ(CON,*) semantics: every call consumes exactly one record (row), and
# list-directed reads of an array `(ARR(J),J=1,N)` consume ONE row holding N
# comma-separated values -- NOT N separate rows. (This tripped up the first
# draft of the Python validation harness; flagging it here so it doesn't trip
# up the next person editing this file.)
# ------------------------------------------------------------------------------
mutable struct CSVCursor
    rows::Vector{Vector{String}}
    pos::Int
end

function CSVCursor(path::AbstractString)
    rows = Vector{Vector{String}}()
    for line in eachline(path)
        fields = split(line, ',')
        # Trim trailing empty padding columns (w2_con.csv pads every row to
        # ~250 columns with trailing commas -- cosmetic Excel-export artifact).
        last_nonempty = findlast(f -> !isempty(strip(f)), fields)
        fields = last_nonempty === nothing ? String[] : fields[1:last_nonempty]
        push!(rows, [strip(f, [' ', '"']) for f in fields])
    end
    return CSVCursor(rows, 1)
end

"""Skip n records without reading -- mirrors a bare `READ (CON,*)` discard line."""
function skip!(c::CSVCursor, n::Int=1)
    c.pos += n
    return nothing
end

"""Read one record (row), returning its fields as strings."""
function read_row!(c::CSVCursor)
    row = c.rows[c.pos]
    c.pos += 1
    return row
end

"""
    read_list!(c, n)

Mirrors `READ (CON,*) (ARR(J), J=1,n)`: ONE row, first n comma-separated
fields. Pads with empty strings if the row is short (shouldn't happen on a
well-formed control file, but fails loudly via debug output rather than a
bounds error if it does).
"""
function read_list!(c::CSVCursor, n::Int)
    row = read_row!(c)
    if length(row) < n
        @warn "read_list!: row $(c.pos-1) has $(length(row)) fields, expected $n -- padding with empty strings" row
        row = vcat(row, fill("", n - length(row)))
    end
    return row[1:n]
end

"""Read a single scalar field from the next row (first field only)."""
read_scalar!(c::CSVCursor) = read_row!(c)[1]

onoff(s::AbstractString) = uppercase(strip(s)) == "ON"

"""
    expect_row!(c, expected)

Reads one row and checks its first `length(expected)` fields (case-insensitive,
whitespace-trimmed) match `expected`. This replaces a blind `skip!()` over the
control file's label rows -- w2_con.csv carries a label row before almost
every data row (sometimes the real field names, e.g. "NWB, NBR, IMX, KMX,
NPROC, CLOSEC"; sometimes just a generic "WB1,WB2,...,WB10" / "BR1,...,BR10"
cardinality marker). Either way, checking it catches cursor drift -- a
miscounted skip!() upstream, or a control file laid out slightly differently
than expected -- right where it happens, instead of producing silently wrong
numbers three sections later. This was Kunal's suggestion (read 27/06):
rather than add an Excel-reading dependency, the same self-checking benefit
is already sitting in the CSV's own label rows, just discarded until now.

Mismatches `@warn` rather than `error()` by default -- some label rows are
genuinely generic (WBn/BRn) and a strict failure there would be noise, not
signal. Pass `strict=true` to escalate to a hard error (recommended once
this has been run cleanly against a few more real control files).
"""
function expect_row!(c::CSVCursor, expected::Vector{String}; strict::Bool=false)
    row = read_row!(c)
    n = length(expected)
    got = [uppercase(strip(f)) for f in (length(row) >= n ? row[1:n] : vcat(row, fill("", n - length(row))))]
    want = uppercase.(strip.(expected))
    if got != want
        msg = "InputReader: label mismatch at row $(c.pos-1) -- expected $want, got $got"
        strict ? error(msg) : @warn msg
    end
    return row
end

# ------------------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------------------
"""
    read_control_file(path; debug=true)

Parses a CE-QUAL-W2 `w2_con.csv` file (the .csv branch of input.F90's
read logic) through the grid/waterbody definition block. Returns
`(g::W2Global, geom::W2Geometry, tc::W2TimeControl)`.

`debug=true` prints each section's parsed values as it goes -- this is the
"debug statements" the person asked for, meant to be diffed against the
known-good Detroit values during testing, not left on for production runs.
"""
function read_control_file(path::AbstractString; debug::Bool=true)
    c = CSVCursor(path)
    g = W2Global()
    tc = W2TimeControl()

    debug && @info "InputReader: opened $path ($(length(c.rows)) rows)"

    # --- Title block: 3 header/label lines, then 10 title lines ---
    skip!(c, 3)
    titles = [read_scalar!(c) for _ in 1:10]
    debug && @info "InputReader: title[1] = $(titles[1])"

    # --- Dimension block (input.F90 lines 40-49, .csv branch) ---
    skip!(c, 1); expect_row!(c, ["NWB", "NBR", "IMX", "KMX", "NPROC", "CLOSEC"])
    dims1 = read_list!(c, 6)
    g.NWB, g.NBR, g.IMX, g.KMX, g.NPROC = parse.(Int, dims1[1:5])
    g.CLOSEC = onoff(dims1[6])

    skip!(c, 1); expect_row!(c, ["NTR", "NST", "NIW", "NWD", "NGT", "NSP", "NPI", "NPU"])
    dims2 = read_list!(c, 8)
    g.NTR, g.NST, g.NIW, g.NWD, g.NGT, g.NSP, g.NPI, g.NPU = parse.(Int, dims2)

    skip!(c, 1); expect_row!(c, ["NGC", "NSS", "NAL", "NEP", "NBOD", "NMC", "NZP"])
    dims3 = read_list!(c, 7)
    g.NGC, g.NSS, g.NAL, g.NEP, g.NBOD, g.NMC, g.NZP = parse.(Int, dims3)

    # NOTE: the control file's own label calls this field "NDAY" -- the Fortran
    # *variable* is NOD (input.F90), so g.NOD below is correct, but don't expect
    # to see "NOD" in the file itself; that would be the wrong thing to assert.
    skip!(c, 1); expect_row!(c, ["NDAY", "SELECTC", "HABTATC", "ENVIRPC", "AERATEC", "INITUWL", "ORGCC", "SED_DIAG", "DZMAX"])
    dims4 = read_list!(c, 9)
    g.NOD = parse(Int, dims4[1])
    # dims4[2:8] = SELECTC, HABTATC, ENVIRPC, AERATEC, INITUWL, ORGCC, SED_DIAG
    # dims4[9]   = DZMAX -- not stored as a separate field per current State.jl
    # design (DZMAX already lives on W2Global with its Fortran default).
    dzmax = tryparse(Float64, dims4[9])
    if dzmax !== nothing && dzmax != 0.0
        g.DZMAX = dzmax
    end

    if g.NPROC == 0
        g.NPROC = 1   # mirrors input.F90 line 55: `if(NPROC == 0) NPROC=1`
    end

    debug && @info "InputReader: dimensions" g.NWB g.NBR g.IMX g.KMX g.NPROC g.CLOSEC g.NTR g.NST g.NIW g.NWD g.NGT g.NSP g.NPI g.NPU g.NGC g.NSS g.NAL g.NEP g.NBOD g.NMC g.NZP g.NOD g.DZMAX

    # --- Derived constituent numbering (mirrors input.F90 lines 68-148) ---
    # NOT computed here. input.F90 builds a running constituent index map
    # (NTDS=1, NGCS=2, NGCE=NGCS+NGC-1, NSSS=NGCE+1, ... through NATE, with a
    # branch on ORGC_CALC and an NBOD>0 stoichiometry expansion) that this
    # phase deliberately does not replicate -- g.NCT stays at its default (0)
    # until Tier 1 (water-quality.f90 / wqconstituents.F90) needs it and it
    # can be built and cross-checked against that file directly, rather than
    # guessed at here from the dimension counts alone.
    debug && @info "InputReader: g.NCT left at default (0) -- constituent numbering deferred to Tier 1"

    # --- Time control block (input.F90 lines 717-737, .csv branch) ---
    skip!(c, 1); expect_row!(c, ["TMSTRT", "TMEND", "YEAR"])
    tcrow = read_list!(c, 3)
    tc.TMSTRT = parse(Float64, tcrow[1])
    tc.TMEND = parse(Float64, tcrow[2])
    tc.YEAR = parse(Int, tcrow[3])

    skip!(c, 1); expect_row!(c, ["NDLT", "DLTMIN", "DLTINTER"])
    ndltrow = read_list!(c, 3)
    tc.NDLT = parse(Int, ndltrow[1])
    tc.DLTMIN = parse(Float64, ndltrow[2])
    tc.DLTINTER = onoff(ndltrow[3])

    skip!(c, 1); expect_row!(c, fill("DLTD", tc.NDLT)); tc.DLTD = parse.(Float64, read_list!(c, tc.NDLT))
    skip!(c, 1); expect_row!(c, fill("DLTMAX", tc.NDLT)); tc.DLTMAX = parse.(Float64, read_list!(c, tc.NDLT))
    skip!(c, 1); expect_row!(c, fill("DLTF", tc.NDLT)); tc.DLTF = parse.(Float64, read_list!(c, tc.NDLT))

    skip!(c, 1); expect_row!(c, ["WB1"])   # generic per-waterbody marker, not a real field name -- see expect_row! docstring
    tc.VISC = onoff.(read_list!(c, g.NWB))
    tc.CELC = onoff.(read_list!(c, g.NWB))
    tc.DLTADD = onoff.(read_list!(c, g.NWB))

    debug && @info "InputReader: time control" tc.TMSTRT tc.TMEND tc.YEAR tc.NDLT tc.DLTMIN tc.DLTINTER tc.DLTD tc.DLTMAX tc.DLTF

    # --- Grid definition block, per branch (input.F90 lines 740-748) ---
    geom = W2Geometry()
    skip!(c, 1); expect_row!(c, ["BR1"])   # generic per-branch marker
    g.US  = parse.(Int, read_list!(c, g.NBR))
    g.DS  = parse.(Int, read_list!(c, g.NBR))
    g.UHS = parse.(Int, read_list!(c, g.NBR))
    g.DHS = parse.(Int, read_list!(c, g.NBR))
    geom.NL = parse.(Int, read_list!(c, g.NBR))
    geom.SLOPE = parse.(Float64, read_list!(c, g.NBR))
    geom.SLOPEC = parse.(Float64, read_list!(c, g.NBR))

    debug && @info "InputReader: grid definition (per branch)" g.US g.DS g.UHS g.DHS geom.NL geom.SLOPE geom.SLOPEC

    # --- Waterbody definition block, per waterbody (input.F90 lines 753-758) ---
    skip!(c, 1); expect_row!(c, ["WB1"])   # generic per-waterbody marker
    geom.LAT    = parse.(Float64, read_list!(c, g.NWB))
    geom.LONGIT = parse.(Float64, read_list!(c, g.NWB))
    geom.ELBOT  = parse.(Float64, read_list!(c, g.NWB))
    g.BS   = parse.(Int, read_list!(c, g.NWB))
    g.BE   = parse.(Int, read_list!(c, g.NWB))
    g.JBDN = parse.(Int, read_list!(c, g.NWB))

    debug && @info "InputReader: waterbody definition" geom.LAT geom.LONGIT geom.ELBOT g.BS g.BE g.JBDN

    # --- Initial condition block, per waterbody (input.F90 line 703) ---
    skip!(c, 1); expect_row!(c, ["WB1"])   # generic per-waterbody marker
    geom.T2I    = parse.(Float64, read_list!(c, g.NWB))
    geom.ICEI   = parse.(Float64, read_list!(c, g.NWB))
    geom.WTYPEC = strip.(read_list!(c, g.NWB))
    geom.GRIDC  = strip.(read_list!(c, g.NWB))

    debug && @info "InputReader: initial conditions" geom.T2I geom.ICEI geom.WTYPEC geom.GRIDC
    debug && @info "InputReader: stopped after initial condition block, row $(c.pos) of $(length(c.rows)) -- BTHFN not yet reached (see module docstring)"

    return g, geom, tc
end

end # module InputReader
