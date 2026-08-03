# ==============================================================================
# IO/BathymetryReader.jl
#
# Reads a single waterbody's bathymetry file (the file named by BTHFN(JW) in
# w2_con.csv -- e.g. the Detroit example's ".\InputFiles\bth1.csv"), as read
# by input.F90 lines 2251-2292.
#
# STATUS: implemented from direct line-by-line translation of input.F90, but
# UNVALIDATED -- no sample bathymetry file was available to test against when
# this was written. Do not trust this against a real file until it's been
# run against one. Flagging exactly what's untested rather than presenting
# this with the same confidence as InputReader.jl's grid-definition block
# (which WAS independently checked against real Detroit data).
#
# Fortran detects two file formats by peeking at the first character:
#   '$' prefix -> "new" format: list-directed reads (READ(BTH,*)), one
#                 logical row per quantity (DLX, ELWS, PHI0, FRIC), then KMX
#                 rows of "H(K,JW)  B(K,I=US-1..DS+1)".
#   otherwise  -> "old" format: fixed-width FORMATTED reads, '(//(10F8.0))'
#                 -- skip 2 lines, then up to 10 values of 8 characters each
#                 per physical line, continuing across as many lines as
#                 needed. This is NOT comma-delimited despite the .csv
#                 extension some example files use -- it's column-positional.
#
# Both branches populate the same target arrays: DLX, ELWS, PHI0, FRIC (all
# indexed over the waterbody's full segment range, US(BS(JW))-1 : DS(BE(JW))+1
# -- i.e. including the one-segment boundary pad on each end) and B(K,I) / H(K,JW).
# ==============================================================================

module BathymetryReader

export read_bathymetry!

"""
    read_bathymetry!(geom, g, path, jw)

Reads the bathymetry file at `path` for waterbody `jw` (1-indexed) into
`geom` (a `W2Geometry`) and `g` (a `W2Global`, for `H`). Requires `g.US`,
`g.DS`, `g.BS`, `g.BE` to already be populated (i.e. call this AFTER
`InputReader.read_control_file`).

Segment range follows the Fortran convention exactly: `US(BS(jw))-1` through
`DS(BE(jw))+1` -- one boundary segment of padding on each end of the
waterbody's branch range.
"""
function read_bathymetry!(geom, g, path::AbstractString, jw::Int; debug::Bool=true)
    i_lo = g.US[g.BS[jw]] - 1
    i_hi = g.DS[g.BE[jw]] + 1
    n_seg = i_hi - i_lo + 1
    kmx = g.KMX

    debug && @info "BathymetryReader: waterbody $jw, segments $i_lo:$i_hi ($n_seg cells), $kmx layers, file=$path"

    if length(geom.DLX) < g.IMX || size(geom.B, 2) < g.IMX
        error("BathymetryReader: geom.DLX/geom.B are not sized to IMX=$(g.IMX) yet -- " *
              "allocate W2Geometry's per-segment/per-layer arrays before calling read_bathymetry!. " *
              "This phase doesn't include that allocate! step yet (see InputReader.jl scope note).")
    end

    lines = readlines(path)
    li = Ref(1)   # line cursor

    nextline!() = (l = lines[li[]]; li[] += 1; l)

    first_char = isempty(strip(lines[1])) ? "" : string(strip(lines[1])[1])

    # --- pull n floating values, however many physical lines that takes,
    #     where each line holds up to `per_line` fixed-width 8-char fields ---
    function read_fixed_floats!(n::Int; per_line::Int=10, width::Int=8)
        vals = Float64[]
        while length(vals) < n
            line = nextline!()
            nchunks = min(per_line, cld(length(line), width))
            for k in 0:(nchunks-1)
                lo = k*width + 1
                hi = min(lo + width - 1, length(line))
                lo > length(line) && break
                chunk = strip(line[lo:hi])
                isempty(chunk) && continue
                v = tryparse(Float64, chunk)
                v === nothing && continue
                push!(vals, v)
                length(vals) >= n && break
            end
        end
        return vals[1:n]
    end

    DLX = zeros(Float64, n_seg)
    ELWS = zeros(Float64, n_seg)
    PHI0 = zeros(Float64, n_seg)
    FRIC = zeros(Float64, n_seg)
    H = zeros(Float64, kmx)
    B = zeros(Float64, kmx, n_seg)

    if first_char == "\$"
        debug && @info "BathymetryReader: detected NEW ('\$') format"
        li[] = 1
        nextline!()             # the '$' marker/title line itself (e.g. "$Detroit Simplified Grid")
        nextline!()             # discarded by a bare READ(BTH,*) -- on real Detroit data this is a
                                 # segment-number header row (",1,2,3,...,31,,"), NOT actually blank,
                                 # but it's a list-directed read with no target variables either way,
                                 # so it's correctly discarded regardless of content. (Confirmed against
                                 # the real bth1.csv -- this comment was wrong about the content, even
                                 # though the discard logic itself was right.)
        DLX  .= parse.(Float64, split(nextline!(), ',')[2:end][1:n_seg])
        ELWS .= parse.(Float64, split(nextline!(), ',')[2:end][1:n_seg])
        PHI0 .= parse.(Float64, split(nextline!(), ',')[2:end][1:n_seg])
        FRIC .= parse.(Float64, split(nextline!(), ',')[2:end][1:n_seg])
        nextline!()             # discarded; real data has a "LAYERH...K" label row here
        for k in 1:kmx
            fields = split(nextline!(), ',')
            H[k] = parse(Float64, fields[1])
            B[k, :] .= parse.(Float64, fields[2:end][1:n_seg])
        end
    else
        debug && @info "BathymetryReader: detected OLD (fixed-width) format"
        li[] = 1
        nextline!(); nextline!()              # leading // skip
        DLX  .= read_fixed_floats!(n_seg)
        nextline!(); nextline!()
        ELWS .= read_fixed_floats!(n_seg)
        nextline!(); nextline!()
        PHI0 .= read_fixed_floats!(n_seg)
        nextline!(); nextline!()
        FRIC .= read_fixed_floats!(n_seg)
        nextline!(); nextline!()
        H .= read_fixed_floats!(kmx)
        for i in 1:n_seg
            nextline!(); nextline!()
            B[:, i] .= read_fixed_floats!(kmx)
        end
    end

    if i_lo < 1
        error("BathymetryReader: computed i_lo=$i_lo < 1 for waterbody $jw (US(BS(jw))=$(g.US[g.BS[jw]])) -- " *
              "this happens when a branch's upstream segment is 1, and hasn't been worked through yet " *
              "(Fortran's 0-indexed boundary-cell convention doesn't have an obvious 1-based Julia equivalent " *
              "in that case). Flagging rather than guessing at an offset.")
    end

    # Fill into the shared arrays at the right segment offset, for the fields
    # that actually have a home in W2Geometry today. Fortran's loop index
    # (US(BS(JW))-1 : DS(BE(JW))+1) already lines up directly with a 1-based
    # Julia Vector{}(1:IMX) -- no extra offset needed (verified against the
    # Detroit numbers: US(BS(1))-1 = US(1)-1 = 2-1 = 1, the natural first
    # index of a 1:IMX array, not 0).
    seg_range = i_lo:i_hi
    geom.DLX[seg_range]  .= DLX
    geom.ELWS[seg_range] .= ELWS
    geom.B[:, seg_range] .= B

    # PHI0, FRIC, and H(K,JW) don't have a struct field yet:
    #  - PHI0, FRIC are GLOBAL-module, per-segment (IMX) in the Fortran source
    #    -- not yet added to W2Global, since nothing else has needed them.
    #  - H(K,JW) is per-waterbody (KMX x NWB), and W2Geometry's H/H1/H2 fields
    #    are currently typed as per-segment (KMX x IMX) -- a real mismatch,
    #    not just a missing field. Needs a deliberate decision, not a quick
    #    patch, so it's returned directly instead of forced into the wrong
    #    shape.
    # Returning all four below so nothing is silently dropped; add the struct
    # fields once a second caller actually needs them persisted.

    debug && @info "BathymetryReader: read DLX[1:3]=$(DLX[1:min(3,end)]) H[1:3]=$(H[1:min(3,end)])"

    return (; DLX, ELWS, PHI0, FRIC, H, B)
end

end # module BathymetryReader
