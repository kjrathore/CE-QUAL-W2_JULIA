# ==============================================================================
# IO/OutputWriter.jl
#
# First-cut TSR (time-series) CSV writer, ported in spirit (not literally) from
# outputinitw2tools.F90:984-1057 -- one CSV file per output segment, named
# `<base_name>_seg<N>.csv`, matching the real per-segment TSR file convention
# (`TSRFN1(1:L1-1)//'_seg'//SEGNUM//'.'//ext`).
#
# SCOPE -- REDUCED, MATCHING Hydrodynamics/FreeSurface.jl's reduced-physics
# first cut (CLAUDE.md "MVP hydrodynamic run", step 7/8). The real TSR header
# is `JDAY,DLT(s),ELWS(m),T2(C),U(ms-1),Q(m3s-1),SRON(Wm-2),EXT(m-1),...`
# (outputinitw2tools.F90:1003) plus many more kinetics columns. This writer
# only emits JDAY, DLT(s), ELWS(m), U(ms-1) -- the columns this port actually
# computes today. T2/Q/SRON/EXT/etc are NOT written as fabricated zeros or
# placeholders; they're simply absent from the header, since nothing in this
# codebase computes them yet (no temperature transport, no boundary-flow IO,
# no meteorology). Add columns here only once the corresponding physics
# module exists -- don't pad the header ahead of the data.
# ==============================================================================

module OutputWriter

export TSRWriter, open_tsr_files, write_tsr_row!, close_tsr_files!

"""
    TSRWriter

Holds one open output stream per requested segment. `segments[i]` corresponds
to `streams[i]`.
"""
struct TSRWriter
    segments::Vector{Int}
    streams::Vector{IO}
end

"""
    open_tsr_files(output_dir, base_name, segments) -> TSRWriter

Opens one CSV file per segment in `segments` under `output_dir`, named
`<base_name>_seg<segment>.csv`, and writes the header row. Creates
`output_dir` if it doesn't exist.
"""
function open_tsr_files(output_dir::AbstractString, base_name::AbstractString, segments::Vector{Int})
    isdir(output_dir) || mkpath(output_dir)
    streams = IO[]
    for seg in segments
        path = joinpath(output_dir, "$(base_name)_seg$(seg).csv")
        io = open(path, "w")
        println(io, "JDAY,DLT(s),ELWS(m),U(ms-1)")
        push!(streams, io)
    end
    return TSRWriter(segments, streams)
end

"""
    write_tsr_row!(w, g, geom, jday, dlt)

Writes one row per open segment: JDAY, DLT(s), ELWS(m), and U(ms-1) at that
segment's current top active layer (`g.KTI[seg]`, per-segment since
`Core/InitGeometry.jl`'s `compute_water_surface_layer!` -- not the
waterbody-wide `g.KTWB`, since a segment's own top layer is what a TSR file
reports).
"""
function write_tsr_row!(w::TSRWriter, g, geom, jday::Real, dlt::Real)
    for (io, seg) in zip(w.streams, w.segments)
        kt = g.KTI[seg]
        u = g.U[kt, seg]
        println(io, "$(jday),$(dlt),$(geom.ELWS[seg]),$(u)")
    end
    return w
end

"""
    close_tsr_files!(w)

Closes every open stream. Call once at the end of a simulation run.
"""
function close_tsr_files!(w::TSRWriter)
    foreach(close, w.streams)
    return nothing
end

end # module OutputWriter
