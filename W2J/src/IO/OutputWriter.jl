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
# emits JDAY, DLT(s), ELWS(m), U(ms-1), and (optionally, `include_temp=true`)
# T2(C) -- the columns this port actually computes. Q/SRON/EXT/etc are NOT
# written as fabricated zeros or placeholders; they're simply absent from the
# header, since nothing in this codebase computes them yet (no boundary flow
# reported per-segment, no meteorology). Add columns here only once the
# corresponding physics module exists -- don't pad the header ahead of the
# data. T2(C) was added 2026-08-23 once `Hydrodynamics/Transport.jl`'s
# `temperature_transport!` had a real per-timestep driver to call it from
# (`Simulation.jl`'s `run_forced_simulation!`) -- reads `g.T1[kt,seg]`, the
# just-solved new-timestep temperature (matches the real TSR's own "T2(C)"
# naming, which is confusingly the CURRENT temperature by the time output
# happens, not literally this port's `g.T1`/`g.HYD[:,:,4]` naming).
# ==============================================================================

module OutputWriter

export TSRWriter, open_tsr_files, write_tsr_row!, close_tsr_files!

"""
    TSRWriter

Holds one open output stream per requested segment. `segments[i]` corresponds
to `streams[i]`. `include_temp` controls whether `write_tsr_row!` also reads
and writes `g.T1` (see module docstring).
"""
struct TSRWriter
    segments::Vector{Int}
    streams::Vector{IO}
    include_temp::Bool
end

"""
    open_tsr_files(output_dir, base_name, segments; include_temp=false) -> TSRWriter

Opens one CSV file per segment in `segments` under `output_dir`, named
`<base_name>_seg<segment>.csv`, and writes the header row. Creates
`output_dir` if it doesn't exist. `include_temp=true` adds a `T2(C)` column
-- pass this only when the caller actually runs `temperature_transport!`
each step (see module docstring).
"""
function open_tsr_files(output_dir::AbstractString, base_name::AbstractString, segments::Vector{Int}; include_temp::Bool=false)
    isdir(output_dir) || mkpath(output_dir)
    streams = IO[]
    header = include_temp ? "JDAY,DLT(s),ELWS(m),T2(C),U(ms-1)" : "JDAY,DLT(s),ELWS(m),U(ms-1)"
    for seg in segments
        path = joinpath(output_dir, "$(base_name)_seg$(seg).csv")
        io = open(path, "w")
        println(io, header)
        push!(streams, io)
    end
    return TSRWriter(segments, streams, include_temp)
end

"""
    segment_waterbody(g, seg) -> jw

Finds which waterbody `seg` belongs to by membership in that waterbody's
full segment range (`US(BS(jw))-1 : DS(BE(jw))+1`, the same boundary-padded
range `IO/BathymetryReader.jl`/`Core/InitGeometry.jl` use elsewhere). `NWB`
is small (1-4 in every case this port has been validated against), so a
linear scan is cheap -- no need to thread `Core/Grid.jl`'s `BranchNetwork`
through just for this lookup.
"""
function segment_waterbody(g, seg::Int)
    for jw in 1:g.NWB
        lo = g.US[g.BS[jw]] - 1
        hi = g.DS[g.BE[jw]] + 1
        lo <= seg <= hi && return jw
    end
    error("segment_waterbody: segment $seg not found in any waterbody's range")
end

"""
    write_tsr_row!(w, g, geom, jday, dlt)

Writes one row per open segment: JDAY, DLT(s), ELWS(m), (T2(C) if
`w.include_temp`), and U(ms-1) at the WATERBODY's top active layer
`g.KTWB[jw]` (confirmed against real Fortran, `outputa2w2tools.F90`'s TSR
row write uses `U(K,I)`/`T1(K,I)` with `K` set from `KTWB(JW)`, e.g. line
180's `DO K=KTWB(JW),KB(I)` default case -- NOT `g.KTI[seg]`, the bathymetry
sub-layer index, which is a different quantity entirely used for computing
partial-top-layer width `BI`, not for indexing `U`/`T1`. FOUND AS A REAL
BUG, not a hypothetical: an early `run_forced_simulation!` validation run
against real DET data (2026-08-23) showed `U(ms-1)` reading `0.0` for every
row despite real QIN/QOT forcing -- traced to this exact mismatch (`g.KTI`
on Detroit/DET's bathymetry is frequently a K index BELOW the true top
active layer, an unrelated/mostly-empty row of `U`).
"""
function write_tsr_row!(w::TSRWriter, g, geom, jday::Real, dlt::Real)
    for (io, seg) in zip(w.streams, w.segments)
        kt = g.KTWB[segment_waterbody(g, seg)]
        u = g.U[kt, seg]
        if w.include_temp
            t2 = g.T1[kt, seg]
            println(io, "$(jday),$(dlt),$(geom.ELWS[seg]),$(t2),$(u)")
        else
            println(io, "$(jday),$(dlt),$(geom.ELWS[seg]),$(u)")
        end
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
