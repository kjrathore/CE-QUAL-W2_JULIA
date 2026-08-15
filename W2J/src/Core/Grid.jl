# ==============================================================================
# Core/Grid.jl
#
# Static branch/waterbody connectivity, built once from US/DS/UHS/DHS/BS/BE
# (already read by IO/InputReader.jl) into direct lookups instead of the
# `findfirst(jjb -> UHS(JB) in US(jjb):DS(jjb), ...)` search pattern that
# recurs at the top of nearly every ENTRY in waterbody.f90 (UPSTREAM_VELOCITY,
# UPSTREAM_WATERBODY, DOWNSTREAM_WATERBODY, UPSTREAM_BRANCH, DOWNSTREAM_BRANCH,
# UPSTREAM_FLOW, DOWNSTREAM_FLOW, UPSTREAM_CONSTITUENT, DOWNSTREAM_CONSTITUENT
# -- confirmed by reading waterbody.f90:36-42, the same "find JJB, then find
# JJW" boilerplate at the start of UPSTREAM_VELOCITY).
#
# SCOPE: this file owns the STATIC topology query (which branch/waterbody a
# given branch connects to). The actual per-timestep BOUNDARY VALUE
# computation those ENTRY points do (interpolating velocity/elevation/
# constituent profiles across a junction) is Hydrodynamics/Waterbody.jl's
# job, not this file's -- still a stub, not started this session.
#
# Already has one real caller: Core/InitGeometry.jl's `compute_bottom_layers!`
# duplicated this exact lookup ad hoc (to find JBUH/JBDH for the cross-branch
# KB-matching logic in init-geom.F90:221-250) before this file existed --
# refactored to call `build_branch_network` instead once it did.
# ==============================================================================

"""
    BranchNetwork

Per-branch connectivity, computed once by `build_branch_network`. Sentinel
`0` means "no connected branch" -- an external boundary (UHS/DHS == 0 or
-1), matching UQ_EXTERNAL/DQ_EXTERNAL/UH_EXTERNAL/DH_EXTERNAL in
`Core/InitGeometry.jl`.

Field names deliberately don't mirror Fortran variable names 1:1 (unlike
`W2Global`/`W2Geometry`) -- `JBUH`/`JBDH`/`JWUH`/`JWDH` are exactly the
per-branch arrays this struct's `upstream_branch`/`downstream_branch`/
`upstream_waterbody`/`downstream_waterbody` fields hold, so `g.JBUH` (a
plain field also still written directly by `compute_bottom_layers!` for
call sites that expect it there) and `net.upstream_branch` are the same
data by two names during the transition -- see that function's docstring.
"""
struct BranchNetwork
    segment_range::Vector{UnitRange{Int}}    # per branch JB: US[JB]:DS[JB]
    waterbody_of_branch::Vector{Int}         # per branch JB: which JW owns it
    segment_to_branch::Vector{Int}           # per segment I: which JB owns it (0 for the 1-segment boundary pad on each branch end, which belongs to no single branch)

    upstream_branch::Vector{Int}             # per branch JB: JJB s.t. UHS(JB) in US(JJB):DS(JJB), or 0
    upstream_waterbody::Vector{Int}
    downstream_branch::Vector{Int}           # per branch JB: JJB s.t. DHS(JB) in US(JJB):DS(JJB), or 0
    downstream_waterbody::Vector{Int}
end

"""
    build_branch_network(g)

Ports the branch-linkage lookup at init-geom.F90:221-250 (and repeated at
the top of every waterbody.f90 ENTRY) into a one-time precomputation. `g`
needs US/DS/UHS/DHS/BS/BE already populated (i.e. call after
`InputReader.read_control_file`).
"""
function build_branch_network(g)
    nbr = g.NBR
    segment_range = [g.US[jb]:g.DS[jb] for jb in 1:nbr]

    waterbody_of_branch = zeros(Int, nbr)
    for jw in 1:g.NWB, jb in g.BS[jw]:g.BE[jw]
        waterbody_of_branch[jb] = jw
    end

    segment_to_branch = zeros(Int, g.IMX)
    for jb in 1:nbr, i in g.US[jb]:g.DS[jb]
        segment_to_branch[i] = jb
    end

    upstream_branch = zeros(Int, nbr)
    upstream_waterbody = zeros(Int, nbr)
    downstream_branch = zeros(Int, nbr)
    downstream_waterbody = zeros(Int, nbr)

    find_branch(seg) = findfirst(jjb -> seg >= g.US[jjb] && seg <= g.DS[jjb], 1:nbr)

    for jb in 1:nbr
        if g.UHS[jb] > 0
            jjb = find_branch(g.UHS[jb])
            if jjb !== nothing
                upstream_branch[jb] = jjb
                upstream_waterbody[jb] = waterbody_of_branch[jjb]
            end
        end
        if g.DHS[jb] > 0
            jjb = find_branch(g.DHS[jb])
            if jjb !== nothing
                downstream_branch[jb] = jjb
                downstream_waterbody[jb] = waterbody_of_branch[jjb]
            end
        end
    end

    return BranchNetwork(segment_range, waterbody_of_branch, segment_to_branch,
                          upstream_branch, upstream_waterbody,
                          downstream_branch, downstream_waterbody)
end

"""
    branch_processing_order(g, net, jw)

Returns waterbody `jw`'s branches in an order safe for any solve whose
per-branch boundary condition draws on a *connected* branch's just-updated
state within the same timestep -- confirmed necessary for the free-surface
elevation solve (w2_4_win.f90:896-1050): branch JB's implicit tridiagonal
system uses `Z(UHS(JB))` (if `UH_INTERNAL(JB)`) and/or `Z(DHS(JB))` (if
`DH_INTERNAL(JB)`) as boundary terms, both of which belong to a *different*
branch and must already be freshly solved. This is a REAL sequential
dependency, unlike the confirmed-parallel TKE/momentum TRIDIAG -- see
CLAUDE.md "Open questions".

GENERALITY (the actual point of this function): naively looping branches
`BS[jw]:BE[jw]` in stored numeric order happens to work for Detroit only
because its branches were numbered in downstream-to-upstream order in the
control file -- nothing in the Fortran source *requires* that numbering,
and a different reservoir's control file has no reason to follow it. This
computes a real topological order via Kahn's algorithm on the dependency
graph {JB depends on `net.upstream_branch[JB]` if `UH_INTERNAL`, and on
`net.downstream_branch[JB]` if `DH_INTERNAL`, restricted to dependencies
within the same waterbody}, so it's correct for any branch numbering.
Errors loudly (does not guess) if the dependency graph has a cycle or an
unresolvable branch -- a genuinely malformed or cross-waterbody-dependent
control file, not something to silently paper over.

Sanity check against Detroit (own numbering happens to already be
topological, so this is a real cross-check, not a tautology): branches
2-4 each depend only on branch 1 (`net.downstream_branch == [0,1,1,1]`),
so Kahn's algorithm must produce `[1,2,3,4]` -- see test/runtests.jl.
"""
function branch_processing_order(g, net, jw)
    return reduce(vcat, branch_processing_tiers(g, net, jw); init=Int[])
end

"""
    branch_processing_tiers(g, net, jw) -> Vector{Vector{Int}}

Same dependency graph and Kahn's-algorithm construction as
`branch_processing_order` (see that docstring for the full derivation and
the free-surface-solve sequential-dependency proof), but returns the
round-by-round "ready" sets *without* flattening them into one order.

Every branch within one tier (inner vector) depends only on branches in
*earlier* tiers, never on another branch in the same tier -- by
construction of Kahn's algorithm, a branch only enters the `ready` set once
every dependency is already in a prior tier. This is the actual point of
this function, for `Hydrodynamics/FreeSurface.jl`'s parallel-processing
phase (CLAUDE.md "MVP hydrodynamic run", Pillar 1): branches in the same
tier can run under `Threads.@threads` with no data race, since each
branch's implicit tridiagonal solve only ever reads a *different* branch's
state (`net.upstream_branch`/`net.downstream_branch`, both necessarily in
an earlier tier) and only ever writes its own segment range plus its own
one-segment boundary pad -- never another branch's.

Sanity check against Detroit: branch 1 has no internal dependency
(`net.downstream_branch[1] == 0`, external at the dam) so tier 1 is `[1]`;
branches 2-4 each depend only on branch 1, so tier 2 is `[2,3,4]` -- three
branches genuinely safe to process concurrently, not just a coincidence of
Detroit's specific topology (see the reversed-numbering synthetic test in
`test/runtests.jl` for the generality proof, same as
`branch_processing_order`).
"""
function branch_processing_tiers(g, net, jw)
    branches = collect(g.BS[jw]:g.BE[jw])
    inset(jb) = jb in branches

    depends_on = Dict{Int,Vector{Int}}()
    for jb in branches
        deps = Int[]
        g.UH_INTERNAL[jb] && net.upstream_branch[jb] != 0 && inset(net.upstream_branch[jb]) &&
            push!(deps, net.upstream_branch[jb])
        g.DH_INTERNAL[jb] && net.downstream_branch[jb] != 0 && inset(net.downstream_branch[jb]) &&
            push!(deps, net.downstream_branch[jb])
        depends_on[jb] = deps
    end

    tiers = Vector{Vector{Int}}()
    resolved = Set{Int}()
    remaining = Set(branches)
    while !isempty(remaining)
        ready = sort([jb for jb in remaining if all(d -> d in resolved, depends_on[jb])])
        if isempty(ready)
            error("branch_processing_tiers: cyclic or unresolvable branch dependency in waterbody $jw among $(sort(collect(remaining)))")
        end
        push!(tiers, ready)
        union!(resolved, ready)
        setdiff!(remaining, ready)
    end
    return tiers
end
