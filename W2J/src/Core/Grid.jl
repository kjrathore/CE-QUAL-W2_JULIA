"""
Grid / branch-network topology.

Originating Fortran context: waterbody.f90's ENTRY points
(UPSTREAM_BRANCH, DOWNSTREAM_BRANCH, UPSTREAM_WATERBODY,
DOWNSTREAM_WATERBODY) define how segments connect across branch
boundaries — confirmed real coupling, not just USE-graph noise (see
project research notes, call-graph analysis of waterbody.f90).

This is exactly the piece flagged as the structural reason W2J is NOT
built on Oceananigans (Decision Log #1): a branched network of connected
channels is not a first-class concept in that framework's grid model.
It needs to be modeled directly here.

TODO: design BranchNetwork type once waterbody.f90's ENTRY points are
fully decomposed (Tier 0). Needs to represent: segment ranges per branch,
upstream/downstream branch connectivity, and junction behavior.
"""

# struct BranchNetwork
#     n_branches::Int
#     segment_range::Vector{UnitRange{Int}}   # which segments belong to each branch
#     upstream_branch::Vector{Int}            # connectivity — which branch feeds into this one, if any
#     downstream_branch::Vector{Int}
# end
