"""
Waterbody boundary conditions — branch and segment connectivity.

Originating Fortran context: waterbody.f90 (real call-graph hub, 6 distinct
callers — see project research notes). ENTRY points to decompose:
UPSTREAM_VELOCITY, UPSTREAM_WATERBODY, DOWNSTREAM_WATERBODY,
UPSTREAM_BRANCH, DOWNSTREAM_BRANCH, UPSTREAM_FLOW, DOWNSTREAM_FLOW,
UPSTREAM_CONSTITUENT, DOWNSTREAM_CONSTITUENT.

Status: Tier 0 — not yet started. Priority: HIGH (porting-order rationale —
see README, Decision Log #4 — this is one of the three files, alongside
transport.f90 and water-quality.f90, where ENTRY-shared state must be
mapped before porting begins).

TODO:
  1. Pull local-variable declarations for the WATERBODY host subroutine
     (same exercise as performed for water-quality.f90's KINETICS) to
     confirm whether these entries share write/read state the same way
     TEMPERATURE_RATES did, or whether they're more independent.
  2. Define explicit functions per entry, each accepting Core.Grid's
     BranchNetwork (once defined) rather than re-deriving connectivity.
"""

# function upstream_branch(grid, branch_id) end
# function downstream_branch(grid, branch_id) end
# function upstream_constituent(state, grid, branch_id) end
# function downstream_constituent(state, grid, branch_id) end
