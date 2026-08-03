"""
Hydraulic structures — gates, spillways, pipes, and selective withdrawal.

Originating Fortran context: gate-spill-pipe.f90 (8 distinct callers —
second-highest call-graph hub in the codebase) and withdrawal.f90.
ENTRY points to decompose: PIPE_FLOW, OPEN_CHANNEL (+ dealloc variants)
in gate-spill-pipe.f90; DOWNSTREAM_WITHDRAWAL, LATERAL_WITHDRAWAL
(+ _ESTIMATE variants), Set_Flow_Fracs(2) in withdrawal.f90.

Status: Tier 1 — not yet started.

IMPORTANT — this is the cluster of physics flagged during the
Oceananigans-vs-from-scratch evaluation (Decision Log #1) as having no
off-the-shelf equivalent: gate/spillway/pipe head-discharge formulas and
the withdrawal-zone-of-influence calculation are dam/reservoir-specific
engineering submodels with nothing comparable in ocean-modeling
frameworks. This needs to be written as original physics regardless of
transport-engine choice — there is no library to lean on here.

TODO:
  1. Port the withdrawal-zone calculation (determines which vertical
     layers around an outlet contribute flow, based on outlet geometry
     and density stratification).
  2. Port head-discharge formulas for gates/spillways/pipes.
  3. Confirm how withdrawal volume removal couples to the free-surface
     elevation solve (see project research notes open question on the
     free-surface solve's sequential-dependency structure — not yet
     resolved).
"""

# function gate_flow(head, opening, ...) end
# function lateral_withdrawal(state, grid, structure) end
