"""
Turbulence closure — TKE (turbulent kinetic energy) vertical mixing.

Originating Fortran context: az.f90. Contains 4 near-identical TRIDIAG
call sites (lines 149, 310, 381, 424 in source), one per TKEBC(JW)
boundary-condition branch. CONFIRMED Case A (independent per-segment
systems, called inside `DO I` loop) — see project research notes.

Status: Tier 2 — not yet started.

TODO:
  1. Consolidate the 4 near-duplicate TRIDIAG blocks into one parallel
     loop with boundary-condition handling as a parameter/dispatch,
     rather than 4 separate Julia functions (flagged explicitly during
     the call-graph investigation — avoid reproducing the duplication).
  2. Wire into Solvers/Tridiagonal.jl's solve_vertical_systems!.
"""

# function calculate_tke!(arch, state, grid, tke_bc::Int) end
