"""
Transport — advection-diffusion of momentum, temperature, and constituents.

Originating Fortran context: transport.f90. The single most-depended-upon
physics file in the codebase (9 distinct callers — confirmed real call-graph
hub, the highest of any file — see project research notes). ENTRY points
to decompose: INTERPOLATION_MULTIPLIERS, HORIZONTAL_MULTIPLIERS(1),
VERTICAL_MULTIPLIERS(1), DEALLOCATE_TRANSPORT.

Status: Tier 0 — not yet started. HIGHEST priority hydrodynamics file
(see README "Base Module / Walking-Skeleton Plan" — this is the solver
pattern, not file, that's foundational: TRIDIAG usage here should call
Solvers/Tridiagonal.jl's solve_vertical_systems!, not a reimplementation).

TODO:
  1. Map INTERPOLATION_MULTIPLIERS / HORIZONTAL_MULTIPLIERS / VERTICAL_
     MULTIPLIERS shared state, same exercise as performed for
     water-quality.f90's TEMPERATURE_RATES (trace writes vs. reads).
  2. Confirm TRIDIAG call site here is also Case A (independent per-segment)
     — likely, given az.f90 and w2_4_win.f90 both confirmed Case A, but
     not yet explicitly checked for this file.
"""

# function transport!(arch, state, grid) end
