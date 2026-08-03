# W2J — CE-QUAL-W2 in Julia

**W2J** is the package name for the from-scratch Julia reimplementation of CE-QUAL-W2
discussed in the project research notes (`README.md` at the project root — the porting
analysis, Fortran call-graph audit, ENTRY-decomposition findings, and decision log live
there; this file is about the Julia codebase itself).

**Status: 🚧 architecture scaffolding only. No physics implemented yet.**
Syntax has not been runtime-verified (no Julia available in the environment this was
drafted in) — review carefully before first `instantiate`/`test`.

---

## The Three Pillars (restated briefly — full rationale in the project research notes)

1. **Parallel processing** — multi-core speedup via per-cell kinetics and per-segment
   tridiagonal solves (confirmed embarrassingly parallel for both).
2. **Differentiable programming** — gradient-based calibration via `Enzyme.jl`.
3. **Architectural flexibility** — kinetics decoupled from the hydrodynamics/grid layer.

---

## Design Principles (how the pillars become code structure)

- **No global mutable state.** Everything threaded explicitly through structs
  (`Core/State.jl`), replacing `w2modules.F90` + the `ENTRY`-shared-state pattern that
  has no Julia equivalent.
- **Kinetics functions are spatially agnostic.** Local cell state in, rate terms out —
  no `(K,I)` indexing assumptions baked into the water-quality math itself. This is what
  keeps Pillar 3 alive without committing to building 3D support now.
- **Shared numerical primitives are ported once.** The per-column tridiagonal solver
  (`Solvers/Tridiagonal.jl`) is implemented a single time and reused everywhere `TRIDIAG`
  was called in the original (currently: the TKE closure, the momentum solve, and likely
  constituent transport) — not reimplemented per call site as in the Fortran source.
- **Architecture (CPU/GPU) is an explicit, swappable argument**, not baked into kernel
  logic — modeled on Oceananigans.jl's `CPU`/`GPU{D}` pattern, even though W2J does not
  depend on Oceananigans (see project research notes, Decision Log #1, for why).

---

## Package Structure

```
W2J/
├── Project.toml
├── src/
│   ├── W2J.jl                          — module entry point
│   ├── Core/
│   │   ├── Architecture.jl             — CPU/GPU dispatch type (implemented)
│   │   ├── State.jl                    — shared-state struct design (skeleton)
│   │   └── Grid.jl                     — branch-network topology (skeleton)
│   ├── Solvers/
│   │   └── Tridiagonal.jl              — per-column Thomas algorithm (IMPLEMENTED — see status table)
│   ├── Hydrodynamics/
│   │   ├── Waterbody.jl                — branch/segment boundary conditions (stub)
│   │   ├── Transport.jl                — advection-diffusion (stub)
│   │   ├── Turbulence.jl                — TKE closure (stub)
│   │   └── Structures.jl               — gates/spillways/pipes/withdrawal (stub)
│   ├── WaterQuality/
│   │   ├── RateMultipliers.jl          — temperature rate multipliers (skeleton, formulas confirmed)
│   │   ├── Kinetics.jl                 — ENTRY decomposition orchestration (stub)
│   │   └── Constituents/
│   │       ├── Nutrients.jl
│   │       ├── Organics.jl
│   │       ├── Algae.jl
│   │       └── Sediment.jl
│   ├── IO/
│   │   ├── InputReader.jl              — stub
│   │   └── OutputWriter.jl             — stub
│   └── Simulation.jl                   — driver loop / walking-skeleton orchestrator (stub)
├── test/
│   ├── runtests.jl
│   └── reference_cases/                — empty; needs a Fortran reference run to compare against
└── docs/                               — empty
```

---

## Module → Fortran Source Map

| W2J module | Originating Fortran file(s) | Porting tier | Status |
|---|---|---|---|
| `Core/State.jl` | `w2modules.F90` | Foundation | Skeleton |
| `Core/Grid.jl` | `waterbody.f90` (branch connectivity entries) | Foundation | Skeleton |
| `Solvers/Tridiagonal.jl` | `TRIDIAG` (called from `az.f90`, `w2_4_win.f90`) | Foundation | **Implemented, unvalidated** |
| `Hydrodynamics/Waterbody.jl` | `waterbody.f90` | Tier 0 | Stub |
| `Hydrodynamics/Transport.jl` | `transport.f90` | Tier 0 | Stub |
| `Hydrodynamics/Turbulence.jl` | `az.f90` | Tier 2 | Stub |
| `Hydrodynamics/Structures.jl` | `gate-spill-pipe.f90`, `withdrawal.f90` | Tier 1 | Stub |
| `WaterQuality/RateMultipliers.jl` | `water-quality.f90` (`TEMPERATURE_RATES` entry) | Tier 1 | Skeleton, formulas confirmed |
| `WaterQuality/Kinetics.jl` + `Constituents/` | `water-quality.f90` (~70 entries) | Tier 1 | Stub |
| `IO/InputReader.jl` | `input.F90` | Tier 2 | Stub |
| `IO/OutputWriter.jl` | `output.f90`, `outputa2w2tools.F90`, `outputinitw2tools.F90` | — | Stub |
| `Simulation.jl` | `w2_4_win.f90` | Tier 3 (file) / Foundation (skeleton role) | Stub |

---

## What's Actually Implemented vs. What's a Placeholder

Be precise about this so nothing gets mistaken for working code:

- **`Solvers/Tridiagonal.jl`** — real Thomas algorithm + threaded batch dispatch, written
  out fully. **Not yet validated** against `TRIDIAG`'s actual output for a known input —
  do that before relying on it.
- **`WaterQuality/RateMultipliers.jl`** — formulas are direct transcriptions from the
  confirmed Fortran trace (not guessed), but the function itself is a stub because
  `LAM1`'s own derivation was never traced (flagged honestly in the file).
- **Everything else** — docstrings only, describing what each piece corresponds to,
  its confirmed/unconfirmed findings so far, and what needs to happen before real code
  goes in. Treat these as a map, not an implementation.

---

## Build / Test

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. -e "using Pkg; Pkg.test()"
```

(Untested in this environment — no Julia runtime available where this was drafted.)

---

## Progress Tracker

- [x] Fortran call-graph / dependency audit complete (see project research notes)
- [x] `ENTRY`-decomposition risk mapped for `water-quality.f90`, `transport.f90`,
      `waterbody.f90`, `heat-exchange.f90`, `time-varying-data.f90`,
      `gate-spill-pipe.f90`, `withdrawal.f90`
- [x] `TRIDIAG` confirmed Case A (independent per-segment) at both call sites
- [x] Package architecture scaffolded (this commit)
- [ ] `thomas_solve!` validated against a known tridiagonal test case
- [ ] `LAM1` derivation traced; `RateMultipliers.jl` implemented for real
- [ ] First walking-skeleton run (simplified physics, single unbranched reservoir)
- [ ] First Fortran reference-output comparison
- [ ] Cross-group kinetics coupling investigation resumed (paused mid-trace — see
      project research notes Open Questions, and the note in `Kinetics.jl`)
- [ ] Free-surface elevation solve sequential-dependency check (separate from the
      already-confirmed-parallel TKE/momentum solve)
- [ ] I/O reference mapping against v4.5 Excel doc (if available) + source files
