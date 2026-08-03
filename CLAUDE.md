# W2J — Claude Code Context File

This file is the continuity layer from a prior Claude.ai session where the
architecture was designed and the IO/geometry base module was built. Read this
before touching any file. Do not read the entire Fortran source tree at startup.

---

## What this project is

A from-scratch Julia reimplementation of CE-QUAL-W2 (a 2D laterally-averaged
hydrodynamic and water quality model used by USACE for reservoir management).
Fortran source: `w2_source_5-21-2026/` (~48 files, ~45-48k LOC).
Test case: Detroit Reservoir example in `DetroitReservoir/`.

Three pillars drive every design decision:
1. **Parallel processing** — `Threads.@threads` over segments/cells (confirmed embarrassingly parallel for both `TRIDIAG` solves and per-cell kinetics).
2. **Differentiable programming** — `Enzyme.jl` for gradient-based calibration. Requires type-stable code and no global mutable state.
3. **Architectural flexibility** — kinetics module spatially agnostic (cell state in, rate terms out, no `(K,I)` hardcoding), so it can plug into a future transport engine without rewriting ecology math.

---

## Repo layout expected

```
W2J/
├── CLAUDE.md                        ← this file
├── README.md                        ← Fortran analysis: call-graph, ENTRY map, decision log
├── W2J_README.md                    ← Julia package design and progress tracker
├── w2_source_5-21-2026/             ← original Fortran source (read-only reference)
├── DetroitReservoir/
│   ├── w2_con.csv                   ← primary test control file
│   ├── InputFiles/bth1.csv          ← bathymetry (Detroit, $-format, confirmed)
│   └── w2_con_DetroitReservoir4_5.xlsm
└── src/
    ├── W2J.jl
    ├── Core/
    │   ├── State.jl                 ← UPDATED this session (see status below)
    │   ├── Grid.jl                  ← stub
    │   └── Architecture.jl          ← implemented
    ├── Solvers/
    │   └── Tridiagonal.jl           ← implemented, unvalidated
    ├── IO/
    │   ├── InputReader.jl           ← IMPLEMENTED this session (see status below)
    │   └── BathymetryReader.jl      ← IMPLEMENTED this session (see status below)
    ├── Plotting/
    │   └── LongitudinalProfile.jl   ← IMPLEMENTED this session, unvalidated
    ├── Hydrodynamics/               ← all stubs
    ├── WaterQuality/                ← all stubs
    └── Simulation.jl                ← stub
```

---

## Status of each file built this session

### `Core/State.jl` — skeleton, designed
Mirrors `w2modules.F90` MODULE GLOBAL / GEOMC / NAMESC. Key structs:
- `W2Global` — all dimension scalars (NWB, NBR, IMX, KMX, NPROC, ...) + connectivity arrays (US, DS, UHS, DHS, BS, BE, JBDN, KB, ...)
- `W2Geometry` — per-segment/per-layer geometry (DLX, ELWS, B, H, SLOPE, SLOPEC, ...) + LAT, LONGIT, ELBOT, NL added this session
- `W2TimeControl` — TMSTRT, TMEND, YEAR, NDLT, DLTD/DLTMAX/DLTF, VISC/CELC/DLTADD
- `W2Names` — constituent naming tables (stub)

Field names are kept **identical to Fortran** throughout — this is a hard rule, not a preference. The w2_con.csv label rows use the same names; users cross-reference W2 documentation directly.

**NPROC**: declared in Fortran as "INACTIVE" but we repurpose it to bound `Threads.@threads` loops. No new config field needed — existing control files already carry it.

**Not yet resolved**: `H(K,JW)` is per-waterbody (KMX×NWB) in Fortran but `W2Geometry.H` was typed KMX×IMX — a real mismatch flagged in BathymetryReader.jl. Do not silently patch; decide the right shape when `init-geom.F90` is ported.

### `IO/InputReader.jl` — implemented and validated against real Detroit data
Parses `w2_con.csv` through Phase A (dimensions → time control → grid/waterbody definition). Key design decisions:

- **Two-phase read**: reads dimension scalars first, then arrays sized from them — mirrors Fortran's ALLOCATE-then-READ pattern exactly.
- **One row = one Fortran READ(CON,*) call**: for arrays, all N values sit on ONE comma-separated row (not one-value-per-row). This is the most common mistake when reading W2 CSV files — `(ARR(J),J=1,NWB)` consumes one row of NWB values.
- **`expect_row!()` assertions**: every `skip!(c,2)` was replaced with `skip!(c,1)` + `expect_row!(c, [...])` which reads the label row and checks its field names. Catches cursor drift immediately rather than producing silently wrong numbers. All 12 assertions confirmed `[OK]` against real Detroit file.
- **CSV format only**: the legacy `.npt` fixed-width branch of `input.F90` is not ported. Every control file in this project is `.csv`.

**Scope boundary (important)**: Phase A stops after the grid/waterbody definition block (row 62 of 919 in the Detroit file). `BTHFN` (the bathymetry filename) is at row 868 — everything between is structures, withdrawals, the full constituent block, output control, kinetics rates. These are Tier 1 and have not been counted yet. Until Tier 1 is built, `BTHFN` is passed explicitly to `read_bathymetry!` rather than read from the control file.

**Validated values for Detroit**:
```
NWB=1  NBR=4  IMX=31  KMX=117  NPROC=1  CLOSEC=OFF
NTR=2  NST=1  NIW=0   NWD=2    NGT=0    NSP=0  NPI=0  NPU=0
NGC=1  NSS=2  NAL=0   NEP=0    NBOD=0   NMC=0  NZP=1
NOD=400  TMSTRT=1  TMEND=365  YEAR=2002  NDLT=1  DLTMAX=1200
US=[2,14,22,28]  DS=[11,19,25,30]  UHS=[0,0,0,0]  DHS=[0,9,10,11]
BS=1  BE=4  LAT=45.7299  LONGIT=122.177  ELBOT=364  JBDN=1
```

### `IO/BathymetryReader.jl` — implemented, validated against real bth1.csv
Reads the bathymetry file for one waterbody into `W2Geometry` (DLX, ELWS, B) and returns (DLX, ELWS, PHI0, FRIC, H, B) as a named tuple.

**`bth1.csv` format confirmed**: `$`-prefix (new format). Layout:
- Line 1: `$Detroit Simplified Grid,...` (title, discarded)
- Line 2: `,1,2,...,31,,` (segment-number header, discarded)
- Line 3: `DLX,0,844.9,1086.3,...` (label in col 0, then N values)
- Line 4: `ELWS,441.09,...`
- Line 5: `PHI0,0,2.382,...`
- Line 6: `FRIC,0.04,...`
- Line 7: `LAYERH,...,K,...` (discarded)
- Lines 8–124: `H_k, B_k_1, B_k_2, ..., B_k_31, k` (layer k, K index in last col, discarded)

Boundary segments (1, 12-13, 20-21, 26-27, 31) correctly read as zero-width.
Active layer counts per branch: BR1 segments 2-11 (63–110 layers), BR2 14-19 (63–86), BR3 22-25 (68–93), BR4 28-30 (77–97).

**Known gaps**: PHI0, FRIC, and H(K,JW) don't yet have struct homes. Returned in the named tuple; wire them into structs when `init-geom.F90` is ported. The guard against `i_lo < 1` will fire if US(BS(jw)) == 1 — that edge case hasn't been worked through.

### `Plotting/LongitudinalProfile.jl` — implemented, untested
Depends on `Plots.jl`. Takes `(g, geom, bathy_by_wb::Dict{Int,NamedTuple})`, draws side-view profiles per branch. Bottom elevation is approximated (active layer count × avg layer height) — not the real `KB`/`EL` computation from `init-geom.F90`. Good enough to catch bathymetry read errors, not good enough for a report figure.

**First thing to do in Claude Code**: wire these three together and run against Detroit data.

---

## Key Fortran analysis facts (do not re-derive, treat as confirmed)

**Call-graph centrality** (distinct callers):
`transport.f90`=9, `gate-spill-pipe.f90`=8, `waterbody.f90`=6, `systdg.f90`/`water-quality.f90`=4

**USE-dependency graph**: 43 nodes, 59 edges. `w2modules.F90` has 42 outgoing edges — literally every other file depends on it.

**ENTRY counts** (files with shared-SAVE state requiring decomposition):
- `water-quality.f90`: 60 ENTRY points + a bare `SAVE` at line 53 (every local var persists). Highest-risk file.
- `transport.f90`: 6 entries (`INTERPOLATION_MULTIPLIERS`, `HORIZONTAL_MULTIPLIERS{1}`, `VERTICAL_MULTIPLIERS{1}`, `DEALLOCATE_TRANSPORT`)
- `waterbody.f90`: 10 entries (9 named + DEALLOCATE)
- `heat-exchange.f90`: 3 entries (SHORT_WAVE_RADIATION, EQUILIBRIUM_TEMPERATURE, SURFACE_TERMS)
- `gate-spill-pipe.f90`: 4 entries (PIPE_FLOW, DEALLOCATE_PIPE_FLOW, OPEN_CHANNEL, DEALLOCATE_OPEN_CHANNEL)
- `withdrawal.f90`: 7 entries
- `time-varying-data.f90`: 3 entries (READ_INPUT_DATA, INTERPOLATE_INPUTS, DEALLOCATE_TIME_VARYING_DATA) + 36 `SLEEPQQ` calls (Windows real-time file polling for RESSIM coupling — replace with `sleep()` + mtime polling, or skip for standalone calibration runs)

**TRIDIAG**: Case A (independent per-segment) confirmed at both call sites (`az.f90` line 149/310/381/424, `w2_4_win.f90` line 1325). Safe for `Threads.@threads` with no algorithm changes. AT/VT/CT/DT arrays are per-(K,I) with no cross-I reach.

**T1/T2 array swap**: Fortran uses POINTER swapping to avoid copying the "old timestep" arrays. Julia doesn't need this — rebinding is free. Do not port the POINTER pattern; just use two plain arrays and swap references in the timestep loop.

**`input.F90` two-phase read pattern**: (1) read 6 dimension groups → (2) one giant ALLOCATE block (~150 lines) → (3) continue reading into now-allocated arrays. The Julia equivalent is `W2Global()` empty constructor + a separate `allocate!(g)` function, then continue reading. Phase A only covers steps (1) and the start of (3) for the geometry fields.

**`w2_con.csv` row convention**: label row (blank + field names) precedes every data row. Per-array reads consume ONE row of N comma-separated values, not N rows. All label rows confirmed via `expect_row!()` assertions. The file calls the `NOD` field `NDAY` in its own label — the Fortran variable is `NOD`, the file label is `NDAY`. Assertion uses `NDAY`.

---

## Architecture decisions (do not re-litigate without a strong reason)

| # | Decision |
|---|---|
| 1 | Build hydrodynamics from scratch, not on Oceananigans.jl (no hydraulic structures, no branch-network topology, wrong domain) |
| 2 | Kinetics module spatially agnostic — cell state in, rate terms out |
| 3 | `w2modules.F90` shared state → Julia structs (W2Global, W2Geometry, etc.), not globals |
| 4 | ENTRY subroutines → explicit functions with explicit state passed in |
| 5 | TRIDIAG → one shared `thomas_solve!` in `Solvers/Tridiagonal.jl`, not reimplemented per call site |
| 6 | `*TRM` arrays and `DO1/DO2/DO3` → explicit struct passed into constituent functions |
| 7 | Windows-specific code (`SLEEPQQ`, `DFWINTY`, `screen_output_intel.f90`) → replace, don't port |
| 8 | `NPROC` field repurposed to bound `Threads.@threads` thread count (was "INACTIVE" in Fortran) |
| 9 | CSV format only — `.npt` fixed-width branch not ported |
| 10 | Assertion-based label checking in IO reader (expect_row!) instead of blind positional skips |

---

## Open questions (do not guess at these — investigate from source before deciding)

- **`H(K,JW)` shape mismatch**: per-waterbody (KMX×NWB) in Fortran vs. per-segment in current W2Geometry. Resolve when porting `init-geom.F90`.
- **Free-surface elevation solve**: separate from the confirmed-parallel TKE/momentum TRIDIAG. Check `w2_4_win.f90` for sequential along-branch dependency before assuming it parallelizes.
- **Cross-group kinetics coupling**: `Kinetic_rates` mass-balance assembly — is there real same-pass coupling between Nutrients/Algae/Sediment groups? Paused mid-trace.
- **`LAM1` derivation**: needed before `RateMultipliers.jl` can be implemented for real. Trace from `TEMPERATURE_RATES` entry in `water-quality.f90`.
- **Phase B+ of InputReader**: everything between row 62 and the BTHFN row (868) in `w2_con.csv` — structures, withdrawals, full constituent block, output control, kinetics rates. Tier 1 work, not yet counted.

---

## Immediate next task — DONE (2026-08-03)

The first end-to-end test is wired up and runtime-verified against Julia 1.11.3 and
real Detroit data. `InputReader.allocate_geometry!(g, geom)` now exists
(`W2J/src/IO/InputReader.jl`) and sizes `geom.DLX`/`geom.ELWS` to `IMX` and `geom.B` to
`KMX×IMX` — call it between `read_control_file` and `read_bathymetry!`. `W2J.jl` now
`include`s `IO/BathymetryReader.jl` and `Plotting/LongitudinalProfile.jl` (both were
missing from the module before). `test/runtests.jl` has a regression test covering the
whole chain (11/11 passing); see `W2J_README.md` Progress Tracker for the one real bug
that surfaced (a parse-before-slice ordering bug in `BathymetryReader.jl`'s `$`-format
row reader) and the Project.toml/test-environment gotchas (placeholder UUID, `Test`
stdlib UUID) if this needs to be re-instantiated from scratch elsewhere.

**Next actual next task**: pick the next stub to implement — likely `Solvers/Tridiagonal.jl`'s
unit test (pure math, no W2 dependencies) or continuing Tier 1 `InputReader.jl` coverage
(structures/withdrawals/constituent block between row 62 and `BTHFN`).

---

## Working discipline (enforce this)

- **Never read the whole Fortran source tree at once.** Read one Fortran file at a time, only when it's the immediate source for the Julia function being written.
- **Validate before writing Julia.** For any new parsing logic, write and run a Python/Julia one-liner against real Detroit data first. If the numbers don't match the validated values table above, fix the logic before writing the Julia function.
- **Ask for the Fortran file; don't assume.** If implementing a function that translates `foo.f90` and it hasn't been read this session, ask for it rather than guessing from memory of similar files.
- **Don't write stubs speculatively.** Every file in the package should either be implemented+validated, or explicitly marked stub with a one-line docstring saying what it's waiting for.
- **Keep the progress tracker in `W2J_README.md` current.** After any session where something moves from stub → implemented, update the tracker before ending.
