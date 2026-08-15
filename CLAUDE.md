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

- ~~**`H(K,JW)` shape mismatch**~~ RESOLVED (2026-08-11): `geom.H` is now `KMX×IMX`, broadcast per waterbody across its segments (`BathymetryReader.jl` fills it; `allocate_geometry!` sizes it), matching every other geometry array instead of the raw Fortran's `KMX×NWB`. See `Core/InitGeometry.jl` module docstring.
- ~~**Free-surface elevation solve**~~ RESOLVED (2026-08-11): CONFIRMED sequential across connected branches, NOT parallel like the TKE/momentum TRIDIAG. Traced `w2_4_win.f90:896-1050` directly: branch JB's implicit tridiagonal system uses `Z(UHS(JB))`/`Z(DHS(JB))` (if `UH_INTERNAL`/`DH_INTERNAL`) as boundary terms, both belonging to a *different*, already-processed branch within the same timestep. `Core/Grid.jl`'s `branch_processing_order` (Kahn's-algorithm topological sort, not a numbering assumption -- validated against both Detroit and a synthetic reversed-numbering topology) now computes the correct per-waterbody processing order for any control file.
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

## MVP hydrodynamic run — scoping (2026-08-11), and progress against it

User asked what's needed for a minimal viable hydrodynamic run with base constituents.
Dependency-ordered critical path (each step blocks the next):

1. ~~`init-geom.F90` port~~ **DONE** (`Core/InitGeometry.jl`) — `EL`/`KTI`/`KTWB`/`KB`/`VOL`/
   `DEPTHB`/`DEPTHM`/areas, resolves the `H(K,JW)` shape question. Validated against real
   Detroit data (see `W2J_README.md` "What's Actually Implemented"). Known deferrals
   documented in the file's own module docstring — read that before building on top of
   it, especially the **approximated** (not faithfully ported) cross-branch boundary-
   width interpolation, since Detroit's 4-branches-in-1-waterbody topology does exercise it.
2. ~~`thomas_solve!` validation~~ **DONE** — cross-checked against a literal `TRIDIAG`
   transcription and Julia's dense solver.
3. ~~`Core/Grid.jl` branch connectivity~~ **DONE** — `BranchNetwork` + `build_branch_network`,
   the static branch-to-branch/branch-to-waterbody lookup (ported from the pattern at the
   top of every `waterbody.f90` ENTRY, e.g. `UPSTREAM_VELOCITY` at `waterbody.f90:36-42`,
   not the full runtime boundary-value logic those ENTRYs do — that's still
   `Hydrodynamics/Waterbody.jl`, a stub). Validated against real Detroit topology.
   `Core/InitGeometry.jl`'s `compute_bottom_layers!` refactored to consume it instead of
   its own duplicated `findfirst` search — `init_geometry!` now returns `(g, geom, net)`.
   ALSO added `branch_processing_order` (Kahn's-algorithm topological sort over branch
   dependencies within a waterbody) once tracing step 4 below showed it was needed —
   see that entry.
4. Free-surface elevation solve (continuity) + momentum solve (horizontal velocity, uses
   `TRIDIAG`) — **DONE, 2026-08-13** (`Hydrodynamics/FreeSurface.jl`), as REDUCED
   PHYSICS by deliberate choice (confirmed with user, 2026-08-12 — see that file's
   module docstring for the exact real-vs-stubbed line). Real: `RHO`/`P`/`GRAV`/`HPG`,
   the implicit free-surface tridiagonal solve with real branch sequencing via
   `branch_processing_order`, the explicit velocity update. Stubbed+flagged: `SB`/`ST`
   (needs `Turbulence.jl` + meteorology), `ADMX`/`ADMZ`/`DM` (real formulas exist, zero
   whenever `U=0`, port before trusting nonzero inflow), dam-flow branch case
   (unreachable for Detroit), implicit vertical-eddy-viscosity correction (needs
   `Turbulence.jl`). **Bug found+fixed along the way**: `input.F90:2265,2275`'s
   `H2(:,I) = H(:,JW)` baseline copy (every layer, not just top-active `KT`) was never
   ported into `BathymetryReader.jl`, leaving `geom.H2` zero below `KT` and causing
   `0.0/0.0` NaN in the velocity update's `UXBR/H2` term. Fixed in
   `Core/InitGeometry.jl`'s `allocate_init_geometry!` (`geom.H2 = copy(geom.H)`, since
   `geom.H` is already populated by `BathymetryReader.jl` before that function runs —
   cleaner than moving `H2`'s allocation earlier). Also revealed the earlier `VOL` test
   was too weak (finite-and-non-negative passes trivially for all-zero) —
   strengthened. Validated (`test/runtests.jl`): zero-flow sanity check against real
   Detroit data holds ELWS/U stable to `< 1e-9` over 5 timesteps.
5. Meteorology reading + surface heat exchange — needed because temperature couples into
   momentum via density. NOT started; also blocked on Tier-1 IO (met file isn't read yet).
6. Temperature transport (advection-diffusion, uses `TRIDIAG`) + minimal turbulence
   closure. NOT started.
7. `Simulation.jl` time-stepping driver tying 3-6 together — **DONE, 2026-08-13**
   (first-cut only: `run_zero_flow_sanity_check!`, ONE fixed `dlt = tc.DLTMAX[1]` for
   the whole run, no real adaptive-timestep logic, no boundary-condition IO, no
   kinetics). `IO/OutputWriter.jl` also **DONE, first-cut**: real per-segment TSR CSVs
   (`outputinitw2tools.F90:984-1057` naming convention) with only the columns this port
   actually computes (`JDAY,DLT(s),ELWS(m),U(ms-1)`) — other real TSR columns
   (`T2`/`Q`/`SRON`/`EXT`/...) deliberately absent, not fabricated as zeros. Validated
   end-to-end against real Detroit data: 5-step run writes real CSVs showing stable
   ELWS and zero U. **This is the "first cut running and getting TSR outputs for
   Detroit" milestone the user asked for (2026-08-13).**
8. Tier-1 boundary-condition IO — actual inflow/outflow/withdrawal/tributary time-series
   file reading (only the control-file *pointers* to these files are read today, via
   `tools/xlsx_to_yaml` into `w2_config.yaml` — `InputReader.jl` itself doesn't read them
   into `W2Global`/`W2Geometry` yet). NOT started.
9. Base constituent transport (e.g. TDS as a passive tracer) — reuses the transport
   scheme from step 6. NOT started.

Steps 1-4 and 7 (first-cut only) done. **Next: parallel processing** (per user's
explicit "then next is to prepare for parallel processing" instruction, 2026-08-13) —
not yet started. After that, real adaptive timestep + Tier-1 boundary IO + non-zero
forcing (ADMX/ADMZ/DM, SB/ST) are needed before `Simulation.jl` can run anything beyond
the zero-flow sanity check.

Context for how step 4 (density prerequisite) landed, from the prior session:
`Hydrodynamics/Density.jl`'s `density()` (the equation of state feeding the free-surface
solve's `BHRHO` term) is implemented and validated against known physical reference
values (fresh water's 4degC density maximum, 0/20/25degC points). Getting
`FRESH_WATER(JW)`/`SALT_WATER(JW)` right required extending `InputReader.jl` to read the
initial-condition block (`T2I`/`ICEI`/`WTYPEC`/`GRIDC`) and discovering those flags also
depend on the global `CONSTITUENTS` flag (from `CCC`) and `CAC(NTDS)` (Tier 1, not read
into structs yet) — not just `WTYPEC` alone. Rather than default those two to `true`,
`compute_water_type_flags!` takes them as required explicit keyword arguments (see
`W2J_README.md` "What's Actually Implemented" for the full reasoning). Also resolved
`TRAPEZOIDAL(JW)` from "assumed false" to actually read from `GRIDC`, with
`compute_areas_volumes!` now erroring loudly (test-proven, not just asserted) instead of
silently computing wrong RECT-formula areas for a TRAP-gridded waterbody.

**Parallel processing — DONE (first pass), 2026-08-14.** Per the user's explicit
"then next is to prepare for parallel processing" instruction. Two changes:

1. `Core/Grid.jl` gained `branch_processing_tiers(g, net, jw)` — the same Kahn's-
   algorithm dependency graph `branch_processing_order` already used, but returning the
   round-by-round "ready" sets *without* flattening them, so branches with no same-tier
   dependency on each other are identifiable as a group. `branch_processing_order` is
   now just `reduce(vcat, branch_processing_tiers(...))` — one graph construction, not
   two. Validated (`test/runtests.jl`): Detroit's tiers are `[[1],[2,3,4]]` (branch 1 has
   no internal dependency; branches 2-4 each depend only on branch 1) — three branches
   genuinely safe to run concurrently, not a coincidence of Detroit's specific topology
   (also validated against the existing synthetic reversed-numbering topology, tiers
   `[[4],[1,2,3]]`).
2. `Hydrodynamics/FreeSurface.jl`: `solve_free_surface!`'s per-branch body was pulled out
   into `solve_branch_free_surface!` and is now called under `Threads.@threads for jb in
   tier` for each tier from `branch_processing_tiers` — tiers still run strictly in
   order (each tier's `Threads.@threads` loop is an implicit barrier), so a later tier's
   cross-branch boundary read always sees the earlier tier's fully-updated `Z`/`ELWS`;
   only branches *within* one tier run concurrently, which is safe because each only
   reads a different (already-resolved, earlier-tier) branch's state and only ever
   writes its own segment range plus its own boundary pad. `NWB` (waterbodies) is also
   threaded — independent by construction (disjoint segment ranges). Every per-column
   computation (`compute_density_field!`, `compute_pressure_field!`,
   `compute_gravity_term!`, `compute_pressure_gradient!`, `update_velocities!`) is
   threaded over segments `i` too — each column only ever writes its own data (a couple
   read a neighbor column, none write one), confirmed embarrassingly parallel, matching
   the project's own `TRIDIAG` case-A analysis (Pillar 1).

**Validated, not just "should be correct"**: full test suite (484/484) passes at both
`--threads=1` and `--threads=4`. Beyond that, ran `hydrodynamic_step!` for 20 timesteps
against real Detroit data at both thread counts and confirmed **bit-identical** `ELWS`/
`U` output (`441.090000000000089` etc., matching to the last printed digit) — proof
there's no data race, not just that the looser sanity-check tests still pass.

**Regression found + fixed same day: threading was actually SLOWER for Detroit.**
User asked "have we tested the timestamps for parallel vs simple" — hadn't; only
correctness had been checked. Benchmarked `hydrodynamic_step!` wall-clock time
(500-step loop, real Detroit data): **0.136ms/step at 1 thread, climbing to
0.484ms/step at 8 threads** — parallel got *worse* as thread count went up, not better.
Root cause: Detroit's grid is small (IMX=31, at most ~31 columns / 1-4 branches per
loop), so `Threads.@threads`'s own spawn/scheduling overhead exceeds the actual
per-column work (a few dozen flops). Naively threading every loop (the first pass
above) was correct but not yet a win at this problem size.

**Fix**: `Core/Parallel.jl` (new file) — `parallel_foreach(f, range; threshold=...)`
picks serial vs. `Threads.@threads` execution PER CALL based on `length(range)` against
`PARALLEL_THRESHOLD`. Threshold empirically benchmarked (not guessed): a synthetic
per-column workload shaped like `update_velocities!`'s heaviest inner-K loop
(~100 layers) at `--threads=4` showed serial winning up to N~64-96, threaded starting to
win at N=128 (1.33x) and improving to 3.6x by N=4096 — set `PARALLEL_THRESHOLD = 128`
(clearly past the crossover, not sitting at the noisy breakeven point; documented as
machine-dependent in the file's own docstring, not a universal constant). Every
`Threads.@threads` site in `Hydrodynamics/FreeSurface.jl` (`compute_density_field!`,
`compute_pressure_field!`, `compute_gravity_term!`, `compute_pressure_gradient!`,
`update_velocities!`, and `solve_free_surface!`'s waterbody/branch-tier loops) now goes
through `parallel_foreach` instead. Since Detroit's loops never reach 128 columns, this
makes Detroit run every one of them serially — matching or beating the original
1-thread baseline at every thread count re-tested (post-fix: ~0.10-0.20ms/step
regardless of `--threads`, vs. the pre-fix 0.484ms/step blowup at 8 threads) — while a
future larger reservoir's bigger grid would cross the threshold and thread
automatically, no code change needed either way. This is the actual point of Pillar 1
generalizing beyond Detroit (same "GENERALITY REMINDER" discipline as
`branch_processing_order`/`branch_processing_tiers`).

**Also tested**: a full simulated year for Detroit (`TMSTRT=1` to `TMEND=365`,
`DLTMAX=1200s` → 26,208 steps) — runs in ~2.7s, and stays stable long-term (max `ELWS`
drift `1.16e-10` over the whole run, essentially floating-point noise, not slow
divergence) at both 1 and 4 threads. `test/runtests.jl` gained a 2000-step long-run
stability test (not the full 26,208, to keep the suite fast, but two orders of
magnitude more steps than the original 5-step smoke test) plus a
`Core/Parallel: parallel_foreach` unit test that exercises both the serial and threaded
code paths deterministically via the `threshold` override (not dependent on incidental
range lengths or `Threads.nthreads()` at test time). Full suite now 494/494.

**Next — allocation/GC overhead investigation (not started, current focus, 2026-08-14).**
User's original motivation for the Julia port was speedup — threading alone hasn't
delivered that at Detroit's scale (see above), so before assuming a large-grid
threading demo would be the next lever, the more likely culprit at THIS scale needs
checking first: `Hydrodynamics/FreeSurface.jl`'s `compute_density_field!`/
`compute_pressure_field!`/`compute_gravity_term!`/`compute_pressure_gradient!` each
allocate a fresh `zeros(Float64, kmx, imx)` (~3,600 elements for Detroit) EVERY
timestep just to re-zero and refill it, and `solve_branch_free_surface!` allocates six
more `IMX`-length arrays per branch per step — classic avoidable GC pressure in a hot
per-timestep loop, and unlike the threading threshold, not something whose payoff
depends on grid size. Plan (agreed with user, not yet executed): profile first
(`@allocated`/`--track-allocation`) to confirm this is actually the dominant cost
before changing anything — same "validate before implementing" discipline as
everything else in this project — then preallocate reusable buffers instead of
guessing.

Also still open: real adaptive timestep (DLTF/DLTMIN/DLTD breakpoints), Tier-1
boundary-condition IO, and non-zero forcing (ADMX/ADMZ/DM, SB/ST via `Turbulence.jl` +
meteorology) are needed before `Simulation.jl`/`FreeSurface.jl` can run anything beyond
the zero-flow sanity check — none of these started yet. A larger-grid parallel-speedup
demonstration (a synthetic reservoir big enough to cross `PARALLEL_THRESHOLD`) would
also be worth doing before claiming Pillar 1 delivers real wins, not just avoids losses
— not done yet, only the isolated per-column-workload benchmark that derived the
threshold itself.

**GENERALITY REMINDER (user, 2026-08-11)**: this model must generalize to other
reservoirs later, not just Detroit — do not hardcode values/assumptions beyond what the
Fortran source itself hardcodes. Concrete example from this session:
`branch_processing_order` could have naively looped `BS[jw]:BE[jw]` in stored numeric
order and it would have "worked" for Detroit, since Detroit's branches happen to be
numbered downstream-to-upstream — but nothing in the Fortran source requires that
numbering. Built a real topological sort instead and validated it against a synthetic
*reversed*-numbering topology (not just Detroit) specifically to catch this class of
bug — see `test/runtests.jl` "Core/Grid: branch_processing_order". Audited `W2J/src/`
for hardcoded Detroit-specific values (`364`, `441.09`, `117`, `NWB==1`, etc.) before
starting this work; found none outside of comments/docstrings and `test/` fixtures
(where Detroit-specific expected values are correct and expected). Keep applying this
check before every future module.

---

## Working discipline (enforce this)

- **Never read the whole Fortran source tree at once.** Read one Fortran file at a time, only when it's the immediate source for the Julia function being written.
- **Validate before writing Julia.** For any new parsing logic, write and run a Python/Julia one-liner against real Detroit data first. If the numbers don't match the validated values table above, fix the logic before writing the Julia function.
- **Ask for the Fortran file; don't assume.** If implementing a function that translates `foo.f90` and it hasn't been read this session, ask for it rather than guessing from memory of similar files.
- **Don't write stubs speculatively.** Every file in the package should either be implemented+validated, or explicitly marked stub with a one-line docstring saying what it's waiting for.
- **Keep the progress tracker in `W2J_README.md` current.** After any session where something moves from stub → implemented, update the tracker before ending.
