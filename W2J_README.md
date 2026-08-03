# W2J — CE-QUAL-W2 in Julia

**W2J** is the package name for the from-scratch Julia reimplementation of CE-QUAL-W2
discussed in the project research notes (`README.md` at the project root — the porting
analysis, Fortran call-graph audit, ENTRY-decomposition findings, and decision log live
there; this file is about the Julia codebase itself).

**Status: 🚧 architecture scaffolding + IO base module. No physics implemented yet.**
`Pkg.instantiate()` / `Pkg.test()` now runtime-verified against Julia 1.11.3 — the
IO/geometry base stack (`InputReader` → `allocate_geometry!` → `BathymetryReader` →
`LongitudinalProfile`) runs end-to-end against real Detroit Reservoir data (11/11 tests
passing, see Progress Tracker).

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
├── tools/
│   ├── xlsx_to_yaml/                   — standalone Excel->YAML converter (own Project.toml: XLSX.jl only)
│   │   └── convert.jl                  — IMPLEMENTED + runtime-verified, see "Excel -> YAML conversion" below
│   └── review_config/                  — standalone w2_config.yaml validator (own Project.toml: YAML.jl only)
│       └── review_config.jl            — IMPLEMENTED + runtime-verified, see "Preprocessor checks ported" below
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
| `IO/InputReader.jl` | `input.F90` | Tier 2 | **Implemented + validated (Phase A + `allocate_geometry!`)** |
| `IO/BathymetryReader.jl` | `input.F90` (bathymetry read, `$`-format) | Tier 2 | **Implemented + validated** against real `bth1.csv` |
| `Plotting/LongitudinalProfile.jl` | n/a (debugging tool) | — | **Implemented + runtime-verified** end-to-end against Detroit data |
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

Verified working on Julia 1.11.3. Notes for anyone re-running this from scratch:
- `Project.toml`'s package UUID must be a real generated UUID, not `00000000-...` —
  Julia's loader silently refuses to resolve the all-zero placeholder.
- `test/Project.toml` needs an explicit `Test` stdlib entry (with its *real* UUID,
  `8dfed614-e22c-5e08-85e1-65c5234f0b40` — guessing this one is an easy way to lose an
  hour) for `Pkg.test()`'s sandboxed test environment to resolve.
- `Plots` is a real dependency now (`Plotting/LongitudinalProfile.jl`) — first
  `instantiate`/`test` will compile the full GR/Qt6 stack, budget a few minutes.

---

## Excel -> YAML conversion (`tools/xlsx_to_yaml/convert.jl`)

The real source of truth for a CE-QUAL-W2 run is the `.xlsm` control-file workbook, not
`w2_con.csv` -- the CSV is exported *from* the Excel sheet, and `w2_con.csv`'s
row-position-only format (no self-describing labels) is exactly why `IO/InputReader.jl`
could only safely cover Phase A (rows 1-62 of 919) without risking a silent row-count
mismatch. The Excel workbook's "w2_con.csv" sheet, by contrast, already carries real
structure in columns B (a section title or "VARNAME - description" label) and C onward
(either data or column headers) -- column A is pure commentary the workbook itself
flags as not written to the real CSV.

`tools/xlsx_to_yaml/convert.jl` exploits that structure to emit a single, human-editable,
commented `w2_config.yaml`. It is a **standalone script**, not part of the `W2J` package
-- it has its own `Project.toml` (just `XLSX.jl`) so the runtime package doesn't carry a
spreadsheet-reading dependency. Run it with:

```
cd W2J/tools/xlsx_to_yaml
julia --project=. -e "using Pkg; Pkg.instantiate()"   # first time only
julia --project=. convert.jl <path-to-.xlsm> <output.yaml>   # both args optional, default to Detroit
```

**How it classifies rows** (see the module docstring in `convert.jl` for full detail):
a row is a header ONLY when it immediately follows a blank-separator row (column B
missing, column C a lone space) -- confirmed as a hard invariant across the entire
Detroit workbook. From there: a repeated token ("DLTD","DLTD",...) is a Fortran array: a
run of generic per-instance tokens ("WB1","WB2",...,"PIPE1",...) is an entity block,
consumed field-row-by-field-row (key from column B, values from column C on) until the
next blank separator; anything else is disambiguated by *lookahead* -- exactly one
following data row means a dimension-scalar block (header names N fields, that one row
supplies N values 1:1), more than one means a vertical field block (same field-row-by-
field-row consumption, just not entity-shaped). Anything the script can't confidently
classify is preserved verbatim under `_raw_unclassified` (row number + cells) rather
than dropped or guessed at.

**Validated against the real Detroit workbook**: dimensions, `US`/`DS`/`BS`/`BE`/`JBDN`,
branch grid, waterbody location, structures, withdrawals, and every file-name block
round-trip correctly through a real YAML parser (`YAML.jl`) and match the values in
"Validated values for Detroit" in `CLAUDE.md`.

**`CST` / `CST DERI` constituent tables**: given a dedicated special case
(`emit_cst_block!`), since they don't fit the general column-B-is-the-key convention --
`CST`'s row identity sits in column B as a bare row *index integer* (constituent code is
column C), while `CST DERI`'s column B is blank entirely (its short name is the first
column-C token, same alignment as its header row -- these two tables are NOT shaped the
same way, confirmed by direct inspection, not assumed). Per constituent, keyed by its
**short name** (`CNAME2`/`CDNAME2`, e.g. `TDS`), with `CNAME`/`CDNAME` (the long name) kept
as `long_name`, `CAC` as `active`, `CMULT`/`CDMULT` as `multiplier`, and the remaining
columns grouped and sliced to the *actual* `NWB`/`NBR`/`NTR` counts (captured earlier in
the same conversion pass, from the "GRID/NPROC/CLOSE DIALOG BOX" and "IN/OUTFLOW"
dimension blocks) rather than the header's own placeholder-column width -- confirmed
against `W2manual455_Part3_InputOutputFiles_rev1.pdf` pp.156-157 (CST) and p.139 (CST
DERI). This also caught a real bug in the general converter: a stray bare-integer
annotation row (row 400 in the Detroit workbook, left over from an Excel formula) was
not a valid blank-separator OR a valid field row, and was silently merging the entire
CST table into the *preceding* entity block until `ends_block` was added to treat any
non-blank, non-string column B as a block terminator too.

**Known remaining limitation**: a handful of rows in the periphyton/epiphyton-group
block (`EPIC`/`EPRC`, ~15 rows) repeat the same extracted key in an *interleaved* pattern
(EPIC, EPRC, EPIC, EPRC, ...) rather than consecutive runs, which the "merge consecutive
same-key rows into an array" logic doesn't catch -- these currently overwrite each other
under a duplicate YAML key (YAML.jl warns but doesn't error; last value wins). Not yet
fixed -- flagging rather than guessing at a fix without checking the manual first.

`W2J`'s `InputReader.jl` does NOT yet read `w2_config.yaml` -- this script only produces
the file. Pointing `InputReader.jl` at YAML instead of `w2_con.csv` is future work.

**Other sections with prior-value-dependent column counts**: entity blocks (WBn/BRn/
PIPEn/...) are safe by construction -- `trimmed()` only reads as many columns as the
workbook actually populated, which already matches the real NWB/NBR/etc. counts (Excel
authored it that way), so there's no separate count to get wrong. CST/CST DERI needed
special-casing specifically because their row *identity* (not just column count) doesn't
follow the column-B-is-a-label convention every other block uses -- that's what to check
for before trusting a newly-encountered section, not column count by itself.

---

## Preprocessor checks ported (`tools/review_config/review_config.jl`)

The reference CE-QUAL-W2 preprocessor, `v455/preprocessor_422/preprocessor_4.22/pre_ivf_422.f90`
(13,725 lines), validates a control file before a run. It is NOT a small set of checks --
reconnaissance (via `Grep`, never a full read of the file) found **~1,174 distinct
validation sites** (`CALL ERRORS`/`CALL WARNINGS`, each paired with a `WRITE(ERR,...)`/
`WRITE(WRN,...)` message), spread across nearly every input section: dimensions,
branch/grid geometry, bathymetry, structures, withdrawals, gates, pipes, meteorology,
kinetics rate coefficients, output control, and more. Porting all of them is a
multi-session effort, not a single pass.

`tools/review_config/review_config.jl` (own `Project.toml`, just `YAML.jl`) ports a
**deliberately scoped first subset** -- only the sections `InputReader.jl` /
`xlsx_to_yaml/convert.jl` already produce data for. It reads `w2_config.yaml` and runs:

- **Dimensions/NPROC/CLOSEC** (pre_ivf_422.f90:3499-3512): NPROC range, CLOSEC ON/OFF.
- **Branch & grid geometry** (2794-2820, 3306-3335, 3489-3494, 8417-8427): BS/BE sum vs
  NBR, US(1)=2, inter-branch inactive-segment spacing, DS(last)+1=IMX, SLOPE>=SLOPEC,
  UHS/DHS validity, JBDN within its own waterbody's branch range.
- **CST constituent activation** (10383-10443): CAC/CPRWBC must be ON/OFF, initial
  concentration < -2 is an error, initial concentration = 0 while active is a warning
  (exempting age/residence-time tracers -- see the exemption-logic note in the file;
  this caught a real bug where the exemption checked only the short name and missed
  `Gen1`/"Water age, days", which the Fortran source's own equivalent check would exempt).

Run with `julia --project=. review_config.jl [path-to-w2_config.yaml]`. Against the real
Detroit `w2_config.yaml`: **0 errors, 1 warning** (`ISS1` starts at zero concentration
while active -- a legitimate, expected warning, not a bug).

**A clean run of this script does not mean the control file is fully valid** -- roughly
1,150 more checks from the reference preprocessor are not ported yet. Each check in
`review_config.jl` cites its source line range in `pre_ivf_422.f90` so it can be
independently verified or used as a starting point to add more. Extending scope
(structures/withdrawals, meteorology, kinetics rates, output control, ...) is future work.

---

## Progress Tracker

- [x] Fortran call-graph / dependency audit complete (see project research notes)
- [x] `ENTRY`-decomposition risk mapped for `water-quality.f90`, `transport.f90`,
      `waterbody.f90`, `heat-exchange.f90`, `time-varying-data.f90`,
      `gate-spill-pipe.f90`, `withdrawal.f90`
- [x] `TRIDIAG` confirmed Case A (independent per-segment) at both call sites
- [x] Package architecture scaffolded (this commit)
- [x] IO/geometry base module wired end-to-end and runtime-verified: `read_control_file`
      → `allocate_geometry!` → `read_bathymetry!` → `plot_longitudinal_profile`, all
      against real Detroit Reservoir data, with a regression test in `test/runtests.jl`
      (11/11 passing). Fixed a real bug in `BathymetryReader.jl` along the way: the
      `$`-format row parser was calling `parse.(Float64, ...)` on the *whole*
      comma-split row (including trailing blank padding columns) before slicing to
      `n_seg`, instead of slicing first — blew up on every real Detroit row.
- [ ] `thomas_solve!` validated against a known tridiagonal test case
- [ ] `LAM1` derivation traced; `RateMultipliers.jl` implemented for real
- [ ] First walking-skeleton run (simplified physics, single unbranched reservoir)
- [ ] First Fortran reference-output comparison
- [ ] Cross-group kinetics coupling investigation resumed (paused mid-trace — see
      project research notes Open Questions, and the note in `Kinetics.jl`)
- [ ] Free-surface elevation solve sequential-dependency check (separate from the
      already-confirmed-parallel TKE/momentum solve)
- [x] I/O reference mapping against v4.5 Excel doc: `tools/xlsx_to_yaml/convert.jl` reads
      the real `.xlsm` control-file workbook (the source of truth `w2_con.csv` is
      exported from) and emits a commented `w2_config.yaml`, validated against known
      Detroit values -- see "Excel -> YAML conversion" above for scope/limitations.
      `InputReader.jl` itself still reads `w2_con.csv`, not this YAML -- retargeting it
      is separate future work.
- [x] `LongitudinalProfile.jl` plot layout changed to a single-row, whole-grid
      (x = segment index 1:IMX) plot with one shared legend, output now written to
      `DetroitReservoir/Outputs/` instead of the repo root; `test/runtests.jl` updated
      to save there
