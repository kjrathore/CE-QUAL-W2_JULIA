# W2J — CE-QUAL-W2 in Julia

**W2J** is the package name for the from-scratch Julia reimplementation of CE-QUAL-W2
discussed in the project research notes (`README.md` at the project root — the porting
analysis, Fortran call-graph audit, ENTRY-decomposition findings, and decision log live
there; this file is about the Julia codebase itself).

**Status: 🚧 first-cut reduced-physics hydrodynamic run working end-to-end.**
`Pkg.instantiate()` / `Pkg.test()` runtime-verified against Julia 1.11.3 (481/481 tests
passing). The IO/geometry base stack (`InputReader` → `allocate_geometry!` →
`BathymetryReader` → `init_geometry!`) AND a reduced-physics free-surface/momentum
solve (`Hydrodynamics/FreeSurface.jl`) AND a time-stepping driver
(`Simulation.run_zero_flow_sanity_check!`) all run end-to-end against real Detroit
Reservoir data, producing real per-segment TSR CSV output (`IO/OutputWriter.jl`) that
shows a stable water surface under a zero-boundary-flow sanity check — see Progress
Tracker. Parallel processing is the explicitly-agreed next phase, not yet started.

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
│   │   ├── State.jl                    — shared-state struct design (IMPLEMENTED — 130/50 fields, W2Global/W2Geometry)
│   │   ├── Grid.jl                     — branch-network topology (IMPLEMENTED + validated — see status table)
│   │   └── InitGeometry.jl             — geometry finalization: EL/KB/KTI/VOL/... (IMPLEMENTED — see status table)
│   ├── Solvers/
│   │   └── Tridiagonal.jl              — per-column Thomas algorithm (IMPLEMENTED + VALIDATED — see status table)
│   ├── Hydrodynamics/
│   │   ├── Density.jl                  — water density equation of state (IMPLEMENTED + validated)
│   │   ├── FreeSurface.jl              — free-surface elevation + momentum solve (IMPLEMENTED, REDUCED PHYSICS — see status table)
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
│   │   ├── InputReader.jl              — IMPLEMENTED + validated (Phase A + allocate_geometry!)
│   │   ├── BathymetryReader.jl         — IMPLEMENTED + validated against real bth1.csv
│   │   └── OutputWriter.jl             — IMPLEMENTED, first-cut TSR CSV writer (reduced columns — see status table)
│   └── Simulation.jl                   — IMPLEMENTED, first-cut fixed-DLT zero-flow driver (see status table)
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
| `Core/State.jl` | `w2modules.F90` | Foundation | **Implemented** (grown as later modules need fields) |
| `Core/Grid.jl` | `waterbody.f90:36-42` (branch-linkage lookup pattern, not the full ENTRY points) | Foundation | **Implemented + validated** against real Detroit topology |
| `Core/InitGeometry.jl` | `init-geom.F90` + `init.F90:205-249,415-418` | Foundation | **Implemented + validated** against real Detroit data (see below) |
| `Solvers/Tridiagonal.jl` | `TRIDIAG` (called from `az.f90`, `w2_4_win.f90`) | Foundation | **Implemented + validated** (200 randomized cases + 1 fixed case, cross-checked against a literal `TRIDIAG` transcription and Julia's dense solver) |
| `Hydrodynamics/Density.jl` | `density.f90` | Foundation | **Implemented + validated** against known physical reference values |
| `Hydrodynamics/FreeSurface.jl` | `w2_4_win.f90` (free-surface/momentum block, :879-1327) | Foundation | **Implemented, REDUCED PHYSICS, validated** — see "What's Actually Implemented" below |
| `Hydrodynamics/Waterbody.jl` | `waterbody.f90` | Tier 0 | Stub |
| `Hydrodynamics/Transport.jl` | `transport.f90` | Tier 0 | Stub |
| `Hydrodynamics/Turbulence.jl` | `az.f90` | Tier 2 | Stub |
| `Hydrodynamics/Structures.jl` | `gate-spill-pipe.f90`, `withdrawal.f90` | Tier 1 | Stub |
| `WaterQuality/RateMultipliers.jl` | `water-quality.f90` (`TEMPERATURE_RATES` entry) | Tier 1 | Skeleton, formulas confirmed |
| `WaterQuality/Kinetics.jl` + `Constituents/` | `water-quality.f90` (~70 entries) | Tier 1 | Stub |
| `IO/InputReader.jl` | `input.F90` | Tier 2 | **Implemented + validated** (Phase A + `allocate_geometry!` + initial-condition block: `T2I`/`ICEI`/`WTYPEC`/`GRIDC`) |
| `IO/BathymetryReader.jl` | `input.F90` (bathymetry read, `$`-format) | Tier 2 | **Implemented + validated** against real `bth1.csv` |
| `Plotting/LongitudinalProfile.jl` | n/a (debugging tool) | — | **Implemented + runtime-verified** end-to-end against Detroit data |
| `IO/OutputWriter.jl` | `output.f90`, `outputa2w2tools.F90`, `outputinitw2tools.F90` | — | **Implemented, first-cut** — per-segment TSR CSV, reduced columns (JDAY/DLT/ELWS/U only — see below) |
| `Simulation.jl` | `w2_4_win.f90` | Tier 3 (file) / Foundation (skeleton role) | **Implemented, first-cut** — fixed-DLT zero-flow driver (`run_zero_flow_sanity_check!`), see below |

---

## What's Actually Implemented vs. What's a Placeholder

Be precise about this so nothing gets mistaken for working code:

- **`Solvers/Tridiagonal.jl`** — real Thomas algorithm + threaded batch dispatch, written
  out fully and **validated**: cross-checked against Julia's dense `\` solver AND a
  literal transcription of the real `TRIDIAG` subroutine (transport.f90:572-593) over
  200 randomized diagonally-dominant systems plus one hand-verified fixed case (see
  `test/runtests.jl`).
- **`Core/InitGeometry.jl`** — ports `init-geom.F90` (641 lines): layer elevations
  `EL(K,I)`, top/bottom active layers `KTI`/`KTWB`/`KB`, cell volumes/depths/areas, and
  the branch boundary-condition flags it depends on from `init.F90`. **Validated**
  against real Detroit data: `EL[KMX,:] == ELBOT`, `KB` matches the documented range,
  volumes are finite and non-negative, and the model's own "upstream active segment"
  logic (`CUS`) correctly detects that branch 1's segments 2-5 are dry at the Jan-1
  drawdown pool while segments 6-11 near the dam stay wet -- a real, physically-sensible
  edge case this session initially mistook for a bug before checking directly (see
  `Core/InitGeometry.jl`'s module docstring for the full deferral list: RESTART_IN,
  TRAPEZOIDAL grids, `w2_constriction.csv`, and -- the one real accuracy gap, not just a
  skip -- an **approximated** cross-branch boundary-width interpolation that copies the
  adjacent branch's column instead of the real elevation-matched interpolation).
  This also resolves the `H(K,JW)` shape question flagged in `CLAUDE.md`: `geom.H` is
  now populated as `KMX x IMX` (broadcast per waterbody across its segments), matching
  every other geometry array, rather than the raw Fortran's `KMX x NWB`.
- **`Core/Grid.jl`** — `BranchNetwork` + `build_branch_network`, a one-time
  precomputation of the branch-to-branch/branch-to-waterbody connectivity lookup that
  otherwise gets re-derived via `findfirst` at the top of nearly every `waterbody.f90`
  ENTRY point (confirmed by reading `UPSTREAM_VELOCITY`, `waterbody.f90:36-42` -- the
  same "find JJB, then find JJW" pattern). `Core/InitGeometry.jl`'s
  `compute_bottom_layers!` was refactored to use it instead of its own inline copy of
  the same search. **Validated** against real Detroit topology: branches 2-4 correctly
  resolve to branch 1 as their downstream connection at segments 9/10/11, branch 1
  correctly resolves to no downstream branch (external, at the dam), and no branch has
  an internal upstream connection (`UHS` all 0 for Detroit). Scope note: this file owns
  the *static* topology query only -- the actual per-timestep boundary-value
  interpolation those `waterbody.f90` ENTRY points do (velocity/elevation/constituent
  profiles across a junction) is `Hydrodynamics/Waterbody.jl`'s job, still a stub.
  ALSO added `branch_processing_order` -- a Kahn's-algorithm topological sort answering
  the "Free-surface elevation solve" open question in `CLAUDE.md`: tracing
  `w2_4_win.f90:896-1050` confirmed the free-surface solve is sequential across
  connected branches (branch JB's tridiagonal boundary term pulls `Z(UHS(JB))`/
  `Z(DHS(JB))` from a *different*, already-solved branch within the same timestep),
  not parallel like the TKE/momentum `TRIDIAG`. **Validated for generality, not just
  Detroit**: Detroit's branches happen to already be numbered downstream-to-upstream,
  so matching `[1,2,3,4]` alone wouldn't prove the algorithm is doing real topological
  work rather than just returning stored order -- also tested against a synthetic
  *reversed*-numbering topology (branch 4 is the true downstream/terminal branch;
  branches 1-3 all connect into it) and confirmed branch 4 is still correctly ordered
  first despite its highest number (see `test/runtests.jl`).
- **`Hydrodynamics/Density.jl`** — `density()`, the water equation of state (`density.f90`,
  a clean self-contained 27-line function, chosen as the next concrete step after tracing
  showed the free-surface+momentum solve itself isn't separable -- see `CLAUDE.md` "MVP
  hydrodynamic run"). Ported with EXPLICIT `fresh_water`/`salt_water`/`susp_solids`
  boolean arguments instead of the real Fortran's implicit `USE GLOBAL, ONLY:JW` shared-
  state read (Decision Log #4). **Validated** against known physical reference points, not
  just internal consistency: fresh water's density maximum at 4 degC (~999.97 kg/m^3),
  and values at 0/20/25 degC, all matching standard textbook figures to 3 decimal places.
  Getting the real `FRESH_WATER(JW)`/`SALT_WATER(JW)` flags right required tracing one
  level deeper than expected: they depend on `WTYPEC` (now read by `InputReader.jl`'s new
  initial-condition block) AND the global `CONSTITUENTS` flag (from `CCC`) AND
  `CAC(NTDS)` (is TDS itself active) -- both Tier 1, not read into structs yet. Rather
  than default those two to `true` (which would silently diverge from the real model for
  any control file with water quality off or TDS untracked), `Core/InitGeometry.jl`'s
  `compute_water_type_flags!` takes them as REQUIRED explicit keyword arguments, and is
  deliberately NOT called automatically by `init_geometry!` for that reason.
  Also (same investigation): `Core/InitGeometry.jl`'s previously-assumed-false
  `TRAPEZOIDAL(JW)` is now read for real (`GRIDC`, no Tier-1 dependency, unlike
  `FRESH_WATER`/`SALT_WATER`) and `compute_areas_volumes!` now ERRORS loudly (proven by a
  test that flips Detroit's own real, non-TRAP flag to force the guard to fire) instead of
  silently computing wrong RECT-formula areas for a TRAP-gridded waterbody.
- **`Hydrodynamics/FreeSurface.jl`** — traced directly from `w2_4_win.f90` (an
  `INTEGER FUNCTION`, not a separable `SUBROUTINE` like `init-geom.F90`), confirming the
  free-surface elevation solve is genuinely sequential across connected branches (see
  `Core/Grid.jl`'s `branch_processing_order` above). Confirmed with the user
  (2026-08-12) as **REDUCED PHYSICS BY DELIBERATE CHOICE**, not a faithful full port:
  REAL formulas for density (`RHO`), hydrostatic pressure (`P`), gravity/channel-slope
  (`GRAV`), the horizontal pressure gradient (`HPG`, though HDG/HPG are merged into one
  computation rather than kept as the two separate old/new-geometry passes the real
  source uses), the implicit free-surface tridiagonal solve itself (including real
  branch-to-branch sequencing and `DH_INTERNAL`/`UH_INTERNAL` boundary coupling), and
  the explicit velocity update. STUBBED AND FLAGGED (not silently omitted): bottom/wind
  shear (`SB`/`ST`, needs `Turbulence.jl` + meteorology IO), advection/dispersion of
  momentum (`ADMX`/`ADMZ`/`DM`, real formulas exist but are 0 whenever `U=0` anyway —
  port before trusting nonzero-inflow runs), the dam-flow/reciprocal-head-flow branch
  case (unreachable for Detroit), and the implicit vertical-eddy-viscosity correction
  step (needs `Turbulence.jl`). **Validated**: a real bug was found and fixed along the
  way — `input.F90:2265,2275`'s `H2(:,I) = H(:,JW)` baseline copy (every layer, not just
  the top active layer) was never ported into `BathymetryReader.jl`, leaving `geom.H2`
  zero below the top layer and causing `0.0/0.0` NaNs in the velocity update; fixed in
  `Core/InitGeometry.jl`'s `allocate_init_geometry!` (`geom.H2 = copy(geom.H)`, since
  `geom.H` is already populated by the time that function runs). This also revealed the
  earlier `VOL` test was too weak (finite-and-non-negative passes trivially for
  all-zero) — strengthened to assert nonzero volume below the top layer. Post-fix,
  `test/runtests.jl`'s zero-flow sanity check confirms water surface elevation and
  velocity stay bit-for-bit stable (`< 1e-9` drift) over 5 timesteps against real
  Detroit data — a real test of the tridiagonal assembly and branch sequencing (a bug
  in branch ordering or the `BHRHO`/`A`/`V`/`C` coefficients would show up as drift or
  NaN), not a tautology.
- **`Simulation.jl` / `IO/OutputWriter.jl`** — `run_zero_flow_sanity_check!` ties
  `init_geometry!` → `hydrodynamic_step!` → per-step TSR CSV output into one runnable
  driver. First-cut scope, matching `FreeSurface.jl`: ONE fixed `dlt = tc.DLTMAX[1]` for
  the whole run (the real adaptive DLTF/DLTMIN/DLTD breakpoint machinery isn't ported —
  fine for a reduced-physics zero-flow run, wrong once real forcing is added), no
  boundary-condition IO (Tier 1, not built), no kinetics/constituent transport.
  `OutputWriter.jl`'s TSR files follow the real per-segment naming convention
  (`outputinitw2tools.F90:984-1057`, `<base>_seg<N>.csv`) but only emit the columns this
  port actually computes (`JDAY,DLT(s),ELWS(m),U(ms-1)`) — `T2`/`Q`/`SRON`/`EXT`/etc are
  real TSR columns, deliberately absent rather than written as fabricated zeros.
  **Validated** end-to-end against real Detroit data (`test/runtests.jl`): runs 5 steps,
  writes real CSVs for one representative segment per branch, confirms stable ELWS and
  zero U in the written output itself (not just in-memory state).
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

**FIXED**: the periphyton/epiphyton-group block (`EPIC`/`EPRC`) had rows repeating the
same extracted key in an *interleaved* pattern (EPIC, EPRC, EPIC, EPRC, ...), which the
original "merge only *consecutive* same-key rows" logic didn't catch. Investigating (via
`review_config.jl` surfacing the warnings during a normal run, not a dedicated audit)
found this was one instance of a broader pattern, not a one-off: `SPILLWAYS`' upstream
and downstream SPECIFY-elevation fields are BOTH labeled `ETUSP`/`EBUSP` in the source
workbook (rows 201 and 206 -- the downstream one should logically be `ETDSP` to match
the `KTDSP`/`KBDSP` naming two rows later, but the workbook itself has the mislabel, not
just our extraction), and `GENERIC CONSTITUENT`'s `CGS` is reused for two unrelated
fields ("Settling rate" and "Gas transfer saturation concentration", with `CGLDK`/`CGKLF`
sitting between the two `CGS` rows). `emit_vertical_fields!` now groups by key across the
*whole* block (preserving first-occurrence order) instead of only merging adjacent rows
-- catches all three cases in one pass, confirmed by inspecting each one directly rather
than assumed. Re-running `review_config.jl` against the regenerated `w2_config.yaml`
shows **zero** duplicate-key warnings now (down from 3 known instances, ~30-something
affected rows total).

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
- **Outlet structures** (8896-8969): NSTR=0 with an external downstream boundary is a
  warning, NSTR>0 without one is an error, KTSTR/KBSTR ordering, SINKC must be LINE/
  POINT, WSTR>0 when SINKC=LINE.
- **Pipes** (8976-9119): WPI/DLXPI/FPI/FMINPI range checks, PUPIC/PDPIC/LATPIC/DYNPIPE
  enum validity, KTUPI/KTDPI vs KBUPI/KBDPI ordering (downstream fields only checked
  when IDPI != 0, matching the Fortran guard).
- **Withdrawals** (9931-9968): IWD must not sit on a branch boundary segment, KTWD>=2,
  KTWD<=KBWD.
- **Tributaries** (9976-10010): PTRC enum validity, ITR must fall within some branch's
  active segment range, ELTRT>=ELTRB when PTRC=SPECIFY.
- **CST constituent activation** (10383-10443): CAC/CPRWBC must be ON/OFF, initial
  concentration < -2 is an error, initial concentration = 0 while active is a warning
  (exempting age/residence-time tracers -- see the exemption-logic note in the file;
  this caught a real bug where the exemption checked only the short name and missed
  `Gen1`/"Water age, days", which the Fortran source's own equivalent check would exempt).

**Grid-elevation-dependent checks are deliberately skipped** -- a large fraction of the
structures/pipes/withdrawals/tributaries checks in the Fortran source compare against
`EL(K,I)` (computed layer elevation) and `KB(I)` (bottom active layer), both of which
come from the geometry solve in `init-geom.F90`, not yet ported to W2J (see CLAUDE.md
"Open questions" -- the `H(K,JW)` shape mismatch blocks this). Only checks that need
nothing beyond the raw control-file values are ported; each function's docstring in
`review_config.jl` says exactly which checks were left out and why.

Run with `julia --project=. review_config.jl [path-to-w2_config.yaml] [path-to-log]`.
Prints the report to stdout AND writes the identical text to a log file (default
`<yaml_path's directory>/Outputs/reviewed.log`, e.g. `DetroitReservoir/Outputs/
reviewed.log` -- matching the project's existing `Outputs/` convention from
`LongitudinalProfile.jl`) so the run is a durable, revisitable artifact, not just
terminal scrollback. Against the real Detroit `w2_config.yaml`: **0 errors, 1 warning**
(`ISS1` starts at zero concentration while active -- a legitimate, expected warning, not
a bug). Sanity-checked by mutating a copy of the config in-memory (bad `SINKC`,
`KTWD`>`KBWD`, out-of-grid `ITR`) and confirming all three fire correctly -- a clean run
isn't proof the checks aren't dead code.

Also fixed a real bug found while expanding scope: YAML `null` parses to Julia `nothing`,
and `string(nothing)` is the literal text `"nothing"` -- unguarded string checks (e.g. on
`PUPIC`/`DYNPIPE` for the `PIPES` section, all null when `NPI=0`) were treating that as
real, invalid data and raising false-positive errors. Fixed via a single `str_or_empty`
helper used by every string-valued check now, rather than patching each site ad hoc.

**A clean run of this script does not mean the control file is fully valid** -- roughly
1,100 more checks from the reference preprocessor (gates, meteorology, kinetics rate
coefficients, output control, and the grid-elevation-dependent checks noted above) are
not ported yet. Each check in `review_config.jl` cites its source line range in
`pre_ivf_422.f90` so it can be independently verified or used as a starting point to add
more.

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
- [x] `thomas_solve!` validated against a known tridiagonal test case: cross-checked
      against a literal transcription of the real `TRIDIAG` subroutine
      (transport.f90:572-593) and Julia's dense solver, 200 randomized cases + 1 fixed
      case (see `test/runtests.jl`)
- [x] `init-geom.F90` ported (`Core/InitGeometry.jl`): `EL`/`KTI`/`KTWB`/`KB`/`VOL`/
      `DEPTHB`/`DEPTHM`/areas, plus the `init.F90` boundary-condition flags and
      `ALPHA`/`SINA`/`COSA` it depends on. Validated against real Detroit data. Resolves
      the `H(K,JW)` shape question below. Known deferrals (RESTART_IN, trapezoidal
      grids, `w2_constriction.csv`, JBTR/JBWD, and an **approximated** cross-branch
      boundary-width interpolation) are documented in the file's module docstring --
      see "What's Actually Implemented" above.
- [x] `Core/Grid.jl` branch connectivity: `BranchNetwork` + `build_branch_network`,
      validated against real Detroit topology, and `Core/InitGeometry.jl` refactored to
      use it instead of its own duplicated lookup. Scope: static topology only --
      `Hydrodynamics/Waterbody.jl` (the per-timestep boundary-value interpolation) is
      still a stub, see "What's Actually Implemented" above.
- [x] Free-surface elevation solve sequential-dependency question: RESOLVED, confirmed
      sequential (not parallel like TKE/momentum). `Core/Grid.jl`'s
      `branch_processing_order` computes the correct per-waterbody topological order for
      any control file (validated against a synthetic reversed-numbering topology, not
      just Detroit's already-convenient numbering) -- see "What's Actually Implemented".
      The free-surface solve itself is not yet implemented (see Progress Tracker /
      CLAUDE.md "MVP hydrodynamic run" for why it's not a single-file port like
      `init-geom.F90` was, and the recommended decomposition).
- [x] `Hydrodynamics/Density.jl` (`density()`, water equation of state) implemented +
      validated against known physical reference values (4degC density maximum, 0/20/25
      degC points). `InputReader.jl` extended to read the initial-condition block
      (`T2I`/`ICEI`/`WTYPEC`/`GRIDC`); `Core/InitGeometry.jl`'s `TRAPEZOIDAL(JW)` now
      read for real instead of assumed false, and `compute_areas_volumes!` errors loudly
      (test-proven) rather than silently computing wrong areas for a TRAP grid.
      `FRESH_WATER`/`SALT_WATER` require explicit `constituents`/`tds_active` args
      (their real Fortran dependencies, not yet read into structs) rather than being
      defaulted to `true` -- see "What's Actually Implemented" above.
- [x] Free-surface + momentum solve implemented (`Hydrodynamics/FreeSurface.jl`),
      REDUCED PHYSICS by deliberate choice (confirmed with user, 2026-08-12): real
      density/pressure/gravity/pressure-gradient/tridiagonal-solve/branch-sequencing/
      velocity-update; shear/advection/dispersion/dam-flow/vertical-eddy-viscosity
      explicitly zeroed and flagged (need `Turbulence.jl` + meteorology IO, not built).
      Found + fixed a real bug along the way: `input.F90:2265,2275`'s `H2(:,I)=H(:,JW)`
      baseline copy was never ported, causing `0.0/0.0` NaNs in the velocity update --
      fixed in `Core/InitGeometry.jl`. See "What's Actually Implemented" above.
- [x] `Simulation.jl` first-cut driver (`run_zero_flow_sanity_check!`) + `IO/
      OutputWriter.jl` first-cut TSR CSV writer implemented and validated end-to-end
      against real Detroit data: 5-step zero-boundary-flow run produces real per-segment
      TSR CSVs showing a stable water surface (`< 1e-9` drift) and zero velocity. Fixed
      `dlt`, no boundary-condition IO, no kinetics -- see "What's Actually Implemented".
      **This is the "first cut running and getting TSR outputs for Detroit" milestone
      the user asked for (2026-08-13).** Parallel processing is the next agreed phase,
      not yet started.
- [ ] `LAM1` derivation traced; `RateMultipliers.jl` implemented for real
- [ ] Parallel processing (`Threads.@threads` over segments/cells) -- NEXT STEP, per
      user's explicit "then next is to prepare for parallel processing" (2026-08-13)
- [ ] Real adaptive timestep (DLTF/DLTMIN/DLTD breakpoints), boundary-condition IO
      (Tier 1: inflow/outflow/withdrawal/tributary time series), and non-zero-flow
      forcing (ADMX/ADMZ/DM, SB/ST via `Turbulence.jl` + meteorology) -- needed before
      `Simulation.jl`/`FreeSurface.jl` can run anything beyond the zero-flow sanity check
- [ ] First Fortran reference-output comparison
- [ ] Cross-group kinetics coupling investigation resumed (paused mid-trace — see
      project research notes Open Questions, and the note in `Kinetics.jl`)
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
