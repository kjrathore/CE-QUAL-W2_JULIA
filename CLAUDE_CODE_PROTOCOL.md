# Claude Code Operating Protocol for W2J

Instructions for working efficiently with Claude Code on this project.
Read this once, then keep it open as a reference.

---

## Setup (do this once)

1. Place `CLAUDE.md` in the root of your W2J repo.
   Claude Code reads this automatically at session start — it's the context transfer.

2. Update `W2J_README.md` progress tracker to reflect current state before your first
   Claude Code session (the one in the project is stale — the IO files built in the
   Claude.ai session are not reflected there yet).

3. Confirm your Julia environment has `Plots.jl` installed:
   ```
   julia --project=. -e "using Pkg; Pkg.add(\"Plots\")"
   ```

---

## How to start each session

Open Claude Code and give it one of these openers depending on what you're doing:

**Starting a new feature:**
> "Read CLAUDE.md, then read `src/IO/InputReader.jl`. I want to implement [specific thing]. The relevant Fortran source is `w2_source/input.F90`. Read that file next and then propose an approach before writing any code."

**Fixing a bug from a failed run:**
> "Read CLAUDE.md. I ran [command] and got this error: [paste error]. Read only the specific file the error points to."

**Continuing from last session:**
> "Read CLAUDE.md and `W2J_README.md`. What's the next unchecked item in the progress tracker?"

**Never start with:**
> "Read all the source files and tell me what to do next."
That's the token disaster pattern.

---

## Token-efficient patterns

### The one-file rule
Claude Code should read **one Fortran file per Julia function being written**, not the whole source tree. The call-graph and ENTRY map are already in `CLAUDE.md` — there's no need to re-derive them.

Tell Claude Code explicitly:
> "Do not read any Fortran file I haven't asked you to read."

### Validate before implementing
Before writing a new reader or parser, ask Claude Code to write a
10-line validation script first:
> "Before writing the Julia function, write a quick Julia/Python script
> that reads [filename] and prints the first few parsed values. Run it.
> If the output matches the validated values in CLAUDE.md, then write the function."

This costs one small script run instead of an implement-run-fix-run cycle
that might take 5 turns.

### Batching related fields
When implementing a new section of `InputReader.jl` (e.g. the initial
conditions block), ask Claude Code to trace the relevant portion of
`input.F90` first, write out the expected field order as a comment block,
get your confirmation, then write the code. One read of the Fortran file,
one code write, done.

### Stopping before stubs
If Claude Code is about to write a stub for something that won't be needed
for several weeks, stop it:
> "Don't write the stub. Just add a TODO comment and move on."
Stubs cost tokens to write and tokens to re-read later; a comment costs nothing.

---

## What to do when Claude Code reads too much

If you see it starting to read multiple Fortran files or the entire `src/`
directory, interrupt:
> "Stop. Read only [specific file] and nothing else. Tell me what you found
> in that file before reading anything else."

Claude Code respects explicit scope constraints when stated upfront. The
issue is usually not that it ignores them — it's that they weren't stated.

---

## Session-ending checklist

Before closing Claude Code:
1. Run the test suite: `julia --project=. -e "using Pkg; Pkg.test()"`
2. If any new file moved from stub → implemented, update the progress
   tracker in `W2J_README.md`.
3. If any new architectural decision was made, add it to the decision log
   in `README.md` and the relevant section of `CLAUDE.md`.
4. Commit. Each session should produce one clean commit with a message
   like `IO: implement Phase B initial-conditions block, validated against Detroit`.

---

## The right order to build things

From `CLAUDE.md`'s immediate next task, working forward:

1. **This week**: First end-to-end run — `read_control_file` → `read_bathymetry!` → `plot_longitudinal_profile`. Fix whatever breaks. The allocation gap for `W2Geometry` arrays (DLX, ELWS, B not yet sized) will be the first error.

2. **Then**: Phase B of `InputReader.jl` — trace `input.F90` lines 762–840 (initial conditions, boundary conditions, meteorology, ice, turbulence closure). Same assertion-based approach as Phase A.

3. **Then**: `init-geom.F90` → `Core/InitGeom.jl` — this resolves the `H(K,JW)` shape question and gives real `KB`/`EL` for the longitudinal profile plot.

4. **Then**: `thomas_solve!` validation against a known tridiagonal test case — write a 3×3 system by hand, check the output.

5. **Then**: Phase C of `InputReader.jl` — constituent block (the `NCT` count that's currently left at 0).

Do not skip ahead to kinetics or transport until the geometry/IO stack produces
correct numbers for Detroit.

---

## Quick reference: Detroit validated values

If any of these change after a parser modification, something is wrong:
```
NWB=1  NBR=4  IMX=31  KMX=117  NPROC=1
US=[2,14,22,28]  DS=[11,19,25,30]
BS=1   BE=4   LAT=45.7299   LONGIT=122.177   ELBOT=364
TMSTRT=1  TMEND=365  YEAR=2002  DLTMAX=1200
bth1: DLX[2]=844.9  ELWS[1]=441.09  B[60,2]=84.07
```
