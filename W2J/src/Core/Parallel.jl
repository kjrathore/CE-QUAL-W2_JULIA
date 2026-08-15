# ==============================================================================
# Core/Parallel.jl
#
# Threading utility used throughout Hydrodynamics/FreeSurface.jl (and future
# per-cell modules -- Turbulence.jl, Transport.jl, Kinetics.jl, all confirmed
# embarrassingly parallel per Pillar 1's `TRIDIAG` case-A analysis).
#
# THE PROBLEM THIS SOLVES: naively wrapping every per-column loop in
# `Threads.@threads` (the first-pass parallel-processing implementation,
# 2026-08-14) made Detroit's zero-flow sanity check SLOWER, not faster --
# measured 0.136ms/step at 1 thread vs 0.484ms/step at 8 threads. Detroit's
# grid is small (IMX=31, at most ~31 columns in any one loop), so
# `Threads.@threads`'s per-call thread-spawn/scheduling overhead exceeds the
# actual per-column work (a few dozen flops). This is not a Detroit-only
# concern -- W2J is explicitly meant to generalize to other reservoirs
# (CLAUDE.md "GENERALITY REMINDER"), some much larger than Detroit, where
# the same per-column loops WILL be worth parallelizing. `parallel_foreach`
# picks serial or threaded execution PER CALL based on the actual loop
# length, so both cases are handled by the same code path without the
# caller needing to know in advance which one it will get.
# ==============================================================================

"""
    PARALLEL_THRESHOLD

Minimum loop length (segment/column/branch count) at which `parallel_foreach`
switches from a plain serial loop to `Threads.@threads`. Empirically
determined (not guessed) by benchmarking a synthetic per-column workload
shaped like `Hydrodynamics/FreeSurface.jl`'s heaviest per-column loop
(`update_velocities!`'s ~100-layer inner-K loop) at `--threads=4`: serial
beats threaded up to N~64-96, threaded starts winning around N=128 (1.3x
measured) and keeps improving (3.6x measured by N=4096). This is
machine/OS/Julia-version dependent -- a different crossover point on
different hardware is expected, not a bug -- but 128 is a reasonable,
conservative default (clearly past the measured crossover, not sitting at
the noisy breakeven point).
"""
const PARALLEL_THRESHOLD = 128

"""
    parallel_foreach(f, range; threshold=PARALLEL_THRESHOLD)

Calls `f(i)` for every `i in range`, using `Threads.@threads` only when
`length(range) >= threshold` AND more than one thread is actually available
(`Threads.nthreads() > 1` -- skips the `Threads.@threads` machinery's own
overhead entirely on a single-threaded Julia process, where it can only
ever add cost). Below that, a plain `for` loop -- see `PARALLEL_THRESHOLD`'s
docstring for why this isn't just "always thread."

`threshold` defaults to `PARALLEL_THRESHOLD` and exists as an override
purely so `test/runtests.jl` can exercise both the serial and threaded code
paths deterministically regardless of how many segments the test's control
file happens to have, rather than depending on incidental range lengths.

Every call site using this must still satisfy the same no-data-race
requirement `Threads.@threads` itself required (each `i` writes only its
own data) -- this function changes WHEN to thread, not what's safe to.
"""
function parallel_foreach(f, range; threshold::Int=PARALLEL_THRESHOLD)
    if Threads.nthreads() > 1 && length(range) >= threshold
        Threads.@threads for i in range
            f(i)
        end
    else
        for i in range
            f(i)
        end
    end
    return nothing
end
