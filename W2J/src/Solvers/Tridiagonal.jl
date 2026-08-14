"""
Shared tridiagonal (Thomas algorithm) solver.

Originating Fortran context: TRIDIAG, called from az.f90 (4 call sites,
TKE turbulence closure) and w2_4_win.f90 (1 call site, implicit momentum
solve). CONFIRMED Case A in both: each call operates on a single segment's
own vertical column only (sliced at fixed I, spanning KT to KB(I)/KBMIN(I)),
called once per segment inside a `DO I` loop — see project research notes,
"Checking the bottleneck" investigation. No algorithm redesign needed
(no cyclic reduction / Spike / PCR required) — confirmed independent
per-segment systems, so the win is parallelizing the *batch*, not the
algorithm itself.

Ported ONCE here (Decision Log #5) and reused by every hydrodynamics
module that needs it (Turbulence.jl, Transport.jl, and the momentum solve
in whatever replaces w2_4_win.f90's main loop) — rather than reimplemented
per call site as in the original Fortran.
"""

"""
    thomas_solve!(a, b, c, d, x)

Solve a single tridiagonal system Ax = d using the Thomas algorithm.
`a` = sub-diagonal, `b` = diagonal, `c` = super-diagonal, `d` = RHS.
Operates on one column — this function itself stays serial (correct,
since a single Thomas elimination has a genuine sequential dependency);
the parallelism comes from calling this independently across many
columns at once (see `solve_vertical_systems!` below).

VALIDATED (test/runtests.jl, "Solvers/Tridiagonal: thomas_solve!"): cross-
checked against Julia's dense `\\` solver and against a literal transcription
of `TRIDIAG` (transport.f90:572-593) over 200 randomized diagonally-dominant
systems, plus one hand-verified fixed case. The two recurrences are
algebraically equivalent (this one normalizes incrementally; TRIDIAG keeps
values unnormalized until the final back-substitution division) — confirmed
by direct comparison, not just by inspection.
"""
function thomas_solve!(a::AbstractVector, b::AbstractVector, c::AbstractVector,
                        d::AbstractVector, x::AbstractVector)
    n = length(b)
    # Forward elimination
    cp = similar(b); dp = similar(d)
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]
    @inbounds for i in 2:n
        m = b[i] - a[i] * cp[i-1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i-1]) / m
    end
    # Back substitution
    x[n] = dp[n]
    @inbounds for i in (n-1):-1:1
        x[i] = dp[i] - cp[i] * x[i+1]
    end
    return x
end

"""
    solve_vertical_systems!(arch, U, AT, VT, CT, DT, KT, KBMIN, segments)

Batch entry point — solves one independent tridiagonal system per segment
in `segments`, dispatched on `arch` (see Core/Architecture.jl). This is
the function that actually replaces the `DO I ... CALL TRIDIAG(...) ...
END DO` pattern confirmed in both az.f90 and w2_4_win.f90.

TODO: implement GPU dispatch path once CPU path is validated (likely via
CUSPARSE's batched solver — gtsv2StridedBatch — rather than a hand-rolled
KernelAbstractions kernel, given this maps directly onto "many small
independent tridiagonal systems", the documented use case for that
routine — see project research notes).
"""
function solve_vertical_systems!(::CPU, U, AT, VT, CT, DT, KT, KBMIN, segments)
    Threads.@threads for i in segments
        kb = KBMIN[i]
        @views thomas_solve!(AT[KT:kb, i], VT[KT:kb, i], CT[KT:kb, i], DT[KT:kb, i], U[KT:kb, i])
    end
    return U
end

function solve_vertical_systems!(::GPU, args...)
    error("GPU path not yet implemented — see TODO above. Validate CPU path first.")
end
