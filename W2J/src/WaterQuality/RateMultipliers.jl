"""
Temperature-dependent rate multipliers and oxygen-limitation switches.

Originating Fortran context: water-quality.f90's TEMPERATURE_RATES entry
(writes) and the DO1/DO2/DO3 computation immediately following it. This is
the MOST THOROUGHLY CONFIRMED piece of the whole port — traced line-by-line
against the actual source (lines 230-252, 294-296) rather than inferred from
naming. Single-writer / many-reader, independent per (K,I) cell — see
project research notes, "trm_usage.txt" trace and the constituent-workflow
diagram. This is also the namesake / motivating example for Core/State.jl's
RateMultipliers struct.

Confirmed exact formulas (direct transcription from source):
    NH4TRM(K,I) = LAM1/(1.0+LAM1-NH4K1(JW))
    NO3TRM(K,I) = LAM1/(1.0+LAM1-NO3K1(JW))
    OMTRM(K,I)  = LAM1/(1.0+LAM1-OMK1(JW))
    SODTRM(K,I) = LAM1/(1.0+LAM1-SODK1(JW))
    ATRM(K,I,JA) = ATRMR(K,I,JA) * ATRMF(K,I,JA)
    ETRM(K,I,JE) = ETRMR(K,I,JE) * ETRMF(K,I,JE)
    DO1(K,I) = O2(K,I) / (O2(K,I) + KDO)
    DO2(K,I) = 1.0 - DO1(K,I)
    DO3(K,I) = (1.0 + SIGN(1.0, O2(K,I) - 1.0e-10)) * 0.5

NOT YET CONFIRMED — flagged honestly rather than guessed: where LAM1 itself
comes from. It's used here but its own derivation (presumably an Arrhenius-
style temperature term computed earlier in the same host subroutine) was
never traced. TODO before implementing: grep water-quality.f90 for `LAM1 =`
to find its defining expression before porting this function for real.

Note the DO3 switch's SIGN-based formulation is deliberately smooth/
algebraic rather than a hard IF/THEN branch — exactly the pattern Pillar 2
(differentiability) wants more of. Worth preferring this style consistently
when porting the remaining ~70 ENTRY points in Kinetics.jl.
"""

"""
    rate_multipliers!(out::RateMultipliers, T, O2, rate_constants, K_range, I_range)

Compute all rate multipliers for the given cells. Independent per (K,I) —
safe to call inside a `Threads.@threads` loop over I, same pattern as
Solvers/Tridiagonal.jl's batch solve.

TODO: resolve LAM1's derivation (see module docstring) before implementing.
Placeholder signature only.
"""
function rate_multipliers!(out, T, O2, rate_constants, K_range, I_range)
    error("Not yet implemented — LAM1 derivation must be traced first (see TODO in module docstring).")
end
