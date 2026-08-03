"""
Architecture abstraction — determines where kernels execute (CPU threads vs GPU).

Design note: this pattern is deliberately modeled on Oceananigans.jl's
Architectures.jl (CPU / GPU{D} types passed explicitly into solver functions)
even though W2J does not depend on Oceananigans (see Decision Log #1 in the
project research notes). Keeping this as an explicit, swappable argument now
is what lets Solvers/Tridiagonal.jl (and later kernels) run on CPU threads
today and GPU later without restructuring call sites.

Originating Fortran context: N/A — no equivalent concept exists in the
original codebase, which is single-threaded throughout (Pillar 1 driver).
"""

abstract type AbstractArchitecture end

"Run kernels using `Threads.@threads` across available CPU threads."
struct CPU <: AbstractArchitecture end

"Run kernels on GPU. Placeholder — not implemented until CPU threading path
is validated against Fortran reference output (see README Progress Tracker)."
struct GPU <: AbstractArchitecture end

# TODO: once KernelAbstractions.jl is added as a dependency, replace the
# Threads.@threads calls in Solvers/Tridiagonal.jl with @kernel functions
# dispatched on `arch::AbstractArchitecture`, matching this pattern:
#   loop_kernel!(arch::CPU, ...) = Threads.@threads for i in ... ; ... ; end
#   loop_kernel!(arch::GPU, ...) = <KernelAbstractions @kernel launch>
