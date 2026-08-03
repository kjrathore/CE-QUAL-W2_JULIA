"""
Water quality kinetics — the ~70 ENTRY decomposition.

Originating Fortran context: water-quality.f90's KINETICS subroutine.
SINGLE HIGHEST-RISK FILE in the codebase (see README Decision Log #4,
project research notes). ~70 ENTRY points sharing implicit SAVE'd state;
Julia has no equivalent construct, so each entry becomes an explicit
function, grouped here by physical category (matching the constituent-
workflow diagram already produced):

    Nutrients  — AMMONIUM, NITRATE, PHOSPHORUS, ...
    Organics   — LABILE_DOM, REFRACTORY_DOM, LABILE_POM, REFRACTORY_POM, ...
    Algae      — ALGAE, EPIPHYTON, ZOOPLANKTON, ...
    Sediment   — SEDIMENT, SEDIMENT1, SEDIMENT2, SEDIMENTP/N/C, ...

Each group lives in its own file under Constituents/ (see includes below).

OPEN QUESTION (not yet resolved — see README Open Questions): whether
Kinetic_rates / the mass-balance assembly step has genuine same-pass
cross-group coupling (e.g. DO consuming terms from Algae, Nutrients, and
Sediment groups simultaneously) or whether groups are fully independent
per timestep. Investigation paused mid-trace (was checking O2 update
lines) — pick back up before finalizing whether these Constituent
functions can run in any order / in parallel with each other, or whether
there's a required sequencing beyond "RateMultipliers first."

Design intent (Pillar 3): every function below must accept local cell
state and return rate terms only — no (K,I) indexing assumptions baked
into the kinetics math itself. Spatial layout is the caller's concern.
"""

include("Constituents/Nutrients.jl")
include("Constituents/Organics.jl")
include("Constituents/Algae.jl")
include("Constituents/Sediment.jl")

# function kinetic_rates(nutrients, organics, algae, sediment, rate_mult::RateMultipliers)
#     # TODO: this is the aggregation/join point flagged in the open question
#     # above. Do not implement until cross-group coupling is confirmed.
# end
