"""
Nutrient kinetics — NH4, NO3, PO4 and related entries.

Originating ENTRY points (water-quality.f90): AMMONIUM, NITRATE, PHOSPHORUS,
DISSOLVED_SILICA, PARTICULATE_SILICA. Confirmed read-only consumers of
RateMultipliers (NH4TRM, NO3TRM via the trm_usage.txt trace — e.g. line 316:
NH4D(K,I) = NH4TRM(K,I)*NH4DK(JW)*NH4(K,I)*DO1(K,I)).

TODO: not yet ported. Port after RateMultipliers.jl's LAM1 question is
resolved.
"""

# function ammonium(nh4, rate_mult::RateMultipliers, params) end
# function nitrate(no3, rate_mult::RateMultipliers, params) end
