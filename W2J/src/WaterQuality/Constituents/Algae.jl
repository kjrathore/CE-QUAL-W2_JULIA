"""
Algae and epiphyton kinetics.

Originating ENTRY points (water-quality.f90): ALGAE, EPIPHYTON,
ZOOPLANKTON, MACROPHYTE. Confirmed read-only consumers of ATRM/ETRM
(e.g. lines 686-689: AGR/ARR/AER computed from ATRM(K,I,JA)).

NOTE: this is one of the two most likely candidates for genuine cross-
group coupling with Nutrients/Sediment (algal growth consumes nutrient
pools) — relevant to the open question flagged in Kinetics.jl. Check
here first when resuming that investigation.

TODO: not yet ported.
"""
