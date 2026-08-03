"""
Sediment oxygen demand and sediment decay kinetics.

Originating ENTRY points (water-quality.f90): SEDIMENT, SEDIMENT1,
SEDIMENT2, SEDIMENTP, SEDIMENTN, SEDIMENTC. Confirmed read-only consumers
of SODTRM (e.g. lines 297-309 in trm_usage.txt trace).

Distinct from the larger CEMA sediment DIAGENESIS modules (Diagenesis
Sediment Model 03.f90, Diagenesis Sediment Flux Model 05.f90, etc.) —
those are a separate, more detailed optional submodel and a later Tier 0
porting target in their own right, not part of this file's scope.

TODO: not yet ported.
"""
