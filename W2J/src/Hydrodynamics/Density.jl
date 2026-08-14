# ==============================================================================
# Hydrodynamics/Density.jl
#
# Ports DENSITY (w2_source_5-21-2026/density.f90, 27 lines), the water-density
# equation of state used throughout the hydrodynamic core (the free-surface
# solve's BHRHO term at w2_4_win.f90:904, buoyancy/momentum terms, RHO(K,I)
# feeding GRAV/ADMZ) and water quality (dissolved oxygen saturation, etc).
# Chosen as the next concrete implementation target after Core/Grid.jl's
# `branch_processing_order`, per the recommended decomposition of "step 4"
# (free-surface + momentum solve) in CLAUDE.md "MVP hydrodynamic run" -- a
# clean, isolated, well-defined formula, unlike the free-surface solve
# itself, which is not separable (see that file's tracing notes).
#
# EXPLICIT STATE, NOT IMPLICIT GLOBAL (Decision Log #4): the real Fortran
# function reads FRESH_WATER(JW)/SALT_WATER(JW)/SUSP_SOLIDS via
# `USE GLOBAL, ONLY:JW` -- JW is set by the CALLER as a shared module
# variable before invoking DENSITY, the same ENTRY-style implicit-state
# pattern flagged project-wide as something to NOT port literally. Ported
# here as ordinary explicit boolean arguments instead.
# ==============================================================================

"""
    density(T, TDS, SS, fresh_water::Bool, salt_water::Bool, susp_solids::Bool)

density.f90:5-27 -- water density (kg/m^3) as a function of temperature `T`
(deg C), total dissolved solids `TDS` (g/m^3), and suspended solids `SS`
(g/m^3), with separate polynomial fits below and above 0 deg C.

`fresh_water`/`salt_water`/`susp_solids` are explicit booleans, not derived
from any global state here -- see `Core/InitGeometry.jl`'s
`compute_water_type_flags!` for how `fresh_water`/`salt_water` are actually
derived (and why it's NOT simply "WTYPEC==FRESH" -- there's a real, easy to
silently drop, dependency on whether constituents/TDS tracking are even
active). `susp_solids` mirrors Fortran's `SUSP_SOLIDS` global, itself just
`NSS > 0` (call sites should pass `g.NSS > 0`).

VALIDATED (test/runtests.jl, "Hydrodynamics/Density: density"): matches
known physical reference points -- fresh water's density maximum at 4 deg C
(~999.97 kg/m^3), density at 20 deg C (~998.2 kg/m^3), and 0 deg C
(~999.84 kg/m^3) -- not just internal self-consistency.
"""
function density(T, TDS, SS, fresh_water::Bool, salt_water::Bool, susp_solids::Bool)
    if T <= 0.0
        rho = 0.842594
        susp_solids && (rho += 6.2e-4 * SS)
        fresh_water && (rho += TDS * 8.221e-4)
        salt_water && (rho += TDS * 0.824493 + (-5.72466e-3) * TDS^1.5 + 4.8314e-4 * TDS * TDS)
    else
        rho = ((((6.536332e-9 * T - 1.120083e-6) * T + 1.001685e-4) * T - 9.09529e-3) * T + 6.793952e-2) * T + 0.842594
        susp_solids && (rho += 6.2e-4 * SS)
        fresh_water && (rho += TDS * ((4.99e-8 * T - 3.87e-6) * T + 8.221e-4))
        salt_water && (rho += TDS * ((((5.3875e-9 * T - 8.2467e-7) * T + 7.6438e-5) * T - 4.0899e-3) * T + 0.824493) +
                               ((-1.6546e-6 * T + 1.0227e-4) * T - 5.72466e-3) * TDS^1.5 + 4.8314e-4 * TDS * TDS)
    end
    return rho + 999.0
end
