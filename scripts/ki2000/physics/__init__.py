"""
K&I 2000 Physics Package

First-principles heating and cooling rates for ISM thermal equilibrium,
based on Koyama & Inutsuka (2000, ApJ 532, 980).

Modules:
    heating - Heating processes (photoelectric, cosmic ray, X-ray, H2)
    cooling - Cooling processes (CII, OI, Ly-alpha, H2, CO, gas-grain)
    equilibrium - Thermal equilibrium solver
"""

from .heating import (
    photoelectric_heating,
    cosmic_ray_heating,
    xray_heating,
    h2_formation_heating,
    total_heating_rate,
    heating_breakdown
)

from .cooling import (
    cii_cooling,
    oi_cooling,
    lyman_alpha_cooling,
    h2_cooling,
    co_cooling,
    gas_grain_cooling,
    total_cooling_rate,
    cooling_breakdown
)

from .equilibrium import (
    ThermalEquilibrium,
    ChemicalEquilibrium,
    reproduce_ki2000_panel_c
)

__all__ = [
    # Heating
    'photoelectric_heating',
    'cosmic_ray_heating', 
    'xray_heating',
    'h2_formation_heating',
    'total_heating_rate',
    'heating_breakdown',
    # Cooling
    'cii_cooling',
    'oi_cooling',
    'lyman_alpha_cooling',
    'h2_cooling',
    'co_cooling',
    'gas_grain_cooling',
    'total_cooling_rate',
    'cooling_breakdown',
    # Equilibrium
    'ThermalEquilibrium',
    'ChemicalEquilibrium',
    'reproduce_ki2000_panel_c'
]
