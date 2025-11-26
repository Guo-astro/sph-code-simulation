"""
BNS Ejecta Models Package

Provides tools for constructing initial conditions for BNS merger ejecta simulations:
- Loading Radice+2018 GR simulation data
- Analytic fits (Hotokezaka, Bauswein, etc.)
- Homologous expansion mapping
- SRGSPH particle sampling
"""

from .radice_loader import RadiceDataLoader
from .homologous_expansion import HomologousExpander
from .analytic_profiles import (
    HotokezakaProfile,
    BausweinProfile,
    PowerLawProfile,
)
from .particle_sampler import ParticleSampler

__all__ = [
    'RadiceDataLoader',
    'HomologousExpander',
    'HotokezakaProfile',
    'BausweinProfile',
    'PowerLawProfile',
    'ParticleSampler',
]

__version__ = '0.1.0'
