"""
Relativistic Riemann Solver

A Python package for solving relativistic Riemann problems in special 
relativistic hydrodynamics. Based on Marti and Mueller, J. Fluid Mech., (1994).
"""

from .solver import RelativisiticRiemannSolver
from .test_cases import (
    test_case_sr_sod,
    test_case_blast_wave,
    test_case_relativistic_shock,
    test_case_two_shocks,
)

__version__ = "0.1.0"
__all__ = [
    "RelativisiticRiemannSolver",
    "test_case_sr_sod",
    "test_case_blast_wave",
    "test_case_relativistic_shock",
    "test_case_two_shocks",
]
