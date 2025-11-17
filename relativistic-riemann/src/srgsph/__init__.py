"""
SRGSPH: Special Relativistic Godunov Smoothed Particle Hydrodynamics

Implementation based on:
Kitajima, K., Inutsuka, S., & Seno, I. (2025)
"Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver"
arXiv:2510.18251v1
"""

__version__ = "0.1.0"

from .particle import Particle
from .kernel import GaussianKernel
from .simulator import SRGSPH1D

__all__ = ["Particle", "GaussianKernel", "SRGSPH1D"]
