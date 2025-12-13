"""
Analysis modules for BNS cocoon shock breakout simulations.

Modules:
- shock_tracker: Track shock position and Lorentz factor
- optical_depth: Compute optical depth profiles
- observables: Calculate observable quantities (E_rad, T_obs, L_peak)
"""

from .shock_tracker import ShockTracker
from .optical_depth import OpticalDepthCalculator
from .observables import BreakoutObservables

__all__ = ['ShockTracker', 'OpticalDepthCalculator', 'BreakoutObservables']
