#!/usr/bin/env python3
"""
Optical Depth Calculator for Shock Breakout

Computes optical depth τ(r) from a given density profile:
    τ(r) = ∫_r^∞ κ ρ(r') dr'

Breakout occurs when τ = 1 (or c/v for relativistic case).
"""

import argparse
from pathlib import Path
from typing import Tuple, Optional
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class BreakoutState:
    """State at shock breakout."""
    time: float           # Breakout time
    r_breakout: float     # Breakout radius
    v_breakout: float     # Breakout velocity [c]
    gamma_breakout: float # Breakout Lorentz factor
    tau_at_shock: float   # Optical depth at shock
    
    # Shell properties
    internal_energy: float  # Internal energy in shocked shell
    shell_volume: float     # Volume of shocked shell
    shell_mass: float       # Mass of shocked shell


class OpticalDepthCalculator:
    """
    Calculate optical depth and detect shock breakout.
    
    Usage:
        calc = OpticalDepthCalculator(opacity=0.2)
        tau, r = calc.compute_optical_depth(r, rho)
        breakout = calc.find_breakout(r, rho, v_r, r_shock)
    """
    
    def __init__(self, opacity: float = 0.2, relativistic: bool = True):
        """
        Initialize calculator.
        
        Args:
            opacity: Opacity κ [cm²/g]
            relativistic: Use relativistic breakout condition τ = c/v
        """
        self.opacity = opacity
        self.relativistic = relativistic
    
    def compute_optical_depth(self, r: np.ndarray, 
                              rho: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute optical depth from infinity inward.
        
        τ(r) = ∫_r^∞ κ ρ(r') dr'
        
        Args:
            r: Radial positions (increasing)
            rho: Density at each radius
            
        Returns:
            tau: Optical depth at each radius
            r: Corresponding radii
        """
        # Integrate from outside in
        # τ(r) = κ ∫_r^∞ ρ dr' ≈ κ Σ_{i>r_idx} ρ_i Δr_i
        
        dr = np.gradient(r)
        integrand = self.opacity * rho * dr
        
        # Cumulative sum from the end
        tau = np.cumsum(integrand[::-1])[::-1]
        
        return tau, r
    
    def find_photosphere(self, r: np.ndarray, tau: np.ndarray) -> Optional[float]:
        """
        Find photosphere radius where τ = 1.
        
        Args:
            r: Radial positions
            tau: Optical depth at each radius
            
        Returns:
            Photosphere radius or None if τ < 1 everywhere
        """
        if np.all(tau < 1):
            return None
        
        # Linear interpolation
        idx = np.searchsorted(tau[::-1], 1.0)
        idx = len(tau) - 1 - idx
        
        if idx <= 0 or idx >= len(r) - 1:
            return r[idx]
        
        # Interpolate
        r_phot = np.interp(1.0, [tau[idx+1], tau[idx]], [r[idx+1], r[idx]])
        
        return r_phot
    
    def compute_breakout_tau(self, v_shock: float) -> float:
        """
        Compute optical depth at breakout for relativistic shock.
        
        For relativistic shocks, breakout occurs at τ ≈ c/v.
        For non-relativistic, τ ≈ 1.
        
        Args:
            v_shock: Shock velocity [c]
            
        Returns:
            Critical optical depth for breakout
        """
        if self.relativistic and v_shock > 0.01:
            return 1.0 / v_shock  # τ_bo = c/v
        else:
            return 1.0
    
    def is_breakout(self, tau_at_shock: float, v_shock: float) -> bool:
        """Check if breakout condition is satisfied."""
        tau_crit = self.compute_breakout_tau(v_shock)
        return tau_at_shock <= tau_crit
    
    def find_breakout_from_history(self, shock_history: pd.DataFrame,
                                   density_profiles: dict) -> Optional[BreakoutState]:
        """
        Find breakout from shock history and density profiles.
        
        Args:
            shock_history: DataFrame with columns [time, r_shock, v_shock, ...]
            density_profiles: Dict of {time: (r, rho)} arrays
            
        Returns:
            BreakoutState at first breakout or None
        """
        for _, row in shock_history.iterrows():
            time = row['time']
            r_shock = row['r_shock']
            v_shock = row['v_shock']
            gamma_shock = row['gamma_shock']
            
            if time not in density_profiles:
                continue
            
            r, rho = density_profiles[time]
            tau, _ = self.compute_optical_depth(r, rho)
            
            # Find τ at shock
            idx_shock = np.argmin(np.abs(r - r_shock))
            tau_at_shock = tau[idx_shock]
            
            if self.is_breakout(tau_at_shock, v_shock):
                # Estimate shell properties
                # Shell extends from shock to some inner radius
                mask = r <= r_shock
                shell_mass = np.trapz(4 * np.pi * r[mask]**2 * rho[mask], r[mask])
                shell_volume = (4/3) * np.pi * r_shock**3
                
                # Internal energy (rough estimate)
                # For strong shock: e_int ≈ (Γ-1) M c² for relativistic
                internal_energy = shell_mass * (gamma_shock - 1)  # in M c² units
                
                return BreakoutState(
                    time=time,
                    r_breakout=r_shock,
                    v_breakout=v_shock,
                    gamma_breakout=gamma_shock,
                    tau_at_shock=tau_at_shock,
                    internal_energy=internal_energy,
                    shell_volume=shell_volume,
                    shell_mass=shell_mass,
                )
        
        return None


def main():
    parser = argparse.ArgumentParser(description='Compute optical depth and breakout')
    parser.add_argument('shock_history', help='Shock history CSV file')
    parser.add_argument('results_dir', help='Directory with snapshots')
    parser.add_argument('--opacity', type=float, default=0.2,
                       help='Opacity κ [cm²/g]')
    parser.add_argument('--relativistic', action='store_true',
                       help='Use relativistic breakout condition')
    parser.add_argument('--output', '-o', default='breakout.json',
                       help='Output file')
    
    args = parser.parse_args()
    
    # Load shock history
    shock_history = pd.read_csv(args.shock_history)
    
    calc = OpticalDepthCalculator(
        opacity=args.opacity,
        relativistic=args.relativistic,
    )
    
    print(f"Loaded shock history with {len(shock_history)} entries")
    print(f"Using opacity κ = {args.opacity} cm²/g")
    print(f"Relativistic: {args.relativistic}")
    
    # For now, just compute optical depth at final snapshot
    # Full implementation would load all snapshots
    print("\nNote: Full breakout detection requires loading all snapshots.")
    print("This is a placeholder - implement density_profiles loading.")


if __name__ == '__main__':
    main()
