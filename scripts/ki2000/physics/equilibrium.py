#!/usr/bin/env python3
"""
Thermal equilibrium solver for ISM gas.

Solves for the equilibrium temperature where heating = cooling,
reproducing K&I 2000 Figure 1c from first principles.

Usage:
    python equilibrium.py [--N_H 1e19] [--output equilibrium.png]
    
    # Or as a module:
    from ki2000.physics.equilibrium import ThermalEquilibrium
    eq = ThermalEquilibrium(N_H=1e19)
    T, x_e = eq.solve_equilibrium(n=1.0)
"""

import numpy as np
from scipy.optimize import brentq, fsolve
from typing import Union, Tuple, Optional
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from physics.heating import (
    total_heating_rate, heating_breakdown,
    photoelectric_heating, cosmic_ray_heating, xray_heating,
    h2_formation_heating, h2_photodissociation_heating, photoelectric_cooling
)
from physics.cooling import (
    total_cooling_rate, cooling_breakdown,
    cii_cooling, oi_cooling, lyman_alpha_cooling,
    h2_cooling, co_cooling, gas_grain_cooling
)

# Physical constants
K_B = 1.380649e-16  # Boltzmann constant [erg/K]


class ChemicalEquilibrium:
    """
    Simplified chemical equilibrium for electron fraction.
    
    For a full treatment, see ki2000_chemistry_ode.py.
    This provides approximate x_e(n, T) for thermal equilibrium.
    """
    
    def __init__(
        self,
        zeta_CR: float = 1.8e-17,
        x_C: float = 3e-4,
        x_O: float = 4.6e-4
    ):
        """
        Parameters
        ----------
        zeta_CR : float
            Cosmic ray ionization rate [s^-1]
        x_C : float
            Carbon abundance (singly ionized in diffuse ISM)
        x_O : float
            Oxygen abundance (neutral via charge exchange)
        """
        self.zeta_CR = zeta_CR
        self.x_C = x_C
        self.x_O = x_O
    
    def electron_fraction(
        self,
        n: Union[float, np.ndarray],
        T: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """
        Calculate equilibrium electron fraction.
        
        In diffuse ISM, electrons come from:
        - C+ (carbon has IP < 13.6 eV, so ionized by UV)
        - Cosmic ray ionization of H
        - Residual ionization from recombination balance
        
        Parameters
        ----------
        n : float or array
            Hydrogen nucleus density [cm^-3]
        T : float or array
            Gas temperature [K]
        
        Returns
        -------
        float or array
            Electron fraction x_e = n_e/n
        """
        T_safe = np.maximum(T, 10.0)
        n_safe = np.maximum(n, 1e-10)
        
        # At low density: x_e ~ x_C (carbon fully ionized)
        # At high density: recombination dominates
        
        # Recombination coefficient (case B, approximate)
        alpha_B = 2.6e-13 * (T_safe / 1e4)**(-0.7)  # cm^3/s
        
        # Ionization rate (cosmic ray + residual UV)
        ionization_rate = self.zeta_CR
        
        # Equilibrium: ionization = recombination
        # zeta * n_H = alpha * n_e^2
        # x_e = sqrt(zeta / (alpha * n))
        x_e_CR = np.sqrt(ionization_rate / (alpha_B * n_safe + 1e-30))
        
        # Add carbon contribution (C+ provides floor)
        x_e = np.sqrt(x_e_CR**2 + self.x_C**2)
        
        # Limit to physically reasonable range
        x_e = np.clip(x_e, 1e-8, 1.0)
        
        return x_e
    
    def h2_fraction(
        self,
        n: Union[float, np.ndarray],
        T: Union[float, np.ndarray],
        N_H: float = 1e19
    ) -> Union[float, np.ndarray]:
        """
        Approximate equilibrium H2 fraction.
        
        Parameters
        ----------
        n : float or array
            Hydrogen nucleus density [cm^-3]
        T : float or array
            Gas temperature [K]
        N_H : float
            Shielding column density [cm^-2]
        
        Returns
        -------
        float or array
            H2 fraction x_H2 = n_H2/n
        """
        T_safe = np.maximum(T, 10.0)
        
        # Formation rate coefficient on grains
        S = 1.0 / (1.0 + 0.04 * np.sqrt(T_safe + 8.0) + 2e-3 * T_safe + 8e-6 * T_safe**2)
        R = 6e-17 * (T_safe / 300.0)**0.5 * S  # cm^3/s
        
        # Photodissociation rate (with self-shielding)
        # R_pd ~ 3e-11 * G0 * f_shield [s^-1]
        G0 = 1.7
        # Approximate self-shielding based on column
        N_H2_approx = N_H * 0.01  # Very rough
        f_shield = np.minimum(1.0, (N_H2_approx / 1e14)**(-0.75))
        R_pd = 3e-11 * G0 * f_shield
        
        # Equilibrium: R * n * x_H = R_pd * x_H2
        # x_H2 / x_H = R * n / R_pd
        ratio = R * n / (R_pd + 1e-30)
        x_H2 = ratio / (1.0 + 2.0 * ratio)
        
        # Limit range
        x_H2 = np.clip(x_H2, 1e-10, 0.5)
        
        return x_H2


class ThermalEquilibrium:
    """
    Solve for thermal equilibrium temperature.
    """
    
    def __init__(
        self,
        N_H: float = 1e19,
        G0: float = 1.7,
        zeta_CR: float = 1.8e-17
    ):
        """
        Parameters
        ----------
        N_H : float
            Shielding column density [cm^-2]
        G0 : float
            FUV field strength in Habing units
        zeta_CR : float
            Cosmic ray ionization rate [s^-1]
        """
        self.N_H = N_H
        self.G0 = G0
        self.zeta_CR = zeta_CR
        self.chem = ChemicalEquilibrium(zeta_CR=zeta_CR)
    
    def net_heating(
        self,
        T: float,
        n: float,
        x_e: Optional[float] = None,
        x_H2: Optional[float] = None
    ) -> float:
        """
        Net heating rate (heating - cooling).
        
        Zero at thermal equilibrium.
        """
        if x_e is None:
            x_e = self.chem.electron_fraction(n, T)
        if x_H2 is None:
            x_H2 = self.chem.h2_fraction(n, T, self.N_H)
        
        x_H = 1.0 - 2.0 * x_H2
        
        gamma = total_heating_rate(n, T, x_e, x_H, x_H2, self.G0, self.N_H, self.zeta_CR)
        lambd = total_cooling_rate(n, T, x_e, x_H, x_H2)
        
        return gamma - lambd
    
    def solve_equilibrium(
        self,
        n: float,
        T_range: Tuple[float, float] = (10.0, 2e4)
    ) -> Tuple[float, float, float]:
        """
        Find equilibrium temperature for given density.
        
        Parameters
        ----------
        n : float
            Hydrogen nucleus density [cm^-3]
        T_range : tuple
            Temperature search range [K]
        
        Returns
        -------
        T_eq : float
            Equilibrium temperature [K]
        x_e : float
            Electron fraction at equilibrium
        x_H2 : float
            H2 fraction at equilibrium
        """
        T_min, T_max = T_range
        
        # Check if solution exists in range
        net_min = self.net_heating(T_min, n)
        net_max = self.net_heating(T_max, n)
        
        # In thermally unstable region, there may be multiple solutions
        # Try to find at least one
        try:
            T_eq = brentq(lambda T: self.net_heating(T, n), T_min, T_max)
        except ValueError:
            # No zero crossing - return boundary value
            if abs(net_min) < abs(net_max):
                T_eq = T_min
            else:
                T_eq = T_max
        
        x_e = self.chem.electron_fraction(n, T_eq)
        x_H2 = self.chem.h2_fraction(n, T_eq, self.N_H)
        
        return T_eq, x_e, x_H2
    
    def equilibrium_curve(
        self,
        log_n_range: Tuple[float, float] = (-2, 6),
        n_points: int = 100
    ) -> dict:
        """
        Calculate equilibrium curve T(n) and heating/cooling rates.
        
        Returns arrays for plotting panel C of K&I 2000 Figure 1.
        
        Parameters
        ----------
        log_n_range : tuple
            Range of log10(n) to compute
        n_points : int
            Number of density points
        
        Returns
        -------
        dict with keys:
            'log_n': log10 of density
            'n': density [cm^-3]
            'T': equilibrium temperature [K]
            'x_e': electron fraction
            'x_H2': H2 fraction
            'heating': dict of heating rates
            'cooling': dict of cooling rates
        """
        log_n = np.linspace(log_n_range[0], log_n_range[1], n_points)
        n = 10**log_n
        
        T_eq = np.zeros(n_points)
        x_e = np.zeros(n_points)
        x_H2 = np.zeros(n_points)
        
        for i, ni in enumerate(n):
            T_eq[i], x_e[i], x_H2[i] = self.solve_equilibrium(ni)
        
        x_H = 1.0 - 2.0 * x_H2
        
        # Calculate individual rates at equilibrium
        heating = heating_breakdown(n, T_eq, x_e, x_H, x_H2, self.G0, self.N_H, self.zeta_CR)
        cooling = cooling_breakdown(n, T_eq, x_e, x_H, x_H2)
        
        return {
            'log_n': log_n,
            'n': n,
            'T': T_eq,
            'x_e': x_e,
            'x_H2': x_H2,
            'heating': heating,
            'cooling': cooling
        }


def reproduce_ki2000_panel_c(output_path: str = None, N_H: float = 1e19):
    """
    Reproduce K&I 2000 Figure 1c from first principles.
    
    Parameters
    ----------
    output_path : str, optional
        Path to save figure. If None, display interactively.
    N_H : float
        Shielding column density [cm^-2]
    """
    import matplotlib.pyplot as plt
    
    eq = ThermalEquilibrium(N_H=N_H)
    result = eq.equilibrium_curve(log_n_range=(-2, 6), n_points=200)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Heating rates (dashed lines)
    heating = result['heating']
    log_n = result['log_n']
    
    heating_processes = [
        ('PE', 'Photoelectric (PE)', 'C0'),
        ('CR', 'Cosmic Ray (CR)', 'C1'),
        ('XR', 'X-ray (XR)', 'C2'),
        ('H2_form', 'H₂ formation (H₂)', 'C3'),
    ]
    
    for key, label, color in heating_processes:
        rate = heating[key]
        valid = rate > 1e-35
        if np.any(valid):
            ax.semilogy(log_n[valid], rate[valid], '--', color=color, label=f'{label} (heat)', lw=1.5)
    
    # Cooling rates (solid lines)
    cooling = result['cooling']
    
    cooling_processes = [
        ('CII', '[CII] 158μm', 'C4'),
        ('OI', '[OI] 63μm', 'C5'),
        ('Lya', 'Ly-α', 'C6'),
        ('CO', 'CO', 'C7'),
        ('GR', 'Gas-Grain', 'C8'),
    ]
    
    for key, label, color in cooling_processes:
        rate = cooling[key]
        valid = rate > 1e-35
        if np.any(valid):
            ax.semilogy(log_n[valid], rate[valid], '-', color=color, label=f'{label} (cool)', lw=1.5)
    
    # Total rates
    ax.semilogy(log_n, heating['total'], 'k--', lw=2, label='Total heating')
    ax.semilogy(log_n, cooling['total'], 'k-', lw=2, label='Total cooling')
    
    ax.set_xlabel('log n [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Rate per H nucleus [erg s$^{-1}$]', fontsize=12)
    ax.set_title(f'Heating & Cooling at Thermal Equilibrium (N_H = {N_H:.0e} cm$^{{-2}}$)', fontsize=14)
    ax.set_xlim(-2, 6)
    ax.set_ylim(1e-28, 1e-23)
    ax.legend(loc='upper right', fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved figure to {output_path}")
    else:
        plt.show()
    
    plt.close()
    
    return result


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Compute thermal equilibrium for ISM')
    parser.add_argument('--N_H', type=float, default=1e19,
                       help='Shielding column density [cm^-2]')
    parser.add_argument('--output', type=str, default=None,
                       help='Output figure path')
    args = parser.parse_args()
    
    output = args.output
    if output is None:
        output = f'/Users/guo/Downloads/sphcode/scripts/ki2000/physics/panel_c_N{int(np.log10(args.N_H))}_firstprinciples.png'
    
    result = reproduce_ki2000_panel_c(output_path=output, N_H=args.N_H)
    
    print(f"\nEquilibrium summary (N_H = {args.N_H:.0e}):")
    print(f"  log(n) range: {result['log_n'][0]:.1f} to {result['log_n'][-1]:.1f}")
    print(f"  T range: {result['T'].min():.0f} to {result['T'].max():.0f} K")
    print(f"  x_e range: {result['x_e'].min():.2e} to {result['x_e'].max():.2e}")
