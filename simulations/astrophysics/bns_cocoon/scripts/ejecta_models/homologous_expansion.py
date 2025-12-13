#!/usr/bin/env python3
"""
Homologous Expansion Module

Converts velocity distributions to spatial distributions assuming homologous expansion:
    r = v × t₀

where t₀ is a reference time after the merger when expansion becomes approximately
free (typically ~1 second for dynamical ejecta).
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass
class SpatialProfile:
    """Container for spatial density/velocity profile."""
    # 1D radial profile
    r: np.ndarray              # Radial positions [code units]
    rho: np.ndarray            # Density at each radius [code units]
    v_r: np.ndarray            # Radial velocity at each radius [c]
    
    # 2D axisymmetric profile (optional)
    theta: Optional[np.ndarray] = None     # Polar angle bins [rad]
    rho_2d: Optional[np.ndarray] = None    # Density ρ(r,θ) [code units]
    v_r_2d: Optional[np.ndarray] = None    # Radial velocity v(r,θ) [c]
    
    # Domain
    r_min: float = 0.0
    r_max: float = 1.0
    
    # Reference time
    t0: float = 1.0            # Reference time [s]
    
    # Units
    length_unit: float = 1.0   # Length unit [cm]
    density_unit: float = 1.0  # Density unit [g/cm³]


class HomologousExpander:
    """
    Convert velocity-based ejecta profiles to spatial density distributions.
    
    Assumes homologous expansion where each mass element at velocity v is
    located at radius r = v × t₀ at reference time t₀.
    
    Usage:
        expander = HomologousExpander(t0=1.0)  # t0 in seconds
        spatial = expander.expand_1d(v_inf, dM_dv)
    """
    
    # Physical constants
    C_CGS = 2.998e10      # Speed of light [cm/s]
    M_SUN_CGS = 1.989e33  # Solar mass [g]
    
    def __init__(self, t0: float = 1.0, length_unit: float = 1e10):
        """
        Initialize homologous expander.
        
        Args:
            t0: Reference time after merger [seconds]
            length_unit: Code length unit [cm], default 10^10 cm ≈ 0.1 AU
        """
        self.t0 = t0
        self.length_unit = length_unit
        self.density_unit = self.M_SUN_CGS / length_unit**3
        
    def v_to_r(self, v: np.ndarray) -> np.ndarray:
        """
        Convert velocity to radius at reference time.
        
        Args:
            v: Velocity [c]
            
        Returns:
            r: Radius [code units]
        """
        r_cm = v * self.C_CGS * self.t0
        return r_cm / self.length_unit
    
    def expand_1d(self, v_inf: np.ndarray, dM_dv: np.ndarray) -> SpatialProfile:
        """
        Convert 1D velocity distribution to radial density profile.
        
        Under homologous expansion:
            r = v × t₀
            dr = t₀ × dv
            ρ(r) × 4πr² dr = dM/dv × dv
            ρ(r) = (dM/dv) / (4πr² t₀)
        
        Args:
            v_inf: Velocity bins [c]
            dM_dv: Differential mass dM/dv [M_sun/c]
            
        Returns:
            SpatialProfile with radial density and velocity
        """
        # Convert velocities to radii
        r = self.v_to_r(v_inf)
        
        # Avoid division by zero
        r_safe = np.maximum(r, 1e-10)
        
        # Compute density
        # dM/dv = ρ × 4πr² × dr/dv = ρ × 4πr² × t₀ × c
        # ρ = (dM/dv) / (4πr² × t₀ × c)
        rho_cgs = dM_dv * self.M_SUN_CGS / (4 * np.pi * (r_safe * self.length_unit)**2 
                                            * self.t0 * self.C_CGS)
        rho = rho_cgs / self.density_unit
        
        # Radial velocity equals v (in c)
        v_r = v_inf.copy()
        
        return SpatialProfile(
            r=r,
            rho=rho,
            v_r=v_r,
            r_min=r.min(),
            r_max=r.max(),
            t0=self.t0,
            length_unit=self.length_unit,
            density_unit=self.density_unit,
        )
    
    def expand_2d(self, v_inf: np.ndarray, cos_theta: np.ndarray,
                  dM_dv_dOmega: np.ndarray) -> SpatialProfile:
        """
        Convert 2D velocity-angle distribution to axisymmetric density profile.
        
        Under homologous expansion:
            r = v × t₀
            ρ(r,θ) = (dM/dv/dΩ) / (r² × t₀)
        
        Args:
            v_inf: Velocity bins [c], shape (n_v,)
            cos_theta: Angular bins cos(θ), shape (n_theta,)
            dM_dv_dOmega: Differential mass dM/(dv dΩ) [M_sun/c/sr], shape (n_theta, n_v)
            
        Returns:
            SpatialProfile with 2D density and velocity
        """
        # Convert velocities to radii
        r = self.v_to_r(v_inf)
        theta = np.arccos(cos_theta)
        
        # Avoid division by zero
        r_safe = np.maximum(r, 1e-10)
        
        # Compute 2D density
        # dM/(dv dΩ) = ρ × r² × dr/dv = ρ × r² × t₀ × c
        # ρ(r,θ) = (dM/dv/dΩ) / (r² × t₀ × c)
        rho_2d_cgs = np.zeros((len(theta), len(r)))
        for i_theta in range(len(theta)):
            for i_r in range(len(r)):
                rho_2d_cgs[i_theta, i_r] = (dM_dv_dOmega[i_theta, i_r] * self.M_SUN_CGS 
                                           / ((r_safe[i_r] * self.length_unit)**2 
                                              * self.t0 * self.C_CGS))
        
        rho_2d = rho_2d_cgs / self.density_unit
        
        # Radial velocity at each (r, θ) - in homologous expansion, v_r = r/t₀
        # But we store velocity in units of c
        v_r_2d = np.outer(np.ones_like(theta), v_inf)
        
        # Also compute angle-averaged 1D profile
        dOmega = 2 * np.pi * np.abs(np.gradient(cos_theta))
        rho_1d = np.sum(rho_2d * dOmega[:, np.newaxis], axis=0) / (4 * np.pi)
        v_r_1d = v_inf.copy()
        
        return SpatialProfile(
            r=r,
            rho=rho_1d,
            v_r=v_r_1d,
            theta=theta,
            rho_2d=rho_2d,
            v_r_2d=v_r_2d,
            r_min=r.min(),
            r_max=r.max(),
            t0=self.t0,
            length_unit=self.length_unit,
            density_unit=self.density_unit,
        )
    
    def add_external_medium(self, profile: SpatialProfile, 
                            n_ism: float = 1e-3,
                            r_extend: float = 10.0) -> SpatialProfile:
        """
        Extend profile with external interstellar medium.
        
        Args:
            profile: Input spatial profile
            n_ism: ISM number density [cm⁻³]
            r_extend: Factor to extend radial domain
            
        Returns:
            Extended SpatialProfile
        """
        # ISM density (assuming proton mass)
        m_p = 1.67e-24  # g
        rho_ism_cgs = n_ism * m_p
        rho_ism = rho_ism_cgs / self.density_unit
        
        # Extend radial grid
        r_max_new = profile.r_max * r_extend
        n_extend = int(len(profile.r) * (r_extend - 1))
        r_extend_arr = np.linspace(profile.r_max, r_max_new, n_extend)
        
        # Combine
        r_new = np.concatenate([profile.r, r_extend_arr])
        rho_new = np.concatenate([profile.rho, np.ones(n_extend) * rho_ism])
        v_r_new = np.concatenate([profile.v_r, np.zeros(n_extend)])
        
        return SpatialProfile(
            r=r_new,
            rho=rho_new,
            v_r=v_r_new,
            r_min=profile.r_min,
            r_max=r_max_new,
            t0=profile.t0,
            length_unit=profile.length_unit,
            density_unit=profile.density_unit,
        )


def test_homologous_expansion():
    """Test the homologous expansion module."""
    import matplotlib.pyplot as plt
    
    # Create synthetic velocity distribution (power law)
    v_inf = np.linspace(0.05, 0.5, 50)
    # dM/dv ∝ v^(-2) with total mass ~0.01 M_sun
    dM_dv = 0.001 * (v_inf / 0.2)**(-2)
    
    # Expand to spatial profile
    expander = HomologousExpander(t0=1.0, length_unit=1e10)
    profile = expander.expand_1d(v_inf, dM_dv)
    
    # Add ISM
    profile_with_ism = expander.add_external_medium(profile, n_ism=1e-3)
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Velocity distribution
    ax = axes[0, 0]
    ax.semilogy(v_inf, dM_dv)
    ax.set_xlabel('Velocity [c]')
    ax.set_ylabel('dM/dv [M☉/c]')
    ax.set_title('Input Velocity Distribution')
    
    # Density profile
    ax = axes[0, 1]
    ax.loglog(profile.r, profile.rho)
    ax.set_xlabel('Radius [code units]')
    ax.set_ylabel('Density [code units]')
    ax.set_title(f'Density Profile (t₀ = {expander.t0} s)')
    
    # Velocity profile
    ax = axes[1, 0]
    ax.plot(profile.r, profile.v_r)
    ax.set_xlabel('Radius [code units]')
    ax.set_ylabel('Velocity [c]')
    ax.set_title('Velocity Profile')
    
    # Extended profile with ISM
    ax = axes[1, 1]
    ax.loglog(profile_with_ism.r, profile_with_ism.rho)
    ax.axhline(1e-3 * 1.67e-24 / expander.density_unit, 
               color='r', linestyle='--', label='ISM')
    ax.set_xlabel('Radius [code units]')
    ax.set_ylabel('Density [code units]')
    ax.set_title('Extended Profile with ISM')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('homologous_test.png', dpi=150)
    print("Saved: homologous_test.png")
    
    # Print summary
    print(f"\n=== Homologous Expansion Test ===")
    print(f"Reference time: {expander.t0} s")
    print(f"Length unit: {expander.length_unit:.2e} cm")
    print(f"Density unit: {expander.density_unit:.2e} g/cm³")
    print(f"Radial range: {profile.r.min():.3f} - {profile.r.max():.3f} code units")
    print(f"Radial range: {profile.r.min() * expander.length_unit:.3e} - "
          f"{profile.r.max() * expander.length_unit:.3e} cm")


if __name__ == '__main__':
    test_homologous_expansion()
