#!/usr/bin/env python3
"""
Analytic Ejecta Profile Models

Provides parametric models for BNS merger ejecta based on fits to GR simulations:
- Hotokezaka et al. (2013, 2015): Velocity distribution
- Bauswein, Goriely & Janka (2013): Angular distribution
- Generic power-law models

These can be used when full numerical data are not needed or for parameter surveys.
"""

import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple, Optional


@dataclass
class EjectaParameters:
    """Parameters defining an ejecta model."""
    total_mass: float = 0.01        # Total ejecta mass [M_sun]
    avg_velocity: float = 0.2       # Mass-weighted average velocity [c]
    max_velocity: float = 0.5       # Maximum (cutoff) velocity [c]
    min_velocity: float = 0.05      # Minimum velocity [c]
    velocity_power: float = -2.0    # Power-law index for dM/dv ∝ v^n
    angular_conc: float = 2.0       # Angular concentration parameter
    polar_opening: float = 60.0     # Polar opening angle [degrees]


class EjectaProfile(ABC):
    """Abstract base class for ejecta profiles."""
    
    @abstractmethod
    def velocity_distribution(self, v: np.ndarray) -> np.ndarray:
        """
        Compute differential mass distribution dM/dv.
        
        Args:
            v: Velocity array [c]
            
        Returns:
            dM/dv [M_sun/c]
        """
        pass
    
    @abstractmethod
    def angular_distribution(self, theta: np.ndarray) -> np.ndarray:
        """
        Compute angular mass distribution dM/dΩ.
        
        Args:
            theta: Polar angle array [radians]
            
        Returns:
            dM/dΩ [M_sun/sr]
        """
        pass
    
    def velocity_angle_distribution(self, v: np.ndarray, 
                                   theta: np.ndarray) -> np.ndarray:
        """
        Compute 2D velocity-angle distribution dM/(dv dΩ).
        
        Default: assumes separable distribution.
        
        Args:
            v: Velocity array [c]
            theta: Polar angle array [radians]
            
        Returns:
            dM/(dv dΩ) [M_sun/c/sr], shape (n_theta, n_v)
        """
        dM_dv = self.velocity_distribution(v)
        dM_dOmega = self.angular_distribution(theta)
        
        # Normalize angular distribution to integrate to 1
        dcos = -np.gradient(np.cos(theta))
        norm = 2 * np.pi * np.sum(dM_dOmega * dcos)
        dM_dOmega_norm = dM_dOmega / norm
        
        return np.outer(dM_dOmega_norm, dM_dv)


class HotokezakaProfile(EjectaProfile):
    """
    Ejecta profile based on Hotokezaka et al. (2013, 2015).
    
    Key features:
    - Kinetic energy distribution: E(≥β) ∝ β^(-0.5)
    - Typical average velocity: β ≈ 0.2
    - Fast tail up to β ≈ 0.4-0.8
    
    This translates to a velocity distribution:
        dM/dβ ∝ β^(-2) × exp(-β/β_max)
    """
    
    def __init__(self, params: Optional[EjectaParameters] = None):
        """
        Initialize Hotokezaka profile.
        
        Args:
            params: EjectaParameters, uses defaults if None
        """
        self.params = params or EjectaParameters()
    
    def velocity_distribution(self, v: np.ndarray) -> np.ndarray:
        """
        Compute dM/dv following Hotokezaka et al. fits.
        
        Form: dM/dv ∝ v^(-2) × exp(-v/v_max) for v > v_min
        """
        p = self.params
        
        # Power-law with exponential cutoff
        dM_dv = np.zeros_like(v)
        mask = (v >= p.min_velocity) & (v <= p.max_velocity)
        
        # v^(-2) behavior from E(≥β) ∝ β^(-0.5)
        # E = (Γ-1)Mc² ≈ 0.5 × v² × M for v << c
        # dE/dv = v × dM/dv
        # E(≥v) = ∫_v^∞ v' dM/dv' dv' ∝ v^(-0.5)
        # => dM/dv ∝ v^(-2.5) approximately
        # We use v^(-2) for simplicity
        
        dM_dv[mask] = (v[mask] / p.avg_velocity)**p.velocity_power
        dM_dv[mask] *= np.exp(-v[mask] / p.max_velocity)
        
        # Normalize to total mass
        dv = np.gradient(v)
        total = np.sum(dM_dv * dv)
        if total > 0:
            dM_dv *= p.total_mass / total
        
        return dM_dv
    
    def angular_distribution(self, theta: np.ndarray) -> np.ndarray:
        """
        Compute dM/dΩ - concentration near equatorial plane.
        
        Form: dM/dΩ ∝ 1 + A × sin^n(θ)
        Ejecta concentrated near θ = π/2 (equator)
        """
        p = self.params
        
        # Angular concentration
        A = p.angular_conc
        n = 4  # Controls sharpness
        
        # sin^n(θ) gives concentration at equator
        dM_dOmega = 1.0 + A * np.sin(theta)**n
        
        return dM_dOmega


class BausweinProfile(EjectaProfile):
    """
    Ejecta profile based on Bauswein, Goriely & Janka (2013).
    
    Key features:
    - Ejecta concentrated within ~60° from orbital plane
    - More detailed angular structure
    - EOS and mass-ratio dependence
    """
    
    def __init__(self, params: Optional[EjectaParameters] = None):
        """
        Initialize Bauswein profile.
        
        Args:
            params: EjectaParameters, uses defaults if None
        """
        self.params = params or EjectaParameters()
    
    def velocity_distribution(self, v: np.ndarray) -> np.ndarray:
        """
        Compute dM/dv with Bauswein-style distribution.
        """
        p = self.params
        
        dM_dv = np.zeros_like(v)
        mask = (v >= p.min_velocity) & (v <= p.max_velocity)
        
        # Power-law with softer cutoff
        dM_dv[mask] = (v[mask] / p.avg_velocity)**p.velocity_power
        
        # Smooth high-velocity cutoff
        sigma_v = 0.05  # Width of cutoff
        dM_dv[mask] *= 0.5 * (1 - np.tanh((v[mask] - p.max_velocity) / sigma_v))
        
        # Normalize
        dv = np.gradient(v)
        total = np.sum(dM_dv * dv)
        if total > 0:
            dM_dv *= p.total_mass / total
        
        return dM_dv
    
    def angular_distribution(self, theta: np.ndarray) -> np.ndarray:
        """
        Compute dM/dΩ with piecewise angular structure.
        
        Bauswein+2013 find ejecta concentrated within θ_open from equator.
        """
        p = self.params
        
        theta_open = np.radians(p.polar_opening)
        theta_equator = np.pi / 2
        
        # Distance from equatorial plane
        delta_theta = np.abs(theta - theta_equator)
        
        # Gaussian profile centered on equator
        sigma_theta = theta_open / 2
        dM_dOmega = np.exp(-delta_theta**2 / (2 * sigma_theta**2))
        
        return dM_dOmega


class PowerLawProfile(EjectaProfile):
    """
    Simple power-law ejecta profile for parameter studies.
    
    Allows easy variation of velocity and angular distribution indices.
    """
    
    def __init__(self, params: Optional[EjectaParameters] = None,
                 velocity_index: float = -2.0,
                 angular_index: float = 4.0):
        """
        Initialize power-law profile.
        
        Args:
            params: EjectaParameters
            velocity_index: Power-law index for dM/dv ∝ v^n
            angular_index: Power for sin^n(θ) angular distribution
        """
        self.params = params or EjectaParameters()
        self.velocity_index = velocity_index
        self.angular_index = angular_index
    
    def velocity_distribution(self, v: np.ndarray) -> np.ndarray:
        """Pure power-law velocity distribution with cutoffs."""
        p = self.params
        
        dM_dv = np.zeros_like(v)
        mask = (v >= p.min_velocity) & (v <= p.max_velocity)
        
        dM_dv[mask] = (v[mask] / p.avg_velocity)**self.velocity_index
        
        # Normalize
        dv = np.gradient(v)
        total = np.sum(dM_dv * dv)
        if total > 0:
            dM_dv *= p.total_mass / total
        
        return dM_dv
    
    def angular_distribution(self, theta: np.ndarray) -> np.ndarray:
        """Sin^n angular distribution."""
        return np.sin(theta)**self.angular_index


class FastTailProfile(EjectaProfile):
    """
    Two-component ejecta: bulk + fast tail.
    
    Motivated by observations of GW170817 suggesting a fast (v > 0.4c)
    component in addition to bulk ejecta at v ~ 0.2c.
    """
    
    def __init__(self, 
                 bulk_mass: float = 0.01,
                 bulk_velocity: float = 0.2,
                 tail_mass: float = 1e-4,
                 tail_velocity: float = 0.6,
                 tail_index: float = -3.0):
        """
        Initialize two-component profile.
        
        Args:
            bulk_mass: Mass in bulk component [M_sun]
            bulk_velocity: Characteristic velocity of bulk [c]
            tail_mass: Mass in fast tail [M_sun]
            tail_velocity: Characteristic velocity of tail [c]
            tail_index: Power-law index for tail
        """
        self.bulk_mass = bulk_mass
        self.bulk_velocity = bulk_velocity
        self.tail_mass = tail_mass
        self.tail_velocity = tail_velocity
        self.tail_index = tail_index
        
        # Total parameters
        self.params = EjectaParameters(
            total_mass=bulk_mass + tail_mass,
            avg_velocity=bulk_velocity,
            max_velocity=0.95,  # Relativistic limit
            min_velocity=0.05,
        )
    
    def velocity_distribution(self, v: np.ndarray) -> np.ndarray:
        """Two-component velocity distribution."""
        # Bulk component - Gaussian-like
        dM_dv_bulk = self.bulk_mass * np.exp(-(v - self.bulk_velocity)**2 
                                              / (2 * 0.05**2))
        dM_dv_bulk /= np.sqrt(2 * np.pi) * 0.05
        
        # Fast tail - power law
        dM_dv_tail = np.zeros_like(v)
        mask = v >= self.tail_velocity
        if np.any(mask):
            dM_dv_tail[mask] = (v[mask] / self.tail_velocity)**self.tail_index
            # Normalize tail
            dv = np.gradient(v)
            tail_total = np.sum(dM_dv_tail * dv)
            if tail_total > 0:
                dM_dv_tail *= self.tail_mass / tail_total
        
        return dM_dv_bulk + dM_dv_tail
    
    def angular_distribution(self, theta: np.ndarray) -> np.ndarray:
        """Use equatorial concentration."""
        return 1.0 + 2.0 * np.sin(theta)**4


def create_profile(profile_type: str, **kwargs) -> EjectaProfile:
    """
    Factory function to create ejecta profiles.
    
    Args:
        profile_type: One of 'hotokezaka', 'bauswein', 'powerlaw', 'fasttail'
        **kwargs: Parameters to pass to profile constructor
        
    Returns:
        EjectaProfile instance
    """
    profiles = {
        'hotokezaka': HotokezakaProfile,
        'bauswein': BausweinProfile,
        'powerlaw': PowerLawProfile,
        'fasttail': FastTailProfile,
    }
    
    if profile_type.lower() not in profiles:
        raise ValueError(f"Unknown profile type: {profile_type}. "
                        f"Available: {list(profiles.keys())}")
    
    return profiles[profile_type.lower()](**kwargs)


def test_profiles():
    """Test and plot different ejecta profiles."""
    import matplotlib.pyplot as plt
    
    # Velocity grid
    v = np.linspace(0.01, 0.8, 100)
    theta = np.linspace(0, np.pi, 50)
    
    # Create profiles
    params = EjectaParameters(total_mass=0.01, avg_velocity=0.2)
    profiles = {
        'Hotokezaka': HotokezakaProfile(params),
        'Bauswein': BausweinProfile(params),
        'Power-law': PowerLawProfile(params),
        'Fast tail': FastTailProfile(),
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Velocity distributions
    ax = axes[0, 0]
    for name, prof in profiles.items():
        dM_dv = prof.velocity_distribution(v)
        ax.semilogy(v, dM_dv, label=name)
    ax.set_xlabel('Velocity [c]')
    ax.set_ylabel('dM/dv [M☉/c]')
    ax.set_title('Velocity Distributions')
    ax.legend()
    ax.set_xlim(0, 0.8)
    ax.grid(True, alpha=0.3)
    
    # Angular distributions
    ax = axes[0, 1]
    for name, prof in profiles.items():
        dM_dOmega = prof.angular_distribution(theta)
        ax.plot(np.degrees(theta), dM_dOmega, label=name)
    ax.set_xlabel('Polar angle θ [degrees]')
    ax.set_ylabel('dM/dΩ [arb. units]')
    ax.set_title('Angular Distributions')
    ax.axvline(90, color='k', linestyle='--', alpha=0.3, label='Equator')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2D distribution example (Hotokezaka)
    ax = axes[1, 0]
    prof = profiles['Hotokezaka']
    dM_2d = prof.velocity_angle_distribution(v, theta)
    im = ax.pcolormesh(v, np.degrees(theta), np.log10(dM_2d + 1e-10), 
                       cmap='viridis', shading='auto')
    ax.set_xlabel('Velocity [c]')
    ax.set_ylabel('Polar angle [degrees]')
    ax.set_title('Hotokezaka 2D Distribution\nlog₁₀(dM/dv/dΩ)')
    plt.colorbar(im, ax=ax)
    
    # Integrated quantities
    ax = axes[1, 1]
    ax.axis('off')
    text = "Profile Summary:\n\n"
    for name, prof in profiles.items():
        dM_dv = prof.velocity_distribution(v)
        dv = np.gradient(v)
        M_tot = np.sum(dM_dv * dv)
        v_avg = np.sum(v * dM_dv * dv) / M_tot if M_tot > 0 else 0
        text += f"{name}:\n"
        text += f"  Total mass: {M_tot:.4f} M☉\n"
        text += f"  Avg velocity: {v_avg:.3f} c\n\n"
    ax.text(0.1, 0.9, text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig('analytic_profiles_test.png', dpi=150)
    print("Saved: analytic_profiles_test.png")


if __name__ == '__main__':
    test_profiles()
