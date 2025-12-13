#!/usr/bin/env python3
"""
SRGSPH Particle Sampler

Samples SPH particles from ejecta density profiles to create initial conditions
for SR-GSPH simulations. Handles both 1D (spherical) and 2D (axisymmetric) cases.
"""

import numpy as np
from typing import Tuple, Optional, List, Dict
from dataclasses import dataclass
import json


@dataclass
class Particle:
    """Single SPH particle."""
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    mass: float
    density: float
    pressure: float
    energy: float


@dataclass
class ParticleDistribution:
    """Collection of SPH particles."""
    positions: np.ndarray    # Shape (N, 3)
    velocities: np.ndarray   # Shape (N, 3)
    masses: np.ndarray       # Shape (N,)
    densities: np.ndarray    # Shape (N,)
    pressures: np.ndarray    # Shape (N,)
    energies: np.ndarray     # Shape (N,)
    
    @property
    def n_particles(self) -> int:
        return len(self.masses)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'positions': self.positions.tolist(),
            'velocities': self.velocities.tolist(),
            'masses': self.masses.tolist(),
            'densities': self.densities.tolist(),
            'pressures': self.pressures.tolist(),
            'energies': self.energies.tolist(),
            'n_particles': self.n_particles,
        }


class ParticleSampler:
    """
    Sample SPH particles from density profiles.
    
    Supports:
    - 1D spherical: Radial shells with angular averaging
    - 2D axisymmetric: Radial-angular grid sampling
    
    Usage:
        sampler = ParticleSampler(n_particles=10000)
        particles = sampler.sample_1d(r, rho, v_r, gamma=2.0)
    """
    
    def __init__(self, n_particles: int = 10000, 
                 gamma: float = 2.0,
                 random_seed: int = 42):
        """
        Initialize particle sampler.
        
        Args:
            n_particles: Target number of particles
            gamma: Adiabatic index for EOS
            random_seed: Random seed for reproducibility
        """
        self.n_particles = n_particles
        self.gamma = gamma
        self.rng = np.random.default_rng(random_seed)
    
    def sample_1d(self, r: np.ndarray, rho: np.ndarray, 
                  v_r: np.ndarray, 
                  polytropic_K: float = 1.0) -> ParticleDistribution:
        """
        Sample particles for 1D spherical profile.
        
        Distributes particles in radial shells with uniform angular distribution.
        
        Args:
            r: Radial positions [code units]
            rho: Density at each radius [code units]
            v_r: Radial velocity at each radius [c]
            polytropic_K: Polytropic constant for P = K * rho^gamma
            
        Returns:
            ParticleDistribution
        """
        # Compute mass in each shell
        r_edges = self._compute_edges(r)
        shell_volumes = (4/3) * np.pi * (r_edges[1:]**3 - r_edges[:-1]**3)
        shell_masses = rho * shell_volumes
        total_mass = np.sum(shell_masses)
        
        # Distribute particles proportional to mass
        if total_mass <= 0:
            raise ValueError("Total mass must be positive")
        
        n_per_shell = np.round(self.n_particles * shell_masses / total_mass).astype(int)
        n_per_shell = np.maximum(n_per_shell, 0)
        
        # Ensure we have the right total (adjust largest shell)
        diff = self.n_particles - np.sum(n_per_shell)
        idx_max = np.argmax(n_per_shell)
        n_per_shell[idx_max] += diff
        
        # Sample particles
        positions = []
        velocities = []
        densities = []
        
        for i, (r_in, r_out, n_shell) in enumerate(zip(r_edges[:-1], r_edges[1:], n_per_shell)):
            if n_shell <= 0:
                continue
            
            # Sample radii uniformly in volume
            u = self.rng.uniform(0, 1, n_shell)
            r_sample = (r_in**3 + u * (r_out**3 - r_in**3))**(1/3)
            
            # Sample angles uniformly on sphere
            cos_theta = self.rng.uniform(-1, 1, n_shell)
            phi = self.rng.uniform(0, 2*np.pi, n_shell)
            sin_theta = np.sqrt(1 - cos_theta**2)
            
            # Convert to Cartesian
            x = r_sample * sin_theta * np.cos(phi)
            y = r_sample * sin_theta * np.sin(phi)
            z = r_sample * cos_theta
            positions.append(np.column_stack([x, y, z]))
            
            # Radial velocity (interpolate)
            v_shell = np.interp(r_sample, r, v_r)
            vx = v_shell * sin_theta * np.cos(phi)
            vy = v_shell * sin_theta * np.sin(phi)
            vz = v_shell * cos_theta
            velocities.append(np.column_stack([vx, vy, vz]))
            
            # Density (interpolate)
            rho_shell = np.interp(r_sample, r, rho)
            densities.append(rho_shell)
        
        # Combine
        positions = np.vstack(positions)
        velocities = np.vstack(velocities)
        densities = np.concatenate(densities)
        
        # Compute derived quantities
        n_total = len(densities)
        masses = np.full(n_total, total_mass / n_total)
        
        # Pressure from polytropic EOS
        pressures = polytropic_K * densities**self.gamma
        
        # Specific internal energy
        energies = pressures / (densities * (self.gamma - 1))
        
        return ParticleDistribution(
            positions=positions,
            velocities=velocities,
            masses=masses,
            densities=densities,
            pressures=pressures,
            energies=energies,
        )
    
    def sample_2d(self, r: np.ndarray, theta: np.ndarray,
                  rho_2d: np.ndarray, v_r_2d: np.ndarray,
                  polytropic_K: float = 1.0) -> ParticleDistribution:
        """
        Sample particles for 2D axisymmetric profile.
        
        Args:
            r: Radial positions [code units], shape (n_r,)
            theta: Polar angles [rad], shape (n_theta,)
            rho_2d: Density ρ(θ,r) [code units], shape (n_theta, n_r)
            v_r_2d: Radial velocity v(θ,r) [c], shape (n_theta, n_r)
            polytropic_K: Polytropic constant
            
        Returns:
            ParticleDistribution
        """
        n_r = len(r)
        n_theta = len(theta)
        
        # Compute edges
        r_edges = self._compute_edges(r)
        theta_edges = self._compute_edges(theta)
        
        # Compute mass in each cell
        cell_masses = np.zeros((n_theta, n_r))
        for i_theta in range(n_theta):
            th_in, th_out = theta_edges[i_theta], theta_edges[i_theta + 1]
            dOmega = 2 * np.pi * np.abs(np.cos(th_in) - np.cos(th_out))
            
            for i_r in range(n_r):
                r_in, r_out = r_edges[i_r], r_edges[i_r + 1]
                dV = dOmega * (r_out**3 - r_in**3) / 3
                cell_masses[i_theta, i_r] = rho_2d[i_theta, i_r] * dV
        
        total_mass = np.sum(cell_masses)
        if total_mass <= 0:
            raise ValueError("Total mass must be positive")
        
        # Flatten and compute particles per cell
        flat_masses = cell_masses.flatten()
        n_per_cell = np.round(self.n_particles * flat_masses / total_mass).astype(int)
        n_per_cell = np.maximum(n_per_cell, 0)
        
        # Adjust total
        diff = self.n_particles - np.sum(n_per_cell)
        if diff != 0 and np.sum(n_per_cell) > 0:
            idx_max = np.argmax(n_per_cell)
            n_per_cell[idx_max] = max(0, n_per_cell[idx_max] + diff)
        
        # Sample particles
        positions = []
        velocities = []
        densities = []
        
        idx = 0
        for i_theta in range(n_theta):
            th_in, th_out = theta_edges[i_theta], theta_edges[i_theta + 1]
            
            for i_r in range(n_r):
                n_cell = n_per_cell[idx]
                idx += 1
                
                if n_cell <= 0:
                    continue
                
                r_in, r_out = r_edges[i_r], r_edges[i_r + 1]
                
                # Sample radii (uniform in volume)
                u = self.rng.uniform(0, 1, n_cell)
                r_sample = (r_in**3 + u * (r_out**3 - r_in**3))**(1/3)
                
                # Sample theta (uniform in cos(theta))
                cos_in, cos_out = np.cos(th_in), np.cos(th_out)
                cos_sample = self.rng.uniform(min(cos_in, cos_out), 
                                              max(cos_in, cos_out), n_cell)
                theta_sample = np.arccos(cos_sample)
                sin_theta = np.sin(theta_sample)
                
                # Sample phi uniformly
                phi = self.rng.uniform(0, 2*np.pi, n_cell)
                
                # Convert to Cartesian
                x = r_sample * sin_theta * np.cos(phi)
                y = r_sample * sin_theta * np.sin(phi)
                z = r_sample * cos_sample
                positions.append(np.column_stack([x, y, z]))
                
                # Velocity (use cell value)
                v_r_cell = v_r_2d[i_theta, i_r]
                vx = v_r_cell * sin_theta * np.cos(phi)
                vy = v_r_cell * sin_theta * np.sin(phi)
                vz = v_r_cell * cos_sample
                velocities.append(np.column_stack([vx, vy, vz]))
                
                # Density
                rho_cell = rho_2d[i_theta, i_r]
                densities.append(np.full(n_cell, rho_cell))
        
        if not positions:
            raise ValueError("No particles sampled - check input profiles")
        
        # Combine
        positions = np.vstack(positions)
        velocities = np.vstack(velocities)
        densities = np.concatenate(densities)
        
        n_total = len(densities)
        masses = np.full(n_total, total_mass / n_total)
        pressures = polytropic_K * densities**self.gamma
        energies = pressures / (densities * (self.gamma - 1))
        
        return ParticleDistribution(
            positions=positions,
            velocities=velocities,
            masses=masses,
            densities=densities,
            pressures=pressures,
            energies=energies,
        )
    
    def add_cocoon(self, particles: ParticleDistribution,
                   cocoon_energy: float,
                   cocoon_opening_angle: float,
                   cocoon_velocity: float,
                   n_cocoon_particles: int = 1000) -> ParticleDistribution:
        """
        Add cocoon/jet energy injection region.
        
        Args:
            particles: Existing particle distribution
            cocoon_energy: Total cocoon energy [code units]
            cocoon_opening_angle: Half-opening angle [degrees]
            cocoon_velocity: Initial cocoon velocity [c]
            n_cocoon_particles: Number of cocoon particles
            
        Returns:
            Modified ParticleDistribution
        """
        theta_open = np.radians(cocoon_opening_angle)
        
        # Sample cocoon particles near polar axis
        cos_theta = self.rng.uniform(np.cos(theta_open), 1.0, n_cocoon_particles)
        theta = np.arccos(cos_theta)
        phi = self.rng.uniform(0, 2*np.pi, n_cocoon_particles)
        
        # Place near center (small radius)
        r_cocoon = self.rng.uniform(0.01, 0.1, n_cocoon_particles)
        
        sin_theta = np.sin(theta)
        x = r_cocoon * sin_theta * np.cos(phi)
        y = r_cocoon * sin_theta * np.sin(phi)
        z = r_cocoon * cos_theta
        
        # Radial velocity
        vx = cocoon_velocity * sin_theta * np.cos(phi)
        vy = cocoon_velocity * sin_theta * np.sin(phi)
        vz = cocoon_velocity * cos_theta
        
        # Cocoon properties
        cocoon_mass = n_cocoon_particles * np.mean(particles.masses)
        cocoon_density = np.mean(particles.densities)  # placeholder
        cocoon_energy_per_particle = cocoon_energy / n_cocoon_particles
        
        # Combine with existing particles
        new_positions = np.vstack([particles.positions, 
                                   np.column_stack([x, y, z])])
        new_velocities = np.vstack([particles.velocities,
                                    np.column_stack([vx, vy, vz])])
        new_masses = np.concatenate([particles.masses,
                                     np.full(n_cocoon_particles, cocoon_mass / n_cocoon_particles)])
        new_densities = np.concatenate([particles.densities,
                                        np.full(n_cocoon_particles, cocoon_density)])
        new_pressures = np.concatenate([particles.pressures,
                                        np.full(n_cocoon_particles, cocoon_density * 0.1)])  # hot
        new_energies = np.concatenate([particles.energies,
                                       np.full(n_cocoon_particles, cocoon_energy_per_particle)])
        
        return ParticleDistribution(
            positions=new_positions,
            velocities=new_velocities,
            masses=new_masses,
            densities=new_densities,
            pressures=new_pressures,
            energies=new_energies,
        )
    
    def _compute_edges(self, x: np.ndarray) -> np.ndarray:
        """Compute bin edges from centers."""
        dx = np.diff(x)
        edges = np.zeros(len(x) + 1)
        edges[1:-1] = x[:-1] + dx / 2
        edges[0] = x[0] - dx[0] / 2
        edges[-1] = x[-1] + dx[-1] / 2
        return edges
    
    def save_json(self, particles: ParticleDistribution, filepath: str):
        """Save particles to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(particles.to_dict(), f, indent=2)
        print(f"Saved {particles.n_particles} particles to {filepath}")


def test_particle_sampler():
    """Test the particle sampler."""
    import matplotlib.pyplot as plt
    
    # Create test profile
    r = np.linspace(0.1, 2.0, 50)
    rho = 1.0 * np.exp(-(r - 0.5)**2 / 0.1)  # Gaussian shell
    v_r = 0.2 * np.ones_like(r)  # Constant velocity
    
    # Sample particles
    sampler = ParticleSampler(n_particles=5000, gamma=2.0)
    particles = sampler.sample_1d(r, rho, v_r)
    
    print(f"Sampled {particles.n_particles} particles")
    print(f"Total mass: {np.sum(particles.masses):.4f}")
    print(f"Position range: {np.linalg.norm(particles.positions, axis=1).min():.3f} - "
          f"{np.linalg.norm(particles.positions, axis=1).max():.3f}")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # XY projection
    ax = axes[0, 0]
    ax.scatter(particles.positions[:, 0], particles.positions[:, 1], 
               s=1, c=particles.densities, cmap='viridis', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('XY Projection')
    ax.set_aspect('equal')
    
    # XZ projection
    ax = axes[0, 1]
    ax.scatter(particles.positions[:, 0], particles.positions[:, 2], 
               s=1, c=particles.densities, cmap='viridis', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title('XZ Projection')
    ax.set_aspect('equal')
    
    # Radial histogram
    ax = axes[1, 0]
    r_particles = np.linalg.norm(particles.positions, axis=1)
    ax.hist(r_particles, bins=50, alpha=0.7)
    ax.axvline(0.5, color='r', linestyle='--', label='Profile center')
    ax.set_xlabel('Radius')
    ax.set_ylabel('Count')
    ax.set_title('Radial Distribution')
    ax.legend()
    
    # Velocity histogram
    ax = axes[1, 1]
    v_particles = np.linalg.norm(particles.velocities, axis=1)
    ax.hist(v_particles, bins=50, alpha=0.7)
    ax.set_xlabel('Velocity [c]')
    ax.set_ylabel('Count')
    ax.set_title('Velocity Distribution')
    
    plt.tight_layout()
    plt.savefig('particle_sampler_test.png', dpi=150)
    print("Saved: particle_sampler_test.png")


if __name__ == '__main__':
    test_particle_sampler()
