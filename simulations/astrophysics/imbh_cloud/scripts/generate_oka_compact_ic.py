#!/usr/bin/env python3
"""
Generate Oka et al. (2017) style compact cloud initial conditions.

This creates a Gaussian cloud matching the exact Oka paper parameters:
- Gaussian radial distribution with σ = 0.2 pc
- Isotropic velocity dispersion of 1.43 km/s
- Total cloud mass: 1000 M☉
- Cloud centered at origin (to be translated by config file)

Unlike the Lane-Emden hydrostatic sphere (R=1.13 pc), this is a compact
Gaussian cloud suitable for exact comparison with Oka et al. N-body results.

Usage:
    python generate_oka_compact_ic.py --output results/IC/oka_compact/snapshot_0000.csv
    python generate_oka_compact_ic.py --particles 1000  # Exact Oka (test particles)
    python generate_oka_compact_ic.py --particles 10000  # SPH resolution

Reference:
    Oka et al. (2017) Nature Astronomy, Methods section:
    "The cloud was assumed to be consisted of 1000 point masses of 1 M☉ each,
     with velocity dispersion 1.43 km s⁻¹ and Gaussian radial dispersion 0.2 pc."
"""

import numpy as np
import argparse
from pathlib import Path
from datetime import datetime


# =============================================================================
# OKA PAPER PARAMETERS
# =============================================================================

# Gaussian cloud parameters (from Oka et al. 2017 Methods)
SIGMA_RADIAL_PC = 0.2       # Gaussian radial dispersion [pc]
VELOCITY_DISPERSION_KMS = 1.43  # Isotropic velocity dispersion [km/s]
TOTAL_MASS_MSUN = 1000      # Total cloud mass [M☉]

# Default particle count (Oka used 1000, but SPH needs more for accuracy)
DEFAULT_PARTICLES = 10000

# Thermal state (cold molecular gas, ~10-50 K)
DEFAULT_TEMPERATURE_K = 20  # Assumed initial temperature


# =============================================================================
# UNIT SYSTEM (matching imbh_cloud simulations)
# =============================================================================
# Code units: length=1pc, mass=1000M☉, velocity=1km/s
# This gives G ≈ 1 naturally

# Physical constants in code units
GAMMA = 5/3  # Adiabatic index
MU = 2.4     # Mean molecular weight (molecular H2 with He)
K_B_CGS = 1.38e-16  # erg/K
M_H_CGS = 1.67e-24  # g

# Unit conversions
PC_TO_CM = 3.086e18
MSUN_TO_G = 1.989e33
KMS_TO_CMS = 1e5


def gaussian_3d_sample(n_particles: int, sigma: float, seed: int = None) -> np.ndarray:
    """
    Generate 3D positions from isotropic Gaussian distribution.

    Each coordinate is drawn from N(0, σ), giving a 3D Gaussian cloud.

    Args:
        n_particles: Number of particles
        sigma: Gaussian standard deviation (same in all directions)
        seed: Random seed for reproducibility

    Returns:
        positions: (n_particles, 3) array of positions
    """
    if seed is not None:
        np.random.seed(seed)

    # Draw each coordinate independently from N(0, σ)
    x = np.random.normal(0, sigma, n_particles)
    y = np.random.normal(0, sigma, n_particles)
    z = np.random.normal(0, sigma, n_particles)

    return np.column_stack([x, y, z])


def maxwell_boltzmann_velocities(n_particles: int, sigma_v: float, seed: int = None) -> np.ndarray:
    """
    Generate 3D velocities from Maxwell-Boltzmann distribution.

    Each velocity component is drawn from N(0, σ_v), giving isotropic dispersion.

    Args:
        n_particles: Number of particles
        sigma_v: 1D velocity dispersion (same in all directions)
        seed: Random seed for reproducibility

    Returns:
        velocities: (n_particles, 3) array of velocities
    """
    if seed is not None:
        np.random.seed(seed + 1000)  # Different seed from positions

    # Draw each component from N(0, σ_v)
    vx = np.random.normal(0, sigma_v, n_particles)
    vy = np.random.normal(0, sigma_v, n_particles)
    vz = np.random.normal(0, sigma_v, n_particles)

    # Remove bulk motion (center of mass velocity = 0)
    vx -= np.mean(vx)
    vy -= np.mean(vy)
    vz -= np.mean(vz)

    return np.column_stack([vx, vy, vz])


def estimate_local_density(positions: np.ndarray, masses: np.ndarray, n_neighbors: int = 50) -> np.ndarray:
    """
    Estimate local SPH density using nearest neighbor method.

    For Gaussian distribution, we use an analytical estimate:
    ρ(r) ≈ ρ_0 * exp(-r²/(2σ²)) where ρ_0 = M_total / (2πσ²)^(3/2)

    Args:
        positions: (N, 3) particle positions
        masses: (N,) particle masses
        n_neighbors: Number of neighbors for smoothing length estimate

    Returns:
        densities: (N,) local density estimates
    """
    n_particles = len(positions)
    total_mass = np.sum(masses)

    # Estimate sigma from positions
    sigma = np.std(positions)

    # Central density for Gaussian: ρ_0 = M / (2πσ²)^(3/2)
    rho_0 = total_mass / (2 * np.pi * sigma**2)**1.5

    # Density at each particle position
    r2 = np.sum(positions**2, axis=1)
    densities = rho_0 * np.exp(-r2 / (2 * sigma**2))

    # Ensure minimum density (avoid numerical issues)
    densities = np.maximum(densities, rho_0 * 1e-4)

    return densities


def estimate_smoothing_length(densities: np.ndarray, masses: np.ndarray, n_neighbors: int = 50, eta: float = 1.2) -> np.ndarray:
    """
    Estimate smoothing length from density.

    h = η * (m / ρ)^(1/3) * (n_neighbors)^(1/3)

    For standard SPH with 50 neighbors and η ≈ 1.2.
    """
    return eta * (n_neighbors * masses / densities)**(1/3)


def compute_thermal_state(densities: np.ndarray, temperature_K: float = 20) -> tuple:
    """
    Compute pressure and internal energy for fixed temperature.

    For ideal gas: P = ρ k_B T / (μ m_H)
                   u = P / ((γ-1) ρ) = k_B T / ((γ-1) μ m_H)

    In code units where velocity = 1 km/s:
    u [code] = u [km²/s²] = T [K] * k_B / (μ m_H) / (km/s)²
    """
    # Sound speed squared in cgs: c_s² = γ P / ρ = γ k_B T / (μ m_H)
    cs2_cgs = GAMMA * K_B_CGS * temperature_K / (MU * M_H_CGS)

    # Convert to code units (km/s)²
    cs2_code = cs2_cgs / KMS_TO_CMS**2

    # Internal energy: u = P / ((γ-1)ρ) = c_s² / (γ(γ-1))
    ene = cs2_code / (GAMMA * (GAMMA - 1))

    # Pressure: P = (γ-1) ρ u
    pres = (GAMMA - 1) * densities * ene

    # Sound speed
    sound = np.sqrt(cs2_code) * np.ones_like(densities)

    return pres, np.full_like(densities, ene), sound


def generate_oka_ic(
    n_particles: int = DEFAULT_PARTICLES,
    sigma_pc: float = SIGMA_RADIAL_PC,
    vel_disp_kms: float = VELOCITY_DISPERSION_KMS,
    total_mass_msun: float = TOTAL_MASS_MSUN,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    seed: int = 42
) -> dict:
    """
    Generate complete Oka-style compact cloud initial conditions.

    Returns a dictionary with all particle data ready for CSV output.
    """
    print("=" * 70)
    print("  Generating Oka et al. (2017) Compact Cloud IC")
    print("=" * 70)
    print()

    # Generate positions (Gaussian distribution)
    print(f"Generating {n_particles} particles with Gaussian distribution...")
    print(f"  σ_radial = {sigma_pc:.2f} pc")
    positions = gaussian_3d_sample(n_particles, sigma_pc, seed=seed)

    # Statistics
    r = np.sqrt(np.sum(positions**2, axis=1))
    print(f"  r_mean = {np.mean(r):.3f} pc")
    print(f"  r_max = {np.max(r):.3f} pc")
    print(f"  r_90% = {np.percentile(r, 90):.3f} pc")
    print()

    # Generate velocities (Maxwell-Boltzmann)
    print(f"Generating velocities with isotropic dispersion...")
    print(f"  σ_v = {vel_disp_kms:.2f} km/s")
    velocities = maxwell_boltzmann_velocities(n_particles, vel_disp_kms, seed=seed)

    # Velocity statistics
    v = np.sqrt(np.sum(velocities**2, axis=1))
    print(f"  |v|_mean = {np.mean(v):.3f} km/s")
    print(f"  |v|_max = {np.max(v):.3f} km/s")
    print()

    # Particle masses (equal mass particles)
    # Total mass in code units: 1000 M☉ = 1 code mass
    mass_code = total_mass_msun / 1000  # 1 code mass = 1000 M☉
    masses = np.full(n_particles, mass_code / n_particles)
    print(f"Particle masses:")
    print(f"  Total mass = {total_mass_msun:.0f} M☉ = {mass_code:.4f} code units")
    print(f"  Particle mass = {total_mass_msun/n_particles:.4f} M☉")
    print()

    # Estimate densities
    print("Computing density field...")
    densities = estimate_local_density(positions, masses)
    print(f"  ρ_center = {np.max(densities):.4e} code units")
    print(f"  ρ_edge = {np.min(densities):.4e} code units")
    print()

    # Compute smoothing lengths
    print("Computing smoothing lengths...")
    sml = estimate_smoothing_length(densities, masses)
    print(f"  h_mean = {np.mean(sml):.4f} pc")
    print(f"  h_range = [{np.min(sml):.4f}, {np.max(sml):.4f}] pc")
    print()

    # Compute thermal state
    print(f"Computing thermal state (T = {temperature_K} K)...")
    pres, ene, sound = compute_thermal_state(densities, temperature_K)
    print(f"  c_s = {sound[0]:.3f} km/s")
    print(f"  u = {ene[0]:.4e} (km/s)²")
    print()

    # Compile results
    data = {
        'id': np.arange(n_particles),
        'pos_x': positions[:, 0],
        'pos_y': positions[:, 1],
        'pos_z': positions[:, 2],
        'vel_x': velocities[:, 0],
        'vel_y': velocities[:, 1],
        'vel_z': velocities[:, 2],
        'acc_x': np.zeros(n_particles),
        'acc_y': np.zeros(n_particles),
        'acc_z': np.zeros(n_particles),
        'mass': masses,
        'dens': densities,
        'pres': pres,
        'ene': ene,
        'sml': sml,
        'sound': sound,
        'alpha': np.ones(n_particles),
        'balsara': np.ones(n_particles),
        'gradh': np.ones(n_particles),
        'phi': np.zeros(n_particles),
        'grav_acc_x': np.zeros(n_particles),
        'grav_acc_y': np.zeros(n_particles),
        'grav_acc_z': np.zeros(n_particles),
        'neighbor': np.full(n_particles, 50),
        'is_ghost': np.zeros(n_particles, dtype=int),
    }

    # Summary statistics
    print("=" * 70)
    print("  IC Summary")
    print("=" * 70)
    print(f"  Particles: {n_particles}")
    print(f"  Cloud size: σ = {sigma_pc:.2f} pc (vs Oka σ=0.2 pc)")
    print(f"  Velocity dispersion: {vel_disp_kms:.2f} km/s (Oka: 1.43 km/s)")
    print(f"  Total mass: {total_mass_msun:.0f} M☉")
    print(f"  Temperature: {temperature_K:.0f} K")
    print()

    return data


def write_csv(data: dict, output_path: Path, n_particles: int, sigma_pc: float, vel_disp_kms: float):
    """Write IC data to CSV file with SPH-compatible header."""

    timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")

    # Compute energies
    kinetic = 0.5 * np.sum(data['mass'] * (data['vel_x']**2 + data['vel_y']**2 + data['vel_z']**2))
    thermal = np.sum(data['mass'] * data['ene'])

    header = f"""# SPH Simulation Output - CSV Format v1
# Timestamp: {timestamp}
# Simulation: Oka et al. (2017) Compact Cloud IC
#
# === Unit System ===
# Type: IMBH_ENCOUNTER
# Length: 1 pc
# Mass: 1000 M☉
# Velocity: 1 km/s
# Time: 0.978 Myr
#
# === Oka Paper Parameters ===
# Gaussian radial dispersion: σ = {sigma_pc:.2f} pc
# Velocity dispersion: {vel_disp_kms:.2f} km/s (isotropic)
# Total mass: {np.sum(data['mass'])*1000:.0f} M☉
#
# === Simulation State ===
# Step: 0
# Time (code): 0
# Time (physical): 0 Myr
# Particle Count: {n_particles}
#
# === Physics Parameters ===
# Gamma: 1.666666667
# G: 1
# Neighbor Number: 50
# SPH Type: Godunov SPH
# Kernel: Wendland C4
# Gravity: enabled
# Balsara Switch: enabled
# Time-Dependent AV: disabled
#
# === Energy ===
# Kinetic: {kinetic:.6e} [code units]
# Thermal: {thermal:.6e} [code units]
# Potential: 0 [code units]
# Total: {kinetic + thermal:.6e} [code units]
#
# === Columns ===
# id: Particle ID
# pos_x, pos_y, pos_z: Position [pc]
# vel_x, vel_y, vel_z: Velocity [km/s]
# acc_x, acc_y, acc_z: Acceleration [code units]
# mass: Particle mass [1000 M☉]
# dens: Density [code units]
# pres: Pressure [code units]
# ene: Specific internal energy [(km/s)²]
# sml: Smoothing length [pc]
# sound: Sound speed [km/s]
# alpha: Artificial viscosity coefficient
# balsara: Balsara factor
# gradh: grad-h term factor
# phi: Gravitational potential [code units]
# grav_acc_x, grav_acc_y, grav_acc_z: Gravitational acceleration [code units]
# neighbor: Number of neighbors
# is_ghost: Ghost particle flag (0=real, 1=ghost)
#
"""

    # Column order
    columns = ['id', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z',
               'acc_x', 'acc_y', 'acc_z', 'mass', 'dens', 'pres', 'ene', 'sml',
               'sound', 'alpha', 'balsara', 'gradh', 'phi',
               'grav_acc_x', 'grav_acc_y', 'grav_acc_z', 'neighbor', 'is_ghost']

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write(header)
        f.write(','.join(columns) + '\n')

        for i in range(n_particles):
            row = []
            for col in columns:
                val = data[col][i]
                if col == 'id' or col == 'neighbor' or col == 'is_ghost':
                    row.append(str(int(val)))
                else:
                    row.append(f"{val}")
            f.write(','.join(row) + '\n')

    print(f"  Written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Oka et al. (2017) compact cloud initial conditions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_oka_compact_ic.py --particles 1000   # Exact Oka parameters
  python generate_oka_compact_ic.py --particles 10000  # Higher SPH resolution
  python generate_oka_compact_ic.py --sigma 0.1        # Even more compact cloud
        """
    )

    parser.add_argument('--output', '-o', type=str,
                        default='simulations/astrophysics/imbh_cloud/results/IC/oka_compact/snapshot_0000.csv',
                        help='Output CSV file path')
    parser.add_argument('--particles', '-n', type=int, default=DEFAULT_PARTICLES,
                        help=f'Number of particles (default: {DEFAULT_PARTICLES})')
    parser.add_argument('--sigma', type=float, default=SIGMA_RADIAL_PC,
                        help=f'Gaussian radial dispersion in pc (default: {SIGMA_RADIAL_PC})')
    parser.add_argument('--velocity-dispersion', type=float, default=VELOCITY_DISPERSION_KMS,
                        help=f'Velocity dispersion in km/s (default: {VELOCITY_DISPERSION_KMS})')
    parser.add_argument('--mass', type=float, default=TOTAL_MASS_MSUN,
                        help=f'Total cloud mass in M☉ (default: {TOTAL_MASS_MSUN})')
    parser.add_argument('--temperature', type=float, default=DEFAULT_TEMPERATURE_K,
                        help=f'Initial temperature in K (default: {DEFAULT_TEMPERATURE_K})')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')

    args = parser.parse_args()

    # Generate IC
    data = generate_oka_ic(
        n_particles=args.particles,
        sigma_pc=args.sigma,
        vel_disp_kms=args.velocity_dispersion,
        total_mass_msun=args.mass,
        temperature_K=args.temperature,
        seed=args.seed
    )

    # Write output
    output_path = Path(args.output)
    print(f"Writing IC to: {output_path}")
    write_csv(data, output_path, args.particles, args.sigma, args.velocity_dispersion)

    print()
    print("=" * 70)
    print("  Done! Use this IC with oka_compact.json config")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Create matching config file:")
    print("     simulations/astrophysics/imbh_cloud/config/presets/categories/CAT_OKA/A_61k/oka_compact.json")
    print("  2. Update resumeFromSnapshot to point to this IC")
    print("  3. Run simulation:")
    print("     make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT_OKA_A_oka_compact")


if __name__ == '__main__':
    main()
