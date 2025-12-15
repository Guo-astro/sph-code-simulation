#!/usr/bin/env python3
"""
Generate Initial Conditions for ISM-Equilibrium Cloud

Creates a self-gravitating cloud in hydrostatic equilibrium where:
- Pressure follows the ISM thermal equilibrium curve P(ρ) from Inoue & Inutsuka (2008)
- Each shell has T = T_eq(n) from the cooling/heating balance

This ensures the cloud remains stable when ISM cooling is applied.

Usage:
    python generate_ism_equilibrium_ic.py --n_particles 64000 --output cloud_ism_eq.csv
"""

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import argparse

# Physical constants (CGS)
k_B = 1.38e-16      # erg/K
m_p = 1.67e-24      # g
m_n = 1.27 * m_p    # mean particle mass
gamma = 5.0/3.0
Gamma_0 = 2e-26     # erg/s (I&I 2008 heating rate)

# Unit conversions (Galactic units)
density_to_n = 31.9      # code density -> n_H [cm^-3]
u_to_cgs = 1e10          # code energy -> erg/g
G_code = 1.0             # gravitational constant in code units


def lambda_over_gamma(T):
    """I&I 2008 cooling coefficient ratio Lambda/Gamma"""
    T = np.maximum(T, 10.0)
    term1 = 1e7 * np.exp(-114800.0 / (T + 1000.0))
    term2 = 0.014 * np.sqrt(T) * np.exp(-92.0 / T)
    return term1 + term2


def T_equilibrium(n_H):
    """Find equilibrium temperature for given number density n_H [cm^-3]"""
    # At equilibrium: Gamma = n * Lambda, so Lambda/Gamma = 1/n
    # Solve for T where lambda_over_gamma(T) = 1/n

    if np.isscalar(n_H):
        n_H = np.array([n_H])

    T_eq = np.zeros_like(n_H)

    for i, n in enumerate(n_H):
        if n <= 0:
            T_eq[i] = 1e5
            continue

        target = 1.0 / n

        # Binary search for T
        T_lo, T_hi = 10.0, 1e6
        for _ in range(100):
            T_mid = np.sqrt(T_lo * T_hi)
            ratio = lambda_over_gamma(T_mid)
            if ratio > target:
                T_lo = T_mid
            else:
                T_hi = T_mid
            if abs(T_hi - T_lo) / T_mid < 1e-6:
                break
        T_eq[i] = T_mid

    return T_eq if len(T_eq) > 1 else T_eq[0]


def u_equilibrium_code(rho_code):
    """Get equilibrium specific internal energy in code units"""
    n_H = rho_code * density_to_n
    T_eq = T_equilibrium(n_H)
    u_cgs = k_B * T_eq / ((gamma - 1) * m_n)
    return u_cgs / u_to_cgs


def P_equilibrium_code(rho_code):
    """Get equilibrium pressure in code units: P = rho * (gamma-1) * u"""
    u_eq = u_equilibrium_code(rho_code)
    return rho_code * (gamma - 1) * u_eq


def dPdrho_equilibrium(rho_code):
    """Derivative dP/drho for ISM equilibrium"""
    eps = rho_code * 1e-6
    P1 = P_equilibrium_code(rho_code - eps)
    P2 = P_equilibrium_code(rho_code + eps)
    return (P2 - P1) / (2 * eps)


def hydrostatic_ode(y, r, rho_c):
    """
    ODE for hydrostatic equilibrium with ISM P(rho)

    y = [rho, M]  where M is enclosed mass

    dP/dr = -rho * G * M / r^2
    dM/dr = 4 * pi * r^2 * rho

    Using chain rule: drho/dr = (dP/dr) / (dP/drho)
    """
    rho, M = y

    if rho <= 1e-10 or r <= 1e-10:
        return [0.0, 0.0]

    # Gravitational acceleration
    g = G_code * M / (r * r)

    # Pressure gradient from hydrostatic equilibrium
    dPdr = -rho * g

    # Convert to density gradient using EOS
    dPdrho = dPdrho_equilibrium(rho)
    if abs(dPdrho) < 1e-20:
        drhodr = 0.0
    else:
        drhodr = dPdr / dPdrho

    # Mass continuity
    dMdr = 4 * np.pi * r * r * rho

    return [drhodr, dMdr]


def generate_density_profile(rho_c, r_max=1.0, n_points=1000):
    """
    Integrate hydrostatic equilibrium to get rho(r) profile

    Args:
        rho_c: Central density (code units)
        r_max: Maximum radius to integrate to
        n_points: Number of radial points

    Returns:
        r_arr, rho_arr, M_arr
    """
    # Initial conditions at small r
    r0 = 1e-6
    M0 = 4/3 * np.pi * r0**3 * rho_c
    y0 = [rho_c, M0]

    # Radial grid
    r_arr = np.linspace(r0, r_max, n_points)

    # Integrate
    solution = odeint(hydrostatic_ode, y0, r_arr, args=(rho_c,))
    rho_arr = solution[:, 0]
    M_arr = solution[:, 1]

    # Ensure positivity
    rho_arr = np.maximum(rho_arr, 1e-10)

    return r_arr, rho_arr, M_arr


def sample_particles_from_profile(r_arr, rho_arr, n_particles):
    """
    Sample particle positions from density profile using rejection sampling
    """
    # Create interpolator for rho(r)
    rho_interp = interp1d(r_arr, rho_arr, kind='linear', fill_value=1e-10, bounds_error=False)

    # Find cloud radius (where density drops significantly)
    rho_edge = rho_arr[0] * 1e-4
    r_cloud = r_arr[rho_arr > rho_edge][-1] if any(rho_arr > rho_edge) else r_arr[-1]

    print(f"Cloud radius: {r_cloud:.4f}")
    print(f"Central density: {rho_arr[0]:.4f}")
    print(f"Edge density: {rho_arr[rho_arr > rho_edge][-1] if any(rho_arr > rho_edge) else rho_arr[-1]:.6f}")

    # Compute mass profile for importance sampling
    r_sample = np.linspace(0, r_cloud, 1000)
    rho_sample = rho_interp(r_sample)

    # CDF of mass: M(<r) = integral of 4*pi*r^2*rho dr
    dM = 4 * np.pi * r_sample**2 * rho_sample
    M_cumulative = np.cumsum(dM) * (r_sample[1] - r_sample[0])
    M_cumulative /= M_cumulative[-1]  # Normalize

    # Inverse CDF sampling for radii
    M_interp = interp1d(M_cumulative, r_sample, kind='linear')
    u = np.random.uniform(0, 1, n_particles)
    r_particles = M_interp(u)

    # Random angles for spherical distribution
    theta = np.arccos(2 * np.random.uniform(0, 1, n_particles) - 1)
    phi = np.random.uniform(0, 2*np.pi, n_particles)

    # Convert to Cartesian
    x = r_particles * np.sin(theta) * np.cos(phi)
    y = r_particles * np.sin(theta) * np.sin(phi)
    z = r_particles * np.cos(theta)

    # Get density at each particle position
    rho_particles = rho_interp(r_particles)

    return x, y, z, rho_particles


def main():
    parser = argparse.ArgumentParser(description='Generate ISM-equilibrium cloud IC')
    parser.add_argument('--n_particles', type=int, default=64000, help='Number of particles')
    parser.add_argument('--rho_c', type=float, default=1.0, help='Central density (code units)')
    parser.add_argument('--r_max', type=float, default=1.0, help='Maximum radius')
    parser.add_argument('--output', type=str, default='cloud_ism_eq.csv', help='Output file')
    parser.add_argument('--mass', type=float, default=1.0, help='Total cloud mass (code units)')
    args = parser.parse_args()

    print("="*60)
    print("Generating ISM-Equilibrium Cloud Initial Conditions")
    print("="*60)

    # Check equilibrium temperatures
    print("\nISM Equilibrium Curve (I&I 2008):")
    for n in [0.1, 0.5, 1, 5, 10, 30, 100]:
        T = T_equilibrium(n)
        print(f"  n_H = {n:5.1f} cm^-3: T_eq = {T:7.0f} K")

    # Generate density profile
    print(f"\nGenerating density profile (rho_c = {args.rho_c})...")
    r_arr, rho_arr, M_arr = generate_density_profile(args.rho_c, args.r_max)

    # Sample particles
    print(f"\nSampling {args.n_particles} particles...")
    x, y, z, rho = sample_particles_from_profile(r_arr, rho_arr, args.n_particles)

    # Compute particle mass (equal mass particles)
    total_mass = args.mass
    particle_mass = total_mass / args.n_particles

    # Get equilibrium energy for each particle
    print("Computing equilibrium energies...")
    u_eq = u_equilibrium_code(rho)

    # Pressure
    pres = rho * (gamma - 1) * u_eq

    # Check temperatures
    n_H = rho * density_to_n
    T = T_equilibrium(n_H)
    print(f"\nTemperature range: {T.min():.0f} - {T.max():.0f} K")
    print(f"Density range: {rho.min():.4f} - {rho.max():.4f} (code)")
    print(f"n_H range: {n_H.min():.2f} - {n_H.max():.2f} cm^-3")

    # Create DataFrame
    df = pd.DataFrame({
        'id': np.arange(args.n_particles),
        'pos_x': x,
        'pos_y': y,
        'pos_z': z,
        'vel_x': np.zeros(args.n_particles),
        'vel_y': np.zeros(args.n_particles),
        'vel_z': np.zeros(args.n_particles),
        'acc_x': np.zeros(args.n_particles),
        'acc_y': np.zeros(args.n_particles),
        'acc_z': np.zeros(args.n_particles),
        'mass': np.full(args.n_particles, particle_mass),
        'dens': rho,
        'pres': pres,
        'ene': u_eq,
        'sml': np.zeros(args.n_particles),  # Will be computed by SPH
        'sound': np.sqrt(gamma * pres / rho),
        'gradh': np.ones(args.n_particles),
        'omega': np.ones(args.n_particles),
        'grav_acc_x': np.zeros(args.n_particles),
        'grav_acc_y': np.zeros(args.n_particles),
        'grav_acc_z': np.zeros(args.n_particles),
        'grav_pot': np.zeros(args.n_particles),
        'divv': np.zeros(args.n_particles),
        'neighbors': np.zeros(args.n_particles, dtype=int),
        'type': np.zeros(args.n_particles, dtype=int),
    })

    # Write header and data
    output_path = args.output
    with open(output_path, 'w') as f:
        f.write("# SPH Simulation Output - CSV Format v1\n")
        f.write(f"# ISM-Equilibrium Cloud IC\n")
        f.write(f"# n_particles: {args.n_particles}\n")
        f.write(f"# rho_c: {args.rho_c}\n")
        f.write(f"# total_mass: {total_mass}\n")
        f.write(f"# T_range: {T.min():.0f}-{T.max():.0f} K\n")
        f.write("#\n")
        f.write("# === Unit System ===\n")
        f.write("# type: IMBH_ENCOUNTER\n")
        f.write("# length_pc: 1.0\n")
        f.write("# mass_1e3Msun: 1.0\n")
        f.write("# velocity_kms: 1.0\n")
        f.write("#\n")

    df.to_csv(output_path, mode='a', index=False)
    print(f"\nWrote {args.n_particles} particles to {output_path}")

    print("\n" + "="*60)
    print("To use this IC with ISM cooling, create a config with:")
    print('  "resumeFromSnapshot": "<path_to_this_file>"')
    print('  "enableCooling": true')
    print("="*60)


if __name__ == '__main__':
    main()
