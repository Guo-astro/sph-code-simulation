#!/usr/bin/env python3
"""
Generate initial conditions for BNS cocoon shock breakout simulations.

This script creates CSV files with SPH particle data that can be loaded
by the SRGSPH simulation code.

The IC represents:
1. Inner cocoon: Hot, relativistic plasma from jet interaction (Γ ~ 2-10)
2. Outer ejecta: Cold, mildly relativistic BNS merger ejecta (v ~ 0.1-0.3c)
3. Transition: Shock interface between cocoon and ejecta
"""

import os
import sys
import argparse
from pathlib import Path
from dataclasses import dataclass
import numpy as np


@dataclass
class CocoonProfile:
    """Physical parameters for cocoon+ejecta system."""
    
    # Cocoon parameters
    E_cocoon: float = 1e49        # Cocoon energy [erg]
    M_cocoon: float = 1e-5        # Cocoon mass [Msun]
    gamma_cocoon: float = 5.0     # Cocoon Lorentz factor
    r_cocoon: float = 0.1         # Cocoon radius (code units, fraction of domain)
    
    # Ejecta parameters
    M_ejecta: float = 1e-3        # Ejecta mass [Msun]
    v_min: float = 0.05           # Minimum ejecta velocity (c)
    v_max: float = 0.30           # Maximum ejecta velocity (c)
    power_law: float = -2.0       # dM/dv ∝ v^power_law
    
    # Numerical parameters
    n_particles: int = 10000
    domain_size: float = 1.0
    gamma_eos: float = 4.0/3.0    # Relativistic gas γ = 4/3
    
    # Code units: c = 1, G = 1
    c: float = 1.0


def generate_cocoon_ejecta_1d(profile: CocoonProfile) -> dict:
    """
    Generate 1D cocoon + ejecta initial conditions.
    
    Parameters
    ----------
    profile : CocoonProfile
        Physical parameters for the system.
    
    Returns
    -------
    dict with particle arrays:
        pos, vel, dens, pres, mass, sml
    """
    n = profile.n_particles
    n_cocoon = int(0.2 * n)  # 20% in cocoon
    n_ejecta = n - n_cocoon
    
    # === Cocoon particles (inner, hot, fast) ===
    r_cocoon = np.linspace(0.001, profile.r_cocoon, n_cocoon)
    
    # Cocoon velocity from Lorentz factor (roughly uniform)
    v_cocoon = profile.c * np.sqrt(1 - 1/profile.gamma_cocoon**2)
    vel_cocoon = v_cocoon * np.ones(n_cocoon)
    
    # Hot cocoon: high pressure, moderate density
    # P ~ E / V for relativistic plasma
    V_cocoon = 4/3 * np.pi * profile.r_cocoon**3
    P_cocoon = profile.E_cocoon / (3 * V_cocoon)  # P ~ E/(3V) for relativistic
    
    # Scale to code units (normalize)
    P_cocoon_code = 1.0  # High pressure normalized
    
    # Density from mass conservation
    rho_cocoon = profile.M_cocoon / V_cocoon
    rho_cocoon_code = 1.0  # Normalized
    
    dens_cocoon = rho_cocoon_code * np.ones(n_cocoon)
    pres_cocoon = P_cocoon_code * np.ones(n_cocoon)
    
    # === Ejecta particles (outer, cold, slow) ===
    # Velocity follows v_min to v_max in shell structure
    r_inner = profile.r_cocoon
    r_outer = profile.domain_size
    
    # Use homologous expansion: r = v * t0
    # We set up at t0 = r_outer / v_max
    t0 = r_outer / profile.v_max
    
    # Velocity distribution: power law
    # dM/dv ∝ v^α => M(>v) ∝ v^(α+1)
    u_random = np.linspace(0, 1, n_ejecta)  # Uniform in mass shells
    
    # Invert CDF for power law
    alpha = profile.power_law
    if alpha != -1:
        # M(<v) = M_total * [(v/v_min)^(α+1) - 1] / [(v_max/v_min)^(α+1) - 1]
        v_ratio = profile.v_max / profile.v_min
        v_scale = (1 + u_random * (v_ratio**(alpha+1) - 1)) ** (1/(alpha+1))
        vel_ejecta = profile.v_min * v_scale
    else:
        # Logarithmic distribution
        vel_ejecta = profile.v_min * np.exp(u_random * np.log(profile.v_max/profile.v_min))
    
    # Position from velocity (homologous)
    r_ejecta = vel_ejecta * t0
    
    # Shift to be outside cocoon
    r_ejecta = r_inner + (r_ejecta - r_ejecta.min()) * (r_outer - r_inner) / (r_ejecta.max() - r_ejecta.min() + 1e-10)
    
    # Ejecta density: steep power law
    # ρ ∝ r^-3 roughly for expanding ejecta
    dens_ejecta = 0.01 * (r_inner / r_ejecta) ** 3
    dens_ejecta = np.clip(dens_ejecta, 1e-8, 1.0)
    
    # Cold ejecta: low pressure
    pres_ejecta = 1e-4 * dens_ejecta  # P << ρc²
    
    # === Combine ===
    pos = np.concatenate([r_cocoon, r_ejecta])
    vel = np.concatenate([vel_cocoon, vel_ejecta])
    dens = np.concatenate([dens_cocoon, dens_ejecta])
    pres = np.concatenate([pres_cocoon, pres_ejecta])
    
    # Sort by position
    idx = np.argsort(pos)
    pos = pos[idx]
    vel = vel[idx]
    dens = dens[idx]
    pres = pres[idx]
    
    # Mass per particle (equal mass SPH)
    # In 1D spherical: dV = 4πr²dr
    dr = np.diff(pos, prepend=pos[0]/2, append=pos[-1]*1.1)
    dr = 0.5 * (dr[:-1] + dr[1:])  # Average spacing
    
    # Volume element
    vol = 4 * np.pi * pos**2 * dr
    mass = dens * vol
    
    # Normalize to total mass
    mass = mass / mass.sum() * (profile.M_cocoon + profile.M_ejecta)
    
    # Smoothing length from particle spacing
    sml = 2.0 * dr  # h ~ 2 * particle spacing
    
    return {
        'pos': pos,
        'vel': vel,
        'dens': dens,
        'pres': pres,
        'mass': mass,
        'sml': sml
    }


def compute_sr_quantities(pos, vel, dens, pres, mass, gamma_eos=4.0/3.0, c=1.0):
    """
    Compute special relativistic quantities from primitives.
    
    Parameters
    ----------
    pos, vel, dens, pres, mass : array
        Primitive variables.
    gamma_eos : float
        Adiabatic index.
    c : float
        Speed of light.
    
    Returns
    -------
    dict with SR quantities:
        gamma_lor, N, S, e, enthalpy, nu, ene
    """
    # Lorentz factor
    v = np.abs(vel)
    gamma_lor = 1.0 / np.sqrt(1 - (v/c)**2)
    
    # Lab-frame density: N = γ * n (where n = rest-frame density)
    # We interpret input dens as rest-frame density
    n = dens
    N = gamma_lor * n
    
    # Specific internal energy: u = P / [(γ-1) * n]
    u = pres / ((gamma_eos - 1) * n + 1e-30)
    
    # Specific enthalpy: H = 1 + u/c² + P/(n*c²)
    enthalpy = 1.0 + u/c**2 + pres/(n * c**2 + 1e-30)
    
    # Canonical momentum per baryon: S = γ * H * v
    S = gamma_lor * enthalpy * vel
    
    # Canonical energy per baryon: e = γ*H - P/(N*c²)
    e = gamma_lor * enthalpy - pres / (N * c**2 + 1e-30)
    
    # Baryon number per particle: ν = m_particle (constant)
    nu = mass
    
    return {
        'gamma_lor': gamma_lor,
        'N': N,
        'S': S,
        'e': e,
        'enthalpy': enthalpy,
        'nu': nu,
        'ene': u
    }


def write_ic_csv(filepath: Path, data: dict, sr_data: dict):
    """
    Write initial conditions to CSV file.
    
    Format compatible with SRGSPH code.
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    n = len(data['pos'])
    
    with open(filepath, 'w') as f:
        # Header
        f.write("# BNS Cocoon Initial Conditions\n")
        f.write("# Generated by generate_cocoon_ic.py\n")
        f.write(f"# N = {n} particles\n")
        f.write("#\n")
        
        # Column headers
        columns = [
            'id', 'pos_x', 'vel_x', 'mass', 'dens', 'pres', 'ene', 'sml',
            'gamma_lor', 'N', 'S_x', 'e_canon', 'enthalpy', 'nu'
        ]
        f.write(','.join(columns) + '\n')
        
        # Data
        for i in range(n):
            row = [
                i,
                data['pos'][i],
                data['vel'][i],
                data['mass'][i],
                data['dens'][i],
                data['pres'][i],
                sr_data['ene'][i],
                data['sml'][i],
                sr_data['gamma_lor'][i],
                sr_data['N'][i],
                sr_data['S'][i],
                sr_data['e'][i],
                sr_data['enthalpy'][i],
                sr_data['nu'][i]
            ]
            f.write(','.join(f'{x:.16e}' for x in row) + '\n')
    
    print(f"Wrote {n} particles to {filepath}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate BNS cocoon initial conditions for SRGSPH"
    )
    
    parser.add_argument(
        "-o", "--output", type=Path,
        default=Path("simulations/astrophysics/bns_cocoon/data/ic/cocoon_1d_ic.csv"),
        help="Output CSV file path"
    )
    parser.add_argument(
        "-n", "--n-particles", type=int, default=10000,
        help="Number of SPH particles"
    )
    parser.add_argument(
        "--gamma-cocoon", type=float, default=5.0,
        help="Cocoon Lorentz factor"
    )
    parser.add_argument(
        "--r-cocoon", type=float, default=0.1,
        help="Cocoon radius (fraction of domain)"
    )
    parser.add_argument(
        "--v-min", type=float, default=0.05,
        help="Minimum ejecta velocity (units of c)"
    )
    parser.add_argument(
        "--v-max", type=float, default=0.30,
        help="Maximum ejecta velocity (units of c)"
    )
    parser.add_argument(
        "--gamma-eos", type=float, default=1.333333,
        help="Adiabatic index (4/3 for relativistic gas)"
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Plot the initial conditions"
    )
    
    args = parser.parse_args()
    
    # Create profile
    profile = CocoonProfile(
        n_particles=args.n_particles,
        gamma_cocoon=args.gamma_cocoon,
        r_cocoon=args.r_cocoon,
        v_min=args.v_min,
        v_max=args.v_max,
        gamma_eos=args.gamma_eos
    )
    
    print("=" * 60)
    print("BNS Cocoon IC Generator")
    print("=" * 60)
    print(f"N particles:      {profile.n_particles}")
    print(f"Cocoon Γ:         {profile.gamma_cocoon}")
    print(f"Cocoon radius:    {profile.r_cocoon}")
    print(f"Ejecta v_min:     {profile.v_min}c")
    print(f"Ejecta v_max:     {profile.v_max}c")
    print(f"EOS γ:            {profile.gamma_eos}")
    print("=" * 60)
    
    # Generate primitive variables
    data = generate_cocoon_ejecta_1d(profile)
    
    # Compute SR quantities
    sr_data = compute_sr_quantities(
        data['pos'], data['vel'], data['dens'], data['pres'], data['mass'],
        gamma_eos=profile.gamma_eos, c=profile.c
    )
    
    # Write output
    write_ic_csv(args.output, data, sr_data)
    
    # Plot if requested
    if args.plot:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        ax = axes[0, 0]
        ax.semilogy(data['pos'], data['dens'], 'b-', lw=1.5)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Density')
        ax.set_title('Density Profile')
        ax.grid(True, alpha=0.3)
        
        ax = axes[0, 1]
        ax.semilogy(data['pos'], data['pres'], 'r-', lw=1.5)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Pressure')
        ax.set_title('Pressure Profile')
        ax.grid(True, alpha=0.3)
        
        ax = axes[1, 0]
        ax.plot(data['pos'], data['vel'], 'g-', lw=1.5)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Velocity [c]')
        ax.set_title('Velocity Profile')
        ax.grid(True, alpha=0.3)
        
        ax = axes[1, 1]
        ax.semilogy(data['pos'], sr_data['gamma_lor'], 'm-', lw=1.5)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Lorentz factor Γ')
        ax.set_title('Lorentz Factor Profile')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        plot_file = args.output.with_suffix('.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {plot_file}")
        
        plt.show()


if __name__ == "__main__":
    main()
