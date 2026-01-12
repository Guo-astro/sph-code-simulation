#!/usr/bin/env python3
"""
Plot radial profiles of temperature, density, and cooling rate
for IMBH cloud Phase 2.5 simulation with ISM cooling.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Physical constants
k_B = 1.38e-16      # Boltzmann constant [erg/K]
m_p = 1.67e-24      # Proton mass [g]
mu = 1.27           # Mean molecular weight for neutral H
gamma = 5/3         # Adiabatic index

# Unit conversions from config
density_to_n_H = 31.86    # cm^-3 per code density
u_to_cgs = 1e10           # erg/g per code energy
t_to_cgs = 3.086e13       # s per code time

# Inoue & Inutsuka (2008) cooling parameters
GAMMA_HEAT = 2.0e-26      # erg s^-1 (heating rate per particle)

def temperature_from_u(u_code):
    """Convert specific internal energy to temperature [K]."""
    u_cgs = u_code * u_to_cgs
    T = u_cgs * mu * m_p * (gamma - 1) / k_B
    return T

def cooling_rate(n_H, T):
    """
    Inoue & Inutsuka (2008) cooling function.
    Returns net cooling rate per volume [erg cm^-3 s^-1].
    Positive = cooling, negative = heating.
    """
    # Lambda/Gamma ratio (Eq. 8-9)
    Lambda_over_Gamma = (1e7 * np.exp(-114800 / (T + 1000)) + 
                         0.014 * np.sqrt(T) * np.exp(-92 / T))
    
    # Net cooling rate: rho_n * L = n_n * (-Gamma + n_n * Lambda)
    # = n_H * (-Gamma + n_H * Gamma * Lambda/Gamma)
    # = n_H * Gamma * (-1 + n_H * Lambda/Gamma)
    net_cooling = n_H * GAMMA_HEAT * (-1 + n_H * Lambda_over_Gamma)
    
    return net_cooling

def load_snapshot(filepath):
    """Load snapshot CSV, skipping comment lines."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('id,'):
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError(f"Could not find header in {filepath}")
    
    df = pd.read_csv(filepath, skiprows=header_idx)
    return df

def compute_radial_profiles(df, n_bins=50):
    """Compute radial profiles from particle data."""
    # Filter out ghost particles
    real = df[df['is_ghost'] == 0].copy()
    
    # Compute radius (2D for cylindrical cloud)
    real['r'] = np.sqrt(real['pos_x']**2 + real['pos_y']**2)
    
    # Compute derived quantities
    # Internal energy: u = pres / (dens * (gamma-1))
    real['u'] = real['pres'] / (real['dens'] * (gamma - 1))
    real['T'] = temperature_from_u(real['u'])
    real['n_H'] = real['dens'] * density_to_n_H
    real['cool_rate'] = cooling_rate(real['n_H'], real['T'])
    
    # Radial binning
    r_max = real['r'].max()
    r_bins = np.linspace(0, r_max, n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    # Compute mean values in each bin
    profiles = {
        'r': r_centers,
        'dens': np.zeros(n_bins),
        'T': np.zeros(n_bins),
        'n_H': np.zeros(n_bins),
        'cool_rate': np.zeros(n_bins),
        'count': np.zeros(n_bins)
    }
    
    real['r_bin'] = np.digitize(real['r'], r_bins) - 1
    real['r_bin'] = real['r_bin'].clip(0, n_bins - 1)
    
    for i in range(n_bins):
        mask = real['r_bin'] == i
        if mask.sum() > 0:
            profiles['dens'][i] = real.loc[mask, 'dens'].mean()
            profiles['T'][i] = real.loc[mask, 'T'].mean()
            profiles['n_H'][i] = real.loc[mask, 'n_H'].mean()
            profiles['cool_rate'][i] = real.loc[mask, 'cool_rate'].mean()
            profiles['count'][i] = mask.sum()
    
    return profiles, real

def plot_profiles(profiles, real, output_path, snapshot_name):
    """Create multi-panel plot of radial profiles."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    r = profiles['r']
    
    # Panel 1: Density
    ax1 = axes[0, 0]
    ax1.semilogy(r, profiles['dens'], 'b-', lw=2, label='Code density')
    ax1.set_xlabel('Radius [code units]')
    ax1.set_ylabel('Density [code units]')
    ax1.set_title('Radial Density Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Add n_H on secondary axis
    ax1b = ax1.twinx()
    ax1b.semilogy(r, profiles['n_H'], 'g--', lw=2, alpha=0.7, label='n_H')
    ax1b.set_ylabel('n_H [cm⁻³]', color='g')
    ax1b.tick_params(axis='y', labelcolor='g')
    
    # Panel 2: Temperature
    ax2 = axes[0, 1]
    ax2.semilogy(r, profiles['T'], 'r-', lw=2)
    ax2.axhline(y=6000, color='orange', ls='--', label='WNM equilibrium (~6000 K)')
    ax2.axhline(y=100, color='cyan', ls='--', label='CNM equilibrium (~100 K)')
    ax2.set_xlabel('Radius [code units]')
    ax2.set_ylabel('Temperature [K]')
    ax2.set_title('Radial Temperature Profile')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Panel 3: Cooling rate
    ax3 = axes[1, 0]
    cool = profiles['cool_rate']
    # Plot positive (cooling) and negative (heating) separately
    pos_mask = cool > 0
    neg_mask = cool < 0
    
    if pos_mask.any():
        ax3.semilogy(r[pos_mask], cool[pos_mask], 'b-', lw=2, label='Cooling')
    if neg_mask.any():
        ax3.semilogy(r[neg_mask], -cool[neg_mask], 'r--', lw=2, label='Heating')
    
    ax3.set_xlabel('Radius [code units]')
    ax3.set_ylabel('|Cooling rate| [erg cm⁻³ s⁻¹]')
    ax3.set_title('Radial Cooling/Heating Rate')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Panel 4: Scatter plot T vs n_H with equilibrium curve
    ax4 = axes[1, 1]
    
    # Sample particles for scatter plot
    sample = real.sample(min(5000, len(real)))
    scatter = ax4.scatter(sample['n_H'], sample['T'], c=sample['r'], 
                          cmap='viridis', s=1, alpha=0.5)
    plt.colorbar(scatter, ax=ax4, label='Radius')
    
    # Plot equilibrium curve
    n_eq = np.logspace(-1, 2.5, 100)
    T_eq = []
    for n in n_eq:
        # Find T where cooling = 0
        # -1 + n * Lambda/Gamma = 0 => Lambda/Gamma = 1/n
        # Solve numerically
        from scipy.optimize import brentq
        def eq_func(T):
            LG = 1e7 * np.exp(-114800/(T+1000)) + 0.014*np.sqrt(T)*np.exp(-92/T)
            return LG - 1/n
        try:
            T_sol = brentq(eq_func, 10, 20000)
            T_eq.append(T_sol)
        except:
            T_eq.append(np.nan)
    
    ax4.loglog(n_eq, T_eq, 'k-', lw=2, label='Thermal equilibrium')
    ax4.set_xlabel('n_H [cm⁻³]')
    ax4.set_ylabel('Temperature [K]')
    ax4.set_title('Phase Diagram (T vs n_H)')
    ax4.set_xlim(0.1, 300)
    ax4.set_ylim(1, 20000)
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    # Summary stats
    T_mean = real['T'].mean()
    T_min, T_max = real['T'].min(), real['T'].max()
    n_mean = real['n_H'].mean()
    
    fig.suptitle(f'Phase 2.5 ISM Cooling - {snapshot_name}\n'
                 f'T: {T_min:.1f} - {T_max:.1f} K (mean: {T_mean:.1f} K), '
                 f'n_H mean: {n_mean:.1f} cm⁻³', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()
    
    return T_mean, T_min, T_max, n_mean

def main():
    results_dir = Path('/Users/guo-opt-p148/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/phase2.5_cooling_gamma53')
    
    # Find latest snapshot
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    
    if not snapshots:
        print("No snapshots found!")
        sys.exit(1)
    
    latest = snapshots[-1]
    print(f"Loading: {latest.name}")
    
    df = load_snapshot(latest)
    print(f"Loaded {len(df)} particles ({(df['is_ghost']==0).sum()} real)")
    
    profiles, real = compute_radial_profiles(df)
    
    output_path = results_dir / 'radial_profiles.png'
    T_mean, T_min, T_max, n_mean = plot_profiles(profiles, real, output_path, latest.name)
    
    print(f"\nSummary:")
    print(f"  Temperature: {T_min:.1f} - {T_max:.1f} K (mean: {T_mean:.1f} K)")
    print(f"  n_H mean: {n_mean:.1f} cm⁻³")
    print(f"  Cooling regime: {'Net cooling' if T_mean > 1000 else 'Near equilibrium'}")

if __name__ == '__main__':
    main()
