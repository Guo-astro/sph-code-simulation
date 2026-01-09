#!/usr/bin/env python3
"""
Plot density profile from SPH simulation vs analytical solution.

Analytical profile: Modified isothermal sphere (BE-like)
    rho(r) = rho_c / (1 + (r/r_c)^2)

Physical parameters (from config):
    M_cloud = 40 M_sun
    T = 15 K
    n_center = 1000 cm^-3
    P_ext/k_B = 3000 K cm^-3
    R_cloud = 0.791 pc
    r_core = 0.396 pc
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Physical parameters from simulation
R_CLOUD = 0.791197  # pc
R_CORE = 0.395599   # pc
RHO_CENTER = 57.5846  # M_sun/pc^3 (code units)
RHO_EDGE = 11.5169    # M_sun/pc^3

# Number density conversion: rho [M_sun/pc^3] to n [cm^-3]
# n = rho / (mu * m_H) where mu = 2.33, m_H in M_sun
M_H_MSUN = 8.4e-58  # hydrogen mass in M_sun
MU = 2.33
RHO_TO_N = 1.0 / (MU * M_H_MSUN) * (3.086e18)**3 / 1e6  # approximate conversion


def analytical_profile(r, rho_c, r_c):
    """Modified isothermal sphere: rho = rho_c / (1 + (r/r_c)^2)"""
    return rho_c / (1 + (r / r_c)**2)


def load_snapshot(filepath):
    """Load snapshot CSV file, skipping header comments."""
    # Count comment lines
    with open(filepath, 'r') as f:
        skip_rows = 0
        for line in f:
            if line.startswith('#'):
                skip_rows += 1
            else:
                break

    df = pd.read_csv(filepath, skiprows=skip_rows)
    return df


def compute_radial_profile(df, n_bins=50):
    """Compute radial density profile from particle data."""
    # Calculate radius for each particle
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    rho = df['dens']

    # Filter out ghost particles if present
    if 'is_ghost' in df.columns:
        mask = df['is_ghost'] == 0
        r = r[mask]
        rho = rho[mask]

    # Create radial bins
    r_max = r.max() * 1.1
    bins = np.linspace(0, r_max, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # Compute mean density in each bin
    rho_mean = np.zeros(n_bins)
    rho_std = np.zeros(n_bins)
    counts = np.zeros(n_bins)

    for i in range(n_bins):
        mask = (r >= bins[i]) & (r < bins[i+1])
        if mask.sum() > 0:
            rho_mean[i] = rho[mask].mean()
            rho_std[i] = rho[mask].std()
            counts[i] = mask.sum()

    return bin_centers, rho_mean, rho_std, counts


def main():
    # Directories
    script_dir = Path(__file__).parent.parent
    relax_dir = script_dir / "results" / "be_sphere_relaxation"
    hydro_dir = script_dir / "results" / "be_sphere_hydrostatic"

    # Find available snapshots
    relax_snapshots = sorted(relax_dir.glob("snapshot_*.csv"))
    hydro_snapshots = sorted(hydro_dir.glob("snapshot_*.csv"))

    print(f"Found {len(relax_snapshots)} relaxation snapshots")
    print(f"Found {len(hydro_snapshots)} hydrostatic snapshots")

    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Analytical solution
    r_analytic = np.linspace(0.001, R_CLOUD * 1.2, 500)
    rho_analytic = analytical_profile(r_analytic, RHO_CENTER, R_CORE)

    # Plot 1: Relaxation results (initial, middle, final)
    ax1 = axes[0]
    if relax_snapshots:
        snapshots_to_plot = [
            relax_snapshots[0],
            relax_snapshots[len(relax_snapshots)//2],
            relax_snapshots[-1]
        ]
        colors = ['blue', 'orange', 'green']
        labels = ['Initial', 'Middle', 'Final']

        for snap, color, label in zip(snapshots_to_plot, colors, labels):
            df = load_snapshot(snap)
            r_bins, rho_mean, rho_std, counts = compute_radial_profile(df)

            # Filter bins with data
            mask = counts > 0
            ax1.scatter(r_bins[mask], rho_mean[mask], s=20, alpha=0.7,
                       color=color, label=f'{label} (SPH)')

    ax1.plot(r_analytic, rho_analytic, 'k-', lw=2, label='Analytical')
    ax1.axvline(R_CLOUD, color='gray', linestyle='--', alpha=0.5, label=f'R_cloud = {R_CLOUD:.2f} pc')
    ax1.axvline(R_CORE, color='gray', linestyle=':', alpha=0.5, label=f'r_core = {R_CORE:.2f} pc')

    ax1.set_xlabel('Radius [pc]', fontsize=12)
    ax1.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
    ax1.set_title('Phase 1: Relaxation', fontsize=14)
    ax1.legend(fontsize=9)
    ax1.set_xlim(0, R_CLOUD * 1.3)
    ax1.set_ylim(0, RHO_CENTER * 1.1)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Hydrostatic test results
    ax2 = axes[1]
    if hydro_snapshots:
        snapshots_to_plot = [
            hydro_snapshots[0],
            hydro_snapshots[len(hydro_snapshots)//2],
            hydro_snapshots[-1]
        ]
        colors = ['blue', 'orange', 'green']
        labels = ['t=0', 't=0.5 Myr', 't=1.0 Myr']

        for snap, color, label in zip(snapshots_to_plot, colors, labels):
            df = load_snapshot(snap)
            r_bins, rho_mean, rho_std, counts = compute_radial_profile(df)

            mask = counts > 0
            ax2.scatter(r_bins[mask], rho_mean[mask], s=20, alpha=0.7,
                       color=color, label=f'{label} (SPH)')

    ax2.plot(r_analytic, rho_analytic, 'k-', lw=2, label='Analytical')
    ax2.axvline(R_CLOUD, color='gray', linestyle='--', alpha=0.5, label=f'R_cloud = {R_CLOUD:.2f} pc')
    ax2.axvline(R_CORE, color='gray', linestyle=':', alpha=0.5, label=f'r_core = {R_CORE:.2f} pc')

    ax2.set_xlabel('Radius [pc]', fontsize=12)
    ax2.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
    ax2.set_title('Phase 2: Hydrostatic Test (with self-gravity)', fontsize=14)
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, R_CLOUD * 1.3)
    ax2.set_ylim(0, RHO_CENTER * 1.1)
    ax2.grid(True, alpha=0.3)

    plt.suptitle('Bonnor-Ebert Sphere: Density Profile Comparison\n' +
                 r'$\rho(r) = \rho_c / (1 + (r/r_c)^2)$, ' +
                 f'$\\rho_c$ = {RHO_CENTER:.1f} M$_\\odot$/pc$^3$, $r_c$ = {R_CORE:.2f} pc',
                 fontsize=12)

    plt.tight_layout()

    # Save figure
    output_path = script_dir / "results" / "density_profile_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_path}")

    # Also create a log-scale version
    fig2, ax3 = plt.subplots(figsize=(10, 7))

    # Final relaxation snapshot
    if relax_snapshots:
        df = load_snapshot(relax_snapshots[-1])
        r_bins, rho_mean, rho_std, counts = compute_radial_profile(df, n_bins=40)
        mask = counts > 0
        ax3.scatter(r_bins[mask], rho_mean[mask], s=40, alpha=0.8,
                   color='blue', label='Relaxed (Phase 1 final)', zorder=3)
        ax3.errorbar(r_bins[mask], rho_mean[mask], yerr=rho_std[mask],
                    fmt='none', color='blue', alpha=0.3, capsize=2)

    # Final hydrostatic snapshot
    if hydro_snapshots:
        df = load_snapshot(hydro_snapshots[-1])
        r_bins, rho_mean, rho_std, counts = compute_radial_profile(df, n_bins=40)
        mask = counts > 0
        ax3.scatter(r_bins[mask], rho_mean[mask], s=40, alpha=0.8,
                   color='green', marker='s', label='Hydrostatic (Phase 2 final)', zorder=3)

    ax3.plot(r_analytic, rho_analytic, 'k-', lw=2.5, label='Analytical: ' +
             r'$\rho = \rho_c/(1+(r/r_c)^2)$', zorder=2)

    ax3.axvline(R_CLOUD, color='red', linestyle='--', alpha=0.7, lw=1.5,
                label=f'$R_{{cloud}}$ = {R_CLOUD:.2f} pc')
    ax3.axvline(R_CORE, color='orange', linestyle=':', alpha=0.7, lw=1.5,
                label=f'$r_{{core}}$ = {R_CORE:.2f} pc')
    ax3.axhline(RHO_EDGE, color='purple', linestyle='-.', alpha=0.5, lw=1,
                label=f'$\\rho_{{edge}}$ = {RHO_EDGE:.1f} M$_\\odot$/pc$^3$')

    ax3.set_xlabel('Radius [pc]', fontsize=14)
    ax3.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=14)
    ax3.set_title('Bonnor-Ebert Sphere Density Profile\n' +
                  f'M = 40 M$_\\odot$, T = 15 K, $n_c$ = 1000 cm$^{{-3}}$',
                  fontsize=14)
    ax3.legend(fontsize=10, loc='upper right')
    ax3.set_xlim(0, R_CLOUD * 1.4)
    ax3.set_ylim(0, RHO_CENTER * 1.15)
    ax3.grid(True, alpha=0.3)

    # Add secondary y-axis for number density
    ax3_twin = ax3.twinx()
    n_min, n_max = 0, RHO_CENTER * 1.15 / (MU * 1.67e-24) * (2e33) / (3.086e18)**3
    # Approximate: n [cm^-3] ~ rho [M_sun/pc^3] * 17.4
    ax3_twin.set_ylim(0, RHO_CENTER * 1.15 * 17.4)
    ax3_twin.set_ylabel(r'Number Density [cm$^{-3}$]', fontsize=14)

    plt.tight_layout()

    output_path2 = script_dir / "results" / "density_profile_final.png"
    plt.savefig(output_path2, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path2}")


if __name__ == "__main__":
    main()
