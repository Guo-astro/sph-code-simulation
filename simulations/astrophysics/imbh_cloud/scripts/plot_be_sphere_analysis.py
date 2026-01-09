#!/usr/bin/env python3
"""
Comprehensive Bonnor-Ebert sphere analysis with TRUE analytical profile overlay.

Solves the isothermal Lane-Emden equation:
    (1/xi^2) d/dxi(xi^2 dpsi/dxi) = exp(-psi)

where:
    xi = r / r_0        (dimensionless radius)
    r_0 = c_s / sqrt(4*pi*G*rho_c)
    psi = -ln(rho/rho_c)
    rho = rho_c * exp(-psi)

Option A Parameters (Atomic CNM on K&I 2000 curve):
    mu = 1.27 (atomic HI + He)
    T = 9 K
    n_center = 980 cm^-3
    n_edge = 70 cm^-3
    P_ext = 630 K cm^-3
    R_cloud = 1.2 pc
    xi_s = 6.45 (critical BE sphere)
    M_cloud = 40 M_sun
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp
import json
import argparse

# Physical constants
K_B_CGS = 1.380649e-16       # erg/K
M_PROTON_CGS = 1.6726219e-24 # g
MSUN_CGS = 1.989e33          # g
PC_CGS = 3.086e18            # cm
KMS_CGS = 1.0e5              # cm/s
G_CGS = 6.674e-8             # cm^3/g/s^2

# Option A default parameters
DEFAULT_PARAMS = {
    'mu': 1.27,
    'T_cloud': 9.0,
    'n_center': 980.0,
    'xi_s': 6.45,
    'M_cloud': 40.0,
}


def solve_lane_emden(xi_max, n_points=1000):
    """
    Solve isothermal Lane-Emden equation using RK45.

    Returns:
        xi: Dimensionless radius array
        psi: Dimensionless potential array
        dpsi: dpsi/dxi array
    """
    def equations(xi, y):
        psi, dpsi = y
        if xi < 1e-10:
            return [0, 0]
        d2psi = np.exp(-psi) - 2 * dpsi / xi
        return [dpsi, d2psi]

    # Initial conditions from series expansion
    xi0 = 1e-6
    psi0 = xi0**2 / 6
    dpsi0 = xi0 / 3

    xi_eval = np.linspace(xi0, xi_max, n_points)

    sol = solve_ivp(equations, [xi0, xi_max], [psi0, dpsi0],
                    t_eval=xi_eval, method='RK45', dense_output=True)

    return sol.t, sol.y[0], sol.y[1]


def compute_be_profile(params, n_points=500):
    """
    Compute physical BE profile from parameters.

    Returns:
        r: Radius array [pc]
        rho: Density array [M_sun/pc^3]
        n: Number density array [cm^-3]
        M_enc: Enclosed mass array [M_sun]
    """
    mu = params.get('mu', DEFAULT_PARAMS['mu'])
    T = params.get('T_cloud', DEFAULT_PARAMS['T_cloud'])
    n_center = params.get('n_center', DEFAULT_PARAMS['n_center'])
    xi_s = params.get('xi_s', DEFAULT_PARAMS['xi_s'])

    # Sound speed
    c_s_cgs = np.sqrt(K_B_CGS * T / (mu * M_PROTON_CGS))  # cm/s
    c_s = c_s_cgs / KMS_CGS  # km/s

    # Central density: n_center [cm^-3] -> rho_c [M_sun/pc^3]
    rho_c_cgs = n_center * mu * M_PROTON_CGS  # g/cm^3
    rho_c = rho_c_cgs * PC_CGS**3 / MSUN_CGS  # M_sun/pc^3

    # Scale length: r_0 = c_s / sqrt(4*pi*G*rho_c)
    # In code units: [pc, M_sun, km/s, Myr]
    G_code = 0.00430091  # pc^3 / (M_sun * Myr^2)
    r_0 = c_s / np.sqrt(4 * np.pi * G_code * rho_c)  # pc

    # Solve Lane-Emden
    xi, psi, dpsi = solve_lane_emden(xi_s * 1.1, n_points)

    # Physical quantities
    r = xi * r_0
    rho = rho_c * np.exp(-psi)

    # Number density
    density_to_n = (MSUN_CGS / PC_CGS**3) / (mu * M_PROTON_CGS)
    n = rho * density_to_n

    # Enclosed mass: M(r) = 4*pi*rho_c*r_0^3 * xi^2 * dpsi
    M_enc = 4 * np.pi * rho_c * r_0**3 * xi**2 * dpsi

    # Truncate at xi_s
    mask = xi <= xi_s

    return {
        'r': r[mask],
        'rho': rho[mask],
        'n': n[mask],
        'M_enc': M_enc[mask],
        'r_0': r_0,
        'R_cloud': xi_s * r_0,
        'rho_c': rho_c,
        'n_center': n_center,
        'c_s': c_s,
        'xi_s': xi_s,
        'mu': mu,
        'T': T,
    }


def load_snapshot(filepath):
    """Load snapshot CSV file, skipping header comments."""
    with open(filepath, 'r') as f:
        skip_rows = 0
        for line in f:
            if line.startswith('#'):
                skip_rows += 1
            else:
                break

    df = pd.read_csv(filepath, skiprows=skip_rows)
    return df


def compute_radial_profile(df, n_bins=50, exclude_ghosts=True):
    """Compute radial density profile from particle data."""
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    rho = df['dens']

    if exclude_ghosts and 'is_ghost' in df.columns:
        mask = df['is_ghost'] == 0
        r = r[mask]
        rho = rho[mask]

    r_max = r.max() * 1.05
    bins = np.linspace(0, r_max, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    rho_mean = np.zeros(n_bins)
    rho_std = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        mask = (r >= bins[i]) & (r < bins[i+1])
        if mask.sum() > 0:
            rho_mean[i] = rho[mask].mean()
            rho_std[i] = rho[mask].std()
            counts[i] = mask.sum()

    return bin_centers, rho_mean, rho_std, counts


def compute_velocity_dispersion(df, exclude_ghosts=True):
    """Compute radial velocity dispersion profile."""
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    vx, vy, vz = df['vel_x'], df['vel_y'], df['vel_z']

    if exclude_ghosts and 'is_ghost' in df.columns:
        mask = df['is_ghost'] == 0
        r, vx, vy, vz = r[mask], vx[mask], vy[mask], vz[mask]

    # Total velocity magnitude
    v_tot = np.sqrt(vx**2 + vy**2 + vz**2)

    # Global velocity dispersion
    sigma_global = np.std(v_tot)

    # Radial profile
    n_bins = 20
    r_max = r.max() * 1.05
    bins = np.linspace(0, r_max, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    sigma_r = np.zeros(n_bins)

    for i in range(n_bins):
        mask = (r >= bins[i]) & (r < bins[i+1])
        if mask.sum() > 5:
            sigma_r[i] = np.std(v_tot[mask])

    return bin_centers, sigma_r, sigma_global


def plot_density_comparison(results_dir, params, output_name='density_profile.png'):
    """
    Create density profile comparison plot with analytical BE curve.
    """
    results_path = Path(results_dir)
    snapshots = sorted(results_path.glob('snapshot_*.csv'))

    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return

    print(f"Found {len(snapshots)} snapshots")

    # Compute analytical profile
    be = compute_be_profile(params)

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: Density vs radius (linear)
    ax1 = axes[0, 0]
    ax1.plot(be['r'], be['rho'], 'k-', lw=2.5, label='Analytical BE', zorder=10)

    colors = plt.cm.viridis(np.linspace(0, 1, min(5, len(snapshots))))
    indices = np.linspace(0, len(snapshots)-1, min(5, len(snapshots)), dtype=int)

    for idx, color in zip(indices, colors):
        snap = snapshots[idx]
        df = load_snapshot(snap)
        r_bins, rho_mean, rho_std, counts = compute_radial_profile(df, n_bins=40)
        mask = counts > 3

        step_num = int(snap.stem.split('_')[-1])
        ax1.scatter(r_bins[mask], rho_mean[mask], s=30, alpha=0.7, color=color,
                   label=f'Step {step_num}')

    ax1.axvline(be['R_cloud'], color='red', ls='--', alpha=0.6, label=f"R = {be['R_cloud']:.2f} pc")
    ax1.axvline(be['r_0'], color='orange', ls=':', alpha=0.6, label=f"r_0 = {be['r_0']:.2f} pc")
    ax1.set_xlabel('Radius [pc]', fontsize=12)
    ax1.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
    ax1.set_title('Density Profile', fontsize=14)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_xlim(0, be['R_cloud'] * 1.3)
    ax1.set_ylim(0, be['rho_c'] * 1.1)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Density vs radius (log scale)
    ax2 = axes[0, 1]
    ax2.semilogy(be['r'], be['rho'], 'k-', lw=2.5, label='Analytical BE', zorder=10)

    for idx, color in zip(indices, colors):
        snap = snapshots[idx]
        df = load_snapshot(snap)
        r_bins, rho_mean, rho_std, counts = compute_radial_profile(df, n_bins=40)
        mask = (counts > 3) & (rho_mean > 0)

        step_num = int(snap.stem.split('_')[-1])
        ax2.scatter(r_bins[mask], rho_mean[mask], s=30, alpha=0.7, color=color,
                   label=f'Step {step_num}')

    ax2.axvline(be['R_cloud'], color='red', ls='--', alpha=0.6)
    ax2.set_xlabel('Radius [pc]', fontsize=12)
    ax2.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=12)
    ax2.set_title('Density Profile (log scale)', fontsize=14)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.set_xlim(0, be['R_cloud'] * 1.3)
    ax2.grid(True, alpha=0.3)

    # Plot 3: Number density
    ax3 = axes[1, 0]
    ax3.plot(be['r'], be['n'], 'k-', lw=2.5, label='Analytical BE', zorder=10)

    # Convert simulation density to number density
    mu = params.get('mu', DEFAULT_PARAMS['mu'])
    density_to_n = (MSUN_CGS / PC_CGS**3) / (mu * M_PROTON_CGS)

    for idx, color in zip(indices, colors):
        snap = snapshots[idx]
        df = load_snapshot(snap)
        r_bins, rho_mean, rho_std, counts = compute_radial_profile(df, n_bins=40)
        mask = counts > 3

        step_num = int(snap.stem.split('_')[-1])
        ax3.scatter(r_bins[mask], rho_mean[mask] * density_to_n, s=30, alpha=0.7,
                   color=color, label=f'Step {step_num}')

    ax3.axhline(params.get('n_center', 980), color='green', ls='-.', alpha=0.6,
                label=f"n_c = {params.get('n_center', 980):.0f} cm$^{{-3}}$")
    ax3.axhline(70, color='purple', ls='-.', alpha=0.6, label=r"n$_{edge}$ = 70 cm$^{-3}$")
    ax3.set_xlabel('Radius [pc]', fontsize=12)
    ax3.set_ylabel(r'Number Density [cm$^{-3}$]', fontsize=12)
    ax3.set_title('Number Density Profile', fontsize=14)
    ax3.legend(fontsize=9, loc='upper right')
    ax3.set_xlim(0, be['R_cloud'] * 1.3)
    ax3.set_ylim(0, be['n_center'] * 1.1)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Velocity dispersion evolution
    ax4 = axes[1, 1]

    steps = []
    sigma_vals = []

    for snap in snapshots:
        df = load_snapshot(snap)
        _, _, sigma_global = compute_velocity_dispersion(df)
        step_num = int(snap.stem.split('_')[-1])
        steps.append(step_num)
        sigma_vals.append(sigma_global)

    ax4.plot(steps, sigma_vals, 'b-o', markersize=4, label='Velocity dispersion')
    ax4.axhline(be['c_s'], color='red', ls='--', alpha=0.6,
                label=f"c_s = {be['c_s']:.3f} km/s (thermal)")
    ax4.set_xlabel('Relaxation Step', fontsize=12)
    ax4.set_ylabel(r'$\sigma_v$ [km/s]', fontsize=12)
    ax4.set_title('Velocity Dispersion Evolution', fontsize=14)
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, max(sigma_vals) * 1.5 if sigma_vals else 1)

    # Title
    plt.suptitle(
        f"Bonnor-Ebert Sphere Analysis (Option A: Atomic CNM)\n"
        f"M = {params.get('M_cloud', 40):.0f} M$_\\odot$, "
        f"T = {params.get('T_cloud', 9):.0f} K, "
        f"$\\mu$ = {params.get('mu', 1.27):.2f}, "
        f"$\\xi_s$ = {params.get('xi_s', 6.45):.2f}",
        fontsize=14
    )

    plt.tight_layout()

    output_path = results_path / output_name
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")

    return fig


def plot_final_comparison(results_dir, params, output_name='final_comparison.png'):
    """
    Create single-panel final comparison with error bars.
    """
    results_path = Path(results_dir)
    snapshots = sorted(results_path.glob('snapshot_*.csv'))

    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return

    # Compute analytical profile
    be = compute_be_profile(params)

    # Load final snapshot
    final_snap = snapshots[-1]
    df = load_snapshot(final_snap)
    r_bins, rho_mean, rho_std, counts = compute_radial_profile(df, n_bins=50)
    mask = counts > 3

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))

    # Analytical
    ax.plot(be['r'], be['rho'], 'k-', lw=2.5, label='Analytical BE solution', zorder=10)
    ax.fill_between(be['r'], be['rho'] * 0.95, be['rho'] * 1.05,
                    color='gray', alpha=0.2, label='$\\pm$5% band')

    # Simulation
    ax.errorbar(r_bins[mask], rho_mean[mask], yerr=rho_std[mask],
               fmt='o', markersize=6, color='blue', capsize=3,
               label='SPH simulation (final)', zorder=5)

    # Markers
    ax.axvline(be['R_cloud'], color='red', ls='--', lw=1.5, alpha=0.7,
               label=f"R$_{{cloud}}$ = {be['R_cloud']:.2f} pc")
    ax.axvline(be['r_0'], color='orange', ls=':', lw=1.5, alpha=0.7,
               label=f"r$_0$ = {be['r_0']:.3f} pc")

    ax.set_xlabel('Radius [pc]', fontsize=14)
    ax.set_ylabel(r'Density [M$_\odot$/pc$^3$]', fontsize=14)
    ax.set_title(
        f"Final Relaxed Bonnor-Ebert Sphere\n"
        f"M = {params.get('M_cloud', 40):.0f} M$_\\odot$, "
        f"T = {params.get('T_cloud', 9):.0f} K, "
        f"$\\mu$ = {params.get('mu', 1.27):.2f}",
        fontsize=14
    )
    ax.legend(fontsize=11, loc='upper right')
    ax.set_xlim(0, be['R_cloud'] * 1.3)
    ax.set_ylim(0, be['rho_c'] * 1.15)
    ax.grid(True, alpha=0.3)

    # Add secondary axis for number density
    mu = params.get('mu', DEFAULT_PARAMS['mu'])
    density_to_n = (MSUN_CGS / PC_CGS**3) / (mu * M_PROTON_CGS)
    ax2 = ax.twinx()
    ax2.set_ylim(0, be['rho_c'] * 1.15 * density_to_n)
    ax2.set_ylabel(r'Number Density [cm$^{-3}$]', fontsize=14)

    plt.tight_layout()

    output_path = results_path / output_name
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")

    # Print statistics
    print("\n=== Final State Statistics ===")
    print(f"Analytical: rho_c = {be['rho_c']:.2f} M_sun/pc^3, n_c = {be['n_center']:.0f} cm^-3")

    # Find central bin
    central_mask = r_bins < 0.1
    if central_mask.any() and counts[central_mask].sum() > 0:
        rho_sim_center = rho_mean[central_mask & (counts > 0)].mean()
        n_sim_center = rho_sim_center * density_to_n
        print(f"Simulation: rho_c ~ {rho_sim_center:.2f} M_sun/pc^3, n_c ~ {n_sim_center:.0f} cm^-3")
        print(f"Relative error: {abs(rho_sim_center - be['rho_c'])/be['rho_c']*100:.1f}%")

    return fig


def main():
    parser = argparse.ArgumentParser(description='Bonnor-Ebert sphere analysis')
    parser.add_argument('--results', '-r', type=str,
                        default='results/phase1_relaxation',
                        help='Results directory')
    parser.add_argument('--config', '-c', type=str,
                        help='Config JSON file (optional)')
    args = parser.parse_args()

    script_dir = Path(__file__).parent.parent
    results_dir = script_dir / args.results

    # Load parameters from config or use defaults
    params = DEFAULT_PARAMS.copy()
    if args.config:
        config_path = Path(args.config)
        if config_path.exists():
            with open(config_path) as f:
                config = json.load(f)
                for key in DEFAULT_PARAMS:
                    if key in config:
                        params[key] = config[key]

    print("=== Bonnor-Ebert Sphere Analysis ===")
    print(f"Results directory: {results_dir}")
    print(f"Parameters: mu={params['mu']}, T={params['T_cloud']}K, "
          f"n_c={params['n_center']} cm^-3, xi_s={params['xi_s']}")

    # Generate plots
    plot_density_comparison(results_dir, params)
    plot_final_comparison(results_dir, params)


if __name__ == "__main__":
    main()
