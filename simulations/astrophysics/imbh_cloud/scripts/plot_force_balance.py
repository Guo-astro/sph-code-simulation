#!/usr/bin/python3
"""Plot radial profiles of density, gravitational acceleration, and hydrostatic acceleration.

This script analyzes SPH simulation snapshots to verify force balance in
self-gravitating gas clouds (e.g., Bonnor-Ebert spheres).

Forces computed:
- Gravity: g_r = -dφ/dr (from potential gradient)
- Hydro: a_hydro = -(1/ρ) dP/dr (from pressure gradient)
- For equilibrium: g_r + a_hydro ≈ 0

Usage:
    python plot_force_balance.py <snapshot.h5> [--output plot.png]
    python plot_force_balance.py <result_dir> --gif  # Create animation
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

def load_hdf5(filename):
    """Load HDF5 snapshot file."""
    if not HAS_H5PY:
        raise ImportError("h5py required for HDF5 files")

    data = {}
    with h5py.File(filename, 'r') as f:
        field_map = {
            'pos_x': 'particles/pos_x',
            'pos_y': 'particles/pos_y',
            'pos_z': 'particles/pos_z',
            'acc_x': 'particles/acc_x',
            'acc_y': 'particles/acc_y',
            'acc_z': 'particles/acc_z',
            'dens': 'particles/dens',
            'mass': 'particles/mass',
            'pres': 'particles/pres',
            'phi': 'particles/phi',
            'sml': 'particles/sml',
        }

        for key, path in field_map.items():
            if path in f:
                data[key] = np.array(f[path])

        if 'metadata' in f and 'time' in f['metadata']:
            data['_time'] = float(f['metadata/time'][()])
        else:
            data['_time'] = 0.0

    return data


def load_config(result_dir):
    """Try to load config from result directory or parent."""
    search_paths = [
        Path(result_dir) / "config.json",
        Path(result_dir).parent / "config.json",
    ]

    for config_file in Path(result_dir).parents:
        search_paths.append(config_file / "config.json")
        if "results" in str(config_file):
            break

    for p in search_paths:
        if p.exists():
            with open(p) as f:
                return json.load(f)

    # Search for JSON file with same base name as directory
    base_name = Path(result_dir).name
    for ext_dir in [Path(result_dir).parent.parent / "config" / "presets"]:
        if ext_dir.exists():
            for json_file in ext_dir.rglob(f"*{base_name}*.json"):
                with open(json_file) as f:
                    return json.load(f)

    return None


def solve_lane_emden(xi_s, n_points=1000):
    """Solve Lane-Emden equation for isothermal case using RK4."""
    def f(xi, psi, dpsi):
        """Return d²ψ/dξ² = exp(-ψ) - (2/ξ) dψ/dξ"""
        if xi < 1e-10:
            return 0.0
        return np.exp(-psi) - (2.0 / xi) * dpsi

    # RK4 integration
    xi = np.linspace(1e-10, xi_s, n_points)
    h = xi[1] - xi[0]
    psi = np.zeros(n_points)
    dpsi = np.zeros(n_points)

    # Initial conditions at xi ~ 0
    psi[0] = 0.0
    dpsi[0] = 0.0

    for i in range(n_points - 1):
        xi_i = xi[i]
        psi_i = psi[i]
        dpsi_i = dpsi[i]

        # RK4 for coupled ODEs: dψ/dξ = dpsi, d(dpsi)/dξ = f
        k1_psi = dpsi_i
        k1_dpsi = f(xi_i, psi_i, dpsi_i)

        k2_psi = dpsi_i + 0.5 * h * k1_dpsi
        k2_dpsi = f(xi_i + 0.5 * h, psi_i + 0.5 * h * k1_psi, dpsi_i + 0.5 * h * k1_dpsi)

        k3_psi = dpsi_i + 0.5 * h * k2_dpsi
        k3_dpsi = f(xi_i + 0.5 * h, psi_i + 0.5 * h * k2_psi, dpsi_i + 0.5 * h * k2_dpsi)

        k4_psi = dpsi_i + h * k3_dpsi
        k4_dpsi = f(xi_i + h, psi_i + h * k3_psi, dpsi_i + h * k3_dpsi)

        psi[i + 1] = psi_i + (h / 6.0) * (k1_psi + 2 * k2_psi + 2 * k3_psi + k4_psi)
        dpsi[i + 1] = dpsi_i + (h / 6.0) * (k1_dpsi + 2 * k2_dpsi + 2 * k3_dpsi + k4_dpsi)

    return xi, psi, dpsi


def compute_analytic_forces(config, r_points):
    """Compute analytical force profiles for Bonnor-Ebert sphere.

    Returns gravitational and hydrostatic accelerations (radial components).
    Positive = outward, Negative = inward.
    """
    xi_s = config.get('xi_s', 6.5)
    n_center = config.get('n_center', 500.0)
    T_cloud = config.get('T_cloud', 35.0)
    mu = config.get('mu', 1.27)
    G = config.get('G', 0.00430091)

    # Physical constants in code units
    k_B = 0.831446  # (km/s)^2 * amu / K
    m_H = 0.8378    # Msun (hydrogen mass in code mass unit)

    # Sound speed
    c_s = np.sqrt(k_B * T_cloud / (mu * m_H))  # km/s

    # Central density in code units (Msun/pc^3)
    m_H_cgs = 1.6726e-24  # g
    Msun = 1.989e33  # g
    pc = 3.086e18  # cm
    rho_center = n_center * mu * m_H_cgs * (pc**3) / Msun

    # Scale radius
    r0 = c_s / np.sqrt(4 * np.pi * G * rho_center)

    # Solve Lane-Emden
    xi, psi, dpsi = solve_lane_emden(xi_s, n_points=2000)
    R_cloud = xi_s * r0

    # Interpolate to r_points
    xi_interp = r_points / r0
    xi_interp = np.clip(xi_interp, xi[0], xi[-1])

    psi_interp = np.interp(xi_interp, xi, psi)
    dpsi_interp = np.interp(xi_interp, xi, dpsi)

    # Density profile: rho(r) = rho_c * exp(-psi)
    rho_analytic = rho_center * np.exp(-psi_interp)

    # Gravitational acceleration (radial, inward is negative)
    # From Lane-Emden: g_r = -c_s^2 * (1/r) * dpsi/d(ln xi) = -c_s^2 * dpsi/dr
    # But dpsi/dr = dpsi/dxi * dxi/dr = dpsi/dxi / r0
    # So g_r = -c_s^2 * dpsi/dxi / r0
    g_grav = np.zeros_like(r_points)
    mask = xi_interp > 1e-6
    g_grav[mask] = -c_s**2 * dpsi_interp[mask] / r0  # negative = inward

    # Hydrostatic acceleration: a_hydro = -(1/rho) dP/dr = -c_s^2 * (d ln rho/dr)
    # For isothermal: P = rho * c_s^2, so dP/dr = c_s^2 * drho/dr
    # a_hydro = -c_s^2 * (drho/dr) / rho = -c_s^2 * d(ln rho)/dr
    # Since rho = rho_c * exp(-psi), d(ln rho)/dr = -dpsi/dr = -dpsi/dxi / r0
    # So a_hydro = c_s^2 * dpsi/dxi / r0
    a_hydro = np.zeros_like(r_points)
    a_hydro[mask] = c_s**2 * dpsi_interp[mask] / r0  # positive = outward

    return {
        'r': r_points,
        'rho': rho_analytic,
        'g_grav': g_grav,  # gravity (negative = inward)
        'a_hydro': a_hydro,  # hydro (positive = outward)
        'R_cloud': R_cloud,
        'c_s': c_s,
        'rho_center': rho_center,
    }


def compute_radial_profiles(data, n_bins=100):
    """Compute binned radial profiles from particle data."""
    # Compute radius
    x = data.get('pos_x', np.zeros(1))
    y = data.get('pos_y', np.zeros(1))
    z = data.get('pos_z', np.zeros(1))
    r = np.sqrt(x**2 + y**2 + z**2)

    # Radial component of total acceleration
    acc_x = data.get('acc_x', np.zeros_like(x))
    acc_y = data.get('acc_y', np.zeros_like(y))
    acc_z = data.get('acc_z', np.zeros_like(z))
    acc_mag = np.sqrt(acc_x**2 + acc_y**2 + acc_z**2)

    # Radial acceleration (positive = outward)
    r_safe = np.maximum(r, 1e-10)
    acc_r = (acc_x * x + acc_y * y + acc_z * z) / r_safe

    # Bin by radius
    r_max = np.percentile(r, 99)
    r_bins = np.linspace(0, r_max, n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

    # Compute binned profiles
    profiles = {
        'r': r_centers,
        'dens': np.zeros(n_bins),
        'pres': np.zeros(n_bins),
        'phi': np.zeros(n_bins),
        'acc_r': np.zeros(n_bins),
        'acc_mag': np.zeros(n_bins),
        'count': np.zeros(n_bins, dtype=int),
    }

    indices = np.digitize(r, r_bins) - 1
    indices = np.clip(indices, 0, n_bins - 1)

    for i in range(n_bins):
        mask = indices == i
        if np.sum(mask) > 0:
            profiles['count'][i] = np.sum(mask)
            profiles['dens'][i] = np.mean(data.get('dens', np.zeros_like(x))[mask])
            profiles['pres'][i] = np.mean(data.get('pres', np.zeros_like(x))[mask])
            profiles['phi'][i] = np.mean(data.get('phi', np.zeros_like(x))[mask])
            profiles['acc_r'][i] = np.mean(acc_r[mask])
            profiles['acc_mag'][i] = np.mean(acc_mag[mask])

    # Compute numerical gradients for forces
    dr = r_centers[1] - r_centers[0] if len(r_centers) > 1 else 1.0

    # Gravity from potential: g = -dφ/dr
    profiles['g_grav'] = -np.gradient(profiles['phi'], dr)

    # Hydro from pressure: a_hydro = -(1/ρ) dP/dr
    rho_safe = np.maximum(profiles['dens'], 1e-10)
    profiles['a_hydro'] = -np.gradient(profiles['pres'], dr) / rho_safe

    return profiles


def plot_force_balance(filename, output_file=None, config=None, show=True):
    """Create force balance plot."""
    data = load_hdf5(filename)
    time = data.get('_time', 0.0)

    profiles = compute_radial_profiles(data)

    # Try to load config
    if config is None:
        result_dir = str(Path(filename).parent)
        config = load_config(result_dir)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: Density profile
    ax1 = axes[0, 0]
    ax1.plot(profiles['r'], profiles['dens'], 'b-', linewidth=2, label='SPH')

    if config:
        analytic = compute_analytic_forces(config, profiles['r'])
        ax1.plot(profiles['r'], analytic['rho'], 'r--', linewidth=2, label='Analytic BE')
        ax1.axvline(analytic['R_cloud'], color='gray', linestyle=':', label=f"R_cloud={analytic['R_cloud']:.2f}")

    ax1.set_xlabel('Radius [pc]')
    ax1.set_ylabel('Density [M☉/pc³]')
    ax1.set_title(f'Density Profile (t={time:.2f})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Panel 2: Gravitational potential
    ax2 = axes[0, 1]
    ax2.plot(profiles['r'], profiles['phi'], 'g-', linewidth=2, label='φ (potential)')
    ax2.set_xlabel('Radius [pc]')
    ax2.set_ylabel('Gravitational Potential φ')
    ax2.set_title('Gravitational Potential')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Panel 3: Accelerations
    ax3 = axes[1, 0]
    ax3.plot(profiles['r'], profiles['g_grav'], 'b-', linewidth=2, label='Gravity (g = -dφ/dr)')
    ax3.plot(profiles['r'], profiles['a_hydro'], 'r-', linewidth=2, label='Hydro (a = -(1/ρ)dP/dr)')
    ax3.plot(profiles['r'], profiles['acc_r'], 'k--', linewidth=2, label='Total (radial)')

    if config:
        ax3.plot(profiles['r'], analytic['g_grav'], 'b:', linewidth=1.5, label='Gravity (analytic)')
        ax3.plot(profiles['r'], analytic['a_hydro'], 'r:', linewidth=1.5, label='Hydro (analytic)')

    ax3.axhline(0, color='gray', linestyle='-', alpha=0.5)
    ax3.set_xlabel('Radius [pc]')
    ax3.set_ylabel('Acceleration [(km/s)²/pc]')
    ax3.set_title('Force Balance: Gravity vs Hydrostatic')
    ax3.legend(loc='best', fontsize=8)
    ax3.grid(True, alpha=0.3)

    # Panel 4: Net acceleration (deviation from equilibrium)
    ax4 = axes[1, 1]
    net_acc = profiles['g_grav'] + profiles['a_hydro']
    ax4.plot(profiles['r'], net_acc, 'purple', linewidth=2, label='Net (g + a_hydro)')
    ax4.plot(profiles['r'], profiles['acc_r'], 'k--', linewidth=2, label='Total (from sim)')
    ax4.axhline(0, color='gray', linestyle='-', alpha=0.5)
    ax4.set_xlabel('Radius [pc]')
    ax4.set_ylabel('Net Acceleration')
    ax4.set_title('Deviation from Hydrostatic Equilibrium')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150)
        print(f"Saved: {output_file}")

    if show:
        plt.show()
    else:
        plt.close()

    return profiles


def create_animation(result_dir, output_file='force_balance.gif', config=None):
    """Create GIF animation of force evolution."""
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except ImportError:
        print("matplotlib animation support required for GIF creation")
        return

    # Find snapshots
    result_path = Path(result_dir)
    h5_files = sorted(result_path.glob("snapshot_*.h5"))

    if not h5_files:
        print(f"No HDF5 snapshots found in {result_dir}")
        return

    print(f"Found {len(h5_files)} snapshots")

    # Load config
    if config is None:
        config = load_config(result_dir)

    # Setup figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    def update(frame_idx):
        filename = h5_files[frame_idx]

        for ax in axes.flat:
            ax.clear()

        try:
            data = load_hdf5(str(filename))
            time = data.get('_time', 0.0)
            profiles = compute_radial_profiles(data)

            # Panel 1: Density
            ax1 = axes[0, 0]
            ax1.plot(profiles['r'], profiles['dens'], 'b-', linewidth=2, label='SPH')
            if config:
                analytic = compute_analytic_forces(config, profiles['r'])
                ax1.plot(profiles['r'], analytic['rho'], 'r--', linewidth=2, label='Analytic')
            ax1.set_xlabel('Radius [pc]')
            ax1.set_ylabel('Density [M☉/pc³]')
            ax1.set_title(f'Density Profile (t={time:.2f})')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            ax1.set_ylim(0, 20)

            # Panel 2: Potential
            ax2 = axes[0, 1]
            ax2.plot(profiles['r'], profiles['phi'], 'g-', linewidth=2)
            ax2.set_xlabel('Radius [pc]')
            ax2.set_ylabel('φ')
            ax2.set_title('Gravitational Potential')
            ax2.grid(True, alpha=0.3)

            # Panel 3: Accelerations
            ax3 = axes[1, 0]
            ax3.plot(profiles['r'], profiles['g_grav'], 'b-', linewidth=2, label='Gravity')
            ax3.plot(profiles['r'], profiles['a_hydro'], 'r-', linewidth=2, label='Hydro')
            ax3.axhline(0, color='gray', linestyle='-', alpha=0.5)
            ax3.set_xlabel('Radius [pc]')
            ax3.set_ylabel('Acceleration')
            ax3.set_title('Force Components')
            ax3.legend()
            ax3.grid(True, alpha=0.3)

            # Panel 4: Net
            ax4 = axes[1, 1]
            net_acc = profiles['g_grav'] + profiles['a_hydro']
            ax4.plot(profiles['r'], net_acc, 'purple', linewidth=2, label='Net')
            ax4.plot(profiles['r'], profiles['acc_r'], 'k--', linewidth=2, label='Total (sim)')
            ax4.axhline(0, color='gray', linestyle='-', alpha=0.5)
            ax4.set_xlabel('Radius [pc]')
            ax4.set_ylabel('Net Acceleration')
            ax4.set_title('Equilibrium Deviation')
            ax4.legend()
            ax4.grid(True, alpha=0.3)

        except Exception as e:
            print(f"Error processing {filename}: {e}")

        plt.tight_layout()
        print(f"\rFrame {frame_idx + 1}/{len(h5_files)}", end='', flush=True)

    anim = FuncAnimation(fig, update, frames=len(h5_files), interval=200)

    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer)
    print(f"\nSaved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Plot force balance from SPH snapshots')
    parser.add_argument('input', help='Snapshot file or result directory')
    parser.add_argument('--output', '-o', help='Output image file')
    parser.add_argument('--gif', action='store_true', help='Create GIF animation')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')

    args = parser.parse_args()

    if args.gif:
        output = args.output or 'force_balance.gif'
        create_animation(args.input, output)
    else:
        output = args.output or 'force_balance.png'
        plot_force_balance(args.input, output, show=not args.no_show)


if __name__ == '__main__':
    main()
