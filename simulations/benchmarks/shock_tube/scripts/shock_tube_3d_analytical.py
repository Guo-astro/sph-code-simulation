#!/usr/bin/env python3
"""
3D Shock Tube Analytical Solution Comparison

Compares SPH simulation results with the exact Riemann solution for a shock tube.
The analytical solution is computed in 1D (x-direction) and averaged over y and z
directions for comparison with 3D simulation data.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path


def exact_sod_solution(x, t, gamma=1.4, x0=0.5):
    """
    Exact solution for the Sod shock tube problem.

    Initial conditions:
        Left state (x < x0):  rho=1.0, p=1.0, u=0.0
        Right state (x > x0): rho=0.125, p=0.1, u=0.0
    """
    rho_L, p_L, u_L = 1.0, 1.0, 0.0
    rho_R, p_R, u_R = 0.125, 0.1, 0.0

    c_L = np.sqrt(gamma * p_L / rho_L)
    c_R = np.sqrt(gamma * p_R / rho_R)

    # Exact star region values for standard Sod problem (gamma=1.4)
    p_star = 0.30313
    u_star = 0.92745
    rho_star_L = 0.42632
    rho_star_R = 0.26557

    c_star = np.sqrt(gamma * p_star / rho_star_L)

    # Wave positions at time t
    x_head = x0 - c_L * t
    x_tail = x0 + (u_star - c_star) * t
    x_contact = x0 + u_star * t
    shock_speed = rho_star_R * u_star / (rho_star_R - rho_R)
    x_shock = x0 + shock_speed * t

    n = len(x)
    rho = np.zeros(n)
    u = np.zeros(n)
    p = np.zeros(n)

    for i in range(n):
        xi = x[i]

        if xi < x_head:
            rho[i], u[i], p[i] = rho_L, u_L, p_L
        elif xi < x_tail:
            u[i] = 2.0 / (gamma + 1) * (c_L + (xi - x0) / t)
            c = c_L - (gamma - 1) / 2 * u[i]
            rho[i] = rho_L * (c / c_L) ** (2.0 / (gamma - 1))
            p[i] = p_L * (c / c_L) ** (2.0 * gamma / (gamma - 1))
        elif xi < x_contact:
            rho[i], u[i], p[i] = rho_star_L, u_star, p_star
        elif xi < x_shock:
            rho[i], u[i], p[i] = rho_star_R, u_star, p_star
        else:
            rho[i], u[i], p[i] = rho_R, u_R, p_R

    e = p / ((gamma - 1) * rho)

    wave_info = {
        'x_head': x_head,
        'x_tail': x_tail,
        'x_contact': x_contact,
        'x_shock': x_shock,
        'shock_speed': shock_speed,
        'p_star': p_star,
        'u_star': u_star,
    }

    return rho, u, p, e, wave_info


def load_snapshot(filename):
    """Load SPH snapshot CSV file."""
    data = {}
    metadata = {}

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break

        f.seek(0)
        lines = [l for l in f.readlines() if not l.startswith('#')]

        if len(lines) < 2:
            raise ValueError("No data in file")

        header = lines[0].strip().split(',')

        for col in header:
            data[col] = []

        for line in lines[1:]:
            values = line.strip().split(',')
            for i, col in enumerate(header):
                try:
                    data[col].append(float(values[i]))
                except (ValueError, IndexError):
                    pass

        for key in data:
            data[key] = np.array(data[key])

    data['metadata'] = metadata
    return data


def plot_shock_tube_3d_comparison(snapshot_file, output_file=None, show_plot=True):
    """
    Create comparison plots between SPH and analytical solution.
    For 3D data, particles are averaged in y and z directions for each x-position bin.
    """
    print(f"Loading {snapshot_file}...")
    data = load_snapshot(snapshot_file)

    try:
        time_str = data['metadata'].get('Time (code)', '0.1')
        time = float(time_str)
    except:
        print("Warning: Could not extract time from metadata, using t=0.1")
        time = 0.1

    try:
        gamma_str = data['metadata'].get('Gamma', '1.4')
        gamma = float(gamma_str)
    except:
        gamma = 1.4

    print(f"Time: {time:.4f}, Gamma: {gamma}")

    x = data['pos_x']
    rho = data['dens']
    vx = data['vel_x']
    pres = data['pres']
    ene = data['ene']

    print(f"Total particles: {len(x)}")

    # Bin data in x-direction (average over y and z)
    x_min, x_max = 0.0, 1.0
    n_bins = 100
    x_bins = np.linspace(x_min, x_max, n_bins + 1)
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    rho_avg = np.zeros(n_bins)
    vx_avg = np.zeros(n_bins)
    pres_avg = np.zeros(n_bins)
    ene_avg = np.zeros(n_bins)
    counts = np.zeros(n_bins)

    for i in range(len(x)):
        bin_idx = int((x[i] - x_min) / (x_max - x_min) * n_bins)
        if 0 <= bin_idx < n_bins:
            rho_avg[bin_idx] += rho[i]
            vx_avg[bin_idx] += vx[i]
            pres_avg[bin_idx] += pres[i]
            ene_avg[bin_idx] += ene[i]
            counts[bin_idx] += 1

    mask = counts > 0
    rho_avg[mask] /= counts[mask]
    vx_avg[mask] /= counts[mask]
    pres_avg[mask] /= counts[mask]
    ene_avg[mask] /= counts[mask]

    # Analytical solution
    x_exact = np.linspace(x_min, x_max, 1000)
    rho_exact, u_exact, p_exact, e_exact, wave_info = exact_sod_solution(
        x_exact, time, gamma=gamma, x0=0.5
    )

    print(f"Exact shock position: x = {wave_info['x_shock']:.4f}")
    print(f"Exact shock speed: {wave_info['shock_speed']:.4f}")
    print(f"Exact contact position: x = {wave_info['x_contact']:.4f}")

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'3D Shock Tube: GSPH HLL (t = {time:.4f})\n'
                 f'Exact shock at x={wave_info["x_shock"]:.3f}, '
                 f'contact at x={wave_info["x_contact"]:.3f}',
                 fontsize=14, fontweight='bold')

    # Density
    ax = axes[0, 0]
    ax.plot(x_exact, rho_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], rho_avg[mask], c='red', s=20, alpha=0.6, label='SPH (yz-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5, label='Exact shock')
    ax.axvline(wave_info['x_contact'], color='orange', ls=':', alpha=0.5, label='Contact')
    ax.set_xlabel('x')
    ax.set_ylabel('Density')
    ax.set_title('Density Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Velocity
    ax = axes[0, 1]
    ax.plot(x_exact, u_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], vx_avg[mask], c='blue', s=20, alpha=0.6, label='SPH (yz-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5)
    ax.axvline(wave_info['x_contact'], color='orange', ls=':', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('Velocity')
    ax.set_title('Velocity Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Pressure
    ax = axes[1, 0]
    ax.plot(x_exact, p_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], pres_avg[mask], c='green', s=20, alpha=0.6, label='SPH (yz-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5)
    ax.axvline(wave_info['x_contact'], color='orange', ls=':', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('Pressure')
    ax.set_title('Pressure Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Internal Energy
    ax = axes[1, 1]
    ax.plot(x_exact, e_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], ene_avg[mask], c='purple', s=20, alpha=0.6, label='SPH (yz-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5)
    ax.axvline(wave_info['x_contact'], color='orange', ls=':', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('Specific Internal Energy')
    ax.set_title('Internal Energy Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")

    if show_plot:
        plt.show()
    else:
        plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Plot 3D shock tube with analytical solution')
    parser.add_argument('snapshot', help='Path to snapshot CSV file')
    parser.add_argument('-o', '--output', help='Output image file')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')

    args = parser.parse_args()

    plot_shock_tube_3d_comparison(
        args.snapshot,
        output_file=args.output,
        show_plot=not args.no_show
    )


if __name__ == '__main__':
    main()
