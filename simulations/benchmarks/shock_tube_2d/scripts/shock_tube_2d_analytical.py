#!/usr/bin/env python3
"""
2D Shock Tube Analytical Solution Comparison

Compares SPH simulation results with the exact Riemann solution for a shock tube.
The analytical solution is computed in 1D (x-direction) and averaged over the y-direction
for comparison with 2D simulation data.

Usage:
    python3 shock_tube_2d_analytical.py <snapshot.csv> -o <output.png>
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

    Uses known exact values for the standard Sod problem with gamma=1.4.

    Returns: rho, u, p, e (density, velocity, pressure, specific internal energy)
    """
    # Initial states
    rho_L, p_L, u_L = 1.0, 1.0, 0.0
    rho_R, p_R, u_R = 0.125, 0.1, 0.0

    # Sound speeds
    c_L = np.sqrt(gamma * p_L / rho_L)  # = 1.1832
    c_R = np.sqrt(gamma * p_R / rho_R)  # = 1.0583

    # Exact star region values for standard Sod problem (gamma=1.4)
    # These are the exact values from solving the Riemann problem
    p_star = 0.30313
    u_star = 0.92745
    rho_star_L = 0.42632  # Density left of contact (in star region)
    rho_star_R = 0.26557  # Density right of contact (in star region)

    # Sound speed in star region (left side)
    c_star = np.sqrt(gamma * p_star / rho_star_L)

    # Wave positions at time t
    # 1. Rarefaction head: moves at -c_L into left state
    x_head = x0 - c_L * t

    # 2. Rarefaction tail: moves at u_star - c_star
    x_tail = x0 + (u_star - c_star) * t

    # 3. Contact discontinuity: moves at u_star
    x_contact = x0 + u_star * t

    # 4. Shock wave speed from Rankine-Hugoniot relations
    #    Mass conservation: rho_R * S = rho_star_R * (S - u_star)
    #    Solving: S = rho_star_R * u_star / (rho_star_R - rho_R)
    #    Or equivalently: S = c_R * sqrt((gamma+1)/(2*gamma) * (p_star/p_R) + (gamma-1)/(2*gamma))
    #    NOTE: Do NOT add u_star - that was a bug!
    shock_speed = rho_star_R * u_star / (rho_star_R - rho_R)  # ≈ 1.75
    x_shock = x0 + shock_speed * t

    # Initialize solution arrays
    n = len(x)
    rho = np.zeros(n)
    u = np.zeros(n)
    p = np.zeros(n)

    for i in range(n):
        xi = x[i]

        if xi < x_head:
            # Region 1: Undisturbed left state
            rho[i], u[i], p[i] = rho_L, u_L, p_L
        elif xi < x_tail:
            # Region 2: Rarefaction fan (self-similar solution)
            # u = 2/(gamma+1) * (c_L + (gamma-1)/2 * u_L + (x-x0)/t)
            # c = c_L + (gamma-1)/2 * (u_L - u)
            u[i] = 2.0 / (gamma + 1) * (c_L + (xi - x0) / t)
            c = c_L - (gamma - 1) / 2 * u[i]
            rho[i] = rho_L * (c / c_L) ** (2.0 / (gamma - 1))
            p[i] = p_L * (c / c_L) ** (2.0 * gamma / (gamma - 1))
        elif xi < x_contact:
            # Region 3: Star region (left of contact)
            rho[i], u[i], p[i] = rho_star_L, u_star, p_star
        elif xi < x_shock:
            # Region 4: Star region (right of contact)
            rho[i], u[i], p[i] = rho_star_R, u_star, p_star
        else:
            # Region 5: Undisturbed right state
            rho[i], u[i], p[i] = rho_R, u_R, p_R

    # Specific internal energy: e = p / ((gamma - 1) * rho)
    e = p / ((gamma - 1) * rho)

    # Return wave positions for annotation
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
        # Read metadata
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break

        # Read data
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

        # Convert to numpy arrays
        for key in data:
            data[key] = np.array(data[key])

    data['metadata'] = metadata
    return data


def plot_shock_tube_2d_comparison(snapshot_file, output_file=None, show_plot=True):
    """
    Create comparison plots between SPH and analytical solution.

    For 2D data, particles are averaged in y-direction for each x-position bin.
    """
    # Load snapshot
    print(f"Loading {snapshot_file}...")
    data = load_snapshot(snapshot_file)

    # Extract simulation time from metadata
    try:
        time_str = data['metadata'].get('Time (code)', '0.2')
        time = float(time_str)
    except:
        print("Warning: Could not extract time from metadata, using t=0.2")
        time = 0.2

    # Extract gamma
    try:
        gamma_str = data['metadata'].get('Gamma', '1.4')
        gamma = float(gamma_str)
    except:
        gamma = 1.4

    print(f"Time: {time:.4f}, Gamma: {gamma}")

    # Get particle data
    x = data['pos_x']
    y = data['pos_y']
    rho = data['dens']
    vx = data['vel_x']
    pres = data['pres']
    ene = data['ene']

    # Bin data in x-direction (average over y)
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

    # Average
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

    # Estimate SPH shock position
    sph_shock_x = None
    for i in range(len(x_centers)):
        if x_centers[i] > 0.55 and mask[i] and rho_avg[i] < 0.2:
            sph_shock_x = x_centers[i]
            break

    if sph_shock_x:
        sph_shock_speed = (sph_shock_x - 0.5) / time
        print(f"SPH shock position: x = {sph_shock_x:.4f}")
        print(f"SPH shock speed: {sph_shock_speed:.4f}")
        print(f"SPH/Exact shock speed ratio: {sph_shock_speed/wave_info['shock_speed']*100:.1f}%")

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'2D Shock Tube Comparison (t = {time:.4f})\n'
                 f'Exact shock at x={wave_info["x_shock"]:.3f}, '
                 f'contact at x={wave_info["x_contact"]:.3f}',
                 fontsize=14, fontweight='bold')

    # Density
    ax = axes[0, 0]
    ax.plot(x_exact, rho_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], rho_avg[mask], c='red', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5, label=f'Exact shock')
    if sph_shock_x:
        ax.axvline(sph_shock_x, color='blue', ls=':', alpha=0.5, label=f'SPH shock')
    ax.set_xlabel('x')
    ax.set_ylabel('Density')
    ax.set_title('Density Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Velocity
    ax = axes[0, 1]
    ax.plot(x_exact, u_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], vx_avg[mask], c='blue', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('Velocity')
    ax.set_title('Velocity Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Pressure
    ax = axes[1, 0]
    ax.plot(x_exact, p_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], pres_avg[mask], c='green', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('Pressure')
    ax.set_title('Pressure Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Internal Energy
    ax = axes[1, 1]
    ax.plot(x_exact, e_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], ene_avg[mask], c='purple', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.axvline(wave_info['x_shock'], color='green', ls='--', alpha=0.5)
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

    parser = argparse.ArgumentParser(description='Plot 2D shock tube with analytical solution')
    parser.add_argument('snapshot', help='Path to snapshot CSV file')
    parser.add_argument('-o', '--output', help='Output image file')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')

    args = parser.parse_args()

    plot_shock_tube_2d_comparison(
        args.snapshot,
        output_file=args.output,
        show_plot=not args.no_show
    )


if __name__ == '__main__':
    main()
