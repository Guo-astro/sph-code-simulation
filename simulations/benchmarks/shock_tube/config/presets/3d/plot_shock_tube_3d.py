#!/usr/bin/env python3
"""Plot 3D shock tube results with analytical Sod solution overlay.

Shows RAW particle data (no binning) with analytic solution.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys


def read_csv_sph(filename):
    """Read SPH CSV file, skipping comment header."""
    skip_rows = 0
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('id,'):
                skip_rows = i
                break
    data = np.genfromtxt(filename, delimiter=',', names=True, skip_header=skip_rows)
    return data


def sod_exact_solution(x, t, gamma=1.4, x0=0.5):
    """Exact solution for the Sod shock tube problem.

    Initial conditions:
    - Left (x < 0.5):  rho=1.0, P=1.0, u=0
    - Right (x > 0.5): rho=0.125, P=0.1, u=0
    """
    rho_L, P_L, u_L = 1.0, 1.0, 0.0
    rho_R, P_R, u_R = 0.125, 0.1, 0.0
    c_L = np.sqrt(gamma * P_L / rho_L)

    # Star region values (pre-computed for standard Sod)
    P_star = 0.30313
    u_star = 0.92745
    rho_star_L = rho_L * (P_star / P_L) ** (1.0 / gamma)
    rho_star_R = rho_R * ((P_star / P_R + (gamma - 1) / (gamma + 1)) /
                          ((gamma - 1) / (gamma + 1) * P_star / P_R + 1))
    c_star_L = np.sqrt(gamma * P_star / rho_star_L)

    # Wave positions
    x_head = x0 - c_L * t
    x_tail = x0 + (u_star - c_star_L) * t
    x_contact = x0 + u_star * t
    W_s = (rho_star_R * u_star) / (rho_star_R - rho_R)
    x_shock = x0 + W_s * t

    rho, u, P = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi < x_head:
            rho[i], u[i], P[i] = rho_L, u_L, P_L
        elif xi < x_tail:
            u[i] = (2 / (gamma + 1)) * ((xi - x0) / t + c_L)
            c = c_L - (gamma - 1) / 2 * u[i]
            rho[i] = rho_L * (c / c_L) ** (2 / (gamma - 1))
            P[i] = P_L * (c / c_L) ** (2 * gamma / (gamma - 1))
        elif xi < x_contact:
            rho[i], u[i], P[i] = rho_star_L, u_star, P_star
        elif xi < x_shock:
            rho[i], u[i], P[i] = rho_star_R, u_star, P_star
        else:
            rho[i], u[i], P[i] = rho_R, u_R, P_R
    return rho, u, P


def main():
    # Configuration - updated for Morton/iterative traversal test
    data_dir = Path("simulations/benchmarks/shock_tube/results/gsph_3d_wendland_n120")

    # Find latest snapshot
    csv_files = sorted(data_dir.glob("snapshot_*.csv"))
    if not csv_files:
        print(f"No CSV files found in {data_dir}")
        print("Run the simulation first!")
        sys.exit(1)

    csv_file = csv_files[-1]  # Use latest

    # Extract time from filename (snapshot_NNNN.csv -> t = NNNN * output_time)
    snap_num = int(csv_file.stem.split('_')[1])
    t = snap_num * 0.01  # outputTime from config

    print(f"Reading: {csv_file}")
    print(f"Simulation time: t = {t:.3f}")

    data = read_csv_sph(csv_file)
    x, vx, rho, pres = data['pos_x'], data['vel_x'], data['dens'], data['pres']

    print(f"Total particles: {len(x)}")
    print(f"x range: [{x.min():.3f}, {x.max():.3f}]")
    print(f"rho range: [{rho.min():.4f}, {rho.max():.4f}]")

    # Analytical solution at simulation time
    x_exact = np.linspace(0, 1, 1000)
    rho_exact, u_exact, P_exact = sod_exact_solution(x_exact, t)
    gamma = 1.4

    # Plot - RAW particle data (no binning)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Density - RAW particles
    ax = axes[0, 0]
    ax.plot(x_exact, rho_exact, 'r-', lw=2, label='Exact solution', zorder=10)
    ax.scatter(x, rho, s=1, c='blue', alpha=0.5, label=f'Particles (N={len(x)})')
    ax.set_xlabel('x')
    ax.set_ylabel('Density')
    ax.set_title('Density (raw particle data)')
    ax.legend()
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)

    # Velocity - RAW particles
    ax = axes[0, 1]
    ax.plot(x_exact, u_exact, 'r-', lw=2, label='Exact solution', zorder=10)
    ax.scatter(x, vx, s=1, c='blue', alpha=0.5, label=f'Particles (N={len(x)})')
    ax.set_xlabel('x')
    ax.set_ylabel('Velocity')
    ax.set_title('Velocity (raw particle data)')
    ax.legend()
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)

    # Pressure - RAW particles
    ax = axes[1, 0]
    ax.plot(x_exact, P_exact, 'r-', lw=2, label='Exact solution', zorder=10)
    ax.scatter(x, pres, s=1, c='blue', alpha=0.5, label=f'Particles (N={len(x)})')
    ax.set_xlabel('x')
    ax.set_ylabel('Pressure')
    ax.set_title('Pressure (raw particle data)')
    ax.legend()
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)

    # Internal energy - RAW particles
    ax = axes[1, 1]
    e_exact = P_exact / ((gamma - 1) * rho_exact)
    e_particles = pres / ((gamma - 1) * rho)
    ax.plot(x_exact, e_exact, 'r-', lw=2, label='Exact solution', zorder=10)
    ax.scatter(x, e_particles, s=1, c='blue', alpha=0.5, label=f'Particles (N={len(x)})')
    ax.set_xlabel('x')
    ax.set_ylabel('Internal Energy')
    ax.set_title('Internal Energy (raw particle data)')
    ax.legend()
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'3D Sod Shock Tube - GSPH HLL Wendland N=120 (t={t:.2f})\n'
                 f'Morton ordering + Iterative traversal enabled', fontsize=14)
    plt.tight_layout()

    output_file = data_dir / f"shock_tube_3d_raw_t{t:.2f}.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.show()


if __name__ == "__main__":
    main()
