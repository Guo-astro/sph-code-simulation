#!/usr/bin/env python3
"""Plot 3D shock tube results with analytical Sod solution overlay."""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


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
    """Exact solution for the Sod shock tube problem."""
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
    # Configuration - use HLL results
    data_dir = Path("simulations/benchmarks/shock_tube/results/gsph_3d_hll")
    csv_file = data_dir / "snapshot_0020.csv"
    print(f"Reading: {csv_file}")

    data = read_csv_sph(csv_file)
    x, vx, rho, pres = data['pos_x'], data['vel_x'], data['dens'], data['pres']

    print(f"Total particles: {len(x)}")
    print(f"x range: [{x.min():.3f}, {x.max():.3f}]")

    # Bin particles by x and compute averages
    x_bins = np.linspace(0, 1, 101)
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    rho_avg, vx_avg, pres_avg = [], [], []

    for i in range(len(x_bins) - 1):
        mask = (x >= x_bins[i]) & (x < x_bins[i+1])
        if np.sum(mask) > 0:
            rho_avg.append(np.mean(rho[mask]))
            vx_avg.append(np.mean(vx[mask]))
            pres_avg.append(np.mean(pres[mask]))
        else:
            rho_avg.append(np.nan)
            vx_avg.append(np.nan)
            pres_avg.append(np.nan)

    rho_avg, vx_avg, pres_avg = np.array(rho_avg), np.array(vx_avg), np.array(pres_avg)

    # Analytical solution
    x_exact = np.linspace(0, 1, 1000)
    rho_exact, u_exact, P_exact = sod_exact_solution(x_exact, 0.2)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Density
    ax = axes[0, 0]
    ax.plot(x_exact, rho_exact, 'k-', lw=2, label='Exact')
    ax.scatter(x, rho, s=3, c='blue', alpha=0.3, label='Particles')
    ax.plot(x_centers, rho_avg, 'r-', lw=1.5, label='Binned avg')
    ax.set_xlabel('x'); ax.set_ylabel('Density'); ax.set_title('Density')
    ax.legend(); ax.set_xlim(0, 1); ax.grid(True, alpha=0.3)

    # Velocity
    ax = axes[0, 1]
    ax.plot(x_exact, u_exact, 'k-', lw=2, label='Exact')
    ax.scatter(x, vx, s=3, c='blue', alpha=0.3, label='Particles')
    ax.plot(x_centers, vx_avg, 'r-', lw=1.5, label='Binned avg')
    ax.set_xlabel('x'); ax.set_ylabel('Velocity'); ax.set_title('Velocity')
    ax.legend(); ax.set_xlim(0, 1); ax.grid(True, alpha=0.3)

    # Pressure
    ax = axes[1, 0]
    ax.plot(x_exact, P_exact, 'k-', lw=2, label='Exact')
    ax.scatter(x, pres, s=3, c='blue', alpha=0.3, label='Particles')
    ax.plot(x_centers, pres_avg, 'r-', lw=1.5, label='Binned avg')
    ax.set_xlabel('x'); ax.set_ylabel('Pressure'); ax.set_title('Pressure')
    ax.legend(); ax.set_xlim(0, 1); ax.grid(True, alpha=0.3)

    # Internal energy
    ax = axes[1, 1]
    gamma = 1.4
    e_exact = P_exact / ((gamma - 1) * rho_exact)
    e_particles = pres / ((gamma - 1) * rho)
    e_avg = pres_avg / ((gamma - 1) * rho_avg)
    ax.plot(x_exact, e_exact, 'k-', lw=2, label='Exact')
    ax.scatter(x, e_particles, s=3, c='blue', alpha=0.3, label='Particles')
    ax.plot(x_centers, e_avg, 'r-', lw=1.5, label='Binned avg')
    ax.set_xlabel('x'); ax.set_ylabel('Internal Energy'); ax.set_title('Internal Energy')
    ax.legend(); ax.set_xlim(0, 1); ax.grid(True, alpha=0.3)

    plt.suptitle('3D Sod Shock Tube - GSPH HLL (t=0.2, periodic BC)', fontsize=14)
    plt.tight_layout()

    output_file = data_dir / "shock_tube_3d_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.show()


if __name__ == "__main__":
    main()
