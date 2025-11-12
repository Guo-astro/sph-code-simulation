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

def exact_riemann_solver(x, t, gamma=1.4, x0=0.5, rho_L=1.0, rho_R=0.125, 
                         p_L=1.0, p_R=0.1, u_L=0.0, u_R=0.0):
    """
    Exact Riemann solution for Sod shock tube problem.
    
    Initial conditions:
        Left state (x < x0):  rho=1.0, p=1.0, u=0.0
        Right state (x > x0): rho=0.125, p=0.1, u=0.0
    
    Returns: rho, u, p, e (density, velocity, pressure, specific internal energy)
    """
    # Sound speeds
    c_L = np.sqrt(gamma * p_L / rho_L)
    c_R = np.sqrt(gamma * p_R / rho_R)
    
    # Solve for pressure in star region (using iterative method)
    # For simplicity, use approximate formula
    p_star = 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L) * (rho_L + rho_R) * (c_L + c_R)
    p_star = max(p_star, 0.001 * min(p_L, p_R))  # Ensure positive
    
    # Better approximation via Newton iteration
    for _ in range(10):
        f_L = (p_star / p_L - 1.0) * np.sqrt(rho_L / (gamma * ((gamma + 1) / 2 * p_star + (gamma - 1) / 2 * p_L)))
        f_R = (p_star / p_R - 1.0) * np.sqrt(rho_R / (gamma * ((gamma + 1) / 2 * p_star + (gamma - 1) / 2 * p_R)))
        f = f_L + f_R + (u_R - u_L)
        
        df_L = np.sqrt(rho_L / (gamma * ((gamma + 1) / 2 * p_star + (gamma - 1) / 2 * p_L))) * \
               (1.0 / p_L - 0.5 * (gamma + 1) * (p_star + p_L * (gamma - 1) / (gamma + 1)) / \
                (p_star * ((gamma + 1) / 2 * p_star + (gamma - 1) / 2 * p_L)))
        df_R = np.sqrt(rho_R / (gamma * ((gamma + 1) / 2 * p_star + (gamma - 1) / 2 * p_R))) * \
               (1.0 / p_R - 0.5 * (gamma + 1) * (p_star + p_R * (gamma - 1) / (gamma + 1)) / \
                (p_star * ((gamma + 1) / 2 * p_star + (gamma - 1) / 2 * p_R)))
        df = df_L + df_R
        
        p_star = p_star - f / df
        p_star = max(p_star, 0.001 * min(p_L, p_R))
    
    # Velocity in star region
    u_star = 0.5 * (u_L + u_R) + 0.5 * (f_L - f_R)
    
    # Densities in star regions
    rho_star_L = rho_L * (p_star / p_L + (gamma - 1) / (gamma + 1)) / \
                 ((gamma - 1) / (gamma + 1) * p_star / p_L + 1.0)
    rho_star_R = rho_R * (p_star / p_R) ** (1.0 / gamma)
    
    # Wave speeds
    S_L = u_L - c_L * np.sqrt((gamma + 1) / (2 * gamma) * p_star / p_L + (gamma - 1) / (2 * gamma))
    S_R = u_R + c_R * np.sqrt((gamma + 1) / (2 * gamma) * p_star / p_R + (gamma - 1) / (2 * gamma))
    S_contact = u_star
    S_head_fan = u_L - c_L
    c_star_L = c_L * (p_star / p_L) ** ((gamma - 1) / (2 * gamma))
    S_tail_fan = u_star - c_star_L
    
    # Initialize solution arrays
    n = len(x)
    rho = np.zeros(n)
    u = np.zeros(n)
    p = np.zeros(n)
    
    # Position relative to initial discontinuity
    x_rel = (x - x0) / t if t > 0 else x - x0
    
    for i in range(n):
        xi = x_rel[i] if t > 0 else (1e10 if x[i] > x0 else -1e10)
        
        if xi < S_L:
            # Left state
            rho[i] = rho_L
            u[i] = u_L
            p[i] = p_L
        elif xi < S_head_fan:
            # Left state (before rarefaction)
            rho[i] = rho_L
            u[i] = u_L
            p[i] = p_L
        elif xi < S_tail_fan:
            # Rarefaction fan
            u[i] = 2.0 / (gamma + 1) * (c_L + (gamma - 1) / 2 * u_L + xi)
            c = c_L + (gamma - 1) / 2 * (u_L - xi)
            rho[i] = rho_L * (c / c_L) ** (2.0 / (gamma - 1))
            p[i] = p_L * (c / c_L) ** (2.0 * gamma / (gamma - 1))
        elif xi < S_contact:
            # Left star state
            rho[i] = rho_star_L
            u[i] = u_star
            p[i] = p_star
        elif xi < S_R:
            # Right star state
            rho[i] = rho_star_R
            u[i] = u_star
            p[i] = p_star
        else:
            # Right state
            rho[i] = rho_R
            u[i] = u_R
            p[i] = p_R
    
    # Specific internal energy: e = p / ((gamma - 1) * rho)
    e = p / ((gamma - 1) * rho)
    
    return rho, u, p, e


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
    rho_exact, u_exact, p_exact, e_exact = exact_riemann_solver(
        x_exact, time, gamma=gamma, x0=0.5,
        rho_L=1.0, rho_R=0.125, p_L=1.0, p_R=0.1
    )
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'2D Shock Tube Comparison (t = {time:.4f})', fontsize=14, fontweight='bold')
    
    # Density
    ax = axes[0, 0]
    ax.plot(x_exact, rho_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], rho_avg[mask], c='red', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.set_xlabel('x')
    ax.set_ylabel('Density')
    ax.set_title('Density Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Velocity
    ax = axes[0, 1]
    ax.plot(x_exact, u_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], vx_avg[mask], c='blue', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.set_xlabel('x')
    ax.set_ylabel('Velocity')
    ax.set_title('Velocity Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Pressure
    ax = axes[1, 0]
    ax.plot(x_exact, p_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], pres_avg[mask], c='green', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.set_xlabel('x')
    ax.set_ylabel('Pressure')
    ax.set_title('Pressure Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Internal Energy
    ax = axes[1, 1]
    ax.plot(x_exact, e_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
    ax.scatter(x_centers[mask], ene_avg[mask], c='purple', s=20, alpha=0.6, label='SPH (y-averaged)', zorder=5)
    ax.set_xlabel('x')
    ax.set_ylabel('Specific Internal Energy')
    ax.set_title('Internal Energy Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
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
