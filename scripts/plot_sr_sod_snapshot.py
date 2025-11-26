#!/usr/bin/env python3
"""
Plot SR Sod shock tube snapshot
Shows all primitive variables with analytical solution overlay
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_snapshot(filename):
    """Load CSV snapshot and extract particle data (excluding ghost particles)"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.strip().startswith('id'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 22:  # Updated to include is_ghost column
                try:
                    row = [float(x) for x in parts]
                    # Filter out ghost particles (is_ghost is column 21, 0-indexed)
                    if len(row) > 21 and row[21] == 0:  # is_ghost == 0 means real particle
                        data.append(row)
                except ValueError:
                    continue

    data = np.array(data)

    result = {
        'pos_x': data[:, 1],
        'vel_x': data[:, 4],
        'nu': data[:, 10],
        'dens': data[:, 11],    # rest-frame density n
        'pres': data[:, 12],
        'ene': data[:, 13],
        'sml': data[:, 14],
        'sound': data[:, 15]
    }

    # Sort by position
    idx = np.argsort(result['pos_x'])
    for key in result:
        result[key] = result[key][idx]

    return result

def compute_derived(data):
    """Compute derived quantities"""
    c = 1.0
    v = data['vel_x']
    gamma = 1.0 / np.sqrt(1.0 - v**2 / c**2)
    N = gamma * data['dens']  # Lab-frame density

    return gamma, N

def get_analytical_solution(t, gamma_eos=5.0/3.0):
    """
    Get analytical solution for SR Sod test
    Using exact Riemann solver
    """
    # Initial conditions (Kitajima et al.)
    P_L, n_L = 1.0, 1.0
    P_R, n_R = 0.1, 0.125

    # For now, return placeholder
    # TODO: Implement exact Riemann solver or load from file
    return None

def plot_snapshot(data, gamma_lor, N, output_file, time=0.0):
    """Create comprehensive plot of snapshot"""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    x = data['pos_x']

    # Plot 1: Rest-frame density n
    ax = axes[0, 0]
    ax.plot(x, data['dens'], 'ko-', markersize=3, linewidth=1.5, label='SPH')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Rest-frame Density n')
    ax.set_title(f'Rest-Frame Density (t={time:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Lab-frame density N = γn
    ax = axes[0, 1]
    ax.plot(x, N, 'ro-', markersize=3, linewidth=1.5, label='SPH: N = γn')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Lab-frame Density N')
    ax.set_title(f'Lab-Frame Density (t={time:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Pressure
    ax = axes[0, 2]
    ax.plot(x, data['pres'], 'bo-', markersize=3, linewidth=1.5, label='SPH')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Pressure P')
    ax.set_title(f'Pressure (t={time:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Velocity
    ax = axes[1, 0]
    ax.plot(x, data['vel_x'], 'go-', markersize=3, linewidth=1.5, label='SPH')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Velocity v')
    ax.set_title(f'Velocity (t={time:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 5: Lorentz factor
    ax = axes[1, 1]
    ax.plot(x, gamma_lor, 'mo-', markersize=3, linewidth=1.5, label='γ')
    ax.axhline(1, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Lorentz Factor γ')
    ax.set_title(f'Lorentz Factor (t={time:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 6: Specific internal energy
    ax = axes[1, 2]
    ax.plot(x, data['ene'], 'co-', markersize=3, linewidth=1.5, label='SPH')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Internal Energy u')
    ax.set_title(f'Internal Energy (t={time:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_sr_sod_snapshot.py <snapshot.csv> [output.png]")
        sys.exit(1)

    snapshot_file = sys.argv[1]

    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = snapshot_file.replace('.csv', '.png')

    # Extract time from filename or file
    time = 0.0
    with open(snapshot_file, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                time = float(line.split(':')[1].strip())
                break

    data = load_snapshot(snapshot_file)
    gamma_lor, N = compute_derived(data)

    plot_snapshot(data, gamma_lor, N, output_file, time)

if __name__ == "__main__":
    main()
