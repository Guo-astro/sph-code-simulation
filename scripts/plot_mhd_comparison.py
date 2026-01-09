#!/usr/bin/env python3
"""
Plot MHD simulation results with analytic/reference solution overlay.

Usage:
    python plot_mhd_comparison.py --simulation results/mhd_shock_tube_1/snapshot_0010.csv \
                                  --reference reference_1.csv \
                                  --output comparison.png
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import re


def read_simulation_csv(filepath):
    """Read simulation snapshot CSV file."""
    data = {}
    header_lines = 0

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines += 1
                # Extract time from header
                if 'Time (code):' in line:
                    match = re.search(r'Time \(code\):\s*([\d.eE+-]+)', line)
                    if match:
                        data['time'] = float(match.group(1))
            elif line.startswith('id,'):
                header_lines += 1
                columns = line.strip().split(',')
                break

    # Read particle data
    raw_data = np.genfromtxt(filepath, delimiter=',', skip_header=header_lines)

    # Map columns (1D case)
    data['x'] = raw_data[:, 1]   # pos_x
    data['vx'] = raw_data[:, 2]  # vel_x
    data['rho'] = raw_data[:, 4] if raw_data.shape[1] > 4 else raw_data[:, 3]  # dens
    data['P'] = raw_data[:, 5] if raw_data.shape[1] > 5 else raw_data[:, 4]    # pres

    # Sort by x position
    idx = np.argsort(data['x'])
    data['x'] = data['x'][idx]
    data['vx'] = data['vx'][idx]
    data['rho'] = data['rho'][idx]
    data['P'] = data['P'][idx]

    return data


def read_reference_csv(filepath):
    """Read reference solution CSV file."""
    data = {}
    header_lines = 0

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines += 1
            elif 'x,rho' in line:
                header_lines += 1
                break

    raw_data = np.genfromtxt(filepath, delimiter=',', skip_header=header_lines)

    data['x'] = raw_data[:, 0]
    data['rho'] = raw_data[:, 1]
    data['vx'] = raw_data[:, 2]
    data['P'] = raw_data[:, 3]
    data['Bx'] = raw_data[:, 4]
    data['By'] = raw_data[:, 5]

    return data


def plot_comparison(sim_data, ref_data, output_file, title="MHD Shock Tube"):
    """Create comparison plot of simulation vs reference."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Density
    ax = axes[0, 0]
    ax.plot(ref_data['x'], ref_data['rho'], 'k-', label='Reference', linewidth=2)
    ax.scatter(sim_data['x'], sim_data['rho'], c='r', s=3, alpha=0.7, label='Simulation')
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\rho$')
    ax.set_title('Density')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Velocity
    ax = axes[0, 1]
    ax.plot(ref_data['x'], ref_data['vx'], 'k-', label='Reference', linewidth=2)
    ax.scatter(sim_data['x'], sim_data['vx'], c='r', s=3, alpha=0.7, label='Simulation')
    ax.set_xlabel('x')
    ax.set_ylabel(r'$v_x$')
    ax.set_title('Velocity')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Pressure
    ax = axes[1, 0]
    ax.plot(ref_data['x'], ref_data['P'], 'k-', label='Reference', linewidth=2)
    ax.scatter(sim_data['x'], sim_data['P'], c='r', s=3, alpha=0.7, label='Simulation')
    ax.set_xlabel('x')
    ax.set_ylabel('P')
    ax.set_title('Pressure')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Magnetic field By
    ax = axes[1, 1]
    ax.plot(ref_data['x'], ref_data['By'], 'k-', label='Reference', linewidth=2)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$B_y$')
    ax.set_title('Magnetic Field By')
    ax.legend()
    ax.grid(True, alpha=0.3)

    time_str = f"t = {sim_data.get('time', 'N/A'):.4f}" if 'time' in sim_data else ""
    fig.suptitle(f'{title} {time_str}', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Comparison plot saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Plot MHD simulation vs reference solution')
    parser.add_argument('--simulation', '-s', type=str, required=True,
                       help='Path to simulation snapshot CSV')
    parser.add_argument('--reference', '-r', type=str, required=True,
                       help='Path to reference solution CSV')
    parser.add_argument('--output', '-o', type=str, default='mhd_comparison.png',
                       help='Output plot file')
    parser.add_argument('--title', '-t', type=str, default='MHD Shock Tube',
                       help='Plot title')

    args = parser.parse_args()

    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        return 1

    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}")
        return 1

    print(f"Reading simulation data from: {args.simulation}")
    sim_data = read_simulation_csv(args.simulation)

    print(f"Reading reference data from: {args.reference}")
    ref_data = read_reference_csv(args.reference)

    plot_comparison(sim_data, ref_data, args.output, args.title)

    return 0


if __name__ == '__main__':
    exit(main())
