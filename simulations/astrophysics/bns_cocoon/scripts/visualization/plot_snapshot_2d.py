#!/usr/bin/env python3
"""
BNS Cocoon 2D Snapshot Visualization

Generates scatter plots and density maps from 2D simulation output.
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm


def load_csv(filename):
    """Load CSV snapshot with header comments."""
    import pandas as pd
    # Read skipping comment lines
    data = pd.read_csv(filename, comment='#')
    return data


def plot_snapshot(data, output_path, snapshot_num, field='dens'):
    """Plot a single snapshot as a scatter plot."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    x = data['pos_x']
    y = data['pos_y']
    
    # Density
    ax = axes[0]
    sc = ax.scatter(x, y, c=data['dens'], s=1, cmap='viridis', 
                    norm=LogNorm(vmin=max(data['dens'].min(), 1e-12), 
                                 vmax=data['dens'].max()))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Density (t = {snapshot_num})')
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax, label=r'$\rho$')
    
    # Pressure
    ax = axes[1]
    pres = np.abs(data['pres']) + 1e-15
    sc = ax.scatter(x, y, c=pres, s=1, cmap='hot',
                    norm=LogNorm(vmin=max(pres.min(), 1e-15), 
                                 vmax=pres.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Pressure')
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax, label=r'$P$')
    
    # Velocity magnitude
    ax = axes[2]
    vx = data['vel_x']
    vy = data['vel_y']
    vmag = np.sqrt(vx**2 + vy**2)
    sc = ax.scatter(x, y, c=vmag, s=1, cmap='coolwarm', vmin=0, vmax=min(1, vmag.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'|v|/c')
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax, label=r'$|v|/c$')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Plot 2D BNS Cocoon snapshots')
    parser.add_argument('--input', '-i', required=True, 
                        help='Input directory with CSV files')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory for plots')
    parser.add_argument('--pattern', default='snapshot_*.csv',
                        help='File pattern to match')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    # Find all snapshot files
    pattern = os.path.join(args.input, args.pattern)
    files = sorted(glob.glob(pattern))
    
    if not files:
        print(f"No files found matching {pattern}")
        sys.exit(1)
    
    print(f"Found {len(files)} snapshots")
    
    for f in files:
        basename = os.path.basename(f)
        num = basename.replace('snapshot_', '').replace('.csv', '')
        output_path = os.path.join(args.output, f'snapshot_{num}.png')
        
        try:
            data = load_csv(f)
            plot_snapshot(data, output_path, num)
        except Exception as e:
            print(f"  Error processing {f}: {e}")
    
    print(f"\nPlots saved to {args.output}")


if __name__ == '__main__':
    main()
