#!/usr/bin/env python3
"""
Create 2D spatial visualization of shock tube showing density, pressure, velocity, and energy fields.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pandas as pd
import argparse
from pathlib import Path

def load_snapshot(filepath):
    """Load a snapshot CSV file."""
    metadata = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    if key == 'Time (code)':
                        metadata['Time'] = value.split()[0] if value else '0.0'
                    elif key == 'Gamma':
                        metadata['Gamma'] = value.split()[0] if value else '1.4'
                    else:
                        metadata[key] = value
            else:
                break
    
    data = pd.read_csv(filepath, comment='#')
    return data, metadata

def plot_2d_spatial(snapshot_file, output_file=None, show_plot=True):
    """
    Create 2D spatial visualization of shock tube.
    
    Parameters:
    -----------
    snapshot_file : str
        Path to SPH snapshot CSV file
    output_file : str, optional
        Path to save the plot
    show_plot : bool
        Whether to display the plot
    """
    # Load data
    data, metadata = load_snapshot(snapshot_file)
    
    # Extract time
    time = float(metadata.get('Time', 0.0))
    
    # Get positions
    x = data['pos_x'].values
    y = data['pos_y'].values
    
    # Get physical quantities
    rho = data['dens'].values
    p = data['pres'].values
    vx = data['vel_x'].values
    vy = data['vel_y'].values
    v_mag = np.sqrt(vx**2 + vy**2)
    e = data['ene'].values
    
    # Create figure with 2x2 layout
    fig, axes = plt.subplots(2, 2, figsize=(16, 8))
    
    # Common parameters
    x_extent = (0.0, 1.0)
    y_extent = (0.0, 0.2)
    point_size = 30
    
    # 1. Density field
    ax = axes[0, 0]
    sc = ax.scatter(x, y, c=rho, s=point_size, cmap='viridis', 
                    vmin=0, vmax=rho.max(), marker='o', edgecolors='none')
    ax.set_xlim(x_extent)
    ax.set_ylim(y_extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Density (t={time:.4f})', fontsize=13, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add shock discontinuity marker
    ax.axvline(x=0.5, color='red', linestyle=':', linewidth=2, alpha=0.5, label='Initial discontinuity')
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('ρ', fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    
    # 2. Pressure field
    ax = axes[0, 1]
    sc = ax.scatter(x, y, c=p, s=point_size, cmap='plasma',
                    vmin=0, vmax=p.max(), marker='o', edgecolors='none')
    ax.set_xlim(x_extent)
    ax.set_ylim(y_extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Pressure (t={time:.4f})', fontsize=13, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add shock discontinuity marker
    ax.axvline(x=0.5, color='red', linestyle=':', linewidth=2, alpha=0.5, label='Initial discontinuity')
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('p', fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    
    # 3. Velocity (x-component)
    ax = axes[1, 0]
    vx_max = max(abs(vx.min()), abs(vx.max()))
    sc = ax.scatter(x, y, c=vx, s=point_size, cmap='coolwarm',
                    vmin=-vx_max, vmax=vx_max, marker='o', edgecolors='none')
    ax.set_xlim(x_extent)
    ax.set_ylim(y_extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'X-Velocity (t={time:.4f})', fontsize=13, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add shock discontinuity marker
    ax.axvline(x=0.5, color='black', linestyle=':', linewidth=2, alpha=0.5, label='Initial discontinuity')
    
    # Add velocity vectors (subsample for clarity)
    skip = max(1, len(x) // 200)  # Show ~200 vectors
    ax.quiver(x[::skip], y[::skip], vx[::skip], vy[::skip],
              color='black', alpha=0.5, scale=5, width=0.002, headwidth=3)
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('vx', fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    
    # 4. Internal Energy
    ax = axes[1, 1]
    sc = ax.scatter(x, y, c=e, s=point_size, cmap='hot',
                    vmin=0, vmax=e.max(), marker='o', edgecolors='none')
    ax.set_xlim(x_extent)
    ax.set_ylim(y_extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Internal Energy (t={time:.4f})', fontsize=13, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add shock discontinuity marker
    ax.axvline(x=0.5, color='cyan', linestyle=':', linewidth=2, alpha=0.5, label='Initial discontinuity')
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('e', fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    
    # Overall title
    fig.suptitle(f'2D Shock Tube Spatial Distribution', fontsize=15, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create 2D spatial plot of shock tube')
    parser.add_argument('snapshot', help='Path to snapshot CSV file')
    parser.add_argument('-o', '--output', help='Output image file')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')
    
    args = parser.parse_args()
    
    plot_2d_spatial(
        args.snapshot,
        output_file=args.output,
        show_plot=not args.no_show
    )

if __name__ == '__main__':
    main()
