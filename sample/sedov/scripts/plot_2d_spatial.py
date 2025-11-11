#!/usr/bin/env python3
"""
Create 2D spatial visualization of Sedov blast wave showing density, pressure, and velocity fields.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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
    Create 2D spatial visualization of Sedov blast wave.
    
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
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Common parameters
    extent = 0.5
    point_size = 50
    
    # 1. Density field
    ax = axes[0, 0]
    sc = ax.scatter(x, y, c=rho, s=point_size, cmap='viridis', 
                    vmin=0, vmax=rho.max())
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Density', fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('ρ', fontsize=12)
    
    # Add circular grid lines
    for r in [0.1, 0.2, 0.3, 0.4]:
        circle = plt.Circle((0, 0), r, fill=False, color='white', 
                           linestyle=':', linewidth=1, alpha=0.5)
        ax.add_patch(circle)
    
    # 2. Pressure field (log scale)
    ax = axes[0, 1]
    sc = ax.scatter(x, y, c=p, s=point_size, cmap='plasma',
                    norm=LogNorm(vmin=max(p.min(), 1e-7), vmax=p.max()))
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Pressure (log scale)', fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('p', fontsize=12)
    
    # Add circular grid lines
    for r in [0.1, 0.2, 0.3, 0.4]:
        circle = plt.Circle((0, 0), r, fill=False, color='white',
                           linestyle=':', linewidth=1, alpha=0.5)
        ax.add_patch(circle)
    
    # 3. Velocity magnitude
    ax = axes[1, 0]
    sc = ax.scatter(x, y, c=v_mag, s=point_size, cmap='coolwarm',
                    vmin=0, vmax=v_mag.max())
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Velocity Magnitude', fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('|v|', fontsize=12)
    
    # Add velocity vectors (subsample for clarity)
    skip = max(1, len(x) // 400)  # Show ~400 vectors
    ax.quiver(x[::skip], y[::skip], vx[::skip], vy[::skip],
              color='black', alpha=0.6, scale=5, width=0.003)
    
    # Add circular grid lines
    for r in [0.1, 0.2, 0.3, 0.4]:
        circle = plt.Circle((0, 0), r, fill=False, color='gray',
                           linestyle=':', linewidth=1, alpha=0.5)
        ax.add_patch(circle)
    
    # 4. Internal Energy (log scale)
    ax = axes[1, 1]
    sc = ax.scatter(x, y, c=e, s=point_size, cmap='hot',
                    norm=LogNorm(vmin=max(e.min(), 1e-7), vmax=e.max()))
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Internal Energy (log scale)', fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('e', fontsize=12)
    
    # Add circular grid lines
    for r in [0.1, 0.2, 0.3, 0.4]:
        circle = plt.Circle((0, 0), r, fill=False, color='white',
                           linestyle=':', linewidth=1, alpha=0.5)
        ax.add_patch(circle)
    
    # Overall title
    fig.suptitle(f'Sedov Blast Wave: 2D Spatial Distribution (t = {time:.4f})',
                 fontsize=15, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create 2D spatial visualization of Sedov blast wave')
    parser.add_argument('snapshot', type=str, help='Snapshot CSV file')
    parser.add_argument('-o', '--output', type=str, help='Output PNG file')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')
    
    args = parser.parse_args()
    
    plot_2d_spatial(args.snapshot, args.output, not args.no_show)

if __name__ == '__main__':
    main()
