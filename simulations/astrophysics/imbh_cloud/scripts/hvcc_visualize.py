#!/usr/bin/env python3
"""
HVCC 10K Isothermal Simulation Visualization
Generates animation and final plots from CSV snapshots
"""

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.animation as animation

def load_snapshot(filename):
    """Load CSV snapshot file"""
    try:
        data = np.genfromtxt(filename, delimiter=',', skip_header=1, 
                            dtype=None, encoding='utf-8', names=True)
        return data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def get_snapshot_files(output_dir):
    """Get sorted list of snapshot files"""
    pattern = os.path.join(output_dir, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))
    # Exclude snapshot_final.csv if it exists (it's a copy)
    files = [f for f in files if 'final' not in f]
    return files

def create_animation(output_dir, output_file="animation.gif"):
    """Create animation from snapshots"""
    files = get_snapshot_files(output_dir)
    if not files:
        print(f"No snapshot files found in {output_dir}")
        return False
    
    print(f"Found {len(files)} snapshots")
    
    # Load first snapshot to get dimensions
    data0 = load_snapshot(files[0])
    if data0 is None:
        return False
    
    # Get field names
    field_names = data0.dtype.names
    
    # Setup figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Determine axis limits from first snapshot
    try:
        mask = data0['is_ghost'] == 0 if 'is_ghost' in field_names else np.ones(len(data0), dtype=bool)
    except:
        mask = np.ones(len(data0), dtype=bool)
    
    x0 = data0['pos_x'][mask] if 'pos_x' in field_names else data0['x'][mask]
    y0 = data0['pos_y'][mask] if 'pos_y' in field_names else data0['y'][mask]
    
    xlim = (x0.min() - 0.05, x0.max() + 0.05)
    ylim = (y0.min() - 0.05, y0.max() + 0.05)
    
    def animate(frame_idx):
        for ax in axes:
            ax.clear()
        
        data = load_snapshot(files[frame_idx])
        if data is None:
            return
        
        field_names = data.dtype.names
        
        # Filter non-ghost particles
        try:
            mask = data['is_ghost'] == 0 if 'is_ghost' in field_names else np.ones(len(data), dtype=bool)
        except:
            mask = np.ones(len(data), dtype=bool)
        
        x = data['pos_x'][mask] if 'pos_x' in field_names else data['x'][mask]
        y = data['pos_y'][mask] if 'pos_y' in field_names else data['y'][mask]
        z = data['pos_z'][mask] if 'pos_z' in field_names else data['z'][mask]
        dens = data['dens'][mask] if 'dens' in field_names else data['density'][mask]
        
        # XY projection
        sc = axes[0].scatter(x, y, c=dens, s=1, cmap='viridis', norm=LogNorm(vmin=dens.min()+0.01, vmax=dens.max()))
        axes[0].set_xlabel('X [pc]')
        axes[0].set_ylabel('Y [pc]')
        axes[0].set_title('XY Projection')
        axes[0].set_aspect('equal')
        axes[0].set_xlim(xlim)
        axes[0].set_ylim(ylim)
        
        # XZ projection  
        axes[1].scatter(x, z, c=dens, s=1, cmap='viridis', norm=LogNorm(vmin=dens.min()+0.01, vmax=dens.max()))
        axes[1].set_xlabel('X [pc]')
        axes[1].set_ylabel('Z [pc]')
        axes[1].set_title('XZ Projection')
        axes[1].set_aspect('equal')
        axes[1].set_xlim(xlim)
        axes[1].set_ylim(ylim)
        
        # Radial density profile
        r = np.sqrt(x**2 + y**2 + z**2)
        axes[2].scatter(r, dens, s=1, alpha=0.5)
        axes[2].set_xlabel('Radius [pc]')
        axes[2].set_ylabel(r'$\rho$ [code]')
        axes[2].set_title('Radial Density Profile')
        axes[2].set_yscale('log')
        axes[2].grid(True, alpha=0.3)
        
        fig.suptitle(f'HVCC 10K - Frame {frame_idx+1}/{len(files)}', fontsize=14)
        plt.tight_layout()
    
    print("Creating animation...")
    anim = animation.FuncAnimation(fig, animate, frames=len(files), interval=200)
    
    output_path = os.path.join(output_dir, output_file)
    try:
        anim.save(output_path, writer='pillow', fps=5)
        print(f"Animation saved: {output_path}")
    except Exception as e:
        print(f"Warning: Could not save animation: {e}")
    plt.close()
    return True

def create_final_plot(output_dir, output_file="final_plot.png"):
    """Create final state plot"""
    files = get_snapshot_files(output_dir)
    if not files:
        print(f"No snapshot files found in {output_dir}")
        return False
    
    # Load final snapshot
    data = load_snapshot(files[-1])
    if data is None:
        return False
    
    field_names = data.dtype.names
    
    # Filter non-ghost particles
    try:
        mask = data['is_ghost'] == 0 if 'is_ghost' in field_names else np.ones(len(data), dtype=bool)
    except:
        mask = np.ones(len(data), dtype=bool)
    
    x = data['pos_x'][mask] if 'pos_x' in field_names else data['x'][mask]
    y = data['pos_y'][mask] if 'pos_y' in field_names else data['y'][mask]
    z = data['pos_z'][mask] if 'pos_z' in field_names else data['z'][mask]
    dens = data['dens'][mask] if 'dens' in field_names else data['density'][mask]
    pres = data['pres'][mask] if 'pres' in field_names else np.ones_like(dens)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # XY projection colored by density
    sc = axes[0, 0].scatter(x, y, c=dens, s=2, cmap='viridis', norm=LogNorm())
    axes[0, 0].set_xlabel('X [pc]')
    axes[0, 0].set_ylabel('Y [pc]')
    axes[0, 0].set_title('XY Projection (Density)')
    axes[0, 0].set_aspect('equal')
    plt.colorbar(sc, ax=axes[0, 0], label=r'$\rho$ [code]')
    
    # XZ projection colored by density
    sc = axes[0, 1].scatter(x, z, c=dens, s=2, cmap='viridis', norm=LogNorm())
    axes[0, 1].set_xlabel('X [pc]')
    axes[0, 1].set_ylabel('Z [pc]')
    axes[0, 1].set_title('XZ Projection (Density)')
    axes[0, 1].set_aspect('equal')
    plt.colorbar(sc, ax=axes[0, 1], label=r'$\rho$ [code]')
    
    # Radial density profile
    r = np.sqrt(x**2 + y**2 + z**2)
    r_bins = np.linspace(0, r.max(), 30)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    dens_mean = []
    dens_std = []
    for i in range(len(r_bins)-1):
        in_bin = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(in_bin) > 0:
            dens_mean.append(np.mean(dens[in_bin]))
            dens_std.append(np.std(dens[in_bin]))
        else:
            dens_mean.append(np.nan)
            dens_std.append(np.nan)
    
    axes[1, 0].errorbar(r_centers, dens_mean, yerr=dens_std, fmt='o-', capsize=3)
    axes[1, 0].set_xlabel('Radius [pc]')
    axes[1, 0].set_ylabel(r'$\rho$ [code units]')
    axes[1, 0].set_title('Radial Density Profile')
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Pressure vs density
    axes[1, 1].scatter(dens, pres, s=1, alpha=0.5)
    axes[1, 1].set_xlabel(r'$\rho$ [code units]')
    axes[1, 1].set_ylabel('P [code units]')
    axes[1, 1].set_title('Pressure vs Density')
    axes[1, 1].set_xscale('log')
    axes[1, 1].set_yscale('log')
    axes[1, 1].grid(True, alpha=0.3)
    
    fig.suptitle(f'HVCC 10K Final State - {len(x)} particles', fontsize=14)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, output_file)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Final plot saved: {output_path}")
    plt.close()
    return True

def main():
    if len(sys.argv) < 2:
        print("Usage: python hvcc_visualize.py <output_dir> [animate|plot|all]")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    action = sys.argv[2] if len(sys.argv) > 2 else "all"
    
    if not os.path.isdir(output_dir):
        print(f"Directory not found: {output_dir}")
        sys.exit(1)
    
    if action in ["animate", "all"]:
        create_animation(output_dir)
    
    if action in ["plot", "all"]:
        create_final_plot(output_dir)

if __name__ == "__main__":
    main()
