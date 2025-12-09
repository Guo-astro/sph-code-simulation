#!/usr/bin/env python3
"""
Generate comparison animation for 2D planar Lane-Emden slab grad-h study.

Creates a 2x2 animated GIF showing all 4 SPH methods side by side:
- GSPH + grad-h
- GSPH - grad-h  
- SSPH + grad-h
- SSPH - grad-h

Uses SSOT module from scripts.shared.lane_emden for Lane-Emden solutions.
"""

import sys
from pathlib import Path

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize
import os
import glob

from scripts.shared.lane_emden import solve_lane_emden_planar as _solve_lane_emden_planar

# Configuration
BASE_DIR = "sample/gradh_study_2d_planar/results"
OUTPUT_DIR = "sample/gradh_study_2d_planar/results/comparison"

# Methods to compare
METHODS = {
    'gsph_gradh': {'label': 'GSPH + grad-h', 'color': '#0066CC'},
    'gsph_nogradh': {'label': 'GSPH - grad-h', 'color': '#9933FF'},
    'ssph_gradh': {'label': 'SSPH + grad-h', 'color': '#CC0000'},
    'ssph_nogradh': {'label': 'SSPH - grad-h', 'color': '#FF6600'},
}

# Physical parameters
RHO_CENTER = 1.0
K = 1.0
G = 1.0
GAMMA = 5.0/3.0


def get_surface_position():
    """Get analytical surface position using SSOT."""
    n = 1.0 / (GAMMA - 1.0)
    alpha_sq = K * (n + 1.0) * RHO_CENTER**(1.0 - n) / (2.0 * np.pi * G)
    alpha = np.sqrt(alpha_sq)
    
    xi_le, theta_le = _solve_lane_emden_planar(n, xi_max=10.0, n_points=10000)
    theta_le = np.maximum(theta_le, 0)  # θ ≥ 0
    
    # Find surface (where theta = 0)
    surface_idx = np.argmax(theta_le <= 0) if np.any(theta_le <= 0) else len(theta_le) - 1
    xi_surface = xi_le[surface_idx]
    
    return alpha * xi_surface


def read_snapshot(directory, snapshot_num):
    """Read a snapshot CSV file."""
    filepath = os.path.join(directory, f"snapshot_{snapshot_num:04d}.csv")
    if not os.path.exists(filepath):
        return None
    return pd.read_csv(filepath, comment='#')


def count_snapshots(directory):
    """Count available snapshots."""
    pattern = os.path.join(directory, "snapshot_*.csv")
    return len(glob.glob(pattern))


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 70)
    print("ANIMATION: 2D Planar Lane-Emden Slab Comparison")
    print("=" * 70)
    
    # Find available methods
    available_methods = []
    for method_name in METHODS:
        method_dir = os.path.join(BASE_DIR, method_name)
        if os.path.exists(method_dir) and count_snapshots(method_dir) > 0:
            available_methods.append(method_name)
            print(f"  Found: {method_name} ({count_snapshots(method_dir)} snapshots)")
    
    if len(available_methods) < 2:
        print("ERROR: Need at least 2 methods with data for comparison")
        sys.exit(1)
    
    # Get number of snapshots (use minimum)
    num_snapshots = min(count_snapshots(os.path.join(BASE_DIR, m)) for m in available_methods)
    print(f"\nUsing {num_snapshots} snapshots")
    
    # Get surface position
    y_surface = get_surface_position()
    print(f"Analytical surface: y = ±{y_surface:.3f}")
    
    # Setup figure
    n_methods = len(available_methods)
    n_cols = 2
    n_rows = (n_methods + 1) // 2
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 5*n_rows))
    if n_methods == 1:
        axes = np.array([[axes]])
    elif n_methods == 2:
        axes = axes.reshape(1, 2)
    axes = axes.flatten()
    
    # Get data range from first snapshot
    first_snap = read_snapshot(os.path.join(BASE_DIR, available_methods[0]), 0)
    x_range = (first_snap['pos_x'].min() - 0.5, first_snap['pos_x'].max() + 0.5)
    y_range = (-y_surface * 1.5, y_surface * 1.5)
    
    # Density colormap range
    vmin, vmax = 0.0, 1.2 * RHO_CENTER
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Output interval
    output_interval = 0.5
    
    # Create scatter plots
    scatters = []
    titles = []
    
    for i, method_name in enumerate(available_methods):
        ax = axes[i]
        snap = read_snapshot(os.path.join(BASE_DIR, method_name), 0)
        
        sc = ax.scatter(snap['pos_x'], snap['pos_y'], c=snap['dens'],
                       cmap='viridis', norm=norm, s=5, alpha=0.8)
        scatters.append(sc)
        
        # Surface lines
        ax.axhline(y=y_surface, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.axhline(y=-y_surface, color='red', linestyle='--', alpha=0.5, linewidth=1)
        
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        title = ax.set_title(f"{METHODS[method_name]['label']}\nt = 0.0")
        titles.append(title)
        ax.set_aspect('equal')
    
    # Hide unused axes
    for i in range(len(available_methods), len(axes)):
        axes[i].axis('off')
    
    # Add colorbar
    fig.colorbar(scatters[0], ax=axes[:len(available_methods)], label='Density', shrink=0.8)
    
    fig.suptitle('2D Planar Lane-Emden Slab: GSPH/SSPH × grad-h Comparison\n(γ=5/3, Wendland kernel, 1D gravity in y)', fontsize=14)
    
    def update(frame):
        for i, method_name in enumerate(available_methods):
            snap = read_snapshot(os.path.join(BASE_DIR, method_name), frame)
            if snap is not None:
                # Update scatter data
                scatters[i].set_offsets(np.c_[snap['pos_x'], snap['pos_y']])
                scatters[i].set_array(snap['dens'])
                
                # Update title
                time = frame * output_interval
                rho_max = snap['dens'].max()
                titles[i].set_text(f"{METHODS[method_name]['label']}\nt = {time:.1f}, ρ_max = {rho_max:.3f}")
        
        return scatters + titles
    
    print("\nGenerating animation...")
    anim = FuncAnimation(fig, update, frames=num_snapshots, interval=100, blit=True)
    
    output_file = os.path.join(OUTPUT_DIR, '4method_comparison.gif')
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer)
    plt.close()
    
    print(f"\n✓ Animation saved: {output_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
