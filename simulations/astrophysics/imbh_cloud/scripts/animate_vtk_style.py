#!/usr/bin/env python3
"""
VTK-Style Animation for IMBH-Cloud Flyby Simulations

Creates professional scientific visualizations with:
- 3D-like scatter plots with proper depth shading
- Scientific colormaps (viridis, plasma, inferno)
- Professional axis styling
- Black hole marker with gravitational influence ring
- Time and statistics annotations

Usage:
    python animate_vtk_style.py <results_dir> <output_gif> <title> [--field velocity|temperature|density]
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
from matplotlib.collections import PathCollection
import matplotlib.animation as animation
import numpy as np
import pandas as pd

# Professional dark style for VTK-like appearance
plt.style.use('dark_background')


def load_snapshot(filepath: str) -> pd.DataFrame:
    """Load a snapshot CSV file with comment header."""
    # Find the header line (starts with 'id,')
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('id,'):
                header_line = i
                break
        else:
            raise ValueError(f"Could not find header line in {filepath}")
    
    df = pd.read_csv(filepath, skiprows=header_line)
    return df


def get_snapshot_files(results_dir: str) -> list:
    """Get sorted list of snapshot files."""
    pattern = os.path.join(results_dir, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))
    return files


def compute_velocity_magnitude(df: pd.DataFrame) -> np.ndarray:
    """Compute velocity magnitude in km/s."""
    vx = df['vel_x'].values
    vy = df['vel_y'].values
    if 'vel_z' in df.columns:
        vz = df['vel_z'].values
        return np.sqrt(vx**2 + vy**2 + vz**2)
    return np.sqrt(vx**2 + vy**2)


def compute_temperature(df: pd.DataFrame, gamma: float = 5/3, mu: float = 1.27) -> np.ndarray:
    """Compute temperature from internal energy."""
    k_B = 1.38e-16  # erg/K
    m_H = 1.67e-24  # g
    u = df['ene'].values  # code units (km/s)^2
    u_cgs = u * 1e10  # erg/g
    T = (gamma - 1) * mu * m_H * u_cgs / k_B
    return T


def compute_density_cgs(df: pd.DataFrame) -> np.ndarray:
    """Convert density to number density (cm^-3)."""
    M_sun = 1.989e33  # g
    pc = 3.086e18  # cm
    mu = 1.27
    m_H = 1.67e-24  # g
    
    rho = df['dens'].values  # M_sun/pc^3
    rho_cgs = rho * M_sun / pc**3  # g/cm^3
    n = rho_cgs / (mu * m_H)  # cm^-3
    return n


def create_vtk_frame(ax, df, field='velocity', time=0, title='', bh_pos=(0,0)):
    """Create a single VTK-style frame."""
    ax.clear()
    
    x = df['pos_x'].values
    y = df['pos_y'].values
    
    # Compute field values
    if field == 'velocity':
        values = compute_velocity_magnitude(df)
        cmap = 'plasma'
        vmin, vmax = 0, 150
        cbar_label = 'Velocity [km/s]'
    elif field == 'temperature':
        values = compute_temperature(df)
        cmap = 'inferno'
        vmin, vmax = 5, 100
        cbar_label = 'Temperature [K]'
    elif field == 'density':
        values = compute_density_cgs(df)
        cmap = 'viridis'
        vmin, vmax = 1, 1000
        cbar_label = 'Number Density [cm⁻³]'
    else:
        values = compute_velocity_magnitude(df)
        cmap = 'plasma'
        vmin, vmax = 0, 150
        cbar_label = 'Velocity [km/s]'
    
    # Use log scale for density
    if field == 'density':
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    # Sort by field value for proper depth rendering (low values behind)
    sort_idx = np.argsort(values)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]
    values_sorted = values[sort_idx]
    
    # Compute alpha based on local density for VTK-like transparency
    alpha_base = 0.6
    
    # Main scatter plot with VTK-style appearance
    scatter = ax.scatter(
        x_sorted, y_sorted,
        c=values_sorted,
        cmap=cmap,
        norm=norm,
        s=3,
        alpha=alpha_base,
        edgecolors='none',
        rasterized=True
    )
    
    # Add black hole marker with glow effect
    for r, alpha in [(0.8, 0.1), (0.5, 0.2), (0.3, 0.3)]:
        circle = Circle(bh_pos, r, color='white', alpha=alpha, zorder=10)
        ax.add_patch(circle)
    
    # Central black hole
    ax.plot(bh_pos[0], bh_pos[1], 'o', color='black', markersize=12, 
            markeredgecolor='white', markeredgewidth=2, zorder=11)
    ax.annotate('IMBH\n10⁵ M☉', bh_pos, 
                textcoords='offset points', xytext=(15, 15),
                fontsize=9, color='white', weight='bold',
                ha='left', va='bottom')
    
    # Gravitational influence radius (Hill sphere approximation)
    r_hill = 2.0  # pc, approximate
    influence_circle = Circle(bh_pos, r_hill, fill=False, 
                              color='cyan', linestyle='--', linewidth=1, alpha=0.5)
    ax.add_patch(influence_circle)
    
    # Styling
    ax.set_xlim(-25, 10)
    ax.set_ylim(-8, 8)
    ax.set_aspect('equal')
    ax.set_xlabel('X [pc]', fontsize=12, color='white')
    ax.set_ylabel('Y [pc]', fontsize=12, color='white')
    
    # Grid with subtle styling
    ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5)
    ax.tick_params(colors='white', labelsize=10)
    
    # Title with time
    ax.set_title(f'{title}\nt = {time:.3f} Myr', fontsize=14, color='white', weight='bold')
    
    # Statistics annotation
    stats_text = (
        f'N = {len(df):,}\n'
        f'{field.title()}: {np.median(values):.1f} (med)\n'
        f'Range: [{values.min():.1f}, {values.max():.1f}]'
    )
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=9, color='white', alpha=0.8,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))
    
    return scatter


def create_vtk_animation(results_dir: str, output_gif: str, title: str, field: str = 'velocity'):
    """Create VTK-style animation."""
    snapshot_files = get_snapshot_files(results_dir)
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Sample frames for faster GIF (every 2nd frame)
    step = 2
    snapshot_files = snapshot_files[::step]
    print(f"Using {len(snapshot_files)} frames (every {step})")
    
    # Load first frame for setup
    df0 = load_snapshot(snapshot_files[0])
    
    # Create figure with dark background
    fig, ax = plt.subplots(figsize=(12, 8), facecolor='black')
    ax.set_facecolor('black')
    
    # Create initial plot
    scatter = create_vtk_frame(ax, df0, field=field, time=0, title=title)
    
    # Add colorbar
    if field == 'velocity':
        cmap = 'plasma'
        vmin, vmax = 0, 150
        cbar_label = 'Velocity [km/s]'
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    elif field == 'temperature':
        cmap = 'inferno'
        vmin, vmax = 5, 100
        cbar_label = 'Temperature [K]'
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    elif field == 'density':
        cmap = 'viridis'
        vmin, vmax = 1, 1000
        cbar_label = 'Number Density [cm⁻³]'
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(cbar_label, fontsize=12, color='white')
    cbar.ax.tick_params(colors='white')
    
    plt.tight_layout()
    
    def update(frame_idx):
        if frame_idx % 10 == 0:
            print(f"  Processing frame {frame_idx + 1}/{len(snapshot_files)}")
        
        df = load_snapshot(snapshot_files[frame_idx])
        time = frame_idx * step * 0.01  # Approximate time based on snapshot interval
        create_vtk_frame(ax, df, field=field, time=time, title=title)
        return []
    
    print("Creating animation...")
    anim = animation.FuncAnimation(
        fig, update,
        frames=len(snapshot_files),
        interval=100,  # 100ms per frame
        blit=False
    )
    
    print(f"Saving to {output_gif}...")
    anim.save(output_gif, writer='pillow', fps=10, dpi=100)
    
    plt.close()
    print(f"✓ Saved {output_gif}")
    
    # Get file size
    size_kb = os.path.getsize(output_gif) / 1024
    print(f"  File size: {size_kb:.0f} KB")


def create_comparison_animation(adiabatic_dir: str, cooling_dir: str, output_gif: str, field: str = 'velocity'):
    """Create side-by-side comparison animation."""
    adiabatic_files = get_snapshot_files(adiabatic_dir)
    cooling_files = get_snapshot_files(cooling_dir)
    
    if not adiabatic_files or not cooling_files:
        print("Missing snapshot files")
        return
    
    # Match frame counts
    n_frames = min(len(adiabatic_files), len(cooling_files))
    print(f"Found {n_frames} matching frames")
    
    # Sample frames
    step = 2
    indices = list(range(0, n_frames, step))
    print(f"Using {len(indices)} frames")
    
    # Field settings
    if field == 'velocity':
        cmap = 'plasma'
        vmin, vmax = 0, 150
        cbar_label = 'Velocity [km/s]'
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    elif field == 'temperature':
        cmap = 'inferno'
        vmin, vmax = 5, 100
        cbar_label = 'Temperature [K]'
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    elif field == 'density':
        cmap = 'viridis'
        vmin, vmax = 1, 1000
        cbar_label = 'Number Density [cm⁻³]'
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(20, 8), facecolor='black')
    
    for ax in axes:
        ax.set_facecolor('black')
    
    # Add shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=axes, shrink=0.8, pad=0.02, location='right')
    cbar.set_label(cbar_label, fontsize=14, color='white')
    cbar.ax.tick_params(colors='white')
    
    def update(i):
        frame_idx = indices[i]
        if i % 10 == 0:
            print(f"  Processing frame {i + 1}/{len(indices)}")
        
        time = frame_idx * 0.01
        
        # Adiabatic
        df_adi = load_snapshot(adiabatic_files[frame_idx])
        create_vtk_frame(axes[0], df_adi, field=field, time=time, 
                        title='Phase 3b: Adiabatic (No Cooling)')
        
        # Cooling
        df_cool = load_snapshot(cooling_files[frame_idx])
        create_vtk_frame(axes[1], df_cool, field=field, time=time,
                        title='Phase 3.5: K&I 2000 Cooling')
        
        return []
    
    plt.tight_layout()
    
    print("Creating comparison animation...")
    anim = animation.FuncAnimation(
        fig, update,
        frames=len(indices),
        interval=100,
        blit=False
    )
    
    print(f"Saving to {output_gif}...")
    anim.save(output_gif, writer='pillow', fps=10, dpi=100)
    
    plt.close()
    print(f"✓ Saved {output_gif}")
    
    size_kb = os.path.getsize(output_gif) / 1024
    print(f"  File size: {size_kb:.0f} KB")


def main():
    parser = argparse.ArgumentParser(description='VTK-Style Animation Generator')
    parser.add_argument('results_dir', help='Results directory with snapshots')
    parser.add_argument('output_gif', help='Output GIF filename')
    parser.add_argument('title', help='Plot title')
    parser.add_argument('--field', choices=['velocity', 'temperature', 'density'],
                       default='velocity', help='Field to visualize')
    parser.add_argument('--compare', help='Second results directory for comparison')
    
    args = parser.parse_args()
    
    if args.compare:
        create_comparison_animation(args.results_dir, args.compare, args.output_gif, args.field)
    else:
        create_vtk_animation(args.results_dir, args.output_gif, args.title, args.field)


if __name__ == '__main__':
    main()
