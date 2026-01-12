#!/usr/bin/env python3
"""
Create GIF animation for Phase 2.5 ISM Cooling simulation.
Shows density evolution with temperature color coding.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import imageio
import sys

# Physical constants for temperature calculation
k_B = 1.38e-16      # Boltzmann constant [erg/K]
m_p = 1.67e-24      # Proton mass [g]
mu = 1.27           # Mean molecular weight
gamma = 5/3         # Adiabatic index
u_to_cgs = 1e10     # erg/g per code energy

def temperature_from_pressure_density(pres, dens):
    """Convert pressure and density to temperature [K]."""
    u = pres / (dens * (gamma - 1))
    u_cgs = u * u_to_cgs
    T = u_cgs * mu * m_p * (gamma - 1) / k_B
    return T

def load_snapshot(filepath):
    """Load snapshot CSV, skipping comment lines."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('id,'):
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError(f"Could not find header in {filepath}")
    
    df = pd.read_csv(filepath, skiprows=header_idx)
    return df

def extract_time(filepath):
    """Extract simulation time from snapshot header."""
    with open(filepath, 'r') as f:
        for line in f:
            if 'time' in line.lower() and '=' in line:
                try:
                    return float(line.split('=')[1].strip())
                except:
                    pass
            if line.startswith('id,'):
                break
    return None

def create_frame(df, output_path, time_val, frame_num, total_frames):
    """Create a single frame showing x-y projection colored by temperature."""
    # Filter real particles
    real = df[df['is_ghost'] == 0].copy()
    
    # Calculate temperature
    real['T'] = temperature_from_pressure_density(real['pres'], real['dens'])
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: x-y projection colored by density
    ax1 = axes[0]
    scatter1 = ax1.scatter(real['pos_x'], real['pos_y'], 
                           c=real['dens'], cmap='viridis',
                           s=1, alpha=0.6, norm=LogNorm(vmin=1, vmax=200))
    ax1.set_xlabel('x [pc]')
    ax1.set_ylabel('y [pc]')
    ax1.set_title(f'Density (code units)')
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.set_aspect('equal')
    plt.colorbar(scatter1, ax=ax1, label='ρ [code]')
    
    # Right panel: x-y projection colored by temperature
    ax2 = axes[1]
    scatter2 = ax2.scatter(real['pos_x'], real['pos_y'], 
                           c=real['T'], cmap='hot',
                           s=1, alpha=0.6, norm=LogNorm(vmin=1, vmax=100))
    ax2.set_xlabel('x [pc]')
    ax2.set_ylabel('y [pc]')
    ax2.set_title(f'Temperature [K]')
    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-1, 1)
    ax2.set_aspect('equal')
    plt.colorbar(scatter2, ax=ax2, label='T [K]')
    
    # Summary stats
    T_mean = real['T'].mean()
    T_min, T_max = real['T'].min(), real['T'].max()
    dens_mean = real['dens'].mean()
    
    time_str = f't = {time_val:.3f}' if time_val is not None else f'Frame {frame_num}'
    fig.suptitle(f'Phase 2.5: ISM Cooling - {time_str} code time\n'
                 f'T: {T_min:.1f}-{T_max:.1f} K (mean: {T_mean:.1f} K), '
                 f'ρ_mean: {dens_mean:.1f}', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()
    
    return T_mean, dens_mean

def main():
    results_dir = Path('/Users/guo-opt-p148/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/phase2.5_cooling_gamma53')
    
    # Find all snapshots
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    
    if len(snapshots) < 2:
        print(f"Need at least 2 snapshots for GIF, found {len(snapshots)}")
        sys.exit(1)
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Create frames directory
    frames_dir = results_dir / 'frames'
    frames_dir.mkdir(exist_ok=True)
    
    # Generate frames
    frame_paths = []
    for i, snap_path in enumerate(snapshots):
        print(f"Processing {snap_path.name} ({i+1}/{len(snapshots)})")
        
        df = load_snapshot(snap_path)
        time_val = extract_time(snap_path)
        
        frame_path = frames_dir / f'frame_{i:04d}.png'
        create_frame(df, frame_path, time_val, i, len(snapshots))
        frame_paths.append(frame_path)
    
    # Create GIF
    print("Creating GIF...")
    gif_path = results_dir / 'phase2.5_cooling.gif'
    
    images = []
    for fp in frame_paths:
        images.append(imageio.imread(fp))
    
    # Duration per frame in seconds (slower for better viewing)
    imageio.mimsave(gif_path, images, duration=0.5, loop=0)
    
    # Get file size
    gif_size = gif_path.stat().st_size / (1024 * 1024)
    print(f"\n✓ GIF created: {gif_path}")
    print(f"  Size: {gif_size:.1f} MB")
    print(f"  Frames: {len(images)}")

if __name__ == '__main__':
    main()
