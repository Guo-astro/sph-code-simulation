#!/usr/bin/env python3
"""
Create GIF animation showing radial profiles of density and temperature
for Phase 2.5 ISM Cooling simulation.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import imageio.v2 as imageio
import sys

# Physical constants for temperature calculation
k_B = 1.38e-16      # Boltzmann constant [erg/K]
m_p = 1.67e-24      # Proton mass [g]
mu = 1.27           # Mean molecular weight
gamma = 5/3         # Adiabatic index
u_to_cgs = 1e10     # erg/g per code energy
density_to_n_H = 31.86  # cm^-3 per code density

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

def compute_particle_data(df):
    """Compute radial distance and derived quantities for each particle."""
    real = df[df['is_ghost'] == 0].copy()
    
    # Compute radius (2D for cylindrical cloud)
    real['r'] = np.sqrt(real['pos_x']**2 + real['pos_y']**2)
    real['T'] = temperature_from_pressure_density(real['pres'], real['dens'])
    real['n_H'] = real['dens'] * density_to_n_H
    
    return real

def create_frame(df, output_path, time_val, frame_num):
    """Create a single frame showing radial profiles with raw particle data."""
    real = compute_particle_data(df)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left panel: Radial density profile (raw particles)
    ax1 = axes[0]
    ax1.scatter(real['r'], real['dens'], c='blue', s=0.5, alpha=0.3)
    ax1.set_yscale('log')
    ax1.set_xlabel('Radius [pc]', fontsize=12)
    ax1.set_ylabel('Density [code units]', fontsize=12)
    ax1.set_xlim(0, 0.8)
    ax1.set_ylim(1, 200)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Radial Density Profile (raw particles)', fontsize=14)
    
    # Right panel: Radial temperature profile (raw particles)
    ax2 = axes[1]
    ax2.scatter(real['r'], real['T'], c='red', s=0.5, alpha=0.3)
    ax2.set_yscale('log')
    ax2.axhline(y=6000, color='orange', ls='--', lw=1.5, label='WNM eq (~6000 K)')
    ax2.axhline(y=100, color='cyan', ls='--', lw=1.5, label='CNM eq (~100 K)')
    ax2.set_xlabel('Radius [pc]', fontsize=12)
    ax2.set_ylabel('Temperature [K]', fontsize=12)
    ax2.set_xlim(0, 0.8)
    ax2.set_ylim(1, 10000)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    ax2.set_title('Radial Temperature Profile (raw particles)', fontsize=14)
    
    # Summary stats
    T_mean = real['T'].mean()
    T_min, T_max = real['T'].min(), real['T'].max()
    dens_mean = real['dens'].mean()
    
    time_str = f't = {time_val:.3f}' if time_val is not None else f'Frame {frame_num}'
    fig.suptitle(f'Phase 2.5: ISM Cooling Radial Profiles - {time_str} code time\n'
                 f'T: {T_min:.1f}-{T_max:.1f} K (mean: {T_mean:.1f} K)', fontsize=13)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()

def main():
    results_dir = Path('/Users/guo-opt-p148/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/phase2.5_cooling_gamma53')
    
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    
    if len(snapshots) < 2:
        print(f"Need at least 2 snapshots for GIF, found {len(snapshots)}")
        sys.exit(1)
    
    print(f"Found {len(snapshots)} snapshots")
    
    frames_dir = results_dir / 'radial_frames'
    frames_dir.mkdir(exist_ok=True)
    
    frame_paths = []
    for i, snap_path in enumerate(snapshots):
        print(f"Processing {snap_path.name} ({i+1}/{len(snapshots)})")
        
        df = load_snapshot(snap_path)
        time_val = extract_time(snap_path)
        
        frame_path = frames_dir / f'radial_frame_{i:04d}.png'
        create_frame(df, frame_path, time_val, i)
        frame_paths.append(frame_path)
    
    print("Creating GIF...")
    gif_path = results_dir / 'phase2.5_radial_profiles.gif'
    
    images = []
    for fp in frame_paths:
        images.append(imageio.imread(fp))
    
    imageio.mimsave(gif_path, images, duration=0.5, loop=0)
    
    gif_size = gif_path.stat().st_size / (1024 * 1024)
    print(f"\n✓ GIF created: {gif_path}")
    print(f"  Size: {gif_size:.1f} MB")
    print(f"  Frames: {len(images)}")

if __name__ == '__main__':
    main()
