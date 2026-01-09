#!/usr/bin/env python3
"""
Plot density profile for the BE sphere relaxation simulation.

This script reads the final snapshot from the phase1_relaxation simulation
and plots the density as a function of radius.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

def read_csv_with_comments(filepath):
    """Read CSV file skipping comment lines starting with #"""
    with open(filepath, 'r') as f:
        # Skip all comment lines
        lines = f.readlines()
        
    # Find the header line (last comment line with column names)
    data_start = 0
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            data_start = i
            break
    
    # Read the actual data
    df = pd.read_csv(filepath, comment='#')
    return df

def main():
    # Path to results
    results_dir = Path(__file__).parent.parent / "results" / "phase1_relaxation"
    
    # Find the final snapshot
    snapshots = sorted(results_dir.glob("snapshot_*.csv"))
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        sys.exit(1)
    
    final_snapshot = snapshots[-1]
    print(f"Reading: {final_snapshot}")
    
    # Read data
    df = read_csv_with_comments(final_snapshot)
    
    # Filter out ghost particles
    df_real = df[df['is_ghost'] == 0].copy()
    print(f"Total particles: {len(df)}, Real particles: {len(df_real)}")
    
    # Calculate radius
    df_real['r'] = np.sqrt(df_real['pos_x']**2 + df_real['pos_y']**2 + df_real['pos_z']**2)
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Density vs radius scatter plot
    ax1 = axes[0, 0]
    scatter = ax1.scatter(df_real['r'], df_real['dens'], 
                          c=df_real['neighbor'], cmap='viridis', 
                          alpha=0.5, s=10)
    ax1.set_xlabel('Radius [code units]', fontsize=12)
    ax1.set_ylabel('Density [code units]', fontsize=12)
    ax1.set_title('Density Profile (colored by neighbor count)', fontsize=14)
    ax1.set_yscale('log')
    plt.colorbar(scatter, ax=ax1, label='Neighbors')
    ax1.grid(True, alpha=0.3)
    
    # 2. Density histogram/binned profile
    ax2 = axes[0, 1]
    # Bin the data by radius
    n_bins = 50
    r_bins = np.linspace(0, df_real['r'].max(), n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    density_mean = []
    density_std = []
    for i in range(n_bins):
        mask = (df_real['r'] >= r_bins[i]) & (df_real['r'] < r_bins[i+1])
        if mask.sum() > 0:
            density_mean.append(df_real.loc[mask, 'dens'].mean())
            density_std.append(df_real.loc[mask, 'dens'].std())
        else:
            density_mean.append(np.nan)
            density_std.append(np.nan)
    
    density_mean = np.array(density_mean)
    density_std = np.array(density_std)
    
    ax2.plot(r_centers, density_mean, 'b-', linewidth=2, label='Mean density')
    ax2.fill_between(r_centers, 
                     density_mean - density_std, 
                     density_mean + density_std,
                     alpha=0.3, label='±1σ')
    ax2.set_xlabel('Radius [code units]', fontsize=12)
    ax2.set_ylabel('Density [code units]', fontsize=12)
    ax2.set_title('Radially Averaged Density Profile', fontsize=14)
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. 2D slice (XY plane, near z=0)
    ax3 = axes[1, 0]
    z_slice_width = df_real['pos_z'].std() * 0.3
    slice_mask = np.abs(df_real['pos_z']) < z_slice_width
    df_slice = df_real[slice_mask]
    
    scatter3 = ax3.scatter(df_slice['pos_x'], df_slice['pos_y'], 
                           c=np.log10(df_slice['dens']), cmap='hot', 
                           alpha=0.7, s=15)
    ax3.set_xlabel('X [code units]', fontsize=12)
    ax3.set_ylabel('Y [code units]', fontsize=12)
    ax3.set_title(f'XY Slice (|z| < {z_slice_width:.3f})', fontsize=14)
    ax3.set_aspect('equal')
    plt.colorbar(scatter3, ax=ax3, label='log₁₀(Density)')
    
    # 4. Density histogram
    ax4 = axes[1, 1]
    ax4.hist(df_real['dens'], bins=50, edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Density [code units]', fontsize=12)
    ax4.set_ylabel('Count', fontsize=12)
    ax4.set_title('Density Distribution', fontsize=14)
    ax4.axvline(df_real['dens'].mean(), color='r', linestyle='--', 
                label=f'Mean: {df_real["dens"].mean():.2f}')
    ax4.axvline(df_real['dens'].median(), color='g', linestyle='--',
                label=f'Median: {df_real["dens"].median():.2f}')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Add overall title with statistics
    stats_text = (f"BE Sphere Relaxation - Final State\n"
                  f"N = {len(df_real)}, "
                  f"ρ_center ≈ {df_real[df_real['r'] < 0.05]['dens'].mean():.1f}, "
                  f"ρ_edge ≈ {df_real[df_real['r'] > df_real['r'].quantile(0.9)]['dens'].mean():.1f}")
    fig.suptitle(stats_text, fontsize=14, y=1.02)
    
    plt.tight_layout()
    
    # Save the figure
    output_path = results_dir / "density_profile.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot to: {output_path}")
    
    # Also display some statistics
    print("\n=== Density Statistics ===")
    print(f"Min density:  {df_real['dens'].min():.4f}")
    print(f"Max density:  {df_real['dens'].max():.4f}")
    print(f"Mean density: {df_real['dens'].mean():.4f}")
    print(f"Median density: {df_real['dens'].median():.4f}")
    print(f"Density contrast (max/min): {df_real['dens'].max() / df_real['dens'].min():.2f}")
    
    # Central vs edge density
    r_center = 0.05
    r_edge_min = df_real['r'].quantile(0.85)
    rho_center = df_real[df_real['r'] < r_center]['dens'].mean()
    rho_edge = df_real[df_real['r'] > r_edge_min]['dens'].mean()
    print(f"\nCentral density (r < {r_center}): {rho_center:.4f}")
    print(f"Edge density (r > {r_edge_min:.3f}): {rho_edge:.4f}")
    print(f"Center/Edge ratio: {rho_center/rho_edge:.2f}")
    
    plt.show()

if __name__ == "__main__":
    main()
