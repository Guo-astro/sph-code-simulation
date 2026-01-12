#!/usr/bin/env python3
"""
Compare Lane-Emden multi-method results with analytical solution.
Plots radial profiles of density, pressure, and velocity for all SPH methods.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from pathlib import Path

def load_analytical_solution(datafile):
    """Load Lane-Emden n=1.5 analytical solution."""
    if not os.path.exists(datafile):
        print(f"Warning: Analytical data file not found: {datafile}")
        return None
    
    data = np.loadtxt(datafile, skiprows=1)
    xi = data[:, 0]
    theta = data[:, 1]
    dtheta_dxi = data[:, 2]
    
    # For n=1.5, gamma=5/3
    gamma = 5.0/3.0
    n = 1.5
    
    # Normalization constants from Lane-Emden (M=1, R=1)
    xi_1 = xi[-1]  # First zero
    M = 1.0
    R = 1.0
    
    # Density profile: ρ(r) = ρ_c * θ^n
    # ρ_c = 1.4301 for M=1, R=1, n=1.5
    rho_c = 1.4301
    
    # Convert to physical radius
    r = xi * R / xi_1
    rho = rho_c * theta**n
    
    # Pressure from polytrope: P = K ρ^γ
    # K = 0.424209 for this setup
    K = 0.424209
    P = K * rho**gamma
    
    return r, rho, P

def load_snapshot(filepath):
    """Load CSV snapshot."""
    # Skip comment lines starting with #
    df = pd.read_csv(filepath, comment='#')
    
    # Calculate radius - handle different column naming conventions
    x_col = 'pos_x' if 'pos_x' in df.columns else 'x'
    y_col = 'pos_y' if 'pos_y' in df.columns else 'y'
    z_col = 'pos_z' if 'pos_z' in df.columns else 'z'
    
    r = np.sqrt(df[x_col]**2 + df[y_col]**2 + df[z_col]**2)
    
    r = np.sqrt(df[x_col]**2 + df[y_col]**2 + df[z_col]**2)
    
    # Handle different column naming conventions
    rho_col = 'dens' if 'dens' in df.columns else ('rho' if 'rho' in df.columns else 'density')
    P_col = 'pres' if 'pres' in df.columns else ('P' if 'P' in df.columns else 'pressure')
    vx_col = 'vel_x' if 'vel_x' in df.columns else 'vx'
    vy_col = 'vel_y' if 'vel_y' in df.columns else 'vy'
    vz_col = 'vel_z' if 'vel_z' in df.columns else 'vz'
    
    return {
        'r': r.values,
        'rho': df[rho_col].values,
        'P': df[P_col].values,
        'vx': df[vx_col].values,
        'vy': df[vy_col].values,
        'vz': df[vz_col].values,
    }

def plot_comparison(comparison_dir, output_dir, snapshot_num='_0064'):
    """Generate comparison plots."""
    
    methods = ['ssph', 'gsph', 'disph', 'gdisph']
    colors = {'ssph': 'blue', 'gsph': 'green', 'disph': 'red', 'gdisph': 'purple'}
    labels = {'ssph': 'SSPH', 'gsph': 'GSPH', 'disph': 'DISPH', 'gdisph': 'GDISPH'}
    
    # Load analytical solution
    analytical_file = 'lane_emden/data/numerical_solutions/3d/n1.5.dat'
    if not os.path.exists(analytical_file):
        analytical_file = 'data/lane_emden/n1.5_3d.dat'
    
    r_ana, rho_ana, P_ana = load_analytical_solution(analytical_file) if os.path.exists(analytical_file) else (None, None, None)
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Load data for each method
    for method in methods:
        snapshot_file = f"{comparison_dir}/{method}/snapshot{snapshot_num}.csv"
        if not os.path.exists(snapshot_file):
            print(f"Warning: Snapshot not found: {snapshot_file}")
            continue
        
        data = load_snapshot(snapshot_file)
        
        # Sort by radius for plotting
        idx = np.argsort(data['r'])
        r = data['r'][idx]
        rho = data['rho'][idx]
        P = data['P'][idx]
        v = np.sqrt(data['vx']**2 + data['vy']**2 + data['vz']**2)[idx]
        
        # Bin data for clearer plotting
        r_bins = np.linspace(0, 1.2, 50)
        rho_binned = np.zeros(len(r_bins)-1)
        P_binned = np.zeros(len(r_bins)-1)
        v_binned = np.zeros(len(r_bins)-1)
        
        for i in range(len(r_bins)-1):
            mask = (r >= r_bins[i]) & (r < r_bins[i+1])
            if np.sum(mask) > 0:
                rho_binned[i] = np.mean(rho[mask])
                P_binned[i] = np.mean(P[mask])
                v_binned[i] = np.mean(v[mask])
        
        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        
        # Plot density
        axes[0].plot(r_centers, rho_binned, 'o-', color=colors[method], 
                    label=labels[method], alpha=0.7, markersize=4)
        
        # Plot pressure
        axes[1].plot(r_centers, P_binned, 'o-', color=colors[method], 
                    label=labels[method], alpha=0.7, markersize=4)
        
        # Plot velocity
        axes[2].plot(r_centers, v_binned, 'o-', color=colors[method], 
                    label=labels[method], alpha=0.7, markersize=4)
    
    # Overlay analytical solution
    if r_ana is not None:
        axes[0].plot(r_ana, rho_ana, 'k--', linewidth=2, label='Analytical', zorder=10)
        axes[1].plot(r_ana, P_ana, 'k--', linewidth=2, label='Analytical', zorder=10)
    
    # Format plots
    axes[0].set_xlabel('Radius', fontsize=12)
    axes[0].set_ylabel('Density', fontsize=12)
    axes[0].set_title('Density Profile', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Radius', fontsize=12)
    axes[1].set_ylabel('Pressure', fontsize=12)
    axes[1].set_title('Pressure Profile', fontsize=14, fontweight='bold')
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_xlabel('Radius', fontsize=12)
    axes[2].set_ylabel('Velocity', fontsize=12)
    axes[2].set_title('Velocity Profile', fontsize=14, fontweight='bold')
    axes[2].legend(fontsize=10)
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = f"{output_dir}/lane_emden_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    
    plt.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: compare_lane_emden.py <comparison_dir> <output_dir>")
        sys.exit(1)
    
    comparison_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("Generating Lane-Emden comparison plots...")
    plot_comparison(comparison_dir, output_dir)
    print("✓ All comparison plots generated successfully!")
