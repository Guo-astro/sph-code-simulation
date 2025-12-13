#!/usr/bin/env python3
"""
Compare GSPH/SSPH × grad-h/no-grad-h methods for 3D cylindrical Lane-Emden.

This script generates:
1. Central density evolution comparison
2. Radial density profile comparison at different times
3. Energy conservation analysis
4. Cross-section profiles (xy-plane)

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
import os
import glob

from scripts.shared.lane_emden import solve_lane_emden_cylindrical as _solve_lane_emden_cylindrical

# Configuration
BASE_DIR = "simulations/stability/gradh_study_3d_cylinder/results"
OUTPUT_DIR = "simulations/stability/gradh_study_3d_cylinder/results/comparison"

# Methods to compare (using high-contrast colors)
METHODS = {
    'gsph_gradh': {'label': 'GSPH + grad-h', 'color': '#0066CC', 'linestyle': '-', 'marker': 'o'},
    'gsph_nogradh': {'label': 'GSPH - grad-h', 'color': '#9933FF', 'linestyle': '--', 'marker': 's'},
    'ssph_gradh': {'label': 'SSPH + grad-h', 'color': '#CC0000', 'linestyle': '-', 'marker': '^'},
    'ssph_nogradh': {'label': 'SSPH - grad-h', 'color': '#FF6600', 'linestyle': '--', 'marker': 'v'},
}

# Physical parameters
RHO_CENTER = 1.0
K = 1.0
G = 1.0
GAMMA = 5.0/3.0  # γ = 5/3 → n = 1.5 polytrope


def get_analytical_profile_cylinder(r_vals, rho_c, K, G, gamma):
    """Get analytical cylindrical Lane-Emden density profile using SSOT."""
    n = 1.0 / (gamma - 1.0)
    
    # Length scale for cylindrical case: α² = K(n+1)ρ_c^(1-n) / (4G)
    alpha_sq = K * (n + 1.0) * rho_c**(1.0 - n) / (4.0 * G)
    alpha = np.sqrt(alpha_sq)
    
    # Solve Lane-Emden using SSOT
    xi_le, theta_le = _solve_lane_emden_cylindrical(n, xi_max=10.0, n_points=10000)
    theta_le = np.maximum(theta_le, 0)  # θ ≥ 0
    
    # Find surface (where theta = 0)
    surface_idx = np.argmax(theta_le <= 0) if np.any(theta_le <= 0) else len(theta_le) - 1
    xi_surface = xi_le[surface_idx]
    
    # Interpolate to get density at given r values
    rho_vals = np.zeros_like(r_vals)
    for i, r in enumerate(r_vals):
        xi = abs(r) / alpha
        if xi <= xi_surface:
            theta = np.interp(xi, xi_le, theta_le)
            rho_vals[i] = rho_c * theta**n
        else:
            rho_vals[i] = 0.0
    
    return rho_vals, alpha * xi_surface


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


def analyze_snapshot_3d_cylinder(snap):
    """Extract key metrics from a 3D cylindrical snapshot."""
    if snap is None:
        return None
    
    # Compute cylindrical radius (xy-plane only)
    r_cyl = np.sqrt(snap['pos_x']**2 + snap['pos_y']**2)
    
    # Central density (r < 0.5)
    mask_center = r_cyl < 0.5
    if mask_center.sum() > 0:
        rho_center = snap[mask_center]['dens'].mean()
    else:
        rho_center = snap['dens'].max()
    
    # Extent in radial direction
    r_max = r_cyl.max()
    
    # Max velocity
    vel_xy = np.sqrt(snap['vel_x']**2 + snap['vel_y']**2)
    vel_mag = np.sqrt(snap['vel_x']**2 + snap['vel_y']**2 + snap['vel_z']**2)
    max_vel = vel_mag.max()
    max_vel_radial = vel_xy.max()
    
    # Max density
    max_dens = snap['dens'].max()
    
    return {
        'rho_center': rho_center,
        'max_dens': max_dens,
        'r_max': r_max,
        'max_vel': max_vel,
        'max_vel_radial': max_vel_radial
    }


def read_energy_file(directory):
    """Read energy.dat file."""
    filepath = os.path.join(directory, "energy.dat")
    if not os.path.exists(filepath):
        return None
    
    data = np.loadtxt(filepath, skiprows=1)
    return {
        'time': data[:, 0],
        'kinetic': data[:, 1],
        'thermal': data[:, 2],
        'potential': data[:, 3],
        'total': data[:, 4]
    }


def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 70)
    print("4-METHOD COMPARISON: GSPH/SSPH × grad-h/no-grad-h")
    print("3D Cylindrical Lane-Emden (γ=5/3, n=1.5)")
    print("with Kernel-Convolved 2D Gravity (xy-plane only)")
    print("=" * 70)
    
    # Get analytical solution
    r_analytical = np.linspace(0, 10, 1000)
    rho_analytical, r_surface = get_analytical_profile_cylinder(r_analytical, RHO_CENTER, K, G, GAMMA)
    print(f"Analytical solution: r_surface = {r_surface:.4f}")
    print(f"Polytropic index n = {1.0/(GAMMA-1.0):.2f}, γ = {GAMMA:.4f}")
    
    # Load all snapshots for each method
    method_data = {}
    num_snapshots = None
    
    for method_name in METHODS:
        method_dir = os.path.join(BASE_DIR, method_name)
        if not os.path.exists(method_dir):
            print(f"  {method_name}: NOT FOUND - skipping")
            continue
            
        n_snaps = count_snapshots(method_dir)
        if n_snaps == 0:
            print(f"  {method_name}: no snapshots - skipping")
            continue
            
        if num_snapshots is None:
            num_snapshots = n_snaps
        print(f"  {method_name}: {n_snaps} snapshots")
        
        snaps = []
        metrics = []
        for i in range(n_snaps):
            snap = read_snapshot(method_dir, i)
            snaps.append(snap)
            metrics.append(analyze_snapshot_3d_cylinder(snap))
        
        energy = read_energy_file(method_dir)
        
        method_data[method_name] = {
            'snaps': snaps,
            'metrics': metrics,
            'energy': energy
        }
    
    if not method_data:
        print("ERROR: No simulation data found!")
        print("Run simulations first with: make gradh_3d_compare_run")
        return
    
    # Estimate output interval from config (default 0.5)
    output_interval = 0.5
    times = np.arange(num_snapshots) * output_interval
    
    print(f"\nTime range: 0 to {times[-1]:.1f}")
    
    # =========================================================================
    # Plot 1: Central density evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        if method_name not in method_data:
            continue
        metrics = method_data[method_name]['metrics']
        rho_center = [m['rho_center'] for m in metrics if m is not None]
        t = times[:len(rho_center)]
        ax.plot(t, rho_center, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=max(1, len(t)//10),
                label=props['label'], linewidth=2, markersize=6)
    
    ax.axhline(y=RHO_CENTER, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label='Analytical ρ_c = 1.0')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Central Density ρ_c', fontsize=12)
    ax.set_title('Central Density Evolution: 3D Cylinder\n(γ=5/3 Polytrope, Kernel-Convolved 2D Gravity in xy)', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'central_density_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/central_density_evolution.png")
    
    # =========================================================================
    # Plot 2: Maximum density evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        if method_name not in method_data:
            continue
        metrics = method_data[method_name]['metrics']
        max_dens = [m['max_dens'] for m in metrics if m is not None]
        t = times[:len(max_dens)]
        ax.plot(t, max_dens, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=max(1, len(t)//10),
                label=props['label'], linewidth=2, markersize=6)
    
    ax.axhline(y=RHO_CENTER, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label='Analytical ρ_c = 1.0')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Maximum Density ρ_max', fontsize=12)
    ax.set_title('Maximum Density Evolution: 3D Cylinder', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'max_density_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/max_density_evolution.png")
    
    # =========================================================================
    # Plot 3: Cylinder radius evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        if method_name not in method_data:
            continue
        metrics = method_data[method_name]['metrics']
        r_max = [m['r_max'] for m in metrics if m is not None]
        t = times[:len(r_max)]
        ax.plot(t, r_max, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=max(1, len(t)//10),
                label=props['label'], linewidth=2, markersize=6)
    
    ax.axhline(y=r_surface, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label=f'Analytical r_surface = {r_surface:.3f}')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Cylinder Radius (xy-plane)', fontsize=12)
    ax.set_title('Cylinder Radius Evolution: 3D Lane-Emden Cylinder', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'radius_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/radius_evolution.png")
    
    # =========================================================================
    # Plot 4: Max velocity evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        if method_name not in method_data:
            continue
        metrics = method_data[method_name]['metrics']
        max_vel = [m['max_vel'] for m in metrics if m is not None]
        t = times[:len(max_vel)]
        ax.plot(t, max_vel, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=max(1, len(t)//10),
                label=props['label'], linewidth=2, markersize=6)
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Maximum |v|', fontsize=12)
    ax.set_title('Maximum Velocity Evolution: 3D Cylinder\n(Indicator of Stability)', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'velocity_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/velocity_evolution.png")
    
    # =========================================================================
    # Plot 5: Energy conservation
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    for method_name, props in METHODS.items():
        if method_name not in method_data:
            continue
        energy = method_data[method_name]['energy']
        if energy is None:
            continue
        
        t = energy['time']
        total_E0 = energy['total'][0]
        
        # Total energy relative error
        ax = axes[0, 0]
        rel_error = np.abs(energy['total'] - total_E0) / np.abs(total_E0) * 100
        ax.plot(t, rel_error, color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.5)
        
        # Kinetic energy
        ax = axes[0, 1]
        ax.plot(t, energy['kinetic'], color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.5)
        
        # Thermal energy
        ax = axes[1, 0]
        ax.plot(t, energy['thermal'], color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.5)
        
        # Potential energy
        ax = axes[1, 1]
        ax.plot(t, energy['potential'], color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.5)
    
    axes[0, 0].set_xlabel('Time')
    axes[0, 0].set_ylabel('|ΔE/E₀| (%)')
    axes[0, 0].set_title('Total Energy Conservation')
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_yscale('log')
    
    axes[0, 1].set_xlabel('Time')
    axes[0, 1].set_ylabel('Kinetic Energy')
    axes[0, 1].set_title('Kinetic Energy Evolution')
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('Time')
    axes[1, 0].set_ylabel('Thermal Energy')
    axes[1, 0].set_title('Thermal Energy Evolution')
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('Time')
    axes[1, 1].set_ylabel('Potential Energy')
    axes[1, 1].set_title('Potential Energy Evolution')
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.suptitle('Energy Conservation: 3D Cylindrical Lane-Emden', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'energy_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/energy_evolution.png")
    
    # =========================================================================
    # Plot 6: Radial density profile at final time (xy-plane cross-section)
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Analytical profile
    ax.plot(r_analytical, rho_analytical, 'k-', linewidth=2, alpha=0.7, label='Analytical')
    
    for method_name, props in METHODS.items():
        if method_name not in method_data:
            continue
        snaps = method_data[method_name]['snaps']
        if snaps[-1] is None:
            continue
        
        snap = snaps[-1]
        # Get cross section near z=0
        mask = np.abs(snap['pos_z']) < 0.5
        r_vals = np.sqrt(snap[mask]['pos_x']**2 + snap[mask]['pos_y']**2).values
        rho_vals = snap[mask]['dens'].values
        
        # Sort by r
        sort_idx = np.argsort(r_vals)
        ax.scatter(r_vals[sort_idx], rho_vals[sort_idx], 
                   color=props['color'], marker=props['marker'],
                   label=props['label'], alpha=0.6, s=20)
    
    ax.set_xlabel('Cylindrical radius r', fontsize=12)
    ax.set_ylabel('Density ρ', fontsize=12)
    ax.set_title('Radial Density Profile at Final Time (cross-section at z≈0)\n3D Cylindrical Lane-Emden', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, r_surface*1.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'density_profiles.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/density_profiles.png")
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
