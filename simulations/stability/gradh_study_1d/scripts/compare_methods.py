#!/usr/bin/env python3
"""
Compare GSPH/SSPH × grad-h/no-grad-h methods for 1D polytropic slab with kernel gravity.

This script generates:
1. Central density evolution comparison
2. Density profile comparison at different times
3. Energy conservation analysis
4. Animation comparing all 4 methods

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
import os
import glob

from scripts.shared.lane_emden import solve_lane_emden_planar

# Configuration
BASE_DIR = "simulations/stability/gradh_study_1d/results/kernel_gravity"
CSMOOTH_DIR = "simulations/stability/gradh_study_1d/results/csmooth_test"
OUTPUT_DIR = "simulations/stability/gradh_study_1d/results/comparison"
RELAXATION_DIR = "simulations/stability/gradh_study_1d/results/relaxation"

# Methods to compare (using high-contrast colors)
# Original 6 methods (C_smooth = 1, default)
METHODS = {
    'gsph_gradh': {'label': 'GSPH + grad-h', 'color': '#0066CC', 'linestyle': '-', 'marker': 'o'},      # Dark blue
    'gsph_nogradh': {'label': 'GSPH - grad-h', 'color': '#9933FF', 'linestyle': '--', 'marker': 's'},   # Purple
    'ssph_gradh': {'label': 'SSPH + grad-h', 'color': '#CC0000', 'linestyle': '-', 'marker': '^'},      # Dark red
    'ssph_nogradh': {'label': 'SSPH - grad-h', 'color': '#FF6600', 'linestyle': '--', 'marker': 'v'},   # Orange
    'gdisph_gradh': {'label': 'GDISPH + grad-h', 'color': '#009933', 'linestyle': '-', 'marker': 'D'},  # Dark green
    'gdisph_nogradh': {'label': 'GDISPH - grad-h', 'color': '#00CCCC', 'linestyle': '--', 'marker': 'p'},  # Cyan
}

# C_smooth = 2 methods (smoother h variation) - GSPH only
CSMOOTH_METHODS = {
    'gsph_gradh_cs2': {'label': 'GSPH + grad-h (C=2)', 'color': '#FF00FF', 'linestyle': '-', 'marker': 'h'},     # Magenta
    'gsph_nogradh_cs2': {'label': 'GSPH - grad-h (C=2)', 'color': '#00FF00', 'linestyle': '--', 'marker': 'H'},  # Lime green
}

# Physical parameters (from polytropic slab setup)
# γ = 5/3 polytrope with polytropic index n = 1/(γ-1) = 1.5
# This is the standard test case for hydrostatic equilibrium
RHO_CENTER = 1.0
K = 1.0
G = 1.0
GAMMA = 5.0/3.0  # γ = 5/3 → n = 1.5 polytrope


def get_analytical_profile(x_vals, rho_c, K, G, gamma):
    """Get analytical Lane-Emden density profile using SSOT."""
    n = 1.0 / (gamma - 1.0)
    
    # Length scale: α² = K(n+1)ρ_c^(1-n) / (2πG)
    alpha_sq = K * (n + 1.0) * rho_c**(1.0 - n) / (2.0 * np.pi * G)
    alpha = np.sqrt(alpha_sq)
    
    # Solve Lane-Emden using SSOT
    xi_le, theta_le = solve_lane_emden_planar(n, xi_max=10.0, n_points=10000)
    theta_le = np.maximum(theta_le, 0)  # θ ≥ 0
    
    # Find surface (where theta = 0)
    surface_idx = np.argmax(theta_le <= 0) if np.any(theta_le <= 0) else len(theta_le) - 1
    xi_surface = xi_le[surface_idx]
    
    # Interpolate to get density at given x values
    rho_vals = np.zeros_like(x_vals)
    for i, x in enumerate(x_vals):
        xi = abs(x) / alpha
        if xi <= xi_surface:
            theta = np.interp(xi, xi_le, theta_le)
            rho_vals[i] = rho_c * theta**n
        else:
            rho_vals[i] = 0.0
    
    return rho_vals, alpha * xi_surface


def compute_analytical_energies(rho_c, K, G, gamma):
    """
    Compute analytical energy components for 1D polytropic slab.
    
    For a polytrope with EOS P = K*rho^gamma:
      - Thermal energy: E_th = 2 * ∫_0^x_surface P/(γ-1) dx
      - Gravitational energy: E_grav = (1/2) ∫ ρ*Φ dx
        where Φ(x) = -2πG ∫ |x - x'| ρ(x') dx' (1D Green's function)
    
    Using Lane-Emden solution: ρ = ρ_c * θ^n, x = α*ξ
    """
    n = 1.0 / (gamma - 1.0)  # n = 1.5 for γ = 5/3
    
    # Length scale: α² = K(n+1)ρ_c^(1/n-1) / (2πG)
    alpha_sq = K * (n + 1.0) * rho_c**(1.0/n - 1.0) / (2.0 * np.pi * G)
    alpha = np.sqrt(alpha_sq)
    
    # Solve Lane-Emden using SSOT
    xi_le, theta_le = solve_lane_emden_planar(n, xi_max=10.0, n_points=10000)
    theta_le = np.maximum(theta_le, 0)  # θ ≥ 0
    
    # Physical coordinates
    x_le = alpha * xi_le
    rho_le = rho_c * theta_le**n
    P_le = K * rho_le**gamma
    
    # Thermal energy: E_th = 2 * ∫_0^x_surface P/(γ-1) dx
    E_thermal = 2.0 * np.trapezoid(P_le / (gamma - 1.0), x_le)
    
    # Gravitational energy using 1D Green's function:
    # Φ(x) = -2πG ∫ |x - x'| ρ(x') dx'
    # E_grav = (1/2) ∫ ρ Φ dx
    #
    # Use subset of points for efficiency
    N = min(200, len(x_le))
    idx = np.linspace(0, len(x_le)-1, N, dtype=int)
    x_sub = x_le[idx]
    rho_sub = rho_le[idx]
    
    # Build full symmetric array (slab from -x_surface to +x_surface)
    x_full = np.concatenate([-x_sub[::-1][:-1], x_sub])
    rho_full = np.concatenate([rho_sub[::-1][:-1], rho_sub])
    
    # Integration weights
    dx = np.diff(x_full)
    dx = np.append(dx, dx[-1])
    
    # Compute |x_i - x_j| matrix
    X_i, X_j = np.meshgrid(x_full, x_full, indexing='ij')
    dist_matrix = np.abs(X_i - X_j)
    
    # Φ_i = -2πG ∫ |x_i - x_j| ρ_j dx_j
    phi_full = -2.0 * np.pi * G * np.sum(dist_matrix * rho_full * dx, axis=1)
    
    # E_grav = (1/2) ∫ ρ Φ dx
    E_gravitational = 0.5 * np.sum(rho_full * phi_full * dx)
    
    # Total energy (kinetic is zero for static equilibrium)
    E_kinetic = 0.0
    E_total = E_thermal + E_gravitational + E_kinetic
    
    return {
        'thermal': E_thermal,
        'gravitational': E_gravitational,
        'kinetic': E_kinetic,
        'total': E_total
    }


def calculate_momentum_from_snapshot(snap):
    """Calculate total momentum from a snapshot."""
    if snap is None:
        return None
    # Total momentum P = Σ m_i * v_i
    return (snap['mass'] * snap['vel_x']).sum()


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


def analyze_snapshot(snap):
    """Extract key metrics from a snapshot."""
    if snap is None:
        return None
    
    # Central density (|x| < 0.5)
    mask_center = np.abs(snap['pos_x']) < 0.5
    if mask_center.sum() > 0:
        rho_center = snap[mask_center]['dens'].mean()
    else:
        rho_center = snap['dens'].max()
    
    # Extent
    x_min = snap['pos_x'].min()
    x_max = snap['pos_x'].max()
    
    # Max velocity
    max_vel = np.abs(snap['vel_x']).max()
    
    # Max density
    max_dens = snap['dens'].max()
    
    return {
        'rho_center': rho_center,
        'max_dens': max_dens,
        'x_min': x_min,
        'x_max': x_max,
        'max_vel': max_vel
    }


def read_energy_file(directory):
    """Read energy.dat file."""
    filepath = os.path.join(directory, "energy.dat")
    if not os.path.exists(filepath):
        return None
    
    # Energy file format: time, kinetic, thermal, potential, total
    try:
        data = np.loadtxt(filepath, skiprows=1)
        if data.ndim < 2 or data.shape[0] == 0:
            return None
        return {
            'time': data[:, 0],
            'kinetic': data[:, 1],
            'thermal': data[:, 2],
            'potential': data[:, 3],
            'total': data[:, 4]
        }
    except Exception:
        return None


def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 70)
    print("8-METHOD COMPARISON: GSPH/SSPH/GDISPH × grad-h/no-grad-h")
    print("with Kernel-Convolved 1D Gravity (+ GSPH C_smooth=2 variants)")
    print("=" * 70)
    
    # Get analytical solution
    x_analytical = np.linspace(-8, 8, 1000)
    rho_analytical, x_surface = get_analytical_profile(x_analytical, RHO_CENTER, K, G, GAMMA)
    print(f"Analytical solution: x_surface = {x_surface:.4f}")
    print(f"Polytropic index n = {1.0/(GAMMA-1.0):.2f}, γ = {GAMMA:.4f}")
    
    # Compute analytical energies
    analytical_E = compute_analytical_energies(RHO_CENTER, K, G, GAMMA)
    print("\nAnalytical energies:")
    print(f"  E_thermal     = {analytical_E['thermal']:.6f}")
    print(f"  E_gravitational = {analytical_E['gravitational']:.6f}")
    print(f"  E_kinetic     = {analytical_E['kinetic']:.6f}")
    print(f"  E_total       = {analytical_E['total']:.6f}")
    
    # Load all snapshots for each method (C_smooth = 1)
    method_data = {}
    num_snapshots = None
    
    print("\n--- Loading C_smooth = 1 methods ---")
    for method_name in METHODS:
        method_dir = os.path.join(BASE_DIR, method_name)
        n_snaps = count_snapshots(method_dir)
        if num_snapshots is None:
            num_snapshots = n_snaps
        print(f"  {method_name}: {n_snaps} snapshots")
        
        snaps = []
        metrics = []
        momentum = []
        for i in range(n_snaps):
            snap = read_snapshot(method_dir, i)
            snaps.append(snap)
            metrics.append(analyze_snapshot(snap))
            momentum.append(calculate_momentum_from_snapshot(snap))
        
        energy = read_energy_file(method_dir)
        
        method_data[method_name] = {
            'snaps': snaps,
            'metrics': metrics,
            'energy': energy,
            'momentum': momentum,
            'c_smooth': 1.0
        }
    
    # Load C_smooth = 2 methods
    csmooth_data = {}
    num_csmooth_snapshots = None
    
    print("\n--- Loading C_smooth = 2 methods ---")
    for method_name in CSMOOTH_METHODS:
        method_dir = os.path.join(CSMOOTH_DIR, method_name)
        n_snaps = count_snapshots(method_dir)
        if n_snaps == 0:
            print(f"  {method_name}: No data found (skipping)")
            continue
        if num_csmooth_snapshots is None:
            num_csmooth_snapshots = n_snaps
        print(f"  {method_name}: {n_snaps} snapshots")
        
        snaps = []
        metrics = []
        momentum = []
        for i in range(n_snaps):
            snap = read_snapshot(method_dir, i)
            snaps.append(snap)
            metrics.append(analyze_snapshot(snap))
            momentum.append(calculate_momentum_from_snapshot(snap))
        
        energy = read_energy_file(method_dir)
        
        csmooth_data[method_name] = {
            'snaps': snaps,
            'metrics': metrics,
            'energy': energy,
            'momentum': momentum,
            'c_smooth': 2.0
        }
    
    # Get time from first snapshot header
    first_snap = method_data['gsph_gradh']['snaps'][-1]
    if first_snap is not None:
        # Output interval: 5.0 for t=0→500 (101 snapshots)
        output_interval = 5.0
        times = np.arange(num_snapshots) * output_interval
    
    # C_smooth times (may be different duration)
    if num_csmooth_snapshots:
        csmooth_times = np.arange(num_csmooth_snapshots) * output_interval
    else:
        csmooth_times = times
    
    print(f"\nTime range (C_smooth=1): 0 to {times[-1]:.1f}")
    if num_csmooth_snapshots:
        print(f"Time range (C_smooth=2): 0 to {csmooth_times[-1]:.1f}")
    
    # =========================================================================
    # Plot 1: Central density evolution (all methods including C_smooth=2)
    # =========================================================================
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot C_smooth = 1 methods
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        rho_center = [m['rho_center'] for m in metrics if m is not None]
        t = times[:len(rho_center)]
        ax.plot(t, rho_center, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    # Plot C_smooth = 2 methods
    for method_name, props in CSMOOTH_METHODS.items():
        if method_name not in csmooth_data:
            continue
        metrics = csmooth_data[method_name]['metrics']
        rho_center = [m['rho_center'] for m in metrics if m is not None]
        t = csmooth_times[:len(rho_center)]
        ax.plot(t, rho_center, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=5, alpha=0.7,
                label=props['label'], linewidth=1.5, markersize=5)
    
    ax.axhline(y=RHO_CENTER, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label='Analytical ρ_c = 1.0')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Central Density ρ_c', fontsize=12)
    ax.set_title('Central Density Evolution: Method Comparison\n(γ=5/3 Polytropic Slab, n=1.5, Kernel-Convolved 1D Gravity)', fontsize=14)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'central_density_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/central_density_evolution.png")
    
    # =========================================================================
    # Plot 2: Maximum density evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot C_smooth=1 methods
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        max_dens = [m['max_dens'] for m in metrics if m is not None]
        t = times[:len(max_dens)]
        ax.plot(t, max_dens, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    # Plot C_smooth=2 methods
    for method_name, props in CSMOOTH_METHODS.items():
        if method_name not in csmooth_data:
            continue
        metrics = csmooth_data[method_name]['metrics']
        max_dens = [m['max_dens'] for m in metrics if m is not None]
        t = csmooth_times[:len(max_dens)]
        ax.plot(t, max_dens, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=5, alpha=0.7,
                label=props['label'], linewidth=1.5, markersize=5)
    
    ax.axhline(y=RHO_CENTER, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label='Analytical ρ_c = 1.0')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Maximum Density ρ_max', fontsize=12)
    ax.set_title('Maximum Density Evolution: Method Comparison\n(Indicator of Core Collapse)', fontsize=14)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'max_density_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/max_density_evolution.png")
    
    # =========================================================================
    # Plot 3: Slab extent evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot C_smooth=1 methods
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        extent = [(m['x_max'] - m['x_min'])/2 for m in metrics if m is not None]
        t = times[:len(extent)]
        ax.plot(t, extent, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    # Plot C_smooth=2 methods
    for method_name, props in CSMOOTH_METHODS.items():
        if method_name not in csmooth_data:
            continue
        metrics = csmooth_data[method_name]['metrics']
        extent = [(m['x_max'] - m['x_min'])/2 for m in metrics if m is not None]
        t = csmooth_times[:len(extent)]
        ax.plot(t, extent, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=5, alpha=0.7,
                label=props['label'], linewidth=1.5, markersize=5)
    
    ax.axhline(y=x_surface, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label=f'Analytical x_surface = {x_surface:.3f}')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Slab Half-Width', fontsize=12)
    ax.set_title('Slab Extent Evolution: Method Comparison', fontsize=14)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'extent_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/extent_evolution.png")
    
    # =========================================================================
    # Plot 4: Max velocity evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot C_smooth=1 methods
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        max_vel = [m['max_vel'] for m in metrics if m is not None]
        t = times[:len(max_vel)]
        ax.plot(t, max_vel, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    # Plot C_smooth=2 methods
    for method_name, props in CSMOOTH_METHODS.items():
        if method_name not in csmooth_data:
            continue
        metrics = csmooth_data[method_name]['metrics']
        max_vel = [m['max_vel'] for m in metrics if m is not None]
        t = csmooth_times[:len(max_vel)]
        ax.plot(t, max_vel, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=5, alpha=0.7,
                label=props['label'], linewidth=1.5, markersize=5)
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Maximum |v_x|', fontsize=12)
    ax.set_title('Maximum Velocity Evolution: Method Comparison\n(Indicator of Stability)', fontsize=14)
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
    
    # Plot C_smooth=1 methods
    for method_name, props in METHODS.items():
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
    
    # Plot C_smooth=2 methods
    for method_name, props in CSMOOTH_METHODS.items():
        if method_name not in csmooth_data:
            continue
        energy = csmooth_data[method_name]['energy']
        if energy is None:
            continue
        
        t = energy['time']
        total_E0 = energy['total'][0]
        
        # Total energy relative error
        ax = axes[0, 0]
        rel_error = np.abs(energy['total'] - total_E0) / np.abs(total_E0) * 100
        ax.plot(t, rel_error, color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.2, alpha=0.7)
        
        # Kinetic energy
        ax = axes[0, 1]
        ax.plot(t, energy['kinetic'], color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.2, alpha=0.7)
        
        # Thermal energy
        ax = axes[1, 0]
        ax.plot(t, energy['thermal'], color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.2, alpha=0.7)
        
        # Potential energy
        ax = axes[1, 1]
        ax.plot(t, energy['potential'], color=props['color'], linestyle=props['linestyle'],
                label=props['label'], linewidth=1.2, alpha=0.7)
    
    axes[0, 0].set_xlabel('Time')
    axes[0, 0].set_ylabel('|ΔE/E₀| (%)')
    axes[0, 0].set_title('Total Energy Conservation\n' + r'$\Delta E = |E(t) - E_0| / |E_0|$')
    axes[0, 0].legend(fontsize=6, ncol=2)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_yscale('log')
    
    # Add analytical reference for kinetic energy (should be 0 for equilibrium)
    axes[0, 1].axhline(y=analytical_E['kinetic'], color='green', linestyle=':', 
                       linewidth=2, alpha=0.7, label=f"Analytical ({analytical_E['kinetic']:.3f})")
    axes[0, 1].set_xlabel('Time')
    axes[0, 1].set_ylabel('Kinetic Energy')
    axes[0, 1].set_title('Kinetic Energy\n' + r'$E_{\mathrm{kin}} = \frac{1}{2}\sum_i m_i v_i^2$')
    axes[0, 1].legend(fontsize=6, ncol=2)
    axes[0, 1].grid(True, alpha=0.3)
    
    # Add analytical reference for thermal energy
    axes[1, 0].axhline(y=analytical_E['thermal'], color='green', linestyle=':', 
                       linewidth=2, alpha=0.7, label=f"Analytical ({analytical_E['thermal']:.4f})")
    axes[1, 0].set_xlabel('Time')
    axes[1, 0].set_ylabel('Thermal Energy')
    axes[1, 0].set_title('Thermal Energy\n' + r'$E_{\mathrm{th}} = \sum_i m_i u_i = \int \frac{P}{\gamma-1} dx$')
    axes[1, 0].legend(fontsize=6, ncol=2)
    axes[1, 0].grid(True, alpha=0.3)
    
    # Add analytical reference for gravitational energy
    axes[1, 1].axhline(y=analytical_E['gravitational'], color='green', linestyle=':', 
                       linewidth=2, alpha=0.7, label=f"Analytical ({analytical_E['gravitational']:.4f})")
    axes[1, 1].set_xlabel('Time')
    axes[1, 1].set_ylabel('Potential Energy')
    axes[1, 1].set_title('Gravitational Energy\n' + r'$E_{\mathrm{grav}} = \frac{1}{2}\sum_i m_i \Phi_i$')
    axes[1, 1].legend(fontsize=6, ncol=2)
    axes[1, 1].grid(True, alpha=0.3)
    
    # Add text box with formulas
    formula_text = (
        r'$E_{\mathrm{total}} = E_{\mathrm{kin}} + E_{\mathrm{th}} + E_{\mathrm{grav}}$' + '\n'
        r'Analytical: $\Phi(x) = -2\pi G \int |x-x^\prime| \rho(x^\prime) dx^\prime$'
    )
    fig.text(0.5, 0.02, formula_text, ha='center', fontsize=10, 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    n_methods = 6 + len([m for m in CSMOOTH_METHODS if m in csmooth_data])
    plt.suptitle(f'Energy Evolution: {n_methods}-Method Comparison (incl. C_smooth=2)\n(Analytical: $E_{{total}}$ = {analytical_E["total"]:.4f})', fontsize=14)
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.savefig(os.path.join(OUTPUT_DIR, 'energy_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/energy_evolution.png")
    
    # =========================================================================
    # Plot 5b: Momentum Conservation
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot C_smooth=1 methods
    for method_name, props in METHODS.items():
        momentum = method_data[method_name]['momentum']
        momentum_vals = [p for p in momentum if p is not None]
        t = times[:len(momentum_vals)]
        ax.plot(t, momentum_vals, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    # Plot C_smooth=2 methods
    for method_name, props in CSMOOTH_METHODS.items():
        if method_name not in csmooth_data:
            continue
        momentum = csmooth_data[method_name]['momentum']
        momentum_vals = [p for p in momentum if p is not None]
        t = csmooth_times[:len(momentum_vals)]
        ax.plot(t, momentum_vals, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=5, alpha=0.7,
                label=props['label'], linewidth=1.5, markersize=5)
    
    # Analytical momentum (should be 0 for symmetric equilibrium)
    ax.axhline(y=0.0, color='green', linestyle=':', linewidth=2, alpha=0.7,
               label='Analytical (P = 0)')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel(r'Total Momentum $P$', fontsize=12)
    
    # Add formula text box
    formula_text = (
        r'$P = \sum_i m_i v_i$' + '\n'
        r'For symmetric system: $P_{\mathrm{analytical}} = 0$' + '\n'
        r'Conservation: $\frac{dP}{dt} = \sum_i F_i^{\mathrm{ext}} = 0$'
    )
    ax.text(0.98, 0.98, formula_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    n_methods = 6 + len([m for m in CSMOOTH_METHODS if m in csmooth_data])
    ax.set_title(f'Momentum Conservation: {n_methods}-Method Comparison (incl. C_smooth=2)\n(Symmetric System: P should remain ≈ 0)', fontsize=14)
    ax.legend(fontsize=8, loc='upper left', ncol=2)
    ax.grid(True, alpha=0.3)
    
    # Use symlog scale to show both positive and negative values near zero
    ax.set_yscale('symlog', linthresh=1e-10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'momentum_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/momentum_evolution.png")
    
    # =========================================================================
    # Plot 6: Density profiles at different times
    # =========================================================================
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    # Snapshot indices: 0, 20, 40, 60, 80, 100 → t = 0, 100, 200, 300, 400, 500
    time_indices = [0, 20, 40, 60, 80, 100]
    
    for idx, ti in enumerate(time_indices):
        ax = axes[idx]
        
        # Analytical solution (static reference)
        ax.fill_between(x_analytical, 0, rho_analytical, alpha=0.15, color='green')
        ax.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5,
                label='Analytical')
        
        # Simulation results for each method
        for method_name, props in METHODS.items():
            snaps = method_data[method_name]['snaps']
            if ti < len(snaps) and snaps[ti] is not None:
                snap = snaps[ti].sort_values('pos_x')
                ax.plot(snap['pos_x'], snap['dens'], 
                        color=props['color'], linestyle=props['linestyle'],
                        label=props['label'], linewidth=1.5, alpha=0.8)
        
        ax.set_xlabel('x', fontsize=10)
        ax.set_ylabel('ρ', fontsize=10)
        ax.set_title(f't = {ti * 5.0:.0f}', fontsize=11)
        if idx == 0:
            ax.legend(fontsize=7, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.5])
    
    plt.suptitle('Density Profiles: 6-Method Comparison (γ=5/3, n=1.5)\n(Green = Analytical Lane-Emden Solution)', 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'density_profiles.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/density_profiles.png")
    
    # =========================================================================
    # Animation: All 8 methods comparison (6 standard + 2 GSPH C_smooth=2)
    # =========================================================================
    print("\nGenerating 8-method animation...")
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    
    # Setup: one panel per method (4 rows × 2 columns)
    # Row 1-3: Standard 6 methods, Row 4: GSPH C_smooth=2 variants
    lines = {}
    
    method_order = ['gsph_gradh', 'gsph_nogradh', 'ssph_gradh', 'ssph_nogradh', 'gdisph_gradh', 'gdisph_nogradh']
    csmooth_order = ['gsph_gradh_cs2', 'gsph_nogradh_cs2']
    
    # Standard methods (rows 0-2)
    for i, method_name in enumerate(method_order):
        ax = axes[i // 2, i % 2]
        props = METHODS[method_name]
        
        # Analytical background
        ax.fill_between(x_analytical, 0, rho_analytical, alpha=0.2, color='green')
        ax.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5)
        
        # Initial line
        snap = method_data[method_name]['snaps'][0]
        if snap is not None:
            snap = snap.sort_values('pos_x')
            line, = ax.plot(snap['pos_x'], snap['dens'], 
                           color=props['color'], linewidth=2)
            lines[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.5])
        ax.set_xlabel('x')
        ax.set_ylabel('ρ')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    # C_smooth=2 methods (row 3)
    for i, method_name in enumerate(csmooth_order):
        ax = axes[3, i]
        props = CSMOOTH_METHODS[method_name]
        
        # Analytical background
        ax.fill_between(x_analytical, 0, rho_analytical, alpha=0.2, color='green')
        ax.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5)
        
        # Initial line (if data exists)
        if method_name in csmooth_data and len(csmooth_data[method_name]['snaps']) > 0:
            snap = csmooth_data[method_name]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line, = ax.plot(snap['pos_x'], snap['dens'], 
                               color=props['color'], linewidth=2)
                lines[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.5])
        ax.set_xlabel('x')
        ax.set_ylabel('ρ')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    time_text = fig.suptitle('8-Method Comparison: t = 0.0', fontsize=14)
    
    def update(frame):
        # Update standard methods
        for method_name in method_order:
            if method_name not in lines:
                continue
            snaps = method_data[method_name]['snaps']
            if frame < len(snaps) and snaps[frame] is not None:
                snap = snaps[frame].sort_values('pos_x')
                lines[method_name].set_data(snap['pos_x'], snap['dens'])
        
        # Update C_smooth methods (use frame % csmooth_num_snapshots for shorter data)
        csmooth_num = len(csmooth_times) if len(csmooth_times) > 0 else 1
        csmooth_frame = min(frame, csmooth_num - 1)
        for method_name in csmooth_order:
            if method_name not in lines or method_name not in csmooth_data:
                continue
            snaps = csmooth_data[method_name]['snaps']
            if csmooth_frame < len(snaps) and snaps[csmooth_frame] is not None:
                snap = snaps[csmooth_frame].sort_values('pos_x')
                lines[method_name].set_data(snap['pos_x'], snap['dens'])
        
        time_text.set_text(f'8-Method Comparison: t = {frame * 5.0:.0f}')
        return list(lines.values()) + [time_text]
    
    anim = FuncAnimation(fig, update, frames=num_snapshots, interval=100, blit=False)
    
    plt.tight_layout()
    anim_path = os.path.join(OUTPUT_DIR, '8method_comparison.gif')
    anim.save(anim_path, writer=PillowWriter(fps=10))
    plt.close()
    print(f"Saved: {anim_path}")
    
    # =========================================================================
    # Animation: Pressure Profile (8 methods)
    # =========================================================================
    print("Generating 8-method pressure animation...")
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    lines_pres = {}
    
    # Compute analytical pressure profile
    P_analytical = K * rho_analytical**GAMMA
    
    # Standard methods (rows 0-2)
    for i, method_name in enumerate(method_order):
        ax = axes[i // 2, i % 2]
        props = METHODS[method_name]
        
        # Analytical background
        ax.fill_between(x_analytical, 0, P_analytical, alpha=0.2, color='green')
        ax.plot(x_analytical, P_analytical, 'g-', linewidth=1.5, alpha=0.5)
        
        snap = method_data[method_name]['snaps'][0]
        if snap is not None:
            snap = snap.sort_values('pos_x')
            line, = ax.plot(snap['pos_x'], snap['pres'], 
                           color=props['color'], linewidth=2)
            lines_pres[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.2])
        ax.set_xlabel('x')
        ax.set_ylabel('P')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    # C_smooth=2 methods (row 3)
    for i, method_name in enumerate(csmooth_order):
        ax = axes[3, i]
        props = CSMOOTH_METHODS[method_name]
        
        # Analytical background
        ax.fill_between(x_analytical, 0, P_analytical, alpha=0.2, color='green')
        ax.plot(x_analytical, P_analytical, 'g-', linewidth=1.5, alpha=0.5)
        
        if method_name in csmooth_data and len(csmooth_data[method_name]['snaps']) > 0:
            snap = csmooth_data[method_name]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line, = ax.plot(snap['pos_x'], snap['pres'], 
                               color=props['color'], linewidth=2)
                lines_pres[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.2])
        ax.set_xlabel('x')
        ax.set_ylabel('P')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    time_text_pres = fig.suptitle('8-Method Pressure Profile: t = 0.0', fontsize=14)
    
    def update_pres(frame):
        # Update standard methods
        for method_name in method_order:
            if method_name not in lines_pres:
                continue
            snaps = method_data[method_name]['snaps']
            if frame < len(snaps) and snaps[frame] is not None:
                snap = snaps[frame].sort_values('pos_x')
                lines_pres[method_name].set_data(snap['pos_x'], snap['pres'])
        
        # Update C_smooth methods
        csmooth_num = len(csmooth_times) if len(csmooth_times) > 0 else 1
        csmooth_frame = min(frame, csmooth_num - 1)
        for method_name in csmooth_order:
            if method_name not in lines_pres or method_name not in csmooth_data:
                continue
            snaps = csmooth_data[method_name]['snaps']
            if csmooth_frame < len(snaps) and snaps[csmooth_frame] is not None:
                snap = snaps[csmooth_frame].sort_values('pos_x')
                lines_pres[method_name].set_data(snap['pos_x'], snap['pres'])
        
        time_text_pres.set_text(f'8-Method Pressure Profile: t = {frame * 5.0:.0f}')
        return list(lines_pres.values()) + [time_text_pres]
    
    anim_pres = FuncAnimation(fig, update_pres, frames=num_snapshots, interval=100, blit=False)
    
    plt.tight_layout()
    anim_path_pres = os.path.join(OUTPUT_DIR, '8method_pressure.gif')
    anim_pres.save(anim_path_pres, writer=PillowWriter(fps=10))
    plt.close()
    print(f"Saved: {anim_path_pres}")
    
    # =========================================================================
    # Animation: Gravitational Acceleration (8 methods)
    # =========================================================================
    print("Generating 8-method gravity force animation...")
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    lines_grav = {}
    
    # Standard methods (rows 0-2)
    for i, method_name in enumerate(method_order):
        ax = axes[i // 2, i % 2]
        props = METHODS[method_name]
        
        snap = method_data[method_name]['snaps'][0]
        if snap is not None:
            snap = snap.sort_values('pos_x')
            line, = ax.plot(snap['pos_x'], snap['grav_acc_x'], 
                           color=props['color'], linewidth=2)
            lines_grav[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([-3, 3])
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('g_x (grav)')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    # C_smooth=2 methods (row 3)
    for i, method_name in enumerate(csmooth_order):
        ax = axes[3, i]
        props = CSMOOTH_METHODS[method_name]
        
        if method_name in csmooth_data and len(csmooth_data[method_name]['snaps']) > 0:
            snap = csmooth_data[method_name]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line, = ax.plot(snap['pos_x'], snap['grav_acc_x'], 
                               color=props['color'], linewidth=2)
                lines_grav[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([-3, 3])
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('g_x (grav)')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    time_text_grav = fig.suptitle('8-Method Gravitational Acceleration: t = 0.0', fontsize=14)
    
    def update_grav(frame):
        # Update standard methods
        for method_name in method_order:
            if method_name not in lines_grav:
                continue
            snaps = method_data[method_name]['snaps']
            if frame < len(snaps) and snaps[frame] is not None:
                snap = snaps[frame].sort_values('pos_x')
                lines_grav[method_name].set_data(snap['pos_x'], snap['grav_acc_x'])
        
        # Update C_smooth methods
        csmooth_num = len(csmooth_times) if len(csmooth_times) > 0 else 1
        csmooth_frame = min(frame, csmooth_num - 1)
        for method_name in csmooth_order:
            if method_name not in lines_grav or method_name not in csmooth_data:
                continue
            snaps = csmooth_data[method_name]['snaps']
            if csmooth_frame < len(snaps) and snaps[csmooth_frame] is not None:
                snap = snaps[csmooth_frame].sort_values('pos_x')
                lines_grav[method_name].set_data(snap['pos_x'], snap['grav_acc_x'])
        
        time_text_grav.set_text(f'8-Method Gravitational Acceleration: t = {frame * 5.0:.0f}')
        return list(lines_grav.values()) + [time_text_grav]
    
    anim_grav = FuncAnimation(fig, update_grav, frames=num_snapshots, interval=100, blit=False)
    
    plt.tight_layout()
    anim_path_grav = os.path.join(OUTPUT_DIR, '8method_gravity.gif')
    anim_grav.save(anim_path_grav, writer=PillowWriter(fps=10))
    plt.close()
    print(f"Saved: {anim_path_grav}")
    
    # =========================================================================
    # Animation: Hydro Acceleration (Pressure Force) (8 methods)
    # =========================================================================
    print("Generating 8-method hydro force animation...")
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    lines_hydro = {}
    
    # Standard methods (rows 0-2)
    for i, method_name in enumerate(method_order):
        ax = axes[i // 2, i % 2]
        props = METHODS[method_name]
        
        snap = method_data[method_name]['snaps'][0]
        if snap is not None:
            snap = snap.sort_values('pos_x')
            # Hydro acceleration = total acceleration - gravitational acceleration
            acc_hydro = snap['acc_x'] - snap['grav_acc_x']
            line, = ax.plot(snap['pos_x'], acc_hydro, 
                           color=props['color'], linewidth=2)
            lines_hydro[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([-3, 3])
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('a_x (hydro = -∇P/ρ)')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    # C_smooth=2 methods (row 3)
    for i, method_name in enumerate(csmooth_order):
        ax = axes[3, i]
        props = CSMOOTH_METHODS[method_name]
        
        if method_name in csmooth_data and len(csmooth_data[method_name]['snaps']) > 0:
            snap = csmooth_data[method_name]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                acc_hydro = snap['acc_x'] - snap['grav_acc_x']
                line, = ax.plot(snap['pos_x'], acc_hydro, 
                               color=props['color'], linewidth=2)
                lines_hydro[method_name] = line
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([-3, 3])
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('a_x (hydro = -∇P/ρ)')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    time_text_hydro = fig.suptitle('8-Method Hydro Acceleration: t = 0.0', fontsize=14)
    
    def update_hydro(frame):
        # Update standard methods
        for method_name in method_order:
            if method_name not in lines_hydro:
                continue
            snaps = method_data[method_name]['snaps']
            if frame < len(snaps) and snaps[frame] is not None:
                snap = snaps[frame].sort_values('pos_x')
                acc_hydro = snap['acc_x'] - snap['grav_acc_x']
                lines_hydro[method_name].set_data(snap['pos_x'], acc_hydro)
        
        # Update C_smooth methods
        csmooth_num = len(csmooth_times) if len(csmooth_times) > 0 else 1
        csmooth_frame = min(frame, csmooth_num - 1)
        for method_name in csmooth_order:
            if method_name not in lines_hydro or method_name not in csmooth_data:
                continue
            snaps = csmooth_data[method_name]['snaps']
            if csmooth_frame < len(snaps) and snaps[csmooth_frame] is not None:
                snap = snaps[csmooth_frame].sort_values('pos_x')
                acc_hydro = snap['acc_x'] - snap['grav_acc_x']
                lines_hydro[method_name].set_data(snap['pos_x'], acc_hydro)
        
        time_text_hydro.set_text(f'8-Method Hydro Acceleration: t = {frame * 5.0:.0f}')
        return list(lines_hydro.values()) + [time_text_hydro]
    
    anim_hydro = FuncAnimation(fig, update_hydro, frames=num_snapshots, interval=100, blit=False)
    
    plt.tight_layout()
    anim_path_hydro = os.path.join(OUTPUT_DIR, '8method_hydro_force.gif')
    anim_hydro.save(anim_path_hydro, writer=PillowWriter(fps=10))
    plt.close()
    print(f"Saved: {anim_path_hydro}")
    
    # =========================================================================
    # Animation: Force Balance (Gravity + Hydro = Total) (8 methods)
    # =========================================================================
    print("Generating 8-method force balance animation...")
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    lines_balance = {}
    
    # Standard methods (rows 0-2)
    for i, method_name in enumerate(method_order):
        ax = axes[i // 2, i % 2]
        props = METHODS[method_name]
        
        snap = method_data[method_name]['snaps'][0]
        if snap is not None:
            snap = snap.sort_values('pos_x')
            # Plot gravity (red), hydro (blue), total (black)
            line_g, = ax.plot(snap['pos_x'], snap['grav_acc_x'], 
                             color='red', linewidth=1.5, alpha=0.7, label='Gravity')
            acc_hydro = snap['acc_x'] - snap['grav_acc_x']
            line_h, = ax.plot(snap['pos_x'], acc_hydro, 
                             color='blue', linewidth=1.5, alpha=0.7, label='Hydro')
            line_t, = ax.plot(snap['pos_x'], snap['acc_x'], 
                             color='black', linewidth=2, label='Total')
            lines_balance[method_name] = (line_g, line_h, line_t)
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([-3, 3])
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('Acceleration')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8, loc='upper right')
    
    # C_smooth=2 methods (row 3)
    for i, method_name in enumerate(csmooth_order):
        ax = axes[3, i]
        props = CSMOOTH_METHODS[method_name]
        
        if method_name in csmooth_data and len(csmooth_data[method_name]['snaps']) > 0:
            snap = csmooth_data[method_name]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line_g, = ax.plot(snap['pos_x'], snap['grav_acc_x'], 
                                 color='red', linewidth=1.5, alpha=0.7, label='Gravity')
                acc_hydro = snap['acc_x'] - snap['grav_acc_x']
                line_h, = ax.plot(snap['pos_x'], acc_hydro, 
                                 color='blue', linewidth=1.5, alpha=0.7, label='Hydro')
                line_t, = ax.plot(snap['pos_x'], snap['acc_x'], 
                                 color='black', linewidth=2, label='Total')
                lines_balance[method_name] = (line_g, line_h, line_t)
        
        ax.set_xlim([-8, 8])
        ax.set_ylim([-3, 3])
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('x')
        ax.set_ylabel('Acceleration')
        ax.set_title(props['label'])
        ax.grid(True, alpha=0.3)
    
    time_text_balance = fig.suptitle('8-Method Force Balance: t = 0.0\n(Red=Gravity, Blue=Hydro, Black=Total)', fontsize=14)
    
    def update_balance(frame):
        # Update standard methods
        for method_name in method_order:
            if method_name not in lines_balance:
                continue
            snaps = method_data[method_name]['snaps']
            if frame < len(snaps) and snaps[frame] is not None:
                snap = snaps[frame].sort_values('pos_x')
                acc_hydro = snap['acc_x'] - snap['grav_acc_x']
                lines_balance[method_name][0].set_data(snap['pos_x'], snap['grav_acc_x'])
                lines_balance[method_name][1].set_data(snap['pos_x'], acc_hydro)
                lines_balance[method_name][2].set_data(snap['pos_x'], snap['acc_x'])
        
        # Update C_smooth methods
        csmooth_num = len(csmooth_times) if len(csmooth_times) > 0 else 1
        csmooth_frame = min(frame, csmooth_num - 1)
        for method_name in csmooth_order:
            if method_name not in lines_balance or method_name not in csmooth_data:
                continue
            snaps = csmooth_data[method_name]['snaps']
            if csmooth_frame < len(snaps) and snaps[csmooth_frame] is not None:
                snap = snaps[csmooth_frame].sort_values('pos_x')
                acc_hydro = snap['acc_x'] - snap['grav_acc_x']
                lines_balance[method_name][0].set_data(snap['pos_x'], snap['grav_acc_x'])
                lines_balance[method_name][1].set_data(snap['pos_x'], acc_hydro)
                lines_balance[method_name][2].set_data(snap['pos_x'], snap['acc_x'])
        
        time_text_balance.set_text(f'8-Method Force Balance: t = {frame * 5.0:.0f}\n(Red=Gravity, Blue=Hydro, Black=Total)')
        return [line for lines in lines_balance.values() for line in lines] + [time_text_balance]
    
    anim_balance = FuncAnimation(fig, update_balance, frames=num_snapshots, interval=100, blit=False)
    
    plt.tight_layout()
    anim_path_balance = os.path.join(OUTPUT_DIR, '8method_force_balance.gif')
    anim_balance.save(anim_path_balance, writer=PillowWriter(fps=10))
    plt.close()
    print(f"Saved: {anim_path_balance}")
    
    # =========================================================================
    # NEW: 4×2 Animated Composite GIF with C_smooth=2 in 4th row
    # Layout:
    #   Row 1: GSPH+gradh (left), GSPH-gradh (right)
    #   Row 2: SSPH+gradh (left), SSPH-gradh (right)
    #   Row 3: GDISPH+gradh (left), GDISPH-gradh (right)
    #   Row 4: C_smooth=2 GSPH+gradh (left), Energy+Momentum summary (right)
    # =========================================================================
    if csmooth_data:
        print("\nGenerating 4×2 animated composite with C_smooth=2...")
        
        fig = plt.figure(figsize=(14, 18))
        
        # Create 4×2 grid with special handling for row 4
        gs = fig.add_gridspec(4, 2, hspace=0.25, wspace=0.2)
        
        # Rows 1-3: 6 method comparison (density profiles)
        method_order_4x2 = ['gsph_gradh', 'gsph_nogradh', 'ssph_gradh', 'ssph_nogradh', 'gdisph_gradh', 'gdisph_nogradh']
        lines_4x2 = {}
        axes_4x2 = []
        
        for i, method_name in enumerate(method_order_4x2):
            row = i // 2
            col = i % 2
            ax = fig.add_subplot(gs[row, col])
            axes_4x2.append(ax)
            props = METHODS[method_name]
            
            # Analytical background
            ax.fill_between(x_analytical, 0, rho_analytical, alpha=0.2, color='green')
            ax.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5)
            
            # Initial line
            snap = method_data[method_name]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line, = ax.plot(snap['pos_x'], snap['dens'], 
                               color=props['color'], linewidth=2)
                lines_4x2[method_name] = line
            
            ax.set_xlim([-8, 8])
            ax.set_ylim([0, 1.5])
            ax.set_xlabel('x')
            ax.set_ylabel('ρ')
            ax.set_title(f'{props["label"]} (C=1)')
            ax.grid(True, alpha=0.3)
        
        # Row 4: Both GSPH C_smooth=2 variants
        # Row 4, Left: GSPH+gradh with C_smooth=2
        ax_cs2_left = fig.add_subplot(gs[3, 0])
        cs2_method_left = 'gsph_gradh_cs2'
        if cs2_method_left in csmooth_data:
            cs2_props = CSMOOTH_METHODS[cs2_method_left]
            ax_cs2_left.fill_between(x_analytical, 0, rho_analytical, alpha=0.2, color='green')
            ax_cs2_left.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5)
            
            snap = csmooth_data[cs2_method_left]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line_cs2_left, = ax_cs2_left.plot(snap['pos_x'], snap['dens'], 
                                       color=cs2_props['color'], linewidth=2)
                lines_4x2[cs2_method_left] = line_cs2_left
            
            ax_cs2_left.set_xlim([-8, 8])
            ax_cs2_left.set_ylim([0, 1.5])
            ax_cs2_left.set_xlabel('x')
            ax_cs2_left.set_ylabel('ρ')
            ax_cs2_left.set_title(f'{cs2_props["label"]}')
            ax_cs2_left.grid(True, alpha=0.3)
        
        # Row 4, Right: GSPH-gradh with C_smooth=2
        ax_cs2_right = fig.add_subplot(gs[3, 1])
        cs2_method_right = 'gsph_nogradh_cs2'
        if cs2_method_right in csmooth_data:
            cs2_props = CSMOOTH_METHODS[cs2_method_right]
            ax_cs2_right.fill_between(x_analytical, 0, rho_analytical, alpha=0.2, color='green')
            ax_cs2_right.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5)
            
            snap = csmooth_data[cs2_method_right]['snaps'][0]
            if snap is not None:
                snap = snap.sort_values('pos_x')
                line_cs2_right, = ax_cs2_right.plot(snap['pos_x'], snap['dens'], 
                                       color=cs2_props['color'], linewidth=2)
                lines_4x2[cs2_method_right] = line_cs2_right
            
            ax_cs2_right.set_xlim([-8, 8])
            ax_cs2_right.set_ylim([0, 1.5])
            ax_cs2_right.set_xlabel('x')
            ax_cs2_right.set_ylabel('ρ')
            ax_cs2_right.set_title(f'{cs2_props["label"]}')
            ax_cs2_right.grid(True, alpha=0.3)
        
        time_text_4x2 = fig.suptitle('4×2 Method Comparison: t = 0.0\n(Rows 1-3: C_smooth=1, Row 4: GSPH C_smooth=2)', fontsize=14)
        
        def update_4x2(frame):
            # Update rows 1-3 (C_smooth=1 methods)
            for method_name in method_order_4x2:
                snaps = method_data[method_name]['snaps']
                if method_name in lines_4x2 and frame < len(snaps) and snaps[frame] is not None:
                    snap = snaps[frame].sort_values('pos_x')
                    lines_4x2[method_name].set_data(snap['pos_x'], snap['dens'])
            
            # Update row 4 (C_smooth=2 methods)
            csmooth_num = len(csmooth_times) if len(csmooth_times) > 0 else 1
            cs2_frame = min(frame, csmooth_num - 1)
            
            for cs2_method in [cs2_method_left, cs2_method_right]:
                if cs2_method in csmooth_data and cs2_method in lines_4x2:
                    cs2_snaps = csmooth_data[cs2_method]['snaps']
                    if cs2_frame < len(cs2_snaps) and cs2_snaps[cs2_frame] is not None:
                        snap = cs2_snaps[cs2_frame].sort_values('pos_x')
                        lines_4x2[cs2_method].set_data(snap['pos_x'], snap['dens'])
            
            current_time = frame * 5.0
            time_text_4x2.set_text(f'4×2 Method Comparison: t = {current_time:.0f}\n(Rows 1-3: C_smooth=1, Row 4: GSPH C_smooth=2)')
            return list(lines_4x2.values()) + [time_text_4x2]
        
        anim_4x2 = FuncAnimation(fig, update_4x2, frames=num_snapshots, interval=100, blit=False)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        anim_path_4x2 = os.path.join(OUTPUT_DIR, '4x2_csmooth_comparison.gif')
        anim_4x2.save(anim_path_4x2, writer=PillowWriter(fps=10))
        plt.close()
        print(f"Saved: {anim_path_4x2}")
    else:
        print("\n⚠️  No C_smooth=2 data available, skipping 4×2 composite animation")
    
    # =========================================================================
    # C_smooth Comparison Plot: Side-by-side h profiles
    # =========================================================================
    if csmooth_data:
        print("\nGenerating C_smooth comparison plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Get final snapshots for comparison
        final_frame = -1  # Last snapshot
        
        # Top row: GSPH+gradh (C=1 vs C=2)
        ax = axes[0, 0]
        ax.fill_between(x_analytical, 0, rho_analytical, alpha=0.15, color='green')
        ax.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5, label='Analytical')
        
        snap_c1 = method_data['gsph_gradh']['snaps'][final_frame]
        if snap_c1 is not None:
            snap_c1 = snap_c1.sort_values('pos_x')
            ax.plot(snap_c1['pos_x'], snap_c1['dens'], 
                   color=METHODS['gsph_gradh']['color'], linewidth=2, 
                   label='GSPH+gradh (C=1)')
        
        if 'gsph_gradh_cs2' in csmooth_data:
            snap_c2 = csmooth_data['gsph_gradh_cs2']['snaps'][final_frame]
            if snap_c2 is not None:
                snap_c2 = snap_c2.sort_values('pos_x')
                ax.plot(snap_c2['pos_x'], snap_c2['dens'], 
                       color=CSMOOTH_METHODS['gsph_gradh_cs2']['color'], linewidth=2, 
                       linestyle='--', label='GSPH+gradh (C=2)')
        
        ax.set_xlabel('x')
        ax.set_ylabel('ρ')
        ax.set_title('GSPH + grad-h: C_smooth Comparison')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.5])
        
        # Top right: GSPH-gradh (C=1 vs C=2)
        ax = axes[0, 1]
        ax.fill_between(x_analytical, 0, rho_analytical, alpha=0.15, color='green')
        ax.plot(x_analytical, rho_analytical, 'g-', linewidth=1.5, alpha=0.5, label='Analytical')
        
        snap_c1 = method_data['gsph_nogradh']['snaps'][final_frame]
        if snap_c1 is not None:
            snap_c1 = snap_c1.sort_values('pos_x')
            ax.plot(snap_c1['pos_x'], snap_c1['dens'], 
                   color=METHODS['gsph_nogradh']['color'], linewidth=2, 
                   label='GSPH-gradh (C=1)')
        
        if 'gsph_nogradh_cs2' in csmooth_data:
            snap_c2 = csmooth_data['gsph_nogradh_cs2']['snaps'][final_frame]
            if snap_c2 is not None:
                snap_c2 = snap_c2.sort_values('pos_x')
                ax.plot(snap_c2['pos_x'], snap_c2['dens'], 
                       color=CSMOOTH_METHODS['gsph_nogradh_cs2']['color'], linewidth=2, 
                       linestyle='--', label='GSPH-gradh (C=2)')
        
        ax.set_xlabel('x')
        ax.set_ylabel('ρ')
        ax.set_title('GSPH - grad-h: C_smooth Comparison')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([-8, 8])
        ax.set_ylim([0, 1.5])
        
        # Bottom left: Smoothing length h comparison
        ax = axes[1, 0]
        
        if snap_c1 is not None and 'sml' in snap_c1.columns:
            snap_c1 = method_data['gsph_gradh']['snaps'][final_frame].sort_values('pos_x')
            ax.plot(snap_c1['pos_x'], snap_c1['sml'], 
                   color=METHODS['gsph_gradh']['color'], linewidth=2, 
                   label='GSPH+gradh (C=1)')
        
        if 'gsph_gradh_cs2' in csmooth_data:
            snap_c2 = csmooth_data['gsph_gradh_cs2']['snaps'][final_frame]
            if snap_c2 is not None and 'sml' in snap_c2.columns:
                snap_c2 = snap_c2.sort_values('pos_x')
                ax.plot(snap_c2['pos_x'], snap_c2['sml'], 
                       color=CSMOOTH_METHODS['gsph_gradh_cs2']['color'], linewidth=2, 
                       linestyle='--', label='GSPH+gradh (C=2)')
        
        ax.set_xlabel('x')
        ax.set_ylabel('Smoothing Length h')
        ax.set_title('Smoothing Length Profile: Effect of C_smooth')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([-8, 8])
        
        # Bottom right: grad-h (Omega) comparison
        ax = axes[1, 1]
        
        snap_c1 = method_data['gsph_gradh']['snaps'][final_frame]
        if snap_c1 is not None and 'gradh' in snap_c1.columns:
            snap_c1 = snap_c1.sort_values('pos_x')
            ax.plot(snap_c1['pos_x'], snap_c1['gradh'], 
                   color=METHODS['gsph_gradh']['color'], linewidth=2, 
                   label='GSPH+gradh (C=1)')
        
        if 'gsph_gradh_cs2' in csmooth_data:
            snap_c2 = csmooth_data['gsph_gradh_cs2']['snaps'][final_frame]
            if snap_c2 is not None and 'gradh' in snap_c2.columns:
                snap_c2 = snap_c2.sort_values('pos_x')
                ax.plot(snap_c2['pos_x'], snap_c2['gradh'], 
                       color=CSMOOTH_METHODS['gsph_gradh_cs2']['color'], linewidth=2, 
                       linestyle='--', label='GSPH+gradh (C=2)')
        
        ax.axhline(y=1.0, color='green', linestyle=':', linewidth=2, alpha=0.7, label='Ω=1 (uniform)')
        ax.set_xlabel('x')
        ax.set_ylabel('Ω (grad-h factor)')
        ax.set_title('Grad-h Factor Profile: Effect of C_smooth')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([-8, 8])
        
        plt.suptitle('C_smooth = 2 Effect on Smoothing Length Adaptation\n(Final snapshot comparison)', fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(os.path.join(OUTPUT_DIR, 'csmooth_effect_comparison.png'), dpi=150)
        plt.close()
        print(f"Saved: {OUTPUT_DIR}/csmooth_effect_comparison.png")
    
    # =========================================================================
    # Summary statistics
    # =========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        initial = metrics[0]
        final = metrics[-1]
        
        rho_change = (final['rho_center'] / initial['rho_center'] - 1) * 100
        extent_change = ((final['x_max'] - final['x_min']) / 
                        (initial['x_max'] - initial['x_min']) - 1) * 100
        
        print(f"\n{props['label']}:")
        print(f"  Central density change: {rho_change:+.2f}%")
        print(f"  Extent change: {extent_change:+.2f}%")
        print(f"  Final max velocity: {final['max_vel']:.4f}")
        print(f"  Final max density: {final['max_dens']:.4f}")
        
        energy = method_data[method_name]['energy']
        if energy is not None:
            E0 = energy['total'][0]
            Ef = energy['total'][-1]
            E_error = np.abs(Ef - E0) / np.abs(E0) * 100
            print(f"  Energy conservation error: {E_error:.4f}%")
    
    # C_smooth = 2 statistics
    if csmooth_data:
        print("\n" + "-" * 70)
        print("C_smooth = 2 Methods:")
        print("-" * 70)
        
        for method_name, props in CSMOOTH_METHODS.items():
            if method_name not in csmooth_data:
                continue
            metrics = csmooth_data[method_name]['metrics']
            initial = metrics[0]
            final = metrics[-1]
            
            if initial is None or final is None:
                continue
            
            rho_change = (final['rho_center'] / initial['rho_center'] - 1) * 100
            extent_change = ((final['x_max'] - final['x_min']) / 
                            (initial['x_max'] - initial['x_min']) - 1) * 100
            
            print(f"\n{props['label']}:")
            print(f"  Central density change: {rho_change:+.2f}%")
            print(f"  Extent change: {extent_change:+.2f}%")
            print(f"  Final max velocity: {final['max_vel']:.4f}")
            print(f"  Final max density: {final['max_dens']:.4f}")
            
            energy = csmooth_data[method_name]['energy']
            if energy is not None:
                E0 = energy['total'][0]
                Ef = energy['total'][-1]
                E_error = np.abs(Ef - E0) / np.abs(E0) * 100
                print(f"  Energy conservation error: {E_error:.4f}%")
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
