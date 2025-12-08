#!/usr/bin/env python3
"""
Compare GSPH/SSPH × grad-h/no-grad-h methods for 1D polytropic slab with kernel gravity.

This script generates:
1. Central density evolution comparison
2. Density profile comparison at different times
3. Energy conservation analysis
4. Animation comparing all 4 methods
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob

# Configuration
BASE_DIR = "sample/gradh_study_1d/results/kernel_gravity"
OUTPUT_DIR = "sample/gradh_study_1d/results/comparison"
RELAXATION_DIR = "sample/gradh_study_1d/results/relaxation"

# Methods to compare (using high-contrast colors)
METHODS = {
    'gsph_gradh': {'label': 'GSPH + grad-h', 'color': '#0066CC', 'linestyle': '-', 'marker': 'o'},      # Dark blue
    'gsph_nogradh': {'label': 'GSPH - grad-h', 'color': '#9933FF', 'linestyle': '--', 'marker': 's'},   # Purple
    'ssph_gradh': {'label': 'SSPH + grad-h', 'color': '#CC0000', 'linestyle': '-', 'marker': '^'},      # Dark red
    'ssph_nogradh': {'label': 'SSPH - grad-h', 'color': '#FF6600', 'linestyle': '--', 'marker': 'v'},   # Orange
}

# Physical parameters (from polytropic slab setup)
# γ = 5/3 polytrope with polytropic index n = 1/(γ-1) = 1.5
# This is the standard test case for hydrostatic equilibrium
RHO_CENTER = 1.0
K = 1.0
G = 1.0
GAMMA = 5.0/3.0  # γ = 5/3 → n = 1.5 polytrope


def solve_lane_emden(n, dxi=1e-4, max_steps=100000):
    """Solve planar Lane-Emden equation: d²θ/dξ² = -θⁿ"""
    xi_arr = [0.0]
    theta_arr = [1.0]
    dtheta_arr = [0.0]
    
    xi = 0.0
    theta = 1.0
    dtheta = 0.0
    
    for _ in range(max_steps):
        if theta <= 0:
            break
            
        # RK4
        def f1(xi, theta, phi): return phi
        def f2(xi, theta, phi): return -theta**n if theta > 0 else 0
        
        k1_theta = dxi * f1(xi, theta, dtheta)
        k1_phi = dxi * f2(xi, theta, dtheta)
        
        k2_theta = dxi * f1(xi + 0.5*dxi, theta + 0.5*k1_theta, dtheta + 0.5*k1_phi)
        k2_phi = dxi * f2(xi + 0.5*dxi, theta + 0.5*k1_theta, dtheta + 0.5*k1_phi)
        
        k3_theta = dxi * f1(xi + 0.5*dxi, theta + 0.5*k2_theta, dtheta + 0.5*k2_phi)
        k3_phi = dxi * f2(xi + 0.5*dxi, theta + 0.5*k2_theta, dtheta + 0.5*k2_phi)
        
        k4_theta = dxi * f1(xi + dxi, theta + k3_theta, dtheta + k3_phi)
        k4_phi = dxi * f2(xi + dxi, theta + k3_theta, dtheta + k3_phi)
        
        xi += dxi
        theta += (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta) / 6.0
        dtheta += (k1_phi + 2*k2_phi + 2*k3_phi + k4_phi) / 6.0
        
        xi_arr.append(xi)
        theta_arr.append(max(0.0, theta))
        dtheta_arr.append(dtheta)
    
    return np.array(xi_arr), np.array(theta_arr)


def get_analytical_profile(x_vals, rho_c, K, G, gamma):
    """Get analytical Lane-Emden density profile."""
    n = 1.0 / (gamma - 1.0)
    
    # Length scale: α² = K(n+1)ρ_c^(1-n) / (2πG)
    alpha_sq = K * (n + 1.0) * rho_c**(1.0 - n) / (2.0 * np.pi * G)
    alpha = np.sqrt(alpha_sq)
    
    # Solve Lane-Emden
    xi_le, theta_le = solve_lane_emden(n)
    xi_surface = xi_le[-1]
    
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
    print("with Kernel-Convolved 1D Gravity")
    print("=" * 70)
    
    # Get analytical solution
    x_analytical = np.linspace(-8, 8, 1000)
    rho_analytical, x_surface = get_analytical_profile(x_analytical, RHO_CENTER, K, G, GAMMA)
    print(f"Analytical solution: x_surface = {x_surface:.4f}")
    print(f"Polytropic index n = {1.0/(GAMMA-1.0):.2f}, γ = {GAMMA:.4f}")
    
    # Load all snapshots for each method
    method_data = {}
    num_snapshots = None
    
    for method_name in METHODS:
        method_dir = os.path.join(BASE_DIR, method_name)
        n_snaps = count_snapshots(method_dir)
        if num_snapshots is None:
            num_snapshots = n_snaps
        print(f"  {method_name}: {n_snaps} snapshots")
        
        snaps = []
        metrics = []
        for i in range(n_snaps):
            snap = read_snapshot(method_dir, i)
            snaps.append(snap)
            metrics.append(analyze_snapshot(snap))
        
        energy = read_energy_file(method_dir)
        
        method_data[method_name] = {
            'snaps': snaps,
            'metrics': metrics,
            'energy': energy
        }
    
    # Get time from first snapshot header
    first_snap = method_data['gsph_gradh']['snaps'][-1]
    if first_snap is not None:
        # Output interval: 5.0 for t=0→500 (101 snapshots)
        output_interval = 5.0
        times = np.arange(num_snapshots) * output_interval
    
    print(f"\nTime range: 0 to {times[-1]:.1f}")
    
    # =========================================================================
    # Plot 1: Central density evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        rho_center = [m['rho_center'] for m in metrics if m is not None]
        t = times[:len(rho_center)]
        ax.plot(t, rho_center, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    ax.axhline(y=RHO_CENTER, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label='Analytical ρ_c = 1.0')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Central Density ρ_c', fontsize=12)
    ax.set_title('Central Density Evolution: 4-Method Comparison\n(γ=5/3 Polytropic Slab, n=1.5, Kernel-Convolved 1D Gravity)', fontsize=14)
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
        metrics = method_data[method_name]['metrics']
        max_dens = [m['max_dens'] for m in metrics if m is not None]
        t = times[:len(max_dens)]
        ax.plot(t, max_dens, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    ax.axhline(y=RHO_CENTER, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label='Analytical ρ_c = 1.0')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Maximum Density ρ_max', fontsize=12)
    ax.set_title('Maximum Density Evolution: 4-Method Comparison\n(Indicator of Core Collapse)', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'max_density_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/max_density_evolution.png")
    
    # =========================================================================
    # Plot 3: Slab extent evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        extent = [(m['x_max'] - m['x_min'])/2 for m in metrics if m is not None]
        t = times[:len(extent)]
        ax.plot(t, extent, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    ax.axhline(y=x_surface, color='green', linestyle=':', alpha=0.7, 
               linewidth=2, label=f'Analytical x_surface = {x_surface:.3f}')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Slab Half-Width', fontsize=12)
    ax.set_title('Slab Extent Evolution: 4-Method Comparison', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'extent_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/extent_evolution.png")
    
    # =========================================================================
    # Plot 4: Max velocity evolution
    # =========================================================================
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for method_name, props in METHODS.items():
        metrics = method_data[method_name]['metrics']
        max_vel = [m['max_vel'] for m in metrics if m is not None]
        t = times[:len(max_vel)]
        ax.plot(t, max_vel, 
                color=props['color'], linestyle=props['linestyle'],
                marker=props['marker'], markevery=10,
                label=props['label'], linewidth=2, markersize=6)
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Maximum |v_x|', fontsize=12)
    ax.set_title('Maximum Velocity Evolution: 4-Method Comparison\n(Indicator of Stability)', fontsize=14)
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
    axes[0, 1].set_title('Kinetic Energy')
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('Time')
    axes[1, 0].set_ylabel('Thermal Energy')
    axes[1, 0].set_title('Thermal Energy')
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('Time')
    axes[1, 1].set_ylabel('Potential Energy')
    axes[1, 1].set_title('Gravitational Potential Energy')
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.suptitle('Energy Evolution: 4-Method Comparison', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'energy_evolution.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/energy_evolution.png")
    
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
    
    plt.suptitle('Density Profiles: 4-Method Comparison (γ=5/3, n=1.5)\n(Green = Analytical Lane-Emden Solution)', 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'density_profiles.png'), dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/density_profiles.png")
    
    # =========================================================================
    # Animation: All 4 methods comparison
    # =========================================================================
    print("\nGenerating animation...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Setup: one panel per method
    lines = {}
    titles = {}
    
    method_order = ['gsph_gradh', 'gsph_nogradh', 'ssph_gradh', 'ssph_nogradh']
    
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
        titles[method_name] = ax
    
    time_text = fig.suptitle('t = 0.0', fontsize=14)
    
    def update(frame):
        for method_name in method_order:
            snaps = method_data[method_name]['snaps']
            if frame < len(snaps) and snaps[frame] is not None:
                snap = snaps[frame].sort_values('pos_x')
                lines[method_name].set_data(snap['pos_x'], snap['dens'])
        time_text.set_text(f't = {frame * 5.0:.0f}')
        return list(lines.values()) + [time_text]
    
    anim = FuncAnimation(fig, update, frames=num_snapshots, interval=100, blit=False)
    
    plt.tight_layout()
    anim_path = os.path.join(OUTPUT_DIR, '4method_comparison.gif')
    anim.save(anim_path, writer=PillowWriter(fps=10))
    plt.close()
    print(f"Saved: {anim_path}")
    
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
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
