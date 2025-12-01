#!/usr/bin/env python3
"""
Visualization script for IMBH-molecular cloud tidal disruption simulations

Generates:
1. Density evolution plots
2. Tidal deformation (axis ratios)
3. Mass loss and bound/unbound fraction
4. Temperature-density phase diagrams
5. Animations of cloud disruption
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse
from pathlib import Path
import glob

# Physics constants
DENSITY_TO_NH = 40.5  # Conversion factor: n_H = ρ [M_☉/pc³] * 40.5 [cm⁻³]

def load_snapshot(filepath):
    """Load SPH snapshot from CSV file"""
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            break
    
    # Load data
    arr = np.loadtxt(filepath, delimiter=',')
    
    data['id'] = arr[:, 0].astype(int)
    data['x'] = arr[:, 1]
    data['y'] = arr[:, 2]
    data['z'] = arr[:, 3]
    data['vx'] = arr[:, 4]
    data['vy'] = arr[:, 5]
    data['vz'] = arr[:, 6]
    data['mass'] = arr[:, 10]
    data['dens'] = arr[:, 11]
    data['pres'] = arr[:, 12]
    data['ene'] = arr[:, 13]
    
    return data

def compute_axis_ratios(x, y, z, mass):
    """Compute cloud axis ratios from moment of inertia tensor"""
    # Center of mass
    x_cm = np.sum(mass * x) / np.sum(mass)
    y_cm = np.sum(mass * y) / np.sum(mass)
    z_cm = np.sum(mass * z) / np.sum(mass)
    
    # Relative positions
    dx = x - x_cm
    dy = y - y_cm
    dz = z - z_cm
    
    # Inertia tensor
    Ixx = np.sum(mass * (dy**2 + dz**2))
    Iyy = np.sum(mass * (dx**2 + dz**2))
    Izz = np.sum(mass * (dx**2 + dy**2))
    Ixy = -np.sum(mass * dx * dy)
    Ixz = -np.sum(mass * dx * dz)
    Iyz = -np.sum(mass * dy * dz)
    
    I = np.array([[Ixx, Ixy, Ixz],
                  [Ixy, Iyy, Iyz],
                  [Ixz, Iyz, Izz]])
    
    # Eigenvalues give principal moments
    eigvals = np.linalg.eigvalsh(I)
    eigvals = np.sort(eigvals)[::-1]  # Sort descending
    
    # Semi-axes: a >= b >= c
    a = np.sqrt(5 * eigvals[0] / np.sum(mass))
    b = np.sqrt(5 * eigvals[1] / np.sum(mass))
    c = np.sqrt(5 * eigvals[2] / np.sum(mass))
    
    return a, b, c

def plot_density_evolution(results_dir, output_dir):
    """Plot density evolution over time"""
    snapshots = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    
    if len(snapshots) == 0:
        print(f"No snapshots found in {results_dir}")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    times = []
    max_densities = []
    mean_densities = []
    
    for snap in snapshots:
        data = load_snapshot(snap)
        time = float(Path(snap).stem.split('_')[1]) * 0.02  # Assuming dt=0.02
        
        times.append(time)
        max_densities.append(np.max(data['dens']))
        mean_densities.append(np.mean(data['dens']))
    
    ax.plot(times, max_densities, 'r-', label='Max density', linewidth=2)
    ax.plot(times, mean_densities, 'b--', label='Mean density', linewidth=2)
    ax.set_xlabel('Time [Myr]', fontsize=14)
    ax.set_ylabel(r'Density [$M_\odot$ pc$^{-3}$]', fontsize=14)
    ax.set_yscale('log')
    ax.legend(fontsize=12)
    ax.grid(alpha=0.3)
    ax.set_title('Cloud Density Evolution During Tidal Disruption', fontsize=16)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/density_evolution.png", dpi=150)
    plt.close()
    print(f"✓ Saved: {output_dir}/density_evolution.png")

def plot_axis_ratios(results_dir, output_dir):
    """Plot cloud deformation (axis ratios) over time"""
    snapshots = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    
    if len(snapshots) == 0:
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    times = []
    a_vals, b_vals, c_vals = [], [], []
    
    for snap in snapshots:
        data = load_snapshot(snap)
        time = float(Path(snap).stem.split('_')[1]) * 0.02
        
        a, b, c = compute_axis_ratios(data['x'], data['y'], data['z'], data['mass'])
        
        times.append(time)
        a_vals.append(a)
        b_vals.append(b)
        c_vals.append(c)
    
    # Plot absolute axes
    ax1.plot(times, a_vals, 'r-', label='a (long axis)', linewidth=2)
    ax1.plot(times, b_vals, 'g-', label='b (intermediate)', linewidth=2)
    ax1.plot(times, c_vals, 'b-', label='c (short axis)', linewidth=2)
    ax1.set_xlabel('Time [Myr]', fontsize=14)
    ax1.set_ylabel('Semi-axis length [pc]', fontsize=14)
    ax1.legend(fontsize=12)
    ax1.grid(alpha=0.3)
    ax1.set_title('Cloud Semi-Axes Evolution', fontsize=16)
    
    # Plot axis ratios
    b_over_a = np.array(b_vals) / np.array(a_vals)
    c_over_a = np.array(c_vals) / np.array(a_vals)
    
    ax2.plot(times, b_over_a, 'g-', label='b/a', linewidth=2)
    ax2.plot(times, c_over_a, 'b-', label='c/a (elongation)', linewidth=2)
    ax2.axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Spherical')
    ax2.set_xlabel('Time [Myr]', fontsize=14)
    ax2.set_ylabel('Axis Ratio', fontsize=14)
    ax2.set_ylim([0, 1.1])
    ax2.legend(fontsize=12)
    ax2.grid(alpha=0.3)
    ax2.set_title('Cloud Deformation (smaller = more elongated)', fontsize=16)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/axis_ratios.png", dpi=150)
    plt.close()
    print(f"✓ Saved: {output_dir}/axis_ratios.png")

def plot_phase_diagram(results_dir, output_dir, snapshot_idx=-1):
    """Plot density-temperature phase diagram"""
    snapshots = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    
    if len(snapshots) == 0:
        return
    
    data = load_snapshot(snapshots[snapshot_idx])
    
    # Convert to number density
    n_H = data['dens'] * DENSITY_TO_NH
    
    # Temperature from ideal gas: T = (γ-1) μ m_H u / k_B
    # In code units, T ∝ u (specific internal energy)
    gamma = 5.0 / 3.0
    u = data['ene']
    # Approximate T [K] from u (needs proper unit conversion)
    T_approx = u * 1e4  # Placeholder scaling
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 2D histogram
    H, xedges, yedges = np.histogram2d(
        np.log10(n_H), np.log10(T_approx),
        bins=50, weights=data['mass']
    )
    
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    im = ax.imshow(H.T, origin='lower', extent=extent, aspect='auto',
                   cmap='viridis', interpolation='bilinear')
    
    ax.set_xlabel(r'log n$_H$ [cm$^{-3}$]', fontsize=14)
    ax.set_ylabel(r'log T [K]', fontsize=14)
    ax.set_title(f'Density-Temperature Phase Diagram (t={snapshot_idx})', fontsize=16)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Mass [M$_\\odot$]', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/phase_diagram.png", dpi=150)
    plt.close()
    print(f"✓ Saved: {output_dir}/phase_diagram.png")

def create_animation(results_dir, output_dir):
    """Create animation of cloud disruption"""
    snapshots = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    
    if len(snapshots) < 10:
        print("Not enough snapshots for animation (need >= 10)")
        return
    
    fig = plt.figure(figsize=(12, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    def update(frame_idx):
        ax1.clear()
        ax2.clear()
        
        snap = snapshots[frame_idx]
        data = load_snapshot(snap)
        time = frame_idx * 0.02
        
        # XY projection
        ax1.scatter(data['x'], data['y'], c=np.log10(data['dens']), 
                   s=1, cmap='viridis', vmin=-2, vmax=2)
        ax1.set_xlabel('x [pc]', fontsize=12)
        ax1.set_ylabel('y [pc]', fontsize=12)
        ax1.set_aspect('equal')
        ax1.set_title(f'XY Plane (t = {time:.2f} Myr)', fontsize=14)
        ax1.set_xlim([-10, 10])
        ax1.set_ylim([-10, 10])
        
        # XZ projection
        ax2.scatter(data['x'], data['z'], c=np.log10(data['dens']), 
                   s=1, cmap='viridis', vmin=-2, vmax=2)
        ax2.set_xlabel('x [pc]', fontsize=12)
        ax2.set_ylabel('z [pc]', fontsize=12)
        ax2.set_aspect('equal')
        ax2.set_title(f'XZ Plane (t = {time:.2f} Myr)', fontsize=14)
        ax2.set_xlim([-10, 10])
        ax2.set_ylim([-10, 10])
    
    anim = FuncAnimation(fig, update, frames=len(snapshots), interval=100)
    anim.save(f"{output_dir}/disruption_animation.gif", writer='pillow', fps=10)
    plt.close()
    print(f"✓ Saved: {output_dir}/disruption_animation.gif")

def main():
    parser = argparse.ArgumentParser(description='Visualize IMBH-cloud tidal disruption')
    parser.add_argument('--results-dir', type=str, required=True,
                       help='Directory containing simulation results')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Output directory for plots (default: results-dir/plots)')
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir) if args.output_dir else results_dir / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("  IMBH-Cloud Visualization")
    print("=" * 60)
    print(f"Results: {results_dir}")
    print(f"Output:  {output_dir}")
    print()
    
    # Generate plots
    print("Generating plots...")
    plot_density_evolution(results_dir, output_dir)
    plot_axis_ratios(results_dir, output_dir)
    plot_phase_diagram(results_dir, output_dir)
    
    print()
    print("Creating animation...")
    animations_dir = results_dir / 'animations'
    animations_dir.mkdir(parents=True, exist_ok=True)
    create_animation(results_dir, animations_dir)
    
    print()
    print("=" * 60)
    print("✓ All visualizations complete!")
    print("=" * 60)

if __name__ == '__main__':
    main()
