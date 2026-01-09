#!/usr/bin/env python3
"""
Monitor isothermal relaxation progress and generate diagnostic plots.

This script:
1. Reads snapshot CSV files from the relaxation simulation
2. Computes density profile errors vs analytical solution
3. Generates diagnostic plots showing convergence
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
RESULTS_DIR = "simulations/astrophysics/imbh_cloud/results/hvcc_10k/isothermal"
PLOT_DIR = "simulations/astrophysics/imbh_cloud/results/hvcc_10k/isothermal/plots"

def read_snapshot(filepath):
    """Read a snapshot CSV file and extract particle data."""
    data = {
        'step': 0,
        'relax_step': 0,
        'relax_total': 0,
        'time': 0.0,
        'particles': []
    }

    with open(filepath, 'r') as f:
        in_header = True
        columns = []

        for line in f:
            line = line.strip()

            if line.startswith('#'):
                # Parse header
                if 'Step:' in line and 'Relaxation' not in line:
                    try:
                        data['step'] = int(line.split(':')[1].strip())
                    except:
                        pass
                elif 'Relaxation Step:' in line:
                    try:
                        parts = line.split(':')[1].strip().split('/')
                        data['relax_step'] = int(parts[0].strip())
                        data['relax_total'] = int(parts[1].strip())
                    except:
                        pass
                elif 'Time (code):' in line:
                    try:
                        data['time'] = float(line.split(':')[1].strip())
                    except:
                        pass
                continue

            if in_header:
                # First non-comment line is column headers
                columns = line.split(',')
                in_header = False
                continue

            # Parse data line
            try:
                values = line.split(',')
                particle = {}
                for i, col in enumerate(columns):
                    if i < len(values):
                        try:
                            particle[col] = float(values[i])
                        except:
                            particle[col] = values[i]
                data['particles'].append(particle)
            except:
                pass

    return data

def analytical_density(r, rho_center, r_c):
    """Modified isothermal sphere density profile."""
    return rho_center / (1.0 + (r / r_c)**2)

def compute_profile_error(data, rho_center, r_c, R_cloud):
    """Compute error between SPH and analytical density profiles."""
    particles = data['particles']

    # Filter real particles (not ghost)
    real_particles = [p for p in particles if p.get('is_ghost', 0) == 0]

    if not real_particles:
        return None, None, None, None

    # Compute radii and densities
    radii = []
    sph_dens = []
    ana_dens = []

    for p in real_particles:
        x = p.get('pos_x', 0)
        y = p.get('pos_y', 0)
        z = p.get('pos_z', 0)
        r = np.sqrt(x**2 + y**2 + z**2)

        if r < R_cloud * 1.1:  # Only particles within cloud
            radii.append(r)
            sph_dens.append(p.get('dens', 0))
            ana_dens.append(analytical_density(r, rho_center, r_c))

    radii = np.array(radii)
    sph_dens = np.array(sph_dens)
    ana_dens = np.array(ana_dens)

    # Compute errors
    errors = np.abs(sph_dens - ana_dens) / ana_dens * 100  # Percent error

    return radii, sph_dens, ana_dens, errors

def create_diagnostic_plots(snapshots, rho_center, r_c, R_cloud, output_dir):
    """Create diagnostic plots for relaxation progress."""
    os.makedirs(output_dir, exist_ok=True)

    # Collect statistics across all snapshots
    steps = []
    mean_errors = []
    max_errors = []
    rms_errors = []

    for snap_file in snapshots:
        data = read_snapshot(snap_file)
        radii, sph_dens, ana_dens, errors = compute_profile_error(
            data, rho_center, r_c, R_cloud
        )

        if radii is None:
            continue

        steps.append(data['relax_step'])
        mean_errors.append(np.mean(errors))
        max_errors.append(np.max(errors))
        rms_errors.append(np.sqrt(np.mean(errors**2)))

    if not steps:
        print("No valid snapshots found")
        return

    # Sort by step
    sort_idx = np.argsort(steps)
    steps = np.array(steps)[sort_idx]
    mean_errors = np.array(mean_errors)[sort_idx]
    max_errors = np.array(max_errors)[sort_idx]
    rms_errors = np.array(rms_errors)[sort_idx]

    # Plot 1: Error convergence
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(steps, mean_errors, 'b-o', label='Mean Error', markersize=4)
    ax.semilogy(steps, rms_errors, 'g-s', label='RMS Error', markersize=4)
    ax.semilogy(steps, max_errors, 'r-^', label='Max Error', markersize=4)
    ax.set_xlabel('Relaxation Step')
    ax.set_ylabel('Density Error (%)')
    ax.set_title('Isothermal Relaxation Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(steps) * 1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'convergence.png'), dpi=150)
    plt.close()

    # Plot 2: Density profile comparison (latest snapshot)
    latest_snap = snapshots[sort_idx[-1]]
    data = read_snapshot(latest_snap)
    radii, sph_dens, ana_dens, errors = compute_profile_error(
        data, rho_center, r_c, R_cloud
    )

    if radii is not None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: Density profile
        ax1 = axes[0]
        # Bin particles by radius
        r_bins = np.linspace(0, R_cloud, 30)
        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

        binned_dens = []
        binned_std = []
        for i in range(len(r_bins) - 1):
            mask = (radii >= r_bins[i]) & (radii < r_bins[i+1])
            if np.sum(mask) > 0:
                binned_dens.append(np.mean(sph_dens[mask]))
                binned_std.append(np.std(sph_dens[mask]))
            else:
                binned_dens.append(np.nan)
                binned_std.append(0)

        binned_dens = np.array(binned_dens)
        binned_std = np.array(binned_std)

        # Analytical profile
        r_ana = np.linspace(0.001, R_cloud, 100)
        rho_ana = analytical_density(r_ana, rho_center, r_c)

        ax1.plot(r_ana, rho_ana, 'b-', linewidth=2, label='Analytical')
        ax1.errorbar(r_centers, binned_dens, yerr=binned_std, fmt='ro',
                     markersize=4, capsize=2, label='SPH (binned)')
        ax1.scatter(radii[::10], sph_dens[::10], s=1, c='gray', alpha=0.3, label='SPH particles')
        ax1.set_xlabel('Radius [pc]')
        ax1.set_ylabel('Density [M_sun/pc^3]')
        ax1.set_title(f'Density Profile (Step {data["relax_step"]})')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Right: Error vs radius
        ax2 = axes[1]
        ax2.scatter(radii, errors, s=1, c='blue', alpha=0.3)

        # Binned error
        binned_err = []
        for i in range(len(r_bins) - 1):
            mask = (radii >= r_bins[i]) & (radii < r_bins[i+1])
            if np.sum(mask) > 0:
                binned_err.append(np.mean(errors[mask]))
            else:
                binned_err.append(np.nan)

        ax2.plot(r_centers, binned_err, 'r-o', linewidth=2, markersize=4, label='Mean error')
        ax2.axhline(y=np.mean(errors), color='g', linestyle='--', label=f'Overall mean: {np.mean(errors):.1f}%')
        ax2.set_xlabel('Radius [pc]')
        ax2.set_ylabel('Density Error (%)')
        ax2.set_title('Error vs Radius')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, min(100, np.percentile(errors, 95) * 1.5))

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'profile_comparison.png'), dpi=150)
        plt.close()

    # Plot 3: Particle distribution (2D projection)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    particles = data['particles']
    real_particles = [p for p in particles if p.get('is_ghost', 0) == 0]
    ghost_particles = [p for p in particles if p.get('is_ghost', 0) == 1]

    x_real = [p.get('pos_x', 0) for p in real_particles]
    y_real = [p.get('pos_y', 0) for p in real_particles]
    z_real = [p.get('pos_z', 0) for p in real_particles]
    dens_real = [p.get('dens', 0) for p in real_particles]

    x_ghost = [p.get('pos_x', 0) for p in ghost_particles]
    y_ghost = [p.get('pos_y', 0) for p in ghost_particles]

    # XY projection
    scatter = axes[0].scatter(x_real, y_real, c=dens_real, s=1, cmap='viridis')
    if ghost_particles:
        axes[0].scatter(x_ghost, y_ghost, c='red', s=1, alpha=0.3, label='Ghost')
    axes[0].set_xlabel('X [pc]')
    axes[0].set_ylabel('Y [pc]')
    axes[0].set_title('XY Projection')
    axes[0].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[0], label='Density')

    # XZ projection
    scatter = axes[1].scatter(x_real, z_real, c=dens_real, s=1, cmap='viridis')
    axes[1].set_xlabel('X [pc]')
    axes[1].set_ylabel('Z [pc]')
    axes[1].set_title('XZ Projection')
    axes[1].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[1], label='Density')

    # YZ projection
    scatter = axes[2].scatter(y_real, z_real, c=dens_real, s=1, cmap='viridis')
    axes[2].set_xlabel('Y [pc]')
    axes[2].set_ylabel('Z [pc]')
    axes[2].set_title('YZ Projection')
    axes[2].set_aspect('equal')
    plt.colorbar(scatter, ax=axes[2], label='Density')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'particle_distribution.png'), dpi=150)
    plt.close()

    # Print summary
    print("\n" + "="*70)
    print("Isothermal Relaxation Monitoring Summary")
    print("="*70)
    print(f"Snapshots analyzed: {len(steps)}")
    print(f"Latest step: {steps[-1]} / {data['relax_total']}")
    print(f"Progress: {100*steps[-1]/data['relax_total']:.1f}%")
    print(f"\nError statistics (latest snapshot):")
    print(f"  Mean error:  {mean_errors[-1]:.2f}%")
    print(f"  RMS error:   {rms_errors[-1]:.2f}%")
    print(f"  Max error:   {max_errors[-1]:.2f}%")
    print(f"\nPlots saved to: {output_dir}")
    print("="*70)

    return {
        'steps': steps,
        'mean_errors': mean_errors,
        'rms_errors': rms_errors,
        'max_errors': max_errors
    }

def main():
    # Profile parameters (from config)
    rho_center = 5758.46  # M_sun/pc^3
    r_c = 0.0852291       # pc
    R_cloud = 0.170458    # pc

    # Find snapshots
    snapshot_pattern = os.path.join(RESULTS_DIR, "snapshot_*.csv")
    snapshots = sorted(glob.glob(snapshot_pattern))

    if not snapshots:
        print(f"No snapshots found in {RESULTS_DIR}")
        print("Waiting for simulation to produce output...")
        return

    print(f"Found {len(snapshots)} snapshots")

    # Create diagnostic plots
    stats = create_diagnostic_plots(snapshots, rho_center, r_c, R_cloud, PLOT_DIR)

    return stats

if __name__ == "__main__":
    main()
