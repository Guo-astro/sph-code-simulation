#!/usr/bin/env python3
"""
Create GIF animation of isothermal cloud relaxation.
Shows XY projection colored by density and radial density profile.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob
import os
from pathlib import Path

# Try to import imageio for GIF creation
try:
    import imageio.v2 as imageio
    HAS_IMAGEIO = True
except ImportError:
    try:
        import imageio
        HAS_IMAGEIO = True
    except ImportError:
        HAS_IMAGEIO = False
        print("Warning: imageio not found, will save individual frames only")

def load_snapshot(filepath):
    """Load snapshot CSV file."""
    # Find header line
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('id,'):
                header_line = i
                break

    # Load data
    data = np.genfromtxt(filepath, delimiter=',', skip_header=header_line, names=True)
    return data

def get_step_from_filename(filepath):
    """Extract step number from filename."""
    basename = os.path.basename(filepath)
    # snapshot_0001.csv -> 0001
    num = basename.replace('snapshot_', '').replace('.csv', '')
    return int(num) * 2000  # Each snapshot is 2000 steps apart

def create_frame(data, step, output_dir, frame_idx, r_c=0.0852291, rho_center=5758.37):
    """Create a single frame showing particle distribution and density profile."""

    # Filter out ghost particles
    mask = data['is_ghost'] == 0
    x = data['pos_x'][mask]
    y = data['pos_y'][mask]
    z = data['pos_z'][mask]
    dens = data['dens'][mask]

    # Calculate radial distance
    r = np.sqrt(x**2 + y**2 + z**2)

    # Analytical profile
    r_analytical = np.linspace(0.001, 0.2, 100)
    rho_analytical = rho_center / (1 + (r_analytical / r_c)**2)

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: XY projection colored by density
    ax1 = axes[0]
    scatter = ax1.scatter(x, y, c=dens, s=1, cmap='viridis',
                         norm=LogNorm(vmin=1000, vmax=10000))
    ax1.set_xlim(-0.25, 0.25)
    ax1.set_ylim(-0.25, 0.25)
    ax1.set_xlabel('x [pc]', fontsize=12)
    ax1.set_ylabel('y [pc]', fontsize=12)
    ax1.set_title(f'XY Projection (Step {step:,})', fontsize=14)
    ax1.set_aspect('equal')
    cbar = plt.colorbar(scatter, ax=ax1, label='Density [code units]')

    # Draw cloud boundary circle
    theta = np.linspace(0, 2*np.pi, 100)
    R_cloud = 0.170458
    ax1.plot(R_cloud * np.cos(theta), R_cloud * np.sin(theta), 'r--',
             linewidth=1.5, label=f'R_cloud = {R_cloud:.3f} pc')
    ax1.legend(loc='upper right', fontsize=10)

    # Right: Radial density profile
    ax2 = axes[1]

    # Bin particles by radius
    r_bins = np.linspace(0, 0.2, 30)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    rho_binned = []
    rho_std = []

    for i in range(len(r_bins) - 1):
        mask_bin = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask_bin) > 0:
            rho_binned.append(np.mean(dens[mask_bin]))
            rho_std.append(np.std(dens[mask_bin]))
        else:
            rho_binned.append(np.nan)
            rho_std.append(np.nan)

    rho_binned = np.array(rho_binned)
    rho_std = np.array(rho_std)

    # Plot analytical profile
    ax2.semilogy(r_analytical, rho_analytical, 'b-', linewidth=2,
                 label='Analytical: ρ/(1+(r/r_c)²)')

    # Plot SPH profile with error bars
    ax2.errorbar(r_centers, rho_binned, yerr=rho_std, fmt='ro',
                 markersize=5, capsize=3, label='SPH particles')

    ax2.set_xlim(0, 0.2)
    ax2.set_ylim(500, 15000)
    ax2.set_xlabel('Radius r [pc]', fontsize=12)
    ax2.set_ylabel('Density ρ [code units]', fontsize=12)
    ax2.set_title(f'Radial Density Profile (Step {step:,})', fontsize=14)
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)

    # Add text annotations
    n_particles = len(x)
    max_acc = np.max(np.sqrt(data['acc_x']**2 + data['acc_y']**2 + data['acc_z']**2))

    # Compute mean error
    rho_expected = rho_center / (1 + (r / r_c)**2)
    error = np.abs(dens - rho_expected) / rho_expected
    mean_error = np.mean(error) * 100

    info_text = f'Particles: {n_particles:,}\nMax Accel: {max_acc:.2f}\nMean Error: {mean_error:.1f}%'
    ax2.text(0.95, 0.05, info_text, transform=ax2.transAxes,
             fontsize=10, verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    # Save frame
    frame_path = output_dir / f'frame_{frame_idx:04d}.png'
    plt.savefig(frame_path, dpi=120, bbox_inches='tight')
    plt.close()

    return str(frame_path)

def main():
    # Paths
    results_dir = Path('/Users/guo/Downloads/sph-simulators/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/hvcc_10k/isothermal')
    output_dir = results_dir / 'plots'
    output_dir.mkdir(exist_ok=True)

    # Find all snapshots
    snapshot_files = sorted(glob.glob(str(results_dir / 'snapshot_*.csv')))

    if not snapshot_files:
        print("No snapshot files found!")
        return

    print(f"Found {len(snapshot_files)} snapshots")

    # Create frames
    frame_paths = []
    for i, filepath in enumerate(snapshot_files):
        step = get_step_from_filename(filepath)
        print(f"Processing snapshot {i+1}/{len(snapshot_files)}: step {step}")

        data = load_snapshot(filepath)
        frame_path = create_frame(data, step, output_dir, i)
        frame_paths.append(frame_path)

    # Create GIF
    if HAS_IMAGEIO and len(frame_paths) > 1:
        gif_path = results_dir / 'relaxation_animation.gif'
        print(f"\nCreating GIF: {gif_path}")

        images = [imageio.imread(fp) for fp in frame_paths]
        # Duration per frame in seconds (slower for relaxation visualization)
        imageio.mimsave(str(gif_path), images, duration=1.0, loop=0)

        print(f"GIF saved to: {gif_path}")
        print(f"Frames saved to: {output_dir}")
    else:
        print(f"\nFrames saved to: {output_dir}")
        if not HAS_IMAGEIO:
            print("Install imageio to create GIF: pip install imageio")

    return str(results_dir / 'relaxation_animation.gif') if HAS_IMAGEIO else None

if __name__ == '__main__':
    main()
