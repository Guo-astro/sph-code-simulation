#!/usr/bin/env python3
"""Visualize isothermal BE relaxation from CSV snapshots."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob
import os
from pathlib import Path
import imageio.v2 as imageio

def read_csv_snapshot(filepath):
    """Read CSV snapshot, skipping comment lines."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find header line (first non-comment line)
    header_idx = 0
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            header_idx = i
            break

    # Parse header
    header = lines[header_idx].strip().split(',')

    # Parse data
    data = []
    for line in lines[header_idx + 1:]:
        if line.strip() and not line.startswith('#'):
            data.append([float(x) for x in line.strip().split(',')])

    data = np.array(data)
    return {col: data[:, i] for i, col in enumerate(header)}

def main():
    results_dir = Path("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/phase1_compact")
    output_dir = results_dir / "visualization"
    output_dir.mkdir(exist_ok=True)

    # Find all CSV snapshots
    snapshots = sorted(glob.glob(str(results_dir / "snapshot_*.csv")))
    print(f"Found {len(snapshots)} CSV snapshots")

    if not snapshots:
        print("No snapshots found!")
        return

    frames = []

    for i, snap_path in enumerate(snapshots):
        print(f"Processing {os.path.basename(snap_path)}...")

        data = read_csv_snapshot(snap_path)
        x = data['pos_x']
        y = data['pos_y']
        z = data['pos_z']
        dens = data['dens']

        # Calculate radius
        r = np.sqrt(x**2 + y**2 + z**2)

        # Create figure with 2 panels
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Panel 1: XY projection colored by density
        ax1 = axes[0]
        scatter = ax1.scatter(x, y, c=dens, s=0.5, cmap='viridis',
                             norm=LogNorm(vmin=max(dens.min(), 0.01), vmax=dens.max()),
                             alpha=0.7)
        ax1.set_xlabel('X [code units]')
        ax1.set_ylabel('Y [code units]')
        ax1.set_title(f'XY Projection - Step {i}')
        ax1.set_aspect('equal')
        ax1.set_xlim(-1.2, 1.2)
        ax1.set_ylim(-1.2, 1.2)
        plt.colorbar(scatter, ax=ax1, label='Density')

        # Panel 2: Radial density profile
        ax2 = axes[1]
        # Bin particles by radius
        r_bins = np.linspace(0, 1.1, 50)
        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        rho_mean = []
        rho_std = []
        for j in range(len(r_bins) - 1):
            mask = (r >= r_bins[j]) & (r < r_bins[j+1])
            if mask.sum() > 0:
                rho_mean.append(np.mean(dens[mask]))
                rho_std.append(np.std(dens[mask]))
            else:
                rho_mean.append(np.nan)
                rho_std.append(np.nan)

        ax2.plot(r_centers, rho_mean, 'b-', linewidth=2, label='SPH mean')
        ax2.fill_between(r_centers,
                        np.array(rho_mean) - np.array(rho_std),
                        np.array(rho_mean) + np.array(rho_std),
                        alpha=0.3, color='blue', label='±1σ')
        ax2.set_xlabel('Radius [code units]')
        ax2.set_ylabel('Density [code units]')
        ax2.set_title(f'Radial Density Profile - Step {i}')
        ax2.set_yscale('log')
        ax2.set_xlim(0, 1.1)
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.suptitle(f'Isothermal Bonnor-Ebert Relaxation - Snapshot {i:04d}', fontsize=14)
        plt.tight_layout()

        # Save frame
        frame_path = output_dir / f"frame_{i:04d}.png"
        plt.savefig(frame_path, dpi=100)
        plt.close()
        frames.append(str(frame_path))

    # Create GIF
    print("Creating GIF...")
    gif_path = output_dir / "phase1_compact_relaxation.gif"
    images = []
    target_shape = None
    for f in frames:
        img = imageio.imread(f)
        if target_shape is None:
            target_shape = img.shape
        # Pad or crop to match target shape
        if img.shape != target_shape:
            # Create canvas with target shape
            canvas = np.zeros(target_shape, dtype=img.dtype)
            h = min(img.shape[0], target_shape[0])
            w = min(img.shape[1], target_shape[1])
            canvas[:h, :w] = img[:h, :w]
            img = canvas
        images.append(img)
    imageio.mimsave(str(gif_path), images, duration=0.5, loop=0)
    print(f"GIF saved: {gif_path}")

    # Clean up individual frames (optional - keep them for now)
    # for f in frames:
    #     os.remove(f)

if __name__ == "__main__":
    main()
