#!/usr/bin/env -S uv run --with numpy --with h5py --with matplotlib --with imageio python
"""Visualize Lane-Emden relaxation snapshots as GIF"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import imageio.v2 as imageio

# Find all snapshots
results_dir = Path("simulations/astrophysics/imbh_cloud/results/lane_emden_1M/phase1_relaxation")
snapshots = sorted(results_dir.glob("snapshot_*.h5"))
print(f"Found {len(snapshots)} snapshots")

if len(snapshots) == 0:
    print("No snapshots found!")
    exit(1)

frames = []
for snap_path in snapshots:
    with h5py.File(snap_path, 'r') as f:
        p = f['particles']
        x = p['pos_x'][:]
        y = p['pos_y'][:]
        z = p['pos_z'][:]
        dens = p['dens'][:]

    # Get step number from filename
    step = int(snap_path.stem.split('_')[1]) * 50

    # Create figure with 2 panels
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: XY slice (|z| < 0.05)
    mask = np.abs(z) < 0.05
    ax1 = axes[0]
    sc = ax1.scatter(x[mask], y[mask], c=dens[mask], s=0.1, cmap='viridis',
                     norm=LogNorm(vmin=0.01, vmax=2.0))
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title(f'XY slice (|z|<0.05) - Step {step}')
    ax1.set_aspect('equal')
    plt.colorbar(sc, ax=ax1, label='Density')

    # Panel 2: Radial density profile
    r = np.sqrt(x**2 + y**2 + z**2)
    ax2 = axes[1]

    # Bin the data
    r_bins = np.linspace(0, 1.05, 100)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    dens_mean = np.zeros(len(r_centers))
    dens_std = np.zeros(len(r_centers))

    for i in range(len(r_centers)):
        mask_bin = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask_bin) > 0:
            dens_mean[i] = np.mean(dens[mask_bin])
            dens_std[i] = np.std(dens[mask_bin])

    ax2.fill_between(r_centers, dens_mean - dens_std, dens_mean + dens_std,
                     alpha=0.3, color='blue', label='±1σ')
    ax2.plot(r_centers, dens_mean, 'b-', lw=2, label='SPH mean')

    # Lane-Emden analytic (approximate)
    xi_1 = 3.65375
    rho_c = 1.4301
    xi = r_centers * xi_1
    theta = np.maximum(1 - xi**2/6, 0)  # Approximate for small xi
    rho_analytic = rho_c * theta**1.5
    ax2.plot(r_centers, rho_analytic, 'r--', lw=2, label='Analytic (approx)')

    ax2.set_xlabel('Radius')
    ax2.set_ylabel('Density')
    ax2.set_title(f'Radial Density Profile - Step {step}')
    ax2.set_xlim(0, 1.05)
    ax2.set_ylim(0, 2.0)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save frame
    frame_path = f"/tmp/lane_emden_frame_{step:04d}.png"
    plt.savefig(frame_path, dpi=100)
    frames.append(imageio.imread(frame_path))
    plt.close()
    print(f"  Created frame for step {step}")

# Create GIF
gif_path = str(results_dir / "relaxation_progress.gif")
imageio.mimsave(gif_path, frames, duration=1.0, loop=0)
print(f"\nGIF saved to: {gif_path}")
