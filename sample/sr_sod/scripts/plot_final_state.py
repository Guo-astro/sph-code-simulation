#!/usr/bin/env python3
"""Plot final state of SR Sod shock tube simulation"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

def read_snapshot(filepath):
    """Read CSV snapshot file"""
    data = []
    header = None

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            if line.startswith('id,'):
                header = line.strip().split(',')
                continue
            data.append(line.strip().split(','))

    if not data or header is None:
        return None

    # Convert to structured array
    particles = []
    for row in data:
        p = {header[i]: float(row[i]) if i > 0 else int(row[i])
             for i in range(len(header))}
        particles.append(p)

    return particles

def main():
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path(__file__).parent.parent / 'results' / 'sharp'

    # Find last snapshot
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    if not snapshots:
        print("No snapshots found!")
        sys.exit(1)

    final_snap = snapshots[-1]
    print(f"Reading final snapshot: {final_snap.name}")

    data = read_snapshot(final_snap)
    if data is None:
        print("Failed to read snapshot")
        sys.exit(1)

    # Extract arrays
    x = np.array([p['pos_x'] for p in data])
    dens = np.array([p['dens'] for p in data])
    pres = np.array([p['pres'] for p in data])
    vel = np.array([p['vel_x'] for p in data])

    # Sort by position
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    dens = dens[sort_idx]
    pres = pres[sort_idx]
    vel = vel[sort_idx]

    # Extract time from filename
    snap_num = int(final_snap.stem.split('_')[-1])
    time = snap_num * 0.002  # dt from config

    # Create plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 8))
    title = f'SR Sod Shock Tube - Final State\nPons et al. (2000) Shock Speed Formula\nt = {time:.2f}'
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Density
    axes[0].plot(x, dens, 'b.', markersize=2, alpha=0.6)
    axes[0].set_ylabel('Density ρ', fontsize=11)
    axes[0].set_ylim(0, 1.2)
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    axes[0].axhline(y=0.125, color='k', linestyle='--', alpha=0.3)
    axes[0].axvline(x=0, color='k', linestyle=':', alpha=0.5)

    # Velocity
    axes[1].plot(x, vel, 'r.', markersize=2, alpha=0.6)
    axes[1].set_ylabel('Velocity v/c', fontsize=11)
    axes[1].set_ylim(-0.005, 0.004)
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    axes[1].axvline(x=0, color='k', linestyle=':', alpha=0.5)

    # Pressure
    axes[2].plot(x, pres, 'g.', markersize=2, alpha=0.6)
    axes[2].set_ylabel('Pressure P', fontsize=11)
    axes[2].set_xlabel('Position x', fontsize=11)
    axes[2].set_ylim(0, 1.1)
    axes[2].grid(True, alpha=0.3)
    axes[2].axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    axes[2].axhline(y=0.1, color='k', linestyle='--', alpha=0.3)
    axes[2].axvline(x=0, color='k', linestyle=':', alpha=0.5)

    for ax in axes:
        ax.set_xlim(-0.6, 0.6)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save
    output_file = results_dir / 'sr_sod_final_pons_formula.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")

    # Check for oscillations in left region
    left_region = (x < -0.1)
    if np.any(left_region):
        left_dens = dens[left_region]
        left_pres = pres[left_region]
        dens_std = np.std(left_dens)
        pres_std = np.std(pres[left_region])
        print(f"\nLeft region statistics (x < -0.1):")
        print(f"  Density: mean={np.mean(left_dens):.6f}, std={dens_std:.6f}")
        print(f"  Pressure: mean={np.mean(left_pres):.6f}, std={pres_std:.6f}")
        if dens_std < 0.001 and pres_std < 0.001:
            print("  ✓ NO significant oscillations detected!")
        else:
            print("  ⚠ Oscillations still present")

if __name__ == '__main__':
    main()
