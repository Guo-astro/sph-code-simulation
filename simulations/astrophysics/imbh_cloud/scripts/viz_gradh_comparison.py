#!/usr/bin/env python3
"""
Grad-h Comparison Visualization
================================

Consolidated visualization script for comparing GSPH methods with and without
grad-h correction. Generates both static plots and animations.

Modes:
------
1. Static plots: Energy evolution, density profiles, stability metrics
2. Animation: Time evolution of density profiles for all methods side-by-side

Methods Compared:
-----------------
- GSPH + grad-h: Full grad-h correction (baseline, expected stable)
- GSPH - grad-h: No grad-h correction (may collapse)
- GSPH + C_smooth=2: Smoother h variation, no grad-h

Usage:
------
    # Static comparison plots
    python3 viz_gradh_comparison.py <results_base_dir> <output_dir> --mode plots

    # Animation
    python3 viz_gradh_comparison.py <results_base_dir> <output_dir> --mode animate

    # Both
    python3 viz_gradh_comparison.py <results_base_dir> <output_dir> --mode all

Example:
--------
    python3 viz_gradh_comparison.py results/gradh_test animations/gradh_test --mode all

Input Directory Structure:
--------------------------
    results_base_dir/
    ├── gsph_gradh/
    │   ├── energy.dat
    │   └── snapshot_*.csv
    ├── gsph_nogradh/
    │   ├── energy.dat
    │   └── snapshot_*.csv
    └── gsph_csmooth2/
        ├── energy.dat
        └── snapshot_*.csv
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Try to import visualization libraries
try:
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available, visualization disabled")

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


# =============================================================================
# CONFIGURATION
# =============================================================================

METHODS = {
    'gsph_gradh': {
        'label': 'GSPH + grad-h',
        'color': '#0066CC',
        'linestyle': '-',
        'expected': 'stable'
    },
    'gsph_nogradh': {
        'label': 'GSPH - grad-h',
        'color': '#CC0000',
        'linestyle': '--',
        'expected': 'may collapse'
    },
    'gsph_csmooth2': {
        'label': 'GSPH + C_smooth=2',
        'color': '#009933',
        'linestyle': '-.',
        'expected': 'test smoother h'
    },
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_energy_data(results_dir: Path, method: str) -> Optional[np.ndarray]:
    """Load energy.dat file for a method."""
    energy_file = results_dir / method / 'energy.dat'
    if not energy_file.exists():
        return None

    try:
        data = np.loadtxt(energy_file, skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"  Warning: Failed to load {energy_file}: {e}")
        return None


def load_snapshot(filepath: Path) -> Optional[np.ndarray]:
    """Load a single snapshot CSV file."""
    if not filepath.exists():
        return None

    try:
        # Find header line
        with open(filepath, 'r') as f:
            skiprows = 0
            for line in f:
                if line.strip().startswith('id,'):
                    break
                skiprows += 1

        if HAS_PANDAS:
            df = pd.read_csv(filepath, skiprows=skiprows)
            return df
        else:
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows + 1)
            return data
    except Exception as e:
        print(f"  Warning: Failed to load {filepath}: {e}")
        return None


def get_snapshot_count(results_dir: Path, method: str) -> int:
    """Get number of snapshots for a method."""
    pattern = results_dir / method / 'snapshot_*.csv'
    return len(list(results_dir.glob(f"{method}/snapshot_*.csv")))


# =============================================================================
# STATIC PLOTS
# =============================================================================

def create_energy_comparison(results_dir: Path, output_dir: Path) -> None:
    """Create energy comparison plots."""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    has_data = False

    for method, props in METHODS.items():
        data = load_energy_data(results_dir, method)
        if data is None or len(data) == 0:
            print(f"  Skipping {method}: no energy data")
            continue

        has_data = True
        t = data[:, 0]
        E_kin = data[:, 1]
        E_th = data[:, 2]
        E_pot = data[:, 3]
        E_total = data[:, 4] if data.shape[1] > 4 else E_kin + E_th + E_pot

        # Total energy
        axes[0, 0].plot(t, E_total, color=props['color'], linestyle=props['linestyle'],
                       label=props['label'], linewidth=2)

        # Kinetic energy
        axes[0, 1].plot(t, E_kin, color=props['color'], linestyle=props['linestyle'],
                       label=props['label'], linewidth=2)

        # Energy drift
        E0 = E_total[0]
        drift = (E_total - E0) / np.abs(E0) * 100
        axes[1, 0].plot(t, drift, color=props['color'], linestyle=props['linestyle'],
                       label=props['label'], linewidth=2)

        # Potential energy
        axes[1, 1].plot(t, E_pot, color=props['color'], linestyle=props['linestyle'],
                       label=props['label'], linewidth=2)

        print(f"  Loaded {method}: {len(t)} time steps, E_drift_final = {drift[-1]:.2f}%")

    if not has_data:
        print("  No energy data found for any method")
        plt.close()
        return

    axes[0, 0].set_xlabel('Time')
    axes[0, 0].set_ylabel('Total Energy')
    axes[0, 0].set_title('Total Energy Evolution')
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].set_xlabel('Time')
    axes[0, 1].set_ylabel('Kinetic Energy')
    axes[0, 1].set_title('Kinetic Energy (Stability Indicator)')
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].set_xlabel('Time')
    axes[1, 0].set_ylabel('Energy Drift (%)')
    axes[1, 0].set_title('Energy Conservation Error')
    axes[1, 0].legend(fontsize=9)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axhline(y=0, color='gray', linestyle=':', alpha=0.5)

    axes[1, 1].set_xlabel('Time')
    axes[1, 1].set_ylabel('Potential Energy')
    axes[1, 1].set_title('Gravitational Potential Energy')
    axes[1, 1].legend(fontsize=9)
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle('Grad-h Comparison: 3 GSPH Variants\n(kernel gravity, self-gravitating sphere)',
                fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    output_path = output_dir / 'gradh_energy_comparison.png'
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved: {output_path}")


def create_density_comparison(results_dir: Path, output_dir: Path, frame_idx: int = -1) -> None:
    """Create density profile comparison at a specific frame (default: final)."""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (method, props) in enumerate(METHODS.items()):
        ax = axes[i]

        # Find snapshots
        snapshots = sorted((results_dir / method).glob('snapshot_*.csv'))
        if not snapshots:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(props['label'])
            continue

        snap_path = snapshots[frame_idx] if frame_idx < len(snapshots) else snapshots[-1]
        snap = load_snapshot(snap_path)

        if snap is None:
            ax.text(0.5, 0.5, 'Failed to load', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(props['label'])
            continue

        # Compute radial profile
        if HAS_PANDAS:
            r = np.sqrt(snap['pos_x']**2 + snap['pos_y']**2 + snap['pos_z']**2)
            rho = snap['dens']
        else:
            r = np.sqrt(snap[:, 1]**2 + snap[:, 2]**2 + snap[:, 3]**2)
            rho = snap[:, 11]

        # Sort by radius
        sort_idx = np.argsort(r)
        r_sorted = np.array(r)[sort_idx]
        rho_sorted = np.array(rho)[sort_idx]

        ax.scatter(r_sorted, rho_sorted, s=1, alpha=0.3, color=props['color'])
        ax.set_xlabel('Radius r')
        ax.set_ylabel('Density ρ')
        ax.set_title(f"{props['label']}\n({props['expected']})")
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 2])
        ax.grid(True, alpha=0.3)

    plt.suptitle(f'Density Profiles at Frame {frame_idx}', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.92])

    output_path = output_dir / 'gradh_density_comparison.png'
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved: {output_path}")


# =============================================================================
# ANIMATION
# =============================================================================

def create_comparison_animation(
    results_dir: Path,
    output_dir: Path,
    fps: int = 10,
    max_frames: int = 100
) -> None:
    """Create animated comparison of all methods."""
    if not HAS_MATPLOTLIB:
        return

    # Count snapshots for each method
    snap_counts = {}
    for method in METHODS:
        snap_counts[method] = get_snapshot_count(results_dir, method)
        print(f"  {method}: {snap_counts[method]} snapshots")

    max_snaps = max(snap_counts.values()) if snap_counts else 0
    if max_snaps == 0:
        print("  No snapshots found, skipping animation")
        return

    n_frames = min(max_snaps, max_frames)

    # Set up figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Initialize plots
    lines = {}
    for i, (method, props) in enumerate(METHODS.items()):
        ax = axes[i]
        line, = ax.plot([], [], 'o', markersize=1, alpha=0.5, color=props['color'])
        lines[method] = (line, ax)

        ax.set_xlabel('Radius r')
        ax.set_ylabel('Density ρ')
        ax.set_title(f"{props['label']}\n({props['expected']})")
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 2])
        ax.grid(True, alpha=0.3)

    time_text = fig.suptitle('t = 0.0', fontsize=14, fontweight='bold')

    def update(frame):
        for method, (line, ax) in lines.items():
            snap_path = results_dir / method / f'snapshot_{frame:04d}.csv'
            snap = load_snapshot(snap_path)

            if snap is not None:
                if HAS_PANDAS:
                    r = np.sqrt(snap['pos_x']**2 + snap['pos_y']**2 + snap['pos_z']**2)
                    rho = snap['dens']
                else:
                    r = np.sqrt(snap[:, 1]**2 + snap[:, 2]**2 + snap[:, 3]**2)
                    rho = snap[:, 11]

                line.set_data(r, rho)

        time_text.set_text(f'Grad-h Comparison: Frame {frame}')
        return list(l[0] for l in lines.values()) + [time_text]

    print(f"  Creating animation with {n_frames} frames...")
    anim = FuncAnimation(fig, update, frames=n_frames, interval=100, blit=False)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    output_path = output_dir / 'gradh_comparison.gif'
    anim.save(output_path, writer=PillowWriter(fps=fps))
    plt.close()

    print(f"  Saved: {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Grad-h comparison visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('results_dir', type=Path,
                        help='Base directory containing method subdirectories')
    parser.add_argument('output_dir', type=Path,
                        help='Output directory for plots and animations')
    parser.add_argument('--mode', choices=['plots', 'animate', 'all'], default='all',
                        help='Visualization mode (default: all)')
    parser.add_argument('--fps', type=int, default=10,
                        help='Animation frames per second (default: 10)')
    parser.add_argument('--max-frames', type=int, default=100,
                        help='Maximum animation frames (default: 100)')

    args = parser.parse_args()

    if not HAS_MATPLOTLIB:
        print("Error: matplotlib is required for visualization")
        print("  Install with: pip install matplotlib")
        return 1

    if not args.results_dir.exists():
        print(f"Error: Results directory not found: {args.results_dir}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGrad-h Comparison Visualization")
    print(f"================================")
    print(f"Results:  {args.results_dir}")
    print(f"Output:   {args.output_dir}")
    print(f"Mode:     {args.mode}")

    if args.mode in ['plots', 'all']:
        print("\nGenerating static plots...")
        create_energy_comparison(args.results_dir, args.output_dir)
        create_density_comparison(args.results_dir, args.output_dir)

    if args.mode in ['animate', 'all']:
        print("\nGenerating animation...")
        create_comparison_animation(
            args.results_dir,
            args.output_dir,
            fps=args.fps,
            max_frames=args.max_frames
        )

    print("\nDone!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
