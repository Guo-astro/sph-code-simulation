#!/usr/bin/env python3
"""
Compare GR-GSPH (Minkowski) with SR-GSPH Results

This script verifies that GR-GSPH reduces exactly to SR-GSPH in flat spacetime.
Both simulations use the Rosswog Test 1 initial conditions:
  Left:  (n, v, P) = (10, 0, 40/3)
  Right: (n, v, P) = (1, 0, 1e-6)

The results should be IDENTICAL (within numerical precision).
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add srrp solver to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / 'docs'))

try:
    from srrp.Solver import Solver as SRRPSolver
    from srrp.State import State
    HAS_SRRP = True

    # Pre-compute wavefan for Rosswog Test 1
    _solver = SRRPSolver()
    _stateL = State(rho=10.0, vx=0.0, vt=0.0, pressure=40.0/3.0)
    _stateR = State(rho=1.0, vx=0.0, vt=0.0, pressure=1e-6)
    _wavefan = _solver.solve(_stateL, _stateR, 5.0/3.0)
    print(f"SRRP solution type: {_solver.solution_type}")
except ImportError as e:
    print(f"Warning: Could not import SRRP solver: {e}")
    HAS_SRRP = False
    _wavefan = None


def load_csv(filename):
    """Load CSV snapshot with header comments."""
    import pandas as pd
    data = pd.read_csv(filename, comment='#')

    # Extract time from header
    time = 0.0
    with open(filename, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                try:
                    time = float(line.split(':')[1].strip())
                except:
                    pass
                break

    return data, time


def compute_exact_solution(x, t):
    """Compute exact SR solution using SRRP solver."""
    if t <= 0 or _wavefan is None:
        n = np.where(x < 0, 10.0, 1.0)
        v = np.zeros_like(x)
        P = np.where(x < 0, 40.0/3.0, 1e-6)
        return n, v, P

    xi = x / t
    state = _wavefan.getState(xi)
    return state.rho, state.vx, state.pressure


def find_matching_snapshot(grgsph_dir, srgsph_dir, target_time=0.35):
    """Find snapshots closest to target time from both simulations."""
    grgsph_files = sorted(glob.glob(os.path.join(grgsph_dir, 'snapshot_*.csv')))
    srgsph_files = sorted(glob.glob(os.path.join(srgsph_dir, 'snapshot_*.csv')))

    if not grgsph_files or not srgsph_files:
        return None, None

    # Find closest to target time
    def find_closest(files, target):
        best_file = None
        best_diff = float('inf')
        for f in files:
            _, t = load_csv(f)
            diff = abs(t - target)
            if diff < best_diff:
                best_diff = diff
                best_file = f
        return best_file

    return find_closest(grgsph_files, target_time), find_closest(srgsph_files, target_time)


def create_comparison_plot(grgsph_dir, srgsph_dir, output_file):
    """Create side-by-side comparison of GR-GSPH and SR-GSPH."""

    grgsph_file, srgsph_file = find_matching_snapshot(grgsph_dir, srgsph_dir, 0.35)

    if grgsph_file is None or srgsph_file is None:
        print("Could not find matching snapshots!")
        return

    # Load data
    grgsph_data, grgsph_time = load_csv(grgsph_file)
    srgsph_data, srgsph_time = load_csv(srgsph_file)

    print(f"GR-GSPH snapshot: t = {grgsph_time:.4f}")
    print(f"SR-GSPH snapshot: t = {srgsph_time:.4f}")

    # Filter out ghosts
    if 'is_ghost' in grgsph_data.columns:
        gr_mask = grgsph_data['is_ghost'] == 0
    else:
        gr_mask = np.ones(len(grgsph_data), dtype=bool)

    if 'is_ghost' in srgsph_data.columns:
        sr_mask = srgsph_data['is_ghost'] == 0
    else:
        sr_mask = np.ones(len(srgsph_data), dtype=bool)

    gr_x = grgsph_data['pos_x'].values[gr_mask]
    gr_n = grgsph_data['dens'].values[gr_mask]
    gr_v = grgsph_data['vel_x'].values[gr_mask]
    gr_P = grgsph_data['pres'].values[gr_mask]

    sr_x = srgsph_data['pos_x'].values[sr_mask]
    sr_n = srgsph_data['dens'].values[sr_mask]
    sr_v = srgsph_data['vel_x'].values[sr_mask]
    sr_P = srgsph_data['pres'].values[sr_mask]

    # Exact solution
    x_exact = np.linspace(-0.5, 0.5, 500)
    n_exact, v_exact, P_exact = compute_exact_solution(x_exact, grgsph_time)

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    fig.suptitle(f'GR-GSPH (Minkowski) vs SR-GSPH Comparison\n'
                 f'Rosswog Test 1, t = {grgsph_time:.4f}', fontsize=14, fontweight='bold')

    # Row 1: Individual results
    # GR-GSPH Density
    axes[0, 0].scatter(gr_x, gr_n, s=2, alpha=0.6, c='blue', label='GR-GSPH')
    axes[0, 0].plot(x_exact, n_exact, 'r-', lw=1.5, label='Exact SR')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('Density n')
    axes[0, 0].set_xlim(-0.55, 0.55)
    axes[0, 0].set_ylim(0, 12)
    axes[0, 0].legend(loc='upper right')
    axes[0, 0].set_title('GR-GSPH (Minkowski) Density')
    axes[0, 0].grid(True, alpha=0.3)

    # SR-GSPH Density
    axes[0, 1].scatter(sr_x, sr_n, s=2, alpha=0.6, c='green', label='SR-GSPH')
    axes[0, 1].plot(x_exact, n_exact, 'r-', lw=1.5, label='Exact SR')
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('Density n')
    axes[0, 1].set_xlim(-0.55, 0.55)
    axes[0, 1].set_ylim(0, 12)
    axes[0, 1].legend(loc='upper right')
    axes[0, 1].set_title('SR-GSPH Density')
    axes[0, 1].grid(True, alpha=0.3)

    # Overlay comparison
    axes[0, 2].scatter(gr_x, gr_n, s=3, alpha=0.6, c='blue', label='GR-GSPH')
    axes[0, 2].scatter(sr_x, sr_n, s=3, alpha=0.6, c='green', label='SR-GSPH', marker='x')
    axes[0, 2].plot(x_exact, n_exact, 'r-', lw=2, label='Exact', alpha=0.8)
    axes[0, 2].set_xlabel('x')
    axes[0, 2].set_ylabel('Density n')
    axes[0, 2].set_xlim(-0.55, 0.55)
    axes[0, 2].set_ylim(0, 12)
    axes[0, 2].legend(loc='upper right')
    axes[0, 2].set_title('OVERLAY: Both Should Be Identical')
    axes[0, 2].grid(True, alpha=0.3)

    # Row 2: Velocity and Pressure comparison
    # Velocity
    axes[1, 0].scatter(gr_x, gr_v, s=2, alpha=0.6, c='blue', label='GR-GSPH')
    axes[1, 0].scatter(sr_x, sr_v, s=2, alpha=0.6, c='green', label='SR-GSPH', marker='x')
    axes[1, 0].plot(x_exact, v_exact, 'r-', lw=2, label='Exact', alpha=0.8)
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('Velocity v')
    axes[1, 0].set_xlim(-0.55, 0.55)
    axes[1, 0].set_ylim(-0.1, 0.9)
    axes[1, 0].legend(loc='upper left')
    axes[1, 0].set_title('Velocity Comparison')
    axes[1, 0].grid(True, alpha=0.3)

    # Pressure
    axes[1, 1].scatter(gr_x, gr_P, s=2, alpha=0.6, c='blue', label='GR-GSPH')
    axes[1, 1].scatter(sr_x, sr_P, s=2, alpha=0.6, c='green', label='SR-GSPH', marker='x')
    axes[1, 1].plot(x_exact, P_exact, 'r-', lw=2, label='Exact', alpha=0.8)
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('Pressure P')
    axes[1, 1].set_xlim(-0.55, 0.55)
    axes[1, 1].set_ylim(-0.5, 16)
    axes[1, 1].legend(loc='upper right')
    axes[1, 1].set_title('Pressure Comparison')
    axes[1, 1].grid(True, alpha=0.3)

    # Difference plot
    # Interpolate SR-GSPH onto GR-GSPH positions for comparison
    from scipy.interpolate import interp1d

    sr_sort = np.argsort(sr_x)
    sr_x_sorted = sr_x[sr_sort]
    sr_n_sorted = sr_n[sr_sort]

    gr_sort = np.argsort(gr_x)
    gr_x_sorted = gr_x[gr_sort]
    gr_n_sorted = gr_n[gr_sort]

    # Only compare in overlap region
    x_min = max(gr_x_sorted.min(), sr_x_sorted.min())
    x_max = min(gr_x_sorted.max(), sr_x_sorted.max())

    x_compare = np.linspace(x_min + 0.01, x_max - 0.01, 100)

    try:
        f_gr = interp1d(gr_x_sorted, gr_n_sorted, kind='linear', bounds_error=True)
        f_sr = interp1d(sr_x_sorted, sr_n_sorted, kind='linear', bounds_error=True)

        n_gr_interp = f_gr(x_compare)
        n_sr_interp = f_sr(x_compare)

        diff = n_gr_interp - n_sr_interp
        rel_diff = diff / np.maximum(np.abs(n_gr_interp), 1e-10) * 100

        axes[1, 2].plot(x_compare, rel_diff, 'b-', lw=1.5)
        axes[1, 2].axhline(0, color='k', linestyle='--', lw=0.5)
        axes[1, 2].set_xlabel('x')
        axes[1, 2].set_ylabel('Relative Difference (%)')
        axes[1, 2].set_xlim(-0.55, 0.55)
        axes[1, 2].set_title(f'Density Difference (max: {np.max(np.abs(rel_diff)):.2f}%)')
        axes[1, 2].grid(True, alpha=0.3)

        print(f"\nDensity comparison:")
        print(f"  Max relative difference: {np.max(np.abs(rel_diff)):.4f}%")
        print(f"  RMS relative difference: {np.sqrt(np.mean(rel_diff**2)):.4f}%")

    except Exception as e:
        print(f"Could not compute difference: {e}")
        axes[1, 2].text(0.5, 0.5, f'Error: {e}', transform=axes[1, 2].transAxes, ha='center')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison to: {output_file}")
    plt.close()


def main():
    script_dir = Path(__file__).parent
    grgsph_dir = script_dir.parent / 'results' / 'minkowski_shock_tube'
    srgsph_dir = script_dir.parent.parent / 'sr_sod' / 'results' / 'rosswog_test1'
    output_file = script_dir.parent / 'results' / 'grgsph_vs_srgsph_comparison.png'

    print("="*70)
    print("  GR-GSPH (Minkowski) vs SR-GSPH Comparison")
    print("="*70)
    print(f"\nGR-GSPH results: {grgsph_dir}")
    print(f"SR-GSPH results: {srgsph_dir}")
    print(f"Output: {output_file}")
    print()

    if not grgsph_dir.exists():
        print(f"Error: GR-GSPH results not found at {grgsph_dir}")
        return 1

    if not srgsph_dir.exists():
        print(f"Error: SR-GSPH results not found at {srgsph_dir}")
        return 1

    create_comparison_plot(str(grgsph_dir), str(srgsph_dir), str(output_file))
    return 0


if __name__ == '__main__':
    sys.exit(main())
