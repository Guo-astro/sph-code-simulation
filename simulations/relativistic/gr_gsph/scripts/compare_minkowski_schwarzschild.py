#!/usr/bin/env python3
"""
Compare GR-GSPH results to verify correctness:
1. Minkowski metric → should match exact SR solution perfectly
2. Schwarzschild metric → should differ from flat solution (GR effects)

This is the key verification: if Minkowski matches exact SR, our code is correct.
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
except ImportError as e:
    print(f"Error: Could not import SRRP solver: {e}")
    sys.exit(1)

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


def compute_exact_sr_solution(r, t, r_disc=6.0, gamma=5.0/3.0):
    """
    Compute exact SR Riemann solution in flat spacetime.
    """
    # Initial conditions (Rosswog Test 1)
    n_L, v_L, P_L = 10.0, 0.0, 40.0/3.0
    n_R, v_R, P_R = 1.0, 0.0, 1e-6

    if t <= 0:
        n = np.where(r < r_disc, n_L, n_R)
        v = np.where(r < r_disc, v_L, v_R)
        P = np.where(r < r_disc, P_L, P_R)
        return n, v, P

    # Solve Riemann problem
    solver = SRRPSolver()
    stateL = State(rho=n_L, vx=v_L, vt=0.0, pressure=P_L)
    stateR = State(rho=n_R, vx=v_R, vt=0.0, pressure=P_R)
    wavefan = solver.solve(stateL, stateR, gamma)

    # Self-similar variable
    xi = (r - r_disc) / t

    n = np.zeros_like(r)
    v = np.zeros_like(r)
    P = np.zeros_like(r)

    for i, xi_val in enumerate(xi):
        try:
            state = wavefan.getState(xi_val)
            n[i] = state.rho
            v[i] = state.vx
            P[i] = state.pressure
        except:
            n[i] = n_L if r[i] < r_disc else n_R
            v[i] = v_L if r[i] < r_disc else v_R
            P[i] = P_L if r[i] < r_disc else P_R

    return n, v, P


def main():
    script_dir = Path(__file__).parent
    # Try minkowski_shock_tube first, fall back to minkowski_shock
    minkowski_dir = script_dir.parent / 'results' / 'minkowski_shock_tube'
    if not minkowski_dir.exists():
        minkowski_dir = script_dir.parent / 'results' / 'minkowski_shock'
    schwarzschild_dir = script_dir.parent / 'results' / 'schwarzschild_shock'

    print("=" * 70)
    print("  GR-GSPH Verification: Minkowski vs Schwarzschild")
    print("=" * 70)

    if not minkowski_dir.exists():
        print(f"Error: Minkowski results not found at {minkowski_dir}")
        return 1
    if not schwarzschild_dir.exists():
        print(f"Error: Schwarzschild results not found at {schwarzschild_dir}")
        return 1

    # Find matching time snapshots
    mink_files = sorted(glob.glob(str(minkowski_dir / 'snapshot_*.csv')))
    schw_files = sorted(glob.glob(str(schwarzschild_dir / 'snapshot_*.csv')))

    # Use last snapshot for comparison
    mink_data, mink_t = load_csv(mink_files[-1])
    schw_data, schw_t = load_csv(schw_files[-1])

    print(f"Minkowski: {len(mink_files)} snapshots, final t = {mink_t:.3f}")
    print(f"Schwarzschild: {len(schw_files)} snapshots, final t = {schw_t:.3f}")

    # Filter ghosts
    if 'is_ghost' in mink_data.columns:
        mink_mask = mink_data['is_ghost'] == 0
        schw_mask = schw_data['is_ghost'] == 0
    else:
        mink_mask = np.ones(len(mink_data), dtype=bool)
        schw_mask = np.ones(len(schw_data), dtype=bool)

    # Extract data
    r_mink = mink_data['pos_x'].values[mink_mask]
    n_mink = mink_data['dens'].values[mink_mask]
    v_mink = mink_data['vel_x'].values[mink_mask]
    P_mink = mink_data['pres'].values[mink_mask]

    r_schw = schw_data['pos_x'].values[schw_mask]
    n_schw = schw_data['dens'].values[schw_mask]
    v_schw = schw_data['vel_x'].values[schw_mask]
    P_schw = schw_data['pres'].values[schw_mask]

    # Sort by radius
    mink_sort = np.argsort(r_mink)
    schw_sort = np.argsort(r_schw)
    r_mink, n_mink, v_mink, P_mink = r_mink[mink_sort], n_mink[mink_sort], v_mink[mink_sort], P_mink[mink_sort]
    r_schw, n_schw, v_schw, P_schw = r_schw[schw_sort], n_schw[schw_sort], v_schw[schw_sort], P_schw[schw_sort]

    # Compute exact solution at Minkowski time
    r_ana = np.linspace(3.0, 15.0, 500)
    n_exact, v_exact, P_exact = compute_exact_sr_solution(r_ana, mink_t, r_disc=6.0)

    # Create comparison figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Top row: Minkowski vs Exact (verification)
    ax = axes[0, 0]
    ax.scatter(r_mink, n_mink, s=3, alpha=0.6, c='blue', label='GR-GSPH (Minkowski)')
    ax.plot(r_ana, n_exact, 'r-', lw=2, label='Exact SR')
    ax.set_xlabel('r')
    ax.set_ylabel('Density')
    ax.set_xlim(3, 15)
    ax.set_ylim(0, 15)
    ax.legend(fontsize=8)
    ax.set_title(f'Minkowski: Density (t={mink_t:.2f})')
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.scatter(r_mink, v_mink, s=3, alpha=0.6, c='blue', label='GR-GSPH (Minkowski)')
    ax.plot(r_ana, v_exact, 'r-', lw=2, label='Exact SR')
    ax.axhline(0, color='k', linestyle='-', lw=0.5)
    ax.set_xlabel('r')
    ax.set_ylabel('Velocity')
    ax.set_xlim(3, 15)
    ax.set_ylim(-0.5, 1.0)
    ax.legend(fontsize=8)
    ax.set_title('Minkowski: Velocity')
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    ax.scatter(r_mink, P_mink, s=3, alpha=0.6, c='blue', label='GR-GSPH (Minkowski)')
    ax.plot(r_ana, P_exact, 'r-', lw=2, label='Exact SR')
    ax.set_xlabel('r')
    ax.set_ylabel('Pressure')
    ax.set_xlim(3, 15)
    ax.set_ylim(-0.5, 20)
    ax.legend(fontsize=8)
    ax.set_title('Minkowski: Pressure')
    ax.grid(True, alpha=0.3)

    # Bottom row: Schwarzschild vs Minkowski (GR effects)
    ax = axes[1, 0]
    ax.scatter(r_schw, n_schw, s=3, alpha=0.6, c='green', label='Schwarzschild')
    ax.scatter(r_mink, n_mink, s=3, alpha=0.4, c='blue', label='Minkowski')
    ax.plot(r_ana, n_exact, 'r--', lw=1, alpha=0.5, label='Exact SR')
    ax.axvline(2.0, color='black', linestyle='--', lw=1.5, alpha=0.5, label='Horizon')
    ax.set_xlabel('r')
    ax.set_ylabel('Density')
    ax.set_xlim(3, 15)
    ax.set_ylim(0, 15)
    ax.legend(fontsize=8)
    ax.set_title(f'Comparison: Density (t={schw_t:.2f})')
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.scatter(r_schw, v_schw, s=3, alpha=0.6, c='green', label='Schwarzschild')
    ax.scatter(r_mink, v_mink, s=3, alpha=0.4, c='blue', label='Minkowski')
    ax.plot(r_ana, v_exact, 'r--', lw=1, alpha=0.5, label='Exact SR')
    ax.axhline(0, color='k', linestyle='-', lw=0.5)
    ax.axvline(2.0, color='black', linestyle='--', lw=1.5, alpha=0.5)
    ax.set_xlabel('r')
    ax.set_ylabel('Velocity')
    ax.set_xlim(3, 15)
    # Auto-scale y to show full velocity range (GR infall can be very negative)
    # Filter out NaN/Inf values
    v_schw_valid = v_schw[np.isfinite(v_schw)]
    v_mink_valid = v_mink[np.isfinite(v_mink)]
    v_exact_valid = v_exact[np.isfinite(v_exact)]
    if len(v_schw_valid) > 0 and len(v_mink_valid) > 0:
        v_min = min(v_schw_valid.min(), v_mink_valid.min(), v_exact_valid.min()) - 0.1
        v_max = max(v_schw_valid.max(), v_mink_valid.max(), v_exact_valid.max()) + 0.1
    else:
        v_min, v_max = -1.0, 1.0
    ax.set_ylim(v_min, v_max)
    ax.legend(fontsize=8)
    ax.set_title('Comparison: Velocity (negative = GR infall)')
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.scatter(r_schw, P_schw, s=3, alpha=0.6, c='green', label='Schwarzschild')
    ax.scatter(r_mink, P_mink, s=3, alpha=0.4, c='blue', label='Minkowski')
    ax.plot(r_ana, P_exact, 'r--', lw=1, alpha=0.5, label='Exact SR')
    ax.axvline(2.0, color='black', linestyle='--', lw=1.5, alpha=0.5)
    ax.set_xlabel('r')
    ax.set_ylabel('Pressure')
    ax.set_xlim(3, 15)
    ax.set_ylim(-0.5, 20)
    ax.legend(fontsize=8)
    ax.set_title('Comparison: Pressure')
    ax.grid(True, alpha=0.3)

    fig.suptitle('GR-GSPH Verification: Minkowski Should Match Exact SR\n'
                 'Top: Minkowski (blue) vs Exact (red) | Bottom: Schwarzschild (green) vs Minkowski (blue)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()

    # Save
    output_file = script_dir.parent / 'results' / 'verification_comparison.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison to {output_file}")

    # Compute L2 error for Minkowski vs Exact
    # Interpolate exact solution to particle positions
    n_exact_interp = np.interp(r_mink, r_ana, n_exact)
    v_exact_interp = np.interp(r_mink, r_ana, v_exact)
    P_exact_interp = np.interp(r_mink, r_ana, P_exact)

    L2_n = np.sqrt(np.mean((n_mink - n_exact_interp)**2))
    L2_v = np.sqrt(np.mean((v_mink - v_exact_interp)**2))
    L2_P = np.sqrt(np.mean((P_mink - P_exact_interp)**2))

    print("\n" + "=" * 50)
    print("  Verification: L2 Error (Minkowski vs Exact SR)")
    print("=" * 50)
    print(f"  Density:  L2 = {L2_n:.6f}")
    print(f"  Velocity: L2 = {L2_v:.6f}")
    print(f"  Pressure: L2 = {L2_P:.6f}")
    print("=" * 50)

    if L2_n < 0.5 and L2_v < 0.1 and L2_P < 1.0:
        print("  ✓ PASS: Minkowski matches Exact SR (small L2 error)")
    else:
        print("  ✗ FAIL: Large deviation from Exact SR")

    plt.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
