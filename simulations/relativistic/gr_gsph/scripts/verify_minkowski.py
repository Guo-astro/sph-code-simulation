#!/usr/bin/env python3
"""
Quick verification: GR-GSPH Minkowski should match exact SR solution.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / 'docs'))

from srrp.Solver import Solver as SRRPSolver
from srrp.State import State


def load_csv(filename):
    data = pd.read_csv(filename, comment='#')
    time = 0.0
    with open(filename, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                time = float(line.split(':')[1].strip())
                break
    return data, time


def main():
    script_dir = Path(__file__).parent

    # Load snapshot at t ≈ 1.0
    mink_file = script_dir.parent / 'results' / 'minkowski_shock' / 'snapshot_0010.csv'
    schw_file = script_dir.parent / 'results' / 'schwarzschild_shock' / 'snapshot_0010.csv'

    mink_data, mink_t = load_csv(mink_file)
    schw_data, schw_t = load_csv(schw_file)

    print(f"Minkowski t = {mink_t:.4f}")
    print(f"Schwarzschild t = {schw_t:.4f}")

    # Filter ghosts
    mink_mask = mink_data['is_ghost'] == 0
    schw_mask = schw_data['is_ghost'] == 0

    r_mink = mink_data['pos_x'].values[mink_mask]
    n_mink = mink_data['dens'].values[mink_mask]
    v_mink = mink_data['vel_x'].values[mink_mask]
    P_mink = mink_data['pres'].values[mink_mask]

    r_schw = schw_data['pos_x'].values[schw_mask]
    n_schw = schw_data['dens'].values[schw_mask]
    v_schw = schw_data['vel_x'].values[schw_mask]
    P_schw = schw_data['pres'].values[schw_mask]

    # Sort
    mink_sort = np.argsort(r_mink)
    schw_sort = np.argsort(r_schw)
    r_mink, n_mink, v_mink, P_mink = r_mink[mink_sort], n_mink[mink_sort], v_mink[mink_sort], P_mink[mink_sort]
    r_schw, n_schw, v_schw, P_schw = r_schw[schw_sort], n_schw[schw_sort], v_schw[schw_sort], P_schw[schw_sort]

    # Exact SR solution
    solver = SRRPSolver()
    stateL = State(rho=10.0, vx=0.0, vt=0.0, pressure=40.0/3.0)
    stateR = State(rho=1.0, vx=0.0, vt=0.0, pressure=1e-6)
    wavefan = solver.solve(stateL, stateR, 5.0/3.0)

    r_ana = np.linspace(3.0, 15.0, 500)
    xi = (r_ana - 6.0) / mink_t

    # getState returns aggregated State with array attributes when given array input
    state_all = wavefan.getState(xi)
    if hasattr(state_all, 'rho') and isinstance(state_all.rho, np.ndarray):
        n_exact = state_all.rho
        v_exact = state_all.vx
        P_exact = state_all.pressure
    else:
        # Fallback: loop over each point
        n_exact = np.zeros_like(xi)
        v_exact = np.zeros_like(xi)
        P_exact = np.zeros_like(xi)
        for i, x in enumerate(xi):
            s = wavefan.getState(np.array([x]))
            n_exact[i] = s.rho[0] if hasattr(s.rho, '__len__') else s.rho
            v_exact[i] = s.vx[0] if hasattr(s.vx, '__len__') else s.vx
            P_exact[i] = s.pressure[0] if hasattr(s.pressure, '__len__') else s.pressure

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Top: Minkowski vs Exact
    ax = axes[0, 0]
    ax.scatter(r_mink, n_mink, s=5, alpha=0.7, c='blue', label='GR-GSPH Minkowski')
    ax.plot(r_ana, n_exact, 'r-', lw=2, label='Exact SR')
    ax.set_xlabel('r'); ax.set_ylabel('Density')
    ax.set_xlim(3, 15); ax.set_ylim(0, 12)
    ax.legend(); ax.set_title(f'Density (t={mink_t:.3f})')
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.scatter(r_mink, v_mink, s=5, alpha=0.7, c='blue', label='GR-GSPH Minkowski')
    ax.plot(r_ana, v_exact, 'r-', lw=2, label='Exact SR')
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel('r'); ax.set_ylabel('Velocity')
    ax.set_xlim(3, 15); ax.set_ylim(-0.2, 1.0)
    ax.legend(); ax.set_title('Velocity')
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    ax.scatter(r_mink, P_mink, s=5, alpha=0.7, c='blue', label='GR-GSPH Minkowski')
    ax.plot(r_ana, P_exact, 'r-', lw=2, label='Exact SR')
    ax.set_xlabel('r'); ax.set_ylabel('Pressure')
    ax.set_xlim(3, 15); ax.set_ylim(-0.5, 15)
    ax.legend(); ax.set_title('Pressure')
    ax.grid(True, alpha=0.3)

    # Bottom: Schwarzschild vs Minkowski (GR effects)
    ax = axes[1, 0]
    ax.scatter(r_schw, n_schw, s=5, alpha=0.7, c='green', label='Schwarzschild')
    ax.scatter(r_mink, n_mink, s=5, alpha=0.5, c='blue', label='Minkowski')
    ax.plot(r_ana, n_exact, 'r--', lw=1, alpha=0.5, label='Exact SR')
    ax.axvline(2.0, color='k', ls='--', lw=1.5, alpha=0.5, label='Horizon r=2M')
    ax.set_xlabel('r'); ax.set_ylabel('Density')
    ax.set_xlim(3, 15); ax.set_ylim(0, 12)
    ax.legend(fontsize=8); ax.set_title(f'Density (t={schw_t:.3f})')
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.scatter(r_schw, v_schw, s=5, alpha=0.7, c='green', label='Schwarzschild')
    ax.scatter(r_mink, v_mink, s=5, alpha=0.5, c='blue', label='Minkowski')
    ax.plot(r_ana, v_exact, 'r--', lw=1, alpha=0.5, label='Exact SR')
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(2.0, color='k', ls='--', lw=1.5, alpha=0.5)
    ax.set_xlabel('r'); ax.set_ylabel('Velocity')
    ax.set_xlim(3, 15); ax.set_ylim(-0.5, 1.0)
    ax.legend(fontsize=8); ax.set_title('Velocity (negative = infall to BH)')
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.scatter(r_schw, P_schw, s=5, alpha=0.7, c='green', label='Schwarzschild')
    ax.scatter(r_mink, P_mink, s=5, alpha=0.5, c='blue', label='Minkowski')
    ax.plot(r_ana, P_exact, 'r--', lw=1, alpha=0.5, label='Exact SR')
    ax.axvline(2.0, color='k', ls='--', lw=1.5, alpha=0.5)
    ax.set_xlabel('r'); ax.set_ylabel('Pressure')
    ax.set_xlim(3, 15); ax.set_ylim(-0.5, 15)
    ax.legend(fontsize=8); ax.set_title('Pressure')
    ax.grid(True, alpha=0.3)

    fig.suptitle('GR-GSPH Verification at t ≈ 1.0\n'
                 'Top: Minkowski (blue) should match Exact SR (red)\n'
                 'Bottom: Schwarzschild (green) differs due to gravity',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    output = script_dir.parent / 'results' / 'verification_t1.png'
    fig.savefig(output, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output}")

    # L2 error
    n_interp = np.interp(r_mink, r_ana, n_exact)
    v_interp = np.interp(r_mink, r_ana, v_exact)
    P_interp = np.interp(r_mink, r_ana, P_exact)

    L2_n = np.sqrt(np.mean((n_mink - n_interp)**2))
    L2_v = np.sqrt(np.mean((v_mink - v_interp)**2))
    L2_P = np.sqrt(np.mean((P_mink - P_interp)**2))

    print("\n" + "=" * 50)
    print("  L2 Error: Minkowski vs Exact SR (t ≈ 1.0)")
    print("=" * 50)
    print(f"  Density:  {L2_n:.6f}")
    print(f"  Velocity: {L2_v:.6f}")
    print(f"  Pressure: {L2_P:.6f}")

    # Check if Schwarzschild differs from Minkowski
    v_schw_interp = np.interp(r_mink, r_schw, v_schw)
    v_diff = np.mean(v_schw_interp - v_mink)
    print(f"\n  Mean velocity difference (Schw - Mink): {v_diff:.6f}")
    if v_diff < -0.01:
        print("  → Schwarzschild has infall toward BH (correct!)")

    plt.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
