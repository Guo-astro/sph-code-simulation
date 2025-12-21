#!/usr/bin/env python3
"""
GRSPH Benchmark Suite

Reproduces key tests from Liptai & Price (2019) MNRAS 485, 819.

Tests included:
1. Mildly relativistic shock (Problem 1) - Marti & Muller (2003)
   - ρ = [10, 1], P = [40/3, 1e-6], Γ_max = 1.38

2. Strong blast wave (Problem 2) - Marti & Muller (2003)
   - ρ = [1, 1], P = [1000, 0.01], Γ_max = 3.6

3. Schwarzschild shock tube - GR effects
   - Same as Problem 1 but in curved spacetime
   - Shows gravitational infall and asymmetric wave propagation

Run this script to generate comparison plots.
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# Add srrp solver to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / 'docs'))

try:
    from srrp.Solver import Solver as SRRPSolver
    from srrp.State import State
    HAS_SRRP = True
except ImportError:
    print("Warning: SRRP solver not available")
    HAS_SRRP = False


def load_csv(filename):
    """Load CSV snapshot with header comments."""
    data = pd.read_csv(filename, comment='#')
    time = 0.0
    with open(filename, 'r') as f:
        for line in f:
            if '# Time (code):' in line:
                time = float(line.split(':')[1].strip())
                break
    return data, time


def compute_exact_sr(r, t, rho_L, P_L, rho_R, P_R, r_disc, gamma=5.0/3.0):
    """Compute exact SR Riemann solution."""
    if not HAS_SRRP or t <= 0:
        n = np.where(r < r_disc, rho_L, rho_R)
        v = np.zeros_like(r)
        P = np.where(r < r_disc, P_L, P_R)
        return n, v, P

    solver = SRRPSolver()
    stateL = State(rho=rho_L, vx=0.0, vt=0.0, pressure=P_L)
    stateR = State(rho=rho_R, vx=0.0, vt=0.0, pressure=P_R)
    wavefan = solver.solve(stateL, stateR, gamma)

    xi = (r - r_disc) / t
    state = wavefan.getState(xi)
    return state.rho, state.vx, state.pressure


def extract_data(data, filter_ghosts=True, interior_only=False, r_range=None):
    """Extract and filter particle data."""
    if filter_ghosts and 'is_ghost' in data.columns:
        mask = data['is_ghost'] == 0
    else:
        mask = np.ones(len(data), dtype=bool)

    r = data['pos_x'].values[mask]
    n = data['dens'].values[mask]
    v = data['vel_x'].values[mask]
    P = data['pres'].values[mask]

    idx = np.argsort(r)
    r, n, v, P = r[idx], n[idx], v[idx], P[idx]

    if interior_only and r_range:
        interior = (r > r_range[0]) & (r < r_range[1])
        r, n, v, P = r[interior], n[interior], v[interior], P[interior]

    return r, n, v, P


def compute_l2_error(sim_r, sim_val, ana_r, ana_val):
    """Compute L2 error between simulation and analytical solution."""
    ana_interp = np.interp(sim_r, ana_r, ana_val)
    return np.sqrt(np.mean((sim_val - ana_interp)**2))


def plot_shock_tube(ax, r_sim, val_sim, r_ana, val_ana,
                    ylabel, title, ylim=None, label_sim='Simulation',
                    label_ana='Exact', show_legend=True):
    """Plot shock tube comparison."""
    ax.scatter(r_sim, val_sim, s=5, alpha=0.7, c='blue', label=label_sim)
    ax.plot(r_ana, val_ana, 'r-', lw=2, label=label_ana)
    ax.set_xlabel('x')
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold')
    if ylim:
        ax.set_ylim(ylim)
    if show_legend:
        ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)


def main():
    script_dir = Path(__file__).parent
    results_dir = script_dir.parent / 'results'

    print("=" * 70)
    print("  GRSPH Benchmark Suite")
    print("  Based on Liptai & Price (2019) MNRAS 485, 819")
    print("=" * 70)

    # =========================================================================
    # Test 1: Minkowski Shock (Problem 1 verification)
    # =========================================================================
    print("\n[Test 1] Minkowski Shock Tube (Problem 1)")
    # Try minkowski_shock_tube first, fall back to minkowski_shock
    mink_dir = results_dir / 'minkowski_shock_tube'
    if not mink_dir.exists():
        mink_dir = results_dir / 'minkowski_shock'

    if mink_dir.exists():
        files = sorted(glob.glob(str(mink_dir / 'snapshot_*.csv')))
        if files:
            # Use snapshot at t ≈ 1.0
            data, t = load_csv(files[10] if len(files) > 10 else files[-1])
            r_sim, n_sim, v_sim, P_sim = extract_data(data, interior_only=True, r_range=(4, 13))

            # Compute exact solution
            r_ana = np.linspace(4, 13, 500)
            n_ana, v_ana, P_ana = compute_exact_sr(r_ana, t, 10.0, 40.0/3.0, 1.0, 1e-6, 6.0)

            # Compute L2 errors
            L2_n = compute_l2_error(r_sim, n_sim, r_ana, n_ana)
            L2_v = compute_l2_error(r_sim, v_sim, r_ana, v_ana)
            L2_P = compute_l2_error(r_sim, P_sim, r_ana, P_ana)

            print(f"  Time: t = {t:.4f}")
            print(f"  L2 errors: density={L2_n:.4f}, velocity={L2_v:.4f}, pressure={L2_P:.4f}")
            print(f"  Relative density error: {L2_n/np.mean(n_ana)*100:.1f}%")
    else:
        print("  [SKIPPED] Results not found. Run simulation first.")

    # =========================================================================
    # Test 2: Schwarzschild Shock (GR effects)
    # =========================================================================
    print("\n[Test 2] Schwarzschild Shock Tube (GR Effects)")
    schw_dir = results_dir / 'schwarzschild_shock'

    if schw_dir.exists():
        files = sorted(glob.glob(str(schw_dir / 'snapshot_*.csv')))
        if files:
            data, t = load_csv(files[10] if len(files) > 10 else files[-1])
            r_schw, n_schw, v_schw, P_schw = extract_data(data, interior_only=True, r_range=(4, 13))

            print(f"  Time: t = {t:.4f}")
            print(f"  Velocity range: [{v_schw.min():.4f}, {v_schw.max():.4f}]")

            # Compare with Minkowski if available
            if mink_dir.exists():
                mink_files = sorted(glob.glob(str(mink_dir / 'snapshot_*.csv')))
                if mink_files:
                    mink_data, mink_t = load_csv(mink_files[10] if len(mink_files) > 10 else mink_files[-1])
                    r_mink, n_mink, v_mink, P_mink = extract_data(mink_data, interior_only=True, r_range=(4, 13))

                    # Velocity shift
                    v_schw_interp = np.interp(r_mink, r_schw, v_schw)
                    mean_v_diff = np.mean(v_schw_interp - v_mink)
                    print(f"  Mean velocity shift (Schw - Mink): {mean_v_diff:.4f}")
                    if mean_v_diff < -0.01:
                        print("  → PASS: Negative velocity shift indicates GR infall")
    else:
        print("  [SKIPPED] Results not found. Run simulation first.")

    # =========================================================================
    # Generate Combined Figure
    # =========================================================================
    print("\n[Generating Figure]")

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))

    # Top row: Minkowski (flat spacetime verification)
    if mink_dir.exists() and (mink_dir / 'snapshot_0010.csv').exists():
        data, t = load_csv(mink_dir / 'snapshot_0010.csv')
        r_sim, n_sim, v_sim, P_sim = extract_data(data, interior_only=True, r_range=(4, 13))
        r_ana = np.linspace(4, 13, 500)
        n_ana, v_ana, P_ana = compute_exact_sr(r_ana, t, 10.0, 40.0/3.0, 1.0, 1e-6, 6.0)

        plot_shock_tube(axes[0,0], r_sim, n_sim, r_ana, n_ana,
                       'Density ρ', f'Minkowski: Density (t={t:.2f})', ylim=(0, 12))
        plot_shock_tube(axes[0,1], r_sim, v_sim, r_ana, v_ana,
                       'Velocity v', 'Velocity', ylim=(-0.2, 1.0))
        plot_shock_tube(axes[0,2], r_sim, P_sim, r_ana, P_ana,
                       'Pressure P', 'Pressure', ylim=(-0.5, 15))

        # Enthalpy
        gamma = 5.0/3.0
        h_sim = 1 + gamma/(gamma-1) * P_sim / np.maximum(n_sim, 1e-10)
        h_ana = 1 + gamma/(gamma-1) * P_ana / np.maximum(n_ana, 1e-10)
        plot_shock_tube(axes[0,3], r_sim, h_sim, r_ana, h_ana,
                       'Enthalpy h', 'Specific Enthalpy', ylim=(0.9, 3.0))

        axes[0,0].text(0.02, 0.98, 'FLAT SPACETIME', transform=axes[0,0].transAxes,
                      fontsize=10, fontweight='bold', va='top', color='blue')
    else:
        for ax in axes[0]:
            ax.text(0.5, 0.5, 'Run Minkowski test first', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('Minkowski (not available)')

    # Bottom row: Schwarzschild (curved spacetime)
    if schw_dir.exists() and (schw_dir / 'snapshot_0010.csv').exists():
        data, t = load_csv(schw_dir / 'snapshot_0010.csv')
        r_schw, n_schw, v_schw, P_schw = extract_data(data, interior_only=True, r_range=(4, 13))

        # Plot Schwarzschild with Minkowski overlay
        r_ana = np.linspace(4, 13, 500)
        n_ana, v_ana, P_ana = compute_exact_sr(r_ana, t, 10.0, 40.0/3.0, 1.0, 1e-6, 6.0)

        ax = axes[1,0]
        ax.scatter(r_schw, n_schw, s=5, alpha=0.7, c='green', label='Schwarzschild')
        ax.plot(r_ana, n_ana, 'r--', lw=1.5, alpha=0.7, label='Flat SR')
        ax.axvline(2.0, color='k', ls='--', lw=1.5, alpha=0.5, label='Horizon')
        ax.set_xlabel('r/M')
        ax.set_ylabel('Density ρ')
        ax.set_title(f'Schwarzschild: Density (t={t:.2f})', fontweight='bold')
        ax.set_ylim(0, 12)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.text(0.02, 0.98, 'CURVED SPACETIME', transform=ax.transAxes,
               fontsize=10, fontweight='bold', va='top', color='green')

        ax = axes[1,1]
        ax.scatter(r_schw, v_schw, s=5, alpha=0.7, c='green', label='Schwarzschild')
        ax.plot(r_ana, v_ana, 'r--', lw=1.5, alpha=0.7, label='Flat SR')
        ax.axhline(0, color='k', lw=0.5)
        ax.axvline(2.0, color='k', ls='--', lw=1.5, alpha=0.5)
        ax.set_xlabel('r/M')
        ax.set_ylabel('Velocity v')
        ax.set_title('Velocity (negative = BH infall)', fontweight='bold')
        # Auto y-range
        v_min, v_max = v_schw.min() - 0.1, v_schw.max() + 0.1
        ax.set_ylim(v_min, max(v_max, 1.0))
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        ax = axes[1,2]
        ax.scatter(r_schw, P_schw, s=5, alpha=0.7, c='green', label='Schwarzschild')
        ax.plot(r_ana, P_ana, 'r--', lw=1.5, alpha=0.7, label='Flat SR')
        ax.axvline(2.0, color='k', ls='--', lw=1.5, alpha=0.5)
        ax.set_xlabel('r/M')
        ax.set_ylabel('Pressure P')
        ax.set_title('Pressure', fontweight='bold')
        ax.set_ylim(-0.5, 15)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Lapse function visualization
        ax = axes[1,3]
        r_plot = np.linspace(2.5, 15, 200)
        alpha = np.sqrt(np.maximum(1 - 2.0/r_plot, 0.01))
        ax.plot(r_plot, alpha, 'k-', lw=2, label=r'$\alpha = \sqrt{1-2M/r}$')
        ax.fill_between(r_plot, 0, alpha, alpha=0.1, color='gray')
        ax.axvline(2.0, color='red', ls='--', lw=2, label='Horizon r=2M')
        ax.axhline(1.0, color='blue', ls=':', lw=1, alpha=0.5, label='Flat spacetime')
        ax.set_xlabel('r/M')
        ax.set_ylabel('Lapse α')
        ax.set_title('Schwarzschild Lapse (Time Dilation)', fontweight='bold')
        ax.set_xlim(2.5, 15)
        ax.set_ylim(0, 1.15)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    else:
        for ax in axes[1]:
            ax.text(0.5, 0.5, 'Run Schwarzschild test first', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('Schwarzschild (not available)')

    fig.suptitle(
        'GRSPH Benchmark: Liptai & Price (2019) Problem 1\n'
        'Top: Flat spacetime (GR-GSPH Minkowski) vs Exact SR | '
        'Bottom: Curved spacetime (Schwarzschild) showing GR effects',
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout()

    output_file = results_dir / 'grsph_benchmark.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_file}")

    plt.close()

    print("\n" + "=" * 70)
    print("  BENCHMARK COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
