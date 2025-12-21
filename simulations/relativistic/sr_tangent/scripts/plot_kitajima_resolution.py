#!/usr/bin/env python3
"""
Create Kitajima-style resolution comparison plot (Figure 10 style)
4 rows x 3 columns, with HLLC overlay and analytical solution

Used by: make -f Makefile.kitajima plot
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

# Add SRRP library path
SRRP_PATH = Path("/Users/guo/Downloads/sph-simulators/docs/papers/sg-gsph/srrp")
sys.path.insert(0, str(SRRP_PATH))
from srrp.Solver import Solver
from srrp.State import State

# Configuration
SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR.parent / 'results' / 'kitajima_resolution'
DT_OUTPUT = 0.03

# Initial conditions (Kitajima Problem 5)
GAMMA = 5.0 / 3.0
P_L, n_L, vx_L, vt_L = 1000.0, 1.0, 0.0, 0.9
P_R, n_R, vx_R, vt_R = 0.01, 1.0, 0.0, 0.9
X_CONTACT = 0.0

# Solve Riemann problem once
_solver = Solver()
_stateL = State(rho=n_L, vx=vx_L, vt=vt_L, pressure=P_L)
_stateR = State(rho=n_R, vx=vx_R, vt=vt_R, pressure=P_R)
_wavefan = _solver.solve(_stateL, _stateR, GAMMA)


def compute_analytical_solution(x_array, t):
    """Compute analytical solution at given positions and time using SRRP library."""
    if t <= 0:
        rho = np.where(x_array < X_CONTACT, n_L, n_R)
        P = np.where(x_array < X_CONTACT, P_L, P_R)
        vx = np.zeros_like(x_array)
        vt = np.full_like(x_array, vt_L)
        return rho, P, vx, vt

    # Self-similar variable
    xi = (x_array - X_CONTACT) / t

    # Get state at each xi from SRRP wavefan
    state = _wavefan.getState(xi)

    return state.rho, state.pressure, state.vx, state.vt


def load_snapshot(output_dir, snap_num=None):
    """Load snapshot from output directory."""
    files = sorted(glob.glob(str(output_dir / 'snapshot_*.csv')))
    if not files:
        return None
    if snap_num is not None:
        # Find specific snapshot
        for f in files:
            if f'snapshot_{snap_num:04d}.csv' in f or f'snapshot_{snap_num:05d}.csv' in f:
                df = pd.read_csv(f, comment='#').dropna(subset=['pos_x'])
                return df
    # Return last snapshot
    f = files[-1]
    df = pd.read_csv(f, comment='#').dropna(subset=['pos_x'])
    return df


def create_kitajima_style_plot():
    """Create 4x3 grid like Kitajima Figure 10 with HLLC overlay and analytical solution."""

    # Resolutions to show
    # All use snapshot 12 (t = 0.36) for consistent comparison
    SNAP_NUM = 12
    resolutions = [
        # (exact_dir, hllc_dir, label, snap_num)
        (RESULTS_DIR / 'exact_800_800', RESULTS_DIR / 'hllc_800_800', "800 + 800", SNAP_NUM),
        (RESULTS_DIR / 'exact_1600_1600', RESULTS_DIR / 'hllc_1600_1600', "1,600 + 1,600", SNAP_NUM),
        (RESULTS_DIR / 'exact_3200_1600', RESULTS_DIR / 'hllc_3200_1600', "3,200 + 1,600", SNAP_NUM),
        (RESULTS_DIR / 'exact_12800_1600', RESULTS_DIR / 'hllc_12800_1600', "12,800 + 1,600", SNAP_NUM),
    ]

    fig, axes = plt.subplots(4, 3, figsize=(12, 14))

    # High contrast colors
    # Exact: dark colors
    EXACT_P_COLOR = 'darkgreen'
    EXACT_N_COLOR = 'darkred'
    EXACT_VX_COLOR = 'darkblue'
    EXACT_VT_COLOR = 'purple'

    # HLLC: bright/light colors
    HLLC_P_COLOR = 'lime'
    HLLC_N_COLOR = 'orange'
    HLLC_VX_COLOR = 'deepskyblue'
    HLLC_VT_COLOR = 'magenta'

    # Analytical solution color
    ANALYTICAL_COLOR = 'black'

    # Compute analytical solution at t = SNAP_NUM * DT_OUTPUT
    t_snapshot = SNAP_NUM * DT_OUTPUT  # = 0.36 for snap 12
    x_analytical = np.linspace(-0.5, 0.5, 1000)
    rho_ana, pres_ana, vx_ana, vt_ana = compute_analytical_solution(x_analytical, t_snapshot)
    print(f"Using snapshot {SNAP_NUM}, t = {t_snapshot:.4f}")

    for row_idx, (exact_dir, hllc_dir, label, snap_num) in enumerate(resolutions):
        # Load exact and HLLC data
        exact_df = load_snapshot(exact_dir, snap_num)
        hllc_df = load_snapshot(hllc_dir, snap_num)

        if exact_df is None:
            print(f"Missing exact data for {label}: {exact_dir}")
            continue

        if hllc_df is None:
            print(f"Missing HLLC data for {label}: {hllc_dir}")

        # Column 1: P/1000 and n/5
        ax = axes[row_idx, 0]

        # Analytical solution - solid black lines (plot first so it's behind)
        ax.plot(x_analytical, pres_ana / 1000, '-', color=ANALYTICAL_COLOR, lw=1.5, alpha=0.8)
        ax.plot(x_analytical, rho_ana / 5, '-', color=ANALYTICAL_COLOR, lw=1.5, alpha=0.8)

        # Exact solver - x markers
        ax.plot(exact_df['pos_x'], exact_df['pres'] / 1000,
                'x', color=EXACT_P_COLOR, ms=3, mew=0.8, alpha=0.9)
        ax.plot(exact_df['pos_x'], exact_df['dens'] / 5,
                'x', color=EXACT_N_COLOR, ms=3, mew=0.8, alpha=0.9)

        # HLLC solver overlay - o markers
        if hllc_df is not None:
            ax.plot(hllc_df['pos_x'], hllc_df['pres'] / 1000,
                    'o', color=HLLC_P_COLOR, ms=4, mfc='none', mew=0.8, alpha=0.7)
            ax.plot(hllc_df['pos_x'], hllc_df['dens'] / 5,
                    'o', color=HLLC_N_COLOR, ms=4, mfc='none', mew=0.8, alpha=0.7)

        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.05, 1.2)
        ax.set_ylabel(label, fontsize=11, fontweight='bold')
        ax.text(-0.45, 1.05, 'P/1000', color=EXACT_P_COLOR, fontsize=10, fontweight='bold')
        ax.text(-0.45, 0.25, 'n/5', color=EXACT_N_COLOR, fontsize=10, fontweight='bold')
        if row_idx == 0:
            ax.set_title('Pressure & Density', fontsize=11)
        if row_idx == 3:
            ax.set_xlabel('x', fontsize=11)
        ax.grid(True, alpha=0.2)

        # Column 2: v^x
        ax = axes[row_idx, 1]

        # Analytical solution - solid black line
        ax.plot(x_analytical, vx_ana, '-', color=ANALYTICAL_COLOR, lw=1.5, alpha=0.8)

        # Exact - x markers
        ax.plot(exact_df['pos_x'], exact_df['vel_x'],
                'x', color=EXACT_VX_COLOR, ms=3, mew=0.8, alpha=0.9)
        # HLLC - o markers
        if hllc_df is not None:
            ax.plot(hllc_df['pos_x'], hllc_df['vel_x'],
                    'o', color=HLLC_VX_COLOR, ms=4, mfc='none', mew=0.8, alpha=0.7)

        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.02, 0.42)
        ax.text(-0.45, 0.36, r'$v^x$', color=EXACT_VX_COLOR, fontsize=12, fontweight='bold')
        if row_idx == 0:
            ax.set_title('Normal Velocity', fontsize=11)
        if row_idx == 3:
            ax.set_xlabel('x', fontsize=11)
        ax.grid(True, alpha=0.2)

        # Column 3: v^t
        ax = axes[row_idx, 2]

        # Analytical solution - solid black line
        ax.plot(x_analytical, vt_ana, '-', color=ANALYTICAL_COLOR, lw=1.5, alpha=0.8)

        # Exact - x markers
        if 'vel_t' in exact_df.columns:
            ax.plot(exact_df['pos_x'], exact_df['vel_t'],
                    'x', color=EXACT_VT_COLOR, ms=3, mew=0.8, alpha=0.9)
        # HLLC - o markers
        if hllc_df is not None and 'vel_t' in hllc_df.columns:
            ax.plot(hllc_df['pos_x'], hllc_df['vel_t'],
                    'o', color=HLLC_VT_COLOR, ms=4, mfc='none', mew=0.8, alpha=0.7)

        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(0.4, 1.02)
        ax.text(-0.45, 0.95, r'$v^t$', color=EXACT_VT_COLOR, fontsize=12, fontweight='bold')
        if row_idx == 0:
            ax.set_title('Tangent Velocity', fontsize=11)
        if row_idx == 3:
            ax.set_xlabel('x', fontsize=11)
        ax.grid(True, alpha=0.2)

    # Add legend with matching markers and colors
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=ANALYTICAL_COLOR, linestyle='-', linewidth=2, label='Analytical'),
        Line2D([0], [0], marker='x', color=EXACT_VX_COLOR, linestyle='None',
               markersize=8, markeredgewidth=1.5, label='Exact Riemann'),
        Line2D([0], [0], marker='o', color=HLLC_VX_COLOR, linestyle='None',
               markersize=8, markerfacecolor='none', markeredgewidth=1.5, label='HLLC'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=3,
               fontsize=11, bbox_to_anchor=(0.5, 0.98))

    plt.suptitle(r'Resolution Study: $v^t_L = v^t_R = 0.9$, $t = ' + f'{t_snapshot:.2f}$' + '\n',
                fontsize=13, fontweight='bold', y=1.01)

    plt.tight_layout()

    # Save to results directory
    save_path = RESULTS_DIR / 'kitajima_style_comparison.png'
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


if __name__ == '__main__':
    create_kitajima_style_plot()
