#!/usr/bin/env python3
"""
Animated Comparison: Exact vs HLLC Riemann Solvers for SR-GSPH
Kitajima et al. (2025) Tangent Velocity Test (v_t = 0.9)

Creates animated GIF overlaying both solvers with analytical solution.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
import glob
import sys

# Add SRRP library path
SRRP_PATH = Path("/Users/guo/Downloads/sph-simulators/docs/papers/sg-gsph/srrp")
sys.path.insert(0, str(SRRP_PATH))
from srrp.Solver import Solver
from srrp.State import State

# Configuration
GAMMA = 5.0 / 3.0
BASE_DIR = Path(__file__).parent.parent / 'results'
EXACT_DIR = BASE_DIR / 'kitajima_vt09'
HLLC_DIR = BASE_DIR / 'kitajima_vt09_hllc'
OUTPUT_DIR = BASE_DIR / 'comparison'

# Initial conditions (Kitajima Problem 5)
P_L, n_L, vx_L, vt_L = 1000.0, 1.0, 0.0, 0.9
P_R, n_R, vx_R, vt_R = 0.01, 1.0, 0.0, 0.9
X_CONTACT = 0.0

# Solve Riemann problem once
_solver = Solver()
_stateL = State(rho=n_L, vx=vx_L, vt=vt_L, pressure=P_L)
_stateR = State(rho=n_R, vx=vx_R, vt=vt_R, pressure=P_R)
_wavefan = _solver.solve(_stateL, _stateR, GAMMA)

# Analytical shock speed
V_SHOCK = 0.4451

def load_snapshot(filepath):
    """Load a CSV snapshot, skipping comment lines"""
    df = pd.read_csv(filepath, comment='#')
    df = df.dropna(subset=['pos_x'])
    return df

def compute_exact_solution(x_array, t):
    """Compute exact solution at given positions and time using SRRP library."""
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

def find_shock_position(x, dens):
    """Find shock position where density drops from high to low"""
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    dens_sorted = dens[sort_idx]

    # Find where density transitions from >3 to <2 (shock front)
    right_mask = x_sorted > 0.02
    x_right = x_sorted[right_mask]
    dens_right = dens_sorted[right_mask]

    if len(dens_right) < 10:
        return None

    # Find last position where density > 2.5
    high_dens_mask = dens_right > 2.5
    if any(high_dens_mask):
        return x_right[high_dens_mask][-1]
    return None

def main():
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Find snapshots
    exact_files = sorted(glob.glob(str(EXACT_DIR / 'snapshot_*.csv')))
    hllc_files = sorted(glob.glob(str(HLLC_DIR / 'snapshot_*.csv')))

    print(f"Found {len(exact_files)} exact snapshots, {len(hllc_files)} HLLC snapshots")

    if not exact_files or not hllc_files:
        print("No snapshots found!")
        return

    # Use common frames, limit to t=0.45 (frame 15, since t = frame * 0.03)
    MAX_TIME = 0.45
    DT_OUTPUT = 0.03
    max_frame = int(MAX_TIME / DT_OUTPUT) + 1  # +1 to include t=0.45
    n_frames = min(len(exact_files), len(hllc_files), max_frame)
    print(f"Creating animation with {n_frames} frames (t=0 to t={MAX_TIME})...")

    # Set up figure - 1x3 layout like Kitajima et al.
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Iterative vs HLLC Riemann Solver: Kitajima et al. (2025) Problem 5 ($v_t = 0.9$)', fontsize=12)

    # X range for analytical solution
    x_exact_line = np.linspace(-0.5, 0.5, 1000)

    # Initialize plot elements
    lines = {}

    # Colors: Blue for Iterative, Red for HLLC, Black for Analytical
    COLOR_ITER = 'blue'
    COLOR_HLLC = 'red'
    COLOR_ANA = 'black'

    # Left panel: P/1000 and n/5
    ax = axes[0]
    # Density (n/5)
    lines['exact_dens'], = ax.plot([], [], 'o', color=COLOR_ITER, ms=2, alpha=0.6, label='Iterative')
    lines['hllc_dens'], = ax.plot([], [], 's', color=COLOR_HLLC, ms=2, alpha=0.6, label='HLLC')
    lines['ana_dens'], = ax.plot([], [], '-', color=COLOR_ANA, lw=1.5, label='Analytical')
    # Pressure (P/1000) - same colors, different line style indicator
    lines['exact_pres'], = ax.plot([], [], 'o', color=COLOR_ITER, ms=2, alpha=0.6)
    lines['hllc_pres'], = ax.plot([], [], 's', color=COLOR_HLLC, ms=2, alpha=0.6)
    lines['ana_pres'], = ax.plot([], [], '-', color=COLOR_ANA, lw=1.5)
    ax.set_xlabel('x', fontsize=12)
    ax.text(0.05, 0.95, r'$P/1000$', transform=ax.transAxes, fontsize=14, color='darkgreen', va='top')
    ax.text(0.05, 0.85, r'$n/5$', transform=ax.transAxes, fontsize=14, color='darkred', va='top')
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 1.2)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Middle panel: v^x
    ax = axes[1]
    lines['exact_velx'], = ax.plot([], [], 'o', color=COLOR_ITER, ms=2, alpha=0.6, label='Iterative')
    lines['hllc_velx'], = ax.plot([], [], 's', color=COLOR_HLLC, ms=2, alpha=0.6, label='HLLC')
    lines['ana_velx'], = ax.plot([], [], '-', color=COLOR_ANA, lw=1.5, label='Analytical')
    ax.set_xlabel('x', fontsize=12)
    ax.text(0.05, 0.95, r'$v^x$', transform=ax.transAxes, fontsize=14, color='darkblue', va='top')
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 0.4)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Right panel: v^t
    ax = axes[2]
    lines['exact_velt'], = ax.plot([], [], 'o', color=COLOR_ITER, ms=2, alpha=0.6, label='Iterative')
    lines['hllc_velt'], = ax.plot([], [], 's', color=COLOR_HLLC, ms=2, alpha=0.6, label='HLLC')
    lines['ana_velt'], = ax.plot([], [], '-', color=COLOR_ANA, lw=1.5, label='Analytical')
    ax.set_xlabel('x', fontsize=12)
    ax.text(0.05, 0.95, r'$v^t$', transform=ax.transAxes, fontsize=14, color='darkcyan', va='top')
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0.4, 1.0)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Text annotation for time and shock info
    time_text = fig.text(0.5, 0.01, '', fontsize=10, ha='center',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])

    def update(frame):
        # Load data
        exact_df = load_snapshot(exact_files[frame])
        hllc_df = load_snapshot(hllc_files[frame])

        snap_num = int(exact_files[frame].split('_')[-1].replace('.csv', ''))
        t = snap_num * 0.03

        # Compute analytical solution
        rho_ex, pres_ex, vel_ex, vt_ex = compute_exact_solution(x_exact_line, t)

        # Update density (scaled by 1/5)
        lines['exact_dens'].set_data(exact_df['pos_x'], exact_df['dens'] / 5.0)
        lines['hllc_dens'].set_data(hllc_df['pos_x'], hllc_df['dens'] / 5.0)
        lines['ana_dens'].set_data(x_exact_line, rho_ex / 5.0)

        # Update pressure (scaled by 1/1000)
        lines['exact_pres'].set_data(exact_df['pos_x'], exact_df['pres'] / 1000.0)
        lines['hllc_pres'].set_data(hllc_df['pos_x'], hllc_df['pres'] / 1000.0)
        lines['ana_pres'].set_data(x_exact_line, pres_ex / 1000.0)

        # Update normal velocity
        lines['exact_velx'].set_data(exact_df['pos_x'], exact_df['vel_x'])
        lines['hllc_velx'].set_data(hllc_df['pos_x'], hllc_df['vel_x'])
        lines['ana_velx'].set_data(x_exact_line, vel_ex)

        # Update tangent velocity
        lines['exact_velt'].set_data(exact_df['pos_x'], exact_df['vel_t'])
        lines['hllc_velt'].set_data(hllc_df['pos_x'], hllc_df['vel_t'])
        lines['ana_velt'].set_data(x_exact_line, vt_ex)

        # Find shock positions
        x_shock_exact = find_shock_position(exact_df['pos_x'].values, exact_df['dens'].values)
        x_shock_hllc = find_shock_position(hllc_df['pos_x'].values, hllc_df['dens'].values)
        x_shock_ana = V_SHOCK * t if t > 0 else 0

        # Build info text
        text_parts = [f"t = {t:.3f}"]
        text_parts.append(f"Analytical: x_shock = {x_shock_ana:.3f}")

        if x_shock_exact and t > 0:
            err_exact = 100 * (x_shock_exact/t - V_SHOCK) / V_SHOCK
            text_parts.append(f"Exact: x = {x_shock_exact:.3f} ({err_exact:+.1f}%)")

        if x_shock_hllc and t > 0:
            err_hllc = 100 * (x_shock_hllc/t - V_SHOCK) / V_SHOCK
            text_parts.append(f"HLLC: x = {x_shock_hllc:.3f} ({err_hllc:+.1f}%)")

        time_text.set_text("  |  ".join(text_parts))

        print(f"  Frame {frame+1}/{n_frames} (t={t:.3f})")
        return list(lines.values()) + [time_text]

    # Create animation
    anim = FuncAnimation(fig, update, frames=n_frames, blit=False, interval=100)

    # Save animation
    output_file = OUTPUT_DIR / 'exact_vs_hllc_animation.gif'
    print(f"\nSaving animation to: {output_file}")
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=100)

    print(f"\nDone! Animation saved to: {output_file}")

if __name__ == '__main__':
    main()
