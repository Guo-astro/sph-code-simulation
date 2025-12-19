#!/usr/bin/env python3
"""
Compare Exact vs HLLC Riemann Solvers for SR-GSPH
Kitajima et al. (2025) Tangent Velocity Test (v_t = 0.9)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

# Add SRRP library path (same as animate_kitajima_vt09.py)
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

# Extract exact solution values
_states = _wavefan.states
EXACT_SOLUTION = {
    'P_star': _states[1].pressure,
    'vx_star': _states[1].vx,
    'vt_L_star': _states[1].vt,
    'vt_R_star': _states[2].vt,
    'rho_L_prime': _states[1].rho,
    'rho_R_prime': _states[2].rho,
}

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

    print(f"\n{'='*70}")
    print(f"  Exact vs HLLC Riemann Solver Comparison")
    print(f"  Kitajima et al. (2025) Problem 5 (v_t = 0.9)")
    print(f"{'='*70}")
    print(f"\nExact Riemann solution:")
    print(f"  P*     = {EXACT_SOLUTION['P_star']:.4f}")
    print(f"  v_x*   = {EXACT_SOLUTION['vx_star']:.4f}")
    print(f"  v_t_L* = {EXACT_SOLUTION['vt_L_star']:.4f}")
    print(f"  v_t_R* = {EXACT_SOLUTION['vt_R_star']:.4f}")

    # Create comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Exact vs HLLC Riemann Solver Comparison\n' +
                 r'Kitajima et al. (2025) Problem 5: $v_t^L = v_t^R = 0.9$', fontsize=14)

    # Use the most recent common snapshot
    n_common = min(len(exact_files), len(hllc_files)) - 1

    exact_df = load_snapshot(exact_files[n_common])
    hllc_df = load_snapshot(hllc_files[n_common])

    snap_num = int(exact_files[n_common].split('_')[-1].replace('.csv', ''))
    t = snap_num * 0.03

    # Exact solution
    x_exact_line = np.linspace(-0.5, 0.5, 1000)
    rho_ex, pres_ex, vel_ex, vt_ex = compute_exact_solution(x_exact_line, t)

    # Get shock positions
    x_shock_exact_sim = find_shock_position(exact_df['pos_x'].values, exact_df['dens'].values)
    x_shock_hllc_sim = find_shock_position(hllc_df['pos_x'].values, hllc_df['dens'].values)

    V_shock_analytical = 0.4451  # From exact solution

    # Plot Density
    ax = axes[0, 0]
    ax.plot(exact_df['pos_x'], exact_df['dens'], 'b.', ms=2, alpha=0.6, label='Exact Solver')
    ax.plot(hllc_df['pos_x'], hllc_df['dens'], 'r.', ms=2, alpha=0.6, label='HLLC Solver')
    ax.plot(x_exact_line, rho_ex, 'k-', lw=1.5, label='Analytical')
    ax.set_xlabel('x')
    ax.set_ylabel('Rest-Frame Density n')
    ax.set_title('Density')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xlim(-0.5, 0.5)
    ax.grid(True, alpha=0.3)

    # Plot Pressure
    ax = axes[0, 1]
    ax.plot(exact_df['pos_x'], exact_df['pres'], 'b.', ms=2, alpha=0.6, label='Exact Solver')
    ax.plot(hllc_df['pos_x'], hllc_df['pres'], 'r.', ms=2, alpha=0.6, label='HLLC Solver')
    ax.plot(x_exact_line, pres_ex, 'k-', lw=1.5, label='Analytical')
    ax.set_xlabel('x')
    ax.set_ylabel('Pressure P')
    ax.set_title('Pressure')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xlim(-0.5, 0.5)
    ax.grid(True, alpha=0.3)

    # Plot Normal Velocity
    ax = axes[1, 0]
    ax.plot(exact_df['pos_x'], exact_df['vel_x'], 'b.', ms=2, alpha=0.6, label='Exact Solver')
    ax.plot(hllc_df['pos_x'], hllc_df['vel_x'], 'r.', ms=2, alpha=0.6, label='HLLC Solver')
    ax.plot(x_exact_line, vel_ex, 'k-', lw=1.5, label='Analytical')
    ax.axhline(y=EXACT_SOLUTION['vx_star'], color='g', ls='--', alpha=0.5,
               label=f"$v_x^*$={EXACT_SOLUTION['vx_star']:.3f}")
    ax.set_xlabel('x')
    ax.set_ylabel(r'Normal Velocity $v^x$')
    ax.set_title('Normal Velocity')
    ax.legend(loc='upper left', fontsize=8)
    ax.set_xlim(-0.5, 0.5)
    ax.grid(True, alpha=0.3)

    # Plot Tangent Velocity
    ax = axes[1, 1]
    ax.plot(exact_df['pos_x'], exact_df['vel_t'], 'b.', ms=2, alpha=0.6, label='Exact Solver')
    ax.plot(hllc_df['pos_x'], hllc_df['vel_t'], 'r.', ms=2, alpha=0.6, label='HLLC Solver')
    ax.plot(x_exact_line, vt_ex, 'k-', lw=1.5, label='Analytical')
    ax.set_xlabel('x')
    ax.set_ylabel(r'Tangent Velocity $v^t$')
    ax.set_title('Tangent Velocity')
    ax.legend(loc='lower right', fontsize=8)
    ax.set_xlim(-0.5, 0.5)
    ax.grid(True, alpha=0.3)

    # Add shock position analysis text
    text_info = f"Time: t = {t:.3f}\n\n"
    text_info += f"Analytical shock:\n  x = {V_shock_analytical*t:.4f}\n  V = {V_shock_analytical:.4f}\n\n"

    if x_shock_exact_sim:
        V_exact = x_shock_exact_sim / t
        err_exact = 100 * (V_exact - V_shock_analytical) / V_shock_analytical
        text_info += f"Exact solver:\n  x = {x_shock_exact_sim:.4f}\n  V = {V_exact:.4f} ({err_exact:+.1f}%)\n\n"

    if x_shock_hllc_sim:
        V_hllc = x_shock_hllc_sim / t
        err_hllc = 100 * (V_hllc - V_shock_analytical) / V_shock_analytical
        text_info += f"HLLC solver:\n  x = {x_shock_hllc_sim:.4f}\n  V = {V_hllc:.4f} ({err_hllc:+.1f}%)"

    fig.text(0.02, 0.02, text_info, fontsize=9, family='monospace',
             verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout(rect=[0, 0.15, 1, 0.95])

    # Save figure
    output_file = OUTPUT_DIR / 'exact_vs_hllc_comparison.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison to: {output_file}")

    # Print shock analysis
    print(f"\n{'='*70}")
    print(f"  Shock Position Analysis at t = {t:.3f}")
    print(f"{'='*70}")
    print(f"  Analytical shock speed: V = {V_shock_analytical:.4f}")
    print(f"  Analytical shock position: x = {V_shock_analytical*t:.4f}")
    if x_shock_exact_sim:
        print(f"  Exact solver shock: x = {x_shock_exact_sim:.4f}, V = {x_shock_exact_sim/t:.4f} ({100*(x_shock_exact_sim/t - V_shock_analytical)/V_shock_analytical:+.1f}%)")
    if x_shock_hllc_sim:
        print(f"  HLLC solver shock: x = {x_shock_hllc_sim:.4f}, V = {x_shock_hllc_sim/t:.4f} ({100*(x_shock_hllc_sim/t - V_shock_analytical)/V_shock_analytical:+.1f}%)")

    print(f"\nDone!")

if __name__ == '__main__':
    main()
