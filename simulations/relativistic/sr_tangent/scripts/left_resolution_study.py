#!/usr/bin/env python3
"""
Left Resolution Study: Testing Kitajima's Hypothesis
Based on Kitajima et al. (2025) Section 3.2.5

From the paper (lines 621-626):
- Baseline: Both sides contain 1600 particles
- Resolution study: Right side fixed at 1600, left side varies: 1600, 3200, 51200
- Initial conditions: v^t_L = 0.9, v^t_R = 0.9
"""

import subprocess
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

# Add SRRP library path
sys.path.insert(0, str(Path("/Users/guo/Downloads/sphcode/scripts")))
try:
    from srrp import SpecialRelativisticRiemannProblem
    HAS_SRRP = True
except ImportError:
    HAS_SRRP = False
    print("Warning: SRRP library not available, analytical solution will not be shown")

# Configuration
SPHCODE_DIR = Path("/Users/guo/Downloads/sphcode")
BUILD_DIR = SPHCODE_DIR / "build"
RESULTS_DIR = SPHCODE_DIR / "simulations/relativistic/sr_tangent/results/kitajima_paper_resolution"

# Analytical values for v^t = 0.9 case
V_SHOCK = 0.4451  # Approximate shock speed
DT_OUTPUT = 0.03

# Initial conditions for v^t = 0.9 case
P_L, N_L, VX_L, VT_L = 1000.0, 1.0, 0.0, 0.9
P_R, N_R, VX_R, VT_R = 0.01, 1.0, 0.0, 0.9
GAMMA = 5.0/3.0

def get_analytical_solution(t, x_array):
    """Compute analytical solution using SRRP library."""
    if not HAS_SRRP:
        return None

    try:
        rp = SpecialRelativisticRiemannProblem(
            rho_L=N_L, P_L=P_L, vx_L=VX_L, vt_L=VT_L,
            rho_R=N_R, P_R=P_R, vx_R=VX_R, vt_R=VT_R,
            gamma=GAMMA
        )

        pres = np.zeros_like(x_array)
        dens = np.zeros_like(x_array)
        vx = np.zeros_like(x_array)
        vt = np.zeros_like(x_array)

        for i, x in enumerate(x_array):
            state = rp.sample(x, t)
            pres[i] = state['P']
            dens[i] = state['rho']
            vx[i] = state['vx']
            vt[i] = state['vt']

        return {'x': x_array, 'pres': pres, 'dens': dens, 'vx': vx, 'vt': vt}
    except Exception as e:
        print(f"Error computing analytical solution: {e}")
        return None

def create_config(N_left, N_right, solver, output_name, target_x_shock=0.2):
    """Create a config file with specified resolution and solver."""
    target_time = target_x_shock / V_SHOCK + 0.01

    config = {
        "_comment": f"Kitajima et al. (2025) - N_left={N_left}, N_right={N_right}, solver={solver}",
        "sample": "sr_tangent_velocity",
        "outputDirectory": str(RESULTS_DIR / output_name),
        "startTime": 0.0,
        "endTime": target_time + 0.01,
        "outputTime": DT_OUTPUT,
        "SPHType": "srgsph",
        "riemannSolver": solver,
        "cSpeed": 1.0,
        "cShock": 3.0,
        "gamma": GAMMA,
        "kernel": "gaussian",
        "iterativeSmoothingLength": True,
        "neighborNumber": 50,
        "N": max(N_left, N_right),
        "N_left": N_left,
        "N_right": N_right,
        "vt_left": VT_L,
        "vt_right": VT_R,
        "useGhostParticles": True,
        "ghostLayers": 6,
        "periodic": False,
        "rangeMax": [0.5],
        "rangeMin": [-0.5],
        "cflSound": 0.15,
        "cflForce": 0.1,
        "output": {"formats": [{"type": "csv", "precision": 16}]}
    }

    config_path = RESULTS_DIR / f"config_{output_name}.json"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / output_name).mkdir(parents=True, exist_ok=True)

    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)

    return config_path

def run_simulation(config_path):
    """Run the SPH simulation."""
    exe_path = BUILD_DIR / "sph"
    cmd = [str(exe_path), str(config_path)]
    print(f"  Running: {config_path.stem}...")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SPHCODE_DIR))
    if result.returncode != 0:
        print(f"  Error: {result.stderr[:200]}")
    return result.returncode == 0

def find_shock_position(x, dens):
    """Find shock position where density drops from high to low."""
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    dens_sorted = dens[sort_idx]

    right_mask = x_sorted > 0.02
    x_right = x_sorted[right_mask]
    dens_right = dens_sorted[right_mask]

    if len(dens_right) < 10:
        return None

    high_dens_mask = dens_right > 2.5
    if any(high_dens_mask):
        return x_right[high_dens_mask][-1]
    return None

def analyze_results(output_dir, target_x=0.2):
    """Analyze simulation results and find the snapshot closest to x_shock=0.2."""
    files = sorted(glob.glob(str(output_dir / 'snapshot_*.csv')))
    if not files:
        return None, None, None

    best_file = None
    best_x_shock = None
    best_t = None

    for f in files:
        snap_num = int(Path(f).stem.split('_')[-1])
        t = snap_num * DT_OUTPUT
        if t < 0.1:
            continue

        df = pd.read_csv(f, comment='#').dropna(subset=['pos_x'])
        x_shock = find_shock_position(df['pos_x'].values, df['dens'].values)

        if x_shock is not None:
            if best_x_shock is None or abs(x_shock - target_x) < abs(best_x_shock - target_x):
                best_file = f
                best_x_shock = x_shock
                best_t = t

    if best_x_shock is None:
        return None, None, None

    V_measured = best_x_shock / best_t
    error = (V_measured - V_SHOCK) / V_SHOCK * 100
    return best_t, best_x_shock, error

def load_snapshot(output_dir, snap_num=None):
    """Load a snapshot from output directory."""
    files = sorted(glob.glob(str(output_dir / 'snapshot_*.csv')))
    if not files:
        return None

    if snap_num is None:
        f = files[-1]  # Use last snapshot
    else:
        f = output_dir / f'snapshot_{snap_num:04d}.csv'
        if not Path(f).exists():
            f = files[-1]

    df = pd.read_csv(f, comment='#').dropna(subset=['pos_x'])
    return df

def plot_solver_comparison(exact_dir, hllc_dir, N_left, N_right, save_path):
    """Create comparison plot between exact and HLLC solvers with analytical solution."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Load data
    exact_df = load_snapshot(exact_dir)
    hllc_df = load_snapshot(hllc_dir)

    if exact_df is None and hllc_df is None:
        print(f"No data found for comparison plot")
        return

    # Get time from last snapshot
    files = sorted(glob.glob(str(exact_dir / 'snapshot_*.csv')))
    if files:
        snap_num = int(Path(files[-1]).stem.split('_')[-1])
        t = snap_num * DT_OUTPUT
    else:
        t = 0.45

    # Compute analytical solution
    x_anal = np.linspace(-0.4, 0.4, 500)
    anal = get_analytical_solution(t, x_anal)

    # Plot settings
    ms = 3  # marker size
    alpha = 0.6

    # === Panel 1: Pressure ===
    ax = axes[0, 0]
    if anal is not None:
        ax.plot(anal['x'], anal['pres'] / 1000, 'k-', lw=2, label='Analytical', zorder=10)
    if exact_df is not None:
        ax.scatter(exact_df['pos_x'], exact_df['pres'] / 1000, s=ms, c='blue', alpha=alpha, label='Exact Riemann')
    if hllc_df is not None:
        ax.scatter(hllc_df['pos_x'], hllc_df['pres'] / 1000, s=ms, c='red', alpha=alpha, label='HLLC')
    ax.set_ylabel('P / 1000', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.05, 1.1)
    ax.set_title('Pressure', fontsize=12)
    ax.legend(fontsize=9, loc='upper right', markerscale=3)
    ax.grid(True, alpha=0.3)

    # === Panel 2: Number Density ===
    ax = axes[0, 1]
    if anal is not None:
        ax.plot(anal['x'], anal['dens'] / 5, 'k-', lw=2, label='Analytical', zorder=10)
    if exact_df is not None:
        ax.scatter(exact_df['pos_x'], exact_df['dens'] / 5, s=ms, c='blue', alpha=alpha, label='Exact Riemann')
    if hllc_df is not None:
        ax.scatter(hllc_df['pos_x'], hllc_df['dens'] / 5, s=ms, c='red', alpha=alpha, label='HLLC')
    ax.set_ylabel('n / 5', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.05, 1.1)
    ax.set_title('Number Density', fontsize=12)
    ax.legend(fontsize=9, loc='upper right', markerscale=3)
    ax.grid(True, alpha=0.3)

    # === Panel 3: Normal Velocity ===
    ax = axes[1, 0]
    if anal is not None:
        ax.plot(anal['x'], anal['vx'], 'k-', lw=2, label='Analytical', zorder=10)
    if exact_df is not None:
        ax.scatter(exact_df['pos_x'], exact_df['vel_x'], s=ms, c='blue', alpha=alpha, label='Exact Riemann')
    if hllc_df is not None:
        ax.scatter(hllc_df['pos_x'], hllc_df['vel_x'], s=ms, c='red', alpha=alpha, label='HLLC')
    ax.set_ylabel(r'$v^x$', fontsize=12)
    ax.set_xlabel('x', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.02, 0.45)
    ax.set_title('Normal Velocity', fontsize=12)
    ax.legend(fontsize=9, loc='upper left', markerscale=3)
    ax.grid(True, alpha=0.3)

    # === Panel 4: Tangent Velocity ===
    ax = axes[1, 1]
    if anal is not None:
        ax.plot(anal['x'], anal['vt'], 'k-', lw=2, label='Analytical', zorder=10)
    if exact_df is not None and 'vel_t' in exact_df.columns:
        ax.scatter(exact_df['pos_x'], exact_df['vel_t'], s=ms, c='blue', alpha=alpha, label='Exact Riemann')
    if hllc_df is not None and 'vel_t' in hllc_df.columns:
        ax.scatter(hllc_df['pos_x'], hllc_df['vel_t'], s=ms, c='red', alpha=alpha, label='HLLC')
    ax.set_ylabel(r'$v^t$', fontsize=12)
    ax.set_xlabel('x', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(0.4, 1.02)
    ax.set_title('Tangent Velocity', fontsize=12)
    ax.legend(fontsize=9, loc='lower left', markerscale=3)
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Exact vs HLLC Riemann Solver Comparison\n'
                 f'$v^t_L = v^t_R = 0.9$, N_left={N_left}, N_right={N_right}, t={t:.3f}',
                 fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {save_path}")
    plt.close()

def plot_resolution_study(results_list, save_path):
    """Create resolution study comparison plot."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Get time from first result
    t = 0.42  # Approximate

    # Compute analytical solution
    x_anal = np.linspace(-0.4, 0.4, 500)
    anal = get_analytical_solution(t, x_anal)

    # Color map for different resolutions
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(results_list)))

    for idx, (output_dir, label, color) in enumerate(zip(
            [r['output_dir'] for r in results_list],
            [r['label'] for r in results_list],
            colors)):

        df = load_snapshot(output_dir)
        if df is None:
            continue

        ms = 2
        alpha = 0.6

        # Pressure
        axes[0, 0].scatter(df['pos_x'], df['pres'] / 1000, s=ms, c=[color], alpha=alpha, label=label)

        # Density
        axes[0, 1].scatter(df['pos_x'], df['dens'] / 5, s=ms, c=[color], alpha=alpha, label=label)

        # v^x
        axes[1, 0].scatter(df['pos_x'], df['vel_x'], s=ms, c=[color], alpha=alpha, label=label)

        # v^t
        if 'vel_t' in df.columns:
            axes[1, 1].scatter(df['pos_x'], df['vel_t'], s=ms, c=[color], alpha=alpha, label=label)

    # Add analytical solution
    if anal is not None:
        axes[0, 0].plot(anal['x'], anal['pres'] / 1000, 'k-', lw=2, label='Analytical', zorder=10)
        axes[0, 1].plot(anal['x'], anal['dens'] / 5, 'k-', lw=2, label='Analytical', zorder=10)
        axes[1, 0].plot(anal['x'], anal['vx'], 'k-', lw=2, label='Analytical', zorder=10)
        axes[1, 1].plot(anal['x'], anal['vt'], 'k-', lw=2, label='Analytical', zorder=10)

    # Formatting
    axes[0, 0].set_ylabel('P / 1000', fontsize=12)
    axes[0, 0].set_title('Pressure', fontsize=12)
    axes[0, 0].set_xlim(-0.4, 0.4)
    axes[0, 0].set_ylim(-0.05, 1.1)
    axes[0, 0].legend(fontsize=8, loc='upper right', markerscale=4)

    axes[0, 1].set_ylabel('n / 5', fontsize=12)
    axes[0, 1].set_title('Number Density', fontsize=12)
    axes[0, 1].set_xlim(-0.4, 0.4)
    axes[0, 1].set_ylim(-0.05, 1.1)
    axes[0, 1].legend(fontsize=8, loc='upper right', markerscale=4)

    axes[1, 0].set_ylabel(r'$v^x$', fontsize=12)
    axes[1, 0].set_xlabel('x', fontsize=12)
    axes[1, 0].set_title('Normal Velocity', fontsize=12)
    axes[1, 0].set_xlim(-0.4, 0.4)
    axes[1, 0].set_ylim(-0.02, 0.45)
    axes[1, 0].legend(fontsize=8, loc='upper left', markerscale=4)

    axes[1, 1].set_ylabel(r'$v^t$', fontsize=12)
    axes[1, 1].set_xlabel('x', fontsize=12)
    axes[1, 1].set_title('Tangent Velocity', fontsize=12)
    axes[1, 1].set_xlim(-0.4, 0.4)
    axes[1, 1].set_ylim(0.4, 1.02)
    axes[1, 1].legend(fontsize=8, loc='lower left', markerscale=4)

    for ax in axes.flat:
        ax.grid(True, alpha=0.3)

    plt.suptitle(f'Kitajima Resolution Study: Effect of Rarefaction Side Resolution\n'
                 f'$v^t_L = v^t_R = 0.9$, N_right fixed at 1600',
                 fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {save_path}")
    plt.close()

def main():
    print("="*70)
    print("Kitajima et al. (2025) Resolution Study")
    print("Section 3.2.5: Riemann Problem 5 - Tangential Velocity")
    print("="*70)
    print("\nFrom paper: 'particles on the right-hand side is fixed at 1,600,")
    print("while the number of particles on the left-hand side varied from")
    print("1,600, 3,200, to 51,200.'")
    print("="*70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Define test cases - both exact and HLLC solvers
    resolutions = [
        (1600, 1600, "1600+1600 (baseline)"),
        (3200, 1600, "3200+1600 (2x left)"),
        (6400, 1600, "6400+1600 (4x left)"),
        (12800, 1600, "12800+1600 (8x left)"),
    ]

    results_exact = []
    results_hllc = []

    for n_left, n_right, desc in resolutions:
        print(f"\n--- {desc} ---")

        # Run exact solver
        name_exact = f"exact_{n_left}_{n_right}"
        config_exact = create_config(n_left, n_right, "exact", name_exact)
        if run_simulation(config_exact):
            output_dir = RESULTS_DIR / name_exact
            t, x_shock, err = analyze_results(output_dir)
            results_exact.append({
                'name': name_exact,
                'label': f'Exact {desc}',
                'N_left': n_left,
                'N_right': n_right,
                'output_dir': output_dir,
                'error': err
            })
            if err is not None:
                print(f"  Exact: t={t:.3f}, x_shock={x_shock:.4f}, error={err:+.2f}%")

        # Run HLLC solver
        name_hllc = f"hllc_{n_left}_{n_right}"
        config_hllc = create_config(n_left, n_right, "hllc", name_hllc)
        if run_simulation(config_hllc):
            output_dir = RESULTS_DIR / name_hllc
            t, x_shock, err = analyze_results(output_dir)
            results_hllc.append({
                'name': name_hllc,
                'label': f'HLLC {desc}',
                'N_left': n_left,
                'N_right': n_right,
                'output_dir': output_dir,
                'error': err
            })
            if err is not None:
                print(f"  HLLC:  t={t:.3f}, x_shock={x_shock:.4f}, error={err:+.2f}%")

    # Print summary
    print("\n" + "="*70)
    print("Results Summary")
    print("="*70)
    print(f"\n{'Resolution':<20} | {'Exact Error':>12} | {'HLLC Error':>12}")
    print("-"*50)
    for exact, hllc in zip(results_exact, results_hllc):
        exact_err = f"{exact['error']:+.2f}%" if exact['error'] is not None else "N/A"
        hllc_err = f"{hllc['error']:+.2f}%" if hllc['error'] is not None else "N/A"
        print(f"{exact['N_left']}+{exact['N_right']:<14} | {exact_err:>12} | {hllc_err:>12}")

    # Create plots
    if results_exact:
        # Resolution study plot (exact solver only)
        plot_resolution_study(results_exact, RESULTS_DIR / 'resolution_study_exact.png')

        # Resolution study plot (HLLC solver only)
        plot_resolution_study(results_hllc, RESULTS_DIR / 'resolution_study_hllc.png')

        # Solver comparison at baseline resolution
        plot_solver_comparison(
            RESULTS_DIR / "exact_1600_1600",
            RESULTS_DIR / "hllc_1600_1600",
            1600, 1600,
            RESULTS_DIR / 'solver_comparison_1600.png'
        )

        # Solver comparison at high resolution
        plot_solver_comparison(
            RESULTS_DIR / "exact_12800_1600",
            RESULTS_DIR / "hllc_12800_1600",
            12800, 1600,
            RESULTS_DIR / 'solver_comparison_12800.png'
        )

        # Error vs resolution plot
        fig, ax = plt.subplots(figsize=(10, 6))

        exact_n = [r['N_left'] for r in results_exact if r['error'] is not None]
        exact_err = [r['error'] for r in results_exact if r['error'] is not None]
        hllc_n = [r['N_left'] for r in results_hllc if r['error'] is not None]
        hllc_err = [r['error'] for r in results_hllc if r['error'] is not None]

        ax.semilogx(exact_n, exact_err, 'bo-', ms=10, lw=2, label='Exact Riemann Solver')
        ax.semilogx(hllc_n, hllc_err, 'rs-', ms=10, lw=2, label='HLLC Solver')
        ax.axhline(y=0, color='g', ls='--', lw=2, label='Analytical (0%)')

        ax.set_xlabel('N_left (Rarefaction Side Particles)', fontsize=12)
        ax.set_ylabel('Shock Speed Error (%)', fontsize=12)
        ax.set_title('Shock Speed Error vs Left Resolution\n'
                    'N_right fixed at 1600, $v^t = 0.9$', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / 'error_vs_resolution.png', dpi=150, bbox_inches='tight')
        print(f"Plot saved: {RESULTS_DIR / 'error_vs_resolution.png'}")

    print(f"\nResults saved to: {RESULTS_DIR}")

if __name__ == '__main__':
    main()
