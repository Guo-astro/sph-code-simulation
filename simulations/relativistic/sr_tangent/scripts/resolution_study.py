#!/usr/bin/env python3
"""
Resolution Study: Exact vs HLLC Riemann Solvers
Kitajima et al. (2025) Problem 5 (v_t = 0.9)

Compares shock position error at different resolutions.
"""

import subprocess
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys
import shutil

# Add SRRP library path
SRRP_PATH = Path("/Users/guo/Downloads/sph-simulators/docs/papers/sg-gsph/srrp")
sys.path.insert(0, str(SRRP_PATH))
from srrp.Solver import Solver
from srrp.State import State

# Configuration
BASE_DIR = Path(__file__).parent.parent
CONFIG_DIR = BASE_DIR / 'config' / 'presets'
RESULTS_DIR = BASE_DIR / 'results'
SPHCODE_DIR = Path("/Users/guo/Downloads/sphcode")
BUILD_DIR = SPHCODE_DIR / "build"

# Analytical shock speed
V_SHOCK = 0.4451
TARGET_TIME = 0.45
DT_OUTPUT = 0.03

# Resolutions to test
RESOLUTIONS = [400, 800, 1600, 3200]

def create_config(base_config_path, N, solver_type, output_dir):
    """Create a config file with specified resolution."""
    with open(base_config_path, 'r') as f:
        config = json.load(f)

    config['N'] = N
    config['outputDirectory'] = str(output_dir)
    config['riemannSolver'] = solver_type
    config['endTime'] = TARGET_TIME + 0.01  # Slightly longer to ensure we get t=0.45

    # Save to temp config
    temp_config = RESULTS_DIR / f'temp_{solver_type}_N{N}.json'
    with open(temp_config, 'w') as f:
        json.dump(config, f, indent=2)

    return temp_config

def run_simulation(config_path):
    """Run the SPH simulation."""
    exe_path = BUILD_DIR / "sph"
    if not exe_path.exists():
        print(f"Error: {exe_path} not found")
        return False

    cmd = [str(exe_path), str(config_path)]
    print(f"  Running: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SPHCODE_DIR))
    if result.returncode != 0:
        print(f"  Error: {result.stderr[:500]}")
        return False
    return True

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

def analyze_results(output_dir, t_target=TARGET_TIME):
    """Analyze simulation results at target time."""
    snap_num = int(t_target / DT_OUTPUT)
    snapshot_file = output_dir / f'snapshot_{snap_num:05d}.csv'

    if not snapshot_file.exists():
        # Try to find closest snapshot
        files = sorted(glob.glob(str(output_dir / 'snapshot_*.csv')))
        if not files:
            return None, None
        snapshot_file = Path(files[-1])
        snap_num = int(snapshot_file.stem.split('_')[-1])

    t = snap_num * DT_OUTPUT
    df = pd.read_csv(snapshot_file, comment='#')
    df = df.dropna(subset=['pos_x'])

    x_shock = find_shock_position(df['pos_x'].values, df['dens'].values)
    if x_shock is None:
        return t, None

    V_measured = x_shock / t
    error = (V_measured - V_SHOCK) / V_SHOCK * 100
    return t, error

def main():
    print("="*70)
    print("Resolution Study: Exact vs HLLC Riemann Solvers")
    print("Kitajima et al. (2025) Problem 5 (v_t = 0.9)")
    print("="*70)

    # Base configs
    exact_base = CONFIG_DIR / 'sr_tangent_kitajima_vt09.json'
    hllc_base = CONFIG_DIR / 'sr_tangent_kitajima_vt09_hllc.json'

    if not exact_base.exists():
        print(f"Error: Base config not found: {exact_base}")
        return

    results = {'N': [], 'exact_error': [], 'hllc_error': []}

    for N in RESOLUTIONS:
        print(f"\n{'='*50}")
        print(f"Resolution N = {N}")
        print(f"{'='*50}")

        # Run exact solver
        print(f"\nRunning Exact (Iterative) solver...")
        exact_output = RESULTS_DIR / f'resolution_study' / f'exact_N{N}'
        exact_output.mkdir(parents=True, exist_ok=True)
        exact_config = create_config(exact_base, N, 'exact', exact_output)

        if run_simulation(exact_config):
            t, exact_err = analyze_results(exact_output)
            print(f"  Exact solver: error = {exact_err:+.2f}% at t = {t:.3f}")
        else:
            exact_err = None
            print(f"  Exact solver: FAILED")

        # Run HLLC solver
        print(f"\nRunning HLLC solver...")
        hllc_output = RESULTS_DIR / f'resolution_study' / f'hllc_N{N}'
        hllc_output.mkdir(parents=True, exist_ok=True)
        hllc_config = create_config(hllc_base, N, 'hllc', hllc_output)

        if run_simulation(hllc_config):
            t, hllc_err = analyze_results(hllc_output)
            print(f"  HLLC solver: error = {hllc_err:+.2f}% at t = {t:.3f}")
        else:
            hllc_err = None
            print(f"  HLLC solver: FAILED")

        results['N'].append(N)
        results['exact_error'].append(exact_err)
        results['hllc_error'].append(hllc_err)

        # Cleanup temp configs
        exact_config.unlink(missing_ok=True)
        hllc_config.unlink(missing_ok=True)

    # Plot results
    print("\n" + "="*70)
    print("Results Summary")
    print("="*70)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    valid_N = []
    valid_exact = []
    valid_hllc = []

    for i, N in enumerate(results['N']):
        if results['exact_error'][i] is not None:
            valid_N.append(N)
            valid_exact.append(results['exact_error'][i])
        if results['hllc_error'][i] is not None:
            if N not in valid_N:
                valid_N.append(N)
            valid_hllc.append(results['hllc_error'][i])

    print(f"\n{'N':>6} | {'Exact Error':>12} | {'HLLC Error':>12}")
    print("-"*36)
    for i, N in enumerate(results['N']):
        exact_str = f"{results['exact_error'][i]:+.2f}%" if results['exact_error'][i] else "N/A"
        hllc_str = f"{results['hllc_error'][i]:+.2f}%" if results['hllc_error'][i] else "N/A"
        print(f"{N:>6} | {exact_str:>12} | {hllc_str:>12}")

    # Plot
    ax.semilogx(results['N'], results['exact_error'], 'bo-', ms=10, lw=2, label='Iterative (Exact)')
    ax.semilogx(results['N'], results['hllc_error'], 'rs-', ms=10, lw=2, label='HLLC')
    ax.axhline(y=0, color='k', ls='--', alpha=0.5, label='Analytical')

    ax.set_xlabel('Number of Particles N', fontsize=12)
    ax.set_ylabel('Shock Speed Error (%)', fontsize=12)
    ax.set_title(f'Resolution Dependence: Shock Speed Error at t = {TARGET_TIME}\n'
                 f'Kitajima et al. (2025) Problem 5 ($v_t = 0.9$)', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(RESOLUTIONS)
    ax.set_xticklabels([str(N) for N in RESOLUTIONS])

    plt.tight_layout()
    output_file = RESULTS_DIR / 'resolution_study' / 'resolution_dependence.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved plot to: {output_file}")

    print("\nDone!")

if __name__ == '__main__':
    main()
