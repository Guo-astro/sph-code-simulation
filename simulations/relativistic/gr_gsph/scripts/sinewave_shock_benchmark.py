#!/usr/bin/env python3
"""
Sinusoidally Perturbed Shock Tube Benchmark

Reproduces Figure 13 from Liptai & Price (2019) MNRAS 485, 819.

This test verifies that density perturbations (sine waves) are correctly
advected through shock fronts without loss of accuracy.

Initial conditions:
    Left:  (rho, P) = (10, 40/3)          x < 0
    Right: (rho, P) = (2 + 0.3*sin(50x), 5)   x > 0

The sine wave should pass through the shock interface and be advected
with the post-shock flow while maintaining its shape.

Reference: Del Zanna & Bucciantini (2002), Rosswog (2010)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import glob
import pandas as pd

# Add SRRP solver path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / 'docs'))

try:
    from srrp.Solver import Solver as SRRPSolver
    from srrp.State import State
    HAS_SRRP = True
except ImportError:
    print("Warning: SRRP solver not available")
    HAS_SRRP = False


def compute_sine_wave_exact(x, t, gamma=5.0/3.0):
    """
    Compute the expected density profile for sinusoidally perturbed shock.

    Since the shock compresses the right state, we compute solutions for
    the density extremes (rho_min = 1.7, rho_max = 2.3) and show them as
    envelopes.
    """
    if not HAS_SRRP or t <= 0:
        # Initial conditions
        rho = np.where(x < 0, 10.0, 2.0 + 0.3 * np.sin(50 * x))
        v = np.zeros_like(x)
        P = np.where(x < 0, 40.0/3.0, 5.0)
        return rho, v, P, None, None

    solver = SRRPSolver()

    # Compute solution for density extremes
    stateL = State(rho=10.0, vx=0.0, vt=0.0, pressure=40.0/3.0)

    # Minimum density case (rho = 1.7)
    stateR_min = State(rho=1.7, vx=0.0, vt=0.0, pressure=5.0)
    wavefan_min = solver.solve(stateL, stateR_min, gamma)

    # Maximum density case (rho = 2.3)
    stateR_max = State(rho=2.3, vx=0.0, vt=0.0, pressure=5.0)
    wavefan_max = solver.solve(stateL, stateR_max, gamma)

    xi = x / t
    state_min = wavefan_min.getState(xi)
    state_max = wavefan_max.getState(xi)

    # Return main solution (mean) plus envelopes
    rho_min = state_min.rho
    rho_max = state_max.rho

    return None, None, None, rho_min, rho_max


def plot_sinewave_shock(output_dir, sim_data=None, t=0.35):
    """
    Plot sinusoidally perturbed shock tube results.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    x = np.linspace(-0.5, 0.5, 1000)

    # Get exact solution envelopes
    _, _, _, rho_min, rho_max = compute_sine_wave_exact(x, t)

    # Initial sine wave profile (for reference)
    x_right = x[x > 0]
    rho_init_right = 2.0 + 0.3 * np.sin(50 * x_right)

    # Left panel: Full shock tube
    ax = axes[0]

    if rho_min is not None:
        ax.fill_between(x, rho_min, rho_max, alpha=0.3, color='green',
                        label='Analytic envelope (min/max)')
        ax.plot(x, rho_min, 'g--', lw=1.5, alpha=0.7)
        ax.plot(x, rho_max, 'g--', lw=1.5, alpha=0.7)

    if sim_data is not None:
        ax.scatter(sim_data['x'], sim_data['rho'], s=3, c='blue',
                   alpha=0.7, label='C++ SPH Simulation')

    ax.axvline(0, color='k', ls=':', lw=1, alpha=0.5, label='Initial interface')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'Sinusoidally Perturbed Shock Tube (t = {t})', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.55, 0.55)
    ax.set_ylim(0, 12)

    # Right panel: Zoom on sine wave region
    ax = axes[1]

    if rho_min is not None:
        # Focus on region where sine wave exists
        mask = (x > 0.1) & (x < 0.4)
        ax.fill_between(x[mask], rho_min[mask], rho_max[mask], alpha=0.3,
                        color='green', label='Analytic envelope')
        ax.plot(x[mask], rho_min[mask], 'g--', lw=1.5, alpha=0.7)
        ax.plot(x[mask], rho_max[mask], 'g--', lw=1.5, alpha=0.7)

    if sim_data is not None:
        mask = (sim_data['x'] > 0.1) & (sim_data['x'] < 0.4)
        ax.scatter(sim_data['x'][mask], sim_data['rho'][mask], s=10, c='blue',
                   alpha=0.7, label='C++ SPH Simulation')

    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Zoom: Advected Sine Wave', fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    fig.suptitle('Sinusoidally Perturbed SR Shock: Liptai & Price (2019) Fig 13\n'
                 'Initial: Left (10, 40/3), Right (2+0.3sin(50x), 5)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()

    output_file = output_dir / 'sinewave_shock.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_initial_conditions(output_dir):
    """Plot the initial conditions for the sine wave test."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    x = np.linspace(-0.5, 0.5, 1000)

    # Initial density
    rho = np.where(x < 0, 10.0, 2.0 + 0.3 * np.sin(50 * x))
    P = np.where(x < 0, 40.0/3.0, 5.0)

    ax = axes[0]
    ax.plot(x, rho, 'b-', lw=2, label='Initial density')
    ax.axhline(10, color='r', ls='--', lw=1, alpha=0.5, label='Left state')
    ax.axhline(2.0, color='g', ls='--', lw=1, alpha=0.5, label='Right mean')
    ax.axvline(0, color='k', ls=':', lw=1.5, alpha=0.7, label='Interface')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Initial Density Profile', fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.55, 0.55)
    ax.set_ylim(0, 12)

    ax = axes[1]
    ax.plot(x, P, 'b-', lw=2, label='Initial pressure')
    ax.axvline(0, color='k', ls=':', lw=1.5, alpha=0.7, label='Interface')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('Pressure', fontsize=12)
    ax.set_title('Initial Pressure Profile', fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.55, 0.55)

    fig.suptitle('Sinusoidally Perturbed Shock Tube: Initial Conditions\n'
                 'Left: (rho, P) = (10, 40/3), Right: (rho, P) = (2 + 0.3sin(50x), 5)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()

    output_file = output_dir / 'sinewave_initial.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def load_simulation_data(sim_dir, target_time=0.35):
    """
    Load C++ simulation snapshot closest to target time.

    Returns dict with 'x', 'rho', 'v', 'P' arrays, or None if not found.
    """
    if not sim_dir.exists():
        print(f"  Simulation directory not found: {sim_dir}")
        return None

    # Find all snapshot files
    files = sorted(glob.glob(str(sim_dir / 'snapshot_*.csv')))
    if not files:
        print(f"  No snapshot files found in {sim_dir}")
        return None

    # Find snapshot closest to target time
    best_file = None
    best_time_diff = float('inf')
    best_time = 0

    for f in files:
        # Read time from header
        with open(f, 'r') as fp:
            for line in fp:
                if '# Time (code):' in line:
                    try:
                        t = float(line.split(':')[1].strip())
                        if abs(t - target_time) < best_time_diff:
                            best_time_diff = abs(t - target_time)
                            best_file = f
                            best_time = t
                    except:
                        pass
                    break

    if best_file is None:
        print(f"  Could not find suitable snapshot")
        return None

    print(f"  Loading: {Path(best_file).name} (t = {best_time:.4f})")

    # Load data
    data = pd.read_csv(best_file, comment='#')

    # Filter out ghost particles
    if 'is_ghost' in data.columns:
        mask = data['is_ghost'] == 0
        data = data[mask]

    # Extract arrays
    result = {
        'x': data['pos_x'].values,
        'time': best_time
    }

    # Map column names - rest-frame density
    if 'dens' in data.columns:
        # For SR-GSPH, dens is the lab-frame density N = gamma * n
        # Need to recover rest-frame density n = N / gamma
        if 'gamma_lor' in data.columns:
            result['rho'] = data['dens'].values / data['gamma_lor'].values
        else:
            result['rho'] = data['dens'].values
    elif 'rho' in data.columns:
        result['rho'] = data['rho'].values

    if 'vel_x' in data.columns:
        result['v'] = data['vel_x'].values
    if 'pres' in data.columns:
        result['P'] = data['pres'].values

    return result


def main():
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("  Sinusoidally Perturbed Shock Tube Benchmark")
    print("  Based on Liptai & Price (2019) MNRAS 485, 819 - Fig 13")
    print("=" * 70)

    # Try to find C++ simulation results
    project_root = script_dir.parent.parent.parent.parent
    sim_dir = project_root / 'simulations' / 'relativistic' / 'sr_sod' / 'results' / 'sinewave_shock'

    print("\n[1] Loading C++ Simulation Data")
    sim_data = load_simulation_data(sim_dir, target_time=0.35)

    if sim_data is None:
        print("  No simulation data found. Run the simulation first:")
        print("  ./build/sph simulations/relativistic/sr_sod/config/presets/sr_sinewave_shock.json")
    else:
        print(f"  Loaded {len(sim_data['x'])} particles")

    print("\n[2] Initial Conditions")
    plot_initial_conditions(output_dir)

    print("\n[3] Shock Tube Solution (C++ + analytic envelopes)")
    plot_sinewave_shock(output_dir, sim_data=sim_data, t=sim_data['time'] if sim_data else 0.35)

    print("\n" + "=" * 70)
    print("  SINE WAVE BENCHMARK COMPLETE")
    print("=" * 70)

    if sim_data is not None:
        print(f"\nBlue points: C++ SR-GSPH simulation ({len(sim_data['x'])} particles)")
        print("Green envelope: Analytic solution for density extremes (rho=1.7 and 2.3)")

    print("""
NOTES:
======
This test verifies that density perturbations are correctly advected
through shock fronts.

Initial conditions:
  Left:  (rho, P) = (10, 40/3)             for x < 0
  Right: (rho, P) = (2 + 0.3*sin(50x), 5)  for x > 0

The sine wave should:
1. Pass through the shock interface
2. Be advected with post-shock flow
3. Maintain its shape (no loss in accuracy)

This is a stringent test for numerical diffusion and shock capturing.
""")


if __name__ == '__main__':
    main()
