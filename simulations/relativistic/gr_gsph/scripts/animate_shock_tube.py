#!/usr/bin/env python3
"""
GR-GSPH Shock Tube Animation with Analytical Solution Overlay

Uses the SRRP (Special Relativistic Riemann Problem) solver for exact solution.
For Minkowski (flat) spacetime, GR-GSPH should match SR-GSPH exactly.

Based on Rosswog (2010) Test 1: Standard Relativistic Shock Tube
- Left:  (n, v, P) = (10, 0, 40/3)
- Right: (n, v, P) = (1, 0, 10^-6)
- gamma = 5/3

Author: Generated for GR-GSPH validation
"""

import argparse
import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Add srrp solver to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'docs'))

try:
    from srrp.Solver import Solver as SRRPSolver
    from srrp.State import State
    HAS_SRRP = True
except ImportError as e:
    print(f"Warning: Could not import SRRP solver: {e}")
    HAS_SRRP = False

# Pre-compute wavefan for the Rosswog Test 1 at module load
_wavefan = None
_solver = None
if HAS_SRRP:
    try:
        # Rosswog Test 1 initial conditions
        _solver = SRRPSolver()
        _stateL = State(rho=10.0, vx=0.0, vt=0.0, pressure=40.0/3.0)
        _stateR = State(rho=1.0, vx=0.0, vt=0.0, pressure=1e-6)
        _wavefan = _solver.solve(_stateL, _stateR, 5.0/3.0)
        print(f"SRRP solution type: {_solver.solution_type}")
    except Exception as e:
        print(f"Warning: SRRP solver initialization failed: {e}")
        _wavefan = None


def load_csv(filename):
    """Load CSV snapshot with header comments."""
    import pandas as pd
    data = pd.read_csv(filename, comment='#')
    return data


def get_time_from_csv(filename):
    """Extract time from CSV header."""
    with open(filename, 'r') as f:
        for line in f:
            if '# time =' in line.lower() or '# Time (code):' in line:
                try:
                    return float(line.split('=')[1].strip().split()[0])
                except:
                    try:
                        return float(line.split(':')[1].strip())
                    except:
                        pass
    return 0.0


def compute_exact_sr_solution(x, t, gamma=5.0/3.0, test_type='rosswog1'):
    """
    Compute exact SR Riemann solution using SRRP solver.

    Test types:
    - 'rosswog1': Left (n=10, v=0, P=40/3), Right (n=1, v=0, P=1e-6)
    - 'sod': Left (n=1, v=0, P=1), Right (n=0.125, v=0, P=0.1)
    """
    # Setup initial conditions
    if test_type == 'rosswog1':
        n_L, v_L, P_L = 10.0, 0.0, 40.0/3.0
        n_R, v_R, P_R = 1.0, 0.0, 1e-6
    else:
        n_L, v_L, P_L = 1.0, 0.0, 1.0
        n_R, v_R, P_R = 0.125, 0.0, 0.1

    if t <= 0:
        n = np.where(x < 0, n_L, n_R)
        v = np.where(x < 0, v_L, v_R)
        P = np.where(x < 0, P_L, P_R)
        return n, v, P

    # Use pre-computed wavefan for rosswog1 test
    if test_type == 'rosswog1' and _wavefan is not None:
        # Self-similar variable xi = x/t
        xi = x / t
        state = _wavefan.getState(xi)
        return state.rho, state.vx, state.pressure

    # Fallback or other test types: solve on the fly
    if not HAS_SRRP or _solver is None:
        n = np.where(x < 0, n_L, n_R)
        v = np.where(x < 0, v_L, v_R)
        P = np.where(x < 0, P_L, P_R)
        return n, v, P

    try:
        solver = SRRPSolver()
        stateL = State(rho=n_L, vx=v_L, vt=0.0, pressure=P_L)
        stateR = State(rho=n_R, vx=v_R, vt=0.0, pressure=P_R)
        wavefan = solver.solve(stateL, stateR, gamma)

        # Self-similar variable
        xi = x / t
        state = wavefan.getState(xi)
        return state.rho, state.vx, state.pressure

    except Exception as e:
        print(f"SRRP solver failed: {e}")
        n = np.where(x < 0, n_L, n_R)
        v = np.where(x < 0, v_L, v_R)
        P = np.where(x < 0, P_L, P_R)
        return n, v, P


def make_animation(input_dir, output_file, gamma=5.0/3.0, test_type='rosswog1'):
    """Create animation with analytical overlay."""

    # Find all snapshot files
    pattern = 'snapshot_*.csv'
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        print(f"No files found matching {pattern}")
        sys.exit(1)

    print(f"Found {len(files)} snapshots")
    print(f"SRRP solver available: {HAS_SRRP}")

    # Get times
    times = [get_time_from_csv(f) for f in files]

    # Setup figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Determine axis limits based on test type
    if test_type == 'rosswog1':
        n_max = 12.0
        P_max = 16.0
        v_max = 0.9
        u_max = 5.0
    else:
        n_max = 1.2
        P_max = 1.2
        v_max = 1.0
        u_max = 3.0

    def update(frame):
        for ax in axes.flat:
            ax.clear()

        data = load_csv(files[frame])
        x = data['pos_x'].values
        t = times[frame]

        # Filter out ghost particles if column exists
        if 'is_ghost' in data.columns:
            mask = data['is_ghost'] == 0
            x = x[mask]
            n_sim = data['dens'].values[mask]
            v_sim = data['vel_x'].values[mask]
            P_sim = data['pres'].values[mask]
        else:
            n_sim = data['dens'].values
            v_sim = data['vel_x'].values
            P_sim = data['pres'].values

        # Compute analytical solution
        x_ana = np.linspace(-0.5, 0.5, 500)
        n_ana, v_ana, P_ana = compute_exact_sr_solution(x_ana, t, gamma, test_type)

        # Plot density
        ax = axes[0, 0]
        ax.scatter(x, n_sim, s=2, alpha=0.6, c='blue', label='GR-GSPH')
        ax.plot(x_ana, n_ana, 'r-', lw=1.5, label='Exact SR')
        ax.set_xlabel('x')
        ax.set_ylabel('Density n')
        ax.set_xlim(-0.55, 0.55)
        ax.set_ylim(0, n_max)
        ax.legend(loc='upper right')
        ax.set_title('Density')
        ax.grid(True, alpha=0.3)

        # Plot velocity
        ax = axes[0, 1]
        ax.scatter(x, v_sim, s=2, alpha=0.6, c='blue', label='GR-GSPH')
        ax.plot(x_ana, v_ana, 'r-', lw=1.5, label='Exact SR')
        ax.set_xlabel('x')
        ax.set_ylabel('Velocity v')
        ax.set_xlim(-0.55, 0.55)
        ax.set_ylim(-0.1, v_max)
        ax.legend(loc='upper left')
        ax.set_title('Velocity')
        ax.grid(True, alpha=0.3)

        # Plot pressure
        ax = axes[1, 0]
        ax.scatter(x, P_sim, s=2, alpha=0.6, c='blue', label='GR-GSPH')
        ax.plot(x_ana, P_ana, 'r-', lw=1.5, label='Exact SR')
        ax.set_xlabel('x')
        ax.set_ylabel('Pressure P')
        ax.set_xlim(-0.55, 0.55)
        ax.set_ylim(-0.5, P_max)
        ax.legend(loc='upper right')
        ax.set_title('Pressure')
        ax.grid(True, alpha=0.3)

        # Plot internal energy
        ax = axes[1, 1]
        u_sim = P_sim / ((gamma - 1.0) * np.maximum(n_sim, 1e-15))
        u_ana = P_ana / ((gamma - 1.0) * np.maximum(n_ana, 1e-15))
        ax.scatter(x, u_sim, s=2, alpha=0.6, c='blue', label='GR-GSPH')
        ax.plot(x_ana, u_ana, 'r-', lw=1.5, label='Exact SR')
        ax.set_xlabel('x')
        ax.set_ylabel('Internal Energy u')
        ax.set_xlim(-0.55, 0.55)
        ax.set_ylim(-0.5, u_max)
        ax.legend(loc='upper right')
        ax.set_title('Internal Energy')
        ax.grid(True, alpha=0.3)

        fig.suptitle(f'GR-GSPH Shock Tube (Minkowski Spacetime) - t = {t:.4f}', fontsize=14)
        plt.tight_layout()

        return axes.flat

    print("Generating animation...")
    anim = FuncAnimation(fig, update, frames=len(files),
                         interval=150, blit=False)

    # Save as GIF
    writer = PillowWriter(fps=7)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"Saved animation to {output_file}")

    # Also save final frame as PNG
    final_png = output_file.replace('.gif', '_final.png')
    update(len(files) - 1)
    fig.savefig(final_png, dpi=150)
    print(f"Saved final frame to {final_png}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Animate GR-GSPH shock tube')
    parser.add_argument('--input', '-i',
                        default='../results/minkowski_shock_tube',
                        help='Input directory with CSV snapshots')
    parser.add_argument('--output', '-o',
                        default='gr_gsph_shock_tube.gif',
                        help='Output animation file')
    parser.add_argument('--gamma', type=float, default=5.0/3.0,
                        help='Adiabatic index')
    parser.add_argument('--test', choices=['rosswog1', 'sod'], default='rosswog1',
                        help='Test case type')

    args = parser.parse_args()

    make_animation(args.input, args.output, args.gamma, args.test)


if __name__ == '__main__':
    main()
