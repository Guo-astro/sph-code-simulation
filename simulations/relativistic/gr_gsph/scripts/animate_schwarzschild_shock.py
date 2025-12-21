#!/usr/bin/env python3
"""
Animate GR-GSPH Schwarzschild Radial Shock Tube with Analytical Overlay

Tests gravitational source terms in curved spacetime.
Based on Liptai & Price (2019) GRSPH benchmarks.

Key physics to observe:
1. Gravitational acceleration towards BH (inward drift of particles)
2. Asymmetric wave propagation (outgoing waves redshifted, ingoing blueshifted)
3. Comparison with flat spacetime (Minkowski) shows GR effects

The analytical overlay uses the SRRP (Special Relativistic Riemann Problem) solver
for the FLAT spacetime solution. The difference between simulation and analytical
solution shows the effect of curved spacetime (gravitational source terms).
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

# Add srrp solver to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / 'docs'))

try:
    from srrp.Solver import Solver as SRRPSolver
    from srrp.State import State
    HAS_SRRP = True

    # Pre-compute wavefan for Rosswog Test 1 (flat spacetime)
    _solver = SRRPSolver()
    _stateL = State(rho=10.0, vx=0.0, vt=0.0, pressure=40.0/3.0)
    _stateR = State(rho=1.0, vx=0.0, vt=0.0, pressure=1e-6)
    _wavefan = _solver.solve(_stateL, _stateR, 5.0/3.0)
    print(f"SRRP solver loaded, solution type: {_solver.solution_type}")
except ImportError as e:
    print(f"Warning: Could not import SRRP solver: {e}")
    print("Analytical overlay will use initial conditions only.")
    HAS_SRRP = False
    _wavefan = None


def compute_flat_sr_solution(r, t, r_disc=6.0):
    """
    Compute exact SR Riemann solution in FLAT spacetime.

    This serves as a reference to show the effect of curved spacetime.
    In Schwarzschild, waves propagate differently due to gravitational effects.

    Args:
        r: radial positions (array)
        t: time
        r_disc: initial discontinuity location

    Returns:
        n, v, P: density, velocity, pressure arrays
    """
    # Initial conditions
    n_L, v_L, P_L = 10.0, 0.0, 40.0/3.0
    n_R, v_R, P_R = 1.0, 0.0, 1e-6

    if t <= 0 or _wavefan is None:
        n = np.where(r < r_disc, n_L, n_R)
        v = np.where(r < r_disc, v_L, v_R)
        P = np.where(r < r_disc, P_L, P_R)
        return n, v, P

    # Self-similar variable: xi = (r - r_disc) / t
    # This assumes flat spacetime wave propagation
    xi = (r - r_disc) / t

    try:
        state = _wavefan.getState(xi)
        return state.rho, state.vx, state.pressure
    except:
        n = np.where(r < r_disc, n_L, n_R)
        v = np.where(r < r_disc, v_L, v_R)
        P = np.where(r < r_disc, P_L, P_R)
        return n, v, P


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


def compute_schwarzschild_lapse(r, M=1.0):
    """Compute lapse function alpha = sqrt(1 - 2M/r)."""
    alpha_sq = np.maximum(1.0 - 2.0 * M / r, 0.01)
    return np.sqrt(alpha_sq)


def make_animation(input_dir, output_file, M=1.0, gamma=5.0/3.0):
    """Create animation of Schwarzschild shock tube."""

    # Find all snapshot files
    files = sorted(glob.glob(os.path.join(input_dir, 'snapshot_*.csv')))
    if not files:
        print(f"No files found in {input_dir}")
        sys.exit(1)

    print(f"Found {len(files)} snapshots")

    # Load all snapshots
    snapshots = []
    times = []
    for f in files:
        data, t = load_csv(f)
        snapshots.append(data)
        times.append(t)

    # Setup figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('GR-GSPH Schwarzschild Radial Shock Tube (M=1)', fontsize=14, fontweight='bold')

    # Schwarzschild radius
    r_s = 2.0 * M

    # Analytical solution grid
    r_ana = np.linspace(3.0, 15.0, 500)
    r_disc = 6.0

    def update(frame):
        for ax in axes.flat:
            ax.clear()

        data = snapshots[frame]
        t = times[frame]

        # Filter out ghosts
        if 'is_ghost' in data.columns:
            mask = data['is_ghost'] == 0
        else:
            mask = np.ones(len(data), dtype=bool)

        r = data['pos_x'].values[mask]  # Radial coordinate
        n = data['dens'].values[mask]   # Rest-frame density
        v = data['vel_x'].values[mask]  # Radial velocity
        P = data['pres'].values[mask]   # Pressure

        # Compute lapse at each particle position
        alpha = compute_schwarzschild_lapse(r, M)

        # Sort by radius
        sort_idx = np.argsort(r)
        r = r[sort_idx]
        n = n[sort_idx]
        v = v[sort_idx]
        P = P[sort_idx]
        alpha = alpha[sort_idx]

        # Compute flat spacetime analytical solution (for comparison)
        n_flat, v_flat, P_flat = compute_flat_sr_solution(r_ana, t, r_disc)

        # Plot density with analytical overlay
        ax = axes[0, 0]
        ax.scatter(r, n, s=3, alpha=0.6, c='blue', label='C++ GR-GSPH (Schwarzschild)')
        ax.plot(r_ana, n_flat, 'r-', lw=2.5, alpha=0.8, label='Analytic (Flat SR)')
        ax.axvline(r_s, color='black', linestyle='--', lw=1.5, alpha=0.5, label=f'Horizon r={r_s}M')
        ax.axvline(r_disc, color='green', linestyle=':', lw=1, alpha=0.5, label='Initial disc.')
        ax.set_xlabel('r/M')
        ax.set_ylabel('Rest Density n')
        ax.set_xlim(2.5, 16)
        ax.set_ylim(0, 15)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        ax.set_title('Density')
        ax.grid(True, alpha=0.3)

        # Plot velocity with analytical overlay
        ax = axes[0, 1]
        ax.scatter(r, v, s=3, alpha=0.6, c='blue', label='C++ GR-GSPH (Schwarzschild)')
        ax.plot(r_ana, v_flat, 'r-', lw=2.5, alpha=0.8, label='Analytic (Flat SR)')
        ax.axvline(r_s, color='black', linestyle='--', lw=1.5, alpha=0.5)
        ax.axhline(0, color='k', linestyle='-', lw=0.5)
        ax.set_xlabel('r/M')
        ax.set_ylabel('Radial Velocity v')
        ax.set_xlim(2.5, 16)
        # Auto-scale y to show full velocity range (GR infall can be very negative)
        v_valid = v[np.isfinite(v)]
        v_flat_valid = v_flat[np.isfinite(v_flat)]
        if len(v_valid) > 0 and len(v_flat_valid) > 0:
            v_min = min(v_valid.min(), v_flat_valid.min()) - 0.1
            v_max = max(v_valid.max(), v_flat_valid.max()) + 0.1
        else:
            v_min, v_max = -1.0, 1.0
        ax.set_ylim(v_min, max(v_max, 1.0))
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        ax.set_title('Radial Velocity (negative = GR infall)')
        ax.grid(True, alpha=0.3)

        # Plot pressure with analytical overlay
        ax = axes[1, 0]
        ax.scatter(r, P, s=3, alpha=0.6, c='blue', label='C++ GR-GSPH (Schwarzschild)')
        ax.plot(r_ana, P_flat, 'r-', lw=2.5, alpha=0.8, label='Analytic (Flat SR)')
        ax.axvline(r_s, color='black', linestyle='--', lw=1.5, alpha=0.5)
        ax.set_xlabel('r/M')
        ax.set_ylabel('Pressure P')
        ax.set_xlim(2.5, 16)
        ax.set_ylim(-0.5, 20)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        ax.set_title('Pressure')
        ax.grid(True, alpha=0.3)

        # Plot lapse and GR effects explanation
        ax = axes[1, 1]
        r_plot = np.linspace(2.5, 16, 200)
        alpha_plot = compute_schwarzschild_lapse(r_plot, M)
        ax.plot(r_plot, alpha_plot, 'k-', lw=2, label=r'Lapse $\alpha = \sqrt{1-2M/r}$')
        ax.fill_between(r_plot, 0, alpha_plot, alpha=0.1, color='gray')
        ax.axvline(r_s, color='red', linestyle='--', lw=2, label=f'Horizon r={r_s}M')
        ax.axhline(1.0, color='blue', linestyle=':', lw=1, alpha=0.5, label='Flat spacetime')

        # Annotate GR effects
        ax.annotate('Gravitational\nredshift', xy=(4, 0.5), fontsize=9,
                   ha='center', color='darkred')
        ax.annotate('Flat spacetime\n(no redshift)', xy=(12, 0.95), fontsize=9,
                   ha='center', color='blue')

        ax.set_xlabel('r/M')
        ax.set_ylabel(r'Lapse $\alpha$ (time dilation factor)')
        ax.set_xlim(2.5, 16)
        ax.set_ylim(0, 1.15)
        ax.legend(loc='lower right', fontsize=8)
        ax.set_title('Schwarzschild Spacetime Curvature')
        ax.grid(True, alpha=0.3)

        fig.suptitle(f'GR-GSPH Schwarzschild Shock Tube vs Flat SR Solution\n'
                     f'M = 1, t = {t:.3f} (Blue = curved spacetime, Red = flat spacetime)',
                     fontsize=12, fontweight='bold')
        plt.tight_layout()

        return axes.flat

    print("Generating animation...")
    anim = FuncAnimation(fig, update, frames=len(files), interval=150, blit=False)

    # Save as GIF
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"Saved animation to {output_file}")

    # Save final frame
    final_png = output_file.replace('.gif', '_final.png')
    update(len(files) - 1)
    fig.savefig(final_png, dpi=150, bbox_inches='tight')
    print(f"Saved final frame to {final_png}")

    plt.close()


def main():
    script_dir = Path(__file__).parent
    input_dir = script_dir.parent / 'results' / 'schwarzschild_shock'
    output_file = input_dir / 'schwarzschild_shock.gif'

    print("="*70)
    print("  GR-GSPH Schwarzschild Radial Shock Tube Animation")
    print("="*70)

    if not input_dir.exists():
        print(f"Error: Results not found at {input_dir}")
        print("Run the simulation first:")
        print("  ./build/sph simulations/relativistic/gr_gsph/config/presets/gr_schwarzschild_shock.json")
        return 1

    make_animation(str(input_dir), str(output_file))
    return 0


if __name__ == '__main__':
    sys.exit(main())
