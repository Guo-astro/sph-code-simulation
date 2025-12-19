#!/usr/bin/env python3
"""
Animate SR-GSPH Riemann Problem 5 (Kitajima et al. 2025) with Exact Solution from SRRP

v_t = 0.9 case (Figure 12)

Initial conditions:
  (P_L, n_L, v^x_L, v^t_L) = (1000.0, 1.0, 0, 0.9)
  (P_R, n_R, v^x_R, v^t_R) = (0.01, 1.0, 0, 0.9)

Uses SRRP (Special Relativistic Riemann Problem) library for exact solution
with proper characteristic-based rarefaction fan computation.
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

# Add SRRP library to path
SRRP_PATH = Path("/Users/guo/Downloads/sph-simulators/docs/papers/sg-gsph/srrp")
sys.path.insert(0, str(SRRP_PATH))

from srrp.Solver import Solver
from srrp.State import State

# Configuration
RESULTS_DIR = Path(__file__).parent.parent / "results" / "kitajima_vt09"
OUTPUT_DIR = RESULTS_DIR

# Physical parameters
GAMMA_EOS = 5.0/3.0

# Initial conditions for v_t = 0.9
P_L = 1000.0
P_R = 0.01
RHO_L = 1.0
RHO_R = 1.0
VX_L = 0.0
VX_R = 0.0
VT_L = 0.9
VT_R = 0.9

# Contact discontinuity initial position
X_CONTACT = 0.0

# Solve the Riemann problem once at module load
_solver = Solver()
_stateL = State(rho=RHO_L, vx=VX_L, vt=VT_L, pressure=P_L)
_stateR = State(rho=RHO_R, vx=VX_R, vt=VT_R, pressure=P_R)
_wavefan = _solver.solve(_stateL, _stateR, GAMMA_EOS)

# Extract exact values for display
_states = _wavefan.states
EXACT_SOLUTION = {
    'P_star': _states[1].pressure,
    'vx_star': _states[1].vx,
    'vt_L_star': _states[1].vt,  # Left of contact
    'vt_R_star': _states[2].vt,  # Right of contact (different due to shock!)
    'rho_L_prime': _states[1].rho,
    'rho_R_prime': _states[2].rho,
    'solution_type': _solver.solution_type,
}

print(f"SRRP Solution type: {EXACT_SOLUTION['solution_type']}")
print(f"  P*       = {EXACT_SOLUTION['P_star']:.4f}")
print(f"  vx*      = {EXACT_SOLUTION['vx_star']:.4f}")
print(f"  vt_L*    = {EXACT_SOLUTION['vt_L_star']:.4f}")
print(f"  vt_R*    = {EXACT_SOLUTION['vt_R_star']:.4f}")
print(f"  rho_L'   = {EXACT_SOLUTION['rho_L_prime']:.4f}")
print(f"  rho_R'   = {EXACT_SOLUTION['rho_R_prime']:.4f}")


def compute_exact_solution(x_array, t):
    """
    Compute exact solution at given positions and time using SRRP library.
    Uses self-similar coordinate xi = x/t and the wavefan structure.
    """
    if t <= 0:
        rho = np.where(x_array < X_CONTACT, RHO_L, RHO_R)
        P = np.where(x_array < X_CONTACT, P_L, P_R)
        vx = np.zeros_like(x_array)
        vt = np.full_like(x_array, VT_L)
        return rho, P, vx, vt

    # Self-similar variable
    xi = (x_array - X_CONTACT) / t

    # Get state at each xi from SRRP wavefan
    state = _wavefan.getState(xi)

    return state.rho, state.pressure, state.vx, state.vt


def read_csv_snapshot(filepath):
    """Read a CSV snapshot file."""
    data = {'time': 0.0, 'step': 0}

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('# Time (code):'):
                data['time'] = float(line.split(':')[1].strip())
            elif line.startswith('# Step:'):
                data['step'] = int(line.split(':')[1].strip())
            elif not line.startswith('#'):
                break

    header = ['id', 'pos_x', 'vel_x', 'vel_t', 'acc_x', 'mass', 'dens',
             'pres', 'ene', 'sml', 'sound', 'alpha', 'balsara',
             'gradh', 'phi', 'grav_acc_x', 'neighbor', 'is_ghost']

    df_data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = line.strip().split(',')
            if values and values[0]:
                try:
                    df_data.append([float(v) if v else 0.0 for v in values])
                except ValueError:
                    continue

    if not df_data:
        return None

    arr = np.array(df_data)
    result = {}
    for i, col in enumerate(header):
        if i < arr.shape[1]:
            result[col] = arr[:, i]

    result['time'] = data['time']
    result['step'] = data['step']
    return result


def find_shock_position(x, pres, threshold=0.5):
    """Find approximate shock position from pressure jump."""
    for i in range(len(x) - 1):
        if pres[i] < threshold and pres[i+1] > threshold:
            return x[i]
    return None


def create_animation():
    """Create animation with exact solution overlay."""

    snapshot_files = sorted(glob.glob(str(RESULTS_DIR / "snapshot_*.csv")))
    if not snapshot_files:
        print(f"No snapshots found in {RESULTS_DIR}")
        return

    print(f"Found {len(snapshot_files)} snapshots")

    snapshots = []
    for f in snapshot_files:
        data = read_csv_snapshot(f)
        if data is not None:
            snapshots.append(data)

    n_snapshots = len(snapshots)
    if n_snapshots == 0:
        print("No valid snapshots found")
        return

    print(f"\nCreating animation with {n_snapshots} frames...")

    # Set up figure - Kitajima Figure 12 style (linear scale, P/1000, n/5)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(r'SR-GSPH Riemann Problem 5 (Kitajima et al. 2025, Fig 12)' + '\n' +
                 r'$v^t_L = v^t_R = 0.9$, $P_L = 1000$, $P_R = 0.01$',
                 fontsize=14, fontweight='bold')

    # Initialize lines
    num_lines = []
    for ax in axes.flatten():
        line, = ax.plot([], [], 'b.', markersize=2, alpha=0.6, label='Numerical (SRGSPH)')
        num_lines.append(line)

    exact_lines = []
    for ax in axes.flatten():
        line, = ax.plot([], [], 'r-', linewidth=2, alpha=0.8, label='Exact (Pons2000)')
        exact_lines.append(line)

    # Configure axes - LINEAR scale like Kitajima Figure 12
    # Density (n/5 units in Kitajima)
    axes[0, 0].set_ylabel(r'Rest Density $n$', fontsize=12)
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_xlim(-0.5, 0.5)
    axes[0, 0].set_ylim(0, 6)  # Linear scale
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title('Rest-Frame Density')
    axes[0, 0].axhline(EXACT_SOLUTION['rho_R_prime'], color='green', linestyle='--',
                       alpha=0.5, label=f"$n'_R$={EXACT_SOLUTION['rho_R_prime']:.2f}")
    axes[0, 0].legend(loc='upper right', fontsize=9)

    # Pressure (P/1000 units in Kitajima, so show full scale)
    axes[0, 1].set_ylabel(r'Pressure $P$', fontsize=12)
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_xlim(-0.5, 0.5)
    axes[0, 1].set_ylim(0, 1100)  # Linear scale
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_title('Pressure')
    axes[0, 1].axhline(EXACT_SOLUTION['P_star'], color='green', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$P^*$={EXACT_SOLUTION['P_star']:.3f}")
    axes[0, 1].legend(loc='upper right', fontsize=9)

    # Normal velocity
    axes[1, 0].set_ylabel(r'Normal Velocity $v^x/c$', fontsize=12)
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_xlim(-0.5, 0.5)
    axes[1, 0].set_ylim(-0.1, 0.5)
    axes[1, 0].axhline(0, color='k', linestyle='--', linewidth=0.5)
    axes[1, 0].axhline(EXACT_SOLUTION['vx_star'], color='green', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$v^x_*$={EXACT_SOLUTION['vx_star']:.3f}")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_title('Normal Velocity')
    axes[1, 0].legend(loc='upper left', fontsize=9)

    # Tangent velocity
    axes[1, 1].set_ylabel(r'Tangent Velocity $v^t/c$', fontsize=12)
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_xlim(-0.5, 0.5)
    axes[1, 1].set_ylim(0.7, 1.0)  # Adjusted to show vt_R* = 0.77
    axes[1, 1].axhline(VT_L, color='orange', linestyle=':', linewidth=1,
                       alpha=0.7, label=r'$v^t_{init}=0.9$')
    axes[1, 1].axhline(EXACT_SOLUTION['vt_L_star'], color='green', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$v^t_{{L*}}$={EXACT_SOLUTION['vt_L_star']:.3f}")
    axes[1, 1].axhline(EXACT_SOLUTION['vt_R_star'], color='purple', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$v^t_{{R*}}$={EXACT_SOLUTION['vt_R_star']:.3f}")
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_title('Tangent Velocity')
    axes[1, 1].legend(loc='lower right', fontsize=9)

    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12, fontweight='bold',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    x_exact = np.linspace(-0.5, 0.5, 1000)

    def init():
        for line in num_lines + exact_lines:
            line.set_data([], [])
        time_text.set_text('')
        return num_lines + exact_lines + [time_text]

    def animate(frame):
        snap = snapshots[frame]
        t = snap['time']

        mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
        x = snap['pos_x'][mask]
        dens = snap['dens'][mask]
        pres = snap['pres'][mask]
        vel_x = snap['vel_x'][mask]
        vel_t = snap['vel_t'][mask]

        sort_idx = np.argsort(x)
        x = x[sort_idx]
        dens = dens[sort_idx]
        pres = pres[sort_idx]
        vel_x = vel_x[sort_idx]
        vel_t = vel_t[sort_idx]

        rho_ex, P_ex, vx_ex, vt_ex = compute_exact_solution(x_exact, t)

        num_lines[0].set_data(x, dens)
        num_lines[1].set_data(x, pres)
        num_lines[2].set_data(x, vel_x)
        num_lines[3].set_data(x, vel_t)

        exact_lines[0].set_data(x_exact, rho_ex)
        exact_lines[1].set_data(x_exact, P_ex)
        exact_lines[2].set_data(x_exact, vx_ex)
        exact_lines[3].set_data(x_exact, vt_ex)

        shock_pos = find_shock_position(x, pres)
        shock_str = f", shock at x={shock_pos:.3f}" if shock_pos else ""
        time_text.set_text(f't = {t:.4f} (frame {frame+1}/{n_snapshots}){shock_str}')

        return num_lines + exact_lines + [time_text]

    anim = FuncAnimation(fig, animate, init_func=init, frames=n_snapshots,
                         interval=200, blit=True)

    output_path = OUTPUT_DIR / "kitajima_vt09_with_exact.gif"
    print(f"Saving animation to {output_path}...")
    anim.save(output_path, writer=PillowWriter(fps=5))
    print("Animation saved!")

    plt.close(fig)

    # Create final frame when shock is at x=0.2
    for snap in snapshots:
        mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
        x = snap['pos_x'][mask]
        pres = snap['pres'][mask]
        sort_idx = np.argsort(x)
        shock_pos = find_shock_position(x[sort_idx], pres[sort_idx])
        if shock_pos and abs(shock_pos - 0.2) < 0.02:
            create_final_plot(snap, x_exact)
            break


def create_final_plot(snap, x_exact):
    """Create static plot at specific time."""

    t = snap['time']

    mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
    x = snap['pos_x'][mask]
    dens = snap['dens'][mask]
    pres = snap['pres'][mask]
    vel_x = snap['vel_x'][mask]
    vel_t = snap['vel_t'][mask]

    sort_idx = np.argsort(x)
    x = x[sort_idx]
    dens = dens[sort_idx]
    pres = pres[sort_idx]
    vel_x = vel_x[sort_idx]
    vel_t = vel_t[sort_idx]

    rho_ex, P_ex, vx_ex, vt_ex = compute_exact_solution(x_exact, t)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'SR-GSPH Riemann Problem 5 at t = {t:.4f} (shock at x=0.2)\n' +
                 r'Kitajima et al. (2025) Figure 12: $v^t_L = v^t_R = 0.9$',
                 fontsize=14, fontweight='bold')

    # Density
    axes[0, 0].plot(x, dens, 'b.', markersize=3, alpha=0.6, label='Numerical')
    axes[0, 0].plot(x_exact, rho_ex, 'r-', linewidth=2, alpha=0.8, label='Exact')
    axes[0, 0].axhline(EXACT_SOLUTION['rho_R_prime'], color='green', linestyle='--',
                       alpha=0.5, label=f"$n'_R$={EXACT_SOLUTION['rho_R_prime']:.2f}")
    axes[0, 0].set_ylabel(r'Rest Density $n$', fontsize=12)
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_xlim(-0.5, 0.5)
    axes[0, 0].set_ylim(0, 6)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(loc='upper right', fontsize=9)
    axes[0, 0].set_title('Rest-Frame Density')

    # Pressure
    axes[0, 1].plot(x, pres, 'b.', markersize=3, alpha=0.6, label='Numerical')
    axes[0, 1].plot(x_exact, P_ex, 'r-', linewidth=2, alpha=0.8, label='Exact')
    axes[0, 1].axhline(EXACT_SOLUTION['P_star'], color='green', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$P^*$={EXACT_SOLUTION['P_star']:.3f}")
    axes[0, 1].set_ylabel(r'Pressure $P$', fontsize=12)
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_xlim(-0.5, 0.5)
    axes[0, 1].set_ylim(0, 1100)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(loc='upper right', fontsize=9)
    axes[0, 1].set_title('Pressure')

    # Normal velocity
    axes[1, 0].plot(x, vel_x, 'b.', markersize=3, alpha=0.6, label='Numerical')
    axes[1, 0].plot(x_exact, vx_ex, 'r-', linewidth=2, alpha=0.8, label='Exact')
    axes[1, 0].axhline(0, color='k', linestyle='--', linewidth=0.5)
    axes[1, 0].axhline(EXACT_SOLUTION['vx_star'], color='green', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$v^x_*$={EXACT_SOLUTION['vx_star']:.3f}")
    axes[1, 0].set_ylabel(r'Normal Velocity $v^x/c$', fontsize=12)
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_xlim(-0.5, 0.5)
    axes[1, 0].set_ylim(-0.1, 0.5)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(loc='upper left', fontsize=9)
    axes[1, 0].set_title('Normal Velocity')

    # Tangent velocity
    axes[1, 1].plot(x, vel_t, 'b.', markersize=3, alpha=0.6, label='Numerical')
    axes[1, 1].plot(x_exact, vt_ex, 'r-', linewidth=2, alpha=0.8, label='Exact')
    axes[1, 1].axhline(VT_L, color='orange', linestyle=':', linewidth=1,
                       alpha=0.7, label=r'$v^t_{init}=0.9$')
    axes[1, 1].axhline(EXACT_SOLUTION['vt_L_star'], color='green', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$v^t_{{L*}}$={EXACT_SOLUTION['vt_L_star']:.3f}")
    axes[1, 1].axhline(EXACT_SOLUTION['vt_R_star'], color='purple', linestyle='--',
                       linewidth=1.5, alpha=0.7, label=f"$v^t_{{R*}}$={EXACT_SOLUTION['vt_R_star']:.3f}")
    axes[1, 1].set_ylabel(r'Tangent Velocity $v^t/c$', fontsize=12)
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_xlim(-0.5, 0.5)
    axes[1, 1].set_ylim(0.7, 1.0)  # Adjusted to show vt_R* = 0.77
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(loc='lower right', fontsize=9)
    axes[1, 1].set_title('Tangent Velocity')

    plt.tight_layout()

    output_path = OUTPUT_DIR / "kitajima_figure12_final.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved final plot to {output_path}")

    plt.close(fig)


if __name__ == "__main__":
    print("="*70)
    print("  SR-GSPH Riemann Problem 5 Animation (v_t = 0.9)")
    print("  Kitajima et al. (2025) Section 3.2.1 Figure 12")
    print("="*70)
    print(f"\nExact solution values (SRRP solver):")
    print(f"  Solution type: {EXACT_SOLUTION['solution_type']}")
    print(f"  P*       = {EXACT_SOLUTION['P_star']:.4f}")
    print(f"  v_x*     = {EXACT_SOLUTION['vx_star']:.4f}")
    print(f"  v_t_L*   = {EXACT_SOLUTION['vt_L_star']:.4f}")
    print(f"  v_t_R*   = {EXACT_SOLUTION['vt_R_star']:.4f}")
    print(f"  n'_L     = {EXACT_SOLUTION['rho_L_prime']:.4f}")
    print(f"  n'_R     = {EXACT_SOLUTION['rho_R_prime']:.4f}")
    print("="*70)

    if not RESULTS_DIR.exists():
        print(f"Error: Results directory not found: {RESULTS_DIR}")
        print("Run the simulation first:")
        print("  ./build/sph simulations/relativistic/sr_tangent/config/presets/sr_tangent_kitajima_vt09.json")
        sys.exit(1)

    create_animation()
    print("\nDone!")
