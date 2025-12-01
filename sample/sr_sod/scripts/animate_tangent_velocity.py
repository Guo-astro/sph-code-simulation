#!/usr/bin/env python3
"""
Animate SR Tangent Velocity Shock Tube (Kitajima et al. 2025 Section 3.1.5)

This test demonstrates the effect of tangential velocity on relativistic shocks.
The key physics: v^t affects the solution through the Lorentz factor:
  γ = 1/√(1 - (v_x² + v_t²)/c²)

Creates GIF/MP4 animation showing density, pressure, v_x, and v_t profiles
with exact Riemann solution overlay.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import sys
import os
import glob
import re

# Add docs directory to path for relativistic_riemann_solver_tangent
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'docs'))
from relativistic_riemann_solver_tangent import RelativisticRiemannSolverTangent


def load_snapshot(filename):
    """Load CSV snapshot and extract particle data (excluding ghost particles)"""
    data = []
    time = 0.0
    col_indices = {}
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if '# Time (code):' in line:
                    time = float(line.split(':')[1].strip())
                continue
            parts = line.strip().split(',')
            # Parse header to get column indices
            if parts[0] == 'id':
                for i, name in enumerate(parts):
                    col_indices[name.strip()] = i
                continue
            if len(parts) < 6:
                continue
            try:
                row = [float(x) for x in parts]
                # Filter out ghost particles
                ghost_col = col_indices.get('is_ghost', len(row))
                if ghost_col < len(row) and row[ghost_col] == 0:
                    data.append(row)
                elif ghost_col >= len(row):
                    data.append(row)
            except ValueError:
                continue

    data = np.array(data)
    
    # Use column indices from header
    pos_x_col = col_indices.get('pos_x', 1)
    vel_x_col = col_indices.get('vel_x', 2)
    vel_t_col = col_indices.get('vel_t', None)  # New column for tangent velocity
    dens_col = col_indices.get('dens', 5)
    pres_col = col_indices.get('pres', 6)
    
    result = {
        'time': time,
        'pos_x': data[:, pos_x_col],
        'vel_x': data[:, vel_x_col],
        'vel_t': data[:, vel_t_col] if vel_t_col is not None and vel_t_col < data.shape[1] else np.zeros(len(data)),
        'dens': data[:, dens_col],
        'pres': data[:, pres_col],
    }
    
    # Compute Lorentz factor
    v_sq = result['vel_x']**2 + result['vel_t']**2
    result['gamma'] = 1.0 / np.sqrt(np.maximum(1.0 - v_sq, 1e-10))
    
    # Sort by position
    idx = np.argsort(result['pos_x'])
    for key in ['pos_x', 'vel_x', 'vel_t', 'dens', 'pres', 'gamma']:
        result[key] = result[key][idx]
    
    return result


def detect_tangent_velocity(snapshot_dir):
    """Detect tangent velocity from directory name"""
    dirname = os.path.basename(snapshot_dir)
    
    # Parse patterns like tangent_vt09_vt09, tangent_vt0_vt0, tangent_vt099_vt099
    match = re.search(r'vt(\d+)_vt(\d+)', dirname)
    if match:
        # Handle special cases: 0 -> 0.0, 09 -> 0.9, 099 -> 0.99, 9 -> 0.9
        def parse_vt(s):
            if s == '0':
                return 0.0
            elif s == '09' or s == '9':
                return 0.9
            elif s == '099' or s == '99':
                return 0.99
            else:
                # Generic: treat as decimal after 0.
                return float(f"0.{s}")
        
        vt_l = parse_vt(match.group(1))
        vt_r = parse_vt(match.group(2))
        return vt_l, vt_r
    
    # Default
    return 0.9, 0.9


def get_exact_tangent_solution(t, vt_l, vt_r, gamma=5.0/3.0, n_points=500):
    """
    Compute exact Riemann solution for tangent velocity test at time t.
    
    Uses the Pons et al. (2000) solver that correctly accounts for tangent
    velocity effects on both rarefactions and shocks.
    
    The tangent velocity test uses strong blast initial conditions:
      Left:  (P, rho, v^x, v^t) = (1000, 1.0, 0, vt_l)
      Right: (P, rho, v^x, v^t) = (0.01, 1.0, 0, vt_r)
    
    Key physics from Pons et al. (2000):
    - K = h * W * v^t is conserved across waves (h = enthalpy, W = Lorentz factor)
    - The characteristic speed depends on tangent velocity through Eq. (6)
    - This significantly affects the solution when vt is large
    
    Args:
        t: time
        vt_l, vt_r: tangential velocities left and right
        gamma: adiabatic index (default 5/3)
        n_points: number of points for solution
        
    Returns:
        dict with keys: x, pres, dens, vel_x, vel_t
    """
    # Strong blast initial conditions (Section 3.1.5)
    P_L, rho_L, vx_L = 1000.0, 1.0, 0.0
    P_R, rho_R, vx_R = 0.01, 1.0, 0.0
    
    if t <= 0:
        # Initial condition
        x = np.linspace(-0.5, 0.5, n_points)
        pres = np.where(x < 0, P_L, P_R)
        dens = np.where(x < 0, rho_L, rho_R)
        vel_x = np.where(x < 0, vx_L, vx_R)
        vel_t = np.where(x < 0, vt_l, vt_r)
        return {'x': x, 'pres': pres, 'dens': dens, 'vel_x': vel_x, 'vel_t': vel_t}
    
    # Use the Pons et al. (2000) solver with tangent velocity support
    solver = RelativisticRiemannSolverTangent(gamma)
    result = solver.solve(P_L, rho_L, vx_L, vt_l, P_R, rho_R, vx_R, vt_r, t, n_points)
    
    return {
        'x': result['x'],
        'pres': result['pres'],
        'dens': result['dens'],
        'vel_x': result['vel_x'],
        'vel_t': result['vel_t']
    }


def create_animation(snapshot_dir, output_file, skip=1, fps=10):
    """Create animation for tangent velocity shock tube with exact solution overlay"""
    
    # Detect tangent velocities
    vt_l, vt_r = detect_tangent_velocity(snapshot_dir)
    gamma_l = 1.0 / np.sqrt(1.0 - vt_l**2)
    gamma_r = 1.0 / np.sqrt(1.0 - vt_r**2)
    
    print(f"Tangent Velocity Test (Section 3.1.5)")
    print(f"  v^t_L = {vt_l}, γ_L = {gamma_l:.4f}")
    print(f"  v^t_R = {vt_r}, γ_R = {gamma_r:.4f}")
    
    # Find all snapshots
    pattern = os.path.join(snapshot_dir, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))
    
    if not files:
        print(f"No snapshot files found in {snapshot_dir}")
        return
    
    # Apply skip
    files = files[::skip]
    n_frames = len(files)
    
    print(f"Found {len(files)} snapshots (skip={skip})")
    print(f"Creating animation with {n_frames} frames...")
    
    # Load first snapshot to determine axis limits
    data0 = load_snapshot(files[0])
    
    # Strong blast initial conditions: P_L=1000, P_R=0.01
    dens_max = 8.0
    pres_max = 1100.0
    vel_max = 1.0
    gamma_max = max(gamma_l, gamma_r) * 1.5
    
    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Density plot
    ax1 = axes[0, 0]
    line_dens, = ax1.plot([], [], 'b.', markersize=2, alpha=0.6, label='SPH')
    line_dens_exact, = ax1.plot([], [], 'r-', linewidth=2, label='Exact')
    ax1.set_xlim(-0.55, 0.55)
    ax1.set_ylim(0, dens_max)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('Density (n)', fontsize=12)
    ax1.set_title('Rest-Frame Density', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Pressure plot
    ax2 = axes[0, 1]
    line_pres, = ax2.plot([], [], 'r.', markersize=2, alpha=0.6, label='SPH')
    line_pres_exact, = ax2.plot([], [], 'k-', linewidth=2, label='Exact')
    ax2.set_xlim(-0.55, 0.55)
    ax2.set_ylim(0, pres_max)
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('Pressure (P)', fontsize=12)
    ax2.set_title('Pressure', fontsize=14)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Normal velocity plot
    ax3 = axes[1, 0]
    line_vx, = ax3.plot([], [], 'g.', markersize=2, alpha=0.6, label='SPH v^x')
    line_vx_exact, = ax3.plot([], [], 'k-', linewidth=2, label='Exact v^x')
    ax3.set_xlim(-0.55, 0.55)
    ax3.set_ylim(-0.1, vel_max)
    ax3.set_xlabel('x', fontsize=12)
    ax3.set_ylabel('Velocity', fontsize=12)
    ax3.set_title('Normal Velocity v^x', fontsize=14)
    ax3.legend(loc='upper left')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Tangent velocity plot
    ax4 = axes[1, 1]
    line_vt, = ax4.plot([], [], 'm.', markersize=2, alpha=0.6, label='SPH v^t')
    line_vt_exact, = ax4.plot([], [], 'k-', linewidth=2, label='Exact v^t')
    ax4.axhline(y=vt_l, color='orange', linestyle=':', alpha=0.7, label=f'Initial v^t_L={vt_l}')
    if vt_l != vt_r:
        ax4.axhline(y=vt_r, color='cyan', linestyle=':', alpha=0.7, label=f'Initial v^t_R={vt_r}')
    ax4.set_xlim(-0.55, 0.55)
    ax4.set_ylim(-0.1, vel_max)
    ax4.set_xlabel('x', fontsize=12)
    ax4.set_ylabel('Tangent Velocity v^t', fontsize=12)
    ax4.set_title('Tangent Velocity v^t', fontsize=14)
    ax4.legend(loc='upper right', fontsize=9)
    ax4.grid(True, alpha=0.3)
    
    # Title
    title = fig.suptitle('', fontsize=14, fontweight='bold')
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    def init():
        line_dens.set_data([], [])
        line_dens_exact.set_data([], [])
        line_pres.set_data([], [])
        line_pres_exact.set_data([], [])
        line_vx.set_data([], [])
        line_vx_exact.set_data([], [])
        line_vt.set_data([], [])
        line_vt_exact.set_data([], [])
        title.set_text('')
        return (line_dens, line_dens_exact, line_pres, line_pres_exact, 
                line_vx, line_vx_exact, line_vt, line_vt_exact, title)
    
    def update(frame):
        data = load_snapshot(files[frame])
        t = data['time']
        
        # Get exact solution
        exact = get_exact_tangent_solution(t, vt_l, vt_r)
        
        # Update SPH plots
        line_dens.set_data(data['pos_x'], data['dens'])
        line_pres.set_data(data['pos_x'], data['pres'])
        line_vx.set_data(data['pos_x'], data['vel_x'])
        line_vt.set_data(data['pos_x'], data['vel_t'])
        
        # Update exact solution plots
        line_dens_exact.set_data(exact['x'], exact['dens'])
        line_pres_exact.set_data(exact['x'], exact['pres'])
        line_vx_exact.set_data(exact['x'], exact['vel_x'])
        line_vt_exact.set_data(exact['x'], exact['vel_t'])
        
        title.set_text(f'Tangent Velocity Test (v^t_L={vt_l}, v^t_R={vt_r}) - t = {t:.4f}')
        
        return (line_dens, line_dens_exact, line_pres, line_pres_exact,
                line_vx, line_vx_exact, line_vt, line_vt_exact, title)
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=n_frames, init_func=init, 
                         blit=True, interval=1000//fps)
    
    # Save animation
    print(f"Saving to {output_file}...")
    
    if output_file.endswith('.gif'):
        writer = PillowWriter(fps=fps)
        anim.save(output_file, writer=writer)
    elif output_file.endswith('.mp4'):
        writer = FFMpegWriter(fps=fps, bitrate=2000)
        anim.save(output_file, writer=writer)
    else:
        print("Warning: Unknown format, saving as GIF")
        writer = PillowWriter(fps=fps)
        anim.save(output_file + '.gif', writer=writer)
    
    print(f"Animation saved: {output_file}")
    plt.close()


def create_comparison_animation(results_dirs, output_file, skip=1, fps=10):
    """
    Create side-by-side comparison animation for multiple tangent velocity tests
    with exact solution overlay.
    
    Args:
        results_dirs: list of result directories to compare
        output_file: output animation file
        skip: skip every N frames
        fps: frames per second
    """
    n_cases = len(results_dirs)
    
    # Load metadata for each case
    cases = []
    for rdir in results_dirs:
        if not os.path.isdir(rdir):
            print(f"Warning: {rdir} not found, skipping")
            continue
        vt_l, vt_r = detect_tangent_velocity(rdir)
        files = sorted(glob.glob(os.path.join(rdir, "snapshot_*.csv")))
        if files:
            cases.append({
                'dir': rdir,
                'vt_l': vt_l,
                'vt_r': vt_r,
                'files': files[::skip],
                'label': f'v^t={vt_l}'
            })
    
    if not cases:
        print("No valid result directories found")
        return
    
    n_cases = len(cases)
    n_frames = min(len(c['files']) for c in cases)
    
    print(f"Comparing {n_cases} cases with {n_frames} frames each")
    
    # Create figure
    fig, axes = plt.subplots(2, n_cases, figsize=(5*n_cases, 8))
    if n_cases == 1:
        axes = axes.reshape(2, 1)
    
    lines_dens = []
    lines_dens_exact = []
    lines_pres = []
    lines_pres_exact = []
    
    for i, case in enumerate(cases):
        # Density row
        ax1 = axes[0, i]
        line_d, = ax1.plot([], [], 'b.', markersize=2, alpha=0.6, label='SPH')
        line_d_ex, = ax1.plot([], [], 'r-', linewidth=2, label='Exact')
        lines_dens.append(line_d)
        lines_dens_exact.append(line_d_ex)
        ax1.set_xlim(-0.55, 0.55)
        ax1.set_ylim(0, 8.0)
        ax1.set_xlabel('x')
        ax1.set_ylabel('Density (n)')
        ax1.set_title(f"v^t_L={case['vt_l']}, v^t_R={case['vt_r']}")
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(True, alpha=0.3)
        
        # Pressure row
        ax2 = axes[1, i]
        line_p, = ax2.plot([], [], 'r.', markersize=2, alpha=0.6, label='SPH')
        line_p_ex, = ax2.plot([], [], 'k-', linewidth=2, label='Exact')
        lines_pres.append(line_p)
        lines_pres_exact.append(line_p_ex)
        ax2.set_xlim(-0.55, 0.55)
        ax2.set_ylim(0, 1100.0)
        ax2.set_xlabel('x')
        ax2.set_ylabel('Pressure (P)')
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, alpha=0.3)
    
    title = fig.suptitle('', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    def init():
        for ld, ld_ex, lp, lp_ex in zip(lines_dens, lines_dens_exact, lines_pres, lines_pres_exact):
            ld.set_data([], [])
            ld_ex.set_data([], [])
            lp.set_data([], [])
            lp_ex.set_data([], [])
        title.set_text('')
        return lines_dens + lines_dens_exact + lines_pres + lines_pres_exact + [title]
    
    def update(frame):
        t = 0
        for i, case in enumerate(cases):
            data = load_snapshot(case['files'][frame])
            t = data['time']
            
            # Get exact solution for this case
            exact = get_exact_tangent_solution(t, case['vt_l'], case['vt_r'])
            
            # Update SPH data
            lines_dens[i].set_data(data['pos_x'], data['dens'])
            lines_pres[i].set_data(data['pos_x'], data['pres'])
            
            # Update exact solution
            lines_dens_exact[i].set_data(exact['x'], exact['dens'])
            lines_pres_exact[i].set_data(exact['x'], exact['pres'])
        
        title.set_text(f'Tangent Velocity Comparison (Section 3.1.5) - t = {t:.4f}')
        return lines_dens + lines_dens_exact + lines_pres + lines_pres_exact + [title]
    
    anim = FuncAnimation(fig, update, frames=n_frames, init_func=init,
                         blit=True, interval=1000//fps)
    
    print(f"Saving comparison to {output_file}...")
    if output_file.endswith('.gif'):
        writer = PillowWriter(fps=fps)
        anim.save(output_file, writer=writer)
    else:
        writer = FFMpegWriter(fps=fps, bitrate=2000)
        anim.save(output_file, writer=writer)
    
    print(f"Comparison animation saved: {output_file}")
    plt.close()


def main():
    if len(sys.argv) < 2:
        print("Usage: python animate_tangent_velocity.py <snapshot_directory> [output.gif] [skip] [fps]")
        print("       python animate_tangent_velocity.py --compare dir1 dir2 ... [output.gif]")
        print()
        print("Single test animation:")
        print("  python animate_tangent_velocity.py results/tangent_vt09_vt09 tangent.gif")
        print()
        print("Multi-test comparison:")
        print("  python animate_tangent_velocity.py --compare results/tangent_vt0_vt0 results/tangent_vt09_vt09 results/tangent_vt099_vt099 comparison.gif")
        print()
        print("Reference: Kitajima et al. (2025) Section 3.1.5")
        sys.exit(1)
    
    if sys.argv[1] == '--compare':
        # Comparison mode
        args = sys.argv[2:]
        dirs = [a for a in args if not a.endswith('.gif') and not a.endswith('.mp4')]
        output = [a for a in args if a.endswith('.gif') or a.endswith('.mp4')]
        output_file = output[0] if output else 'tangent_comparison.gif'
        create_comparison_animation(dirs, output_file)
    else:
        # Single animation mode
        snapshot_dir = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(snapshot_dir, 'tangent_animation.gif')
        skip = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        fps = int(sys.argv[4]) if len(sys.argv) > 4 else 10
        
        create_animation(snapshot_dir, output_file, skip, fps)


if __name__ == "__main__":
    main()
