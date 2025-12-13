#!/usr/bin/env python3
"""
Animate SR shock tube with exact Riemann solution overlay

Supports all Kitajima et al. (2025) test cases:
- Sod problem (Section 3.1.1)
- Standard blast wave (Section 3.1.2)
- Strong blast wave (Section 3.1.3)
- Ultra-relativistic shock (Section 3.2)

Creates GIF/MP4 animation showing SPH simulation vs analytical solution
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import sys
import os
import glob

# Add docs directory to path for relativistic_riemann_solver
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'docs'))
from relativistic_riemann_solver import RelativisiticRiemannSolver

# Import shared test case definitions (SSOT)
from sr_test_cases import TEST_CASES, detect_test_type, get_initial_conditions


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
    
    # Use column indices from header, with fallbacks for old format
    pos_col = col_indices.get('pos_x', 1)
    vel_col = col_indices.get('vel_x', 2)
    dens_col = col_indices.get('dens', 5)
    pres_col = col_indices.get('pres', 6)
    
    result = {
        'time': time,
        'pos_x': data[:, pos_col],
        'vel_x': data[:, vel_col],
        'dens': data[:, dens_col],    # rest-frame density n
        'pres': data[:, pres_col],
    }
    
    # Sort by position
    idx = np.argsort(result['pos_x'])
    for key in ['pos_x', 'vel_x', 'dens', 'pres']:
        result[key] = result[key][idx]
    
    return result


def get_exact_solution(t, test_type='sod', v_left_override=None, gamma=5.0/3.0, n_points=500):
    """
    Compute exact Riemann solution for any Kitajima test case at time t.
    
    Args:
        t: time
        test_type: 'sod', 'blast_wave', 'strong_blast', or 'ultra_relat'
        v_left_override: override left velocity (for ultra_relat tests)
        gamma: adiabatic index (default 5/3)
        n_points: number of points for solution
        
    Returns:
        x, pres, dens, vel arrays
    """
    if test_type not in TEST_CASES:
        print(f"Warning: Unknown test type '{test_type}', using 'sod'")
        test_type = 'sod'
    
    test = TEST_CASES[test_type]
    pl, rhol, vell = test['left']
    pr, rhor, velr = test['right']
    
    # Override velocity for ultra-relativistic tests
    if v_left_override is not None and test_type == 'ultra_relat':
        vell = v_left_override
    
    if t <= 0:
        # Initial condition
        x = np.linspace(-0.5, 0.5, n_points)
        pres = np.where(x < 0, pl, pr)
        dens = np.where(x < 0, rhol, rhor)
        vel = np.where(x < 0, vell, velr)
        return x, pres, dens, vel
    
    solver = RelativisiticRiemannSolver(gamma)
    solver.set_initial_states(pl, rhol, vell, pr, rhor, velr)
    
    # Solve - returns solution on [0, 1] with discontinuity at 0.5
    x_raw, pres, dens, vel, _ = solver.solve(t, n=n_points)
    
    # Convert to our domain [-0.5, 0.5] with discontinuity at 0
    x = x_raw - 0.5
    
    return x, pres, dens, vel


def create_animation(snapshot_dir, output_file, skip=1, fps=10, test_type=None, v_left=None):
    """Create animation with exact solution overlay"""
    
    # Auto-detect test type if not specified
    if test_type is None:
        test_type, v_left_detected = detect_test_type(snapshot_dir)
        if v_left is None:
            v_left = v_left_detected
    
    print(f"Test type: {test_type}")
    if test_type in TEST_CASES:
        print(f"  {TEST_CASES[test_type]['name']}")
        print(f"  Left:  P={TEST_CASES[test_type]['left'][0]}, n={TEST_CASES[test_type]['left'][1]}, v={v_left if test_type == 'ultra_relat' else TEST_CASES[test_type]['left'][2]}")
        print(f"  Right: P={TEST_CASES[test_type]['right'][0]}, n={TEST_CASES[test_type]['right'][1]}, v={TEST_CASES[test_type]['right'][2]}")
    
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
    
    # Load first snapshot to set up figure and determine axis limits
    data0 = load_snapshot(files[0])
    
    # Determine axis limits based on test type
    if test_type == 'blast_wave':
        dens_max = 12.0
        pres_max = 16.0
        vel_max = 1.0
    elif test_type == 'strong_blast':
        dens_max = 8.0
        pres_max = 1100.0
        vel_max = 1.0
    elif test_type == 'ultra_relat':
        dens_max = 10.0
        pres_max = 50.0
        vel_max = 1.0
    else:  # sod
        dens_max = 1.2
        pres_max = 1.2
        vel_max = 1.0
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Density plot
    ax1 = axes[0]
    line_sph_dens, = ax1.plot([], [], 'b.', markersize=2, alpha=0.6, label='SPH')
    line_exact_dens, = ax1.plot([], [], 'r-', linewidth=2, label='Exact')
    ax1.set_xlim(-0.55, 0.55)
    ax1.set_ylim(0, dens_max)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('Density (n)', fontsize=12)
    ax1.set_title('Density', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Pressure plot
    ax2 = axes[1]
    line_sph_pres, = ax2.plot([], [], 'b.', markersize=2, alpha=0.6, label='SPH')
    line_exact_pres, = ax2.plot([], [], 'r-', linewidth=2, label='Exact')
    ax2.set_xlim(-0.55, 0.55)
    ax2.set_ylim(0, pres_max)
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('Pressure (P)', fontsize=12)
    ax2.set_title('Pressure', fontsize=14)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Velocity plot
    ax3 = axes[2]
    line_sph_vel, = ax3.plot([], [], 'b.', markersize=2, alpha=0.6, label='SPH')
    line_exact_vel, = ax3.plot([], [], 'r-', linewidth=2, label='Exact')
    ax3.set_xlim(-0.55, 0.55)
    ax3.set_ylim(-0.1, vel_max)
    ax3.set_xlabel('x', fontsize=12)
    ax3.set_ylabel('Velocity (v)', fontsize=12)
    ax3.set_title('Velocity', fontsize=14)
    ax3.legend(loc='upper left')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Title with time
    test_title = TEST_CASES.get(test_type, {}).get('name', 'Unknown Test')
    title = fig.suptitle('', fontsize=14, fontweight='bold')
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    def init():
        line_sph_dens.set_data([], [])
        line_exact_dens.set_data([], [])
        line_sph_pres.set_data([], [])
        line_exact_pres.set_data([], [])
        line_sph_vel.set_data([], [])
        line_exact_vel.set_data([], [])
        title.set_text('')
        return line_sph_dens, line_exact_dens, line_sph_pres, line_exact_pres, line_sph_vel, line_exact_vel, title
    
    def update(frame):
        data = load_snapshot(files[frame])
        t = data['time']
        
        # Get exact solution for this test type
        exact_x, exact_p, exact_rho, exact_v = get_exact_solution(t, test_type, v_left)
        
        # Update SPH data
        line_sph_dens.set_data(data['pos_x'], data['dens'])
        line_sph_pres.set_data(data['pos_x'], data['pres'])
        line_sph_vel.set_data(data['pos_x'], data['vel_x'])
        
        # Update exact solution
        line_exact_dens.set_data(exact_x, exact_rho)
        line_exact_pres.set_data(exact_x, exact_p)
        line_exact_vel.set_data(exact_x, exact_v)
        
        title.set_text(f'{test_title} - SRGSPH vs Exact (t = {t:.4f})')
        
        return line_sph_dens, line_exact_dens, line_sph_pres, line_exact_pres, line_sph_vel, line_exact_vel, title
    
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


def main():
    if len(sys.argv) < 2:
        print("Usage: python animate_sr_sod_with_exact.py <snapshot_directory> [output.gif] [skip] [fps] [test_type] [v_left]")
        print()
        print("Example:")
        print("  python animate_sr_sod_with_exact.py results/sod_fig1 animation.gif 2 15")
        print("  python animate_sr_sod_with_exact.py results/blast_wave blast.gif 1 10 blast_wave")
        print("  python animate_sr_sod_with_exact.py results/ultra_v099 ultra.gif 1 10 ultra_relat 0.99")
        print()
        print("Test types (from Kitajima et al. 2025):")
        print("  sod          - Sod shock tube (Section 3.1.1)")
        print("  blast_wave   - Standard relativistic blast wave (Section 3.1.2)")
        print("  strong_blast - Strong relativistic blast wave (Section 3.1.3)")
        print("  ultra_relat  - Ultra-relativistic shock (Section 3.2)")
        print()
        print("If test_type is not specified, it will be auto-detected from directory name.")
        sys.exit(1)
    
    snapshot_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(snapshot_dir, 'animation_with_exact.gif')
    skip = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    fps = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    test_type = sys.argv[5] if len(sys.argv) > 5 else None
    v_left = float(sys.argv[6]) if len(sys.argv) > 6 else None
    
    create_animation(snapshot_dir, output_file, skip, fps, test_type, v_left)


if __name__ == "__main__":
    main()
