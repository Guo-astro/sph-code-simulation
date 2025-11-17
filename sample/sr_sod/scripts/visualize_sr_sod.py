#!/usr/bin/env python3
"""
Visualize SR Sod shock tube results from CSV output.
Creates animated plots of density, velocity, and pressure evolution.
Supports both EXACT and ITERATIVE Riemann solvers.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import glob
import sys
import json

def read_snapshot(filename):
    """Read CSV snapshot file."""
    # Skip header lines (metadata comments) and read data
    # CSV format: id,pos_x,pos_y,pos_z,vel_x,...,dens,pres,...
    data = np.loadtxt(filename, delimiter=',', skiprows=51)
    x = data[:, 1]  # pos_x
    vx = data[:, 4]  # vel_x
    rho = data[:, 11]  # dens
    P = data[:, 12]  # pres
    return x, rho, vx, P

def detect_solver_type(results_dir):
    """Detect which Riemann solver was used from directory name or config"""
    dir_name = Path(results_dir).name.lower()
    
    if 'iterative' in dir_name:
        return 'ITERATIVE', 'van Leer Relativistic'
    elif 'exact' in dir_name:
        return 'EXACT', 'Kitajima'
    else:
        # Check for config file or metadata
        parent = Path(results_dir).parent.parent
        config_files = list(parent.glob('*.json'))
        for cfg in config_files:
            try:
                with open(cfg) as f:
                    data = json.load(f)
                    solver = data.get('riemannSolverSRGSPH', 'EXACT')
                    if solver == 'ITERATIVE':
                        return 'ITERATIVE', 'van Leer Relativistic'
                    else:
                        return 'EXACT', 'Kitajima'
            except:
                pass
        return 'EXACT', 'Kitajima'  # Default

def create_animation(results_dir, output_file='sr_sod_evolution.gif'):
    """Create animation from snapshot files."""
    results_path = Path(results_dir)
    snapshot_files = sorted(results_path.glob('snapshot_*.csv'))
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Detect solver type
    solver_type, solver_name = detect_solver_type(results_dir)
    
    # Read all snapshots
    snapshots = []
    times = []
    for i, fname in enumerate(snapshot_files):
        x, rho, vx, P = read_snapshot(fname)
        snapshots.append((x, rho, vx, P))
        # Extract time from filename or index
        times.append(i * 0.01)  # 0.01 time interval
    
    # Create figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    title = f'SR Sod Shock Tube Evolution (SR-GSPH)\n{solver_type} Riemann Solver ({solver_name})'
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    # Initialize plots
    line1, = ax1.plot([], [], 'b-', lw=2, label='Density')
    line2, = ax2.plot([], [], 'r-', lw=2, label='Velocity')
    line3, = ax3.plot([], [], 'g-', lw=2, label='Pressure')
    
    # Set labels and grids
    ax1.set_ylabel('Density ρ', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_ylabel('Velocity v/c', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    ax3.set_xlabel('Position x', fontsize=12)
    ax3.set_ylabel('Pressure P', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Determine axis limits from all data
    all_x = np.concatenate([s[0] for s in snapshots])
    all_rho = np.concatenate([s[1] for s in snapshots])
    all_vx = np.concatenate([s[2] for s in snapshots])
    all_P = np.concatenate([s[3] for s in snapshots])
    
    x_min, x_max = all_x.min(), all_x.max()
    x_margin = (x_max - x_min) * 0.05
    
    ax1.set_xlim(x_min - x_margin, x_max + x_margin)
    ax2.set_xlim(x_min - x_margin, x_max + x_margin)
    ax3.set_xlim(x_min - x_margin, x_max + x_margin)
    
    ax1.set_ylim(all_rho.min() * 0.95, all_rho.max() * 1.05)
    ax2.set_ylim(all_vx.min() * 1.1, all_vx.max() * 1.1)
    ax3.set_ylim(all_P.min() * 0.95, all_P.max() * 1.05)
    
    # Time text
    time_text = fig.text(0.5, 0.96, '', ha='center', fontsize=12, 
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        time_text.set_text('')
        return line1, line2, line3, time_text
    
    def animate(i):
        x, rho, vx, P = snapshots[i]
        line1.set_data(x, rho)
        line2.set_data(x, vx)
        line3.set_data(x, P)
        time_text.set_text(f't = {times[i]:.3f}')
        return line1, line2, line3, time_text
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=len(snapshots), interval=100,
                                  blit=True, repeat=True)
    
    # Save animation
    output_path = results_path / output_file
    print(f"Saving animation to {output_path}")
    anim.save(str(output_path), writer='pillow', fps=10, dpi=100)
    print(f"✓ Animation saved: {output_path}")
    
    plt.close(fig)

def plot_final_state(results_dir, output_file='sr_sod_final.png'):
    """Plot final state."""
    results_path = Path(results_dir)
    snapshot_files = sorted(results_path.glob('snapshot_*.csv'))
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return
    
    # Detect solver type
    solver_type, solver_name = detect_solver_type(results_dir)
    
    # Read final snapshot
    x, rho, vx, P = read_snapshot(snapshot_files[-1])
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    title = f'SR Sod Shock Tube - Final State (SR-GSPH)\n{solver_type} Riemann Solver ({solver_name})'
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    ax1.plot(x, rho, 'b-', lw=2)
    ax1.set_ylabel('Density ρ', fontsize=12)
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(x, vx, 'r-', lw=2)
    ax2.set_ylabel('Velocity v/c', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    ax3.plot(x, P, 'g-', lw=2)
    ax3.set_xlabel('Position x', fontsize=12)
    ax3.set_ylabel('Pressure P', fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = results_path / output_file
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Final state plot saved: {output_path}")
    plt.close(fig)

def plot_evolution(results_dir, output_file='sr_sod_evolution.png'):
    """Plot evolution at multiple times."""
    results_path = Path(results_dir)
    snapshot_files = sorted(results_path.glob('snapshot_*.csv'))
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return
    
    # Detect solver type
    solver_type, solver_name = detect_solver_type(results_dir)
    
    # Select time slices
    n_times = min(6, len(snapshot_files))
    indices = np.linspace(0, len(snapshot_files)-1, n_times, dtype=int)
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 12))
    title = f'SR Sod Shock Tube Evolution (SR-GSPH)\n{solver_type} Riemann Solver ({solver_name})'
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    colors = plt.cm.viridis(np.linspace(0, 1, n_times))
    
    for idx, color in zip(indices, colors):
        x, rho, vx, P = read_snapshot(snapshot_files[idx])
        t = idx * 0.01
        
        axes[0].plot(x, rho, color=color, lw=1.5, label=f't={t:.2f}')
        axes[1].plot(x, vx, color=color, lw=1.5, label=f't={t:.2f}')
        axes[2].plot(x, P, color=color, lw=1.5, label=f't={t:.2f}')
    
    axes[0].set_ylabel('Density ρ', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc='best', fontsize=9)
    
    axes[1].set_ylabel('Velocity v/c', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc='best', fontsize=9)
    
    axes[2].set_xlabel('Position x', fontsize=12)
    axes[2].set_ylabel('Pressure P', fontsize=12)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc='best', fontsize=9)
    
    plt.tight_layout()
    output_path = results_path / output_file
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Evolution plot saved: {output_path}")
    plt.close(fig)

if __name__ == '__main__':
    # Default results directory
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = 'sample/same_nu/results'
    
    print(f"Visualizing results from: {results_dir}")
    
    # Create all visualizations
    plot_evolution(results_dir)
    plot_final_state(results_dir)
    create_animation(results_dir)
    
    print("\n✓ All visualizations complete!")
