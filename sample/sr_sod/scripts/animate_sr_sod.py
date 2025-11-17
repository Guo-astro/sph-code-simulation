#!/usr/bin/env python3
"""
Animate SR Sod shock tube test results
Generates animation showing density, pressure, and velocity evolution
Supports both EXACT and ITERATIVE Riemann solvers
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import csv
import sys
import json

def read_snapshot(filepath):
    """Read CSV snapshot file, handling comment lines"""
    data = []
    header = None
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('id,'):
                header = line.strip().split(',')
                continue
            data.append(line.strip().split(','))
    
    if not data or header is None:
        return None
    
    # Convert to structured array
    particles = []
    for row in data:
        p = {header[i]: float(row[i]) if i > 0 else int(row[i]) 
             for i in range(len(header))}
        particles.append(p)
    
    return particles

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

def main():
    # Accept results directory as command line argument
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path(__file__).parent.parent / 'results' / 'same_nu'
    
    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        print(f"Usage: {sys.argv[0]} [results_directory]")
        sys.exit(1)
    
    # Find all snapshots
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    print(f"Found {len(snapshots)} snapshots")
    
    if len(snapshots) == 0:
        print("No snapshots found!")
        sys.exit(1)
    
    # Read all snapshots
    all_data = []
    times = []
    
    for snap in snapshots:
        data = read_snapshot(snap)
        if data is None:
            continue
        all_data.append(data)
        
        # Extract time from filename or use index
        snap_num = int(snap.stem.split('_')[-1])
        times.append(snap_num * 0.01)  # Assuming dt=0.01
    
    print(f"Loaded {len(all_data)} snapshots, time range: {times[0]:.3f} to {times[-1]:.3f}")
    
    # Detect solver type
    solver_type, solver_name = detect_solver_type(results_dir)
    print(f"Detected Riemann solver: {solver_type} ({solver_name})")
    
    # Set up the figure
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    title = f'SR Sod Shock Tube (SR-GSPH)\n{solver_type} Riemann Solver ({solver_name})'
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    # Initialize plots
    line_dens, = axes[0].plot([], [], 'b.', markersize=1, alpha=0.6, label='SPH particles')
    line_pres, = axes[1].plot([], [], 'r.', markersize=1, alpha=0.6, label='SPH particles')
    line_vel, = axes[2].plot([], [], 'g.', markersize=1, alpha=0.6, label='SPH particles')
    
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12, fontweight='bold')
    
    # Set axis properties
    axes[0].set_ylabel('Density (n)', fontsize=11)
    axes[0].set_ylim(0, 1.2)
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=1.0, color='k', linestyle='--', alpha=0.3, linewidth=0.5, label='Initial left')
    axes[0].axhline(y=0.125, color='k', linestyle='--', alpha=0.3, linewidth=0.5, label='Initial right')
    axes[0].legend(loc='upper right', fontsize=9)
    
    axes[1].set_ylabel('Pressure (P)', fontsize=11)
    axes[1].set_ylim(0, 1.1)
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=1.0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)
    axes[1].axhline(y=0.1, color='k', linestyle='--', alpha=0.3, linewidth=0.5)
    
    axes[2].set_ylabel('Velocity (v)', fontsize=11)
    axes[2].set_xlabel('Position (x)', fontsize=11)
    axes[2].set_ylim(-0.1, 1.0)
    axes[2].grid(True, alpha=0.3)
    axes[2].axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    for ax in axes:
        ax.set_xlim(-0.5, 0.5)
        ax.axvline(x=0, color='k', linestyle=':', alpha=0.5, linewidth=1, label='Initial interface')
    
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    
    def animate(frame):
        """Animation update function"""
        data = all_data[frame]
        
        # Extract arrays
        x = np.array([p['pos_x'] for p in data])
        dens = np.array([p['dens'] for p in data])
        pres = np.array([p['pres'] for p in data])
        vel = np.array([p['vel_x'] for p in data])
        
        # Sort by position for cleaner plotting
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        dens = dens[sort_idx]
        pres = pres[sort_idx]
        vel = vel[sort_idx]
        
        # Update data
        line_dens.set_data(x, dens)
        line_pres.set_data(x, pres)
        line_vel.set_data(x, vel)
        
        time_text.set_text(f't = {times[frame]:.3f}')
        
        return line_dens, line_pres, line_vel, time_text
    
    # Create animation
    print("Creating animation...")
    anim = animation.FuncAnimation(
        fig, animate, frames=len(all_data),
        interval=100, blit=True, repeat=True
    )
    
    # Save animation
    output_file = results_dir / f'sr_sod_animation_{solver_type.lower()}.gif'
    print(f"Saving animation to {output_file}...")
    
    writer = animation.PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=100)
    
    print(f"✓ Animation saved: {output_file}")
    print(f"  Solver: {solver_type} ({solver_name})")
    print(f"  File size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"  Frames: {len(all_data)}")
    print(f"  Time range: t = {times[0]:.3f} to {times[-1]:.3f}")

if __name__ == '__main__':
    main()
