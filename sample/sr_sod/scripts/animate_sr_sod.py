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
        results_dir = Path(__file__).parent.parent / 'results' / 'sharp'
    
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
    
    # Set up the figure with dark background for better contrast
    fig, axes = plt.subplots(3, 1, figsize=(14, 11))
    fig.patch.set_facecolor('#f8f9fa')
    title = f'SR Sod Shock Tube - Special Relativistic Godunov SPH\n{solver_type} Riemann Solver ({solver_name})'
    fig.suptitle(title, fontsize=16, fontweight='bold', color='#2c3e50')

    # Use vibrant colors with good contrast
    # Modern color palette inspired by scientific visualization
    color_dens = '#e74c3c'  # Vibrant red
    color_pres = '#3498db'  # Bright blue
    color_vel = '#2ecc71'   # Fresh green

    # Initialize plots with BOTH scatter and line plots for clarity
    scatter_dens = axes[0].scatter([], [], s=40, c=color_dens, alpha=0.7,
                                   edgecolors='white', linewidths=0.5, label='SPH particles')
    line_dens, = axes[0].plot([], [], color=color_dens, linewidth=2, alpha=0.8, label='Trend')

    scatter_pres = axes[1].scatter([], [], s=40, c=color_pres, alpha=0.7,
                                   edgecolors='white', linewidths=0.5, label='SPH particles')
    line_pres, = axes[1].plot([], [], color=color_pres, linewidth=2, alpha=0.8, label='Trend')

    scatter_vel = axes[2].scatter([], [], s=40, c=color_vel, alpha=0.7,
                                  edgecolors='white', linewidths=0.5, label='SPH particles')
    line_vel, = axes[2].plot([], [], color=color_vel, linewidth=2, alpha=0.8, label='Trend')

    time_text = fig.text(0.5, 0.96, '', ha='center', fontsize=14, fontweight='bold',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='#34495e', linewidth=2))

    # Set axis properties with better styling
    for i, (ax, ylabel, color) in enumerate(zip(axes,
                                                  ['Rest-Frame Density (n)', 'Pressure (P)', 'Velocity (v/c)'],
                                                  [color_dens, color_pres, color_vel])):
        ax.set_facecolor('#ffffff')
        ax.set_ylabel(ylabel, fontsize=13, fontweight='bold', color=color)
        ax.tick_params(labelsize=11)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8, color='gray')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_color(color)

        # Add subtle background shading for initial regions
        ax.axvspan(-0.5, 0, alpha=0.05, color='blue', label='Initial left' if i == 0 else '')
        ax.axvspan(0, 0.5, alpha=0.05, color='red', label='Initial right' if i == 0 else '')

    # Density panel
    axes[0].set_ylim(0, 1.3)
    axes[0].axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, label='Left IC')
    axes[0].axhline(y=0.125, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, label='Right IC')
    axes[0].legend(loc='upper right', fontsize=10, framealpha=0.9)

    # Pressure panel
    axes[1].set_ylim(0, 1.15)
    axes[1].axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
    axes[1].axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)

    # Velocity panel
    axes[2].set_ylabel('Velocity (v/c)', fontsize=13, fontweight='bold', color=color_vel)
    axes[2].set_xlabel('Position (x)', fontsize=13, fontweight='bold')
    axes[2].set_ylim(-0.5, 0.6)
    axes[2].axhline(y=0, color='black', linestyle='-', alpha=0.6, linewidth=1.5)

    for ax in axes:
        ax.set_xlim(-0.5, 0.5)
        ax.axvline(x=0, color='black', linestyle=':', alpha=0.6, linewidth=2, label='Initial interface')
    
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
        x_sorted = x[sort_idx]
        dens_sorted = dens[sort_idx]
        pres_sorted = pres[sort_idx]
        vel_sorted = vel[sort_idx]

        # Update scatter plots (individual particles)
        scatter_dens.set_offsets(np.c_[x, dens])
        scatter_pres.set_offsets(np.c_[x, pres])
        scatter_vel.set_offsets(np.c_[x, vel])

        # Update line plots (trend lines connecting sorted particles)
        line_dens.set_data(x_sorted, dens_sorted)
        line_pres.set_data(x_sorted, pres_sorted)
        line_vel.set_data(x_sorted, vel_sorted)

        # Update time display with more info
        time_text.set_text(f't = {times[frame]:.4f}  |  Frame {frame+1}/{len(all_data)}')

        return scatter_dens, line_dens, scatter_pres, line_pres, scatter_vel, line_vel, time_text
    
    # Create animation with smoother playback
    print("Creating animation...")
    anim = animation.FuncAnimation(
        fig, animate, frames=len(all_data),
        interval=80, blit=True, repeat=True  # 80ms = ~12.5 fps for smoother viewing
    )

    # Save animation - try MP4 first (faster), fallback to GIF
    try:
        output_file = results_dir / f'sr_sod_animation_{solver_type.lower()}.mp4'
        print(f"Saving animation to {output_file} (MP4)...")
        # Higher quality: increased bitrate and DPI
        writer = animation.FFMpegWriter(fps=12, bitrate=3000, codec='libx264',
                                       extra_args=['-pix_fmt', 'yuv420p'])
        anim.save(output_file, writer=writer, dpi=100)
        print(f"✓ Animation saved: {output_file}")
    except (RuntimeError, FileNotFoundError) as e:
        # FFmpeg not available, use GIF
        output_file = results_dir / f'sr_sod_animation_{solver_type.lower()}.gif'
        print(f"FFmpeg not available ({e}), falling back to GIF...")
        print(f"Saving animation to {output_file} (GIF format)...")
        writer = animation.PillowWriter(fps=10)
        anim.save(output_file, writer=writer, dpi=80)
        print(f"✓ Animation saved: {output_file}")
    print(f"  Solver: {solver_type} ({solver_name})")
    print(f"  File size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"  Frames: {len(all_data)}")
    print(f"  Time range: t = {times[0]:.3f} to {times[-1]:.3f}")

if __name__ == '__main__':
    main()
