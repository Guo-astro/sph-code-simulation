#!/usr/bin/env python3
"""
Animate SR-SOD simulation results.

Creates animations showing the evolution of density, velocity, and pressure profiles
for Special Relativistic Godunov SPH shock tube simulations.

Usage:
    python scripts/animate_sr_sod_results.py [results_dir] [output_file]
    
Examples:
    python scripts/animate_sr_sod_results.py
    python scripts/animate_sr_sod_results.py sample/sr_sod/results/sod_fig1 sr_sod_animation.gif
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec


def load_snapshot(filepath):
    """Load a CSV snapshot file."""
    # Find the header line (starts with 'id,')
    header_line = 0
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('id,'):
                header_line = i
                break
    
    # Read with numpy, skipping metadata lines
    data = np.genfromtxt(filepath, delimiter=',', names=True, 
                         skip_header=header_line, dtype=None, encoding='utf-8')
    return data


def get_snapshot_files(results_dir):
    """Get sorted list of snapshot files."""
    pattern = os.path.join(results_dir, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))
    return files


def extract_time_from_snapshot(filepath):
    """Extract simulation time from the first data row (if available) or from filename."""
    try:
        data = load_snapshot(filepath)
        # The time is usually not in the CSV, so we use the file index
        idx = int(os.path.basename(filepath).split('_')[1].split('.')[0])
        return idx
    except:
        return 0


def create_animation(results_dir, output_file=None, fps=10, dpi=150):
    """
    Create an animation of the SR-SOD results.
    
    Parameters:
    -----------
    results_dir : str
        Directory containing snapshot CSV files
    output_file : str, optional
        Output filename (e.g., 'animation.gif' or 'animation.mp4')
    fps : int
        Frames per second
    dpi : int
        Resolution for output
    """
    # Get snapshot files
    snapshot_files = get_snapshot_files(results_dir)
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return None
    
    print(f"Found {len(snapshot_files)} snapshots in {results_dir}")
    
    # Load first snapshot to set up axes
    data0 = load_snapshot(snapshot_files[0])
    
    # Check available columns
    columns = data0.dtype.names
    print(f"Available columns: {columns}")
    
    # Determine which columns to use
    x_col = 'pos_x' if 'pos_x' in columns else 'x' if 'x' in columns else columns[0]
    
    # Find density, velocity, pressure columns
    dens_col = None
    vel_col = None
    pres_col = None
    
    for col in columns:
        if 'dens' in col.lower() or col == 'rho':
            dens_col = col
        if 'vel' in col.lower() and '_x' in col.lower():
            vel_col = col
        elif 'vel' in col.lower() and vel_col is None:
            vel_col = col
        if 'pres' in col.lower() or col == 'P' or col == 'p':
            pres_col = col
    
    print(f"Using columns: x={x_col}, density={dens_col}, velocity={vel_col}, pressure={pres_col}")
    
    # Set up figure with three subplots
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(3, 1, figure=fig, hspace=0.3)
    
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    
    # Get data ranges from all snapshots for consistent axes
    all_x, all_dens, all_vel, all_pres = [], [], [], []
    
    for f in snapshot_files[::max(1, len(snapshot_files)//10)]:  # Sample 10 files
        data = load_snapshot(f)
        # Filter out ghost particles if is_ghost column exists
        if 'is_ghost' in data.dtype.names:
            mask = data['is_ghost'] == 0
            data = data[mask]
        
        all_x.extend(data[x_col])
        if dens_col:
            all_dens.extend(data[dens_col])
        if vel_col:
            all_vel.extend(data[vel_col])
        if pres_col:
            all_pres.extend(data[pres_col])
    
    x_min, x_max = min(all_x), max(all_x)
    x_margin = (x_max - x_min) * 0.05
    
    # Initialize plots
    line1, = ax1.plot([], [], 'b.', markersize=2, alpha=0.7)
    line2, = ax2.plot([], [], 'r.', markersize=2, alpha=0.7)
    line3, = ax3.plot([], [], 'g.', markersize=2, alpha=0.7)
    
    # Set axis labels and limits
    ax1.set_xlim(x_min - x_margin, x_max + x_margin)
    ax2.set_xlim(x_min - x_margin, x_max + x_margin)
    ax3.set_xlim(x_min - x_margin, x_max + x_margin)
    
    if all_dens:
        dens_min, dens_max = min(all_dens), max(all_dens)
        ax1.set_ylim(dens_min * 0.9, dens_max * 1.1)
    if all_vel:
        vel_min, vel_max = min(all_vel), max(all_vel)
        vel_margin = max(0.1, (vel_max - vel_min) * 0.1)
        ax2.set_ylim(vel_min - vel_margin, vel_max + vel_margin)
    if all_pres:
        pres_min, pres_max = min(all_pres), max(all_pres)
        ax3.set_ylim(pres_min * 0.9, pres_max * 1.1)
    
    ax1.set_ylabel('Density', fontsize=12)
    ax2.set_ylabel('Velocity', fontsize=12)
    ax3.set_ylabel('Pressure', fontsize=12)
    ax3.set_xlabel('Position x', fontsize=12)
    
    ax1.grid(True, alpha=0.3)
    ax2.grid(True, alpha=0.3)
    ax3.grid(True, alpha=0.3)
    
    # Title
    title = fig.suptitle('', fontsize=14)
    
    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        title.set_text('')
        return line1, line2, line3, title
    
    def animate(frame):
        data = load_snapshot(snapshot_files[frame])
        
        # Filter out ghost particles if is_ghost column exists
        if 'is_ghost' in data.dtype.names:
            mask = data['is_ghost'] == 0
            data = data[mask]
        
        x = data[x_col]
        
        # Sort by x for cleaner visualization
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        
        if dens_col:
            dens = data[dens_col][sort_idx]
            line1.set_data(x, dens)
        
        if vel_col:
            vel = data[vel_col][sort_idx]
            line2.set_data(x, vel)
        
        if pres_col:
            pres = data[pres_col][sort_idx]
            line3.set_data(x, pres)
        
        title.set_text(f'SR-SOD Shock Tube - Frame {frame+1}/{len(snapshot_files)}')
        
        return line1, line2, line3, title
    
    # Create animation
    anim = animation.FuncAnimation(
        fig, animate, init_func=init,
        frames=len(snapshot_files), interval=1000//fps,
        blit=True
    )
    
    # Save or show
    if output_file:
        print(f"Saving animation to {output_file}...")
        if output_file.endswith('.gif'):
            writer = animation.PillowWriter(fps=fps)
        elif output_file.endswith('.mp4'):
            writer = animation.FFMpegWriter(fps=fps, bitrate=2000)
        else:
            writer = animation.PillowWriter(fps=fps)
        
        anim.save(output_file, writer=writer, dpi=dpi)
        print(f"Animation saved to {output_file}")
    else:
        plt.show()
    
    plt.close()
    return anim


def load_energy_data(results_dir):
    """Load energy data from energy.dat file."""
    energy_file = os.path.join(results_dir, "energy.dat")
    if not os.path.exists(energy_file):
        return None
    
    try:
        data = np.loadtxt(energy_file, comments='#')
        return {
            'time': data[:, 0],
            'kinetic': data[:, 1],
            'thermal': data[:, 2],
            'potential': data[:, 3],
            'total': data[:, 4]
        }
    except:
        return None


def plot_energy_conservation(results_dir, output_file=None):
    """Plot energy conservation over time."""
    energy = load_energy_data(results_dir)
    if energy is None:
        print(f"No energy data found in {results_dir}")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Energy components
    ax1.plot(energy['time'], energy['kinetic'], 'r-', label='Kinetic', linewidth=1.5)
    ax1.plot(energy['time'], energy['thermal'], 'b-', label='Thermal', linewidth=1.5)
    ax1.plot(energy['time'], energy['total'], 'k-', label='Total', linewidth=2)
    ax1.set_ylabel('Energy', fontsize=12)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Energy Evolution', fontsize=14)
    
    # Conservation (relative error)
    E0 = energy['total'][0]
    relative_error = (energy['total'] - E0) / E0 * 100
    ax2.plot(energy['time'], relative_error, 'k-', linewidth=1.5)
    ax2.set_xlabel('Time', fontsize=12)
    ax2.set_ylabel('Energy Error (%)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    
    # Add final conservation info
    final_error = relative_error[-1]
    ax2.text(0.95, 0.95, f'Final error: {final_error:.2f}%',
             transform=ax2.transAxes, ha='right', va='top',
             fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Energy plot saved to {output_file}")
    else:
        plt.show()
    
    plt.close()


if __name__ == "__main__":
    # Default parameters
    default_results_dir = "sample/sr_sod/results/sod_fig1"
    default_output = "sr_sod_animation.gif"
    
    # Parse arguments
    if len(sys.argv) >= 2:
        results_dir = sys.argv[1]
    else:
        results_dir = default_results_dir
    
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        output_file = default_output
    
    # Create animation
    print(f"Creating animation from {results_dir}...")
    create_animation(results_dir, output_file, fps=10, dpi=150)
    
    # Also create energy conservation plot
    energy_output = output_file.replace('.gif', '_energy.png').replace('.mp4', '_energy.png')
    plot_energy_conservation(results_dir, energy_output)
    
    print("Done!")
