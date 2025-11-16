#!/usr/bin/env python3
"""
Visualization script for SR Sod shock tube results
Creates animations from CSV output files
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import glob
import sys

def load_snapshot(filename):
    """Load a single CSV snapshot"""
    try:
        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        if data.size == 0:
            return None
        
        # Columns: x, y, z, vx, vy, vz, rho, P, u, m, h, gamma_lor, S_mag
        return {
            'x': data[:, 0],
            'rho': data[:, 6],
            'P': data[:, 7],
            'u': data[:, 8],
            'v': np.sqrt(data[:, 3]**2 + data[:, 4]**2 + data[:, 5]**2),
            'gamma': data[:, 11] if data.shape[1] > 11 else np.ones(len(data)),
            'S': data[:, 12] if data.shape[1] > 12 else np.zeros(len(data))
        }
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def create_animation(results_dir):
    """Create animation from CSV files"""
    
    # Find all output files
    csv_files = sorted(glob.glob(f"{results_dir}/output_*.csv"))
    
    if not csv_files:
        print(f"No CSV files found in {results_dir}")
        return
    
    print(f"Found {len(csv_files)} snapshots")
    
    # Load all snapshots
    snapshots = []
    times = []
    for f in csv_files:
        data = load_snapshot(f)
        if data is not None:
            snapshots.append(data)
            # Extract time from filename: output_0.010000.csv
            try:
                time_str = Path(f).stem.split('_')[1]
                times.append(float(time_str))
            except:
                times.append(len(snapshots) * 0.01)
    
    if not snapshots:
        print("No valid snapshots loaded")
        return
    
    print(f"Loaded {len(snapshots)} valid snapshots")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('SR Sod Shock Tube', fontsize=16)
    
    # Subplot titles
    titles = ['Density', 'Pressure', 'Velocity', 'Lorentz Factor']
    ylabels = [r'$\rho$', r'$P$', r'$v/c$', r'$\gamma$']
    variables = ['rho', 'P', 'v', 'gamma']
    
    lines = []
    scatter_plots = []
    
    for ax, title, ylabel in zip(axes.flat, titles, ylabels):
        ax.set_xlabel('x')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        line, = ax.plot([], [], 'b.', markersize=2, alpha=0.6)
        lines.append(line)
    
    # Set axis limits based on all data
    x_min = min([np.min(s['x']) for s in snapshots])
    x_max = max([np.max(s['x']) for s in snapshots])
    
    for ax, var in zip(axes.flat, variables):
        ax.set_xlim(x_min, x_max)
        
        if var == 'rho':
            ax.set_ylim(0, 1.2)
        elif var == 'P':
            ax.set_ylim(0, 1.2)
        elif var == 'v':
            ax.set_ylim(-0.1, 1.0)
        elif var == 'gamma':
            ax.set_ylim(0.9, 2.0)
    
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12)
    
    def animate(frame):
        if frame >= len(snapshots):
            return lines + [time_text]
        
        snapshot = snapshots[frame]
        time = times[frame]
        
        for line, var in zip(lines, variables):
            # Sort by x for better visualization
            idx = np.argsort(snapshot['x'])
            line.set_data(snapshot['x'][idx], snapshot[var][idx])
        
        time_text.set_text(f't = {time:.4f}')
        
        return lines + [time_text]
    
    # Create animation
    anim = animation.FuncAnimation(
        fig, animate, frames=len(snapshots),
        interval=100, blit=True, repeat=True
    )
    
    # Save animation
    output_file = f"{results_dir}/sr_sod_animation.gif"
    print(f"Saving animation to {output_file}...")
    
    try:
        writer = animation.PillowWriter(fps=10)
        anim.save(output_file, writer=writer)
        print(f"✓ Animation saved: {output_file}")
    except Exception as e:
        print(f"Error saving animation: {e}")
        print("Displaying plot instead...")
        plt.show()
    
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = "sample/sr_sod/results"
    
    create_animation(results_dir)
