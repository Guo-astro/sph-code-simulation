#!/usr/bin/env python3
"""
Visualization script for SR-GSPH test results
Plots density, velocity, pressure, and Lorentz factor profiles
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from pathlib import Path

def read_csv_snapshot(filename):
    """Read SPH snapshot from CSV file"""
    try:
        # Skip header lines starting with #
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                if not line.startswith('#'):
                    header_end = i
                    break
        
        # Read data
        data = np.genfromtxt(filename, delimiter=',', names=True, skip_header=header_end)
        return data
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

def plot_sr_profiles(results_dir, output_prefix='sr_profiles'):
    """Plot SR-GSPH simulation profiles"""
    
    # Find all snapshot files
    snapshot_files = sorted(glob.glob(os.path.join(results_dir, 'snapshot_*.csv')))
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Read initial and final snapshots
    initial = read_csv_snapshot(snapshot_files[0])
    final = read_csv_snapshot(snapshot_files[-1])
    
    if initial is None or final is None:
        print("Failed to read snapshots")
        return
    
    # Create figure with subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('SR-GSPH Simulation Results', fontsize=16)
    
    # Plot 1: Density
    ax = axs[0, 0]
    ax.scatter(initial['pos_x'], initial['dens'], s=2, alpha=0.5, label='Initial')
    if len(snapshot_files) > 1:
        ax.scatter(final['pos_x'], final['dens'], s=2, alpha=0.5, label='Final')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Density ρ')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Velocity
    ax = axs[0, 1]
    ax.scatter(initial['pos_x'], initial['vel_x'], s=2, alpha=0.5, label='Initial')
    if len(snapshot_files) > 1:
        ax.scatter(final['pos_x'], final['vel_x'], s=2, alpha=0.5, label='Final')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Velocity v')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Pressure
    ax = axs[1, 0]
    ax.scatter(initial['pos_x'], initial['pres'], s=2, alpha=0.5, label='Initial')
    if len(snapshot_files) > 1:
        ax.scatter(final['pos_x'], final['pres'], s=2, alpha=0.5, label='Final')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Pressure P')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Energy
    ax = axs[1, 1]
    ax.scatter(initial['pos_x'], initial['ene'], s=2, alpha=0.5, label='Initial')
    if len(snapshot_files) > 1:
        ax.scatter(final['pos_x'], final['ene'], s=2, alpha=0.5, label='Final')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Internal Energy u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(results_dir, f'{output_prefix}_profiles.png')
    plt.savefig(output_file, dpi=150)
    print(f"✓ Saved: {output_file}")
    plt.close()
    
    # Create animation frames if multiple snapshots
    if len(snapshot_files) > 2:
        create_animation(snapshot_files, results_dir, output_prefix)
    else:
        print(f"Need at least 3 snapshots for animation (found {len(snapshot_files)})")

def create_animation(snapshot_files, results_dir, output_prefix):
    """Create animation from snapshots"""
    try:
        import matplotlib.animation as animation
        from matplotlib.animation import PillowWriter
        
        print(f"Creating animation from {len(snapshot_files)} frames...")
        
        # First pass: determine global ranges for consistent axes
        all_dens, all_vel, all_pres, all_ene = [], [], [], []
        times = []
        
        for snap_file in snapshot_files:
            data = read_csv_snapshot(snap_file)
            if data is not None:
                all_dens.extend(data['dens'])
                all_vel.extend(data['vel_x'])
                all_pres.extend(data['pres'])
                all_ene.extend(data['ene'])
                # Extract time from filename (snapshot_XXXX.csv corresponds to t = XXXX * dt)
                # We'll read it from the data if available, otherwise use frame index
        
        dens_min, dens_max = min(all_dens), max(all_dens)
        vel_min, vel_max = min(all_vel), max(all_vel)
        pres_min, pres_max = min(all_pres), max(all_pres)
        ene_min, ene_max = min(all_ene), max(all_ene)
        
        # Add 5% margin
        dens_range = dens_max - dens_min
        vel_range = vel_max - vel_min
        ene_range = ene_max - ene_min
        
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        
        def update(frame_idx):
            data = read_csv_snapshot(snapshot_files[frame_idx])
            if data is None:
                return
            
            for ax in axs.flat:
                ax.clear()
            
            # Sort by position for cleaner visualization
            sort_idx = np.argsort(data['pos_x'])
            x_sorted = data['pos_x'][sort_idx]
            
            # Density
            axs[0, 0].scatter(data['pos_x'], data['dens'], s=3, alpha=0.7, c='blue', edgecolors='none')
            axs[0, 0].plot(x_sorted, data['dens'][sort_idx], 'b-', alpha=0.3, linewidth=0.5)
            axs[0, 0].set_ylabel('Density ρ', fontsize=11)
            axs[0, 0].set_ylim(dens_min - 0.05*dens_range, dens_max + 0.05*dens_range)
            axs[0, 0].grid(True, alpha=0.3)
            axs[0, 0].set_xlim(-0.5, 0.5)
            
            # Velocity
            axs[0, 1].scatter(data['pos_x'], data['vel_x'], s=3, alpha=0.7, c='red', edgecolors='none')
            axs[0, 1].plot(x_sorted, data['vel_x'][sort_idx], 'r-', alpha=0.3, linewidth=0.5)
            axs[0, 1].set_ylabel('Velocity v/c', fontsize=11)
            axs[0, 1].set_ylim(vel_min - 0.05*vel_range, vel_max + 0.05*vel_range)
            axs[0, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
            axs[0, 1].grid(True, alpha=0.3)
            axs[0, 1].set_xlim(-0.5, 0.5)
            
            # Pressure (log scale)
            axs[1, 0].scatter(data['pos_x'], data['pres'], s=3, alpha=0.7, c='green', edgecolors='none')
            axs[1, 0].plot(x_sorted, data['pres'][sort_idx], 'g-', alpha=0.3, linewidth=0.5)
            axs[1, 0].set_ylabel('Pressure P', fontsize=11)
            axs[1, 0].set_yscale('log')
            axs[1, 0].set_ylim(pres_min*0.8, pres_max*1.2)
            axs[1, 0].grid(True, alpha=0.3, which='both')
            axs[1, 0].set_xlim(-0.5, 0.5)
            axs[1, 0].set_xlabel('Position x', fontsize=11)
            
            # Internal Energy
            axs[1, 1].scatter(data['pos_x'], data['ene'], s=3, alpha=0.7, c='purple', edgecolors='none')
            axs[1, 1].plot(x_sorted, data['ene'][sort_idx], color='purple', alpha=0.3, linewidth=0.5)
            axs[1, 1].set_ylabel('Internal Energy u', fontsize=11)
            axs[1, 1].set_ylim(ene_min - 0.05*ene_range, ene_max + 0.05*ene_range)
            axs[1, 1].grid(True, alpha=0.3)
            axs[1, 1].set_xlim(-0.5, 0.5)
            axs[1, 1].set_xlabel('Position x', fontsize=11)
            
            # Extract time information (snapshot number * output interval)
            snapshot_num = int(os.path.basename(snapshot_files[frame_idx]).split('_')[1].split('.')[0])
            time_estimate = snapshot_num * 0.01  # Assuming dt_output = 0.01 from JSON
            
            fig.suptitle(f'SR-GSPH Sod Shock Tube - t ≈ {time_estimate:.3f} (Frame {frame_idx+1}/{len(snapshot_files)})', 
                        fontsize=14, fontweight='bold')
            
            plt.tight_layout()
        
        # Create animation
        ani = animation.FuncAnimation(fig, update, frames=len(snapshot_files), 
                                     interval=200, repeat=True, blit=False)
        
        # Save as GIF
        output_file = os.path.join(results_dir, f'{output_prefix}_animation.gif')
        writer = PillowWriter(fps=5, bitrate=1800)
        ani.save(output_file, writer=writer, dpi=100)
        print(f"✓ Animation saved: {output_file}")
        plt.close()
        
    except ImportError as e:
        print(f"Warning: Animation requires Pillow package: {e}")
        print("Install with: pip install pillow")
    except Exception as e:
        print(f"Animation creation failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_sr_profiles.py <results_directory> [output_prefix]")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'sr_profiles'
    
    if not os.path.isdir(results_dir):
        print(f"Error: {results_dir} is not a directory")
        sys.exit(1)
    
    plot_sr_profiles(results_dir, output_prefix)
