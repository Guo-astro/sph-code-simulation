#!/usr/bin/env python3
"""
Enhanced animation script for SR-GSPH Sod shock tube
Creates high-quality animations with multiple visualization options
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter, FFMpegWriter
import glob
import os
import sys
from pathlib import Path

def read_csv_snapshot(filename):
    """Read SPH snapshot from CSV file"""
    try:
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                if not line.startswith('#'):
                    header_end = i
                    break
        
        data = np.genfromtxt(filename, delimiter=',', names=True, skip_header=header_end)
        return data
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

def create_enhanced_animation(results_dir, output_name='sr_sod_animation', format='gif', fps=10):
    """
    Create enhanced animation with better visualization
    
    Parameters:
    -----------
    results_dir : str
        Directory containing snapshot CSV files
    output_name : str
        Base name for output file (without extension)
    format : str
        Output format: 'gif' or 'mp4'
    fps : int
        Frames per second
    """
    
    # Find all snapshot files
    snapshot_files = sorted(glob.glob(os.path.join(results_dir, 'snapshot_*.csv')))
    
    if len(snapshot_files) < 2:
        print(f"Error: Need at least 2 snapshots, found {len(snapshot_files)}")
        return
    
    print(f"Creating animation from {len(snapshot_files)} snapshots...")
    
    # First pass: determine global ranges
    all_data = []
    for snap_file in snapshot_files:
        data = read_csv_snapshot(snap_file)
        if data is not None:
            all_data.append(data)
    
    if not all_data:
        print("Error: Failed to read any snapshots")
        return
    
    # Combine all data to find global ranges
    all_dens = np.concatenate([d['dens'] for d in all_data])
    all_vel = np.concatenate([d['vel_x'] for d in all_data])
    all_pres = np.concatenate([d['pres'] for d in all_data])
    all_ene = np.concatenate([d['ene'] for d in all_data])
    
    # Calculate ranges with margins
    dens_min, dens_max = all_dens.min(), all_dens.max()
    vel_min, vel_max = all_vel.min(), all_vel.max()
    pres_min, pres_max = all_pres.min(), all_pres.max()
    ene_min, ene_max = all_ene.min(), all_ene.max()
    
    # Add margins
    dens_margin = (dens_max - dens_min) * 0.1
    vel_margin = max((vel_max - vel_min) * 0.1, 0.01)  # Minimum margin for velocity
    ene_margin = (ene_max - ene_min) * 0.1
    
    # Create figure
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.25)
    
    ax_dens = fig.add_subplot(gs[0, 0])
    ax_vel = fig.add_subplot(gs[0, 1])
    ax_pres = fig.add_subplot(gs[1, 0])
    ax_ene = fig.add_subplot(gs[1, 1])
    ax_gamma = fig.add_subplot(gs[2, 0])
    ax_phase = fig.add_subplot(gs[2, 1])
    
    def update(frame_idx):
        """Update function for animation"""
        data = all_data[frame_idx]
        
        # Clear all axes
        for ax in [ax_dens, ax_vel, ax_pres, ax_ene, ax_gamma, ax_phase]:
            ax.clear()
        
        # Sort by position for better visualization
        sort_idx = np.argsort(data['pos_x'])
        x = data['pos_x'][sort_idx]
        
        # 1. Density profile
        ax_dens.scatter(data['pos_x'], data['dens'], s=4, alpha=0.7, c='blue', edgecolors='none')
        ax_dens.plot(x, data['dens'][sort_idx], 'b-', alpha=0.3, linewidth=1)
        ax_dens.set_ylabel('Rest-frame density ρ', fontsize=12, fontweight='bold')
        ax_dens.set_ylim(dens_min - dens_margin, dens_max + dens_margin)
        ax_dens.grid(True, alpha=0.3)
        ax_dens.set_xlim(-0.5, 0.5)
        
        # 2. Velocity profile
        ax_vel.scatter(data['pos_x'], data['vel_x'], s=4, alpha=0.7, c='red', edgecolors='none')
        ax_vel.plot(x, data['vel_x'][sort_idx], 'r-', alpha=0.3, linewidth=1)
        ax_vel.set_ylabel('Velocity v/c', fontsize=12, fontweight='bold')
        ax_vel.set_ylim(vel_min - vel_margin, vel_max + vel_margin)
        ax_vel.axhline(y=0, color='k', linestyle='--', alpha=0.4, linewidth=1)
        ax_vel.grid(True, alpha=0.3)
        ax_vel.set_xlim(-0.5, 0.5)
        
        # 3. Pressure profile (log scale)
        ax_pres.scatter(data['pos_x'], data['pres'], s=4, alpha=0.7, c='green', edgecolors='none')
        ax_pres.plot(x, data['pres'][sort_idx], 'g-', alpha=0.3, linewidth=1)
        ax_pres.set_ylabel('Pressure P', fontsize=12, fontweight='bold')
        ax_pres.set_yscale('log')
        ax_pres.set_ylim(pres_min*0.7, pres_max*1.3)
        ax_pres.grid(True, alpha=0.3, which='both')
        ax_pres.set_xlim(-0.5, 0.5)
        ax_pres.set_xlabel('Position x', fontsize=11)
        
        # 4. Internal energy
        ax_ene.scatter(data['pos_x'], data['ene'], s=4, alpha=0.7, c='purple', edgecolors='none')
        ax_ene.plot(x, data['ene'][sort_idx], color='purple', alpha=0.3, linewidth=1)
        ax_ene.set_ylabel('Internal Energy u', fontsize=12, fontweight='bold')
        ax_ene.set_ylim(ene_min - ene_margin, ene_max + ene_margin)
        ax_ene.grid(True, alpha=0.3)
        ax_ene.set_xlim(-0.5, 0.5)
        ax_ene.set_xlabel('Position x', fontsize=11)
        
        # 5. Lorentz factor (if available)
        if 'gamma_lor' in data.dtype.names:
            ax_gamma.scatter(data['pos_x'], data['gamma_lor'], s=4, alpha=0.7, 
                           c='orange', edgecolors='none')
            ax_gamma.plot(x, data['gamma_lor'][sort_idx], 'orange', alpha=0.3, linewidth=1)
            ax_gamma.set_ylabel('Lorentz factor Γ', fontsize=12, fontweight='bold')
            ax_gamma.axhline(y=1, color='k', linestyle='--', alpha=0.4, linewidth=1)
        else:
            # Calculate Lorentz factor from velocity
            gamma_lor = 1.0 / np.sqrt(1.0 - data['vel_x']**2)
            ax_gamma.scatter(data['pos_x'], gamma_lor, s=4, alpha=0.7, 
                           c='orange', edgecolors='none')
            gamma_sorted = 1.0 / np.sqrt(1.0 - data['vel_x'][sort_idx]**2)
            ax_gamma.plot(x, gamma_sorted, 'orange', alpha=0.3, linewidth=1)
            ax_gamma.set_ylabel('Lorentz factor Γ', fontsize=12, fontweight='bold')
            ax_gamma.axhline(y=1, color='k', linestyle='--', alpha=0.4, linewidth=1)
        ax_gamma.grid(True, alpha=0.3)
        ax_gamma.set_xlim(-0.5, 0.5)
        ax_gamma.set_xlabel('Position x', fontsize=11)
        
        # 6. Phase diagram: Pressure vs Density
        ax_phase.scatter(data['dens'], data['pres'], s=4, alpha=0.6, 
                        c=data['pos_x'], cmap='viridis', edgecolors='none')
        ax_phase.set_xlabel('Density ρ', fontsize=11)
        ax_phase.set_ylabel('Pressure P', fontsize=11)
        ax_phase.set_xscale('log')
        ax_phase.set_yscale('log')
        ax_phase.set_title('Phase Diagram (colored by position)', fontsize=10)
        ax_phase.grid(True, alpha=0.3, which='both')
        
        # Extract time from snapshot number
        snapshot_num = int(os.path.basename(snapshot_files[frame_idx]).split('_')[1].split('.')[0])
        time_est = snapshot_num * 0.01  # dt_output from JSON
        
        # Main title
        fig.suptitle(f'SR-GSPH Sod Shock Tube Evolution\nt ≈ {time_est:.3f} | Frame {frame_idx+1}/{len(all_data)}',
                    fontsize=16, fontweight='bold')
    
    # Create animation
    print("Generating frames...")
    ani = animation.FuncAnimation(fig, update, frames=len(all_data), 
                                 interval=1000//fps, repeat=True, blit=False)
    
    # Save animation
    if format.lower() == 'gif':
        output_file = os.path.join(results_dir, f'{output_name}.gif')
        writer = PillowWriter(fps=fps, bitrate=1800)
        ani.save(output_file, writer=writer, dpi=120)
        print(f"✓ GIF animation saved: {output_file}")
        print(f"  Size: {os.path.getsize(output_file) / 1024:.1f} KB")
        
    elif format.lower() == 'mp4':
        output_file = os.path.join(results_dir, f'{output_name}.mp4')
        try:
            writer = FFMpegWriter(fps=fps, bitrate=1800, codec='libx264')
            ani.save(output_file, writer=writer, dpi=120)
            print(f"✓ MP4 animation saved: {output_file}")
            print(f"  Size: {os.path.getsize(output_file) / 1024:.1f} KB")
        except Exception as e:
            print(f"Error creating MP4 (FFmpeg required): {e}")
            print("Falling back to GIF format...")
            output_file = os.path.join(results_dir, f'{output_name}.gif')
            writer = PillowWriter(fps=fps, bitrate=1800)
            ani.save(output_file, writer=writer, dpi=120)
            print(f"✓ GIF animation saved: {output_file}")
    
    plt.close()
    return output_file

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python create_sr_animation.py <results_directory> [output_name] [format] [fps]")
        print("\nArguments:")
        print("  results_directory : Directory containing snapshot_*.csv files")
        print("  output_name       : Output filename without extension (default: 'sr_sod_animation')")
        print("  format            : 'gif' or 'mp4' (default: 'gif')")
        print("  fps               : Frames per second (default: 10)")
        print("\nExample:")
        print("  python create_sr_animation.py sample/sr_sod/results sr_shock gif 10")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_name = sys.argv[2] if len(sys.argv) > 2 else 'sr_sod_animation'
    format_type = sys.argv[3] if len(sys.argv) > 3 else 'gif'
    fps = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    
    if not os.path.isdir(results_dir):
        print(f"Error: {results_dir} is not a directory")
        sys.exit(1)
    
    try:
        output_file = create_enhanced_animation(results_dir, output_name, format_type, fps)
        print(f"\n✅ Animation complete!")
        print(f"📁 Location: {output_file}")
    except Exception as e:
        print(f"\n❌ Animation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
