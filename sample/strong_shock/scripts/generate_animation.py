#!/usr/bin/env python3
"""
Generate Animation for Strong Shock Test Results

Creates an animated GIF from snapshot CSV files.

Usage:
    python3 generate_animation.py RESULTS_DIR [-o OUTPUT] [--fps FPS]
    
Arguments:
    RESULTS_DIR - Directory containing snapshot CSV files
    -o, --output - Output GIF filename (default: strong_shock_animation.gif)
    --fps - Frames per second (default: 5)
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import argparse
from pathlib import Path

# Animation imports
try:
    from matplotlib.animation import FuncAnimation, PillowWriter
    HAS_ANIMATION = True
except ImportError:
    HAS_ANIMATION = False
    print("ERROR: Animation requires pillow")
    print("Install with: pip install pillow")
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate animation from strong shock snapshots')
    parser.add_argument('results_dir', type=str, help='Directory containing snapshot files')
    parser.add_argument('-o', '--output', type=str, default=None, help='Output GIF filename')
    parser.add_argument('--fps', type=int, default=5, help='Frames per second')
    return parser.parse_args()

def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    with open(filename, 'r') as f:
        # Read metadata
        metadata = {}
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        # Read data
        f.seek(0)
        lines = [l for l in f.readlines() if not l.startswith('#')]
        header = lines[0].strip().split(',')
        
        for col_name in header:
            data[col_name] = []
        
        for line in lines[1:]:
            values = line.strip().split(',')
            for i, col_name in enumerate(header):
                try:
                    data[col_name].append(float(values[i]))
                except (ValueError, IndexError):
                    pass
    
    # Convert to numpy arrays
    for key in data:
        data[key] = np.array(data[key])
    
    data['metadata'] = metadata
    return data

def create_animation(snapshot_files, output_file, fps=5):
    """Create animation from snapshot files"""
    
    print(f'Creating animation with {len(snapshot_files)} frames...')
    
    # Load all data
    all_data = []
    for f in snapshot_files:
        all_data.append(load_snapshot(f))
    
    # Determine global ranges for consistent axes
    # Use full simulation domain for x-axis: [-0.5, 0.5] from C++ strong_shock.cpp
    x_range = [-0.5, 0.5]
    
    all_dens = np.concatenate([d['dens'] for d in all_data])
    all_vel = np.concatenate([d['vel_x'] for d in all_data])
    all_pres = np.concatenate([d['pres'] for d in all_data])
    all_ene = np.concatenate([d['ene'] for d in all_data])
    
    # Filter out NaN and Inf values for range calculation
    all_dens = all_dens[np.isfinite(all_dens)]
    all_vel = all_vel[np.isfinite(all_vel)]
    all_pres = all_pres[np.isfinite(all_pres)]
    all_ene = all_ene[np.isfinite(all_ene)]
    
    dens_range = [all_dens.min() * 0.95, all_dens.max() * 1.05]
    vel_range = [all_vel.min() - 0.1, all_vel.max() + 0.1]
    pres_range = [all_pres.min() * 0.95, all_pres.max() * 1.05]
    ene_range = [all_ene.min() * 0.95, all_ene.max() * 1.05]
    
    # Setup figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Initialize empty plots
    line1, = ax1.plot([], [], 'o-', color='#0173B2', linewidth=2, markersize=3)
    line2, = ax2.plot([], [], 'o-', color='#DE8F05', linewidth=2, markersize=3)
    line3, = ax3.plot([], [], 'o-', color='#029E73', linewidth=2, markersize=3)
    line4, = ax4.plot([], [], 'o-', color='#CC78BC', linewidth=2, markersize=3)
    
    # Format axes
    ax1.set_xlim(x_range)
    ax1.set_ylim(dens_range)
    ax1.set_ylabel('Density ρ', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Position x', fontsize=12)
    ax1.set_title('Density Profile', fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':')
    
    ax2.set_xlim(x_range)
    ax2.set_ylim(vel_range)
    ax2.set_ylabel('Velocity u', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Position x', fontsize=12)
    ax2.set_title('Velocity Profile', fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle=':')
    
    ax3.set_xlim(x_range)
    ax3.set_ylim(pres_range)
    ax3.set_ylabel('Pressure P', fontsize=12, fontweight='bold')
    ax3.set_xlabel('Position x', fontsize=12)
    ax3.set_title('Pressure Profile', fontweight='bold')
    ax3.grid(True, alpha=0.3, linestyle=':')
    
    ax4.set_xlim(x_range)
    ax4.set_ylim(ene_range)
    ax4.set_ylabel('Internal Energy e', fontsize=12, fontweight='bold')
    ax4.set_xlabel('Position x', fontsize=12)
    ax4.set_title('Internal Energy Profile', fontweight='bold')
    ax4.grid(True, alpha=0.3, linestyle=':')
    
    # Main title (will be updated per frame)
    title = fig.suptitle('', fontsize=16, fontweight='bold')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    def update(frame):
        """Update function for animation"""
        data = all_data[frame]
        
        # Sort by x position
        x = data['pos_x']
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        
        # Update plots
        line1.set_data(x_sorted, data['dens'][sort_idx])
        line2.set_data(x_sorted, data['vel_x'][sort_idx])
        line3.set_data(x_sorted, data['pres'][sort_idx])
        line4.set_data(x_sorted, data['ene'][sort_idx])
        
        # Update title with time
        time_str = data['metadata'].get('time', f'frame {frame}')
        title.set_text(f'Strong Shock Test (P_ratio = 10,000) | Time: {time_str}')
        
        return line1, line2, line3, line4, title
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(all_data), 
                        interval=1000//fps, blit=True)
    
    # Save
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer)
    plt.close()
    
    print(f'✓ Animation saved: {output_file}')
    print(f'  Frames: {len(all_data)}')
    print(f'  FPS: {fps}')
    print(f'  Duration: {len(all_data)/fps:.1f} seconds')

def main():
    args = parse_args()
    
    # Find snapshot files
    snapshot_files = sorted(glob.glob(f'{args.results_dir}/snapshot_*.csv'))
    
    if not snapshot_files:
        print(f'ERROR: No snapshot files found in {args.results_dir}')
        sys.exit(1)
    
    print('=' * 70)
    print('Strong Shock Animation Generator')
    print('=' * 70)
    print(f'Input directory: {args.results_dir}')
    print(f'Found {len(snapshot_files)} snapshots')
    print()
    
    # Determine output filename
    if args.output is None:
        method_name = os.path.basename(args.results_dir.rstrip('/'))
        output_file = f'{args.results_dir}/{method_name}_animation.gif'
    else:
        output_file = args.output
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    
    # Generate animation
    create_animation(snapshot_files, output_file, args.fps)
    
    print()
    print('=' * 70)
    print('✓ Animation Complete!')
    print('=' * 70)

if __name__ == '__main__':
    main()
