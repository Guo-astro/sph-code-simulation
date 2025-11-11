#!/usr/bin/env python3
"""
Kernel Comparison Animation for Gresho-Chan Vortex
Creates animated comparison between Wendland and Cubic Spline kernels
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import os
import sys
from pathlib import Path
from tqdm import tqdm

def load_snapshot(filepath):
    """Load a snapshot CSV file and extract metadata."""
    metadata = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
    
    data = pd.read_csv(filepath, comment='#')
    return data, metadata

def main():
    base_dir = 'sample/gresho_chan_vortex/results'
    output_file = os.path.join(base_dir, 'kernel_comparison', 'kernel_comparison.gif')
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
    method_labels = ['GSPH', 'SSPH', 'DISPH', 'GDISPH', 'GDISPH+Balsara']
    
    print("=" * 63)
    print("       Kernel Comparison Animation - Wendland vs Cubic Spline")
    print("=" * 63)
    print(f"Base directory: {base_dir}")
    print(f"Output file:    {output_file}")
    print(f"Methods:        {', '.join([m.upper() for m in methods])}")
    print()
    
    # Find common snapshots
    print("Scanning for snapshot files...")
    all_snapshot_numbers = {}
    
    for method in methods:
        wendland_dir = os.path.join(base_dir, f"{method}_wendland")
        cubic_dir = os.path.join(base_dir, method)
        
        wendland_snaps = set()
        cubic_snaps = set()
        
        if os.path.exists(wendland_dir):
            for f in os.listdir(wendland_dir):
                if f.startswith('snapshot_') and f.endswith('.csv'):
                    snap_num = int(f.split('_')[1].split('.')[0])
                    wendland_snaps.add(snap_num)
        
        if os.path.exists(cubic_dir):
            for f in os.listdir(cubic_dir):
                if f.startswith('snapshot_') and f.endswith('.csv'):
                    snap_num = int(f.split('_')[1].split('.')[0])
                    cubic_snaps.add(snap_num)
        
        common = wendland_snaps & cubic_snaps
        all_snapshot_numbers[method] = common
        print(f"  {method.upper():18s}: {len(cubic_snaps)} cubic, {len(wendland_snaps)} wendland, {len(common)} common")
    
    common_numbers = set.intersection(*all_snapshot_numbers.values())
    common_numbers = sorted(list(common_numbers))
    
    if not common_numbers:
        print("\n❌ Error: No common snapshots found!")
        return 1
    
    print(f"\nFound {len(common_numbers)} common snapshot numbers across all methods and kernels")
    
    # Animation frames (every snapshot)
    frame_numbers = common_numbers
    print(f"Animation frames: {len(frame_numbers)} (every 1 snapshots)")
    print()
    
    # Pre-load all data
    print("Pre-loading data...")
    frames_data = []
    
    for snap_num in tqdm(frame_numbers):
        frame_data = {'snap_num': snap_num, 'wendland': {}, 'cubic': {}}
        
        for method in methods:
            # Load Wendland
            w_file = os.path.join(base_dir, f"{method}_wendland", f"snapshot_{snap_num:04d}.csv")
            if os.path.exists(w_file):
                data_w, meta_w = load_snapshot(w_file)
                frame_data['wendland'][method] = (data_w, meta_w)
            
            # Load Cubic
            c_file = os.path.join(base_dir, method, f"snapshot_{snap_num:04d}.csv")
            if os.path.exists(c_file):
                data_c, meta_c = load_snapshot(c_file)
                frame_data['cubic'][method] = (data_c, meta_c)
                frame_data['time'] = float(meta_c.get('Time', 0))
        
        frames_data.append(frame_data)
    
    print(f"Loaded {len(frames_data)} frames")
    print()
    
    # Determine density range for consistent coloring
    all_densities = []
    for frame in frames_data:
        for kernel_type in ['wendland', 'cubic']:
            for method in methods:
                if method in frame[kernel_type]:
                    data, _ = frame[kernel_type][method]
                    if 'dens' in data:
                        all_densities.extend(data['dens'].values)
    
    all_densities = np.array(all_densities)
    log_dens_min = np.log10(np.percentile(all_densities, 1) + 1e-10)
    log_dens_max = np.log10(np.percentile(all_densities, 99) + 1e-10)
    
    print(f"Density range: log₁₀ density from {log_dens_min:.2f} to {log_dens_max:.2f}")
    print()
    
    # Create animation
    fig, axes = plt.subplots(5, 2, figsize=(12, 24))
    
    def init():
        for row in axes:
            for ax in row:
                ax.clear()
        return []
    
    def animate(frame_idx):
        frame = frames_data[frame_idx]
        snap_num = frame['snap_num']
        time = frame.get('time', 0)
        
        fig.suptitle(f'Kernel Comparison: Wendland vs Cubic Spline - t={time:.3f}',
                     fontsize=16, fontweight='bold', y=0.995)
        
        for idx, method in enumerate(methods):
            # Wendland (left column)
            ax = axes[idx, 0]
            ax.clear()
            
            if method in frame['wendland']:
                data_w, _ = frame['wendland'][method]
                x_col = 'x' if 'x' in data_w else 'pos_x'
                y_col = 'y' if 'y' in data_w else 'pos_y'
                
                if x_col in data_w and y_col in data_w and 'dens' in data_w:
                    log_dens = np.log10(data_w['dens'] + 1e-10)
                    ax.scatter(data_w[x_col], data_w[y_col], c=log_dens,
                             cmap='hot', s=15, alpha=0.8, edgecolors='none',
                             vmin=log_dens_min, vmax=log_dens_max)
                    ax.set_ylabel(method_labels[idx], fontsize=11, fontweight='bold',
                                 rotation=0, ha='right', va='center', labelpad=45)
                    ax.set_aspect('equal')
                    ax.set_xlim(-0.5, 0.5)
                    ax.set_ylim(-0.5, 0.5)
                    ax.set_facecolor('black')
                    ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
                    
                    if idx == 0:
                        ax.set_title('Wendland Kernel', fontsize=12, fontweight='bold', pad=15)
                    if idx == len(methods) - 1:
                        ax.set_xlabel('x', fontsize=10)
                    else:
                        ax.set_xticklabels([])
            
            # Cubic Spline (right column)
            ax = axes[idx, 1]
            ax.clear()
            
            if method in frame['cubic']:
                data_c, _ = frame['cubic'][method]
                x_col = 'x' if 'x' in data_c else 'pos_x'
                y_col = 'y' if 'y' in data_c else 'pos_y'
                
                if x_col in data_c and y_col in data_c and 'dens' in data_c:
                    log_dens = np.log10(data_c['dens'] + 1e-10)
                    ax.scatter(data_c[x_col], data_c[y_col], c=log_dens,
                             cmap='hot', s=15, alpha=0.8, edgecolors='none',
                             vmin=log_dens_min, vmax=log_dens_max)
                    ax.set_aspect('equal')
                    ax.set_xlim(-0.5, 0.5)
                    ax.set_ylim(-0.5, 0.5)
                    ax.set_facecolor('black')
                    ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
                    
                    if idx == 0:
                        ax.set_title('Cubic Spline Kernel', fontsize=12, fontweight='bold', pad=15)
                    if idx == len(methods) - 1:
                        ax.set_xlabel('x', fontsize=10)
                    else:
                        ax.set_xticklabels([])
                    
                    ax.set_yticklabels([])
        
        plt.tight_layout(rect=[0, 0, 1, 0.99])
        return []
    
    print("Generating animation...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(frames_data),
                                   interval=100, blit=False, repeat=True)
    
    print("Saving animation...")
    anim.save(output_file, writer='pillow', fps=10, dpi=80)
    plt.close()
    
    print()
    print("=" * 63)
    print("       Animation Complete!")
    print("=" * 63)
    print(f"✓ Saved: {output_file}")
    print("=" * 63)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
