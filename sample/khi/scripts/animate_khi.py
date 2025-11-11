#!/usr/bin/env python3
"""
Create Side-by-Side Animation Comparing Pairing Instability Results

Generates an animated GIF showing time evolution of all five SPH methods simultaneously.

Usage:
    python3 animate_pairing.py [base_dir] [output_file]
    
Arguments:
    base_dir    - Base directory containing gsph/, ssph/, disph/, gdisph/, gdisph_balsara/ subdirs
                  Default: sample/pairing_instability/results
    output_file - Output animation filename
                  Default: sample/pairing_instability/results/comparison/comparison_animation.gif
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
from pathlib import Path

# Try to import animation tools
try:
    from matplotlib.animation import FuncAnimation, PillowWriter
    HAS_ANIMATION = True
except ImportError:
    HAS_ANIMATION = False
    print("ERROR: Animation tools not available")
    print("Install with: pip install pillow")
    sys.exit(1)

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# Configuration
base_dir = sys.argv[1] if len(sys.argv) > 1 else "sample/pairing_instability/results"
output_file = sys.argv[2] if len(sys.argv) > 2 else "sample/pairing_instability/results/comparison/comparison_animation.gif"

# SPH methods to compare
methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
method_labels = ['GSPH (Godunov)', 'SSPH (Standard)', 'DISPH (Density-Independent)', 
                 'GDISPH (Godunov DISPH)', 'GDISPH+Balsara']
method_colors = ['#0173B2', '#DE8F05', '#029E73', '#CC78BC', '#D55E00']

print('=' * 70)
print('Pairing Instability Animation - Multi-Method Comparison')
print('=' * 70)
print(f'Base directory: {base_dir}')
print(f'Output file:    {output_file}')
print(f'Methods:        {", ".join(method_labels)}')
print()

def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    metadata = {}
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
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
    
    for key in data:
        data[key] = np.array(data[key])
    
    return data, metadata

# Scan for snapshot files
print("Scanning for snapshot files...")
method_snapshots = {}
for method in methods:
    method_dir = os.path.join(base_dir, method)
    if not os.path.exists(method_dir):
        print(f"  WARNING: {method.upper()} directory not found")
        continue
    
    snapshots = sorted(glob.glob(os.path.join(method_dir, "snapshot_*.csv")))
    if not snapshots:
        print(f"  WARNING: No snapshots found for {method.upper()}")
        continue
    
    snapshot_numbers = []
    for snap in snapshots:
        basename = os.path.basename(snap)
        num_str = basename.replace('snapshot_', '').replace('.csv', '')
        try:
            snapshot_numbers.append(int(num_str))
        except ValueError:
            pass
    
    method_snapshots[method] = {
        'files': snapshots,
        'numbers': snapshot_numbers
    }
    print(f"  {method.upper():15s}: {len(snapshot_numbers)} snapshots")

if not method_snapshots:
    print("ERROR: No snapshot data found!")
    sys.exit(1)

# Find common snapshots
common_numbers = set(method_snapshots[methods[0]]['numbers'])
for method in methods[1:]:
    if method in method_snapshots:
        common_numbers &= set(method_snapshots[method]['numbers'])

common_numbers = sorted(common_numbers)
print(f"\nFound {len(common_numbers)} common snapshots")

# Animation settings
frame_skip = max(1, len(common_numbers) // 50)  # Target ~50 frames
animation_snapshots = common_numbers[::frame_skip]
print(f"Animation frames: {len(animation_snapshots)} (every {frame_skip} snapshots)\n")

# Pre-load all data
print("Pre-loading data...")
all_frames = []

iterator = tqdm(animation_snapshots) if HAS_TQDM else animation_snapshots

for snap_num in iterator:
    frame_data = {}
    current_time = None
    
    for method in methods:
        if method not in method_snapshots:
            continue
        
        snap_file = os.path.join(base_dir, method, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, metadata = load_snapshot(snap_file)
            frame_data[method] = data
            if current_time is None and 'time' in metadata:
                current_time = metadata['time']
    
    all_frames.append({
        'number': snap_num,
        'time': current_time,
        'data': frame_data
    })

print(f"Loaded {len(all_frames)} frames\n")

# Determine global density range for consistent coloring (use log scale)
all_densities = []
for frame in all_frames:
    for method, data in frame['data'].items():
        if 'dens' in data:
            all_densities.extend(data['dens'])

# Convert to log scale for better dynamic range
all_log_dens = np.log10(np.array(all_densities) + 1e-10)
vmin = np.percentile(all_log_dens, 1)
vmax = np.percentile(all_log_dens, 99)

# Create animation
print("Generating animation...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Pairing Instability - 5 Method Comparison', fontsize=16, fontweight='bold')

# Initialize scatter plots
scatters = {}
for idx, method in enumerate(methods):
    ax = axes.flatten()[idx]
    ax.set_xlabel('x', fontsize=10)
    ax.set_ylabel('y', fontsize=10)
    ax.set_title(method_labels[idx], fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_facecolor('black')  # Black background emphasizes hot spots
    ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
    
    scatters[method] = ax.scatter([], [], c=[], cmap='hot', s=20, alpha=0.8, edgecolors='none',
                                 vmin=vmin, vmax=vmax)
    
    # Add colorbar to first plot only
    if idx == 0:
        cbar = plt.colorbar(scatters[method], ax=ax)
        cbar.set_label('log₁₀(Density)', fontsize=9)

# Info panel
info_ax = axes.flatten()[5]
info_ax.axis('off')
time_text = info_ax.text(0.5, 0.5, '', ha='center', va='center', 
                        fontsize=14, fontweight='bold', transform=info_ax.transAxes)

def update(frame_idx):
    frame = all_frames[frame_idx]
    
    for method in methods:
        if method in frame['data']:
            data = frame['data'][method]
            # Handle both 'x','y' and 'pos_x','pos_y' column naming
            x_col = 'x' if 'x' in data else 'pos_x'
            y_col = 'y' if 'y' in data else 'pos_y'
            if x_col in data and y_col in data and 'dens' in data:
                # Convert density to log scale for better visualization
                dens_array = np.array(data['dens'])
                log_dens = np.log10(dens_array + 1e-10)
                scatters[method].set_offsets(np.c_[data[x_col], data[y_col]])
                scatters[method].set_array(log_dens)
    
    if frame['time']:
        time_text.set_text(f"Time: {frame['time']}\nSnapshot: {frame['number']:04d}")
    else:
        time_text.set_text(f"Snapshot: {frame['number']:04d}")
    
    return list(scatters.values()) + [time_text]

# Create animation
anim = FuncAnimation(fig, update, frames=tqdm(range(len(all_frames))) if HAS_TQDM else range(len(all_frames)),
                    interval=100, blit=False, repeat=True)

# Save animation
writer = PillowWriter(fps=10)
anim.save(output_file, writer=writer, dpi=100)

plt.close()

print()
print('=' * 70)
print('Animation Complete!')
print('=' * 70)
print(f'✓ Saved: {output_file}')
print('=' * 70)
