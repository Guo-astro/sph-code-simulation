#!/usr/bin/env python3
"""
Create Animation Comparing Wendland vs Cubic Spline Kernels

Generates an animated GIF showing time evolution of all five SPH methods
with side-by-side kernel comparison (Wendland left, Cubic Spline right).

Usage:
    python3 animate_kernel_comparison.py [base_dir] [output_file]
    
Arguments:
    base_dir    - Base directory containing method subdirs (gsph/, ssph/, etc.)
                  Default: sample/pairing_instability/results
    output_file - Output animation filename
                  Default: sample/pairing_instability/results/kernel_comparison/kernel_comparison.gif
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
output_file = sys.argv[2] if len(sys.argv) > 2 else "sample/pairing_instability/results/kernel_comparison/kernel_comparison.gif"

# SPH methods to compare
methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
method_labels = ['GSPH\n(Godunov)', 'SSPH\n(Standard)', 'DISPH\n(Density-Indep)', 
                 'GDISPH\n(Godunov DISPH)', 'GDISPH+Balsara']

print('=' * 70)
print('Kernel Comparison Animation - Wendland vs Cubic Spline')
print('=' * 70)
print(f'Base directory: {base_dir}')
print(f'Output file:    {output_file}')
print(f'Methods:        {", ".join([m.upper() for m in methods])}')
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
    # Cubic Spline snapshots
    cubic_dir = os.path.join(base_dir, method)
    cubic_numbers = []
    if os.path.exists(cubic_dir):
        cubic_snaps = sorted(glob.glob(os.path.join(cubic_dir, "snapshot_*.csv")))
        for snap in cubic_snaps:
            basename = os.path.basename(snap)
            num_str = basename.replace('snapshot_', '').replace('.csv', '')
            try:
                cubic_numbers.append(int(num_str))
            except ValueError:
                pass
    
    # Wendland snapshots
    wendland_dir = os.path.join(base_dir, f"{method}_wendland")
    wendland_numbers = []
    if os.path.exists(wendland_dir):
        wendland_snaps = sorted(glob.glob(os.path.join(wendland_dir, "snapshot_*.csv")))
        for snap in wendland_snaps:
            basename = os.path.basename(snap)
            num_str = basename.replace('snapshot_', '').replace('.csv', '')
            try:
                wendland_numbers.append(int(num_str))
            except ValueError:
                pass
    
    # Find common snapshots between both kernels
    common = sorted(set(cubic_numbers) & set(wendland_numbers))
    method_snapshots[method] = common
    print(f"  {method.upper():<18}: {len(cubic_numbers)} cubic, {len(wendland_numbers)} wendland, {len(common)} common")

if not method_snapshots:
    print("ERROR: No snapshot data found!")
    sys.exit(1)

# Find common snapshots across all methods
common_numbers = sorted(set.intersection(*[set(nums) for nums in method_snapshots.values()]))
print(f"\nFound {len(common_numbers)} common snapshot numbers across all methods and kernels")

# Animation settings
frame_skip = max(1, len(common_numbers) // 50)  # Target ~50 frames
animation_snapshots = common_numbers[::frame_skip]
print(f"Animation frames: {len(animation_snapshots)} (every {frame_skip} snapshots)\n")

# Pre-load all data
print("Pre-loading data...")
all_frames = []

iterator = tqdm(animation_snapshots) if HAS_TQDM else animation_snapshots

for snap_num in iterator:
    frame_data = {'wendland': {}, 'cubic': {}}
    current_time = None
    
    for method in methods:
        # Load Wendland
        wendland_file = os.path.join(base_dir, f"{method}_wendland", f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(wendland_file):
            data, metadata = load_snapshot(wendland_file)
            frame_data['wendland'][method] = data
            if current_time is None and 'time' in metadata:
                current_time = metadata['time']
        
        # Load Cubic Spline
        cubic_file = os.path.join(base_dir, method, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(cubic_file):
            data, metadata = load_snapshot(cubic_file)
            frame_data['cubic'][method] = data
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
    for kernel_type in ['wendland', 'cubic']:
        for method, data in frame['data'][kernel_type].items():
            if 'dens' in data:
                all_densities.extend(data['dens'])

# Convert to log scale for better dynamic range
all_log_dens = np.log10(np.array(all_densities) + 1e-10)
vmin = np.percentile(all_log_dens, 1)
vmax = np.percentile(all_log_dens, 99)

print(f"Density range: log₁₀ density from {vmin:.2f} to {vmax:.2f}\n")

# Create animation - 5 rows (methods) × 2 columns (kernels)
print("Generating animation...")

fig, axes = plt.subplots(5, 2, figsize=(12, 24))
fig.suptitle('Kernel Comparison: Wendland (left) vs Cubic Spline (right)', 
             fontsize=16, fontweight='bold', y=0.995)

# Column titles
axes[0, 0].set_title('Wendland Kernel', fontsize=14, fontweight='bold', pad=20)
axes[0, 1].set_title('Cubic Spline Kernel', fontsize=14, fontweight='bold', pad=20)

# Initialize scatter plots
scatters = {}

for row_idx, method in enumerate(methods):
    # Wendland (left column)
    ax_w = axes[row_idx, 0]
    ax_w.set_ylabel(method_labels[row_idx], fontsize=11, fontweight='bold', rotation=0, 
                    ha='right', va='center', labelpad=40)
    ax_w.set_aspect('equal')
    ax_w.set_xlim(-0.5, 0.5)
    ax_w.set_ylim(-0.5, 0.5)
    ax_w.set_facecolor('black')
    ax_w.grid(True, alpha=0.2, color='white', linewidth=0.5)
    if row_idx == 4:
        ax_w.set_xlabel('x', fontsize=10)
    else:
        ax_w.set_xticklabels([])
    
    scatters[('wendland', method)] = ax_w.scatter([], [], c=[], cmap='hot', s=15, 
                                                  alpha=0.8, edgecolors='none',
                                                  vmin=vmin, vmax=vmax)
    
    # Cubic Spline (right column)
    ax_c = axes[row_idx, 1]
    ax_c.set_aspect('equal')
    ax_c.set_xlim(-0.5, 0.5)
    ax_c.set_ylim(-0.5, 0.5)
    ax_c.set_facecolor('black')
    ax_c.grid(True, alpha=0.2, color='white', linewidth=0.5)
    ax_c.set_yticklabels([])
    if row_idx == 4:
        ax_c.set_xlabel('x', fontsize=10)
    else:
        ax_c.set_xticklabels([])
    
    scatters[('cubic', method)] = ax_c.scatter([], [], c=[], cmap='hot', s=15, 
                                               alpha=0.8, edgecolors='none',
                                               vmin=vmin, vmax=vmax)
    
    # Add colorbar to middle row, right plot
    if row_idx == 2:
        cbar = plt.colorbar(scatters[('cubic', method)], ax=ax_c, pad=0.02)
        cbar.set_label('log₁₀(Density)', fontsize=10)

# Add time/frame info at bottom
fig.text(0.5, 0.01, '', ha='center', va='bottom', fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
time_text = fig.texts[0]

plt.tight_layout(rect=[0, 0.02, 1, 0.99])

def update(frame_idx):
    frame = all_frames[frame_idx]
    
    for method in methods:
        # Wendland kernel
        if method in frame['data']['wendland']:
            data = frame['data']['wendland'][method]
            x_col = 'x' if 'x' in data else 'pos_x'
            y_col = 'y' if 'y' in data else 'pos_y'
            if x_col in data and y_col in data and 'dens' in data:
                dens_array = np.array(data['dens'])
                log_dens = np.log10(dens_array + 1e-10)
                scatters[('wendland', method)].set_offsets(np.c_[data[x_col], data[y_col]])
                scatters[('wendland', method)].set_array(log_dens)
        
        # Cubic Spline kernel
        if method in frame['data']['cubic']:
            data = frame['data']['cubic'][method]
            x_col = 'x' if 'x' in data else 'pos_x'
            y_col = 'y' if 'y' in data else 'pos_y'
            if x_col in data and y_col in data and 'dens' in data:
                dens_array = np.array(data['dens'])
                log_dens = np.log10(dens_array + 1e-10)
                scatters[('cubic', method)].set_offsets(np.c_[data[x_col], data[y_col]])
                scatters[('cubic', method)].set_array(log_dens)
    
    if frame['time']:
        time_text.set_text(f"Time: {frame['time']} | Snapshot: {frame['number']:04d}")
    else:
        time_text.set_text(f"Snapshot: {frame['number']:04d}")
    
    return list(scatters.values()) + [time_text]

# Create animation
anim = FuncAnimation(fig, update, frames=tqdm(range(len(all_frames))) if HAS_TQDM else range(len(all_frames)),
                    interval=100, blit=False, repeat=True)

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Save animation
print("Saving animation...")
writer = PillowWriter(fps=10)
anim.save(output_file, writer=writer, dpi=80)

plt.close()

print()
print('=' * 70)
print('Animation Complete!')
print('=' * 70)
print(f'✓ Saved: {output_file}')
print('=' * 70)
