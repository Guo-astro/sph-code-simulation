#!/usr/bin/env python3
"""
Create Side-by-Side Animation Comparing GSPH, SSPH, DISPH, and GDISPH

Generates an animated GIF showing time evolution of all four SPH methods simultaneously.

Usage:
    python3 animate_comparison.py [base_dir] [output_file]
    
Arguments:
    base_dir    - Base directory containing gsph/, ssph/, disph/, gdisph/ subdirs
                  Default: sample/shock_tube/results
    output_file - Output animation filename
                  Default: sample/shock_tube/results/comparison/comparison_animation.gif
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
base_dir = sys.argv[1] if len(sys.argv) > 1 else "sample/shock_tube/results"
output_file = sys.argv[2] if len(sys.argv) > 2 else "sample/shock_tube/results/comparison/comparison_animation.gif"

# SPH methods to compare
methods = ['gsph', 'ssph', 'disph', 'gdisph']
method_labels = ['GSPH (Godunov)', 'SSPH (Standard)', 'DISPH (Density-Independent)', 'GDISPH (Godunov DISPH)']
# Colorblind-friendly palette with high contrast
method_colors = ['#0173B2', '#DE8F05', '#029E73', '#CC78BC']  # Blue, orange, teal, purple
method_styles = ['-', '--', '-.', ':']  # Solid, dashed, dash-dot, dotted
method_markers = ['o', 's', '^', 'D']  # Circle, square, triangle, diamond
method_markevery = 5  # Show marker every N points

print('=' * 70)
print('Shock Tube Animation - Multi-Method Comparison')
print('=' * 70)
print(f'Base directory: {base_dir}')
print(f'Output file:    {output_file}')
print(f'Methods:        {", ".join(method_labels)}')
print()

def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    with open(filename, 'r') as f:
        # Read metadata lines (start with #)
        metadata = {}
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        # Read header and data
        f.seek(0)
        lines = [l for l in f.readlines() if not l.startswith('#')]
        header = lines[0].strip().split(',')
        
        # Parse data
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
    
    # Parse unit system from metadata
    unit_system = metadata.get('unit_system', 'code_units')
    data['unit_system'] = unit_system
    
    return data

def find_snapshots_by_method(base_dir, method):
    """Find all snapshot files for a given method"""
    method_dir = os.path.join(base_dir, method)
    files = sorted(glob.glob(f'{method_dir}/snapshot_*.csv'))
    return files

def extract_snapshot_number(filename):
    """Extract snapshot number from filename"""
    import re
    match = re.search(r'snapshot_(\d+)\.csv', filename)
    return int(match.group(1)) if match else -1

# Collect snapshots for each method
print('Scanning for snapshot files...')
method_snapshots = {}
for method in methods:
    files = find_snapshots_by_method(base_dir, method)
    method_snapshots[method] = files
    print(f'  {method.upper():5s}: {len(files)} snapshots')

if all(len(files) == 0 for files in method_snapshots.values()):
    print('\nERROR: No snapshot files found for any method!')
    print(f'Searched in subdirectories: {", ".join(methods)}')
    sys.exit(1)

print()

# Find common snapshot numbers across all methods
all_snapshot_numbers = []
for method in methods:
    if len(method_snapshots[method]) > 0:
        numbers = [extract_snapshot_number(f) for f in method_snapshots[method]]
        all_snapshot_numbers.append(set(numbers))

if not all_snapshot_numbers:
    print('ERROR: No snapshots found')
    sys.exit(1)

common_snapshots = sorted(list(set.intersection(*all_snapshot_numbers)))

if not common_snapshots:
    print('ERROR: No common snapshot numbers found across all methods')
    print('Each method must have corresponding snapshots (same snapshot numbers)')
    sys.exit(1)

print(f'Found {len(common_snapshots)} common snapshots')

# Determine frame skip for reasonable file size
n_frames = len(common_snapshots)
max_frames = 50  # Limit animation frames
frame_skip = max(1, n_frames // max_frames)
animation_snapshot_nums = common_snapshots[::frame_skip]

print(f'Animation frames: {len(animation_snapshot_nums)} (every {frame_skip} snapshots)')
print()

# Pre-load all data for animation
print('Pre-loading data...')
animation_data = []

pbar_load = tqdm(total=len(animation_snapshot_nums), desc='Loading snapshots') if HAS_TQDM else None

for snap_num in animation_snapshot_nums:
    frame_data = {}
    time_label = None
    
    for method in methods:
        # Find file with this snapshot number
        for filename in method_snapshots[method]:
            if extract_snapshot_number(filename) == snap_num:
                data = load_snapshot(filename)
                frame_data[method] = data
                if time_label is None:
                    time_label = data['metadata'].get('time', f'{snap_num:04d}')
                break
    
    animation_data.append((frame_data, time_label))
    
    if pbar_load:
        pbar_load.update(1)

if pbar_load:
    pbar_load.close()

print(f'Loaded {len(animation_data)} frames')
print()

# Get unit system from first frame
unit_system = 'code_units'
for frame_data, _ in animation_data:
    for data in frame_data.values():
        if 'unit_system' in data:
            unit_system = data['unit_system']
            break
    if unit_system != 'code_units':
        break

def get_unit_labels(unit_system):
    """Get proper unit labels based on unit system"""
    if unit_system == 'solar':
        return {
            'length': 'R☉',
            'velocity': 'km/s',
            'density': 'g/cm³',
            'pressure': 'dyne/cm²',
            'energy': 'erg/g',
            'time': 's'
        }
    elif unit_system == 'cgs':
        return {
            'length': 'cm',
            'velocity': 'cm/s',
            'density': 'g/cm³',
            'pressure': 'dyne/cm²',
            'energy': 'erg/g',
            'time': 's'
        }
    elif unit_system == 'si':
        return {
            'length': 'm',
            'velocity': 'm/s',
            'density': 'kg/m³',
            'pressure': 'Pa',
            'energy': 'J/kg',
            'time': 's'
        }
    else:  # code_units or unknown
        return {
            'length': 'code units',
            'velocity': 'code units',
            'density': 'code units',
            'pressure': 'code units',
            'energy': 'code units',
            'time': 'code units'
        }

units = get_unit_labels(unit_system)

# Create animation
print('Generating animation...')

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))

# Create title with unit system info
title_text = fig.suptitle('', fontsize=18, fontweight='bold', y=0.98)

# Initialize plots
for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=1)

ax1.set_ylabel(f'Density ρ [{units["density"]}]', fontsize=13, fontweight='bold')
ax1.set_title('Density Profile', fontweight='bold', fontsize=14)

ax2.set_ylabel(f'Velocity u [{units["velocity"]}]', fontsize=13, fontweight='bold')
ax2.set_title('Velocity Profile', fontweight='bold', fontsize=14)

ax3.set_ylabel(f'Pressure P [{units["pressure"]}]', fontsize=13, fontweight='bold')
ax3.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
ax3.set_title('Pressure Profile', fontweight='bold', fontsize=14)

ax4.set_ylabel(f'Internal Energy e [{units["energy"]}]', fontsize=13, fontweight='bold')
ax4.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
ax4.set_title('Internal Energy Profile', fontweight='bold', fontsize=14)

# Time text annotation
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes,
                     fontsize=14, fontweight='bold',
                     verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

pbar_anim = tqdm(total=len(animation_data), desc='Rendering frames') if HAS_TQDM else None

def update(frame_num):
    """Update function for animation"""
    data_dict, time_label = animation_data[frame_num]
    
    # Update main title with timestamp and unit system
    title_text.set_text(f'Shock Tube Comparison: GSPH vs SSPH vs DISPH\nTime: {time_label} {units["time"]} | Unit System: {unit_system}')
    
    # Clear axes
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    
    # Plot each method
    for idx, (method, label, color, style, marker) in enumerate(zip(methods, method_labels,
                                                                      method_colors, method_styles,
                                                                      method_markers)):
        if method not in data_dict:
            continue
        
        data = data_dict[method]
        x = data['pos_x']
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        
        # Density
        dens = data['dens'][sort_idx]
        ax1.plot(x_sorted, dens, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=method_markevery,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
        
        # Velocity
        vel = data['vel_x'][sort_idx]
        ax2.plot(x_sorted, vel, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=method_markevery,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
        
        # Pressure
        pres = data['pres'][sort_idx]
        ax3.plot(x_sorted, pres, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=method_markevery,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
        
        # Internal Energy
        ene = data['ene'][sort_idx]
        ax4.plot(x_sorted, ene, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=method_markevery,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
    
    # Re-apply formatting
    ax1.set_ylabel(f'Density ρ [{units["density"]}]', fontsize=13, fontweight='bold')
    ax1.set_title('Density Profile', fontweight='bold', fontsize=14)
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax1.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    ax1.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    
    ax2.set_ylabel(f'Velocity u [{units["velocity"]}]', fontsize=13, fontweight='bold')
    ax2.set_title('Velocity Profile', fontweight='bold', fontsize=14)
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax2.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    ax2.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    
    ax3.set_ylabel(f'Pressure P [{units["pressure"]}]', fontsize=13, fontweight='bold')
    ax3.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    ax3.set_title('Pressure Profile', fontweight='bold', fontsize=14)
    ax3.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax3.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    
    ax4.set_ylabel(f'Internal Energy e [{units["energy"]}]', fontsize=13, fontweight='bold')
    ax4.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    ax4.set_title('Internal Energy Profile', fontweight='bold', fontsize=14)
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax4.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    
    # Update time label (optional since it's in main title now, but keep for emphasis)
    time_text.set_text(f't = {time_label} {units["time"]}')
    
    if pbar_anim:
        pbar_anim.update(1)
    
    return ax1, ax2, ax3, ax4, time_text, title_text

anim = FuncAnimation(fig, update, frames=len(animation_data),
                    interval=100, blit=False, repeat=True)

# Save animation
os.makedirs(os.path.dirname(output_file), exist_ok=True)
writer = PillowWriter(fps=10)
anim.save(output_file, writer=writer, dpi=200)

if pbar_anim:
    pbar_anim.close()

plt.close()

print()
print('=' * 70)
print('Animation Complete!')
print('=' * 70)
print(f'✓ Saved: {output_file}')
print('=' * 70)
