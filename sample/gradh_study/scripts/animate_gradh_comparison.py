#!/usr/bin/env python3
"""
Create Side-by-Side Animation: With vs Without Grad-H Correction

Generates an animated GIF showing the evolution of hydrostatic equilibrium
with and without the grad-h correction term.

Usage:
    python3 animate_gradh_comparison.py --with-gradh <dir1> --no-gradh <dir2> --output <output.gif>
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import sys
import os
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Animate grad-h comparison')
parser.add_argument('--with-gradh', required=True, help='Directory with grad-h results')
parser.add_argument('--no-gradh', required=True, help='Directory without grad-h results')
parser.add_argument('--output', required=True, help='Output GIF file')
parser.add_argument('--fps', type=int, default=10, help='Frames per second')
parser.add_argument('--skip', type=int, default=1, help='Skip every N snapshots')
args = parser.parse_args()

# Configuration
with_gradh_dir = args.with_gradh
no_gradh_dir = args.no_gradh
output_file = args.output

cases = ['with_gradh', 'no_gradh']
case_dirs = [with_gradh_dir, no_gradh_dir]
case_labels = ['With Grad-H (Stable)', 'Without Grad-H (Collapse)']

print('=' * 70)
print('Grad-H Comparison Animation Generator')
print('=' * 70)
print(f'With grad-h:    {with_gradh_dir}')
print(f'Without grad-h: {no_gradh_dir}')
print(f'Output:         {output_file}')
print()


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
        lines = [ln for ln in f.readlines() if not ln.startswith('#')]
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


# Scan for snapshots
print("Scanning for snapshot files...")
case_snapshots = {}

for case, case_dir in zip(cases, case_dirs):
    if not os.path.exists(case_dir):
        print(f"ERROR: Directory not found: {case_dir}")
        sys.exit(1)
    
    snapshots = sorted(glob.glob(os.path.join(case_dir, "snapshot_*.csv")))
    if not snapshots:
        print(f"ERROR: No snapshots found in {case_dir}")
        sys.exit(1)
    
    snapshot_numbers = []
    for snap in snapshots:
        basename = os.path.basename(snap)
        num_str = basename.replace('snapshot_', '').replace('.csv', '')
        try:
            snapshot_numbers.append(int(num_str))
        except ValueError:
            pass
    
    case_snapshots[case] = {
        'files': snapshots,
        'numbers': snapshot_numbers
    }
    print(f"  {case:15s}: {len(snapshot_numbers)} snapshots")

# Find common snapshots
common_numbers = set(case_snapshots['with_gradh']['numbers']) & set(case_snapshots['no_gradh']['numbers'])
common_numbers = sorted(common_numbers)[::args.skip]
print(f"\nUsing {len(common_numbers)} frames (skip={args.skip})")

# Set up figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Grad-H Correction Effect on Hydrostatic Equilibrium', fontsize=14, fontweight='bold')

# Determine global density range
print("Computing density range...")
all_dens = []
for snap_num in common_numbers[:min(10, len(common_numbers))]:
    for case_dir in case_dirs:
        snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, _ = load_snapshot(snap_file)
            if 'dens' in data:
                all_dens.extend(data['dens'])

vmin = np.percentile(all_dens, 1)
vmax = np.percentile(all_dens, 99)
print(f"  Density range: [{vmin:.4f}, {vmax:.4f}]")

# Initialize scatter plots
scatters = []
titles = []
for ax, label in zip(axes, case_labels):
    ax.set_xlim(-0.55, 0.55)
    ax.set_ylim(-0.55, 0.55)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('y', fontsize=11)
    title = ax.set_title(f'{label}\nt = 0.00', fontsize=12)
    titles.append(title)
    
    scatter = ax.scatter([], [], c=[], s=2, cmap='viridis', vmin=vmin, vmax=vmax)
    scatters.append(scatter)

# Add colorbar
cbar = fig.colorbar(scatters[0], ax=axes.ravel().tolist(), label='Density', shrink=0.8)

plt.tight_layout()


def init():
    for scatter in scatters:
        scatter.set_offsets(np.empty((0, 2)))
    return scatters + titles


def animate(frame_idx):
    snap_num = common_numbers[frame_idx]
    
    for i, (case_dir, label) in enumerate(zip(case_dirs, case_labels)):
        snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
        
        if os.path.exists(snap_file):
            data, meta = load_snapshot(snap_file)
            
            x = data.get('pos_x', [])
            y = data.get('pos_y', [])
            dens = data.get('dens', np.ones_like(x))
            
            # Update scatter plot
            offsets = np.column_stack([x, y])
            scatters[i].set_offsets(offsets)
            scatters[i].set_array(dens)
            
            # Update title
            time = float(meta.get('time', 0))
            titles[i].set_text(f'{label}\nt = {time:.2f}')
    
    return scatters + titles


print(f"\nGenerating animation with {len(common_numbers)} frames...")

# Create animation
anim = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=len(common_numbers),
    interval=1000/args.fps,
    blit=True
)

# Save animation
print(f"Saving to {output_file}...")
anim.save(output_file, writer='pillow', fps=args.fps, dpi=100)
print("✓ Animation saved!")

plt.close()
print('=' * 70)
print('Done!')
print('=' * 70)
