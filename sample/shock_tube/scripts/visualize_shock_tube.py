#!/usr/bin/env python3
"""
Visualize 1D Shock Tube (Sod Test) Results

Creates plots and animations showing:
- Density profile
- Velocity profile  
- Pressure profile
- Internal energy profile
- Comparison with analytical solution (if available)

Usage:
    python3 visualize_shock_tube.py [data_dir] [output_prefix]
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
    print("Warning: Animation tools not available")

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# Get data directory and output prefix from command line or use defaults
data_dir = sys.argv[1] if len(sys.argv) > 1 else "sample/shock_tube/results"
output_prefix = sys.argv[2] if len(sys.argv) > 2 else "shock_tube"

# Output directory (same as data directory)
output_dir = data_dir

print('=' * 60)
print('Shock Tube Visualization')
print('=' * 60)
print(f'Data directory:   {data_dir}')
print(f'Output directory: {output_dir}')
print(f'Output prefix:    {output_prefix}')
print()

# Find snapshot files
print('Loading data files...')
files = sorted(glob.glob(f'{data_dir}/snapshot_*.csv'))

if len(files) == 0:
    print('ERROR: No snapshot files found!')
    print(f'Searched in: {data_dir}')
    print('Looking for: snapshot_*.csv')
    sys.exit(1)

print(f'Found {len(files)} CSV snapshot files')
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
    return data

def plot_shock_tube_frame(ax_dens, ax_vel, ax_pres, ax_ene, data, time_label):
    """Plot a single frame of shock tube data"""
    
    # Extract position (x-coordinate)
    x = data['pos_x']
    
    # Sort by position for clean line plots
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    
    # Density
    ax_dens.clear()
    dens = data['dens'][sort_idx]
    ax_dens.plot(x_sorted, dens, 'b-', linewidth=2, label='SPH')
    ax_dens.set_ylabel('Density ρ', fontsize=12, fontweight='bold')
    ax_dens.grid(True, alpha=0.3)
    ax_dens.legend(loc='best')
    ax_dens.set_title('Density Profile', fontweight='bold')
    
    # Velocity
    ax_vel.clear()
    vel = data['vel_x'][sort_idx]
    ax_vel.plot(x_sorted, vel, 'r-', linewidth=2, label='SPH')
    ax_vel.set_ylabel('Velocity u', fontsize=12, fontweight='bold')
    ax_vel.grid(True, alpha=0.3)
    ax_vel.legend(loc='best')
    ax_vel.set_title('Velocity Profile', fontweight='bold')
    
    # Pressure
    ax_pres.clear()
    pres = data['pres'][sort_idx]
    ax_pres.plot(x_sorted, pres, 'g-', linewidth=2, label='SPH')
    ax_pres.set_ylabel('Pressure P', fontsize=12, fontweight='bold')
    ax_pres.set_xlabel('Position x', fontsize=12, fontweight='bold')
    ax_pres.grid(True, alpha=0.3)
    ax_pres.legend(loc='best')
    ax_pres.set_title('Pressure Profile', fontweight='bold')
    
    # Internal Energy
    ax_ene.clear()
    ene = data['ene'][sort_idx]
    ax_ene.plot(x_sorted, ene, 'm-', linewidth=2, label='SPH')
    ax_ene.set_ylabel('Internal Energy e', fontsize=12, fontweight='bold')
    ax_ene.set_xlabel('Position x', fontsize=12, fontweight='bold')
    ax_ene.grid(True, alpha=0.3)
    ax_ene.legend(loc='best')
    ax_ene.set_title('Internal Energy Profile', fontweight='bold')
    
    # Add time label to top plot
    time_text = ax_dens.text(0.02, 0.95, time_label, transform=ax_dens.transAxes,
                             fontsize=14, fontweight='bold',
                             verticalalignment='top',
                             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Create final snapshot plot
print('Creating final snapshot plot...')
final_data = load_snapshot(files[-1])
final_time = final_data['metadata'].get('time', 'unknown')

fig_final, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
fig_final.suptitle(f'Shock Tube Test - Final State (t={final_time})', 
                   fontsize=16, fontweight='bold')

plot_shock_tube_frame(ax1, ax2, ax3, ax4, final_data, f't = {final_time}')

plt.tight_layout()
output_final = os.path.join(output_dir, f'{output_prefix}_final.png')
plt.savefig(output_final, dpi=150, bbox_inches='tight')
print(f'✓ Saved: {output_final}')
plt.close()

# Create animation
if HAS_ANIMATION:
    print()
    print('Creating animation...')
    
    fig_anim, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig_anim.suptitle('Shock Tube Test - Time Evolution', 
                      fontsize=16, fontweight='bold')
    
    # Determine frame skip for reasonable file size
    n_frames = len(files)
    max_frames = 50  # Limit animation frames
    frame_skip = max(1, n_frames // max_frames)
    frame_indices = list(range(0, n_frames, frame_skip))
    
    print(f'Total snapshots: {n_frames}')
    print(f'Animation frames: {len(frame_indices)} (every {frame_skip} snapshots)')
    
    pbar = tqdm(total=len(frame_indices), desc='Generating frames') if HAS_TQDM else None
    
    def update(frame_num):
        file_idx = frame_indices[frame_num]
        data = load_snapshot(files[file_idx])
        time_label = f't = {data["metadata"].get("time", "unknown")}'
        
        plot_shock_tube_frame(ax1, ax2, ax3, ax4, data, time_label)
        
        if pbar:
            pbar.update(1)
        
        return ax1, ax2, ax3, ax4
    
    anim = FuncAnimation(fig_anim, update, frames=len(frame_indices),
                        interval=100, blit=False, repeat=True)
    
    output_anim = os.path.join(output_dir, f'{output_prefix}_animation.gif')
    writer = PillowWriter(fps=10)
    anim.save(output_anim, writer=writer)
    
    if pbar:
        pbar.close()
    
    print(f'✓ Saved: {output_anim}')
    plt.close()
else:
    print()
    print('Skipping animation (matplotlib.animation not available)')

# Create multi-time comparison plot
print()
print('Creating multi-time comparison plot...')

fig_compare, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
fig_compare.suptitle('Shock Tube Test - Time Evolution Comparison', 
                     fontsize=16, fontweight='bold')

# Plot every 5th snapshot or fewer for clarity
compare_indices = [0, len(files)//4, len(files)//2, 3*len(files)//4, -1]
colors = ['blue', 'green', 'orange', 'red', 'purple']
alphas = [0.4, 0.5, 0.6, 0.8, 1.0]

for idx, color, alpha in zip(compare_indices, colors, alphas):
    data = load_snapshot(files[idx])
    time = data['metadata'].get('time', 'unknown')
    
    x = data['pos_x']
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    
    label = f't={time}'
    
    ax1.plot(x_sorted, data['dens'][sort_idx], color=color, alpha=alpha, 
             linewidth=2, label=label)
    ax2.plot(x_sorted, data['vel_x'][sort_idx], color=color, alpha=alpha,
             linewidth=2, label=label)
    ax3.plot(x_sorted, data['pres'][sort_idx], color=color, alpha=alpha,
             linewidth=2, label=label)
    ax4.plot(x_sorted, data['ene'][sort_idx], color=color, alpha=alpha,
             linewidth=2, label=label)

ax1.set_ylabel('Density ρ', fontsize=12, fontweight='bold')
ax1.set_title('Density Evolution', fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best')

ax2.set_ylabel('Velocity u', fontsize=12, fontweight='bold')
ax2.set_title('Velocity Evolution', fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(loc='best')

ax3.set_ylabel('Pressure P', fontsize=12, fontweight='bold')
ax3.set_xlabel('Position x', fontsize=12, fontweight='bold')
ax3.set_title('Pressure Evolution', fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend(loc='best')

ax4.set_ylabel('Internal Energy e', fontsize=12, fontweight='bold')
ax4.set_xlabel('Position x', fontsize=12, fontweight='bold')
ax4.set_title('Energy Evolution', fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.legend(loc='best')

plt.tight_layout()
output_compare = os.path.join(output_dir, f'{output_prefix}_comparison.png')
plt.savefig(output_compare, dpi=150, bbox_inches='tight')
print(f'✓ Saved: {output_compare}')

print()
print('=' * 60)
print('Visualization Complete!')
print('=' * 60)
print('Generated files:')
print(f'  📊 {os.path.basename(output_final)}        - Final state')
print(f'  📊 {os.path.basename(output_compare)}   - Time comparison')
if HAS_ANIMATION:
    print(f'  🎬 {os.path.basename(output_anim)}  - Animation')
print(f'All files saved to: {output_dir}')
print('=' * 60)
