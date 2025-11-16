#!/usr/bin/env python3
"""
Animate GSPH Iterative Riemann Solver Results

Creates animation and plots for the corrected iterative Riemann solver
showing density, velocity, pressure, and internal energy evolution.

Usage:
    python3 animate_iterative.py
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
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

# Data directory
data_dir = "sample/shock_tube/results/gsph_iterative"
output_prefix = "gsph_iterative"

# Output directory (same as data directory)
output_dir = data_dir

print('=' * 80)
print('GSPH Iterative Riemann Solver - Shock Tube Visualization')
print('=' * 80)
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
    exit(1)

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
    ax_dens.plot(x_sorted, dens, 'b-', linewidth=2.5, label='GSPH Iterative')
    ax_dens.set_ylabel('Density ρ', fontsize=13, fontweight='bold')
    ax_dens.grid(True, alpha=0.3, linestyle='--')
    ax_dens.legend(loc='best', fontsize=11)
    ax_dens.set_title('Density Profile', fontweight='bold', fontsize=14)
    ax_dens.set_ylim([0, 1.2])
    
    # Velocity
    ax_vel.clear()
    vel = data['vel_x'][sort_idx]
    ax_vel.plot(x_sorted, vel, 'r-', linewidth=2.5, label='GSPH Iterative')
    ax_vel.set_ylabel('Velocity u', fontsize=13, fontweight='bold')
    ax_vel.grid(True, alpha=0.3, linestyle='--')
    ax_vel.legend(loc='best', fontsize=11)
    ax_vel.set_title('Velocity Profile', fontweight='bold', fontsize=14)
    ax_vel.set_ylim([-0.1, 1.0])
    
    # Pressure
    ax_pres.clear()
    pres = data['pres'][sort_idx]
    ax_pres.plot(x_sorted, pres, 'g-', linewidth=2.5, label='GSPH Iterative')
    ax_pres.set_ylabel('Pressure P', fontsize=13, fontweight='bold')
    ax_pres.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax_pres.grid(True, alpha=0.3, linestyle='--')
    ax_pres.legend(loc='best', fontsize=11)
    ax_pres.set_title('Pressure Profile', fontweight='bold', fontsize=14)
    ax_pres.set_ylim([0, 1.2])
    
    # Internal Energy
    ax_ene.clear()
    ene = data['ene'][sort_idx]
    ax_ene.plot(x_sorted, ene, 'm-', linewidth=2.5, label='GSPH Iterative')
    ax_ene.set_ylabel('Internal Energy e', fontsize=13, fontweight='bold')
    ax_ene.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax_ene.grid(True, alpha=0.3, linestyle='--')
    ax_ene.legend(loc='best', fontsize=11)
    ax_ene.set_title('Internal Energy Profile', fontweight='bold', fontsize=14)
    ax_ene.set_ylim([0, 3.0])
    
    # Add time label to top plot
    time_text = ax_dens.text(0.02, 0.95, time_label, transform=ax_dens.transAxes,
                             fontsize=15, fontweight='bold',
                             verticalalignment='top',
                             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9))

# Create final snapshot plot
print('Creating final snapshot plot...')
final_data = load_snapshot(files[-1])
final_time = final_data['metadata'].get('time', 'unknown')

fig_final, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 13))
fig_final.suptitle(f'GSPH Iterative Solver - Shock Tube Final State (t={final_time})', 
                   fontsize=18, fontweight='bold')

plot_shock_tube_frame(ax1, ax2, ax3, ax4, final_data, f't = {final_time}')

plt.tight_layout()
output_final = os.path.join(output_dir, f'{output_prefix}_final.png')
plt.savefig(output_final, dpi=200, bbox_inches='tight')
print(f'✓ Saved: {output_final}')
plt.close()

# Create animation
if HAS_ANIMATION and len(files) > 1:
    print()
    print('Creating animation...')
    
    fig_anim, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 13))
    fig_anim.suptitle('GSPH Iterative Solver - Shock Tube Time Evolution', 
                      fontsize=18, fontweight='bold')
    
    # Use all available frames
    n_frames = len(files)
    frame_indices = list(range(n_frames))
    
    print(f'Total snapshots: {n_frames}')
    print(f'Animation frames: {len(frame_indices)}')
    
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
                        interval=200, blit=False, repeat=True)
    
    output_anim = os.path.join(output_dir, f'{output_prefix}_animation.gif')
    writer = PillowWriter(fps=5)
    anim.save(output_anim, writer=writer)
    
    if pbar:
        pbar.close()
    
    print(f'✓ Saved: {output_anim}')
    plt.close()
else:
    print()
    if not HAS_ANIMATION:
        print('Skipping animation (matplotlib.animation not available)')
    else:
        print('Skipping animation (only 1 snapshot available)')

# Create multi-time comparison plot
if len(files) > 1:
    print()
    print('Creating multi-time comparison plot...')
    
    fig_compare, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 13))
    fig_compare.suptitle('GSPH Iterative Solver - Time Evolution Comparison', 
                         fontsize=18, fontweight='bold')
    
    # Plot snapshots for comparison
    if len(files) <= 6:
        # Use all snapshots
        compare_indices = list(range(len(files)))
    else:
        # Select evenly spaced snapshots
        compare_indices = [0, len(files)//4, len(files)//2, 3*len(files)//4, -1]
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(compare_indices)))
    alphas = np.linspace(0.5, 1.0, len(compare_indices))
    
    for i, (idx, color, alpha) in enumerate(zip(compare_indices, colors, alphas)):
        data = load_snapshot(files[idx])
        time = data['metadata'].get('time', 'unknown')
        
        x = data['pos_x']
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        
        label = f't={time}'
        
        ax1.plot(x_sorted, data['dens'][sort_idx], color=color, alpha=alpha, 
                 linewidth=2.5, label=label)
        ax2.plot(x_sorted, data['vel_x'][sort_idx], color=color, alpha=alpha,
                 linewidth=2.5, label=label)
        ax3.plot(x_sorted, data['pres'][sort_idx], color=color, alpha=alpha,
                 linewidth=2.5, label=label)
        ax4.plot(x_sorted, data['ene'][sort_idx], color=color, alpha=alpha,
                 linewidth=2.5, label=label)
    
    ax1.set_ylabel('Density ρ', fontsize=13, fontweight='bold')
    ax1.set_title('Density Evolution', fontweight='bold', fontsize=14)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.legend(loc='best', fontsize=10)
    ax1.set_ylim([0, 1.2])
    
    ax2.set_ylabel('Velocity u', fontsize=13, fontweight='bold')
    ax2.set_title('Velocity Evolution', fontweight='bold', fontsize=14)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.legend(loc='best', fontsize=10)
    ax2.set_ylim([-0.1, 1.0])
    
    ax3.set_ylabel('Pressure P', fontsize=13, fontweight='bold')
    ax3.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax3.set_title('Pressure Evolution', fontweight='bold', fontsize=14)
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.legend(loc='best', fontsize=10)
    ax3.set_ylim([0, 1.2])
    
    ax4.set_ylabel('Internal Energy e', fontsize=13, fontweight='bold')
    ax4.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax4.set_title('Energy Evolution', fontweight='bold', fontsize=14)
    ax4.grid(True, alpha=0.3, linestyle='--')
    ax4.legend(loc='best', fontsize=10)
    ax4.set_ylim([0, 3.0])
    
    plt.tight_layout()
    output_compare = os.path.join(output_dir, f'{output_prefix}_comparison.png')
    plt.savefig(output_compare, dpi=200, bbox_inches='tight')
    print(f'✓ Saved: {output_compare}')
    plt.close()

print()
print('=' * 80)
print('✓ Visualization Complete!')
print('=' * 80)
print('Generated files:')
print(f'  📊 {os.path.basename(output_final)}        - Final state plot')
if len(files) > 1:
    print(f'  📊 {os.path.basename(output_compare)}   - Time comparison plot')
if HAS_ANIMATION and len(files) > 1:
    print(f'  🎬 {os.path.basename(output_anim)}  - Animation (GIF)')
print(f'\nAll files saved to: {output_dir}')
print('=' * 80)
