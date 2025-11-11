#!/usr/bin/env python3
"""
Compare Pairing Instability Results from Multiple SPH Methods

Creates side-by-side comparison plots showing GSPH, SSPH, DISPH, GDISPH, and GDISPH+Balsara results.

Usage:
    python3 compare_pairing.py [base_dir] [output_dir]
    
Arguments:
    base_dir    - Base directory containing gsph/, ssph/, disph/, gdisph/, gdisph_balsara/ subdirs
                  Default: sample/pairing_instability/results
    output_dir  - Where to save comparison plots
                  Default: sample/pairing_instability/results/comparison
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
from pathlib import Path

# Configuration
base_dir = sys.argv[1] if len(sys.argv) > 1 else "sample/pairing_instability/results"
output_dir = sys.argv[2] if len(sys.argv) > 2 else "sample/pairing_instability/results/comparison"

# SPH methods to compare
methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
method_labels = ['GSPH (Godunov)', 'SSPH (Standard)', 'DISPH (Density-Independent)', 
                 'GDISPH (Godunov DISPH)', 'GDISPH+Balsara']
# Colorblind-friendly palette with high contrast
method_colors = ['#0173B2', '#DE8F05', '#029E73', '#CC78BC', '#D55E00']  # Blue, orange, teal, purple, red-orange

print('=' * 70)
print('Pairing Instability Multi-Method Comparison (5 Methods)')
print('=' * 70)
print(f'Base directory:   {base_dir}')
print(f'Output directory: {output_dir}')
print(f'Methods:          {", ".join(method_labels)}')
print()

# Create output directory
os.makedirs(output_dir, exist_ok=True)

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
    
    return data, metadata

# Scan for snapshot files
print("Scanning for snapshot files...")
method_snapshots = {}
for method in methods:
    method_dir = os.path.join(base_dir, method)
    if not os.path.exists(method_dir):
        print(f"  WARNING: {method.upper()} directory not found: {method_dir}")
        continue
    
    snapshots = sorted(glob.glob(os.path.join(method_dir, "snapshot_*.csv")))
    if not snapshots:
        print(f"  WARNING: No snapshots found for {method.upper()}")
        continue
    
    # Extract snapshot numbers
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

# Find common snapshot numbers
if not method_snapshots:
    print("ERROR: No snapshot data found!")
    sys.exit(1)

common_numbers = set(method_snapshots[methods[0]]['numbers'])
for method in methods[1:]:
    if method in method_snapshots:
        common_numbers &= set(method_snapshots[method]['numbers'])

common_numbers = sorted(common_numbers)
print(f"\nFound {len(common_numbers)} common snapshot numbers across all methods")

# Select snapshots to plot (initial, 25%, 50%, 75%, final)
if len(common_numbers) >= 5:
    indices = [0, len(common_numbers)//4, len(common_numbers)//2, 3*len(common_numbers)//4, -1]
    plot_snapshots = [common_numbers[i] for i in indices]
else:
    plot_snapshots = common_numbers

print(f"Generating comparison plots for {len(plot_snapshots)} snapshots...\n")

# Generate comparison plots
for snap_num in plot_snapshots:
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Pairing Instability - 5 Method Comparison (Snapshot {snap_num:04d})', 
                 fontsize=16, fontweight='bold')
    
    # Load data for all methods
    method_data = {}
    current_time = None
    
    for method in methods:
        if method not in method_snapshots:
            continue
        
        snap_file = os.path.join(base_dir, method, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, metadata = load_snapshot(snap_file)
            method_data[method] = data
            if current_time is None and 'time' in metadata:
                current_time = metadata['time']
    
    # Plot density fields (2D scatter plots)
    for idx, method in enumerate(methods):
        if method not in method_data:
            continue
        
        ax = axes.flatten()[idx]
        data = method_data[method]
        
        # Create 2D density plot
        # Handle both 'x','y' and 'pos_x','pos_y' column naming
        x_col = 'x' if 'x' in data else 'pos_x'
        y_col = 'y' if 'y' in data else 'pos_y'
        
        if x_col in data and y_col in data and 'dens' in data:
            # Use log scale for better dynamic range and 'hot' colormap for instabilities
            import numpy as np
            dens_array = np.array(data['dens'])
            log_dens = np.log10(dens_array + 1e-10)  # Add small value to avoid log(0)
            
            scatter = ax.scatter(data[x_col], data[y_col], c=log_dens, 
                               cmap='hot', s=25, alpha=0.8, edgecolors='none')
            ax.set_xlabel('x', fontsize=10)
            ax.set_ylabel('y', fontsize=10)
            ax.set_title(method_labels[idx], fontsize=12, fontweight='bold')
            ax.set_aspect('equal')
            ax.set_facecolor('black')  # Black background emphasizes hot spots
            ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
            
            # Add colorbar with log scale label
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('log₁₀(Density)', fontsize=9)
    
    # Use last subplot for legend/info
    ax = axes.flatten()[5]
    ax.axis('off')
    
    # Add time information
    if current_time:
        ax.text(0.5, 0.7, f'Time: {current_time}', 
               ha='center', va='center', fontsize=14, fontweight='bold',
               transform=ax.transAxes)
    
    # Add method legend
    legend_text = '\n'.join([f'{label}' for label in method_labels])
    ax.text(0.5, 0.3, legend_text,
           ha='center', va='center', fontsize=10,
           transform=ax.transAxes,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, f'comparison_snap{snap_num:04d}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'✓ Saved: {os.path.basename(output_file)}')

# Create final state comparison
final_snap = common_numbers[-1]
print(f"\nCreating final state comparison (snapshot {final_snap:04d})...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Pairing Instability - Final State Comparison', 
             fontsize=16, fontweight='bold')

for idx, method in enumerate(methods):
    if method not in method_snapshots:
        continue
    
    snap_file = os.path.join(base_dir, method, f"snapshot_{final_snap:04d}.csv")
    if os.path.exists(snap_file):
        data, metadata = load_snapshot(snap_file)
        
        ax = axes.flatten()[idx]
        
        # Handle both 'x','y' and 'pos_x','pos_y' column naming
        x_col = 'x' if 'x' in data else 'pos_x'
        y_col = 'y' if 'y' in data else 'pos_y'
        
        if x_col in data and y_col in data and 'dens' in data:
            # Use log scale for better dynamic range and 'hot' colormap for instabilities
            import numpy as np
            dens_array = np.array(data['dens'])
            log_dens = np.log10(dens_array + 1e-10)
            
            scatter = ax.scatter(data[x_col], data[y_col], c=log_dens, 
                               cmap='hot', s=25, alpha=0.8, edgecolors='none',
                               vmin=np.percentile(log_dens, 1), 
                               vmax=np.percentile(log_dens, 99))
            ax.set_xlabel('x', fontsize=10)
            ax.set_ylabel('y', fontsize=10)
            ax.set_title(method_labels[idx], fontsize=12, fontweight='bold')
            ax.set_aspect('equal')
            ax.set_facecolor('black')
            ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
            
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('log₁₀(Density)', fontsize=9)

ax = axes.flatten()[5]
ax.axis('off')
if 'time' in metadata:
    ax.text(0.5, 0.5, f'Final Time: {metadata["time"]}', 
           ha='center', va='center', fontsize=14, fontweight='bold',
           transform=ax.transAxes)

plt.tight_layout()
output_file = os.path.join(output_dir, 'comparison_final.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()
print(f'✓ Saved: {os.path.basename(output_file)}')

print()
print('=' * 70)
print('Comparison Complete!')
print('=' * 70)
print(f'Output saved to: {output_dir}')
print('=' * 70)
