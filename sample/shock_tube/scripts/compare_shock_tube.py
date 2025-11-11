#!/usr/bin/env python3
"""
Compare 1D Shock Tube Results from Multiple SPH Methods

Creates side-by-side comparison plots showing GSPH, SSPH, DISPH, GDISPH, and GDISPH+Balsara results.

Usage:
    python3 compare_shock_tube.py [base_dir] [output_dir]
    
Arguments:
    base_dir    - Base directory containing gsph/, ssph/, disph/, gdisph/, gdisph_balsara/ subdirs
                  Default: sample/shock_tube/results
    output_dir  - Where to save comparison plots
                  Default: sample/shock_tube/results/comparison
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
from pathlib import Path

# Configuration
base_dir = sys.argv[1] if len(sys.argv) > 1 else "sample/shock_tube/results"
output_dir = sys.argv[2] if len(sys.argv) > 2 else "sample/shock_tube/results/comparison"

# SPH methods to compare
methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
method_labels = ['GSPH (Godunov)', 'SSPH (Standard)', 'DISPH (Density-Independent)', 
                 'GDISPH (Godunov DISPH)', 'GDISPH+Balsara']
# Colorblind-friendly palette with high contrast
method_colors = ['#0173B2', '#DE8F05', '#029E73', '#CC78BC', '#D55E00']  # Blue, orange, teal, purple, red-orange
method_styles = ['-', '--', '-.', ':', '-']  # Solid, dashed, dash-dot, dotted, solid
method_markers = ['o', 's', '^', 'D', 'v']  # Circle, square, triangle-up, diamond, triangle-down
method_markevery = 5  # Show marker every N points

print('=' * 70)
print('Shock Tube Multi-Method Comparison (5 Methods)')
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

def plot_comparison_frame(data_dict, time_label, output_file):
    """
    Plot comparison of all methods for a single time snapshot
    
    Args:
        data_dict: Dictionary {method: data} for each SPH method
        time_label: Time label string
        output_file: Output filename
    """
    # Get unit system from first available data
    unit_system = 'code_units'
    for data in data_dict.values():
        if 'unit_system' in data:
            unit_system = data['unit_system']
            break
    
    units = get_unit_labels(unit_system)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
    
    # Create comprehensive title with timestamp and unit system
    fig.suptitle(f'Shock Tube Comparison: GSPH vs SSPH vs DISPH\nTime: {time_label} {units["time"]} | Unit System: {unit_system}', 
                 fontsize=18, fontweight='bold', y=0.98)
    
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
    
    # Format density plot
    ax1.set_ylabel(f'Density ρ [{units["density"]}]', fontsize=13, fontweight='bold')
    ax1.set_title('Density Profile', fontweight='bold', fontsize=14)
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax1.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    ax1.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    
    # Format velocity plot
    ax2.set_ylabel(f'Velocity u [{units["velocity"]}]', fontsize=13, fontweight='bold')
    ax2.set_title('Velocity Profile', fontweight='bold', fontsize=14)
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax2.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    ax2.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    
    # Format pressure plot
    ax3.set_ylabel(f'Pressure P [{units["pressure"]}]', fontsize=13, fontweight='bold')
    ax3.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    ax3.set_title('Pressure Profile', fontweight='bold', fontsize=14)
    ax3.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax3.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    
    # Format energy plot
    ax4.set_ylabel(f'Internal Energy e [{units["energy"]}]', fontsize=13, fontweight='bold')
    ax4.set_xlabel(f'Position x [{units["length"]}]', fontsize=13, fontweight='bold')
    ax4.set_title('Internal Energy Profile', fontweight='bold', fontsize=14)
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax4.legend(loc='best', fontsize=12, framealpha=0.95, edgecolor='gray')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    return output_file

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

if all_snapshot_numbers:
    common_snapshots = sorted(list(set.intersection(*all_snapshot_numbers)))
    print(f'Found {len(common_snapshots)} common snapshot numbers across all methods')
else:
    print('WARNING: No common snapshots found, using available snapshots')
    common_snapshots = []

# Generate comparison plots for key snapshots
if common_snapshots:
    # Select snapshots to compare: first, middle, final, and a few in between
    n_common = len(common_snapshots)
    compare_indices = [0]
    if n_common > 4:
        compare_indices.extend([n_common//4, n_common//2, 3*n_common//4])
    if n_common > 1:
        compare_indices.append(n_common - 1)
    
    compare_snapshot_nums = [common_snapshots[i] for i in compare_indices]
    
    print(f'Generating comparison plots for {len(compare_snapshot_nums)} snapshots...')
    print()
    
    for snap_num in compare_snapshot_nums:
        # Load data from each method
        data_dict = {}
        time_label = None
        
        for method in methods:
            # Find file with this snapshot number
            for filename in method_snapshots[method]:
                if extract_snapshot_number(filename) == snap_num:
                    data = load_snapshot(filename)
                    data_dict[method] = data
                    if time_label is None:
                        time_label = data['metadata'].get('time', f'snap{snap_num:04d}')
                    break
        
        if data_dict:
            output_file = os.path.join(output_dir, f'comparison_t{time_label}.png')
            plot_comparison_frame(data_dict, time_label, output_file)
            print(f'✓ Saved: {os.path.basename(output_file)}')

    # Create final comparison plot
    print()
    print('Creating final state comparison...')
    final_snap_num = common_snapshots[-1]
    data_dict = {}
    time_label = None
    
    for method in methods:
        for filename in method_snapshots[method]:
            if extract_snapshot_number(filename) == final_snap_num:
                data = load_snapshot(filename)
                data_dict[method] = data
                if time_label is None:
                    time_label = data['metadata'].get('time', f'snap{final_snap_num:04d}')
                break
    
    if data_dict:
        output_file = os.path.join(output_dir, 'comparison_final.png')
        plot_comparison_frame(data_dict, time_label, output_file)
        print(f'✓ Saved: {os.path.basename(output_file)}')

else:
    print('WARNING: Could not find common snapshots to compare')

print()
print('=' * 70)
print('Comparison Complete!')
print('=' * 70)
print(f'Output saved to: {output_dir}')
print('=' * 70)
