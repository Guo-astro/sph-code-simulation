#!/usr/bin/env python3
"""
Compare 1D Strong Shock Results from Multiple SPH Methods

Creates side-by-side comparison plots showing GSPH, SSPH, DISPH, GDISPH, and GDISPH+Balsara results.

Usage:
    python3 compare_methods.py --results-dir RESULTS_DIR [-o OUTPUT] [--no-show]
    
Arguments:
    --results-dir  - Directory containing method subdirectories (gsph_cubic/, ssph_cubic/, etc.)
    -o, --output   - Output filename (default: comparison_strong_shock.png)
    --no-show      - Don't display plot, just save to file
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Compare strong shock results from multiple SPH methods')
    parser.add_argument('--results-dir', type=str, default='sample/strong_shock/results',
                        help='Base results directory')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output filename')
    parser.add_argument('--snapshot', type=int, default=-1,
                        help='Snapshot number to plot (-1 for final)')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display plot')
    return parser.parse_args()

# SPH methods to compare (Cubic Spline kernel only for 1D)
METHODS = {
    'gsph_cubic': ('GSPH (Godunov)', '#0173B2', '-', 'o'),
    'ssph_cubic': ('SSPH (Standard)', '#DE8F05', '--', 's'),
    'disph_cubic': ('DISPH (Density-Independent)', '#029E73', '-.', '^'),
    'gdisph_cubic': ('GDISPH (Godunov DISPH)', '#CC78BC', ':', 'D'),
    'gdisph_balsara_cubic': ('GDISPH+Balsara', '#D55E00', '-', 'v')
}

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

def find_snapshot_files(method_dir, snapshot_num=-1):
    """Find snapshot file for given snapshot number"""
    files = sorted(glob.glob(f'{method_dir}/snapshot_*.csv'))
    if not files:
        return None
    
    if snapshot_num < 0:
        return files[-1]  # Return final snapshot
    
    # Find specific snapshot
    import re
    for f in files:
        match = re.search(r'snapshot_(\d+)\.csv', f)
        if match and int(match.group(1)) == snapshot_num:
            return f
    
    return None

def plot_comparison(data_dict, time_str, output_file=None, show=True):
    """Plot comparison of all methods"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
    
    fig.suptitle(f'Strong Shock Test (P_ratio = 10,000)\nTime: {time_str} | All Methods with Cubic Spline Kernel', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    # Plot each method
    for method, (label, color, style, marker) in METHODS.items():
        if method not in data_dict:
            continue
        
        data = data_dict[method]
        x = data['pos_x']
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        
        # Density
        dens = data['dens'][sort_idx]
        ax1.plot(x_sorted, dens, color=color, linestyle=style, linewidth=2.5, 
                label=label, alpha=0.9, marker=marker, markevery=10, 
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
        
        # Velocity
        vel = data['vel_x'][sort_idx]
        ax2.plot(x_sorted, vel, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=10,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
        
        # Pressure
        pres = data['pres'][sort_idx]
        ax3.plot(x_sorted, pres, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=10,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
        
        # Internal Energy
        ene = data['ene'][sort_idx]
        ax4.plot(x_sorted, ene, color=color, linestyle=style, linewidth=2.5,
                label=label, alpha=0.9, marker=marker, markevery=10,
                markersize=6, markeredgewidth=1.5, markeredgecolor='white')
    
    # Format plots
    ax1.set_ylabel('Density ρ [code units]', fontsize=13, fontweight='bold')
    ax1.set_xlabel('Position x [code units]', fontsize=13, fontweight='bold')
    ax1.set_title('Density Profile', fontweight='bold', fontsize=14)
    ax1.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax1.legend(loc='best', fontsize=11, framealpha=0.95, edgecolor='gray')
    
    ax2.set_ylabel('Velocity u [code units]', fontsize=13, fontweight='bold')
    ax2.set_xlabel('Position x [code units]', fontsize=13, fontweight='bold')
    ax2.set_title('Velocity Profile', fontweight='bold', fontsize=14)
    ax2.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax2.legend(loc='best', fontsize=11, framealpha=0.95, edgecolor='gray')
    
    ax3.set_ylabel('Pressure P [code units]', fontsize=13, fontweight='bold')
    ax3.set_xlabel('Position x [code units]', fontsize=13, fontweight='bold')
    ax3.set_title('Pressure Profile', fontweight='bold', fontsize=14)
    ax3.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax3.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax3.legend(loc='best', fontsize=11, framealpha=0.95, edgecolor='gray')
    
    ax4.set_ylabel('Internal Energy e [code units]', fontsize=13, fontweight='bold')
    ax4.set_xlabel('Position x [code units]', fontsize=13, fontweight='bold')
    ax4.set_title('Internal Energy Profile', fontweight='bold', fontsize=14)
    ax4.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax4.legend(loc='best', fontsize=11, framealpha=0.95, edgecolor='gray')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if output_file:
        plt.savefig(output_file, dpi=200, bbox_inches='tight')
        print(f'✓ Saved: {output_file}')
    
    if show:
        plt.show()
    else:
        plt.close()

def main():
    args = parse_args()
    
    print('=' * 70)
    print('Strong Shock Multi-Method Comparison')
    print('=' * 70)
    print(f'Results directory: {args.results_dir}')
    print(f'Snapshot:          {args.snapshot if args.snapshot >= 0 else "final"}')
    print()
    
    # Load data from each method
    data_dict = {}
    time_str = None
    
    for method in METHODS.keys():
        method_dir = os.path.join(args.results_dir, method)
        if not os.path.exists(method_dir):
            print(f'⚠️  Skipping {method}: directory not found')
            continue
        
        snapshot_file = find_snapshot_files(method_dir, args.snapshot)
        if snapshot_file is None:
            print(f'⚠️  Skipping {method}: no snapshot found')
            continue
        
        data = load_snapshot(snapshot_file)
        data_dict[method] = data
        
        if time_str is None:
            time_str = data['metadata'].get('time', 'unknown')
        
        print(f'✓ Loaded {method}: {os.path.basename(snapshot_file)}')
    
    if not data_dict:
        print('\n❌ ERROR: No data loaded!')
        sys.exit(1)
    
    print()
    
    # Determine output file
    if args.output is None:
        output_file = os.path.join(args.results_dir, 'comparison', 'strong_shock_comparison.png')
    else:
        output_file = args.output
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Plot comparison
    plot_comparison(data_dict, time_str, output_file, show=not args.no_show)
    
    print()
    print('=' * 70)
    print('✓ Comparison Complete!')
    print('=' * 70)

if __name__ == '__main__':
    main()
