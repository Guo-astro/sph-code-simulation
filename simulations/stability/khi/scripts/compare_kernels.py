#!/usr/bin/env python3
"""
Kernel Comparison Script for Pairing Instability
Compares Wendland vs Cubic Spline kernels across 5 SPH methods
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_kernels.py <base_dir> <output_dir>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 63)
    print("       Kernel Comparison: Wendland vs Cubic Spline")
    print("=" * 63)
    print(f"Base directory:   {base_dir}")
    print(f"Output directory: {output_dir}")
    print(f"SPH Methods:      GSPH, SSPH, DISPH, GDISPH, GDISPH+Balsara")
    print()
    
    # Method definitions
    methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
    method_labels = ['GSPH\n(Godunov)', 'SSPH\n(Standard)', 'DISPH\n(Density-Indep)', 
                     'GDISPH\n(Godunov DISPH)', 'GDISPH+Balsara']
    
    # Find common snapshots
    print("Scanning for snapshot files...")
    method_snapshots = {}
    
    for method in methods:
        # Cubic spline snapshots
        cubic_dir = os.path.join(base_dir, method)
        if os.path.exists(cubic_dir):
            cubic_snaps = sorted([f for f in os.listdir(cubic_dir) if f.startswith('snapshot_') and f.endswith('.csv')])
            cubic_numbers = [int(f.split('_')[1].split('.')[0]) for f in cubic_snaps]
        else:
            cubic_numbers = []
        
        # Wendland snapshots
        wendland_dir = os.path.join(base_dir, f"{method}_wendland")
        if os.path.exists(wendland_dir):
            wendland_snaps = sorted([f for f in os.listdir(wendland_dir) if f.startswith('snapshot_') and f.endswith('.csv')])
            wendland_numbers = [int(f.split('_')[1].split('.')[0]) for f in wendland_snaps]
        else:
            wendland_numbers = []
        
        # Find common snapshots between both kernels
        common = sorted(set(cubic_numbers) & set(wendland_numbers))
        method_snapshots[method] = common
        print(f"  {method.upper():<18}: {len(cubic_numbers)} cubic, {len(wendland_numbers)} wendland, {len(common)} common")
    
    # Find snapshots common across all methods
    common_numbers = sorted(set.intersection(*[set(nums) for nums in method_snapshots.values()]))
    print(f"\nFound {len(common_numbers)} common snapshot numbers across all methods and kernels")
    
    # Select snapshots to plot (first, middle points, last)
    num_plots = min(5, len(common_numbers))
    if len(common_numbers) > 0:
        plot_indices = np.linspace(0, len(common_numbers)-1, num_plots, dtype=int)
        plot_snapshots = [common_numbers[i] for i in plot_indices]
    else:
        print("ERROR: No common snapshots found!")
        sys.exit(1)
    
    print(f"Generating comparison plots for {num_plots} snapshots...")
    print()
    
    # Generate comparison plots for selected snapshots
    for snap_num in plot_snapshots:
        fig, axes = plt.subplots(5, 2, figsize=(12, 24))
        fig.suptitle(f'Kernel Comparison: Wendland vs Cubic Spline - Snapshot {snap_num:04d}', 
                     fontsize=16, fontweight='bold', y=0.995)
        
        for idx, method in enumerate(methods):
            # Load Wendland data
            wendland_file = os.path.join(base_dir, f"{method}_wendland", f"snapshot_{snap_num:04d}.csv")
            if os.path.exists(wendland_file):
                data_w, metadata_w = load_snapshot(wendland_file)
                
                # Wendland plot (left column)
                ax = axes[idx, 0]
                x_col = 'x' if 'x' in data_w else 'pos_x'
                y_col = 'y' if 'y' in data_w else 'pos_y'
                
                if x_col in data_w and y_col in data_w and 'dens' in data_w:
                    log_dens = np.log10(data_w['dens'] + 1e-10)
                    scatter = ax.scatter(data_w[x_col], data_w[y_col], c=log_dens,
                                       cmap='hot', s=15, alpha=0.8, edgecolors='none')
                    ax.set_ylabel(method_labels[idx], fontsize=11, fontweight='bold', rotation=0,
                                 ha='right', va='center', labelpad=45)
                    ax.set_aspect('equal')
                    ax.set_xlim(0.0, 1.0)
                    ax.set_ylim(0.0, 1.0)
                    ax.set_facecolor('black')
                    ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
                    
                    if idx == 0:
                        ax.set_title('Wendland Kernel', fontsize=12, fontweight='bold', pad=15)
                    
                    if idx == len(methods) - 1:
                        ax.set_xlabel('x', fontsize=10)
                    else:
                        ax.set_xticklabels([])
            
            # Load Cubic Spline data
            cubic_file = os.path.join(base_dir, method, f"snapshot_{snap_num:04d}.csv")
            if os.path.exists(cubic_file):
                data_c, metadata_c = load_snapshot(cubic_file)
                
                # Cubic spline plot (right column)
                ax = axes[idx, 1]
                x_col = 'x' if 'x' in data_c else 'pos_x'
                y_col = 'y' if 'y' in data_c else 'pos_y'
                
                if x_col in data_c and y_col in data_c and 'dens' in data_c:
                    log_dens = np.log10(data_c['dens'] + 1e-10)
                    scatter = ax.scatter(data_c[x_col], data_c[y_col], c=log_dens,
                                       cmap='hot', s=15, alpha=0.8, edgecolors='none')
                    ax.set_aspect('equal')
                    ax.set_xlim(0.0, 1.0)
                    ax.set_ylim(0.0, 1.0)
                    ax.set_facecolor('black')
                    ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
                    
                    if idx == 0:
                        ax.set_title('Cubic Spline Kernel', fontsize=12, fontweight='bold', pad=15)
                        cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
                        cbar.set_label('log₁₀(Density)', fontsize=9)
                    
                    if idx == len(methods) - 1:
                        ax.set_xlabel('x', fontsize=10)
                    else:
                        ax.set_xticklabels([])
                    
                    ax.set_yticklabels([])
        
        plt.tight_layout(rect=[0, 0, 1, 0.99])
        output_file = os.path.join(output_dir, f'kernel_compare_snap{snap_num:04d}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'✓ Saved: {os.path.basename(output_file)}')
    
    # Create final state comparison
    final_snap = common_numbers[-1]
    print(f"\nCreating final state comparison (snapshot {final_snap:04d})...")
    
    fig, axes = plt.subplots(5, 2, figsize=(12, 24))
    fig.suptitle('Kernel Comparison: Wendland vs Cubic Spline - Final State', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    for idx, method in enumerate(methods):
        # Wendland
        wendland_file = os.path.join(base_dir, f"{method}_wendland", f"snapshot_{final_snap:04d}.csv")
        if os.path.exists(wendland_file):
            data_w, _ = load_snapshot(wendland_file)
            ax = axes[idx, 0]
            x_col = 'x' if 'x' in data_w else 'pos_x'
            y_col = 'y' if 'y' in data_w else 'pos_y'
            
            if x_col in data_w and y_col in data_w and 'dens' in data_w:
                log_dens = np.log10(data_w['dens'] + 1e-10)
                scatter = ax.scatter(data_w[x_col], data_w[y_col], c=log_dens,
                                   cmap='hot', s=15, alpha=0.8, edgecolors='none',
                                   vmin=np.percentile(log_dens, 1), vmax=np.percentile(log_dens, 99))
                ax.set_ylabel(method_labels[idx], fontsize=11, fontweight='bold', rotation=0,
                             ha='right', va='center', labelpad=45)
                ax.set_aspect('equal')
                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0.0, 1.0)
                ax.set_facecolor('black')
                ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
                
                if idx == 0:
                    ax.set_title('Wendland Kernel', fontsize=12, fontweight='bold', pad=15)
                if idx == len(methods) - 1:
                    ax.set_xlabel('x', fontsize=10)
                else:
                    ax.set_xticklabels([])
        
        # Cubic Spline
        cubic_file = os.path.join(base_dir, method, f"snapshot_{final_snap:04d}.csv")
        if os.path.exists(cubic_file):
            data_c, _ = load_snapshot(cubic_file)
            ax = axes[idx, 1]
            x_col = 'x' if 'x' in data_c else 'pos_x'
            y_col = 'y' if 'y' in data_c else 'pos_y'
            
            if x_col in data_c and y_col in data_c and 'dens' in data_c:
                log_dens = np.log10(data_c['dens'] + 1e-10)
                scatter = ax.scatter(data_c[x_col], data_c[y_col], c=log_dens,
                                   cmap='hot', s=15, alpha=0.8, edgecolors='none',
                                   vmin=np.percentile(log_dens, 1), vmax=np.percentile(log_dens, 99))
                ax.set_aspect('equal')
                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0.0, 1.0)
                ax.set_facecolor('black')
                ax.grid(True, alpha=0.2, color='white', linewidth=0.5)
                
                if idx == 0:
                    ax.set_title('Cubic Spline Kernel', fontsize=12, fontweight='bold', pad=15)
                    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
                    cbar.set_label('log₁₀(Density)', fontsize=9)
                if idx == len(methods) - 1:
                    ax.set_xlabel('x', fontsize=10)
                else:
                    ax.set_xticklabels([])
                
                ax.set_yticklabels([])
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    output_file = os.path.join(output_dir, 'kernel_compare_final.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'✓ Saved: {os.path.basename(output_file)}')
    
    print()
    print("=" * 63)
    print("       Kernel Comparison Complete!")
    print("=" * 63)
    print(f"Output saved to: {output_dir}")
    print("=" * 63)

if __name__ == '__main__':
    main()
