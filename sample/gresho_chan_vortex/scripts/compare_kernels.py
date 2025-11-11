#!/usr/bin/env python3
"""
Kernel Comparison Script for Gresho-Chan Vortex
Compares Wendland vs Cubic Spline kernels across 5 SPH methods
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from pathlib import Path

def load_snapshot(filepath):
    """Load a snapshot CSV file and extract metadata."""
    metadata = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
    
    data = pd.read_csv(filepath, comment='#')
    return data, metadata

def main():
    base_dir = 'sample/gresho_chan_vortex/results'
    output_dir = os.path.join(base_dir, 'kernel_comparison')
    os.makedirs(output_dir, exist_ok=True)
    
    methods = ['gsph', 'ssph', 'disph', 'gdisph', 'gdisph_balsara']
    method_labels = ['GSPH', 'SSPH', 'DISPH', 'GDISPH', 'GDISPH+Balsara']
    
    print("=" * 63)
    print("       Kernel Comparison: Wendland vs Cubic Spline")
    print("=" * 63)
    print(f"Base directory:   {base_dir}")
    print(f"Output directory: {output_dir}")
    print(f"SPH Methods:      {', '.join(method_labels)}")
    print()
    
    # Find common snapshots across all methods and kernels
    print("Scanning for snapshot files...")
    all_snapshot_numbers = {}
    
    for method in methods:
        wendland_dir = os.path.join(base_dir, f"{method}_wendland")
        cubic_dir = os.path.join(base_dir, method)
        
        wendland_snaps = set()
        cubic_snaps = set()
        
        if os.path.exists(wendland_dir):
            for f in os.listdir(wendland_dir):
                if f.startswith('snapshot_') and f.endswith('.csv'):
                    snap_num = int(f.split('_')[1].split('.')[0])
                    wendland_snaps.add(snap_num)
        
        if os.path.exists(cubic_dir):
            for f in os.listdir(cubic_dir):
                if f.startswith('snapshot_') and f.endswith('.csv'):
                    snap_num = int(f.split('_')[1].split('.')[0])
                    cubic_snaps.add(snap_num)
        
        common = wendland_snaps & cubic_snaps
        all_snapshot_numbers[method] = common
        print(f"  {method.upper():18s}: {len(cubic_snaps)} cubic, {len(wendland_snaps)} wendland, {len(common)} common")
    
    # Find snapshots common to all methods
    common_numbers = set.intersection(*all_snapshot_numbers.values())
    common_numbers = sorted(list(common_numbers))
    
    if not common_numbers:
        print("\n❌ Error: No common snapshots found across all methods!")
        return 1
    
    print(f"\nFound {len(common_numbers)} common snapshot numbers across all methods and kernels")
    
    # Select snapshots to plot (5 evenly spaced + final)
    num_snapshots = min(5, len(common_numbers))
    if len(common_numbers) > 1:
        indices = np.linspace(0, len(common_numbers)-2, num_snapshots-1, dtype=int)
        selected_snapshots = [common_numbers[i] for i in indices]
    else:
        selected_snapshots = []
    
    print(f"Generating comparison plots for {num_snapshots} snapshots...")
    print()
    
    for snap_num in selected_snapshots:
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
                    ax.set_xlim(-0.5, 0.5)
                    ax.set_ylim(-0.5, 0.5)
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
                    ax.set_xlim(-0.5, 0.5)
                    ax.set_ylim(-0.5, 0.5)
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
                ax.set_xlim(-0.5, 0.5)
                ax.set_ylim(-0.5, 0.5)
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
                ax.set_xlim(-0.5, 0.5)
                ax.set_ylim(-0.5, 0.5)
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
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
