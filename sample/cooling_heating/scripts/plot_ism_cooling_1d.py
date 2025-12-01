#!/usr/bin/env python3
"""
Visualize 1D ISM cooling test results and compare with Koyama & Inutsuka (2000) curves
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
import sys

def load_snapshot(filename):
    """Load CSV snapshot"""
    df = pd.read_csv(filename, comment='#')
    return df

def load_ki_data(panel='a', var='T', column_density=1e19):
    """Load Koyama-Inutsuka digitized data"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, '..', 'results')
    
    # Map variable names to filenames
    var_map = {
        'T': f'f1a_curve_0.txt' if column_density == 1e19 else f'f1a_curve_1.txt',
        'P': f'f1a_curve_2.txt' if column_density == 1e19 else f'f1a_curve_3.txt',
        'xe': f'f1b_curve_0.txt' if column_density == 1e19 else f'f1b_curve_1.txt',
        'xH2': 'f1b_curve_2.txt',
        'xCO': 'f1b_curve_3.txt',
    }
    
    filepath = os.path.join(data_dir, var_map.get(var, 'f1a_curve_0.txt'))
    if not os.path.exists(filepath):
        return None, None
    
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]  # n_H, value

def plot_snapshot(snapshot_file, output_dir, show=False):
    """Plot single snapshot comparison with K&I curves"""
    df = load_snapshot(snapshot_file)
    
    # Extract data
    n_H = df['density'].values  # Assuming density is in cm^-3
    T = df['temperature'].values if 'temperature' in df.columns else df['pressure'] / df['density']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Temperature vs density
    ax = axes[0]
    ax.loglog(n_H, T, 'o', color='blue', alpha=0.6, label='SPH Simulation')
    
    # Overlay K&I curve
    n_ki, T_ki = load_ki_data('a', 'T', 1e19)
    if n_ki is not None:
        ax.loglog(n_ki, T_ki, 'k-', linewidth=2, label='K&I (2000) N_H=10^19')
    
    ax.set_xlabel(r'Number Density $n_{\rm H}$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel(r'Temperature $T$ [K]', fontsize=12)
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(10, 10000)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Temperature-Density Relation')
    
    # Pressure vs density
    ax = axes[1]
    P = df['pressure'].values
    ax.loglog(n_H, P, 'o', color='red', alpha=0.6, label='SPH Simulation')
    
    # Overlay K&I curve
    n_ki, P_ki = load_ki_data('a', 'P', 1e19)
    if n_ki is not None:
        ax.loglog(n_ki, P_ki, 'k-', linewidth=2, label='K&I (2000) N_H=10^19')
    
    ax.set_xlabel(r'Number Density $n_{\rm H}$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel(r'Pressure $P/k_{\rm B}$ [K cm$^{-3}$]', fontsize=12)
    ax.set_xlim(0.1, 1000)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Pressure-Density Relation')
    
    plt.tight_layout()
    
    # Save
    basename = os.path.basename(snapshot_file).replace('.csv', '.png')
    outfile = os.path.join(output_dir, basename)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved: {outfile}")
    
    if show:
        plt.show()
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_ism_cooling_1d.py <snapshot_file_or_directory>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    
    if os.path.isdir(input_path):
        # Process all snapshots in directory
        snapshot_files = sorted(glob.glob(os.path.join(input_path, 'snapshot_*.csv')))
        output_dir = os.path.join(input_path, 'plots')
        os.makedirs(output_dir, exist_ok=True)
        
        for snapshot_file in snapshot_files:
            print(f"Processing: {snapshot_file}")
            plot_snapshot(snapshot_file, output_dir, show=False)
        
        print(f"\nAll plots saved to: {output_dir}")
    else:
        # Single snapshot
        output_dir = os.path.dirname(input_path)
        plot_snapshot(input_path, output_dir, show=True)

if __name__ == '__main__':
    main()
