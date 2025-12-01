#!/usr/bin/env python3
"""
Plot CNM Thermal Relaxation Results

Creates static plots showing:
1. Temperature evolution over time
2. Final temperature vs K&I equilibrium curve
3. Phase diagram (T vs n)

Usage:
    python3 plot_cnm_relaxation.py [results_dir]
    
Example:
    python3 plot_cnm_relaxation.py sample/cooling_heating/results/cnm_relaxation
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys

def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    metadata = {}
    
    with open(filename, 'r') as f:
        # Read metadata
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        # Read data
        f.seek(0)
        lines = [line for line in f.readlines() if not line.startswith('#')]
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
    
    data['metadata'] = metadata
    return data

def main():
    # Get data directory
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    else:
        data_dir = "sample/cooling_heating/results/cnm_relaxation"
    
    output_dir = os.path.join(data_dir, "plots")
    os.makedirs(output_dir, exist_ok=True)
    
    print('=' * 80)
    print('CNM Thermal Relaxation - Plotting Results')
    print('=' * 80)
    print(f'Data directory:   {data_dir}')
    print(f'Output directory: {output_dir}')
    print()
    
    # Find snapshots
    files = sorted(glob.glob(f'{data_dir}/snapshot_*.csv'))
    
    if len(files) == 0:
        print('ERROR: No snapshot files found!')
        print(f'Run simulation first: make -f sample/cooling_heating/Makefile.cooling cnm_run')
        sys.exit(1)
    
    print(f'Found {len(files)} snapshots')
    
    # Load all snapshots
    snapshots = []
    times = []
    
    for filename in files:
        data = load_snapshot(filename)
        snapshots.append(data)
        
        if 'time' in data['metadata']:
            times.append(float(data['metadata']['time']))
        else:
            idx = int(os.path.basename(filename).split('_')[1].split('.')[0])
            times.append(float(idx))
    
    times = np.array(times)
    
    # Calculate mean temperature over time
    mean_temps = []
    for data in snapshots:
        T = data['pres'] / data['dens']  # T = P / n_H (using actual column names)
        mean_temps.append(np.mean(T))
    
    mean_temps = np.array(mean_temps)
    
    # === PLOT 1: Temperature Evolution ===
    print('Creating temperature evolution plot...')
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(times, mean_temps, 'b-', linewidth=2, label='Mean Temperature')
    ax.axhline(107.0, color='gray', linestyle='--', label='Initial T = 107 K')
    ax.axhline(25.0, color='red', linestyle='--', label='K&I Equilibrium ~25 K (n=10 cm⁻³)')
    
    ax.set_xlabel('Time (code units)', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('CNM Thermal Relaxation - Temperature Evolution', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'temperature_evolution.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'  Saved: {output_file}')
    plt.close()
    
    # === PLOT 2: Phase Diagram (Final State) ===
    print('Creating phase diagram...')
    final_data = snapshots[-1]
    final_T = final_data['pres'] / final_data['dens']
    final_n = final_data['dens']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.scatter(final_n, final_T, c='blue', s=20, alpha=0.6, label='SPH Final State')
    ax.axhline(107.0, color='gray', linestyle='--', alpha=0.5, label='Initial T')
    ax.axhline(25.0, color='red', linestyle='--', alpha=0.5, label='K&I Equilibrium ~25 K')
    ax.axvline(10.0, color='green', linestyle=':', alpha=0.5, label='n = 10 cm⁻³')
    
    ax.set_xlabel('Density (cm⁻³)', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title(f'Phase Diagram - Final State (t={times[-1]:.1f})', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=11)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'phase_diagram.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'  Saved: {output_file}')
    plt.close()
    
    # === Summary ===
    print()
    print('=' * 80)
    print('✓ Plotting Complete!')
    print('=' * 80)
    print(f'Initial temperature: {mean_temps[0]:.1f} K')
    print(f'Final temperature:   {mean_temps[-1]:.1f} K')
    print(f'Cooling: {mean_temps[0] - mean_temps[-1]:.1f} K')
    print()
    print(f'Plots saved to: {output_dir}/')
    print()

if __name__ == '__main__':
    main()
