#!/usr/bin/env python3
"""
Animate CNM Thermal Relaxation Test Results

Creates animation showing thermal evolution of uniform CNM (n=10 cm^-3, T=107 K)
as it relaxes to equilibrium according to Koyama & Inutsuka (2000) cooling.

Usage:
    python3 animate_cnm_relaxation.py [results_dir]
    
Example:
    python3 animate_cnm_relaxation.py simulations/astrophysics/cooling_heating/results/cnm_relaxation
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from pathlib import Path

# Try to import animation tools
try:
    from matplotlib.animation import FuncAnimation, PillowWriter
    HAS_ANIMATION = True
except ImportError:
    HAS_ANIMATION = False
    print("Warning: Animation tools not available. Install pillow: pip install pillow")

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    metadata = {}
    
    with open(filename, 'r') as f:
        # Read metadata lines (start with #)
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
    
    # Add aliases for common names (CSV uses pos_x, pres, dens, ene)
    if 'pos_x' in data:
        data['x'] = data['pos_x']
    
    data['metadata'] = metadata
    return data

def calculate_temperature(data):
    """Calculate temperature from pressure and density: T = P / n_H"""
    if 'pres' in data and 'dens' in data:
        # In code units: P = n_H * T, so T = P / n_H
        return data['pres'] / data['dens']
    else:
        return None

def main():
    # Get data directory from command line or use default
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    else:
        data_dir = "simulations/astrophysics/cooling_heating/results/cnm_relaxation"
    
    output_dir = os.path.join(data_dir, "plots")
    os.makedirs(output_dir, exist_ok=True)
    
    print('=' * 80)
    print('CNM Thermal Relaxation - Animation')
    print('=' * 80)
    print(f'Data directory:   {data_dir}')
    print(f'Output directory: {output_dir}')
    print()
    
    # Find snapshot files
    print('Loading data files...')
    files = sorted(glob.glob(f'{data_dir}/snapshot_*.csv'))
    
    if len(files) == 0:
        print('ERROR: No snapshot files found!')
        print(f'Searched in: {data_dir}')
        print('Run simulation first: make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_run')
        sys.exit(1)
    
    print(f'Found {len(files)} snapshot files')
    print()
    
    # Load all snapshots
    print('Reading snapshots...')
    snapshots = []
    times = []
    
    iterator = tqdm(files) if HAS_TQDM else files
    for filename in iterator:
        data = load_snapshot(filename)
        snapshots.append(data)
        
        # Extract time from metadata
        if 'time' in data['metadata']:
            times.append(float(data['metadata']['time']))
        else:
            # Extract from filename: snapshot_0000.csv -> 0
            idx = int(os.path.basename(filename).split('_')[1].split('.')[0])
            times.append(float(idx))
    
    print(f'Time range: {times[0]:.2f} - {times[-1]:.2f} code units')
    print()
    
    # Create animation
    if not HAS_ANIMATION:
        print('ERROR: Animation not available. Install pillow:')
        print('  pip install pillow')
        sys.exit(1)
    
    print('Creating animation...')
    
    # Set up figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('CNM Thermal Relaxation - Koyama & Inutsuka (2000) Model C10', 
                 fontsize=14, fontweight='bold')
    
    ax_T, ax_P, ax_dens, ax_ene = axes.flatten()
    
    # Calculate temperature for all snapshots
    temperatures = []
    for data in snapshots:
        T = calculate_temperature(data)
        temperatures.append(T)
    
    # Get global ranges for consistent axes
    all_x = snapshots[0]['x']
    T_min = min(np.min(T) for T in temperatures if T is not None)
    T_max = max(np.max(T) for T in temperatures if T is not None)
    P_min = min(np.min(s['pres']) for s in snapshots)
    P_max = max(np.max(s['pres']) for s in snapshots)
    rho_min = min(np.min(s['dens']) for s in snapshots)
    rho_max = max(np.max(s['dens']) for s in snapshots)
    ene_min = min(np.min(s['ene']) for s in snapshots)
    ene_max = max(np.max(s['ene']) for s in snapshots)
    
    # Add some padding
    T_range = T_max - T_min
    P_range = P_max - P_min
    rho_range = rho_max - rho_min
    ene_range = ene_max - ene_min
    
    # Initial conditions from K&I Model C10
    n_init = 10.0  # cm^-3
    T_init = 107.0  # K
    T_eq = 25.0  # Approximate equilibrium temperature at n=10 cm^-3 from K&I Fig 1a
    
    def animate(frame):
        """Update plot for frame"""
        data = snapshots[frame]
        time = times[frame]
        T = temperatures[frame]
        
        # Clear axes
        for ax in axes.flatten():
            ax.clear()
        
        # Temperature
        ax_T.plot(data['x'], T, 'b-', linewidth=2, label='SPH')
        ax_T.axhline(T_init, color='gray', linestyle='--', alpha=0.5, label=f'Initial ({T_init} K)')
        ax_T.axhline(T_eq, color='red', linestyle='--', alpha=0.5, label=f'K&I Equilibrium (~{T_eq} K)')
        ax_T.set_xlabel('x (code units)', fontsize=11)
        ax_T.set_ylabel('Temperature (K)', fontsize=11)
        ax_T.set_title(f'Temperature Evolution (t = {time:.2f})', fontsize=12, fontweight='bold')
        ax_T.set_ylim(T_min - 0.1*T_range, T_max + 0.1*T_range)
        ax_T.grid(True, alpha=0.3)
        ax_T.legend(fontsize=9)
        
        # Pressure
        ax_P.plot(data['x'], data['pres'], 'r-', linewidth=2)
        ax_P.set_xlabel('x (code units)', fontsize=11)
        ax_P.set_ylabel('Pressure (code units)', fontsize=11)
        ax_P.set_title('Pressure', fontsize=12, fontweight='bold')
        ax_P.set_ylim(P_min - 0.1*P_range, P_max + 0.1*P_range)
        ax_P.grid(True, alpha=0.3)
        
        # Density
        ax_dens.plot(data['x'], data['dens'], 'g-', linewidth=2, label='SPH')
        ax_dens.axhline(n_init, color='gray', linestyle='--', alpha=0.5, label=f'Initial ({n_init} cm⁻³)')
        ax_dens.set_xlabel('x (code units)', fontsize=11)
        ax_dens.set_ylabel('Density (cm⁻³)', fontsize=11)
        ax_dens.set_title('Density', fontsize=12, fontweight='bold')
        ax_dens.set_ylim(rho_min - 0.1*rho_range, rho_max + 0.1*rho_range)
        ax_dens.grid(True, alpha=0.3)
        ax_dens.legend(fontsize=9)
        
        # Internal Energy
        ax_ene.plot(data['x'], data['ene'], 'purple', linewidth=2)
        ax_ene.set_xlabel('x (code units)', fontsize=11)
        ax_ene.set_ylabel('Internal Energy (code units)', fontsize=11)
        ax_ene.set_title('Internal Energy', fontsize=12, fontweight='bold')
        ax_ene.set_ylim(ene_min - 0.1*ene_range, ene_max + 0.1*ene_range)
        ax_ene.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return axes.flatten()
    
    # Create animation
    anim = FuncAnimation(fig, animate, frames=len(snapshots), 
                        interval=100, blit=False, repeat=True)
    
    # Save animation
    output_file = os.path.join(output_dir, 'thermal_evolution.gif')
    print(f'Saving animation to: {output_file}')
    
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=100)
    
    print()
    print('=' * 80)
    print('✓ Animation Complete!')
    print('=' * 80)
    print(f'Output: {output_file}')
    print()
    print('View with:')
    print(f'  open {output_file}')
    print()

if __name__ == '__main__':
    main()
