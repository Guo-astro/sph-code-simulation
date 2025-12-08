#!/usr/bin/env python3
"""
Quick visualization of isothermal slab evolution.
Compare grad-h vs no-grad-h simulations.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

def load_csv_snapshot(filepath):
    """Load particle data from CSV snapshot."""
    df = pd.read_csv(filepath, comment='#')
    return {
        'x': df['pos_x'].values,
        'rho': df['dens'].values,
        'vel': df['vel_x'].values if 'vel_x' in df.columns else np.zeros(len(df)),
    }

def plot_comparison(gradh_dir, no_gradh_dir, output_file=None):
    """Compare density profiles between grad-h and no-grad-h runs."""
    
    # Find snapshots
    gradh_snaps = sorted(glob.glob(os.path.join(gradh_dir, 'snapshot_*.csv')))
    no_gradh_snaps = sorted(glob.glob(os.path.join(no_gradh_dir, 'snapshot_*.csv')))
    
    if not gradh_snaps or not no_gradh_snaps:
        print(f"No snapshots found in {gradh_dir} or {no_gradh_dir}")
        return
    
    # Select snapshots to plot (initial, middle, final)
    indices = [0, len(gradh_snaps)//2, -1]
    times = [0, 25, 50]  # Approximate times
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    
    for col, (idx, t) in enumerate(zip(indices, times)):
        # Grad-h
        if idx < len(gradh_snaps):
            data = load_csv_snapshot(gradh_snaps[idx])
            ax = axes[0, col]
            ax.scatter(data['x'], data['rho'], s=2, alpha=0.7, c='blue')
            ax.set_xlabel('x')
            ax.set_ylabel('ρ')
            ax.set_title(f'With grad-h (t ≈ {t})')
            ax.set_ylim([0, 1.5])
            ax.grid(True, alpha=0.3)
        
        # No grad-h
        if idx < len(no_gradh_snaps):
            data = load_csv_snapshot(no_gradh_snaps[idx])
            ax = axes[1, col]
            ax.scatter(data['x'], data['rho'], s=2, alpha=0.7, c='red')
            ax.set_xlabel('x')
            ax.set_ylabel('ρ')
            ax.set_title(f'Without grad-h (t ≈ {t})')
            ax.set_ylim([0, 1.5])
            ax.grid(True, alpha=0.3)
    
    plt.suptitle('1D Isothermal Slab: Grad-h Correction Comparison', fontsize=14)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_file}")
    else:
        plt.show()

def plot_central_density_evolution(gradh_dir, no_gradh_dir, output_file=None):
    """Plot central density vs time for both runs."""
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    for label, results_dir, color in [
        ('With grad-h', gradh_dir, 'blue'),
        ('Without grad-h', no_gradh_dir, 'red')
    ]:
        snaps = sorted(glob.glob(os.path.join(results_dir, 'snapshot_*.csv')))
        
        times = []
        central_dens = []
        
        for snap in snaps:
            data = load_csv_snapshot(snap)
            # Extract time from filename
            snap_num = int(os.path.basename(snap).split('_')[1].split('.')[0])
            times.append(snap_num)  # Use snapshot number as proxy for time
            
            # Central density (particles near x=0)
            center_mask = np.abs(data['x']) < 0.1
            if np.sum(center_mask) > 0:
                central_dens.append(np.mean(data['rho'][center_mask]))
            else:
                central_dens.append(np.max(data['rho']))
        
        ax.plot(times, central_dens, 'o-', color=color, label=label, markersize=3)
    
    ax.set_xlabel('Snapshot number (∝ time)')
    ax.set_ylabel('Central density ρ(x≈0)')
    ax.set_title('Central Density Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_file}")
    else:
        plt.show()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize isothermal slab simulations')
    parser.add_argument('--gradh-dir', default='sample/diffusive_instability/results/slab_gradh',
                        help='Directory with grad-h results')
    parser.add_argument('--no-gradh-dir', default='sample/diffusive_instability/results/slab_no_gradh',
                        help='Directory with no-grad-h results')
    parser.add_argument('--output', '-o', default=None, help='Output file for plots')
    parser.add_argument('--central', action='store_true', help='Plot central density evolution')
    
    args = parser.parse_args()
    
    if args.central:
        plot_central_density_evolution(args.gradh_dir, args.no_gradh_dir, args.output)
    else:
        plot_comparison(args.gradh_dir, args.no_gradh_dir, args.output)
