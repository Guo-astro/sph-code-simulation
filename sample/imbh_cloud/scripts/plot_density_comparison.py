#!/usr/bin/env python3
"""
Plot density profile comparison after relaxation
Shows SPH density vs analytic Lane-Emden solution
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

def load_lane_emden_solution(filepath='data/lane_emden/n1.5_3d.dat'):
    """Load the exact Lane-Emden n=1.5 solution"""
    data = np.loadtxt(filepath, skiprows=4)
    return {
        'xi': data[:, 0],
        'theta': data[:, 1],
        'dtheta': data[:, 2]
    }

def load_snapshot(filepath):
    """Load a single CSV snapshot"""
    # Find the line with "id," header dynamically
    with open(filepath, 'r') as f:
        skiprows = 0
        for line in f:
            if line.startswith('id,'):
                break
            skiprows += 1
    
    # Skip header lines + 1 to get to data
    data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows+1)
    return {
        'pos': data[:, 1:4],
        'rho': data[:, 11],
    }

def main():
    if len(sys.argv) < 2:
        results_dir = 'sample/imbh_cloud/results/lane_emden_50k_relax'
    else:
        results_dir = sys.argv[1]
    
    results_path = Path(results_dir)
    
    # Find latest snapshot
    snapshots = sorted(results_path.glob('snapshot_*.csv'))
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return
    
    final_snapshot = snapshots[-1]
    print(f"Loading final snapshot: {final_snapshot}")
    
    # Load Lane-Emden analytic solution
    le_sol = load_lane_emden_solution()
    
    # Parameters from Lane-Emden setup
    alpha = 0.273691
    rho_c = 1.43009692
    
    # Analytic density profile
    r_analytic = alpha * le_sol['xi']
    rho_analytic = rho_c * le_sol['theta']**1.5
    
    # Load SPH snapshot
    data = load_snapshot(final_snapshot)
    x, y, z = data['pos'][:, 0], data['pos'][:, 1], data['pos'][:, 2]
    rho_sph = data['rho']
    r_sph = np.sqrt(x**2 + y**2 + z**2)
    
    # Measure SPH central density
    central_mask = r_sph < 0.1
    rho_c_sph = rho_sph[central_mask].mean()
    
    print(f"SPH-measured central density: {rho_c_sph:.4f}")
    print(f"Analytic central density: {rho_c:.4f}")
    print(f"Ratio: {rho_c_sph/rho_c:.4f}")
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Density profiles
    ax1.scatter(r_sph, rho_sph, c='blue', s=1, alpha=0.3, label='SPH (corrected formula)')
    ax1.plot(r_analytic, rho_analytic, 'r--', linewidth=2, label=f'Lane-Emden analytic (ρ_c={rho_c:.2f})')
    ax1.set_xlabel('Radius [code units]', fontsize=12)
    ax1.set_ylabel('Density [code units]', fontsize=12)
    ax1.set_title('Radial Density Profile (After Relaxation)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1.0)
    ax1.set_ylim(0, 1.6)
    
    # Plot 2: Relative error
    bins = np.linspace(0, 1.0, 30)
    r_mid = 0.5 * (bins[1:] + bins[:-1])
    
    # Bin SPH data
    rho_binned = []
    for i in range(len(bins)-1):
        mask = (r_sph >= bins[i]) & (r_sph < bins[i+1])
        if mask.sum() > 0:
            rho_binned.append(rho_sph[mask].mean())
        else:
            rho_binned.append(np.nan)
    
    # Interpolate analytic solution
    rho_analytic_interp = np.interp(r_mid, r_analytic, rho_analytic)
    
    # Compute relative errors
    rel_error = 100 * (np.array(rho_binned) - rho_analytic_interp) / rho_analytic_interp
    
    ax2.plot(r_mid, rel_error, 'bo-', label='SPH vs Analytic', linewidth=2, markersize=6)
    ax2.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax2.axhline(-10, color='gray', linestyle=':', alpha=0.3, label='±10%')
    ax2.axhline(10, color='gray', linestyle=':', alpha=0.3)
    ax2.set_xlabel('Radius [code units]', fontsize=12)
    ax2.set_ylabel('Relative Error (%)', fontsize=12)
    ax2.set_title('SPH Density Error vs Analytic Solution', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1.0)
    ax2.set_ylim(-40, 10)
    
    # Add text box with statistics
    rms_error = np.sqrt(np.nanmean(rel_error**2))
    textstr = f'RMS error: {rms_error:.1f}%\nρ_c(SPH)/ρ_c(analytic): {rho_c_sph/rho_c:.3f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax2.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
    
    plt.suptitle(f'Lane-Emden Relaxation: Density Profile Analysis\n{final_snapshot.name}',
                 fontsize=15, fontweight='bold')
    plt.tight_layout()
    
    output_file = results_path / 'density_profile_comparison.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    print(f"  RMS error: {rms_error:.1f}%")
    
if __name__ == '__main__':
    main()
