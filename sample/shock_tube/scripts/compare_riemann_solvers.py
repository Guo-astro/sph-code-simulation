#!/usr/bin/env python3
"""
Compare HLL vs Iterative Riemann solver results for shock tube test
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import sys

def load_snapshot(filepath):
    """Load a snapshot CSV file"""
    df = pd.read_csv(filepath, comment='#')
    return df

def compare_solvers(hll_dir, iterative_dir, output_dir):
    """Compare HLL and Iterative Riemann solver results"""
    
    hll_dir = Path(hll_dir)
    iterative_dir = Path(iterative_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get final snapshots
    hll_files = sorted(hll_dir.glob("snapshot_*.csv"))
    iter_files = sorted(iterative_dir.glob("snapshot_*.csv"))
    
    print(f"HLL snapshots: {len(hll_files)}")
    print(f"Iterative snapshots: {len(iter_files)}")
    
    # Load final states
    hll_final = load_snapshot(hll_files[-1])
    iter_final = load_snapshot(iter_files[-1])
    
    # Rename columns for consistency
    hll_final = hll_final.rename(columns={'pos_x': 'x', 'dens': 'rho', 'pres': 'P'})
    iter_final = iter_final.rename(columns={'pos_x': 'x', 'dens': 'rho', 'pres': 'P'})
    
    # Create comparison plot
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    fig.suptitle('Riemann Solver Comparison: HLL vs Iterative (van Leer 1997)\nSod Shock Tube at t=0.2', 
                 fontsize=14, fontweight='bold')
    
    # Density
    ax = axes[0]
    ax.scatter(hll_final['x'], hll_final['rho'], s=20, alpha=0.6, label='HLL', marker='o')
    ax.scatter(iter_final['x'], iter_final['rho'], s=20, alpha=0.6, label='Iterative', marker='s')
    ax.set_ylabel('Density', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_xlim(0, 1)
    
    # Velocity
    ax = axes[1]
    ax.scatter(hll_final['x'], hll_final['vel_x'], s=20, alpha=0.6, label='HLL', marker='o')
    ax.scatter(iter_final['x'], iter_final['vel_x'], s=20, alpha=0.6, label='Iterative', marker='s')
    ax.set_ylabel('Velocity', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_xlim(0, 1)
    
    # Pressure
    ax = axes[2]
    ax.scatter(hll_final['x'], hll_final['P'], s=20, alpha=0.6, label='HLL', marker='o')
    ax.scatter(iter_final['x'], iter_final['P'], s=20, alpha=0.6, label='Iterative', marker='s')
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Pressure', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_xlim(0, 1)
    
    plt.tight_layout()
    
    output_file = output_dir / "riemann_solver_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved comparison: {output_file}")
    
    # Create difference plot
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    fig.suptitle('Difference: Iterative - HLL\nSod Shock Tube at t=0.2', 
                 fontsize=14, fontweight='bold')
    
    # Interpolate to common grid for difference calculation
    x_common = np.linspace(0, 1, 1000)
    
    # Density difference
    ax = axes[0]
    hll_rho_interp = np.interp(x_common, hll_final['x'].values, hll_final['rho'].values)
    iter_rho_interp = np.interp(x_common, iter_final['x'].values, iter_final['rho'].values)
    diff_rho = iter_rho_interp - hll_rho_interp
    ax.plot(x_common, diff_rho, 'r-', linewidth=1.5)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Δ Density', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    rms_rho = np.sqrt(np.mean(diff_rho**2))
    ax.text(0.02, 0.95, f'RMS: {rms_rho:.6f}', transform=ax.transAxes, 
            fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Velocity difference
    ax = axes[1]
    hll_vx_interp = np.interp(x_common, hll_final['x'].values, hll_final['vel_x'].values)
    iter_vx_interp = np.interp(x_common, iter_final['x'].values, iter_final['vel_x'].values)
    diff_vx = iter_vx_interp - hll_vx_interp
    ax.plot(x_common, diff_vx, 'r-', linewidth=1.5)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Δ Velocity', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    rms_vx = np.sqrt(np.mean(diff_vx**2))
    ax.text(0.02, 0.95, f'RMS: {rms_vx:.6f}', transform=ax.transAxes, 
            fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Pressure difference
    ax = axes[2]
    hll_P_interp = np.interp(x_common, hll_final['x'].values, hll_final['P'].values)
    iter_P_interp = np.interp(x_common, iter_final['x'].values, iter_final['P'].values)
    diff_P = iter_P_interp - hll_P_interp
    ax.plot(x_common, diff_P, 'r-', linewidth=1.5)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Δ Pressure', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    rms_P = np.sqrt(np.mean(diff_P**2))
    ax.text(0.02, 0.95, f'RMS: {rms_P:.6f}', transform=ax.transAxes, 
            fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    output_file = output_dir / "riemann_solver_difference.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved difference plot: {output_file}")
    
    # Print statistics
    print("\n" + "="*60)
    print("Statistical Comparison (Iterative vs HLL)")
    print("="*60)
    print(f"Density RMS difference:  {rms_rho:.8f}")
    print(f"Velocity RMS difference: {rms_vx:.8f}")
    print(f"Pressure RMS difference: {rms_P:.8f}")
    print("="*60)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_riemann_solvers.py <hll_dir> <iterative_dir> [output_dir]")
        sys.exit(1)
    
    hll_dir = sys.argv[1]
    iter_dir = sys.argv[2]
    out_dir = sys.argv[3] if len(sys.argv) > 3 else "."
    
    print("="*60)
    print("Riemann Solver Comparison: HLL vs Iterative")
    print("="*60)
    print(f"HLL directory:       {hll_dir}")
    print(f"Iterative directory: {iter_dir}")
    print(f"Output directory:    {out_dir}")
    print("="*60 + "\n")
    
    compare_solvers(hll_dir, iter_dir, out_dir)
    
    print("\nComparison complete!")
