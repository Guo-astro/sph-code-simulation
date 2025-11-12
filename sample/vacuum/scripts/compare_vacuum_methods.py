#!/usr/bin/env python3
"""
Compare multiple SPH methods for vacuum test
Generates side-by-side comparison plots
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse
import sys


def read_snapshot(filename):
    """Read SPH snapshot from CSV file"""
    # Skip comment lines starting with #
    df = pd.read_csv(filename, comment='#')
    
    # Extract particle data (using new CSV column names)
    x = df['pos_x'].values
    rho = df['dens'].values
    v_x = df['vel_x'].values
    P = df['pres'].values
    u = df['ene'].values
    
    # Sort by position
    idx = np.argsort(x)
    return x[idx], rho[idx], v_x[idx], P[idx], u[idx]


def compute_analytical_solution(x, t, gamma=1.4):
    """
    Compute analytical vacuum test solution
    Two rarefaction waves propagating outward (Toro, Chapter 4)
    """
    x = np.asarray(x)
    rho = np.zeros_like(x)
    v = np.zeros_like(x)
    P = np.zeros_like(x)
    
    # Initial conditions
    rho_L, P_L, v_L = 1.0, 0.4, -2.0
    rho_R, P_R, v_R = 1.0, 0.4, 2.0
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)
    
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    x0 = 0.0
    
    # Two rarefaction waves moving apart
    # Left rarefaction: head moves at v_L - c_L, tail at v_L + c_L
    x_head_L = x0 + (v_L - c_L) * t
    x_tail_L = x0 + v_L * t + c_L * t
    
    # Right rarefaction: tail at v_R - c_R, head at v_R + c_R  
    x_tail_R = x0 + v_R * t - c_R * t
    x_head_R = x0 + (v_R + c_R) * t
    
    # Left uniform region
    mask = x <= x_head_L
    rho[mask] = rho_L
    v[mask] = v_L
    P[mask] = P_L
    
    # Left rarefaction fan
    mask = (x > x_head_L) & (x < x_tail_L)
    if np.any(mask):
        xi = (x[mask] - x0) / t
        c = 2.0 / gp1 * (c_L + gm1 / 2.0 * (v_L - xi))
        v[mask] = 2.0 / gp1 * (gm1 / 2.0 * v_L + c_L + xi)
        rho[mask] = rho_L * (c / c_L) ** (2.0 / gm1)
        P[mask] = P_L * (c / c_L) ** (2.0 * gamma / gm1)
    
    # Central region (between rarefaction waves)
    mask = (x >= x_tail_L) & (x <= x_tail_R)
    if np.any(mask):
        # Check if vacuum condition is met
        vacuum_criterion = v_R - v_L - 2.0 / gm1 * (c_L + c_R)
        if vacuum_criterion > 0:
            # True vacuum - use thermodynamically consistent floor values
            rho_floor = 1e-10
            # Use isentropic relation: P/P_0 = (ρ/ρ_0)^γ
            P[mask] = P_L * (rho_floor / rho_L) ** gamma
            rho[mask] = rho_floor
            v[mask] = 0.0
        else:
            # Intermediate state
            v_star = 0.5 * (v_L + v_R) + (c_L - c_R) / gm1
            v[mask] = v_star
            rho_floor = 1e-10
            P[mask] = P_L * (rho_floor / rho_L) ** gamma
            rho[mask] = rho_floor
    
    # Right rarefaction fan
    mask = (x > x_tail_R) & (x < x_head_R)
    if np.any(mask):
        xi = (x[mask] - x0) / t
        c = 2.0 / gp1 * (c_R - gm1 / 2.0 * (v_R - xi))
        v[mask] = 2.0 / gp1 * (gm1 / 2.0 * v_R - c_R + xi)
        rho[mask] = rho_R * (c / c_R) ** (2.0 / gm1)
        P[mask] = P_R * (c / c_R) ** (2.0 * gamma / gm1)
    
    # Right uniform
    mask = x >= x_head_R
    rho[mask] = rho_R
    v[mask] = v_R
    P[mask] = P_R
    
    e = P / (rho * gm1)
    return rho, v, P, e


def compare_methods(results_dir, methods, snapshot_num, output_file=None, show=True):
    """
    Compare multiple SPH methods
    
    Parameters:
    -----------
    results_dir : Path or str
        Base results directory
    methods : list of str
        Method names (subdirectories)
    snapshot_num : int
        Snapshot number to compare
    output_file : str, optional
        Output plot file
    show : bool
        Whether to display plot
    """
    results_dir = Path(results_dir)
    snapshot_name = f'snapshot_{snapshot_num:04d}.csv'
    
    # Read data for all methods
    data = {}
    for method in methods:
        snapshot_file = results_dir / method / snapshot_name
        if snapshot_file.exists():
            data[method] = read_snapshot(snapshot_file)
        else:
            print(f"Warning: {snapshot_file} not found, skipping {method}")
    
    if not data:
        print("Error: No data found for any method")
        return False
    
    # Assume dt_output = 0.01
    t = snapshot_num * 0.01
    
    # Compute analytical solution
    x_ana = np.linspace(-0.5, 0.5, 1000)
    rho_ana, v_ana, P_ana, e_ana = compute_analytical_solution(x_ana, t)
    
    # Create comparison plot
    n_methods = len(data)
    fig, axes = plt.subplots(4, n_methods, figsize=(5*n_methods, 12))
    
    if n_methods == 1:
        axes = axes.reshape(-1, 1)
    
    for col, (method, (x, rho, v, P, u)) in enumerate(data.items()):
        # Density
        axes[0, col].plot(x_ana, rho_ana, 'k-', label='Analytical', linewidth=2, alpha=0.7)
        axes[0, col].plot(x, rho, 'ro', label='SPH', markersize=2, alpha=0.6)
        axes[0, col].set_ylabel('Density ρ')
        axes[0, col].set_xlim([-0.5, 0.5])
        axes[0, col].set_ylim([-0.1, 1.2])
        axes[0, col].set_title(f'{method.upper()}\nt = {t:.5f}', fontweight='bold')
        axes[0, col].legend()
        axes[0, col].grid(True, alpha=0.3)
        
        # Velocity
        axes[1, col].plot(x_ana, v_ana, 'k-', linewidth=2, alpha=0.7)
        axes[1, col].plot(x, v, 'ro', markersize=2, alpha=0.6)
        axes[1, col].set_ylabel('Velocity v')
        axes[1, col].set_xlim([-0.5, 0.5])
        axes[1, col].set_ylim([-3.0, 3.0])
        axes[1, col].grid(True, alpha=0.3)
        
        # Pressure
        axes[2, col].plot(x_ana, P_ana, 'k-', linewidth=2, alpha=0.7)
        axes[2, col].plot(x, P, 'ro', markersize=2, alpha=0.6)
        axes[2, col].set_ylabel('Pressure P')
        axes[2, col].set_xlim([-0.5, 0.5])
        axes[2, col].set_ylim([-0.05, 0.5])
        axes[2, col].grid(True, alpha=0.3)
        
        # Internal Energy
        axes[3, col].plot(x_ana, e_ana, 'k-', linewidth=2, alpha=0.7)
        axes[3, col].plot(x, u, 'ro', markersize=2, alpha=0.6)
        axes[3, col].set_ylabel('Internal Energy e')
        axes[3, col].set_xlabel('Position x')
        axes[3, col].set_xlim([-0.5, 0.5])
        axes[3, col].set_ylim([0.0, 1.0])
        axes[3, col].grid(True, alpha=0.3)
    
    plt.suptitle('1D Vacuum Test - Multi-Method Comparison', fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Compare multiple SPH methods for vacuum test'
    )
    parser.add_argument('--results-dir', type=str, default='sample/vacuum/results',
                       help='Base results directory')
    parser.add_argument('--methods', nargs='+', 
                       default=['gsph_cubic', 'ssph_cubic', 'disph_cubic', 
                               'gdisph_cubic', 'gdisph_balsara_cubic'],
                       help='Methods to compare')
    parser.add_argument('--snapshot', type=int, default=14,
                       help='Snapshot number to compare (default: 14, final time)')
    parser.add_argument('-o', '--output', type=str,
                       help='Output plot file')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not display plot')
    
    args = parser.parse_args()
    
    success = compare_methods(
        args.results_dir,
        args.methods,
        args.snapshot,
        output_file=args.output,
        show=not args.no_show
    )
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
