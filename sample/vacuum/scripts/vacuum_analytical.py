#!/usr/bin/env python3
"""
1D Vacuum Test Analytical Solution
Exact Riemann solver for vacuum formation problem

Initial conditions:
  Left (x <= 0):  rho=1.0, P=0.4, v=-2.0
  Right (x > 0):  rho=1.0, P=0.4, v=2.0

Physics: Fluids move apart creating vacuum region at x=0
Reference: Yuasa & Mori (2024), Section 4.2.2
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
from pathlib import Path


class VacuumRiemannSolver:
    """Exact Riemann solver for vacuum formation test"""
    
    def __init__(self, gamma=1.4):
        self.gamma = gamma
        
        # Initial conditions
        self.rho_L = 1.0
        self.P_L = 0.4
        self.v_L = -2.0
        
        self.rho_R = 1.0
        self.P_R = 0.4
        self.v_R = 2.0
        
        # Compute sound speeds
        self.c_L = np.sqrt(gamma * self.P_L / self.rho_L)
        self.c_R = np.sqrt(gamma * self.P_R / self.rho_R)
    
    def solve(self, x, t, x0=0.0):
        """
        Solve for density, velocity, pressure, and internal energy
        
        Two rarefaction waves propagating outward from initial discontinuity.
        Based on exact solution of Euler equations (Toro, Chapter 4).
        
        Parameters:
        -----------
        x : array-like
            Spatial positions
        t : float
            Time
        x0 : float
            Initial discontinuity location (default: 0.0)
        
        Returns:
        --------
        rho, v, P, e : arrays
            Density, velocity, pressure, internal energy
        """
        x = np.asarray(x)
        rho = np.zeros_like(x)
        v = np.zeros_like(x)
        P = np.zeros_like(x)
        e = np.zeros_like(x)
        
        gamma = self.gamma
        gm1 = gamma - 1.0
        gp1 = gamma + 1.0
        
        # Two rarefaction waves moving apart
        # Left rarefaction: head moves at v_L - c_L, tail at v_L + c_L
        x_head_L = x0 + (self.v_L - self.c_L) * t
        x_tail_L = x0 + self.v_L * t + self.c_L * t
        
        # Right rarefaction: tail at v_R - c_R, head at v_R + c_R  
        x_tail_R = x0 + self.v_R * t - self.c_R * t
        x_head_R = x0 + (self.v_R + self.c_R) * t
        
        # Left uniform region (unaffected by rarefaction)
        mask_L_uniform = x <= x_head_L
        rho[mask_L_uniform] = self.rho_L
        v[mask_L_uniform] = self.v_L
        P[mask_L_uniform] = self.P_L
        
        # Left rarefaction fan
        # For rarefaction: ρ/ρ_L = (c/c_L)^(2/(γ-1))
        # velocity relation: v - v_L = 2/(γ-1) * (c - c_L)
        # characteristic: dx/dt = v - c
        mask_L_fan = (x > x_head_L) & (x < x_tail_L)
        if np.any(mask_L_fan):
            # In rarefaction fan: x/t = v - c is constant along characteristics
            # For left-moving rarefaction: x/t = v - c
            # Solving: v = 2/(γ+1) * [(γ-1)/2 * v_L + c_L + (x-x0)/t]
            #          c = 2/(γ+1) * [c_L + (γ-1)/2 * (v_L - (x-x0)/t)]
            xi = (x[mask_L_fan] - x0) / t
            c = 2.0 / gp1 * (self.c_L + gm1 / 2.0 * (self.v_L - xi))
            v[mask_L_fan] = 2.0 / gp1 * (gm1 / 2.0 * self.v_L + self.c_L + xi)
            rho[mask_L_fan] = self.rho_L * (c / self.c_L) ** (2.0 / gm1)
            P[mask_L_fan] = self.P_L * (c / self.c_L) ** (2.0 * gamma / gm1)
        
        # Central region (between the two rarefaction waves)
        # If v_R - v_L > 2/(γ-1)*(c_L + c_R), vacuum forms
        # Otherwise, there's a low-density intermediate state
        mask_center = (x >= x_tail_L) & (x <= x_tail_R)
        if np.any(mask_center):
            # Check if vacuum condition is met
            vacuum_criterion = self.v_R - self.v_L - 2.0 / gm1 * (self.c_L + self.c_R)
            if vacuum_criterion > 0:
                # True vacuum - use thermodynamically consistent floor values
                rho_floor = 1e-10
                # Use isentropic relation: P/P_0 = (ρ/ρ_0)^γ to maintain consistency
                P[mask_center] = self.P_L * (rho_floor / self.rho_L) ** gamma
                rho[mask_center] = rho_floor
                v[mask_center] = 0.0
            else:
                # Intermediate state with v_* and P_* = 0 (or very small)
                # The two rarefaction tails meet
                # v_* = 0.5 * (v_L + v_R) + (c_L - c_R) / gm1
                # For symmetric case (c_L = c_R): v_* = (v_L + v_R) / 2
                v_star = 0.5 * (self.v_L + self.v_R) + (self.c_L - self.c_R) / gm1
                v[mask_center] = v_star
                rho_floor = 1e-10
                # Use isentropic relation for thermodynamic consistency
                P[mask_center] = self.P_L * (rho_floor / self.rho_L) ** gamma
                rho[mask_center] = rho_floor
        
        # Right rarefaction fan
        # For right-moving rarefaction: x/t = v + c
        mask_R_fan = (x > x_tail_R) & (x < x_head_R)
        if np.any(mask_R_fan):
            xi = (x[mask_R_fan] - x0) / t
            c = 2.0 / gp1 * (self.c_R - gm1 / 2.0 * (self.v_R - xi))
            v[mask_R_fan] = 2.0 / gp1 * (gm1 / 2.0 * self.v_R - self.c_R + xi)
            rho[mask_R_fan] = self.rho_R * (c / self.c_R) ** (2.0 / gm1)
            P[mask_R_fan] = self.P_R * (c / self.c_R) ** (2.0 * gamma / gm1)
        
        # Right uniform region
        mask_R_uniform = x >= x_head_R
        rho[mask_R_uniform] = self.rho_R
        v[mask_R_uniform] = self.v_R
        P[mask_R_uniform] = self.P_R
        
        # Compute internal energy: e = P / (rho * (gamma - 1))
        e = P / (rho * gm1)
        
        return rho, v, P, e


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


def plot_comparison(snapshot_file, output_file=None, show=True):
    """
    Plot vacuum test with analytical solution overlay
    
    Parameters:
    -----------
    snapshot_file : str
        Path to snapshot CSV file
    output_file : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
    """
    # Read SPH data
    x_sph, rho_sph, v_sph, P_sph, u_sph = read_snapshot(snapshot_file)
    
    # Extract time from filename (assuming format: snapshot_XXXX.csv)
    basename = os.path.basename(snapshot_file)
    snapshot_num = int(basename.split('_')[1].split('.')[0])
    t = snapshot_num * 0.01  # Assuming dt_output = 0.01
    
    # Compute analytical solution
    solver = VacuumRiemannSolver(gamma=1.4)
    x_analytical = np.linspace(-0.5, 0.5, 1000)
    rho_ana, v_ana, P_ana, e_ana = solver.solve(x_analytical, t)
    
    # Create 4-panel plot
    fig, axes = plt.subplots(4, 1, figsize=(10, 12))
    
    # Density
    axes[0].plot(x_analytical, rho_ana, 'k-', label='Analytical', linewidth=2)
    axes[0].plot(x_sph, rho_sph, 'ro', label='SPH', markersize=3, alpha=0.6)
    axes[0].set_ylabel('Density ρ', fontsize=12)
    axes[0].set_xlim([-0.5, 0.5])
    axes[0].set_ylim([-0.1, 1.2])
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title(f'1D Vacuum Test: t = {t:.5f}', fontsize=14, fontweight='bold')
    
    # Velocity
    axes[1].plot(x_analytical, v_ana, 'k-', label='Analytical', linewidth=2)
    axes[1].plot(x_sph, v_sph, 'ro', label='SPH', markersize=3, alpha=0.6)
    axes[1].set_ylabel('Velocity v', fontsize=12)
    axes[1].set_xlim([-0.5, 0.5])
    axes[1].set_ylim([-3.0, 3.0])
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Pressure
    axes[2].plot(x_analytical, P_ana, 'k-', label='Analytical', linewidth=2)
    axes[2].plot(x_sph, P_sph, 'ro', label='SPH', markersize=3, alpha=0.6)
    axes[2].set_ylabel('Pressure P', fontsize=12)
    axes[2].set_xlim([-0.5, 0.5])
    axes[2].set_ylim([-0.05, 0.5])
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    # Internal Energy
    axes[3].plot(x_analytical, e_ana, 'k-', label='Analytical', linewidth=2)
    axes[3].plot(x_sph, u_sph, 'ro', label='SPH', markersize=3, alpha=0.6)
    axes[3].set_ylabel('Internal Energy e', fontsize=12)
    axes[3].set_xlabel('Position x', fontsize=12)
    axes[3].set_xlim([-0.5, 0.5])
    axes[3].set_ylim([0.0, 1.0])
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    if show:
        plt.show()
    else:
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate analytical solution for 1D vacuum test'
    )
    parser.add_argument('snapshot', type=str, help='SPH snapshot CSV file')
    parser.add_argument('-o', '--output', type=str, help='Output plot file')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')
    
    args = parser.parse_args()
    
    # Generate plot
    plot_comparison(
        args.snapshot,
        output_file=args.output,
        show=not args.no_show
    )


if __name__ == '__main__':
    main()
