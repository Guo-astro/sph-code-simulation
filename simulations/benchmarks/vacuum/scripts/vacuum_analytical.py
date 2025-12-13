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
        
        # If t == 0, return initial conditions
        if t <= 0.0:
            rho = np.where(x <= x0, self.rho_L, self.rho_R)
            v = np.where(x <= x0, self.v_L, self.v_R)
            P = np.where(x <= x0, self.P_L, self.P_R)
            e = P / (rho * gm1)
            return rho, v, P, e

        # Self-similar variable
        xi_all = (x - x0) / t

        # Two rarefaction waves moving apart (characteristic speeds in xi = x/t)
        # Initialize to right (default) state to guarantee full coverage
        rho[:] = self.rho_R
        v[:] = self.v_R
        P[:] = self.P_R

        # Check vacuum formation first
        vacuum_criterion = self.v_R - self.v_L - 2.0 / gm1 * (self.c_L + self.c_R)

        if vacuum_criterion > 0:
            # True vacuum: fans expand into center; compute fan edges using original relations
            x_head_L = x0 + (self.v_L - self.c_L) * t
            x_tail_L = x0 + (self.v_L + 2.0 * self.c_L / gm1) * t
            x_tail_R = x0 + (self.v_R - 2.0 * self.c_R / gm1) * t
            x_head_R = x0 + (self.v_R + self.c_R) * t
            # Convert to xi boundaries
            xi_tail_L = (x_tail_L - x0) / t
            xi_tail_R = (x_tail_R - x0) / t

            # Left uniform region
            mask_L_uniform = xi_all <= xi_head_L
            rho[mask_L_uniform] = self.rho_L
            v[mask_L_uniform] = self.v_L
            P[mask_L_uniform] = self.P_L

            # Left rarefaction fan (xi in (xi_head_L, xi_tail_L))
            mask_L_fan = (xi_all > xi_head_L) & (xi_all < xi_tail_L)
            if np.any(mask_L_fan):
                xi = xi_all[mask_L_fan]
                c = 2.0 / gp1 * (self.c_L + gm1 / 2.0 * (self.v_L - xi))
                c = np.maximum(c, 1e-10)
                v[mask_L_fan] = 2.0 / gp1 * (gm1 / 2.0 * self.v_L + self.c_L + xi)
                rho[mask_L_fan] = self.rho_L * (c / self.c_L) ** (2.0 / gm1)
                P[mask_L_fan] = self.P_L * (c / self.c_L) ** (2.0 * gamma / gm1)

            # Central vacuum region between xi_tail_R and xi_tail_L
            mask_center = (xi_all >= xi_tail_R) & (xi_all <= xi_tail_L)
            if np.any(mask_center):
                rho_floor = 1e-10
                P[mask_center] = self.P_L * (rho_floor / self.rho_L) ** gamma
                rho[mask_center] = rho_floor
                v[mask_center] = 0.0

            # Right rarefaction fan
            mask_R_fan = (xi_all > xi_tail_R) & (xi_all < xi_head_R)
            if np.any(mask_R_fan):
                xi = xi_all[mask_R_fan]
                c = 2.0 / gp1 * (self.c_R - gm1 / 2.0 * (self.v_R - xi))
                c = np.maximum(c, 1e-10)
                v[mask_R_fan] = 2.0 / gp1 * (gm1 / 2.0 * self.v_R - self.c_R + xi)
                rho[mask_R_fan] = self.rho_R * (c / self.c_R) ** (2.0 / gm1)
                P[mask_R_fan] = self.P_R * (c / self.c_R) ** (2.0 * gamma / gm1)

            # Right uniform region
            mask_R_uniform = xi_all >= xi_head_R
            rho[mask_R_uniform] = self.rho_R
            v[mask_R_uniform] = self.v_R
            P[mask_R_uniform] = self.P_R
            # Compute internal energy safely
            # Sanity check: ensure every point has been assigned a non-negative density
            zero_count = int((rho == 0.0).sum())
            if zero_count:
                print(f"Warning: {zero_count} analytic points unassigned - filling with right state")
                rho[rho == 0.0] = self.rho_R
                v[rho == 0.0] = self.v_R
                P[rho == 0.0] = self.P_R

            rho_safe = np.maximum(rho, 1e-10)
            e = P / (rho_safe * gm1)
            e = np.maximum(e, 0.0)
            return rho, v, P, e
        else:
            # No vacuum: compute star (intermediate) state using Riemann invariants
            # From Toro Section 4.3.2, Equation 4.53:
            # c_* = (c_L + c_R)/2 - (γ-1)/4 * (v_R - v_L)
            c_star = 0.5 * (self.c_L + self.c_R) - gm1 / 4.0 * (self.v_R - self.v_L)
            c_star = max(c_star, 1e-10)
            # From left rarefaction relation: v_* = v_L + 2/(γ-1)*(c_L - c_*)
            v_star = self.v_L + 2.0 / gm1 * (self.c_L - c_star)
            rho_star = self.rho_L * (c_star / self.c_L) ** (2.0 / gm1)
            P_star = self.P_L * (c_star / self.c_L) ** (2.0 * gamma / gm1)

            # Characteristic boundaries in xi
            xi_tail_L = v_star - c_star
            xi_tail_R = v_star + c_star
            xi_head_L = self.v_L - self.c_L
            xi_head_R = self.v_R + self.c_R

            # Left uniform
            mask_L_uniform = xi_all <= xi_head_L
            rho[mask_L_uniform] = self.rho_L
            v[mask_L_uniform] = self.v_L
            P[mask_L_uniform] = self.P_L

            # Left fan
            mask_L_fan = (xi_all > xi_head_L) & (xi_all < xi_tail_L)
            if np.any(mask_L_fan):
                xi = xi_all[mask_L_fan]
                c = 2.0 / gp1 * (self.c_L + gm1 / 2.0 * (self.v_L - xi))
                c = np.maximum(c, 1e-10)
                v[mask_L_fan] = 2.0 / gp1 * (gm1 / 2.0 * self.v_L + self.c_L + xi)
                rho[mask_L_fan] = self.rho_L * (c / self.c_L) ** (2.0 / gm1)
                P[mask_L_fan] = self.P_L * (c / self.c_L) ** (2.0 * gamma / gm1)

            # Central constant star region
            mask_center = (xi_all >= xi_tail_L) & (xi_all <= xi_tail_R)
            if np.any(mask_center):
                rho[mask_center] = rho_star
                v[mask_center] = v_star
                P[mask_center] = P_star

            # Right fan
            mask_R_fan = (xi_all > xi_tail_R) & (xi_all < xi_head_R)
            if np.any(mask_R_fan):
                xi = xi_all[mask_R_fan]
                c = 2.0 / gp1 * (self.c_R - gm1 / 2.0 * (self.v_R - xi))
                c = np.maximum(c, 1e-10)
                v[mask_R_fan] = 2.0 / gp1 * (gm1 / 2.0 * self.v_R - self.c_R + xi)
                rho[mask_R_fan] = self.rho_R * (c / self.c_R) ** (2.0 / gm1)
                P[mask_R_fan] = self.P_R * (c / self.c_R) ** (2.0 * gamma / gm1)

            # Right uniform
            mask_R_uniform = xi_all >= xi_head_R
            rho[mask_R_uniform] = self.rho_R
            v[mask_R_uniform] = self.v_R
            P[mask_R_uniform] = self.P_R
            # Compute internal energy safely
            # Sanity check: ensure every point has been assigned a non-negative density
            zero_count = int((rho == 0.0).sum())
            if zero_count:
                print(f"Warning: {zero_count} analytic points unassigned - filling with right state")
                rho[rho == 0.0] = self.rho_R
                v[rho == 0.0] = self.v_R
                P[rho == 0.0] = self.P_R

            rho_safe = np.maximum(rho, 1e-10)
            e = P / (rho_safe * gm1)
            e = np.maximum(e, 0.0)
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
