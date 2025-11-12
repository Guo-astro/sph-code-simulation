#!/usr/bin/env python3
"""
Strong Shock Test Analytical Solution Overlay

Exact Riemann solver for strong shock problem with extreme pressure ratio.

Initial conditions:
  Left (x < 0):   rho=1.0, P=1000.0, v=0.0
  Right (x >= 0): rho=1.0, P=0.1, v=0.0
  Pressure ratio: 10,000:1

Physics: Strong shock wave, contact discontinuity, and expansion wave
Reference: Toro (2009) "Riemann Solvers and Numerical Methods for Fluid Dynamics"
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
import sys
from pathlib import Path


class StrongShockRiemannSolver:
    """Exact Riemann solver for strong shock test"""
    
    def __init__(self, gamma=1.4):
        self.gamma = gamma
        
        # Initial conditions (strong shock)
        self.rho_L = 1.0
        self.P_L = 1000.0
        self.v_L = 0.0
        
        self.rho_R = 1.0
        self.P_R = 0.1
        self.v_R = 0.0
        
        # Compute sound speeds
        self.c_L = np.sqrt(gamma * self.P_L / self.rho_L)
        self.c_R = np.sqrt(gamma * self.P_R / self.rho_R)
        
        # Solve for star region (between contact and shock)
        self.P_star, self.v_star = self._solve_star_region()
    
    def _f(self, P, rho_K, P_K, c_K, side='L'):
        """
        Pressure function for exact Riemann solver
        Combines shock and rarefaction formulas
        """
        gamma = self.gamma
        gm1 = gamma - 1.0
        gp1 = gamma + 1.0
        
        if P > P_K:
            # Shock wave
            A = 2.0 / (gp1 * rho_K)
            B = gm1 / gp1 * P_K
            return (P - P_K) * np.sqrt(A / (P + B))
        else:
            # Rarefaction wave
            return 2.0 * c_K / gm1 * ((P / P_K) ** (gm1 / (2.0 * gamma)) - 1.0)
    
    def _fderiv(self, P, rho_K, P_K, c_K, side='L'):
        """Derivative of pressure function"""
        gamma = self.gamma
        gm1 = gamma - 1.0
        gp1 = gamma + 1.0
        
        if P > P_K:
            # Shock wave
            A = 2.0 / (gp1 * rho_K)
            B = gm1 / gp1 * P_K
            return np.sqrt(A / (P + B)) * (1.0 - (P - P_K) / (2.0 * (P + B)))
        else:
            # Rarefaction wave
            return (P / P_K) ** (-gp1 / (2.0 * gamma)) / (rho_K * c_K)
    
    def _solve_star_region(self, tol=1e-6, max_iter=100):
        """
        Iteratively solve for pressure and velocity in star region
        Uses Newton-Raphson method
        """
        # Initial guess (two-rarefaction approximation)
        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        
        P_0 = ((self.c_L + self.c_R - 0.5 * gm1 * (self.v_R - self.v_L)) / 
               (self.c_L / self.P_L ** (gm1 / (2.0 * self.gamma)) + 
                self.c_R / self.P_R ** (gm1 / (2.0 * self.gamma)))) ** (2.0 * self.gamma / gm1)
        
        P = P_0
        
        # Newton-Raphson iteration
        for i in range(max_iter):
            f_L = self._f(P, self.rho_L, self.P_L, self.c_L, 'L')
            f_R = self._f(P, self.rho_R, self.P_R, self.c_R, 'R')
            f = f_L + f_R + (self.v_R - self.v_L)
            
            if abs(f) < tol:
                break
            
            fderiv_L = self._fderiv(P, self.rho_L, self.P_L, self.c_L, 'L')
            fderiv_R = self._fderiv(P, self.rho_R, self.P_R, self.c_R, 'R')
            fderiv = fderiv_L + fderiv_R
            
            P_new = P - f / fderiv
            
            # Ensure positive pressure
            if P_new < 0:
                P_new = tol
            
            if abs(P_new - P) / P < tol:
                P = P_new
                break
            
            P = P_new
        
        # Compute star velocity
        v_star = 0.5 * (self.v_L + self.v_R + self._f(P, self.rho_R, self.P_R, self.c_R, 'R') - 
                        self._f(P, self.rho_L, self.P_L, self.c_L, 'L'))
        
        return P, v_star
    
    def solve(self, x, t, x0=0.0):
        """
        Solve exact Riemann problem at time t
        
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
        
        gamma = self.gamma
        gm1 = gamma - 1.0
        gp1 = gamma + 1.0
        
        P_star = self.P_star
        v_star = self.v_star
        
        # Compute post-shock density (right side - shock wave)
        if P_star > self.P_R:
            # Shock on right
            rho_star_R = self.rho_R * ((P_star / self.P_R + gm1 / gp1) / 
                                        (gm1 / gp1 * P_star / self.P_R + 1.0))
            S_R = self.v_R + self.c_R * np.sqrt(gp1 / (2.0 * gamma) * P_star / self.P_R + 
                                                 gm1 / (2.0 * gamma))
        else:
            # Rarefaction on right
            rho_star_R = self.rho_R * (P_star / self.P_R) ** (1.0 / gamma)
            S_HR = self.v_R + self.c_R  # Head of rarefaction
            c_star_R = self.c_R * (P_star / self.P_R) ** (gm1 / (2.0 * gamma))
            S_TR = v_star + c_star_R  # Tail of rarefaction
        
        # Compute post-expansion density (left side - expansion wave)
        if P_star > self.P_L:
            # Shock on left (unlikely for this problem)
            rho_star_L = self.rho_L * ((P_star / self.P_L + gm1 / gp1) / 
                                        (gm1 / gp1 * P_star / self.P_L + 1.0))
            S_L = self.v_L - self.c_L * np.sqrt(gp1 / (2.0 * gamma) * P_star / self.P_L + 
                                                 gm1 / (2.0 * gamma))
        else:
            # Rarefaction on left
            rho_star_L = self.rho_L * (P_star / self.P_L) ** (1.0 / gamma)
            S_HL = self.v_L - self.c_L  # Head of rarefaction
            c_star_L = self.c_L * (P_star / self.P_L) ** (gm1 / (2.0 * gamma))
            S_TL = v_star - c_star_L  # Tail of rarefaction
        
        # Sample the solution
        for i, xi in enumerate(x):
            s = (xi - x0) / t if t > 0 else 0.0
            
            if s < v_star:
                # Left of contact
                if P_star > self.P_L:
                    # Left shock
                    if s < S_L:
                        rho[i] = self.rho_L
                        v[i] = self.v_L
                        P[i] = self.P_L
                    else:
                        rho[i] = rho_star_L
                        v[i] = v_star
                        P[i] = P_star
                else:
                    # Left rarefaction
                    if s < S_HL:
                        rho[i] = self.rho_L
                        v[i] = self.v_L
                        P[i] = self.P_L
                    elif s < S_TL:
                        # Inside rarefaction fan
                        v[i] = 2.0 / gp1 * (self.c_L + gm1 / 2.0 * self.v_L + s)
                        c = 2.0 / gp1 * (self.c_L + gm1 / 2.0 * (self.v_L - s))
                        rho[i] = self.rho_L * (c / self.c_L) ** (2.0 / gm1)
                        P[i] = self.P_L * (c / self.c_L) ** (2.0 * gamma / gm1)
                    else:
                        rho[i] = rho_star_L
                        v[i] = v_star
                        P[i] = P_star
            else:
                # Right of contact
                if P_star > self.P_R:
                    # Right shock
                    if s < S_R:
                        rho[i] = rho_star_R
                        v[i] = v_star
                        P[i] = P_star
                    else:
                        rho[i] = self.rho_R
                        v[i] = self.v_R
                        P[i] = self.P_R
                else:
                    # Right rarefaction
                    if s < S_TR:
                        rho[i] = rho_star_R
                        v[i] = v_star
                        P[i] = P_star
                    elif s < S_HR:
                        # Inside rarefaction fan
                        v[i] = 2.0 / gp1 * (-self.c_R + gm1 / 2.0 * self.v_R + s)
                        c = 2.0 / gp1 * (self.c_R - gm1 / 2.0 * (self.v_R - s))
                        rho[i] = self.rho_R * (c / self.c_R) ** (2.0 / gm1)
                        P[i] = self.P_R * (c / self.c_R) ** (2.0 * gamma / gm1)
                    else:
                        rho[i] = self.rho_R
                        v[i] = self.v_R
                        P[i] = self.P_R
        
        # Compute internal energy
        e = P / (rho * gm1)
        
        return rho, v, P, e


def load_snapshot(filename):
    """Load SPH snapshot CSV file"""
    data = {}
    with open(filename, 'r') as f:
        # Read metadata
        metadata = {}
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        # Read data
        f.seek(0)
        lines = [l for l in f.readlines() if not l.startswith('#')]
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


def plot_with_analytical(snapshot_file, output_file=None, show=True, gamma=1.4):
    """Plot SPH results with analytical solution overlay"""
    
    # Load SPH data
    data = load_snapshot(snapshot_file)
    
    # Extract time from metadata (try multiple keys)
    metadata = data['metadata']
    time_str = metadata.get('Time (code)', 
                            metadata.get('time', 
                            metadata.get('Time', '0.0')))
    t = float(time_str)
    
    # Get SPH results
    x_sph = data['pos_x']
    sort_idx = np.argsort(x_sph)
    x_sph = x_sph[sort_idx]
    rho_sph = data['dens'][sort_idx]
    v_sph = data['vel_x'][sort_idx]
    P_sph = data['pres'][sort_idx]
    e_sph = data['ene'][sort_idx]
    
    # Compute analytical solution over full simulation domain
    # Domain: [-0.5, 0.5] as per C++ strong_shock.cpp
    solver = StrongShockRiemannSolver(gamma=gamma)
    x_exact = np.linspace(-0.5, 0.5, 1000)
    rho_exact, v_exact, P_exact, e_exact = solver.solve(x_exact, t, x0=0.0)
    
    # Create plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
    
    fig.suptitle(f'Strong Shock Test (P_ratio = 10,000) - SPH vs Exact Solution\nTime: {time_str}', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    # Density
    ax1.plot(x_exact, rho_exact, 'k-', linewidth=2.5, label='Exact', alpha=0.8)
    ax1.plot(x_sph, rho_sph, 'o', color='#0173B2', markersize=4, 
             markeredgewidth=0.5, markeredgecolor='white', label='SPH', alpha=0.7)
    ax1.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax1.set_ylabel('Density ρ', fontsize=13, fontweight='bold')
    ax1.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax1.set_title('Density Profile', fontweight='bold', fontsize=14)
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax1.legend(loc='best', fontsize=12, framealpha=0.95)
    
    # Velocity
    ax2.plot(x_exact, v_exact, 'k-', linewidth=2.5, label='Exact', alpha=0.8)
    ax2.plot(x_sph, v_sph, 'o', color='#DE8F05', markersize=4,
             markeredgewidth=0.5, markeredgecolor='white', label='SPH', alpha=0.7)
    ax2.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax2.set_ylabel('Velocity u', fontsize=13, fontweight='bold')
    ax2.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax2.set_title('Velocity Profile', fontweight='bold', fontsize=14)
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax2.legend(loc='best', fontsize=12, framealpha=0.95)
    
    # Pressure
    ax3.plot(x_exact, P_exact, 'k-', linewidth=2.5, label='Exact', alpha=0.8)
    ax3.plot(x_sph, P_sph, 'o', color='#029E73', markersize=4,
             markeredgewidth=0.5, markeredgecolor='white', label='SPH', alpha=0.7)
    ax3.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax3.set_ylabel('Pressure P', fontsize=13, fontweight='bold')
    ax3.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax3.set_title('Pressure Profile', fontweight='bold', fontsize=14)
    ax3.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax3.legend(loc='best', fontsize=12, framealpha=0.95)
    
    # Internal Energy
    ax4.plot(x_exact, e_exact, 'k-', linewidth=2.5, label='Exact', alpha=0.8)
    ax4.plot(x_sph, e_sph, 'o', color='#CC78BC', markersize=4,
             markeredgewidth=0.5, markeredgecolor='white', label='SPH', alpha=0.7)
    ax4.set_xlim(-0.5, 0.5)  # Domain from C++ strong_shock.cpp
    ax4.set_ylabel('Internal Energy e', fontsize=13, fontweight='bold')
    ax4.set_xlabel('Position x', fontsize=13, fontweight='bold')
    ax4.set_title('Internal Energy Profile', fontweight='bold', fontsize=14)
    ax4.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    ax4.legend(loc='best', fontsize=12, framealpha=0.95)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if output_file:
        plt.savefig(output_file, dpi=200, bbox_inches='tight')
        print(f'✓ Saved: {output_file}')
    
    if show:
        plt.show()
    else:
        plt.close()


def main():
    parser = argparse.ArgumentParser(description='Strong shock analytical solution overlay')
    parser.add_argument('snapshot', type=str, help='Snapshot CSV file')
    parser.add_argument('-o', '--output', type=str, default=None, help='Output filename')
    parser.add_argument('--gamma', type=float, default=1.4, help='Adiabatic index')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.snapshot):
        print(f'ERROR: Snapshot file not found: {args.snapshot}')
        sys.exit(1)
    
    # Determine output filename
    if args.output is None:
        base = os.path.splitext(args.snapshot)[0]
        output_file = f'{base}_analytical.png'
    else:
        output_file = args.output
    
    print('=' * 70)
    print('Strong Shock Analytical Solution Overlay')
    print('=' * 70)
    print(f'Snapshot: {args.snapshot}')
    print(f'Gamma:    {args.gamma}')
    print()
    
    plot_with_analytical(args.snapshot, output_file, show=not args.no_show, gamma=args.gamma)
    
    print()
    print('=' * 70)
    print('✓ Complete!')
    print('=' * 70)


if __name__ == '__main__':
    main()
