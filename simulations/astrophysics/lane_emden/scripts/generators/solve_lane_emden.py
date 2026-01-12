#!/usr/bin/env python3
"""
Lane-Emden Equation Solver
===========================

Solves the Lane-Emden equation for polytropic stellar structure:
    1/ξ² d/dξ(ξ² dθ/dξ) + θ^n = 0
    
with boundary conditions:
    θ(0) = 1, dθ/dξ|_{ξ=0} = 0

Generates high-resolution tables for fast C++ interpolation.

Usage:
    python3 solve_lane_emden.py --n 1.5 --dim 3 --points 10000
    python3 solve_lane_emden.py --n 1.5 --dim 2 --points 10000
"""

import numpy as np
from scipy.integrate import solve_ivp
import argparse
import sys
from pathlib import Path


def lane_emden_ode(xi, y, n, dim):
    """
    Lane-Emden equation in first-order form.
    
    y[0] = θ(ξ)
    y[1] = dθ/dξ
    
    The equation is:
        d²θ/dξ² = -2/ξ dθ/dξ - θ^n     (for 3D sphere)
        d²θ/dξ² = -1/ξ dθ/dξ - θ^n     (for 2D disk)
    
    At ξ=0, apply L'Hôpital's rule to handle singularity:
        lim_{ξ→0} d²θ/dξ² = -1/3 θ^n   (3D)
        lim_{ξ→0} d²θ/dξ² = -1/2 θ^n   (2D)
    """
    theta = y[0]
    dtheta_dxi = y[1]
    
    # Handle singularity at ξ=0
    if xi < 1e-10:
        if dim == 3:
            d2theta_dxi2 = -max(0.0, theta)**n / 3.0
        elif dim == 2:
            d2theta_dxi2 = -max(0.0, theta)**n / 2.0
        else:
            raise ValueError(f"Unsupported dimension: {dim}")
    else:
        # General case - use max(0, theta) to avoid nan from negative**fractional
        theta_safe = max(0.0, theta)
        if dim == 3:
            d2theta_dxi2 = -(2.0 / xi) * dtheta_dxi - theta_safe**n
        elif dim == 2:
            d2theta_dxi2 = -(1.0 / xi) * dtheta_dxi - theta_safe**n
        else:
            raise ValueError(f"Unsupported dimension: {dim}")
    
    return [dtheta_dxi, d2theta_dxi2]


def find_first_zero(xi_array, theta_array):
    """Find first zero crossing of θ(ξ) by linear interpolation."""
    for i in range(len(theta_array) - 1):
        if theta_array[i] > 0 and theta_array[i+1] <= 0:
            # Linear interpolation to find exact zero
            frac = theta_array[i] / (theta_array[i] - theta_array[i+1])
            xi_1 = xi_array[i] + frac * (xi_array[i+1] - xi_array[i])
            return xi_1, i
    return None, None


def solve_lane_emden(n, dim=3, xi_max=30.0, num_points=10000, rtol=1e-10, atol=1e-12):
    """
    Solve Lane-Emden equation using adaptive Runge-Kutta.
    
    Parameters:
    -----------
    n : float
        Polytropic index
    dim : int
        Dimension (2 for disk, 3 for sphere)
    xi_max : float
        Maximum ξ to integrate to (will stop at first zero)
    num_points : int
        Number of output points
    rtol, atol : float
        Relative and absolute tolerances for ODE solver
        
    Returns:
    --------
    xi_array : ndarray
        Dimensionless radius values
    theta_array : ndarray
        θ(ξ) values
    dtheta_array : ndarray
        dθ/dξ values
    xi_1 : float
        First zero of θ (surface radius)
    dtheta_1 : float
        dθ/dξ at first zero
    """
    print(f"Solving Lane-Emden equation for n={n}, dim={dim}...")
    print(f"  Integration range: ξ ∈ [0, {xi_max}]")
    print(f"  Output points: {num_points}")
    print(f"  Tolerance: rtol={rtol}, atol={atol}")
    
    # Initial conditions: θ(0) = 1, dθ/dξ|_{ξ=0} = 0
    y0 = [1.0, 0.0]
    
    # Integration points (dense for accuracy)
    xi_eval = np.linspace(0, xi_max, num_points)
    
    # Event function to stop at first zero
    def theta_zero_event(xi, y):
        return y[0]  # Stop when θ = 0
    
    theta_zero_event.terminal = True
    theta_zero_event.direction = -1  # Only detect zero crossing from positive to negative
    
    # Solve ODE with high accuracy
    print("  Integrating...")
    sol = solve_ivp(
        fun=lambda xi, y: lane_emden_ode(xi, y, n, dim),
        t_span=(0, xi_max),
        y0=y0,
        method='DOP853',  # 8th order Runge-Kutta
        t_eval=xi_eval,
        events=theta_zero_event,
        rtol=rtol,
        atol=atol,
        dense_output=True
    )
    
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    
    xi_array = sol.t
    theta_array = sol.y[0]
    dtheta_array = sol.y[1]
    
    # Find first zero
    xi_1, idx_1 = find_first_zero(xi_array, theta_array)
    
    if xi_1 is None:
        print(f"  ⚠️  Warning: No zero found in range [0, {xi_max}]")
        print(f"     Try increasing xi_max or check polytrope index")
        xi_1 = xi_max
        dtheta_1 = dtheta_array[-1]
    else:
        # Interpolate dθ/dξ at first zero
        dtheta_1 = dtheta_array[idx_1] + \
                   (dtheta_array[idx_1+1] - dtheta_array[idx_1]) * \
                   (xi_1 - xi_array[idx_1]) / (xi_array[idx_1+1] - xi_array[idx_1])
        
        # Truncate arrays at first zero
        idx_end = idx_1 + 1
        xi_array = xi_array[:idx_end]
        theta_array = theta_array[:idx_end]
        dtheta_array = dtheta_array[:idx_end]
        
        print(f"  ✓ First zero found: ξ₁ = {xi_1:.10f}")
        print(f"  ✓ Derivative at zero: dθ/dξ|_{{ξ₁}} = {dtheta_1:.10f}")
        print(f"  ✓ Final table size: {len(xi_array)} points")
    
    return xi_array, theta_array, dtheta_array, xi_1, dtheta_1


def save_table(filename, xi_array, theta_array, dtheta_array, xi_1, dtheta_1, n, dim):
    """
    Save Lane-Emden table to file.
    
    Format:
        # Lane-Emden n=1.5 solution (DIM dimensions)
        # xi_1 = ...
        # dtheta_1 = ...
        # xi  theta  dtheta/dxi
        0.0  1.0  0.0
        ...
    """
    print(f"\nSaving table to {filename}...")
    
    with open(filename, 'w') as f:
        # Header
        f.write(f"# Lane-Emden n={n} solution ({dim}D)\n")
        f.write(f"# Generated with solve_lane_emden.py\n")
        f.write(f"# Points: {len(xi_array)}\n")
        f.write(f"# xi_1 = {xi_1:.15e}\n")
        f.write(f"# dtheta_1 = {dtheta_1:.15e}\n")
        f.write(f"# Format: xi  theta  dtheta/dxi\n")
        
        # Data
        for xi, theta, dtheta in zip(xi_array, theta_array, dtheta_array):
            f.write(f"{xi:.15e}  {theta:.15e}  {dtheta:.15e}\n")
    
    print(f"✓ Saved {len(xi_array)} points")
    
    # Print some verification info
    print(f"\nVerification:")
    print(f"  θ(0) = {theta_array[0]:.10f} (should be 1.0)")
    print(f"  dθ/dξ|_{{ξ=0}} = {dtheta_array[0]:.10e} (should be ~0)")
    print(f"  θ(ξ₁) = {theta_array[-1]:.10e} (should be ~0)")
    print(f"  ξ₁ = {xi_1:.10f}")
    print(f"  |dθ/dξ|_{{ξ₁}}| = {abs(dtheta_1):.10f}")


def verify_known_solutions(n, dim, xi_1, dtheta_1):
    """
    Verify against known analytical solutions.
    
    Known solutions:
    - n=0 (3D): ξ₁ = π, |dθ/dξ|_ξ₁ = 1
    - n=1 (3D): ξ₁ = π, |dθ/dξ|_ξ₁ = 1/π
    - n=5 (3D): ξ₁ = ∞ (no finite radius)
    - n=1.5 (3D): ξ₁ ≈ 3.6538, |dθ/dξ|_ξ₁ ≈ 0.2031
    """
    known = {
        (0, 3): (np.pi, 1.0),
        (1, 3): (np.pi, 1.0/np.pi),
        (1.5, 3): (3.6538, 0.2031),  # Approximate
    }
    
    key = (n, dim)
    if key in known:
        xi_1_known, dtheta_1_known = known[key]
        xi_1_error = abs(xi_1 - xi_1_known) / xi_1_known * 100
        dtheta_1_error = abs(abs(dtheta_1) - dtheta_1_known) / dtheta_1_known * 100
        
        print(f"\nComparison with known solution (n={n}, dim={dim}):")
        print(f"  ξ₁:         {xi_1:.10f} vs {xi_1_known:.10f} (error: {xi_1_error:.6f}%)")
        print(f"  |dθ/dξ|_ξ₁: {abs(dtheta_1):.10f} vs {dtheta_1_known:.10f} (error: {dtheta_1_error:.6f}%)")
        
        if xi_1_error < 0.01 and dtheta_1_error < 0.01:
            print(f"  ✓ Excellent agreement (<0.01% error)")
        elif xi_1_error < 0.1 and dtheta_1_error < 0.1:
            print(f"  ✓ Good agreement (<0.1% error)")
        else:
            print(f"  ⚠️  Warning: Larger than expected error")


def main():
    parser = argparse.ArgumentParser(
        description='Solve Lane-Emden equation and generate interpolation table',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate 3D n=1.5 table (sphere)
  python3 solve_lane_emden.py --n 1.5 --dim 3 --points 10000
  
  # Generate 2D n=1.5 table (disk)
  python3 solve_lane_emden.py --n 1.5 --dim 2 --points 10000
  
  # High-resolution table for extremely accurate simulations
  python3 solve_lane_emden.py --n 1.5 --dim 3 --points 50000 --rtol 1e-12
        """
    )
    
    parser.add_argument('--n', type=float, required=True,
                        help='Polytropic index (e.g., 1.5 for γ=5/3)')
    parser.add_argument('--dim', type=int, choices=[2, 3], required=True,
                        help='Dimension: 2 (disk) or 3 (sphere)')
    parser.add_argument('--points', type=int, default=10000,
                        help='Number of table points (default: 10000)')
    parser.add_argument('--xi-max', type=float, default=30.0,
                        help='Maximum ξ to integrate (default: 30.0)')
    parser.add_argument('--rtol', type=float, default=1e-10,
                        help='Relative tolerance (default: 1e-10)')
    parser.add_argument('--atol', type=float, default=1e-12,
                        help='Absolute tolerance (default: 1e-12)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output filename (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Solve Lane-Emden equation
    xi_array, theta_array, dtheta_array, xi_1, dtheta_1 = solve_lane_emden(
        n=args.n,
        dim=args.dim,
        xi_max=args.xi_max,
        num_points=args.points,
        rtol=args.rtol,
        atol=args.atol
    )
    
    # Determine output filename
    if args.output is None:
        # Auto-generate filename: n1.5_3d.dat or n1.5_2d.dat
        n_str = str(args.n).replace('.', '_')
        output_file = f"data/lane_emden/n{n_str}_{args.dim}d.dat"
    else:
        output_file = args.output
    
    # Create directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save table
    save_table(output_file, xi_array, theta_array, dtheta_array, xi_1, dtheta_1, args.n, args.dim)
    
    # Verify against known solutions
    verify_known_solutions(args.n, args.dim, xi_1, dtheta_1)
    
    print(f"\n{'='*60}")
    print(f"✓ SUCCESS: Lane-Emden table generated")
    print(f"{'='*60}")
    print(f"Output: {output_file}")
    print(f"Usage in C++: LaneEmdenData::load_from_file(\"{output_file}\")")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
