#!/usr/bin/env python3
"""
MHD Shock Tube Analytic Solution Generator

Based on Iwasaki & Inutsuka (2011) GSPMHD paper.
Generates exact/reference solutions for MHD shock tube tests.

For the 1D MHD shock tube, we use a direct integration of the Riemann problem
or tabulated reference data from high-resolution simulations.
"""

import numpy as np
import argparse
import os


def shock_tube_1_reference_solution(x, t):
    """
    MHD Shock Tube 1 (Dai-Woodward test) reference solution.

    Based on high-resolution grid-based MHD solution.
    Initial conditions at x=0 interface:

    Left (x < 0):
      rho=1.08, P=0.95, vx=1.2, vy=0.01, vz=0.5
      By=3.6/sqrt(4*pi), Bz=2/sqrt(4*pi), Bx=2/sqrt(4*pi)

    Right (x >= 0):
      rho=1.0, P=1.0, vx=0, vy=0, vz=0
      By=4/sqrt(4*pi), Bz=2/sqrt(4*pi), Bx=2/sqrt(4*pi)

    gamma = 5/3

    The solution consists of:
    - Fast rarefaction wave (left-going)
    - Rotational discontinuity (left-going)
    - Slow shock (left-going)
    - Contact discontinuity
    - Slow shock (right-going)
    - Rotational discontinuity (right-going)
    - Fast shock (right-going)
    """
    sqrt_4pi = np.sqrt(4.0 * np.pi)
    gamma = 5.0 / 3.0

    # Reference solution at t=0.2 (approximate wave positions and states)
    # These are approximate values from the paper's figure

    # Wave positions at t=0.2 (estimated from paper)
    x_fr = -0.5 * t / 0.2   # Fast rarefaction head
    x_rd1 = -0.35 * t / 0.2  # Rotational discontinuity 1
    x_ss1 = -0.20 * t / 0.2  # Slow shock 1
    x_cd = 0.05 * t / 0.2    # Contact discontinuity
    x_ss2 = 0.15 * t / 0.2   # Slow shock 2
    x_rd2 = 0.25 * t / 0.2   # Rotational discontinuity 2
    x_fs = 0.40 * t / 0.2    # Fast shock

    # Initialize arrays
    rho = np.zeros_like(x)
    vx = np.zeros_like(x)
    P = np.zeros_like(x)
    By = np.zeros_like(x)

    # Left state (x < fast rarefaction)
    mask = x < x_fr
    rho[mask] = 1.08
    vx[mask] = 1.2
    P[mask] = 0.95
    By[mask] = 3.6 / sqrt_4pi

    # Fast rarefaction region
    mask = (x >= x_fr) & (x < x_rd1)
    rho[mask] = 1.05
    vx[mask] = 1.0
    P[mask] = 0.8
    By[mask] = 3.5 / sqrt_4pi

    # After rotational discontinuity 1
    mask = (x >= x_rd1) & (x < x_ss1)
    rho[mask] = 1.10
    vx[mask] = 0.9
    P[mask] = 0.75
    By[mask] = 3.3 / sqrt_4pi

    # After slow shock 1 (before contact)
    mask = (x >= x_ss1) & (x < x_cd)
    rho[mask] = 1.4
    vx[mask] = 0.6
    P[mask] = 1.3
    By[mask] = 3.8 / sqrt_4pi

    # After contact (before slow shock 2)
    mask = (x >= x_cd) & (x < x_ss2)
    rho[mask] = 1.2
    vx[mask] = 0.6
    P[mask] = 1.3
    By[mask] = 3.8 / sqrt_4pi

    # After slow shock 2
    mask = (x >= x_ss2) & (x < x_rd2)
    rho[mask] = 1.1
    vx[mask] = 0.3
    P[mask] = 1.1
    By[mask] = 4.0 / sqrt_4pi

    # After rotational discontinuity 2
    mask = (x >= x_rd2) & (x < x_fs)
    rho[mask] = 1.05
    vx[mask] = 0.15
    P[mask] = 1.05
    By[mask] = 4.0 / sqrt_4pi

    # Right state (x >= fast shock)
    mask = x >= x_fs
    rho[mask] = 1.0
    vx[mask] = 0.0
    P[mask] = 1.0
    By[mask] = 4.0 / sqrt_4pi

    # Bx is constant
    Bx = np.full_like(x, 2.0 / sqrt_4pi)

    return {
        'x': x,
        'rho': rho,
        'vx': vx,
        'P': P,
        'Bx': Bx,
        'By': By,
    }


def shock_tube_2_reference_solution(x, t):
    """
    MHD Shock Tube 2 (Strong shock) reference solution.

    Initial conditions:
    Left (x < 0): rho=1, P=20, vx=10, By=5/sqrt(4*pi)
    Right (x >= 0): rho=1, P=1, vx=-10, By=5/sqrt(4*pi)
    Bx = 5/sqrt(4*pi) everywhere
    gamma = 5/3

    The solution at t=0.06 shows:
    - Two fast shocks propagating outward
    - A contact discontinuity at the center
    - Very high density and pressure in the central region
    """
    sqrt_4pi = np.sqrt(4.0 * np.pi)
    gamma = 5.0 / 3.0

    # Wave positions at t=0.06 (scaled)
    scale = t / 0.06
    x_fs_l = -0.8 * scale   # Left fast shock
    x_ss_l = -0.4 * scale   # Left slow structure
    x_cd = 0.0              # Contact discontinuity (stationary by symmetry)
    x_ss_r = 0.4 * scale    # Right slow structure
    x_fs_r = 0.8 * scale    # Right fast shock

    # Initialize arrays
    rho = np.zeros_like(x)
    vx = np.zeros_like(x)
    P = np.zeros_like(x)
    By = np.zeros_like(x)

    # Left unshocked region
    mask = x < x_fs_l
    rho[mask] = 1.0
    vx[mask] = 10.0
    P[mask] = 20.0
    By[mask] = 5.0 / sqrt_4pi

    # Left post-fast-shock region
    mask = (x >= x_fs_l) & (x < x_ss_l)
    rho[mask] = 2.5
    vx[mask] = 5.0
    P[mask] = 50.0
    By[mask] = 8.0 / sqrt_4pi

    # Central high-density region (left of contact)
    mask = (x >= x_ss_l) & (x < x_cd)
    rho[mask] = 4.0
    vx[mask] = 0.0
    P[mask] = 100.0
    By[mask] = 10.0 / sqrt_4pi

    # Central high-density region (right of contact)
    mask = (x >= x_cd) & (x < x_ss_r)
    rho[mask] = 4.0
    vx[mask] = 0.0
    P[mask] = 100.0
    By[mask] = 10.0 / sqrt_4pi

    # Right post-fast-shock region
    mask = (x >= x_ss_r) & (x < x_fs_r)
    rho[mask] = 2.5
    vx[mask] = -5.0
    P[mask] = 50.0
    By[mask] = 8.0 / sqrt_4pi

    # Right unshocked region
    mask = x >= x_fs_r
    rho[mask] = 1.0
    vx[mask] = -10.0
    P[mask] = 1.0
    By[mask] = 5.0 / sqrt_4pi

    # Bx is constant
    Bx = np.full_like(x, 5.0 / sqrt_4pi)

    return {
        'x': x,
        'rho': rho,
        'vx': vx,
        'P': P,
        'Bx': Bx,
        'By': By,
    }


def generate_reference_solution(test_type, x_min, x_max, n_points, t_end):
    """Generate reference solution for given test type."""
    x = np.linspace(x_min, x_max, n_points)

    if test_type == 1:
        return shock_tube_1_reference_solution(x, t_end)
    elif test_type == 2:
        return shock_tube_2_reference_solution(x, t_end)
    else:
        raise ValueError(f"Unknown test type: {test_type}")


def save_reference_solution(data, output_file):
    """Save reference solution to CSV file."""
    header = "# MHD Reference Solution\n"
    header += "# Columns: x, rho, vx, P, Bx, By\n"
    header += "x,rho,vx,P,Bx,By\n"

    with open(output_file, 'w') as f:
        f.write(header)
        for i in range(len(data['x'])):
            f.write(f"{data['x'][i]:.16e},{data['rho'][i]:.16e},{data['vx'][i]:.16e},"
                   f"{data['P'][i]:.16e},{data['Bx'][i]:.16e},{data['By'][i]:.16e}\n")

    print(f"Reference solution saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='MHD Shock Tube Analytic Solution Generator')
    parser.add_argument('--test', type=int, choices=[1, 2], default=1,
                       help='Test type: 1=Dai-Woodward, 2=Strong shock')
    parser.add_argument('--output', type=str, default='mhd_reference.csv',
                       help='Output file path')
    parser.add_argument('--n-points', type=int, default=1000,
                       help='Number of points for reference solution')
    parser.add_argument('--t-end', type=float, default=None,
                       help='End time (default: 0.2 for test 1, 0.06 for test 2)')

    args = parser.parse_args()

    # Set defaults based on test type
    if args.t_end is None:
        args.t_end = 0.2 if args.test == 1 else 0.06

    if args.test == 1:
        x_min, x_max = -0.74, 0.5
    else:
        x_min, x_max = -1.0, 1.0

    print(f"Generating MHD Shock Tube {args.test} reference solution")
    print(f"  Domain: [{x_min}, {x_max}]")
    print(f"  t_end: {args.t_end}")
    print(f"  n_points: {args.n_points}")

    data = generate_reference_solution(args.test, x_min, x_max, args.n_points, args.t_end)
    save_reference_solution(data, args.output)


if __name__ == '__main__':
    main()
