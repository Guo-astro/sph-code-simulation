"""
SRGSPH Test Problems

Implements the test problems from Section 3 of the paper:
1. Sod Problem (Section 3.1.1, Eqs. 74-75)
2. Standard Relativistic Blast Wave (Section 3.1.2, Eqs. 76-77)
3. Strong Relativistic Blast Wave (Section 3.1.3, Eqs. 78-79)
"""

import numpy as np
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from srgsph.simulator import SRGSPH1D
from kitajima_solver import KitajimaRiemannSolver


def setup_sod_problem(n_left=3200, n_right=400, gamma_c=5.0/3.0, c=1.0):
    """
    Setup Sod shock tube problem (Section 3.1.1, Eqs. 74-75).
    
    Left state:  (P, n, vx, vt) = (1.0, 1.0, 0, 0)
    Right state: (P, n, vx, vt) = (0.1, 0.125, 0, 0)
    
    Args:
        n_left: Number of particles on left side
        n_right: Number of particles on right side
        gamma_c: Adiabatic index
        c: Speed of light
        
    Returns:
        SRGSPH1D simulator object
    """
    # Initial conditions
    P_left = 1.0
    n_left_density = 1.0
    v_left = 0.0
    
    P_right = 0.1
    n_right_density = 0.125
    v_right = 0.0
    
    # Discontinuity at x = 0.5
    x0 = 0.5
    
    # Create particle positions
    x_left = np.linspace(0.0, x0, n_left, endpoint=False)
    x_right = np.linspace(x0, 1.0, n_right, endpoint=True)
    positions = np.concatenate([x_left, x_right])
    
    # Set velocities, densities, pressures
    n_total = n_left + n_right
    velocities = np.zeros(n_total)
    densities = np.zeros(n_total)
    pressures = np.zeros(n_total)
    
    velocities[:n_left] = v_left
    velocities[n_left:] = v_right
    
    densities[:n_left] = n_left_density
    densities[n_left:] = n_right_density
    
    pressures[:n_left] = P_left
    pressures[n_left:] = P_right
    
    # Create simulator
    sim = SRGSPH1D(gamma_c=gamma_c, c=c, use_variable_h=True, use_muscl=True)
    sim.setup_particles(positions, velocities, densities, pressures)
    
    return sim


def setup_standard_blast_wave(n_left=5000, n_right=500, gamma_c=5.0/3.0, c=1.0):
    """
    Setup standard relativistic blast wave problem (Section 3.1.2, Eqs. 76-77).
    
    Left state:  (P, n, vx, vt) = (40/3, 10.0, 0, 0)
    Right state: (P, n, vx, vt) = (1e-6, 1.0, 0, 0)
    
    Args:
        n_left: Number of particles on left side
        n_right: Number of particles on right side
        gamma_c: Adiabatic index
        c: Speed of light
        
    Returns:
        SRGSPH1D simulator object
    """
    # Initial conditions
    P_left = 40.0 / 3.0
    n_left_density = 10.0
    v_left = 0.0
    
    P_right = 1e-6
    n_right_density = 1.0
    v_right = 0.0
    
    # Discontinuity at x = 0.5
    x0 = 0.5
    
    # Create particle positions
    x_left = np.linspace(0.0, x0, n_left, endpoint=False)
    x_right = np.linspace(x0, 1.0, n_right, endpoint=True)
    positions = np.concatenate([x_left, x_right])
    
    # Set velocities, densities, pressures
    n_total = n_left + n_right
    velocities = np.zeros(n_total)
    densities = np.zeros(n_total)
    pressures = np.zeros(n_total)
    
    velocities[:n_left] = v_left
    velocities[n_left:] = v_right
    
    densities[:n_left] = n_left_density
    densities[n_left:] = n_right_density
    
    pressures[:n_left] = P_left
    pressures[n_left:] = P_right
    
    # Create simulator
    sim = SRGSPH1D(gamma_c=gamma_c, c=c, use_variable_h=True, use_muscl=True)
    sim.setup_particles(positions, velocities, densities, pressures)
    
    return sim


def setup_strong_blast_wave(n_left=900, n_right=900, gamma_c=5.0/3.0, c=1.0):
    """
    Setup strong relativistic blast wave problem (Section 3.1.3, Eqs. 78-79).
    
    Left state:  (P, n, vx, vt) = (1000.0, 1.0, 0, 0)
    Right state: (P, n, vx, vt) = (0.01, 1.0, 0, 0)
    
    Args:
        n_left: Number of particles on left side
        n_right: Number of particles on right side
        gamma_c: Adiabatic index
        c: Speed of light
        
    Returns:
        SRGSPH1D simulator object
    """
    # Initial conditions
    P_left = 1000.0
    n_left_density = 1.0
    v_left = 0.0
    
    P_right = 0.01
    n_right_density = 1.0
    v_right = 0.0
    
    # Discontinuity at x = 0.5
    x0 = 0.5
    
    # Create particle positions
    x_left = np.linspace(0.0, x0, n_left, endpoint=False)
    x_right = np.linspace(x0, 1.0, n_right, endpoint=True)
    positions = np.concatenate([x_left, x_right])
    
    # Set velocities, densities, pressures
    n_total = n_left + n_right
    velocities = np.zeros(n_total)
    densities = np.zeros(n_total)
    pressures = np.zeros(n_total)
    
    velocities[:n_left] = v_left
    velocities[n_left:] = v_right
    
    densities[:n_left] = n_left_density
    densities[n_left:] = n_right_density
    
    pressures[:n_left] = P_left
    pressures[n_left:] = P_right
    
    # Create simulator
    sim = SRGSPH1D(gamma_c=gamma_c, c=c, use_variable_h=True, use_muscl=True)
    sim.setup_particles(positions, velocities, densities, pressures)
    
    return sim


def get_exact_solution(problem_type, t, gamma_c=5.0/3.0, c=1.0):
    """
    Get exact Riemann solution for comparison.
    
    Args:
        problem_type: 'sod', 'standard_blast', or 'strong_blast'
        t: Time
        gamma_c: Adiabatic index
        c: Speed of light
        
    Returns:
        Dictionary with exact solution arrays
    """
    solver = KitajimaRiemannSolver(gamma_c, c)
    
    if problem_type == 'sod':
        P_left, n_left, v_left = 1.0, 1.0, 0.0
        P_right, n_right, v_right = 0.1, 0.125, 0.0
    elif problem_type == 'standard_blast':
        P_left, n_left, v_left = 40.0/3.0, 10.0, 0.0
        P_right, n_right, v_right = 1e-6, 1.0, 0.0
    elif problem_type == 'strong_blast':
        P_left, n_left, v_left = 1000.0, 1.0, 0.0
        P_right, n_right, v_right = 0.01, 1.0, 0.0
    else:
        raise ValueError(f"Unknown problem type: {problem_type}")
    
    solver.set_initial_states(P_left, n_left, v_left, P_right, n_right, v_right)
    x, P, n, N, v, u, gamma, S, e = solver.solve(t, x0=0.5, n_points=400)
    
    return {
        'x': x,
        'P': P,
        'n': n,
        'N': N,
        'v': v,
        'u': u,
        'gamma': gamma,
        'S': S,
        'e': e
    }


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Run SRGSPH test problems')
    parser.add_argument('problem', choices=['sod', 'standard_blast', 'strong_blast'],
                        help='Test problem to run')
    parser.add_argument('--t_end', type=float, default=None,
                        help='End time (default: 0.35 for sod, 0.4 for blasts, 0.16 for strong)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output filename (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Setup problem
    print(f"Setting up {args.problem} problem...")
    
    if args.problem == 'sod':
        sim = setup_sod_problem()
        t_end = args.t_end if args.t_end is not None else 0.35
        output_file = args.output if args.output is not None else 'srgsph_sod.dat'
    elif args.problem == 'standard_blast':
        sim = setup_standard_blast_wave()
        t_end = args.t_end if args.t_end is not None else 0.4
        output_file = args.output if args.output is not None else 'srgsph_standard_blast.dat'
    elif args.problem == 'strong_blast':
        sim = setup_strong_blast_wave()
        t_end = args.t_end if args.t_end is not None else 0.16
        output_file = args.output if args.output is not None else 'srgsph_strong_blast.dat'
    
    # Run simulation
    print(f"Running simulation to t = {t_end}...")
    sim.run(t_end, output_freq=100)
    
    # Save solution
    sim.save_solution(output_file)
    
    print(f"\nDone! Solution saved to {output_file}")
