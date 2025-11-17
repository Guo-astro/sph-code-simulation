"""
Example: Running SRGSPH Sod Problem

This script demonstrates how to use the SRGSPH implementation
to solve the relativistic Sod shock tube problem.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend by default
import matplotlib.pyplot as plt
import sys
import os
import argparse

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from srgsph.simulator import SRGSPH1D
from kitajima_solver import KitajimaRiemannSolver


def main(show_plot=False):
    """Run Sod problem and compare with exact solution."""
    
    print("=" * 70)
    print("SRGSPH Example: Relativistic Sod Shock Tube Problem")
    print("=" * 70)
    print()
    
    # Parameters
    gamma_c = 5.0 / 3.0  # Adiabatic index
    c = 1.0              # Speed of light
    t_end = 0.35         # Final time
    
    # Initial conditions (from Eqs. 74-75 in paper)
    P_left = 1.0
    n_left = 1.0
    v_left = 0.0
    
    P_right = 0.1
    n_right = 0.125
    v_right = 0.0
    
    # Particle setup (minimal for fast demo)
    n_particles_left = 50   # Minimal for testing
    n_particles_right = 10  # Minimal for testing
    x0 = 0.5  # Discontinuity position
    
    print("Initial Conditions:")
    print(f"  Left state:  P = {P_left}, n = {n_left}, v = {v_left}")
    print(f"  Right state: P = {P_right}, n = {n_right}, v = {v_right}")
    print(f"  Discontinuity at x = {x0}")
    print(f"  Particles: {n_particles_left} (left) + {n_particles_right} (right)")
    print(f"  Note: Reduced particle count for fast demonstration")
    print()
    
    # Create particle positions
    x_left = np.linspace(0.0, x0, n_particles_left, endpoint=False)
    x_right = np.linspace(x0, 1.0, n_particles_right, endpoint=True)
    positions = np.concatenate([x_left, x_right])
    
    # Set primitive variables
    n_total = n_particles_left + n_particles_right
    velocities = np.concatenate([
        np.full(n_particles_left, v_left),
        np.full(n_particles_right, v_right)
    ])
    
    densities = np.concatenate([
        np.full(n_particles_left, n_left),
        np.full(n_particles_right, n_right)
    ])
    
    pressures = np.concatenate([
        np.full(n_particles_left, P_left),
        np.full(n_particles_right, P_right)
    ])
    
    # Create and setup simulator
    print("Setting up SRGSPH simulator...")
    sim = SRGSPH1D(
        gamma_c=gamma_c,
        c=c,
        C_CFL=0.3,
        eta=1.0,
        C_smooth=2.0,
        C_shock=3.0,
        C_cd=1.0,
        use_variable_h=True,
        use_muscl=False  # Disabled MUSCL for stability
    )
    
    sim.setup_particles(positions, velocities, densities, pressures)
    print(f"  Total particles: {len(sim.particles)}")
    print()
    
    # Run simulation
    print(f"Running simulation to t = {t_end}...")
    print("-" * 70)
    sim.run(t_end, output_freq=100)
    print("-" * 70)
    print()
    
    # Get SRGSPH solution
    sol = sim.get_solution()
    
    # Get exact Riemann solution
    print("Computing exact Riemann solution...")
    riemann = KitajimaRiemannSolver(gamma_c, c)
    riemann.set_initial_states(P_left, n_left, v_left, P_right, n_right, v_right)
    x_exact, P_exact, n_exact, N_exact, v_exact, u_exact, gamma_exact, S_exact, e_exact = \
        riemann.solve(t_end, x0=x0, n_points=400)
    print()
    
    # Plot comparison
    print("Creating comparison plot...")
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'SRGSPH vs Exact: Sod Problem at t = {t_end:.3f}', fontsize=16)
    
    # Sort SRGSPH solution by position
    idx = np.argsort(sol['x'])
    
    # Rest frame density
    ax = axes[0, 0]
    ax.scatter(sol['x'][idx], sol['n'][idx], s=15, alpha=0.6, label='SRGSPH', c='C0')
    ax.plot(x_exact, n_exact, 'k-', linewidth=2, label='Exact', alpha=0.7)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('n (rest frame density)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Pressure
    ax = axes[0, 1]
    ax.scatter(sol['x'][idx], sol['P'][idx], s=15, alpha=0.6, label='SRGSPH', c='C1')
    ax.plot(x_exact, P_exact, 'k-', linewidth=2, label='Exact', alpha=0.7)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('P (pressure)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Velocity
    ax = axes[0, 2]
    ax.scatter(sol['x'][idx], sol['v'][idx], s=15, alpha=0.6, label='SRGSPH', c='C2')
    ax.plot(x_exact, v_exact, 'k-', linewidth=2, label='Exact', alpha=0.7)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('v (velocity)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Lab frame density
    ax = axes[1, 0]
    ax.scatter(sol['x'][idx], sol['N'][idx], s=15, alpha=0.6, label='SRGSPH', c='C3')
    ax.plot(x_exact, N_exact, 'k-', linewidth=2, label='Exact', alpha=0.7)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('N (lab frame density)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Lorentz factor
    ax = axes[1, 1]
    ax.scatter(sol['x'][idx], sol['gamma'][idx], s=15, alpha=0.6, label='SRGSPH', c='C4')
    ax.plot(x_exact, gamma_exact, 'k-', linewidth=2, label='Exact', alpha=0.7)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('γ (Lorentz factor)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Smoothing length
    ax = axes[1, 2]
    ax.scatter(sol['x'][idx], sol['h'][idx], s=15, alpha=0.6, label='h', c='C5')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('h (smoothing length)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    output_file = 'srgsph_sod_example.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    
    # Save data
    data_file = 'srgsph_sod_example.dat'
    sim.save_solution(data_file)
    
    print()
    print("=" * 70)
    print("Example complete!")
    print("=" * 70)
    
    # Show plot only if requested
    if show_plot:
        matplotlib.use('TkAgg')  # Switch to interactive backend
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run SRGSPH Sod shock tube example')
    parser.add_argument('--show-plot', action='store_true', 
                       help='Display interactive plot window (blocks until closed)')
    args = parser.parse_args()
    
    main(show_plot=args.show_plot)
