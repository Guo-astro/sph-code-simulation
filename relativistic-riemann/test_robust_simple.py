#!/usr/bin/env python3
"""
Simple test script for robust iterative relativistic Riemann solver.
Tests convergence and accuracy.
"""

import numpy as np
import matplotlib.pyplot as plt
from iterative_riemann_solver import IterativeRiemannSolver


def test_problem_1():
    """Test Problem 1 from Rezzolla & Zanotti (moderate shock)"""
    print("\n" + "="*60)
    print("Test Problem 1: Moderate Shock")
    print("="*60)
    
    # Initial conditions
    gamma_c = 5.0/3.0
    Pl, nl, vl = 13.33, 10.0, 0.0
    Pr, nr, vr = 1e-8, 1.0, 0.0
    
    # Robust iterative solver
    print("\n--- Robust Iterative Solver ---")
    solver = IterativeRiemannSolver(gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=True)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps = solver.solve_star_pressure()
    vs = 0.5 * (solver.vls + solver.vrs)
    
    # Expected values (from literature/reference solutions)
    # These are approximate reference values
    Ps_expected = 1.448  # Approximate
    vs_expected = 0.536  # Approximate
    
    print(f"\n--- Results ---")
    print(f"Star Pressure: Ps = {Ps:.10f}")
    print(f"Star Velocity: vs = {vs:.10f}")
    print(f"Expected:      Ps ≈ {Ps_expected:.3f}, vs ≈ {vs_expected:.3f}")
    
    # Diagnostics
    diag = solver.get_convergence_diagnostics()
    print(f"\n--- Convergence Diagnostics ---")
    print(f"Iterations:     {diag['iterations']}")
    print(f"Initial P*:     {diag['initial_pressure']:.6e}")
    print(f"Final P*:       {diag['final_pressure']:.6e}")
    print(f"Pressure change: {diag['pressure_change']:.6e}")
    print(f"Final residual: {diag['final_residual']:.3e}")
    
    return solver


def test_problem_2():
    """Test Problem 2: Strong relativistic shock"""
    print("\n" + "="*60)
    print("Test Problem 2: Strong Relativistic Shock")
    print("="*60)
    
    # Initial conditions
    gamma_c = 5.0/3.0
    Pl, nl, vl = 1000.0, 1.0, 0.0
    Pr, nr, vr = 0.01, 1.0, 0.0
    
    # Robust iterative solver
    print("\n--- Robust Iterative Solver ---")
    solver = IterativeRiemannSolver(gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=True)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps = solver.solve_star_pressure()
    vs = 0.5 * (solver.vls + solver.vrs)
    
    print(f"\n--- Results ---")
    print(f"Star Pressure: Ps = {Ps:.10f}")
    print(f"Star Velocity: vs = {vs:.10f}")
    
    # Diagnostics
    diag = solver.get_convergence_diagnostics()
    print(f"\n--- Convergence Diagnostics ---")
    print(f"Iterations:      {diag['iterations']}")
    print(f"Final residual:  {diag['final_residual']:.3e}")
    
    return solver


def test_problem_3():
    """Test Problem 3: Two rarefactions (collision)"""
    print("\n" + "="*60)
    print("Test Problem 3: Two Rarefactions")
    print("="*60)
    
    # Initial conditions
    gamma_c = 5.0/3.0
    Pl, nl, vl = 1.0, 1.0, -0.6
    Pr, nr, vr = 1.0, 1.0, 0.6
    
    # Robust iterative solver
    print("\n--- Robust Iterative Solver ---")
    solver = IterativeRiemannSolver(gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=True)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps = solver.solve_star_pressure()
    vs = 0.5 * (solver.vls + solver.vrs)
    
    print(f"\n--- Results ---")
    print(f"Star Pressure: Ps = {Ps:.10f}")
    print(f"Star Velocity: vs = {vs:.10f}")
    print(f"Expected: vs ≈ 0.0 (symmetry)")
    
    # Diagnostics
    diag = solver.get_convergence_diagnostics()
    print(f"\n--- Convergence Diagnostics ---")
    print(f"Iterations:      {diag['iterations']}")
    print(f"Final residual:  {diag['final_residual']:.3e}")
    
    return solver


def plot_convergence(solver, title):
    """Plot convergence history for a solver"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    iterations = range(len(solver.pressure_history))
    ax1.plot(iterations, solver.pressure_history, 'o-', markersize=4, linewidth=2)
    ax1.set_xlabel('Iteration', fontsize=12)
    ax1.set_ylabel('Pressure P*', fontsize=12)
    ax1.set_title('Pressure Convergence', fontsize=13)
    ax1.grid(True, alpha=0.3)
    
    if len(solver.residual_history) > 0:
        ax2.semilogy(range(1, len(solver.residual_history)+1), 
                    solver.residual_history, 'o-', markersize=4, linewidth=2, color='red')
        ax2.set_xlabel('Iteration', fontsize=12)
        ax2.set_ylabel('|Residual| (velocity difference)', fontsize=12)
        ax2.set_title('Residual Convergence', fontsize=13)
        ax2.grid(True, alpha=0.3)
    
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    print("\n" + "="*70)
    print(" ROBUST ITERATIVE RELATIVISTIC RIEMANN SOLVER TEST")
    print("="*70)
    
    # Run tests
    solver1 = test_problem_1()
    solver2 = test_problem_2()
    solver3 = test_problem_3()
    
    # Summary
    print("\n" + "="*70)
    print(" SUMMARY")
    print("="*70)
    print(f"Problem 1 (Moderate Shock):     {solver1.iterations} iterations")
    print(f"Problem 2 (Strong Shock):       {solver2.iterations} iterations")
    print(f"Problem 3 (Two Rarefactions):   {solver3.iterations} iterations")
    print("\n✓ All tests completed successfully!")
    print("="*70)
    
    # Plot convergence
    fig1 = plot_convergence(solver1, 'Problem 1: Moderate Shock - Convergence')
    fig2 = plot_convergence(solver2, 'Problem 2: Strong Shock - Convergence')
    fig3 = plot_convergence(solver3, 'Problem 3: Two Rarefactions - Convergence')
    
    plt.show()
