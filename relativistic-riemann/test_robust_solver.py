#!/usr/bin/env python3
"""
Test script for robust iterative relativistic Riemann solver.
Compares with Kitajima solver and tests convergence improvements.
"""

import numpy as np
import matplotlib.pyplot as plt
from iterative_riemann_solver import IterativeRiemannSolver
from kitajima_solver import KitajimaRiemannSolver


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
    solver_iter = IterativeRiemannSolver(gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=True)
    solver_iter.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps_iter = solver_iter.solve_star_pressure()
    vs_iter = 0.5 * (solver_iter.vls + solver_iter.vrs)
    
    # Kitajima solver for comparison
    print("\n--- Kitajima Solver (Reference) ---")
    solver_kit = KitajimaRiemannSolver(gamma_c, c=1.0)
    solver_kit.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps_kit, vs_kit = solver_kit.solve_star_state()
    
    # Compare results
    print("\n--- Comparison ---")
    print(f"Star Pressure:")
    print(f"  Iterative: Ps = {Ps_iter:.10f}")
    print(f"  Kitajima:  Ps = {Ps_kit:.10f}")
    print(f"  Difference: ΔPs = {abs(Ps_iter - Ps_kit):.3e} ({abs(Ps_iter - Ps_kit)/Ps_kit*100:.3e}%)")
    
    print(f"Star Velocity:")
    print(f"  Iterative: vs = {vs_iter:.10f}")
    print(f"  Kitajima:  vs = {vs_kit:.10f}")
    print(f"  Difference: Δvs = {abs(vs_iter - vs_kit):.3e}")
    
    # Diagnostics
    diag = solver_iter.get_convergence_diagnostics()
    print(f"\nConvergence Diagnostics:")
    print(f"  Iterations: {diag['iterations']}")
    print(f"  Initial P*: {diag['initial_pressure']:.6e}")
    print(f"  Final P*:   {diag['final_pressure']:.6e}")
    print(f"  Final residual: {diag['final_residual']:.3e}")
    
    return solver_iter, solver_kit


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
    solver_iter = IterativeRiemannSolver(gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=True)
    solver_iter.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps_iter = solver_iter.solve_star_pressure()
    vs_iter = 0.5 * (solver_iter.vls + solver_iter.vrs)
    
    # Kitajima solver
    print("\n--- Kitajima Solver (Reference) ---")
    solver_kit = KitajimaRiemannSolver(gamma_c, c=1.0)
    solver_kit.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps_kit, vs_kit = solver_kit.solve_star_state()
    
    # Compare
    print("\n--- Comparison ---")
    print(f"Star Pressure:")
    print(f"  Iterative: Ps = {Ps_iter:.10f}")
    print(f"  Kitajima:  Ps = {Ps_kit:.10f}")
    print(f"  Difference: ΔPs = {abs(Ps_iter - Ps_kit):.3e} ({abs(Ps_iter - Ps_kit)/Ps_kit*100:.3e}%)")
    
    print(f"Star Velocity:")
    print(f"  Iterative: vs = {vs_iter:.10f}")
    print(f"  Kitajima:  vs = {vs_kit:.10f}")
    print(f"  Difference: Δvs = {abs(vs_iter - vs_kit):.3e}")
    
    return solver_iter, solver_kit


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
    solver_iter = IterativeRiemannSolver(gamma_c, c=1.0, max_iter=100, tol=1e-10, verbose=True)
    solver_iter.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps_iter = solver_iter.solve_star_pressure()
    vs_iter = 0.5 * (solver_iter.vls + solver_iter.vrs)
    
    # Kitajima solver
    print("\n--- Kitajima Solver (Reference) ---")
    solver_kit = KitajimaRiemannSolver(gamma_c, c=1.0)
    solver_kit.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    Ps_kit, vs_kit = solver_kit.solve_star_state()
    
    # Compare
    print("\n--- Comparison ---")
    print(f"Star Pressure:")
    print(f"  Iterative: Ps = {Ps_iter:.10f}")
    print(f"  Kitajima:  Ps = {Ps_kit:.10f}")
    print(f"  Difference: ΔPs = {abs(Ps_iter - Ps_kit):.3e} ({abs(Ps_iter - Ps_kit)/Ps_kit*100:.3e}%)")
    
    print(f"Star Velocity:")
    print(f"  Iterative: vs = {vs_iter:.10f}")
    print(f"  Kitajima:  vs = {vs_kit:.10f}")
    print(f"  Difference: Δvs = {abs(vs_iter - vs_kit):.3e}")
    
    return solver_iter, solver_kit


def plot_convergence_comparison(solvers, labels, title):
    """Plot convergence history for multiple solvers"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    for solver, label in zip(solvers, labels):
        iterations = range(len(solver.pressure_history))
        ax1.plot(iterations, solver.pressure_history, 'o-', label=label, markersize=4)
        
        if len(solver.residual_history) > 0:
            ax2.semilogy(range(1, len(solver.residual_history)+1), 
                        solver.residual_history, 'o-', label=label, markersize=4)
    
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Pressure P*')
    ax1.set_title('Pressure Convergence')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('|Residual| (velocity difference)')
    ax2.set_title('Residual Convergence')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.suptitle(title)
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    print("\n" + "="*60)
    print("ROBUST ITERATIVE RELATIVISTIC RIEMANN SOLVER TEST")
    print("="*60)
    
    # Run tests
    solver1_iter, solver1_kit = test_problem_1()
    solver2_iter, solver2_kit = test_problem_2()
    solver3_iter, solver3_kit = test_problem_3()
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Problem 1 (Moderate Shock):     {solver1_iter.iterations} iterations")
    print(f"Problem 2 (Strong Shock):       {solver2_iter.iterations} iterations")
    print(f"Problem 3 (Two Rarefactions):   {solver3_iter.iterations} iterations")
    
    # Plot convergence
    fig1 = plot_convergence_comparison(
        [solver1_iter], ['Robust Iterative'],
        'Problem 1: Moderate Shock - Convergence'
    )
    
    fig2 = plot_convergence_comparison(
        [solver2_iter], ['Robust Iterative'],
        'Problem 2: Strong Shock - Convergence'
    )
    
    fig3 = plot_convergence_comparison(
        [solver3_iter], ['Robust Iterative'],
        'Problem 3: Two Rarefactions - Convergence'
    )
    
    plt.show()
    
    print("\n✓ All tests completed successfully!")
