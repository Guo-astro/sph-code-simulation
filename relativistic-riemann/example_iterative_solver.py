#!/usr/bin/env python3
"""
Example: Iterative Relativistic Riemann Solver

This script demonstrates the iterative relativistic Riemann solver
and creates animations comparing it with the exact Kitajima solver.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from iterative_riemann_solver import IterativeRiemannSolver
from kitajima_solver import KitajimaRiemannSolver


def create_comparison_animation(solver_exact, solver_iter, 
                                times, output_file='iterative_solver_animation.gif'):
    """
    Create animation comparing exact and iterative solvers over time.
    
    Args:
        solver_exact: KitajimaRiemannSolver instance
        solver_iter: IterativeRiemannSolver instance
        times: Array of time points
        output_file: Output animation filename
    """
    print("\n" + "="*70)
    print("Creating comparison animation...")
    print("="*70)
    
    # Setup figure with 3 rows, 2 columns
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle('Iterative vs Exact Relativistic Riemann Solver', 
                 fontsize=14, fontweight='bold')
    
    # Column titles
    axes[0, 0].set_title('Exact (Kitajima)', fontsize=12, fontweight='bold')
    axes[0, 1].set_title('Iterative (Newton-Raphson)', fontsize=12, fontweight='bold')
    
    # Initialize line plots
    lines_exact = []
    lines_iter = []
    
    # Density plots
    line_e, = axes[0, 0].plot([], [], 'b-', linewidth=2, label='Exact')
    line_i, = axes[0, 1].plot([], [], 'r-', linewidth=2, label='Iterative')
    lines_exact.append(line_e)
    lines_iter.append(line_i)
    axes[0, 0].set_ylabel('Density (n)', fontsize=10)
    axes[0, 0].set_ylim(0, 1.2)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 1].set_ylabel('Density (n)', fontsize=10)
    axes[0, 1].set_ylim(0, 1.2)
    axes[0, 1].grid(True, alpha=0.3)
    
    # Pressure plots
    line_e, = axes[1, 0].plot([], [], 'g-', linewidth=2)
    line_i, = axes[1, 1].plot([], [], 'r-', linewidth=2)
    lines_exact.append(line_e)
    lines_iter.append(line_i)
    axes[1, 0].set_ylabel('Pressure (P)', fontsize=10)
    axes[1, 0].set_ylim(0, 1.1)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 1].set_ylabel('Pressure (P)', fontsize=10)
    axes[1, 1].set_ylim(0, 1.1)
    axes[1, 1].grid(True, alpha=0.3)
    
    # Velocity plots
    line_e, = axes[2, 0].plot([], [], 'm-', linewidth=2)
    line_i, = axes[2, 1].plot([], [], 'r-', linewidth=2)
    lines_exact.append(line_e)
    lines_iter.append(line_i)
    axes[2, 0].set_ylabel('Velocity (v)', fontsize=10)
    axes[2, 0].set_xlabel('Position (x)', fontsize=10)
    axes[2, 0].set_ylim(-0.1, 1.0)
    axes[2, 0].grid(True, alpha=0.3)
    axes[2, 1].set_ylabel('Velocity (v)', fontsize=10)
    axes[2, 1].set_xlabel('Position (x)', fontsize=10)
    axes[2, 1].set_ylim(-0.1, 1.0)
    axes[2, 1].grid(True, alpha=0.3)
    
    # Set x-axis limits for all
    for ax in axes.flat:
        ax.set_xlim(0, 1)
        ax.axvline(x=0.5, color='k', linestyle=':', alpha=0.5, linewidth=1)
    
    # Time text
    time_text = fig.text(0.5, 0.96, '', ha='center', fontsize=12, fontweight='bold')
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    def animate(frame):
        """Animation update function"""
        t = times[frame]
        
        # Solve with both methods
        x_exact, P_exact, n_exact, N_exact, v_exact, u_exact, gamma_exact, S_exact, e_exact = \
            solver_exact.solve(t, x0=0.5, n_points=400)
        
        x_iter, P_iter, n_iter, N_iter, v_iter, u_iter, gamma_iter, S_iter, e_iter = \
            solver_iter.solve(t, x0=0.5, n_points=400)
        
        # Update exact solver plots
        lines_exact[0].set_data(x_exact, n_exact)
        lines_exact[1].set_data(x_exact, P_exact)
        lines_exact[2].set_data(x_exact, v_exact)
        
        # Update iterative solver plots
        lines_iter[0].set_data(x_iter, n_iter)
        lines_iter[1].set_data(x_iter, P_iter)
        lines_iter[2].set_data(x_iter, v_iter)
        
        # Update time text
        time_text.set_text(f't = {t:.3f}  |  Iterations: {solver_iter.iterations}')
        
        return lines_exact + lines_iter + [time_text]
    
    # Create animation
    n_frames = len(times)
    print(f"Generating {n_frames} frames...")
    anim = animation.FuncAnimation(
        fig, animate, frames=n_frames,
        interval=100, blit=True, repeat=True
    )
    
    # Save animation
    print(f"Saving to {output_file}...")
    writer = animation.PillowWriter(fps=10)
    anim.save(output_file, writer=writer, dpi=100)
    
    file_size_mb = Path(output_file).stat().st_size / (1024 * 1024)
    print(f"✓ Animation saved: {output_file}")
    print(f"  File size: {file_size_mb:.2f} MB")
    print(f"  Frames: {n_frames}")
    print(f"  Time range: t = {times[0]:.3f} to {times[-1]:.3f}")
    
    plt.close(fig)


def create_convergence_plot(solver_iter, output_file='iterative_convergence.png'):
    """
    Plot convergence history of iterative solver.
    
    Args:
        solver_iter: IterativeRiemannSolver with convergence history
        output_file: Output filename
    """
    if len(solver_iter.pressure_history) == 0:
        print("No convergence history available")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle('Iterative Solver Convergence', fontsize=14, fontweight='bold')
    
    iterations = np.arange(len(solver_iter.pressure_history))
    
    # Pressure evolution
    ax1.plot(iterations, solver_iter.pressure_history, 'b-o', linewidth=2, markersize=6)
    ax1.set_ylabel('Star Pressure (Ps)', fontsize=11)
    ax1.set_xlabel('Iteration', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=solver_iter.pressure_history[-1], color='r', 
                linestyle='--', alpha=0.5, label='Final value')
    ax1.legend()
    
    # Residual evolution (log scale)
    if len(solver_iter.residual_history) > 0:
        ax2.semilogy(iterations[1:], solver_iter.residual_history, 'r-o', 
                     linewidth=2, markersize=6)
        ax2.set_ylabel('|Residual| (log scale)', fontsize=11)
        ax2.set_xlabel('Iteration', fontsize=11)
        ax2.grid(True, alpha=0.3, which='both')
        ax2.axhline(y=solver_iter.tol, color='g', linestyle='--', 
                    alpha=0.5, label=f'Tolerance = {solver_iter.tol:.1e}')
        ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Convergence plot saved: {output_file}")
    plt.close(fig)


def create_error_analysis(solver_exact, solver_iter, t, output_file='error_analysis.png'):
    """
    Create detailed error analysis between exact and iterative solutions.
    
    Args:
        solver_exact: KitajimaRiemannSolver instance
        solver_iter: IterativeRiemannSolver instance
        t: Time point
        output_file: Output filename
    """
    print(f"\nError analysis at t = {t}...")
    
    # Solve with both methods
    x_exact, P_exact, n_exact, N_exact, v_exact, u_exact, gamma_exact, S_exact, e_exact = \
        solver_exact.solve(t, x0=0.5, n_points=400)
    
    x_iter, P_iter, n_iter, N_iter, v_iter, u_iter, gamma_iter, S_iter, e_iter = \
        solver_iter.solve(t, x0=0.5, n_points=400)
    
    # Compute errors
    n_err = np.abs(n_iter - n_exact)
    P_err = np.abs(P_iter - P_exact)
    v_err = np.abs(v_iter - v_exact)
    
    # L1 and L2 norms
    n_l1 = np.mean(n_err)
    n_l2 = np.sqrt(np.mean(n_err**2))
    P_l1 = np.mean(P_err)
    P_l2 = np.sqrt(np.mean(P_err**2))
    v_l1 = np.mean(v_err)
    v_l2 = np.sqrt(np.mean(v_err**2))
    
    print(f"  Density  - L1 error: {n_l1:.3e}, L2 error: {n_l2:.3e}")
    print(f"  Pressure - L1 error: {P_l1:.3e}, L2 error: {P_l2:.3e}")
    print(f"  Velocity - L1 error: {v_l1:.3e}, L2 error: {v_l2:.3e}")
    
    # Create comparison plot
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle(f'Iterative vs Exact: Error Analysis at t = {t:.3f}', 
                 fontsize=14, fontweight='bold')
    
    # Density
    axes[0, 0].plot(x_exact, n_exact, 'b-', linewidth=2, label='Exact')
    axes[0, 0].plot(x_iter, n_iter, 'r--', linewidth=2, label='Iterative')
    axes[0, 0].set_ylabel('Density (n)', fontsize=10)
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(x_exact, n_err, 'k-', linewidth=2)
    axes[0, 1].set_ylabel('Absolute Error', fontsize=10)
    axes[0, 1].set_title(f'L1={n_l1:.2e}, L2={n_l2:.2e}', fontsize=10)
    axes[0, 1].grid(True, alpha=0.3)
    
    # Pressure
    axes[1, 0].plot(x_exact, P_exact, 'g-', linewidth=2, label='Exact')
    axes[1, 0].plot(x_iter, P_iter, 'r--', linewidth=2, label='Iterative')
    axes[1, 0].set_ylabel('Pressure (P)', fontsize=10)
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].plot(x_exact, P_err, 'k-', linewidth=2)
    axes[1, 1].set_ylabel('Absolute Error', fontsize=10)
    axes[1, 1].set_title(f'L1={P_l1:.2e}, L2={P_l2:.2e}', fontsize=10)
    axes[1, 1].grid(True, alpha=0.3)
    
    # Velocity
    axes[2, 0].plot(x_exact, v_exact, 'm-', linewidth=2, label='Exact')
    axes[2, 0].plot(x_iter, v_iter, 'r--', linewidth=2, label='Iterative')
    axes[2, 0].set_ylabel('Velocity (v)', fontsize=10)
    axes[2, 0].set_xlabel('Position (x)', fontsize=10)
    axes[2, 0].legend()
    axes[2, 0].grid(True, alpha=0.3)
    
    axes[2, 1].plot(x_exact, v_err, 'k-', linewidth=2)
    axes[2, 1].set_ylabel('Absolute Error', fontsize=10)
    axes[2, 1].set_xlabel('Position (x)', fontsize=10)
    axes[2, 1].set_title(f'L1={v_l1:.2e}, L2={v_l2:.2e}', fontsize=10)
    axes[2, 1].grid(True, alpha=0.3)
    
    for ax in axes.flat:
        ax.set_xlim(0, 1)
        ax.axvline(x=0.5, color='k', linestyle=':', alpha=0.5, linewidth=1)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Error analysis saved: {output_file}")
    plt.close(fig)


def main():
    """Run iterative solver example with animations."""
    
    print("\n" + "="*70)
    print("ITERATIVE RELATIVISTIC RIEMANN SOLVER EXAMPLE")
    print("="*70)
    print()
    
    # Problem parameters
    gamma_c = 5.0 / 3.0  # Adiabatic index
    c = 1.0              # Speed of light
    
    # Initial conditions: Classic Sod problem
    P_left = 1.0
    n_left = 1.0
    v_left = 0.0
    
    P_right = 0.1
    n_right = 0.125
    v_right = 0.0
    
    print("Problem Setup:")
    print(f"  Adiabatic index: γ = {gamma_c:.4f}")
    print(f"  Speed of light: c = {c}")
    print(f"  Left state:  P = {P_left}, n = {n_left}, v = {v_left}")
    print(f"  Right state: P = {P_right}, n = {n_right}, v = {v_right}")
    print()
    
    # Initialize solvers
    print("Initializing solvers...")
    solver_exact = KitajimaRiemannSolver(gamma_c, c)
    solver_exact.set_initial_states(P_left, n_left, v_left, P_right, n_right, v_right)
    
    solver_iter = IterativeRiemannSolver(gamma_c, c, max_iter=100, tol=1e-10)
    solver_iter.set_initial_states(P_left, n_left, v_left, P_right, n_right, v_right)
    print("✓ Solvers initialized")
    print()
    
    # Create output directory
    output_dir = Path(__file__).parent / 'iterative_solver_results'
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir}")
    print()
    
    # Test single solution
    t_test = 0.35
    print(f"Testing at t = {t_test}...")
    print("-" * 70)
    
    print("\nExact solver:")
    x_exact, P_exact, n_exact, N_exact, v_exact, u_exact, gamma_exact, S_exact, e_exact = \
        solver_exact.solve(t_test, x0=0.5, n_points=400)
    
    print("\nIterative solver:")
    x_iter, P_iter, n_iter, N_iter, v_iter, u_iter, gamma_iter, S_iter, e_iter = \
        solver_iter.solve(t_test, x0=0.5, n_points=400)
    
    print("-" * 70)
    print()
    
    # Create convergence plot
    create_convergence_plot(solver_iter, 
                           output_dir / 'iterative_convergence.png')
    
    # Create error analysis
    create_error_analysis(solver_exact, solver_iter, t_test,
                         output_dir / 'error_analysis.png')
    
    # Create animation over time
    times = np.linspace(0.05, 0.40, 30)  # 30 frames
    create_comparison_animation(solver_exact, solver_iter, times,
                               output_dir / 'iterative_solver_animation.gif')
    
    # Save final solution data
    solver_iter.save_solution(
        output_dir / f'iterative_solution_t{t_test:.3f}.dat',
        x_iter, P_iter, n_iter, N_iter, v_iter, u_iter, gamma_iter, S_iter, e_iter, t_test
    )
    
    print()
    print("="*70)
    print("EXAMPLE COMPLETE!")
    print("="*70)
    print()
    print("Generated files:")
    print(f"  1. {output_dir / 'iterative_convergence.png'}")
    print(f"  2. {output_dir / 'error_analysis.png'}")
    print(f"  3. {output_dir / 'iterative_solver_animation.gif'}")
    print(f"  4. {output_dir / f'iterative_solution_t{t_test:.3f}.dat'}")
    print()
    print("Summary:")
    print(f"  The iterative solver converged in {solver_iter.iterations} iterations")
    print(f"  Final residual: {solver_iter.residual_history[-1]:.3e}")
    print(f"  The animation shows time evolution from t=0.05 to t=0.40")
    print()


if __name__ == '__main__':
    main()
