"""
Visualization and Animation Tools for SRGSPH

Tools to plot and animate SRGSPH results, comparing with exact Riemann solutions.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_problems import (setup_sod_problem, setup_standard_blast_wave,
                                  setup_strong_blast_wave, get_exact_solution)


def plot_solution(sim, exact_sol=None, filename=None, title="SRGSPH Solution"):
    """
    Plot SRGSPH solution with optional exact solution overlay.
    
    Args:
        sim: SRGSPH1D simulator object
        exact_sol: Dictionary with exact solution (optional)
        filename: Save plot to file (optional)
        title: Plot title
    """
    sol = sim.get_solution()
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"{title} at t = {sol['time']:.4f}", fontsize=16)
    
    # Sort by position for plotting
    idx = np.argsort(sol['x'])
    
    # Plot density in rest frame
    ax = axes[0, 0]
    ax.scatter(sol['x'][idx], sol['n'][idx], s=10, alpha=0.6, label='SRGSPH')
    if exact_sol is not None:
        ax.plot(exact_sol['x'], exact_sol['n'], 'k-', linewidth=2, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('n (rest frame density)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot pressure
    ax = axes[0, 1]
    ax.scatter(sol['x'][idx], sol['P'][idx], s=10, alpha=0.6, label='SRGSPH')
    if exact_sol is not None:
        ax.plot(exact_sol['x'], exact_sol['P'], 'k-', linewidth=2, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('P (pressure)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot velocity
    ax = axes[0, 2]
    ax.scatter(sol['x'][idx], sol['v'][idx], s=10, alpha=0.6, label='SRGSPH')
    if exact_sol is not None:
        ax.plot(exact_sol['x'], exact_sol['v'], 'k-', linewidth=2, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('v (velocity)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot lab frame density
    ax = axes[1, 0]
    ax.scatter(sol['x'][idx], sol['N'][idx], s=10, alpha=0.6, label='SRGSPH')
    if exact_sol is not None:
        ax.plot(exact_sol['x'], exact_sol['N'], 'k-', linewidth=2, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('N (lab frame density)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot Lorentz factor
    ax = axes[1, 1]
    ax.scatter(sol['x'][idx], sol['gamma'][idx], s=10, alpha=0.6, label='SRGSPH')
    if exact_sol is not None:
        ax.plot(exact_sol['x'], exact_sol['gamma'], 'k-', linewidth=2, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('γ (Lorentz factor)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot smoothing length
    ax = axes[1, 2]
    ax.scatter(sol['x'][idx], sol['h'][idx], s=10, alpha=0.6, label='h')
    ax.set_xlabel('x')
    ax.set_ylabel('h (smoothing length)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {filename}")
    else:
        plt.show()
    
    return fig


def animate_simulation(problem_type, t_end, n_frames=50, filename=None):
    """
    Create animation of SRGSPH simulation.
    
    Args:
        problem_type: 'sod', 'standard_blast', or 'strong_blast'
        t_end: End time
        n_frames: Number of frames
        filename: Save animation to file (optional, requires ffmpeg)
    """
    # Setup problem
    print(f"Setting up {problem_type} problem for animation...")
    
    if problem_type == 'sod':
        sim = setup_sod_problem()
    elif problem_type == 'standard_blast':
        sim = setup_standard_blast_wave()
    elif problem_type == 'strong_blast':
        sim = setup_strong_blast_wave()
    else:
        raise ValueError(f"Unknown problem: {problem_type}")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Storage for data
    frames_data = []
    times = []
    
    # Run simulation and collect frames
    dt_frame = t_end / n_frames
    print(f"Running simulation and collecting {n_frames} frames...")
    
    for i in range(n_frames + 1):
        t_target = i * dt_frame
        
        # Advance simulation
        while sim.time < t_target:
            dt = min(sim.compute_timestep(), t_target - sim.time)
            sim.step_forward(dt)
        
        # Collect data
        frames_data.append(sim.get_solution())
        times.append(sim.time)
        
        if i % 10 == 0:
            print(f"  Frame {i}/{n_frames}, t = {sim.time:.4f}")
    
    print("Creating animation...")
    
    # Get exact solution at final time
    exact_sol = get_exact_solution(problem_type, t_end, gamma_c=sim.gamma_c, c=sim.c)
    
    # Animation function
    def update(frame):
        sol = frames_data[frame]
        t = times[frame]
        
        for ax in axes.flat:
            ax.clear()
        
        idx = np.argsort(sol['x'])
        
        # Density
        axes[0, 0].scatter(sol['x'][idx], sol['n'][idx], s=10, alpha=0.6, c='C0')
        if frame == len(frames_data) - 1:
            axes[0, 0].plot(exact_sol['x'], exact_sol['n'], 'k-', linewidth=2, alpha=0.7)
        axes[0, 0].set_xlabel('x')
        axes[0, 0].set_ylabel('n (rest frame density)')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Pressure
        axes[0, 1].scatter(sol['x'][idx], sol['P'][idx], s=10, alpha=0.6, c='C1')
        if frame == len(frames_data) - 1:
            axes[0, 1].plot(exact_sol['x'], exact_sol['P'], 'k-', linewidth=2, alpha=0.7)
        axes[0, 1].set_xlabel('x')
        axes[0, 1].set_ylabel('P (pressure)')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Velocity
        axes[1, 0].scatter(sol['x'][idx], sol['v'][idx], s=10, alpha=0.6, c='C2')
        if frame == len(frames_data) - 1:
            axes[1, 0].plot(exact_sol['x'], exact_sol['v'], 'k-', linewidth=2, alpha=0.7)
        axes[1, 0].set_xlabel('x')
        axes[1, 0].set_ylabel('v (velocity)')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Lab frame density
        axes[1, 1].scatter(sol['x'][idx], sol['N'][idx], s=10, alpha=0.6, c='C3')
        if frame == len(frames_data) - 1:
            axes[1, 1].plot(exact_sol['x'], exact_sol['N'], 'k-', linewidth=2, alpha=0.7)
        axes[1, 1].set_xlabel('x')
        axes[1, 1].set_ylabel('N (lab frame density)')
        axes[1, 1].grid(True, alpha=0.3)
        
        fig.suptitle(f'SRGSPH: {problem_type} at t = {t:.4f}', fontsize=14)
        
        return axes.flat
    
    anim = FuncAnimation(fig, update, frames=len(frames_data), interval=100, blit=False)
    
    if filename:
        print(f"Saving animation to {filename}...")
        anim.save(filename, writer='pillow', fps=10, dpi=100)
        print(f"Animation saved to {filename}")
    else:
        plt.show()
    
    return anim


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize SRGSPH results')
    parser.add_argument('problem', choices=['sod', 'standard_blast', 'strong_blast'],
                        help='Test problem')
    parser.add_argument('--mode', choices=['plot', 'animate'], default='plot',
                        help='Visualization mode')
    parser.add_argument('--t_end', type=float, default=None,
                        help='End time')
    parser.add_argument('--output', type=str, default=None,
                        help='Output filename')
    parser.add_argument('--frames', type=int, default=50,
                        help='Number of frames for animation')
    
    args = parser.parse_args()
    
    # Determine end time
    if args.t_end is None:
        if args.problem == 'sod':
            t_end = 0.35
        elif args.problem == 'standard_blast':
            t_end = 0.4
        else:  # strong_blast
            t_end = 0.16
    else:
        t_end = args.t_end
    
    if args.mode == 'animate':
        # Create animation
        output = args.output if args.output else f'srgsph_{args.problem}.gif'
        animate_simulation(args.problem, t_end, n_frames=args.frames, filename=output)
    else:
        # Run simulation and plot
        print(f"Running {args.problem} problem to t = {t_end}...")
        
        if args.problem == 'sod':
            sim = setup_sod_problem()
        elif args.problem == 'standard_blast':
            sim = setup_standard_blast_wave()
        else:
            sim = setup_strong_blast_wave()
        
        sim.run(t_end)
        
        # Get exact solution
        exact_sol = get_exact_solution(args.problem, t_end, gamma_c=sim.gamma_c, c=sim.c)
        
        # Plot
        output = args.output if args.output else f'srgsph_{args.problem}.png'
        plot_solution(sim, exact_sol=exact_sol, filename=output, 
                     title=f"SRGSPH: {args.problem.replace('_', ' ').title()}")
