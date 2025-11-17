#!/usr/bin/env python3
"""
Animation script for Relativistic Riemann Solver

This script creates animated visualizations of the relativistic Riemann problem solution,
showing the evolution of density, pressure, velocity, and internal energy over time.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
from relativistic_riemann_solver import RelativisiticRiemannSolver


def create_animation(gamma, pl, rhol, vell, pr, rhor, velr, 
                     tmax=0.4, fps=30, duration=5.0,
                     output_file='relativistic_riemann.gif'):
    """
    Create an animation of the relativistic Riemann problem solution.
    
    Args:
        gamma: Adiabatic index
        pl, rhol, vell: Left state (pressure, density, velocity)
        pr, rhor, velr: Right state (pressure, density, velocity)
        tmax: Maximum time to simulate
        fps: Frames per second for animation
        duration: Total duration of animation in seconds
        output_file: Output filename for the animation
    """
    # Initialize solver
    solver = RelativisiticRiemannSolver(gamma)
    solver.set_initial_states(pl, rhol, vell, pr, rhor, velr)
    
    # Setup time array
    n_frames = int(fps * duration)
    times = np.linspace(0.001, tmax, n_frames)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0])  # Density
    ax2 = fig.add_subplot(gs[0, 1])  # Pressure
    ax3 = fig.add_subplot(gs[1, 0])  # Velocity
    ax4 = fig.add_subplot(gs[1, 1])  # Internal Energy
    
    # Compute solution for all times to get y-axis limits
    print("Computing solution range...")
    x_sample, p_sample, rho_sample, vel_sample, u_sample = solver.solve(tmax, n=400)
    
    # Set up axes
    axes = [ax1, ax2, ax3, ax4]
    titles = ['Density ρ', 'Pressure p', 'Velocity v', 'Internal Energy u']
    ylabels = ['ρ', 'p', 'v', 'u']
    data_samples = [rho_sample, p_sample, vel_sample, u_sample]
    
    lines = []
    for ax, title, ylabel, data in zip(axes, titles, ylabels, data_samples):
        ax.set_xlim(0, 1)
        ymin, ymax = data.min(), data.max()
        ymargin = 0.1 * (ymax - ymin) if ymax > ymin else 0.1
        ax.set_ylim(ymin - ymargin, ymax + ymargin)
        ax.set_xlabel('Position x')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0.5, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
        line, = ax.plot([], [], 'b-', linewidth=2)
        lines.append(line)
    
    # Add time text
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=14, weight='bold')
    
    # Add initial condition description
    ic_text = (f'γ={gamma:.2f} | Left: p={pl:.3f}, ρ={rhol:.3f}, v={vell:.3f} | '
               f'Right: p={pr:.3f}, ρ={rhor:.3f}, v={velr:.3f}')
    fig.text(0.5, 0.02, ic_text, ha='center', fontsize=10, style='italic')
    
    def init():
        """Initialize animation"""
        for line in lines:
            line.set_data([], [])
        time_text.set_text('')
        return lines + [time_text]
    
    def animate(frame):
        """Update animation frame"""
        t = times[frame]
        
        # Solve at current time
        x, pressure, density, velocity, internal_energy = solver.solve(t, n=400)
        
        # Update lines
        lines[0].set_data(x, density)
        lines[1].set_data(x, pressure)
        lines[2].set_data(x, velocity)
        lines[3].set_data(x, internal_energy)
        
        # Update time text
        time_text.set_text(f'Time t = {t:.4f}')
        
        if frame % 10 == 0:
            print(f'Frame {frame}/{n_frames} (t={t:.4f})')
        
        return lines + [time_text]
    
    # Create animation
    print(f"Creating animation with {n_frames} frames...")
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=n_frames, interval=1000/fps,
                                   blit=True, repeat=True)
    
    # Save animation
    print(f"Saving to {output_file}...")
    writer = animation.PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"Animation saved to {output_file}")
    
    plt.close()
    
    return anim


def create_static_plot(gamma, pl, rhol, vell, pr, rhor, velr,
                       times=[0.1, 0.2, 0.3, 0.4],
                       output_file='relativistic_riemann_static.png'):
    """
    Create a static multi-time plot of the solution.
    
    Args:
        gamma: Adiabatic index
        pl, rhol, vell: Left state (pressure, density, velocity)
        pr, rhor, velr: Right state (pressure, density, velocity)
        times: List of times to plot
        output_file: Output filename
    """
    # Initialize solver
    solver = RelativisiticRiemannSolver(gamma)
    solver.set_initial_states(pl, rhol, vell, pr, rhor, velr)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Relativistic Riemann Problem: γ={gamma:.2f}', fontsize=16, weight='bold')
    axes = axes.flatten()
    
    titles = ['Density ρ', 'Pressure p', 'Velocity v', 'Internal Energy u']
    ylabels = ['ρ', 'p', 'v', 'u']
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
    
    for t, color in zip(times, colors):
        x, pressure, density, velocity, internal_energy = solver.solve(t, n=400)
        data_arrays = [density, pressure, velocity, internal_energy]
        
        for ax, data, title, ylabel in zip(axes, data_arrays, titles, ylabels):
            ax.plot(x, data, label=f't={t:.2f}', color=color, linewidth=2)
            ax.set_xlabel('Position x')
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            ax.axvline(x=0.5, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
            ax.legend(loc='best')
    
    # Add initial condition description
    ic_text = (f'Initial Conditions: Left: p={pl:.3f}, ρ={rhol:.3f}, v={vell:.3f} | '
               f'Right: p={pr:.3f}, ρ={rhor:.3f}, v={velr:.3f}')
    fig.text(0.5, 0.02, ic_text, ha='center', fontsize=10, style='italic')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Static plot saved to {output_file}")
    plt.close()


def test_case_sr_sod():
    """Standard SR Sod shock tube test"""
    return {
        'name': 'SR Sod Shock Tube',
        'gamma': 1.4,
        'pl': 1.0,
        'rhol': 1.0,
        'vell': 0.0,
        'pr': 0.1,
        'rhor': 0.125,
        'velr': 0.0,
        'tmax': 0.4
    }


def test_case_blast_wave():
    """Strong blast wave test"""
    return {
        'name': 'Blast Wave',
        'gamma': 5.0/3.0,
        'pl': 1000.0,
        'rhol': 1.0,
        'vell': 0.0,
        'pr': 0.01,
        'rhor': 1.0,
        'velr': 0.0,
        'tmax': 0.4
    }


def test_case_relativistic_shock():
    """Relativistic shock test with high velocities"""
    return {
        'name': 'Relativistic Shock',
        'gamma': 4.0/3.0,
        'pl': 10.0,
        'rhol': 1.0,
        'vell': 0.9,  # 0.9c
        'pr': 1.0,
        'rhor': 1.0,
        'velr': 0.0,
        'tmax': 0.5
    }


def test_case_two_shocks():
    """Two shock collision"""
    return {
        'name': 'Two Shocks Collision',
        'gamma': 1.4,
        'pl': 1.0,
        'rhol': 1.0,
        'vell': 0.5,
        'pr': 1.0,
        'rhor': 1.0,
        'velr': -0.5,
        'tmax': 0.6
    }


def main():
    """Main function to generate animations"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Animate Relativistic Riemann Problems')
    parser.add_argument('--test', choices=['sod', 'blast', 'relativistic', 'twoshock', 'all'],
                       default='sod', help='Test case to run')
    parser.add_argument('--output', type=str, default=None, help='Output filename')
    parser.add_argument('--fps', type=int, default=30, help='Frames per second')
    parser.add_argument('--duration', type=float, default=5.0, help='Animation duration (seconds)')
    parser.add_argument('--static', action='store_true', help='Create static plot instead of animation')
    
    args = parser.parse_args()
    
    # Select test cases
    test_cases = {
        'sod': test_case_sr_sod(),
        'blast': test_case_blast_wave(),
        'relativistic': test_case_relativistic_shock(),
        'twoshock': test_case_two_shocks()
    }
    
    if args.test == 'all':
        cases_to_run = list(test_cases.keys())
    else:
        cases_to_run = [args.test]
    
    for case_name in cases_to_run:
        test = test_cases[case_name]
        print(f"\n{'='*60}")
        print(f"Running: {test['name']}")
        print(f"{'='*60}")
        
        if args.static:
            output = args.output or f"relativistic_riemann_{case_name}_static.png"
            create_static_plot(
                test['gamma'], test['pl'], test['rhol'], test['vell'],
                test['pr'], test['rhor'], test['velr'],
                times=np.linspace(0.1, test['tmax'], 4),
                output_file=output
            )
        else:
            output = args.output or f"relativistic_riemann_{case_name}.gif"
            create_animation(
                test['gamma'], test['pl'], test['rhol'], test['vell'],
                test['pr'], test['rhor'], test['velr'],
                tmax=test['tmax'],
                fps=args.fps,
                duration=args.duration,
                output_file=output
            )


if __name__ == "__main__":
    main()
