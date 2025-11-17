#!/usr/bin/env python3
"""
Animation script for Kitajima formulation of relativistic Riemann solver.

Shows evolution of baryon number density N, pressure P, velocity v,
Lorentz factor γ, canonical momentum S, and canonical energy e.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
from kitajima_solver import KitajimaRiemannSolver


def create_animation(gamma_c, Pl, nl, vl, Pr, nr, vr, c=1.0,
                     tmax=0.4, fps=30, duration=5.0,
                     output_file='kitajima_riemann.gif'):
    """
    Create animation of Kitajima formulation solution.
    
    Args:
        gamma_c: Adiabatic index
        Pl, nl, vl: Left state (pressure, rest frame density, velocity)
        Pr, nr, vr: Right state
        c: Speed of light
        tmax: Maximum time
        fps: Frames per second
        duration: Animation duration (seconds)
        output_file: Output filename
    """
    solver = KitajimaRiemannSolver(gamma_c, c)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    
    n_frames = int(fps * duration)
    times = np.linspace(0.001, tmax, n_frames)
    
    # Create figure with 6 panels
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0])  # Lab frame baryon density N
    ax2 = fig.add_subplot(gs[0, 1])  # Pressure P
    ax3 = fig.add_subplot(gs[1, 0])  # Velocity v
    ax4 = fig.add_subplot(gs[1, 1])  # Lorentz factor γ
    ax5 = fig.add_subplot(gs[2, 0])  # Canonical momentum S
    ax6 = fig.add_subplot(gs[2, 1])  # Canonical energy e
    
    # Sample solution for axis limits
    print("Computing solution range...")
    x_s, P_s, n_s, N_s, v_s, u_s, g_s, S_s, e_s = solver.solve(tmax, n_points=400)
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    titles = [
        'Lab Frame Baryon Number Density N = γn',
        'Pressure P',
        'Velocity v',
        'Lorentz Factor γ',
        'Canonical Momentum S = γHv',
        'Canonical Energy e = γH - P/(Nc²)'
    ]
    ylabels = ['N', 'P', 'v/c', 'γ', 'S', 'e']
    samples = [N_s, P_s, v_s/c, g_s, S_s, e_s]
    
    lines = []
    for ax, title, ylabel, data in zip(axes, titles, ylabels, samples):
        ax.set_xlim(0, 1)
        ymin, ymax = data.min(), data.max()
        ymargin = 0.1 * (ymax - ymin) if ymax > ymin else 0.1 * abs(ymax) if ymax != 0 else 0.1
        ax.set_ylim(ymin - ymargin, ymax + ymargin)
        ax.set_xlabel('Position x')
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0.5, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
        line, = ax.plot([], [], 'b-', linewidth=2)
        lines.append(line)
    
    # Title with time and parameters
    time_text = fig.text(0.5, 0.96, '', ha='center', fontsize=14, weight='bold')
    
    # Initial conditions
    ic_text = (f'Kitajima Formulation (Baryon Number Density) | '
               f'γc={gamma_c:.2f}, c={c:.1f} | '
               f'Left: P={Pl:.3f}, n={nl:.3f}, v={vl:.3f} | '
               f'Right: P={Pr:.3f}, n={nr:.3f}, v={vr:.3f}')
    fig.text(0.5, 0.02, ic_text, ha='center', fontsize=9, style='italic')
    
    def init():
        for line in lines:
            line.set_data([], [])
        time_text.set_text('')
        return lines + [time_text]
    
    def animate(frame):
        t = times[frame]
        
        x, P, n, N, v, u, gamma, S, e = solver.solve(t, n_points=400)
        
        lines[0].set_data(x, N)
        lines[1].set_data(x, P)
        lines[2].set_data(x, v/c)
        lines[3].set_data(x, gamma)
        lines[4].set_data(x, S)
        lines[5].set_data(x, e)
        
        time_text.set_text(f'Time t = {t:.4f}')
        
        if frame % 10 == 0:
            print(f'Frame {frame}/{n_frames} (t={t:.4f})')
        
        return lines + [time_text]
    
    print(f"Creating animation with {n_frames} frames...")
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=n_frames, interval=1000/fps,
                                   blit=True, repeat=True)
    
    print(f"Saving to {output_file}...")
    writer = animation.PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"Animation saved!")
    
    plt.close()
    return anim


def create_static_plot(gamma_c, Pl, nl, vl, Pr, nr, vr, c=1.0,
                       times=[0.1, 0.2, 0.3, 0.4],
                       output_file='kitajima_static.png'):
    """Create static multi-time comparison plot."""
    solver = KitajimaRiemannSolver(gamma_c, c)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr)
    
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    axes = [
        fig.add_subplot(gs[0, 0]),  # N
        fig.add_subplot(gs[0, 1]),  # P
        fig.add_subplot(gs[1, 0]),  # v
        fig.add_subplot(gs[1, 1]),  # γ
        fig.add_subplot(gs[2, 0]),  # S
        fig.add_subplot(gs[2, 1]),  # e
    ]
    
    titles = ['Lab Density N', 'Pressure P', 'Velocity v/c', 'Lorentz Factor γ',
              'Canonical Momentum S', 'Canonical Energy e']
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
    
    for t, color in zip(times, colors):
        x, P, n, N, v, u, gamma, S, e = solver.solve(t, n_points=400)
        data = [N, P, v/c, gamma, S, e]
        
        for ax, d, title in zip(axes, data, titles):
            ax.plot(x, d, label=f't={t:.2f}', color=color, linewidth=2)
            ax.set_xlabel('Position x')
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            ax.axvline(x=0.5, color='k', linestyle='--', alpha=0.3)
            ax.legend()
    
    fig.suptitle(f'Kitajima Formulation: γc={gamma_c:.2f}, c={c:.1f}',
                 fontsize=14, weight='bold')
    
    ic_text = (f'Left: P={Pl:.3f}, n={nl:.3f}, v={vl:.3f} | '
               f'Right: P={Pr:.3f}, n={nr:.3f}, v={vr:.3f}')
    fig.text(0.5, 0.02, ic_text, ha='center', fontsize=10, style='italic')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Static plot saved to {output_file}")
    plt.close()


def test_case_sod_kitajima(c=1.0):
    """SR Sod shock tube in Kitajima formulation"""
    return {
        'name': 'SR Sod (Kitajima)',
        'gamma_c': 5.0/3.0,
        'Pl': 1.0,
        'nl': 1.0,
        'vl': 0.0,
        'Pr': 0.1,
        'nr': 0.125,
        'vr': 0.0,
        'c': c,
        'tmax': 0.35
    }


def test_case_blast_kitajima(c=1.0):
    """Blast wave in Kitajima formulation"""
    return {
        'name': 'Blast Wave (Kitajima)',
        'gamma_c': 5.0/3.0,
        'Pl': 1000.0,
        'nl': 1.0,
        'vl': 0.0,
        'Pr': 0.01,
        'nr': 1.0,
        'vr': 0.0,
        'c': c,
        'tmax': 0.4
    }


def test_case_relativistic_kitajima(c=1.0):
    """Ultra-relativistic shock with v=0.9c"""
    return {
        'name': 'Relativistic Shock (Kitajima)',
        'gamma_c': 5.0/3.0,
        'Pl': 1.0,
        'nl': 1.0,
        'vl': 0.9 * c,  # 0.9c
        'Pr': 1.0,
        'nr': 1.0,
        'vr': 0.0,
        'c': c,
        'tmax': 0.3
    }


def test_case_ultra_relativistic_kitajima(c=1.0):
    """Ultra-relativistic with v=0.999c"""
    return {
        'name': 'Ultra-Relativistic (Kitajima)',
        'gamma_c': 5.0/3.0,
        'Pl': 1.0,
        'nl': 1.0,
        'vl': 0.999 * c,  # 0.999c
        'Pr': 1.0,
        'nr': 1.0,
        'vr': 0.0,
        'c': c,
        'tmax': 0.3
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Animate Kitajima formulation relativistic Riemann problems')
    parser.add_argument('--test', choices=['sod', 'blast', 'relativistic', 'ultra', 'all'],
                       default='sod', help='Test case')
    parser.add_argument('--c', type=float, default=1.0,
                       help='Speed of light (default=1.0 for natural units)')
    parser.add_argument('--output', type=str, default=None, help='Output file')
    parser.add_argument('--fps', type=int, default=30, help='Frames per second')
    parser.add_argument('--duration', type=float, default=5.0, help='Duration (seconds)')
    parser.add_argument('--static', action='store_true', help='Create static plot')
    
    args = parser.parse_args()
    
    test_cases = {
        'sod': test_case_sod_kitajima(args.c),
        'blast': test_case_blast_kitajima(args.c),
        'relativistic': test_case_relativistic_kitajima(args.c),
        'ultra': test_case_ultra_relativistic_kitajima(args.c),
    }
    
    if args.test == 'all':
        cases = list(test_cases.keys())
    else:
        cases = [args.test]
    
    for case_name in cases:
        test = test_cases[case_name]
        print(f"\n{'='*60}")
        print(f"Running: {test['name']}")
        print(f"Speed of light c = {test['c']:.1f}")
        print(f"{'='*60}")
        
        if args.static:
            output = args.output or f"kitajima_{case_name}_static.png"
            create_static_plot(
                test['gamma_c'], test['Pl'], test['nl'], test['vl'],
                test['Pr'], test['nr'], test['vr'], test['c'],
                times=np.linspace(0.05, test['tmax'], 4),
                output_file=output
            )
        else:
            output = args.output or f"kitajima_{case_name}.gif"
            create_animation(
                test['gamma_c'], test['Pl'], test['nl'], test['vl'],
                test['Pr'], test['nr'], test['vr'], test['c'],
                tmax=test['tmax'],
                fps=args.fps,
                duration=args.duration,
                output_file=output
            )


if __name__ == "__main__":
    main()
