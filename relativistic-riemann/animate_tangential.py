#!/usr/bin/env python3
"""
Animate tangential velocity effects in relativistic Riemann solver.

Creates animations showing:
1. Comparison of solutions with different tangential velocities
2. Evolution of tangential velocity components
3. Conservation of h*W*v^t across waves
4. Lorentz factor evolution

Based on Pons, Martí & Müller (1999, J. Fluid Mech.)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from kitajima_solver import KitajimaRiemannSolver


def create_comparison_animation():
    """
    Create animation comparing solutions with different tangential velocities.
    """
    print("Creating tangential velocity comparison animation...")
    
    gamma = 1.4
    c = 1.0
    
    # Initial conditions - Sod shock tube
    Pl, nl, vl = 1.0, 1.0, 0.0
    Pr, nr, vr = 0.1, 0.125, 0.0
    
    # Different tangential velocities
    tangential_vels = [0.0, 0.5, 0.9, 0.99]
    colors = ['blue', 'green', 'orange', 'red']
    labels = [f'$v^y = {vy:.2f}$' for vy in tangential_vels]
    
    # Time steps for animation
    t_max = 0.25
    n_frames = 100
    times = np.linspace(0.001, t_max, n_frames)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    ax_p = fig.add_subplot(gs[0, 0])      # Pressure
    ax_rho = fig.add_subplot(gs[0, 1])    # Density
    ax_vx = fig.add_subplot(gs[0, 2])     # Normal velocity
    ax_vy = fig.add_subplot(gs[1, 0])     # Tangential velocity
    ax_gamma = fig.add_subplot(gs[1, 1])  # Lorentz factor
    ax_cons = fig.add_subplot(gs[1, 2])   # Conservation h*W*v^y
    ax_vtot = fig.add_subplot(gs[2, 0])   # Total velocity
    ax_H = fig.add_subplot(gs[2, 1])      # Enthalpy
    ax_e = fig.add_subplot(gs[2, 2])      # Energy per baryon
    
    axes = [ax_p, ax_rho, ax_vx, ax_vy, ax_gamma, ax_cons, ax_vtot, ax_H, ax_e]
    titles = ['Pressure', 'Density', 'Normal Velocity $v^x$', 'Tangential Velocity $v^y$',
              'Lorentz Factor $\\gamma$', 'Conservation $hW v^y$', 
              'Total Velocity $|v|$', 'Enthalpy $h$', 'Energy/Baryon $e$']
    
    for ax, title in zip(axes, titles):
        ax.set_title(title, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel('x', fontsize=9)
    
    # Initialize line collections
    lines = {ax: [] for ax in axes}
    
    def init():
        for ax in axes:
            ax.clear()
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, 1)
        return []
    
    def animate(frame):
        t = times[frame]
        
        # Clear previous lines
        for ax in axes:
            ax.clear()
            ax.grid(True, alpha=0.3)
            ax.set_xlabel('x', fontsize=9)
        
        # Solve for each tangential velocity
        for i, vy in enumerate(tangential_vels):
            solver = KitajimaRiemannSolver(gamma, c)
            solver.set_initial_states(Pl, nl, vl, Pr, nr, vr,
                                     vyl=vy, vzl=0.0, vyr=vy, vzr=0.0)
            
            x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(
                t, x0=0.5, n_points=400)
            
            # Compute derived quantities
            H = 1.0 + u/(c**2) + P/(n * c**2)
            vtot = np.sqrt(v**2 + vy_arr**2 + vz_arr**2)
            conserved = H * gamma_arr * vy_arr
            
            color = colors[i]
            label = labels[i]
            
            # Plot all quantities
            ax_p.plot(x, P, color=color, label=label, linewidth=2)
            ax_rho.plot(x, n, color=color, label=label, linewidth=2)
            ax_vx.plot(x, v, color=color, label=label, linewidth=2)
            ax_vy.plot(x, vy_arr, color=color, label=label, linewidth=2)
            ax_gamma.plot(x, gamma_arr, color=color, label=label, linewidth=2)
            ax_cons.plot(x, conserved, color=color, label=label, linewidth=2)
            ax_vtot.plot(x, vtot, color=color, label=label, linewidth=2)
            ax_H.plot(x, H, color=color, label=label, linewidth=2)
            ax_e.plot(x, e, color=color, label=label, linewidth=2)
        
        # Set titles and limits
        ax_p.set_title(f'Pressure (t={t:.3f})', fontsize=10)
        ax_p.set_ylim(0, 1.2)
        ax_p.legend(fontsize=8, loc='upper right')
        
        ax_rho.set_title(f'Density (t={t:.3f})', fontsize=10)
        ax_rho.set_ylim(0, 1.2)
        
        ax_vx.set_title(f'Normal Velocity $v^x$ (t={t:.3f})', fontsize=10)
        ax_vx.set_ylim(-0.1, 1.0)
        
        ax_vy.set_title(f'Tangential Velocity $v^y$ (t={t:.3f})', fontsize=10)
        ax_vy.set_ylim(-0.1, 1.2)
        
        ax_gamma.set_title(f'Lorentz Factor $\\gamma$ (t={t:.3f})', fontsize=10)
        ax_gamma.set_ylim(0.9, 8)
        ax_gamma.set_yscale('log')
        
        ax_cons.set_title(f'Conservation $hW v^y$ (t={t:.3f})', fontsize=10)
        
        ax_vtot.set_title(f'Total Velocity $|v|$ (t={t:.3f})', fontsize=10)
        ax_vtot.set_ylim(0, 1.1)
        ax_vtot.axhline(y=c, color='black', linestyle='--', alpha=0.5, label='c')
        
        ax_H.set_title(f'Enthalpy $h$ (t={t:.3f})', fontsize=10)
        
        ax_e.set_title(f'Energy/Baryon $e$ (t={t:.3f})', fontsize=10)
        
        fig.suptitle('Tangential Velocity Effects in Relativistic Riemann Problem\n' +
                    'Pons, Martí & Müller (1999) - Kitajima Formulation',
                    fontsize=12, fontweight='bold')
        
        return []
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=n_frames, interval=50, blit=False)
    
    # Save animation
    output_file = 'tangential_velocity_comparison.mp4'
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer='ffmpeg', fps=20, dpi=150)
    print(f"✓ Animation saved: {output_file}")
    
    plt.close()


def create_wave_structure_animation():
    """
    Create animation showing wave structure with tangential velocity.
    """
    print("\nCreating wave structure animation...")
    
    gamma = 1.4
    c = 1.0
    
    # Initial conditions with significant tangential velocity
    Pl, nl, vl = 1.0, 1.0, 0.0
    Pr, nr, vr = 0.1, 0.125, 0.0
    vyl, vyr = 0.9, 0.9
    
    solver = KitajimaRiemannSolver(gamma, c)
    solver.set_initial_states(Pl, nl, vl, Pr, nr, vr,
                             vyl=vyl, vzl=0.0, vyr=vyr, vzr=0.0)
    
    # Time steps
    t_max = 0.25
    n_frames = 100
    times = np.linspace(0.001, t_max, n_frames)
    
    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    
    def animate(frame):
        t = times[frame]
        x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(
            t, x0=0.5, n_points=400)
        
        # Compute derived quantities
        H = 1.0 + u/(c**2) + P/(n * c**2)
        vtot = np.sqrt(v**2 + vy_arr**2)
        conserved_y = H * gamma_arr * vy_arr
        
        for ax in axes:
            ax.clear()
            ax.grid(True, alpha=0.3)
            ax.set_xlabel('x')
        
        # Plot different quantities
        axes[0].plot(x, P, 'b-', linewidth=2)
        axes[0].set_ylabel('Pressure P')
        axes[0].set_title(f'Pressure (t={t:.3f})')
        
        axes[1].plot(x, n, 'r-', linewidth=2)
        axes[1].set_ylabel('Density n')
        axes[1].set_title(f'Density (t={t:.3f})')
        
        axes[2].plot(x, v, 'g-', linewidth=2, label='$v^x$')
        axes[2].plot(x, vy_arr, 'm-', linewidth=2, label='$v^y$')
        axes[2].plot(x, vtot, 'k--', linewidth=1, label='$|v|$')
        axes[2].set_ylabel('Velocity')
        axes[2].set_title(f'Velocities (t={t:.3f})')
        axes[2].legend()
        axes[2].set_ylim(-0.1, 1.1)
        
        axes[3].plot(x, gamma_arr, 'c-', linewidth=2)
        axes[3].set_ylabel('Lorentz factor $\\gamma$')
        axes[3].set_title(f'Lorentz Factor (t={t:.3f})')
        axes[3].set_yscale('log')
        
        axes[4].plot(x, conserved_y, 'orange', linewidth=2)
        axes[4].set_ylabel('$hW v^y$')
        axes[4].set_title(f'Conservation $hW v^y$ (t={t:.3f})')
        
        axes[5].plot(x, H, 'purple', linewidth=2)
        axes[5].set_ylabel('Enthalpy h')
        axes[5].set_title(f'Enthalpy (t={t:.3f})')
        
        # Mark wave positions
        if abs(solver.vshockl) > 1e-10:
            x_left = 0.5 + solver.vshockl * t
            for ax in axes:
                ax.axvline(x=x_left, color='red', linestyle='--', alpha=0.5)
        
        if abs(solver.vshockr) > 1e-10:
            x_right = 0.5 + solver.vshockr * t
            for ax in axes:
                ax.axvline(x=x_right, color='red', linestyle='--', alpha=0.5)
        
        # Contact discontinuity
        x_contact = 0.5 + solver.vls * t
        for ax in axes:
            ax.axvline(x=x_contact, color='blue', linestyle=':', alpha=0.5)
        
        fig.suptitle(f'Wave Structure with Tangential Velocity ($v^y = {vyl}$)\n' +
                    'Pons et al. (1999) - Kitajima Formulation',
                    fontsize=12, fontweight='bold')
        plt.tight_layout()
    
    anim = animation.FuncAnimation(fig, animate, frames=n_frames, 
                                  interval=50, blit=False)
    
    output_file = 'wave_structure_tangential.mp4'
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer='ffmpeg', fps=20, dpi=150)
    print(f"✓ Animation saved: {output_file}")
    
    plt.close()


def create_blast_wave_animation():
    """
    Create animation of relativistic blast wave with tangential velocity.
    """
    print("\nCreating blast wave animation...")
    
    gamma = 5.0/3.0
    c = 1.0
    
    # Initial conditions from Pons et al. Table 1
    Pl, nl, vl = 1000.0, 1.0, 0.0
    Pr, nr, vr = 0.01, 1.0, 0.0
    
    # Test three cases
    cases = [
        (0.0, 0.0, 'No tangential velocity'),
        (0.9, 0.9, 'Moderate tangential velocity ($v^t = 0.9$)'),
        (0.99, 0.99, 'High tangential velocity ($v^t = 0.99$)')
    ]
    
    # Time steps
    t_max = 0.4
    n_frames = 100
    times = np.linspace(0.001, t_max, n_frames)
    
    # Create figure
    fig, axes = plt.subplots(3, 3, figsize=(15, 10))
    
    def animate(frame):
        t = times[frame]
        
        for i, (vyl, vyr, title) in enumerate(cases):
            solver = KitajimaRiemannSolver(gamma, c)
            solver.set_initial_states(Pl, nl, vl, Pr, nr, vr,
                                     vyl=vyl, vzl=0.0, vyr=vyr, vzr=0.0)
            
            x, P, n, N, v, vy_arr, vz_arr, u, gamma_arr, S, e = solver.solve(
                t, x0=0.5, n_points=400)
            
            # Plot pressure
            axes[i, 0].clear()
            axes[i, 0].semilogy(x, P, 'b-', linewidth=2)
            axes[i, 0].set_ylabel('Pressure P')
            axes[i, 0].set_xlabel('x')
            axes[i, 0].grid(True, alpha=0.3)
            axes[i, 0].set_title(f'{title} - Pressure')
            axes[i, 0].set_ylim(1e-3, 1e4)
            
            # Plot density
            axes[i, 1].clear()
            axes[i, 1].semilogy(x, n, 'r-', linewidth=2)
            axes[i, 1].set_ylabel('Density n')
            axes[i, 1].set_xlabel('x')
            axes[i, 1].grid(True, alpha=0.3)
            axes[i, 1].set_title(f'{title} - Density')
            axes[i, 1].set_ylim(1e-3, 1e2)
            
            # Plot velocities
            axes[i, 2].clear()
            axes[i, 2].plot(x, v, 'g-', linewidth=2, label='$v^x$')
            axes[i, 2].plot(x, vy_arr, 'm-', linewidth=2, label='$v^y$')
            axes[i, 2].set_ylabel('Velocity')
            axes[i, 2].set_xlabel('x')
            axes[i, 2].grid(True, alpha=0.3)
            axes[i, 2].set_title(f'{title} - Velocities')
            axes[i, 2].legend()
            axes[i, 2].set_ylim(-0.1, 1.1)
        
        fig.suptitle(f'Relativistic Blast Wave with Tangential Velocities (t={t:.3f})\n' +
                    'Pons et al. (1999) Table 1 - Kitajima Formulation',
                    fontsize=12, fontweight='bold')
        plt.tight_layout()
    
    anim = animation.FuncAnimation(fig, animate, frames=n_frames,
                                  interval=50, blit=False)
    
    output_file = 'blast_wave_tangential.mp4'
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer='ffmpeg', fps=20, dpi=150)
    print(f"✓ Animation saved: {output_file}")
    
    plt.close()


if __name__ == "__main__":
    print("=" * 80)
    print("TANGENTIAL VELOCITY ANIMATIONS")
    print("Based on Pons, Martí & Müller (1999, J. Fluid Mech.)")
    print("=" * 80)
    
    try:
        # Create animations
        create_comparison_animation()
        create_wave_structure_animation()
        create_blast_wave_animation()
        
        print("\n" + "=" * 80)
        print("✓ All animations created successfully!")
        print("=" * 80)
        print("\nGenerated files:")
        print("  1. tangential_velocity_comparison.mp4")
        print("  2. wave_structure_tangential.mp4")
        print("  3. blast_wave_tangential.mp4")
        print("\nThese animations demonstrate:")
        print("  - Tangential velocity effects on wave structure")
        print("  - Conservation of h*W*v^y across waves")
        print("  - Comparison with Pons et al. (1999) results")
        print("=" * 80)
        
    except Exception as e:
        print(f"\n✗ Error creating animations: {e}")
        import traceback
        traceback.print_exc()
