#!/usr/bin/env python3
"""
Plot SR-GSPH Tangent Velocity Test Snapshots

Visualizes the evolution of the relativistic shock tube with tangent velocity
from Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.1.5.

This test demonstrates that 1D SR-GSPH can correctly evolve systems with
transverse (tangent) velocity components, even though the simulation is
spatially 1-dimensional.
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
RESULTS_DIR = Path(__file__).parent.parent / "results" / "tangent_vt099_vt099"
OUTPUT_DIR = RESULTS_DIR  # Save plots in results directory


def read_csv_snapshot(filepath):
    """Read a CSV snapshot file, skipping comment lines."""
    data = {'time': 0.0, 'step': 0}
    
    # First pass: extract metadata from comments
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('# Time (code):'):
                data['time'] = float(line.split(':')[1].strip())
            elif line.startswith('# Step:'):
                data['step'] = int(line.split(':')[1].strip())
            elif not line.startswith('#'):
                break
    
    # Second pass: read the actual data
    # Skip all comment lines
    df_data = []
    header = None
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if header is None:
                header = line.strip().split(',')
            else:
                values = line.strip().split(',')
                df_data.append([float(v) if v else 0.0 for v in values])
    
    if not df_data:
        return None
    
    arr = np.array(df_data)
    result = {col: arr[:, i] for i, col in enumerate(header)}
    result['time'] = data['time']
    result['step'] = data['step']
    return result


def plot_all_snapshots():
    """Create a multi-panel plot showing all snapshots."""
    
    # Find all snapshot files
    snapshot_files = sorted(glob.glob(str(RESULTS_DIR / "snapshot_*.csv")))
    
    if not snapshot_files:
        print(f"No snapshots found in {RESULTS_DIR}")
        return
    
    print(f"Found {len(snapshot_files)} snapshots")
    
    # Read all snapshots
    snapshots = []
    for f in snapshot_files:
        data = read_csv_snapshot(f)
        if data is not None:
            snapshots.append(data)
    
    n_snapshots = len(snapshots)
    if n_snapshots == 0:
        print("No valid snapshots found")
        return
    
    # Create figure with subplots
    fig, axes = plt.subplots(4, 1, figsize=(14, 16), sharex=True)
    fig.suptitle(r'SR-GSPH Tangent Velocity Test ($v^t_L = v^t_R = 0.99$, $\gamma \approx 7.09$)', 
                 fontsize=14, fontweight='bold')
    
    # Color map for time evolution
    colors = plt.cm.viridis(np.linspace(0, 1, n_snapshots))
    
    # Filter out ghost particles (is_ghost == 0 means real particle)
    for i, snap in enumerate(snapshots):
        mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
        
        x = snap['pos_x'][mask]
        vel_x = snap['vel_x'][mask]
        vel_t = snap['vel_t'][mask]
        dens = snap['dens'][mask]
        pres = snap['pres'][mask]
        t = snap['time']
        
        # Sort by position
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        vel_x = vel_x[sort_idx]
        vel_t = vel_t[sort_idx]
        dens = dens[sort_idx]
        pres = pres[sort_idx]
        
        alpha = 0.3 + 0.7 * (i / max(1, n_snapshots - 1))  # Fade earlier snapshots
        
        # Plot density
        axes[0].plot(x, dens, '-', color=colors[i], alpha=alpha, linewidth=0.8,
                     label=f't={t:.4f}' if i % max(1, n_snapshots // 6) == 0 else '')
        
        # Plot pressure
        axes[1].plot(x, pres, '-', color=colors[i], alpha=alpha, linewidth=0.8)
        
        # Plot normal velocity
        axes[2].plot(x, vel_x, '-', color=colors[i], alpha=alpha, linewidth=0.8)
        
        # Plot tangent velocity
        axes[3].plot(x, vel_t, '-', color=colors[i], alpha=alpha, linewidth=0.8)
    
    # Formatting
    axes[0].set_ylabel(r'Density $N = \gamma n$', fontsize=12)
    axes[0].set_yscale('log')
    axes[0].legend(loc='upper right', fontsize=8, ncol=2)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_ylabel(r'Pressure $P$', fontsize=12)
    axes[1].set_yscale('log')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_ylabel(r'Normal Velocity $v^x/c$', fontsize=12)
    axes[2].grid(True, alpha=0.3)
    axes[2].axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    
    axes[3].set_ylabel(r'Tangent Velocity $v^t/c$', fontsize=12)
    axes[3].set_xlabel(r'Position $x$', fontsize=12)
    axes[3].grid(True, alpha=0.3)
    axes[3].axhline(y=0.99, color='r', linestyle='--', linewidth=0.5, label=r'$v^t_{initial}=0.99$')
    axes[3].legend(loc='lower right', fontsize=8)
    
    # Add colorbar for time
    sm = plt.cm.ScalarMappable(cmap='viridis', 
                               norm=plt.Normalize(vmin=snapshots[0]['time'], 
                                                   vmax=snapshots[-1]['time']))
    cbar = fig.colorbar(sm, ax=axes, label='Time (code units)', shrink=0.8)
    
    plt.tight_layout()
    
    # Save
    output_path = OUTPUT_DIR / "tangent_velocity_evolution.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot to {output_path}")
    
    # Also create individual snapshot plots
    plot_individual_snapshots(snapshots)
    
    plt.show()


def plot_individual_snapshots(snapshots):
    """Create a grid of individual snapshots."""
    
    n_snapshots = len(snapshots)
    n_cols = 3
    n_rows = (n_snapshots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    fig.suptitle(r'SR-GSPH Tangent Velocity: Individual Snapshots', fontsize=14, fontweight='bold')
    
    axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes
    
    for i, snap in enumerate(snapshots):
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # Filter ghost particles
        mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
        x = snap['pos_x'][mask]
        dens = snap['dens'][mask]
        pres = snap['pres'][mask]
        vel_x = snap['vel_x'][mask]
        vel_t = snap['vel_t'][mask]
        
        # Sort by position
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        dens = dens[sort_idx]
        pres = pres[sort_idx]
        vel_x = vel_x[sort_idx]
        vel_t = vel_t[sort_idx]
        
        # Plot density on left y-axis
        ax.plot(x, dens, 'b-', linewidth=1, label=r'$N$')
        ax.set_ylabel(r'Density $N$', color='b')
        ax.tick_params(axis='y', labelcolor='b')
        
        # Plot tangent velocity on right y-axis
        ax2 = ax.twinx()
        ax2.plot(x, vel_t, 'r-', linewidth=1, label=r'$v^t$')
        ax2.set_ylabel(r'$v^t/c$', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        ax2.set_ylim([0.98, 1.01])
        
        ax.set_title(f't = {snap["time"]:.4f}', fontsize=10)
        ax.set_xlabel('x')
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for i in range(n_snapshots, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / "tangent_velocity_snapshots_grid.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved grid plot to {output_path}")


def plot_velocity_profiles():
    """Create detailed velocity profile comparison."""
    
    snapshot_files = sorted(glob.glob(str(RESULTS_DIR / "snapshot_*.csv")))
    if not snapshot_files:
        return
    
    # Select a few key snapshots
    indices = [0, len(snapshot_files)//4, len(snapshot_files)//2, 
               3*len(snapshot_files)//4, len(snapshot_files)-1]
    indices = [i for i in indices if i < len(snapshot_files)]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(r'SR-GSPH: Velocity Components Evolution', fontsize=14, fontweight='bold')
    
    colors = plt.cm.plasma(np.linspace(0, 1, len(indices)))
    
    for j, idx in enumerate(indices):
        snap = read_csv_snapshot(snapshot_files[idx])
        if snap is None:
            continue
        
        mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
        x = snap['pos_x'][mask]
        vel_x = snap['vel_x'][mask]
        vel_t = snap['vel_t'][mask]
        
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        vel_x = vel_x[sort_idx]
        vel_t = vel_t[sort_idx]
        
        # Total velocity magnitude
        v_total = np.sqrt(vel_x**2 + vel_t**2)
        
        # Lorentz factor
        gamma = 1.0 / np.sqrt(1.0 - v_total**2)
        
        label = f't={snap["time"]:.4f}'
        
        axes[0, 0].plot(x, vel_x, '-', color=colors[j], linewidth=1.5, label=label)
        axes[0, 1].plot(x, vel_t, '-', color=colors[j], linewidth=1.5, label=label)
        axes[1, 0].plot(x, v_total, '-', color=colors[j], linewidth=1.5, label=label)
        axes[1, 1].plot(x, gamma, '-', color=colors[j], linewidth=1.5, label=label)
    
    axes[0, 0].set_ylabel(r'Normal Velocity $v^x/c$', fontsize=12)
    axes[0, 0].set_xlabel('x')
    axes[0, 0].axhline(0, color='k', linestyle='--', linewidth=0.5)
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title('Normal (x) Velocity Component')
    
    axes[0, 1].set_ylabel(r'Tangent Velocity $v^t/c$', fontsize=12)
    axes[0, 1].set_xlabel('x')
    axes[0, 1].axhline(0.99, color='k', linestyle='--', linewidth=0.5, label='Initial $v^t=0.99$')
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_title('Tangent (perpendicular) Velocity Component')
    
    axes[1, 0].set_ylabel(r'Total Velocity $|v|/c$', fontsize=12)
    axes[1, 0].set_xlabel('x')
    axes[1, 0].axhline(1.0, color='r', linestyle='--', linewidth=0.5, label='Speed of light')
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_title(r'Total Velocity Magnitude $|v| = \sqrt{v_x^2 + v_t^2}$')
    
    axes[1, 1].set_ylabel(r'Lorentz Factor $\gamma$', fontsize=12)
    axes[1, 1].set_xlabel('x')
    axes[1, 1].axhline(7.0888, color='k', linestyle='--', linewidth=0.5, label=r'Initial $\gamma=7.09$')
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_title(r'Lorentz Factor $\gamma = 1/\sqrt{1-|v|^2/c^2}$')
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / "tangent_velocity_profiles.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved velocity profiles to {output_path}")


def create_animation():
    """Create an animation of the shock tube evolution."""
    from matplotlib.animation import FuncAnimation, PillowWriter
    
    snapshot_files = sorted(glob.glob(str(RESULTS_DIR / "snapshot_*.csv")))
    if not snapshot_files:
        print("No snapshots found for animation")
        return
    
    # Read all snapshots
    snapshots = []
    for f in snapshot_files:
        data = read_csv_snapshot(f)
        if data is not None:
            snapshots.append(data)
    
    n_snapshots = len(snapshots)
    print(f"Creating animation with {n_snapshots} frames...")
    
    # Set up the figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(r'SR-GSPH Tangent Velocity Test ($v^t = 0.99$, $\gamma \approx 7.09$)', 
                 fontsize=14, fontweight='bold')
    
    # Initialize empty lines
    lines = []
    for ax in axes.flatten():
        line, = ax.plot([], [], 'b-', linewidth=1.5)
        lines.append(line)
    
    # Configure axes
    axes[0, 0].set_ylabel(r'Density $N$', fontsize=12)
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_xlim(-0.5, 0.5)
    axes[0, 0].set_ylim(0, 10)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title('Density')
    
    axes[0, 1].set_ylabel(r'Pressure $P$', fontsize=12)
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_xlim(-0.5, 0.5)
    axes[0, 1].set_ylim(1e-3, 1e4)
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_title('Pressure')
    
    axes[1, 0].set_ylabel(r'Normal Velocity $v^x/c$', fontsize=12)
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_xlim(-0.5, 0.5)
    axes[1, 0].set_ylim(-0.5, 0.5)
    axes[1, 0].axhline(0, color='k', linestyle='--', linewidth=0.5)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_title('Normal Velocity')
    
    axes[1, 1].set_ylabel(r'Tangent Velocity $v^t/c$', fontsize=12)
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_xlim(-0.5, 0.5)
    axes[1, 1].set_ylim(0.95, 1.02)
    axes[1, 1].axhline(0.99, color='r', linestyle='--', linewidth=0.5, label=r'$v^t_{init}=0.99$')
    axes[1, 1].axhline(1.0, color='k', linestyle=':', linewidth=0.5, label='c')
    axes[1, 1].legend(loc='lower right', fontsize=8)
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_title('Tangent Velocity')
    
    # Time text
    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12, fontweight='bold')
    
    def init():
        for line in lines:
            line.set_data([], [])
        time_text.set_text('')
        return lines + [time_text]
    
    def animate(frame):
        snap = snapshots[frame]
        
        # Filter ghost particles
        mask = snap.get('is_ghost', np.zeros(len(snap['pos_x']))) == 0
        x = snap['pos_x'][mask]
        dens = snap['dens'][mask]
        pres = snap['pres'][mask]
        vel_x = snap['vel_x'][mask]
        vel_t = snap['vel_t'][mask]
        
        # Sort by position
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        dens = dens[sort_idx]
        pres = pres[sort_idx]
        vel_x = vel_x[sort_idx]
        vel_t = vel_t[sort_idx]
        
        lines[0].set_data(x, dens)
        lines[1].set_data(x, pres)
        lines[2].set_data(x, vel_x)
        lines[3].set_data(x, vel_t)
        
        time_text.set_text(f't = {snap["time"]:.5f} (frame {frame+1}/{n_snapshots})')
        
        return lines + [time_text]
    
    anim = FuncAnimation(fig, animate, init_func=init, frames=n_snapshots,
                         interval=300, blit=True)
    
    # Save as GIF
    output_path = OUTPUT_DIR / "tangent_velocity_animation.gif"
    print(f"Saving animation to {output_path}...")
    anim.save(output_path, writer=PillowWriter(fps=3))
    print(f"Animation saved!")
    
    plt.close(fig)


if __name__ == "__main__":
    print("Plotting SR-GSPH Tangent Velocity Test Results")
    print("=" * 50)
    
    if not RESULTS_DIR.exists():
        print(f"Error: Results directory not found: {RESULTS_DIR}")
        sys.exit(1)
    
    plot_all_snapshots()
    plot_velocity_profiles()
    create_animation()
    
    print("\nDone!")
