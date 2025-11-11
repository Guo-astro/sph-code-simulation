#!/usr/bin/env python3
"""
Create comparison animation for Sedov blast wave showing multiple SPH methods
against analytical solution.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import pandas as pd
import argparse
from pathlib import Path
import sys

# Import the analytical solution
from sedov_analytical import SedovSolution, load_snapshot

def create_comparison_animation(results_dir, output_file, methods=None, fps=5):
    """
    Create animation comparing multiple SPH methods with analytical solution.
    
    Parameters:
    -----------
    results_dir : str or Path
        Directory containing method subdirectories (e.g., gsph_wendland/)
    output_file : str or Path
        Output path for the animation GIF
    methods : list of str, optional
        List of method names to include. Default: all available methods
    fps : int
        Frames per second for the animation
    """
    results_dir = Path(results_dir)
    
    # Default methods to compare
    if methods is None:
        methods = ['gsph_wendland', 'ssph_wendland', 'disph_wendland', 'gdisph_wendland']
    
    # Method display names
    method_names = {
        'gsph_wendland': 'GSPH',
        'ssph_wendland': 'SSPH',
        'disph_wendland': 'DISPH',
        'gdisph_wendland': 'GDISPH',
        'gsph_cubic': 'GSPH (Cubic)',
        'ssph_cubic': 'SSPH (Cubic)',
        'disph_cubic': 'DISPH (Cubic)',
        'gdisph_cubic': 'GDISPH (Cubic)',
    }
    
    # Filter to only existing methods
    existing_methods = []
    for method in methods:
        method_dir = results_dir / method
        if method_dir.exists():
            existing_methods.append(method)
        else:
            print(f"Warning: {method} directory not found, skipping")
    
    if not existing_methods:
        print("Error: No method directories found!")
        return
    
    print(f"Creating animation for {len(existing_methods)} methods: {', '.join(existing_methods)}")
    
    # Get list of snapshots from first method
    first_method_dir = results_dir / existing_methods[0]
    snapshots = sorted(first_method_dir.glob('snapshot_*.csv'))
    
    if not snapshots:
        print(f"Error: No snapshots found in {first_method_dir}")
        return
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Setup figure
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    axes = [
        fig.add_subplot(gs[0, 0]),  # Density
        fig.add_subplot(gs[0, 1]),  # Velocity
        fig.add_subplot(gs[1, 0]),  # Pressure
        fig.add_subplot(gs[1, 1]),  # Internal Energy
    ]
    
    # Colors for different methods
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    
    def update(frame):
        """Update function for animation."""
        for ax in axes:
            ax.clear()
        
        snapshot_idx = frame
        snapshot_file = snapshots[snapshot_idx]
        
        # Load analytical solution (same for all methods)
        data_ref, metadata_ref = load_snapshot(snapshot_file)
        time = float(metadata_ref.get('Time', '0.0'))
        gamma = float(metadata_ref.get('Gamma', '1.4'))
        
        # Compute analytical solution
        sedov = SedovSolution(gamma=gamma, E0=1.0, rho0=1.0, nu=2)
        r_analytical, rho_analytical, v_analytical, p_analytical, e_analytical = \
            sedov.solution_at_time(time, n_points=200)
        
        # Plot analytical solution (same on all panels)
        axes[0].plot(r_analytical, rho_analytical, 'k-', linewidth=2, label='Sedov', alpha=0.7, zorder=10)
        axes[1].plot(r_analytical, v_analytical, 'k-', linewidth=2, label='Sedov', alpha=0.7, zorder=10)
        axes[2].plot(r_analytical, p_analytical, 'k-', linewidth=2, label='Sedov', alpha=0.7, zorder=10)
        axes[3].plot(r_analytical, e_analytical, 'k-', linewidth=2, label='Sedov', alpha=0.7, zorder=10)
        
        # Plot each SPH method
        for idx, method in enumerate(existing_methods):
            method_snapshot = results_dir / method / snapshot_file.name
            
            if not method_snapshot.exists():
                continue
            
            # Load SPH data
            data, metadata = load_snapshot(method_snapshot)
            
            # Compute radial distance
            x_col = 'x' if 'x' in data.columns else 'pos_x'
            y_col = 'y' if 'y' in data.columns else 'pos_y'
            r_sph = np.sqrt(data[x_col]**2 + data[y_col]**2)
            
            # Get physical quantities
            rho_sph = data['dens']
            v_x = data.get('vel_x', data.get('vx', 0))
            v_y = data.get('vel_y', data.get('vy', 0))
            v_sph = np.sqrt(v_x**2 + v_y**2)
            p_sph = data['pres']
            e_sph = data['ene']
            
            color = colors[idx % len(colors)]
            method_label = method_names.get(method, method)
            
            # Plot SPH data
            axes[0].scatter(r_sph, rho_sph, c=color, s=10, alpha=0.5, label=method_label)
            axes[1].scatter(r_sph, v_sph, c=color, s=10, alpha=0.5, label=method_label)
            axes[2].scatter(r_sph, p_sph, c=color, s=10, alpha=0.5, label=method_label)
            axes[3].scatter(r_sph, e_sph, c=color, s=10, alpha=0.5, label=method_label)
        
        # Set labels and formatting
        axes[0].set_ylabel('Density', fontsize=12)
        axes[0].set_xlim(0, 0.5)
        axes[0].set_ylim(0, 7)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend(loc='upper right', fontsize=9)
        
        axes[1].set_ylabel('Velocity', fontsize=12)
        axes[1].set_xlim(0, 0.5)
        axes[1].set_ylim(0, 1.5)
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(loc='upper right', fontsize=9)
        
        axes[2].set_xlabel('Radius', fontsize=12)
        axes[2].set_ylabel('Pressure', fontsize=12)
        axes[2].set_xlim(0, 0.5)
        axes[2].set_yscale('log')
        axes[2].set_ylim(1e-6, 10)
        axes[2].grid(True, alpha=0.3, which='both')
        axes[2].legend(loc='upper right', fontsize=9)
        
        axes[3].set_xlabel('Radius', fontsize=12)
        axes[3].set_ylabel('Internal Energy', fontsize=12)
        axes[3].set_xlim(0, 0.5)
        axes[3].set_yscale('log')
        axes[3].set_ylim(1e-6, 10)
        axes[3].grid(True, alpha=0.3, which='both')
        axes[3].legend(loc='upper right', fontsize=9)
        
        # Title
        fig.suptitle(f'Sedov Blast Wave: SPH Methods vs Analytical Solution (t = {time:.4f})',
                     fontsize=14, fontweight='bold', y=0.98)
        
        return axes
    
    # Create animation
    print("Creating animation frames...")
    anim = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                   interval=1000/fps, blit=False, repeat=True)
    
    # Save animation
    print(f"Saving animation to {output_file}...")
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Use pillow writer for GIF
    writer = animation.PillowWriter(fps=fps, metadata=dict(artist='SPH Simulation'))
    anim.save(output_path, writer=writer, dpi=100)
    
    plt.close()
    print(f"✓ Animation saved: {output_file}")
    print(f"  Frames: {len(snapshots)}")
    print(f"  Methods: {len(existing_methods)}")
    print(f"  FPS: {fps}")

def main():
    parser = argparse.ArgumentParser(description='Create Sedov comparison animation')
    parser.add_argument('results_dir', type=str,
                       help='Results directory containing method subdirectories')
    parser.add_argument('output', type=str,
                       help='Output GIF file path')
    parser.add_argument('--methods', type=str, nargs='+',
                       help='Methods to include (default: all Wendland methods)')
    parser.add_argument('--fps', type=int, default=5,
                       help='Frames per second (default: 5)')
    
    args = parser.parse_args()
    
    create_comparison_animation(args.results_dir, args.output, 
                               methods=args.methods, fps=args.fps)

if __name__ == '__main__':
    main()
