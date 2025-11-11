#!/usr/bin/env python3
"""
Generate analytical comparison animations for all SPH variants and kernels.

This script creates comparison plots and animations showing SPH simulation
results against the corrected Sedov analytical solution for:
- SPH types: DISPH, GDISPH, GSPH, SSPH
- Kernels: Wendland C4, Cubic Spline
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path
import sys
import os

# Add parent directory to path to import sedov_analytical
sys.path.insert(0, str(Path(__file__).parent))
from sedov_analytical import SedovSolution, load_snapshot

def create_comparison_plot(snapshot_file, output_file, title_suffix=""):
    """Create a single comparison plot with analytical solution."""
    
    # Load SPH data
    data, metadata = load_snapshot(snapshot_file)
    time = float(metadata.get('Time', 0.0))
    gamma = float(metadata.get('Gamma', '1.4'))
    
    # Compute radial distance
    x_col = 'pos_x'
    y_col = 'pos_y'
    r_sph = np.sqrt(data[x_col]**2 + data[y_col]**2)
    
    # Get physical quantities
    rho_sph = data['dens']
    v_x = data['vel_x']
    v_y = data['vel_y']
    v_sph = np.sqrt(v_x**2 + v_y**2)
    p_sph = data['pres']
    e_sph = data['ene']
    
    # Compute analytical solution
    sedov = SedovSolution(gamma=gamma, E0=1.0, rho0=1.0, nu=2)
    r_anal, rho_anal, v_anal, p_anal, e_anal = sedov.solution_at_time(time, n_points=500)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(f'Sedov Blast Wave - {title_suffix} - t = {time:.4f}', 
                 fontsize=16, fontweight='bold')
    
    # Density
    ax = axes[0, 0]
    ax.scatter(r_sph, rho_sph, s=3, alpha=0.5, label='SPH', color='blue')
    ax.plot(r_anal, rho_anal, 'r-', linewidth=2, label='Analytical')
    ax.set_xlabel('Radius', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Density Profile', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.5)
    
    # Velocity
    ax = axes[0, 1]
    ax.scatter(r_sph, v_sph, s=3, alpha=0.5, label='SPH', color='blue')
    ax.plot(r_anal, v_anal, 'r-', linewidth=2, label='Analytical')
    ax.set_xlabel('Radius', fontsize=12)
    ax.set_ylabel('Velocity', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.5)
    
    # Pressure
    ax = axes[1, 0]
    ax.scatter(r_sph, p_sph, s=3, alpha=0.5, label='SPH', color='blue')
    ax.plot(r_anal, p_anal, 'r-', linewidth=2, label='Analytical')
    ax.set_xlabel('Radius', fontsize=12)
    ax.set_ylabel('Pressure', fontsize=12)
    ax.set_title('Pressure Profile', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xlim(0, 0.5)
    
    # Internal Energy
    ax = axes[1, 1]
    ax.scatter(r_sph, e_sph, s=3, alpha=0.5, label='SPH', color='blue')
    ax.plot(r_anal, e_anal, 'r-', linewidth=2, label='Analytical')
    ax.set_xlabel('Radius', fontsize=12)
    ax.set_ylabel('Specific Internal Energy', fontsize=12)
    ax.set_title('Internal Energy Profile', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xlim(0, 0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    return time

def create_animation(result_dir, output_file, title_suffix=""):
    """Create an animation comparing SPH and analytical solutions over time."""
    
    result_path = Path(result_dir)
    snapshots = sorted(result_path.glob('snapshot_*.csv'))
    
    if len(snapshots) == 0:
        print(f"No snapshots found in {result_dir}")
        return
    
    print(f"Creating animation with {len(snapshots)} snapshots...")
    
    # Load first snapshot to get metadata
    data0, metadata0 = load_snapshot(snapshots[0])
    gamma = float(metadata0.get('Gamma', '1.4'))
    sedov = SedovSolution(gamma=gamma, E0=1.0, rho0=1.0, nu=2)
    
    # Set up the figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(f'Sedov Blast Wave - {title_suffix}', fontsize=16, fontweight='bold')
    
    # Initialize plots
    scatter_plots = []
    line_plots = []
    
    for idx, ax in enumerate(axes.flat):
        scatter = ax.scatter([], [], s=3, alpha=0.5, label='SPH', color='blue')
        line, = ax.plot([], [], 'r-', linewidth=2, label='Analytical')
        scatter_plots.append(scatter)
        line_plots.append(line)
        
        ax.set_xlim(0, 0.5)
        ax.set_xlabel('Radius', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
    
    # Set labels
    axes[0, 0].set_ylabel('Density', fontsize=12)
    axes[0, 0].set_title('Density Profile', fontsize=13, fontweight='bold')
    axes[0, 0].set_ylim(0, 4)
    
    axes[0, 1].set_ylabel('Velocity', fontsize=12)
    axes[0, 1].set_title('Velocity Profile', fontsize=13, fontweight='bold')
    axes[0, 1].set_ylim(0, 2)
    
    axes[1, 0].set_ylabel('Pressure', fontsize=12)
    axes[1, 0].set_title('Pressure Profile', fontsize=13, fontweight='bold')
    axes[1, 0].set_yscale('log')
    axes[1, 0].set_ylim(1e-6, 10)
    
    axes[1, 1].set_ylabel('Specific Internal Energy', fontsize=12)
    axes[1, 1].set_title('Internal Energy Profile', fontsize=13, fontweight='bold')
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_ylim(1e-6, 100)
    
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=14, fontweight='bold')
    
    def init():
        """Initialize animation."""
        for scatter in scatter_plots:
            scatter.set_offsets(np.empty((0, 2)))
        for line in line_plots:
            line.set_data([], [])
        time_text.set_text('')
        return scatter_plots + line_plots + [time_text]
    
    def update(frame):
        """Update animation frame."""
        snapshot_file = snapshots[frame]
        
        # Load SPH data
        data, metadata = load_snapshot(snapshot_file)
        time = float(metadata.get('Time', 0.0))
        
        # Compute radial distance
        r_sph = np.sqrt(data['pos_x']**2 + data['pos_y']**2)
        rho_sph = data['dens'].values
        v_sph = np.sqrt(data['vel_x']**2 + data['vel_y']**2).values
        p_sph = data['pres'].values
        e_sph = data['ene'].values
        
        # Analytical solution
        r_anal, rho_anal, v_anal, p_anal, e_anal = sedov.solution_at_time(time, n_points=500)
        
        # Update scatter plots
        scatter_plots[0].set_offsets(np.column_stack([r_sph, rho_sph]))
        scatter_plots[1].set_offsets(np.column_stack([r_sph, v_sph]))
        scatter_plots[2].set_offsets(np.column_stack([r_sph, p_sph]))
        scatter_plots[3].set_offsets(np.column_stack([r_sph, e_sph]))
        
        # Update line plots
        line_plots[0].set_data(r_anal, rho_anal)
        line_plots[1].set_data(r_anal, v_anal)
        line_plots[2].set_data(r_anal, p_anal)
        line_plots[3].set_data(r_anal, e_anal)
        
        # Update time
        time_text.set_text(f't = {time:.4f}')
        
        return scatter_plots + line_plots + [time_text]
    
    # Create animation
    anim = FuncAnimation(fig, update, init_func=init, frames=len(snapshots),
                        interval=200, blit=True, repeat=True)
    
    # Save animation
    anim.save(output_file, writer='pillow', fps=5, dpi=100)
    plt.close()
    
    print(f'✓ Saved animation: {output_file}')

def main():
    """Generate all comparison plots and animations."""
    
    base_dir = Path(__file__).parent.parent
    results_dir = base_dir / 'results'
    
    # Define all configurations
    sph_types = ['disph', 'gdisph', 'gsph', 'ssph']
    kernels = ['wendland', 'cubic']
    
    print("="*80)
    print("GENERATING SEDOV ANALYTICAL COMPARISONS")
    print("="*80)
    
    # Track which result directories exist
    existing_results = []
    for sph_type in sph_types:
        for kernel in kernels:
            result_name = f"{sph_type}_{kernel}"
            result_path = results_dir / result_name
            
            if result_path.exists():
                existing_results.append((sph_type, kernel, result_name, result_path))
    
    print(f"\nFound {len(existing_results)} result directories:")
    for sph_type, kernel, result_name, _ in existing_results:
        print(f"  - {result_name}")
    
    # Generate animations for each existing result
    animations_dir = results_dir / 'animations'
    animations_dir.mkdir(exist_ok=True)
    
    print(f"\n{'='*80}")
    print("GENERATING ANIMATIONS")
    print(f"{'='*80}\n")
    
    for sph_type, kernel, result_name, result_path in existing_results:
        # Create title
        sph_name = sph_type.upper()
        kernel_name = "Wendland C4" if kernel == "wendland" else "Cubic Spline"
        title = f"{sph_name} - {kernel_name}"
        
        # Create animation
        output_file = animations_dir / f"{result_name}_comparison.gif"
        print(f"Processing {result_name}...")
        
        try:
            create_animation(result_path, output_file, title_suffix=title)
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    # Also create final snapshot comparisons
    print(f"\n{'='*80}")
    print("GENERATING FINAL SNAPSHOT COMPARISONS")
    print(f"{'='*80}\n")
    
    comparison_dir = results_dir / 'comparison'
    comparison_dir.mkdir(exist_ok=True)
    
    for sph_type, kernel, result_name, result_path in existing_results:
        # Find final snapshot
        snapshots = sorted(result_path.glob('snapshot_*.csv'))
        if len(snapshots) == 0:
            continue
        
        final_snapshot = snapshots[-1]
        
        # Create title
        sph_name = sph_type.upper()
        kernel_name = "Wendland C4" if kernel == "wendland" else "Cubic Spline"
        title = f"{sph_name} - {kernel_name}"
        
        # Create comparison plot
        output_file = comparison_dir / f"{result_name}_final.png"
        print(f"Creating final comparison for {result_name}...")
        
        try:
            time = create_comparison_plot(final_snapshot, output_file, title_suffix=title)
            print(f"  ✓ Saved: {output_file} (t={time:.4f})")
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    print(f"\n{'='*80}")
    print("COMPLETE!")
    print(f"{'='*80}")
    print(f"\nAnimations saved to: {animations_dir}")
    print(f"Final comparisons saved to: {comparison_dir}")

if __name__ == '__main__':
    main()
