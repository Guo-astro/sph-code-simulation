#!/usr/bin/env python3
"""
Generate analytical comparison animations for all 2D shock tube SPH variants and kernels.

This script creates comparison plots and animations showing SPH simulation
results against the exact Riemann solution for:
- SPH types: DISPH, GDISPH, GSPH, SSPH, GDISPH+Balsara
- Kernels: Wendland C4, Cubic Spline
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
import sys
import argparse

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))
from shock_tube_2d_analytical import load_snapshot, exact_riemann_solver

def create_animation_for_method(results_dir, method_name, output_dir, fps=5):
    """
    Create animation for a single SPH method.
    
    Parameters:
    -----------
    results_dir : Path
        Base results directory
    method_name : str
        Name of the method (e.g., 'gsph_wendland')
    output_dir : Path
        Output directory for animations
    fps : int
        Frames per second
    """
    method_dir = results_dir / method_name
    if not method_dir.exists():
        print(f"Skipping {method_name} (directory not found)")
        return False
    
    # Find snapshots
    snapshots = sorted(method_dir.glob('snapshot_*.csv'))
    if not snapshots:
        print(f"Skipping {method_name} (no snapshots found)")
        return False
    
    print(f"Creating animation for {method_name} ({len(snapshots)} frames)...")
    
    # Setup figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    def update(frame_idx):
        """Update function for animation."""
        snapshot_file = snapshots[frame_idx]
        
        # Load data
        data = load_snapshot(str(snapshot_file))
        
        # Extract time
        try:
            time = float(data['metadata'].get('Time (code)', '0.0'))
        except:
            time = 0.0
        
        # Extract gamma
        try:
            gamma = float(data['metadata'].get('Gamma', '1.4'))
        except:
            gamma = 1.4
        
        # Get particle data
        x = data['pos_x']
        y = data['pos_y']
        rho = data['dens']
        vx = data['vel_x']
        pres = data['pres']
        ene = data['ene']
        
        # Bin data in x-direction (average over y)
        x_min, x_max = 0.0, 1.0
        n_bins = 100
        x_bins = np.linspace(x_min, x_max, n_bins + 1)
        x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
        
        rho_avg = np.zeros(n_bins)
        vx_avg = np.zeros(n_bins)
        pres_avg = np.zeros(n_bins)
        ene_avg = np.zeros(n_bins)
        counts = np.zeros(n_bins)
        
        for i in range(len(x)):
            bin_idx = int((x[i] - x_min) / (x_max - x_min) * n_bins)
            if 0 <= bin_idx < n_bins:
                rho_avg[bin_idx] += rho[i]
                vx_avg[bin_idx] += vx[i]
                pres_avg[bin_idx] += pres[i]
                ene_avg[bin_idx] += ene[i]
                counts[bin_idx] += 1
        
        # Average
        mask = counts > 0
        rho_avg[mask] /= counts[mask]
        vx_avg[mask] /= counts[mask]
        pres_avg[mask] /= counts[mask]
        ene_avg[mask] /= counts[mask]
        
        # Analytical solution
        x_exact = np.linspace(x_min, x_max, 1000)
        rho_exact, u_exact, p_exact, e_exact = exact_riemann_solver(
            x_exact, time, gamma=gamma, x0=0.5,
            rho_L=1.0, rho_R=0.125, p_L=1.0, p_R=0.1
        )
        
        # Clear axes
        for ax in axes.flat:
            ax.clear()
        
        # Density
        ax = axes[0, 0]
        ax.plot(x_exact, rho_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
        ax.scatter(x_centers[mask], rho_avg[mask], c='red', s=20, alpha=0.6, label='SPH', zorder=5)
        ax.set_xlabel('x')
        ax.set_ylabel('Density')
        ax.set_title('Density Profile')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        # Velocity
        ax = axes[0, 1]
        ax.plot(x_exact, u_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
        ax.scatter(x_centers[mask], vx_avg[mask], c='blue', s=20, alpha=0.6, label='SPH', zorder=5)
        ax.set_xlabel('x')
        ax.set_ylabel('Velocity')
        ax.set_title('Velocity Profile')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        # Pressure
        ax = axes[1, 0]
        ax.plot(x_exact, p_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
        ax.scatter(x_centers[mask], pres_avg[mask], c='green', s=20, alpha=0.6, label='SPH', zorder=5)
        ax.set_xlabel('x')
        ax.set_ylabel('Pressure')
        ax.set_title('Pressure Profile')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        # Internal Energy
        ax = axes[1, 1]
        ax.plot(x_exact, e_exact, 'k-', linewidth=2, label='Analytical', alpha=0.7)
        ax.scatter(x_centers[mask], ene_avg[mask], c='purple', s=20, alpha=0.6, label='SPH', zorder=5)
        ax.set_xlabel('x')
        ax.set_ylabel('Specific Internal Energy')
        ax.set_title('Internal Energy Profile')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        # Overall title
        fig.suptitle(f'2D Shock Tube - {method_name} - t = {time:.4f}', 
                     fontsize=14, fontweight='bold')
        
        plt.tight_layout()
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(snapshots), interval=1000//fps)
    
    # Save animation
    output_file = output_dir / f'{method_name}_animation.gif'
    writer = PillowWriter(fps=fps)
    anim.save(str(output_file), writer=writer, dpi=100)
    plt.close()
    
    print(f"  ✓ Saved: {output_file}")
    return True

def main():
    parser = argparse.ArgumentParser(description='Generate animations for all 2D shock tube methods')
    parser.add_argument('--results-dir', type=Path, default=Path('sample/shock_tube_2d/results'),
                       help='Base results directory')
    parser.add_argument('--fps', type=int, default=5, help='Frames per second')
    
    args = parser.parse_args()
    
    results_dir = args.results_dir
    animations_dir = results_dir / 'animations'
    animations_dir.mkdir(parents=True, exist_ok=True)
    
    # List of all possible methods
    methods = [
        'gsph_wendland', 'ssph_wendland', 'disph_wendland', 'gdisph_wendland', 'gdisph_balsara_wendland',
        'gsph_cubic', 'ssph_cubic', 'disph_cubic', 'gdisph_cubic', 'gdisph_balsara_cubic'
    ]
    
    print("=" * 70)
    print("Generating 2D Shock Tube Animations")
    print("=" * 70)
    print()
    
    created = 0
    for method in methods:
        if create_animation_for_method(results_dir, method, animations_dir, args.fps):
            created += 1
    
    print()
    print("=" * 70)
    print(f"Created {created} animations")
    print(f"Output directory: {animations_dir}")
    print("=" * 70)

if __name__ == '__main__':
    main()
