#!/usr/bin/env python3
"""
2D Neutron Star Merger Visualization

Creates density plots and animations for 2D NS merger simulations.
Based on Shibata & Hotokezaka (2019) setup.

Usage:
    python ns_merger_visualize.py --results-dir results/test_quick --plot
    python ns_merger_visualize.py --results-dir results/test_quick --animate
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation, PillowWriter

# Add scripts directory to path for unit system
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / 'scripts'))

try:
    from units.unit_system import RelativisticUnits
    UNIT_SYSTEM_AVAILABLE = True
except ImportError:
    UNIT_SYSTEM_AVAILABLE = False
    print("Warning: Unit system not available, using code units")


class NSMergerVisualizer:
    """Visualizer for 2D neutron star merger simulations."""
    
    def __init__(self, results_dir: str, config_file: str = None):
        self.results_dir = Path(results_dir).resolve()  # Use absolute path
        self.config_file = config_file
        self.snapshots = []
        self.times = []
        self.units = None
        
        # Load snapshots
        self._load_snapshots()
        
        # Try to load unit system
        if config_file and UNIT_SYSTEM_AVAILABLE:
            import json
            try:
                with open(config_file) as f:
                    config = json.load(f)
                self.units = RelativisticUnits.from_config(config)
            except Exception as e:
                print(f"Warning: Could not load units from config: {e}")
    
    def _load_snapshots(self):
        """Load all snapshot CSV files."""
        pattern = str(self.results_dir / "snapshot_*.csv")
        files = sorted(glob.glob(pattern))
        
        if not files:
            raise FileNotFoundError(f"No snapshot files found in {self.results_dir}")
        
        print(f"Found {len(files)} snapshots")
        
        for f in files:
            df = pd.read_csv(f, comment='#')
            self.snapshots.append(df)
            
            # Extract time from header - look for multiple formats
            time_found = False
            with open(f) as fp:
                for line in fp:
                    if line.startswith('# time='):
                        time = float(line.split('=')[1].strip())
                        self.times.append(time)
                        time_found = True
                        break
                    elif line.startswith('# Time (code):'):
                        time = float(line.split(':')[1].strip())
                        self.times.append(time)
                        time_found = True
                        break
            
            if not time_found:
                # Fallback: use snapshot index
                idx = len(self.times)
                self.times.append(idx * 0.1)  # Approximate
                print(f"Warning: Could not find time in {f}, using index")
        
        print(f"Time range: {self.times[0]:.3f} to {self.times[-1]:.3f}")
    
    def plot_snapshot(self, index: int = -1, save_path: str = None, 
                      show: bool = True, cmap: str = 'magma'):
        """Plot a single snapshot with density coloring."""
        df = self.snapshots[index]
        time = self.times[index]
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        x = df['pos_x'].values
        y = df['pos_y'].values
        density = df['dens'].values
        
        # Log scale for density
        vmin = density[density > 0].min()
        vmax = density.max()
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        
        scatter = ax.scatter(x, y, c=density, s=5, cmap=cmap, norm=norm, 
                            alpha=0.8, edgecolors='none')
        
        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, label='Density')
        
        # Labels
        if self.units:
            ax.set_xlabel(f'x [{self.units.length_label}]')
            ax.set_ylabel(f'y [{self.units.length_label}]')
            time_label = f't = {time:.2f} [{self.units.time_label}]'
        else:
            ax.set_xlabel('x [code units]')
            ax.set_ylabel('y [code units]')
            time_label = f't = {time:.3f}'
        
        ax.set_title(f'2D Neutron Star Merger\n{time_label}')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Domain
        ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
        ax.set_ylim(y.min() - 0.5, y.max() + 0.5)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()
        
        return fig, ax
    
    def create_animation(self, output_path: str = 'ns_merger.gif', 
                        fps: int = 10, cmap: str = 'magma'):
        """Create animation of the merger."""
        print(f"Creating animation with {len(self.snapshots)} frames...")
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Get data ranges from all snapshots
        all_x = np.concatenate([df['pos_x'].values for df in self.snapshots])
        all_y = np.concatenate([df['pos_y'].values for df in self.snapshots])
        all_density = np.concatenate([df['dens'].values for df in self.snapshots])
        
        x_range = (all_x.min() - 0.5, all_x.max() + 0.5)
        y_range = (all_y.min() - 0.5, all_y.max() + 0.5)
        vmin = all_density[all_density > 0].min()
        vmax = all_density.max()
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        
        # Initial plot
        df = self.snapshots[0]
        scatter = ax.scatter(df['pos_x'], df['pos_y'], c=df['dens'], s=5, 
                            cmap=cmap, norm=norm, alpha=0.8, edgecolors='none')
        
        cbar = plt.colorbar(scatter, ax=ax, label='Density')
        
        if self.units:
            ax.set_xlabel(f'x [{self.units.length_label}]')
            ax.set_ylabel(f'y [{self.units.length_label}]')
        else:
            ax.set_xlabel('x [code units]')
            ax.set_ylabel('y [code units]')
        
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        title = ax.set_title('')
        
        def update(frame):
            df = self.snapshots[frame]
            time = self.times[frame]
            
            # Update scatter data
            scatter.set_offsets(np.column_stack([df['pos_x'], df['pos_y']]))
            scatter.set_array(df['dens'])
            
            if self.units:
                title.set_text(f'2D Neutron Star Merger\nt = {time:.2f} [{self.units.time_label}]')
            else:
                title.set_text(f'2D Neutron Star Merger\nt = {time:.3f}')
            
            return scatter, title
        
        anim = FuncAnimation(fig, update, frames=len(self.snapshots),
                            interval=1000/fps, blit=False)
        
        # Save animation
        writer = PillowWriter(fps=fps)
        anim.save(output_path, writer=writer)
        print(f"Animation saved: {output_path}")
        
        plt.close()
    
    def plot_density_evolution(self, save_path: str = None, show: bool = True):
        """Plot density profile evolution along x-axis."""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Select a few representative snapshots
        n_plots = min(8, len(self.snapshots))
        indices = np.linspace(0, len(self.snapshots)-1, n_plots, dtype=int)
        
        colors_list = plt.cm.viridis(np.linspace(0, 1, n_plots))
        
        for i, idx in enumerate(indices):
            df = self.snapshots[idx]
            time = self.times[idx]
            
            # Select particles near y=0
            mask = np.abs(df['pos_y']) < 0.2
            x = df.loc[mask, 'pos_x'].values
            density = df.loc[mask, 'dens'].values
            
            # Sort by x
            order = np.argsort(x)
            x, density = x[order], density[order]
            
            label = f't = {time:.2f}' if self.units else f't = {time:.3f}'
            ax.plot(x, density, 'o-', ms=2, alpha=0.7, color=colors_list[i], label=label)
        
        ax.set_xlabel('x [code units]')
        ax.set_ylabel('Density')
        ax.set_yscale('log')
        ax.set_title('Density Profile Evolution (y ≈ 0)')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_grid(self, n_cols: int = 4, save_path: str = None, show: bool = True):
        """Plot grid of snapshots showing merger evolution."""
        n_total = len(self.snapshots)
        n_plots = min(16, n_total)
        indices = np.linspace(0, n_total-1, n_plots, dtype=int)
        
        n_rows = (n_plots + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
        axes = axes.flatten()
        
        # Get global density range
        all_density = np.concatenate([df['dens'].values for df in self.snapshots])
        vmin = all_density[all_density > 0].min()
        vmax = all_density.max()
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        
        for i, (ax, idx) in enumerate(zip(axes, indices)):
            df = self.snapshots[idx]
            time = self.times[idx]
            
            scatter = ax.scatter(df['pos_x'], df['pos_y'], c=df['dens'], s=1, 
                                cmap='magma', norm=norm, alpha=0.8, edgecolors='none')
            
            ax.set_aspect('equal')
            ax.set_title(f't = {time:.2f}', fontsize=10)
            ax.set_xlim(-6, 6)
            ax.set_ylim(-6, 6)
            ax.tick_params(labelsize=8)
        
        # Hide unused axes
        for ax in axes[len(indices):]:
            ax.set_visible(False)
        
        fig.suptitle('2D Neutron Star Merger Evolution', fontsize=14)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()


def main():
    parser = argparse.ArgumentParser(description='NS Merger Visualization')
    parser.add_argument('--results-dir', '-r', default='results/test_quick',
                       help='Results directory containing snapshots')
    parser.add_argument('--config', '-c', default=None,
                       help='Config JSON file for unit system')
    parser.add_argument('--plot', '-p', action='store_true',
                       help='Plot final snapshot')
    parser.add_argument('--animate', '-a', action='store_true',
                       help='Create animation')
    parser.add_argument('--grid', '-g', action='store_true',
                       help='Plot evolution grid')
    parser.add_argument('--evolution', '-e', action='store_true',
                       help='Plot density evolution')
    parser.add_argument('--output', '-o', default=None,
                       help='Output filename')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not show plots (save only)')
    
    args = parser.parse_args()
    
    # Use current working directory for relative paths
    if not os.path.isabs(args.results_dir):
        results_dir = Path.cwd() / args.results_dir
    else:
        results_dir = Path(args.results_dir)
    
    viz = NSMergerVisualizer(str(results_dir), args.config)
    
    if args.plot:
        output = args.output or str(results_dir / 'ns_merger_final.png')
        viz.plot_snapshot(-1, save_path=output, show=not args.no_show)
    
    if args.animate:
        output = args.output or str(results_dir / 'ns_merger.gif')
        viz.create_animation(output)
    
    if args.grid:
        output = args.output or str(results_dir / 'ns_merger_grid.png')
        viz.plot_grid(save_path=output, show=not args.no_show)
    
    if args.evolution:
        output = args.output or str(results_dir / 'ns_merger_evolution.png')
        viz.plot_density_evolution(save_path=output, show=not args.no_show)
    
    if not any([args.plot, args.animate, args.grid, args.evolution]):
        print("No action specified. Use --plot, --animate, --grid, or --evolution")
        parser.print_help()


if __name__ == '__main__':
    main()
