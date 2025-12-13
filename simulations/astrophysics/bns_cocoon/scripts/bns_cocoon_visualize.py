#!/usr/bin/env python3
"""
Visualization script for BNS cocoon shock breakout simulations.

Generates:
1. Density and Lorentz factor profiles
2. Shock position tracking plots
3. Observable comparison with GRB 170817A
4. Animation of shock propagation
"""

import os
import sys
import argparse
import glob
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.colors as mcolors


# Add parent directory to path for local imports
sys.path.insert(0, str(Path(__file__).parent))
from analysis.shock_tracker import ShockTracker
from analysis.optical_depth import OpticalDepthCalculator
from analysis.observables import BreakoutObservables


@dataclass
class PlotConfig:
    """Configuration for plot styling."""
    figsize: tuple = (10, 8)
    dpi: int = 150
    cmap: str = "viridis"
    style: str = "seaborn-v0_8-whitegrid"
    
    # GRB 170817A reference values
    grb170817a_E_iso: float = 5e46  # erg
    grb170817a_T90: float = 2.0     # s
    grb170817a_Epeak: float = 185   # keV


class CocoonVisualizer:
    """Main visualization class for cocoon shock breakout."""
    
    def __init__(self, results_dir: str | Path, config: PlotConfig | None = None):
        """
        Initialize visualizer.
        
        Parameters
        ----------
        results_dir : str or Path
            Directory containing simulation output CSV files.
        config : PlotConfig, optional
            Plot configuration.
        """
        self.results_dir = Path(results_dir)
        self.config = config or PlotConfig()
        
        # Find snapshot files
        self.snapshots = sorted(glob.glob(str(self.results_dir / "snapshot_*.csv")))
        if not self.snapshots:
            raise ValueError(f"No snapshot files found in {results_dir}")
        
        print(f"Found {len(self.snapshots)} snapshots in {results_dir}")
        
        # Initialize analysis tools
        self.shock_tracker = ShockTracker(snapshot_dir=str(results_dir))
        self.optical_depth = OpticalDepthCalculator()
        self.observables = BreakoutObservables()
    
    def load_snapshot(self, filepath: str) -> tuple[np.ndarray, float]:
        """
        Load a snapshot CSV file.
        
        Returns
        -------
        data : structured ndarray
            Particle data with columns: pos_x, vel_x, dens, pres, etc.
        time : float
            Simulation time.
        """
        # Read header to extract time
        time = 0.0
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    if 'Time' in line:
                        parts = line.split(':')
                        if len(parts) >= 2:
                            try:
                                time = float(parts[1].strip())
                            except ValueError:
                                pass
                else:
                    break
        
        # Load data
        data = np.genfromtxt(filepath, delimiter=',', names=True, comments='#')
        
        return data, time
    
    def plot_radial_profile(
        self,
        snapshot_idx: int = -1,
        quantities: list[str] = ['dens', 'lorentz'],
        output_file: str | None = None,
        show: bool = True
    ):
        """
        Plot radial profiles of density and Lorentz factor.
        
        Parameters
        ----------
        snapshot_idx : int
            Index of snapshot to plot (-1 for last).
        quantities : list of str
            Quantities to plot: 'dens', 'pres', 'lorentz', 'tau'.
        output_file : str, optional
            Save plot to file.
        show : bool
            Show plot interactively.
        """
        filepath = self.snapshots[snapshot_idx]
        data, time = self.load_snapshot(filepath)
        
        # Get radial coordinate
        if 'pos_x' in data.dtype.names:
            r = np.abs(data['pos_x'])
        elif 'pos_r' in data.dtype.names:
            r = data['pos_r']
        else:
            r = np.sqrt(data['pos_x']**2 + data['pos_y']**2)
        
        # Sort by radius
        idx = np.argsort(r)
        r = r[idx]
        
        # Set up figure
        n_plots = len(quantities)
        fig, axes = plt.subplots(n_plots, 1, figsize=(10, 4*n_plots), sharex=True)
        if n_plots == 1:
            axes = [axes]
        
        for ax, qty in zip(axes, quantities):
            if qty == 'dens' and 'dens' in data.dtype.names:
                y = data['dens'][idx]
                ax.semilogy(r, y, 'b-', lw=1.5)
                ax.set_ylabel(r'Density $\rho$ [code units]')
                
            elif qty == 'pres' and 'pres' in data.dtype.names:
                y = data['pres'][idx]
                ax.semilogy(r, y, 'r-', lw=1.5)
                ax.set_ylabel(r'Pressure $P$ [code units]')
                
            elif qty == 'lorentz' and 'vel_x' in data.dtype.names:
                v = np.abs(data['vel_x'][idx])
                gamma = 1.0 / np.sqrt(1 - v**2)
                ax.plot(r, gamma, 'g-', lw=1.5)
                ax.set_ylabel(r'Lorentz factor $\Gamma$')
                ax.set_yscale('log')
                
            elif qty == 'tau':
                # Compute optical depth profile
                dens = data['dens'][idx]
                tau = self.optical_depth.compute_tau_profile(r, dens)
                ax.semilogy(r, tau, 'm-', lw=1.5)
                ax.axhline(1, color='k', ls='--', alpha=0.5, label=r'$\tau = 1$')
                ax.set_ylabel(r'Optical depth $\tau$')
                ax.legend()
            
            ax.grid(True, alpha=0.3)
        
        axes[-1].set_xlabel(r'Radius $r$ [code units]')
        axes[0].set_title(f'Radial profiles at t = {time:.4f}')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_shock_evolution(
        self,
        output_file: str | None = None,
        show: bool = True
    ):
        """
        Plot shock position and Lorentz factor evolution.
        
        Parameters
        ----------
        output_file : str, optional
            Save plot to file.
        show : bool
            Show plot interactively.
        """
        # Track shock through all snapshots
        times, r_shock, gamma_shock, E_int = [], [], [], []
        
        for filepath in self.snapshots:
            data, time = self.load_snapshot(filepath)
            
            # Get radius and velocity
            r = np.abs(data['pos_x'])
            v = np.abs(data['vel_x']) if 'vel_x' in data.dtype.names else np.zeros_like(r)
            dens = data['dens']
            
            # Find shock position (density jump)
            idx_sort = np.argsort(r)
            r_sorted = r[idx_sort]
            dens_sorted = dens[idx_sort]
            
            # Compute density gradient
            grad = np.gradient(np.log10(dens_sorted + 1e-30), r_sorted)
            i_shock = np.argmax(np.abs(grad))
            
            r_sh = r_sorted[i_shock]
            v_sh = np.abs(data['vel_x'][idx_sort][i_shock]) if 'vel_x' in data.dtype.names else 0.1
            gamma_sh = 1.0 / np.sqrt(1 - min(v_sh**2, 0.9999))
            
            times.append(time)
            r_shock.append(r_sh)
            gamma_shock.append(gamma_sh)
        
        times = np.array(times)
        r_shock = np.array(r_shock)
        gamma_shock = np.array(gamma_shock)
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # Shock radius
        ax1.plot(times, r_shock, 'b-o', lw=2, markersize=4)
        ax1.set_ylabel(r'Shock radius $r_{\rm sh}$ [code units]')
        ax1.grid(True, alpha=0.3)
        ax1.set_title('Shock Evolution')
        
        # Lorentz factor
        ax2.semilogy(times, gamma_shock, 'r-o', lw=2, markersize=4)
        ax2.set_ylabel(r'Shock Lorentz factor $\Gamma_{\rm sh}$')
        ax2.set_xlabel('Time [code units]')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_observables(
        self,
        output_file: str | None = None,
        show: bool = True
    ):
        """
        Plot comparison with GRB 170817A observables.
        
        Parameters
        ----------
        output_file : str, optional
            Save plot to file.
        show : bool
            Show plot interactively.
        """
        # Analyze final snapshot for breakout observables
        data, time = self.load_snapshot(self.snapshots[-1])
        
        r = np.abs(data['pos_x'])
        v = np.abs(data['vel_x']) if 'vel_x' in data.dtype.names else np.zeros_like(r)
        dens = data['dens']
        
        # Find shock/breakout
        idx_sort = np.argsort(r)
        r_sorted = r[idx_sort]
        dens_sorted = dens[idx_sort]
        
        # Optical depth
        tau = self.optical_depth.compute_tau_profile(r_sorted, dens_sorted)
        
        # Find breakout radius (τ = 1)
        i_bo = np.searchsorted(-tau, -1)  # Find where τ drops below 1
        if i_bo < len(r_sorted):
            r_bo = r_sorted[i_bo]
            gamma_bo = 1.0 / np.sqrt(1 - min(v[idx_sort][i_bo]**2, 0.9999))
        else:
            r_bo = r_sorted[-1]
            gamma_bo = 1.0
        
        # Compute observables (simplified)
        # In reality, would use full breakout model
        T_obs = self.observables.compute_temperature(gamma_bo, r_bo)
        E_rad = self.observables.compute_radiated_energy(gamma_bo, r_bo)
        
        # Create comparison plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: τ profile
        ax = axes[0, 0]
        ax.semilogy(r_sorted, tau, 'b-', lw=2)
        ax.axhline(1, color='r', ls='--', label=r'$\tau = 1$ (breakout)')
        ax.axvline(r_bo, color='g', ls=':', label=f'$r_{{bo}} = {r_bo:.2e}$')
        ax.set_xlabel('Radius')
        ax.set_ylabel(r'Optical depth $\tau$')
        ax.legend()
        ax.set_title('Optical Depth Profile')
        ax.grid(True, alpha=0.3)
        
        # Panel 2: Lorentz factor profile
        ax = axes[0, 1]
        gamma = 1.0 / np.sqrt(1 - np.minimum(v**2, 0.9999))
        ax.semilogy(r[idx_sort], gamma[idx_sort], 'r-', lw=2)
        ax.axvline(r_bo, color='g', ls=':', label='Breakout')
        ax.axhline(gamma_bo, color='k', ls='--', alpha=0.5)
        ax.set_xlabel('Radius')
        ax.set_ylabel(r'Lorentz factor $\Gamma$')
        ax.legend()
        ax.set_title(f'Lorentz Factor (Γ_bo = {gamma_bo:.2f})')
        ax.grid(True, alpha=0.3)
        
        # Panel 3: E_iso comparison
        ax = axes[1, 0]
        E_model = [E_rad]  # Our model
        E_obs = [self.config.grb170817a_E_iso]
        x = [0, 1]
        colors = ['blue', 'orange']
        labels = ['This simulation', 'GRB 170817A']
        
        bars = ax.bar(x, [E_rad, self.config.grb170817a_E_iso], color=colors)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel(r'$E_{\rm iso}$ [erg]')
        ax.set_yscale('log')
        ax.set_title('Isotropic Energy Comparison')
        
        # Add value labels on bars
        for bar, val in zip(bars, [E_rad, self.config.grb170817a_E_iso]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                   f'{val:.1e}', ha='center', va='bottom')
        
        # Panel 4: Summary text
        ax = axes[1, 1]
        ax.axis('off')
        
        summary = [
            r"$\bf{Breakout\ Parameters:}$",
            f"  Radius: $r_{{bo}} = {r_bo:.2e}$ cm",
            f"  Lorentz factor: $\\Gamma_{{bo}} = {gamma_bo:.2f}$",
            f"  Temperature: $T_{{obs}} = {T_obs:.1f}$ keV",
            f"  Radiated energy: $E_{{rad}} = {E_rad:.2e}$ erg",
            "",
            r"$\bf{GRB\ 170817A:}$",
            f"  $E_{{\\gamma,iso}} = {self.config.grb170817a_E_iso:.1e}$ erg",
            f"  $T_{{90}} = {self.config.grb170817a_T90}$ s",
            f"  $E_{{peak}} = {self.config.grb170817a_Epeak}$ keV",
        ]
        
        ax.text(0.1, 0.9, '\n'.join(summary), transform=ax.transAxes,
               fontsize=12, verticalalignment='top', family='monospace')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def create_animation(
        self,
        output_file: str = "cocoon_evolution.gif",
        fps: int = 10,
        show: bool = False
    ):
        """
        Create animation of shock propagation.
        
        Parameters
        ----------
        output_file : str
            Output GIF filename.
        fps : int
            Frames per second.
        show : bool
            Show animation interactively before saving.
        """
        # Set up figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Load first snapshot for initialization
        data0, time0 = self.load_snapshot(self.snapshots[0])
        r0 = np.abs(data0['pos_x'])
        r_max = r0.max() * 2  # Buffer for expansion
        
        # Initialize lines
        line_dens, = ax1.semilogy([], [], 'b-', lw=2, label='Density')
        line_gamma, = ax2.plot([], [], 'r-', lw=2, label='Lorentz factor')
        
        ax1.set_xlim(0, r_max)
        ax1.set_ylim(1e-10, data0['dens'].max() * 10)
        ax1.set_xlabel('Radius')
        ax1.set_ylabel(r'Density $\rho$')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        ax2.set_xlim(0, r_max)
        ax2.set_ylim(0.9, 10)
        ax2.set_xlabel('Radius')
        ax2.set_ylabel(r'Lorentz factor $\Gamma$')
        ax2.set_yscale('log')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        title = fig.suptitle('')
        
        def init():
            line_dens.set_data([], [])
            line_gamma.set_data([], [])
            return line_dens, line_gamma
        
        def animate(i):
            data, time = self.load_snapshot(self.snapshots[i])
            
            r = np.abs(data['pos_x'])
            idx = np.argsort(r)
            r = r[idx]
            dens = data['dens'][idx]
            
            if 'vel_x' in data.dtype.names:
                v = np.abs(data['vel_x'][idx])
                gamma = 1.0 / np.sqrt(1 - np.minimum(v**2, 0.9999))
            else:
                gamma = np.ones_like(r)
            
            line_dens.set_data(r, dens)
            line_gamma.set_data(r, gamma)
            title.set_text(f'Time = {time:.4f} [code units]')
            
            return line_dens, line_gamma
        
        anim = FuncAnimation(
            fig, animate, init_func=init,
            frames=len(self.snapshots),
            interval=1000//fps, blit=True
        )
        
        # Save
        output_path = self.results_dir / output_file
        print(f"Creating animation with {len(self.snapshots)} frames...")
        anim.save(str(output_path), writer=PillowWriter(fps=fps))
        print(f"Saved: {output_path}")
        
        if show:
            plt.show()
        else:
            plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Visualize BNS cocoon shock breakout simulations"
    )
    
    parser.add_argument(
        "results_dir",
        help="Directory containing simulation snapshot CSV files"
    )
    parser.add_argument(
        "--profile", "-p", action="store_true",
        help="Plot radial profiles"
    )
    parser.add_argument(
        "--shock", "-s", action="store_true",
        help="Plot shock evolution"
    )
    parser.add_argument(
        "--observables", "-o", action="store_true",
        help="Plot GRB 170817A comparison"
    )
    parser.add_argument(
        "--animate", "-a", action="store_true",
        help="Create evolution animation"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Generate all plots and animation"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Output directory for plots (default: results_dir)"
    )
    parser.add_argument(
        "--no-show", action="store_true",
        help="Don't display plots interactively"
    )
    
    args = parser.parse_args()
    
    # Create visualizer
    viz = CocoonVisualizer(args.results_dir)
    
    output_dir = args.output_dir or Path(args.results_dir)
    show = not args.no_show
    
    if args.all or args.profile:
        viz.plot_radial_profile(
            output_file=str(output_dir / "radial_profiles.png"),
            show=show
        )
    
    if args.all or args.shock:
        viz.plot_shock_evolution(
            output_file=str(output_dir / "shock_evolution.png"),
            show=show
        )
    
    if args.all or args.observables:
        viz.plot_observables(
            output_file=str(output_dir / "observables_comparison.png"),
            show=show
        )
    
    if args.all or args.animate:
        viz.create_animation(
            output_file="cocoon_evolution.gif",
            show=show
        )


if __name__ == "__main__":
    main()
