#!/usr/bin/env python3
"""
gradh_animate.py - Grad-h diagnostic animation for IMBH cloud simulations.

Creates animated comparison of SPH methods with Lane-Emden analytic overlay.
Uses SSOT modules from scripts.shared for Lane-Emden solutions and snapshot reading.
"""

import sys
from pathlib import Path

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from scripts.shared.lane_emden import get_analytic_profile, IC_PRESETS
from scripts.shared.snapshot import read_snapshot, get_snapshot_files, get_snapshot_time


def create_comparison_animation(
    result_dirs: list[Path],
    method_names: list[str],
    output_path: Path,
    preset_name: str = "imbh_10k",
    fps: int = 10,
    dpi: int = 100,
    show_analytic: bool = True,
    radial_bins: int = 50,
) -> None:
    """
    Create animated comparison of multiple SPH methods.
    
    Parameters
    ----------
    result_dirs : list[Path]
        List of directories containing simulation results
    method_names : list[str]
        Names for each method (for legend)
    output_path : Path
        Output path for the animation (GIF or MP4)
    preset_name : str
        IC preset name from IC_PRESETS (default: "imbh_10k")
    fps : int
        Frames per second
    dpi : int
        Resolution
    show_analytic : bool
        Whether to show Lane-Emden analytic profile
    radial_bins : int
        Number of radial bins for profile averaging
    """
    # Get preset parameters
    if preset_name not in IC_PRESETS:
        raise ValueError(f"Unknown preset: {preset_name}. Available: {list(IC_PRESETS.keys())}")
    
    preset = IC_PRESETS[preset_name]
    rho_c = preset["rho_c"]
    R = preset["R"]
    
    # Get Lane-Emden analytic profile
    if show_analytic:
        r_analytic, rho_analytic = get_analytic_profile(preset_name, n_points=200)
    
    # Collect snapshots from all directories
    all_snapshots = []
    for result_dir in result_dirs:
        snapshots = get_snapshot_files(result_dir, pattern="snapshot_*.csv")
        all_snapshots.append(sorted(snapshots))
    
    # Find common frame count (minimum across all methods)
    n_frames = min(len(s) for s in all_snapshots)
    if n_frames == 0:
        raise ValueError("No snapshots found in result directories")
    
    print(f"Creating animation with {n_frames} frames from {len(method_names)} methods")
    
    # Set up figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Colors for methods
    colors = plt.cm.tab10(np.linspace(0, 1, len(method_names)))
    
    def animate(frame_idx: int) -> list:
        """Update function for animation."""
        for ax in axes:
            ax.clear()
        
        # Get time from first method's snapshot
        time = get_snapshot_time(all_snapshots[0][frame_idx])
        
        # Left panel: Density vs radius
        ax1 = axes[0]
        
        # Plot analytic profile
        if show_analytic:
            ax1.plot(r_analytic, rho_analytic, 'k-', lw=2, label='Lane-Emden (n=1.5)', zorder=10)
        
        # Plot each method
        for i, (snapshots, name, color) in enumerate(zip(all_snapshots, method_names, colors)):
            if frame_idx >= len(snapshots):
                continue
            
            data = read_snapshot(snapshots[frame_idx])
            
            # Compute radial distance from center
            x = data['x'] - np.mean(data['x'])
            y = data['y'] - np.mean(data['y'])
            z = data['z'] - np.mean(data['z']) if 'z' in data else np.zeros_like(x)
            
            r = np.sqrt(x**2 + y**2 + z**2)
            rho = data['rho']
            
            # Scatter plot with transparency
            ax1.scatter(r, rho, s=1, alpha=0.3, color=color, label=name)
            
            # Compute radial profile (binned average)
            r_bins = np.linspace(0, R * 1.2, radial_bins)
            r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
            rho_profile = np.zeros(len(r_centers))
            
            for j in range(len(r_centers)):
                mask = (r >= r_bins[j]) & (r < r_bins[j + 1])
                if np.sum(mask) > 0:
                    rho_profile[j] = np.mean(rho[mask])
            
            # Plot binned profile
            valid = rho_profile > 0
            ax1.plot(r_centers[valid], rho_profile[valid], '-', lw=2, color=color, alpha=0.8)
        
        ax1.set_xlabel('Radius r', fontsize=12)
        ax1.set_ylabel('Density ρ', fontsize=12)
        ax1.set_title(f'Density Profile (t = {time:.4f})', fontsize=14)
        ax1.set_xlim(0, R * 1.2)
        ax1.set_ylim(0, rho_c * 1.2)
        ax1.legend(loc='upper right', fontsize=9)
        ax1.grid(True, alpha=0.3)
        
        # Right panel: Density error relative to analytic
        ax2 = axes[1]
        
        if show_analytic:
            for i, (snapshots, name, color) in enumerate(zip(all_snapshots, method_names, colors)):
                if frame_idx >= len(snapshots):
                    continue
                
                data = read_snapshot(snapshots[frame_idx])
                
                x = data['x'] - np.mean(data['x'])
                y = data['y'] - np.mean(data['y'])
                z = data['z'] - np.mean(data['z']) if 'z' in data else np.zeros_like(x)
                
                r = np.sqrt(x**2 + y**2 + z**2)
                rho = data['rho']
                
                # Interpolate analytic to particle positions
                rho_exact = np.interp(r, r_analytic, rho_analytic, right=0)
                
                # Compute relative error (avoid division by zero)
                valid = rho_exact > 0.01 * rho_c
                rel_error = np.zeros_like(rho)
                rel_error[valid] = (rho[valid] - rho_exact[valid]) / rho_exact[valid]
                
                ax2.scatter(r[valid], rel_error[valid] * 100, s=1, alpha=0.3, color=color, label=name)
        
        ax2.axhline(0, color='k', ls='--', lw=1)
        ax2.set_xlabel('Radius r', fontsize=12)
        ax2.set_ylabel('Relative Error (%)', fontsize=12)
        ax2.set_title('Density Error vs Lane-Emden', fontsize=14)
        ax2.set_xlim(0, R * 1.2)
        ax2.set_ylim(-50, 50)
        ax2.legend(loc='upper right', fontsize=9)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return axes
    
    # Create animation
    anim = FuncAnimation(fig, animate, frames=n_frames, interval=1000 // fps, blit=False)
    
    # Save animation
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if output_path.suffix.lower() == '.gif':
        writer = PillowWriter(fps=fps)
        anim.save(str(output_path), writer=writer, dpi=dpi)
    else:
        anim.save(str(output_path), fps=fps, dpi=dpi)
    
    plt.close(fig)
    print(f"✓ Animation saved to {output_path}")


def create_single_method_animation(
    result_dir: Path,
    output_path: Path,
    preset_name: str = "imbh_10k",
    fps: int = 10,
    dpi: int = 100,
) -> None:
    """
    Create animation for a single SPH method.
    
    Parameters
    ----------
    result_dir : Path
        Directory containing simulation results
    output_path : Path
        Output path for the animation
    preset_name : str
        IC preset name
    fps : int
        Frames per second
    dpi : int
        Resolution
    """
    create_comparison_animation(
        result_dirs=[result_dir],
        method_names=[result_dir.name],
        output_path=output_path,
        preset_name=preset_name,
        fps=fps,
        dpi=dpi,
    )


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Create grad-h diagnostic animations for IMBH cloud simulations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single method animation
  python gradh_animate.py --results results/gsph --output gsph_anim.gif
  
  # Multi-method comparison
  python gradh_animate.py --results results/gsph results/disph results/gdisph \\
                          --names GSPH DISPH GDISPH --output comparison.gif
  
  # Custom preset and FPS
  python gradh_animate.py --results results/gsph --preset imbh_10k --fps 15
        """,
    )
    
    parser.add_argument(
        "--results", "-r",
        nargs="+",
        type=Path,
        required=True,
        help="Result directories (one or more)",
    )
    parser.add_argument(
        "--names", "-n",
        nargs="+",
        type=str,
        default=None,
        help="Method names for legend (default: directory names)",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("gradh_comparison.gif"),
        help="Output animation path (default: gradh_comparison.gif)",
    )
    parser.add_argument(
        "--preset", "-p",
        type=str,
        default="imbh_10k",
        choices=list(IC_PRESETS.keys()),
        help="IC preset name (default: imbh_10k)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=10,
        help="Frames per second (default: 10)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=100,
        help="Resolution (default: 100)",
    )
    parser.add_argument(
        "--no-analytic",
        action="store_true",
        help="Disable Lane-Emden analytic overlay",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=50,
        help="Number of radial bins (default: 50)",
    )
    
    args = parser.parse_args()
    
    # Validate result directories
    for result_dir in args.results:
        if not result_dir.exists():
            print(f"❌ Error: Result directory not found: {result_dir}")
            sys.exit(1)
    
    # Use directory names if method names not provided
    method_names = args.names if args.names else [d.name for d in args.results]
    
    if len(method_names) != len(args.results):
        print(f"❌ Error: Number of names ({len(method_names)}) must match number of result directories ({len(args.results)})")
        sys.exit(1)
    
    # Create animation
    create_comparison_animation(
        result_dirs=args.results,
        method_names=method_names,
        output_path=args.output,
        preset_name=args.preset,
        fps=args.fps,
        dpi=args.dpi,
        show_analytic=not args.no_analytic,
        radial_bins=args.bins,
    )


if __name__ == "__main__":
    main()
