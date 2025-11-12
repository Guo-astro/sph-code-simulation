#!/usr/bin/env python3
"""
Generate animations for vacuum test results
Creates GIF animations from snapshot comparison plots
"""

import os
import sys
from pathlib import Path
import subprocess
import argparse


def create_animation(results_dir, output_file=None, fps=5):
    """
    Create animation from comparison plots
    
    Parameters:
    -----------
    results_dir : str
        Path to results directory containing plots/ subdirectory
    output_file : str, optional
        Output animation file (default: results_dir/vacuum_animation.gif)
    fps : int
        Frames per second (default: 5)
    """
    results_path = Path(results_dir)
    plots_dir = results_path / 'plots'
    
    if not plots_dir.exists():
        print(f"Error: Plots directory {plots_dir} does not exist")
        print("Run process_all_snapshots.py first to generate plots")
        return False
    
    # Find all comparison plots
    plot_files = sorted(plots_dir.glob('snapshot_*_comparison.png'))
    
    if not plot_files:
        print(f"No comparison plots found in {plots_dir}")
        return False
    
    print(f"Found {len(plot_files)} plots for animation")
    
    # Set output file
    if output_file is None:
        output_file = results_path / 'vacuum_animation.gif'
    else:
        output_file = Path(output_file)
    
    # Create animation using ImageMagick convert
    print(f"Creating animation: {output_file}")
    
    # Build command
    cmd = [
        'convert',
        '-delay', str(100 // fps),  # Delay in 1/100 seconds
        '-loop', '0',  # Loop forever
    ] + [str(f) for f in plot_files] + [str(output_file)]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"✓ Animation saved: {output_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error creating animation: {e}")
        print(f"stderr: {e.stderr.decode()}")
        return False
    except FileNotFoundError:
        print("Error: ImageMagick 'convert' command not found")
        print("Install with: brew install imagemagick")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Generate animation from vacuum test comparison plots'
    )
    parser.add_argument('results_dir', type=str,
                       help='Results directory containing plots/ subdirectory')
    parser.add_argument('-o', '--output', type=str,
                       help='Output animation file (default: results_dir/vacuum_animation.gif)')
    parser.add_argument('--fps', type=int, default=5,
                       help='Frames per second (default: 5)')
    
    args = parser.parse_args()
    
    success = create_animation(args.results_dir, args.output, args.fps)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
