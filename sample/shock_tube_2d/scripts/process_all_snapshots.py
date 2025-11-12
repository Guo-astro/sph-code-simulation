#!/usr/bin/env python3
"""
Process all 2D shock tube snapshots and create comparison plots with analytical solution
"""

import os
import sys
from pathlib import Path
import glob
from shock_tube_2d_analytical import plot_shock_tube_2d_comparison

def main():
    """Process all snapshots in a results directory."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Process all 2D shock tube snapshots')
    parser.add_argument('results_dir', help='Path to results directory containing snapshots')
    parser.add_argument('-o', '--output-dir', help='Output directory for plots (default: results_dir/plots)')
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        print(f"Error: Directory not found: {results_dir}")
        return 1
    
    # Output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = results_dir / 'plots'
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all snapshot files
    snapshots = sorted(results_dir.glob('snapshot_*.csv'))
    
    if not snapshots:
        print(f"No snapshot files found in {results_dir}")
        return 1
    
    print("=" * 70)
    print(f"Processing {len(snapshots)} snapshots from {results_dir}")
    print(f"Output directory: {output_dir}")
    print("=" * 70)
    print()
    
    # Process each snapshot
    for i, snapshot in enumerate(snapshots):
        snapshot_name = snapshot.stem
        output_file = output_dir / f'{snapshot_name}_comparison.png'
        
        print(f"[{i+1}/{len(snapshots)}] Processing {snapshot.name}...")
        
        try:
            plot_shock_tube_2d_comparison(str(snapshot), str(output_file), show_plot=False)
        except Exception as e:
            print(f"  ❌ Error: {e}")
            continue
    
    print()
    print("=" * 70)
    print("Processing complete!")
    print(f"Plots saved to: {output_dir}")
    print("=" * 70)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
