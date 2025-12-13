#!/usr/bin/env python3
"""
Process All Snapshots with Analytical Solution Overlay

Generates analytical comparison plots for all snapshots in a directory.

Usage:
    python3 process_all_snapshots.py RESULTS_DIR [--gamma GAMMA]
"""

import glob
import os
import sys
import argparse
from pathlib import Path

# Import the analytical solver
from strong_shock_analytical import plot_with_analytical


def main():
    parser = argparse.ArgumentParser(description='Process all strong shock snapshots')
    parser.add_argument('results_dir', type=str, help='Directory containing snapshots')
    parser.add_argument('--gamma', type=float, default=1.4, help='Adiabatic index')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results_dir):
        print(f'ERROR: Directory not found: {args.results_dir}')
        sys.exit(1)
    
    # Find all snapshots
    snapshot_files = sorted(glob.glob(f'{args.results_dir}/snapshot_*.csv'))
    
    if not snapshot_files:
        print(f'ERROR: No snapshot files found in {args.results_dir}')
        sys.exit(1)
    
    print('=' * 70)
    print('Processing All Snapshots with Analytical Solution')
    print('=' * 70)
    print(f'Directory: {args.results_dir}')
    print(f'Found {len(snapshot_files)} snapshots')
    print(f'Gamma: {args.gamma}')
    print()
    
    # Create plots directory
    plots_dir = os.path.join(args.results_dir, 'plots_analytical')
    os.makedirs(plots_dir, exist_ok=True)
    
    print(f'Output directory: {plots_dir}')
    print()
    
    # Process each snapshot
    for i, snapshot_file in enumerate(snapshot_files, 1):
        snapshot_name = os.path.basename(snapshot_file)
        output_file = os.path.join(plots_dir, snapshot_name.replace('.csv', '_analytical.png'))
        
        print(f'[{i}/{len(snapshot_files)}] Processing {snapshot_name}...')
        
        try:
            plot_with_analytical(snapshot_file, output_file, show=False, gamma=args.gamma)
        except Exception as e:
            print(f'  ⚠️  Error processing {snapshot_name}: {e}')
            continue
    
    print()
    print('=' * 70)
    print('✓ All Snapshots Processed!')
    print('=' * 70)
    print(f'Output: {plots_dir}/')
    print('=' * 70)


if __name__ == '__main__':
    main()
