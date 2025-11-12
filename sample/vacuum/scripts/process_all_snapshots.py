#!/usr/bin/env python3
"""
Process all snapshots in a directory and generate analytical comparison plots
"""

import os
import sys
import subprocess
from pathlib import Path
import argparse


def process_all_snapshots(results_dir):
    """
    Process all snapshot files in the results directory
    
    Parameters:
    -----------
    results_dir : str
        Path to results directory containing snapshot_*.csv files
    """
    results_path = Path(results_dir)
    
    if not results_path.exists():
        print(f"Error: Directory {results_dir} does not exist")
        return False
    
    # Find all snapshot files
    snapshot_files = sorted(results_path.glob('snapshot_*.csv'))
    
    if not snapshot_files:
        print(f"No snapshot files found in {results_dir}")
        return False
    
    print(f"Found {len(snapshot_files)} snapshots in {results_dir}")
    
    # Create plots directory
    plots_dir = results_path / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    # Get path to analytical script
    script_dir = Path(__file__).parent
    analytical_script = script_dir / 'vacuum_analytical.py'
    
    if not analytical_script.exists():
        print(f"Error: Analytical script not found at {analytical_script}")
        return False
    
    # Process each snapshot
    for i, snapshot_file in enumerate(snapshot_files, 1):
        snapshot_name = snapshot_file.stem
        output_file = plots_dir / f'{snapshot_name}_comparison.png'
        
        print(f"[{i}/{len(snapshot_files)}] Processing {snapshot_file.name}...")
        
        # Run analytical comparison script
        cmd = [
            'python3',
            str(analytical_script),
            str(snapshot_file),
            '-o', str(output_file),
            '--no-show'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"  Error processing {snapshot_file.name}: {e}")
            print(f"  stdout: {e.stdout.decode()}")
            print(f"  stderr: {e.stderr.decode()}")
            continue
    
    print(f"\n✓ All plots saved to {plots_dir}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Process all vacuum test snapshots and generate comparison plots'
    )
    parser.add_argument('results_dir', type=str, 
                       help='Results directory containing snapshot_*.csv files')
    
    args = parser.parse_args()
    
    success = process_all_snapshots(args.results_dir)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
