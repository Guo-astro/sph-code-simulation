#!/usr/bin/env python3
"""
Generate continuation config from an existing simulation config and snapshot.

This script creates a new config file for continuing an interrupted or completed
simulation from a specific snapshot. It:
1. Copies all physics parameters from the original config
2. Sets resumeFromSnapshot to the specified snapshot path
3. Sets resetTimeOnResume to false (continue from snapshot time)
4. Updates endTime to the new target time

Usage:
    python generate_continue_config.py \\
        --original-config config/presets/simulation_10k_b3pc_nocool.json \\
        --snapshot results/imbh_relaxed_2k_b3pc_nocool/snapshot_0100.csv \\
        --end-time 4.0 \\
        --output config/continue/simulation_10k_b3pc_nocool_from_snap0100.json

    # Or with shorthand:
    python generate_continue_config.py -c config/presets/simulation_10k_b3pc_nocool.json \\
        -s results/imbh_relaxed_2k_b3pc_nocool/snapshot_0100.csv -e 4.0

Author: Auto-generated for IMBH cloud simulation workflow
Date: 2024-12-02
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime


def extract_snapshot_number(snapshot_path: str) -> int:
    """Extract snapshot number from filename like snapshot_0100.csv"""
    match = re.search(r'snapshot_(\d+)', snapshot_path)
    if match:
        return int(match.group(1))
    return -1


def estimate_snapshot_time(snapshot_path: str, output_time: float) -> float:
    """Estimate simulation time from snapshot number and output interval"""
    snap_num = extract_snapshot_number(snapshot_path)
    if snap_num >= 0:
        return snap_num * output_time
    return -1.0


def generate_continue_config(
    original_config_path: str,
    snapshot_path: str,
    new_end_time: float,
    output_path: str = None,
    output_directory: str = None,
    description: str = None,
    dry_run: bool = False
) -> dict:
    """
    Generate a continuation config from original config and snapshot.
    
    Args:
        original_config_path: Path to the original simulation config JSON
        snapshot_path: Path to the snapshot file to continue from
        new_end_time: New simulation end time
        output_path: Path for the new config file (auto-generated if None)
        output_directory: Override output directory (uses original if None)
        description: Optional description for the continuation
        dry_run: If True, only print what would be done without writing
    
    Returns:
        The generated config dictionary
    """
    # Load original config
    with open(original_config_path, 'r') as f:
        config = json.load(f)
    
    # Extract info for naming and documentation
    original_name = config.get('name', os.path.basename(original_config_path).replace('.json', ''))
    snap_num = extract_snapshot_number(snapshot_path)
    output_time = config.get('outputTime', 0.02)
    estimated_time = estimate_snapshot_time(snapshot_path, output_time)
    
    # Generate output path if not provided
    if output_path is None:
        base_dir = os.path.dirname(original_config_path)
        continue_dir = os.path.join(os.path.dirname(base_dir), 'continue')
        snap_str = f"snap{snap_num:04d}" if snap_num >= 0 else "snap"
        output_path = os.path.join(continue_dir, f"{original_name}_from_{snap_str}.json")
    
    # Update config for continuation
    config['name'] = f"{original_name}_continue_from_snap{snap_num:04d}" if snap_num >= 0 else f"{original_name}_continue"
    
    # Add/update description
    if description:
        config['description'] = description
    else:
        config['description'] = f"CONTINUE simulation from snapshot_{snap_num:04d}"
    
    # Add continuation metadata
    config['continuation_info'] = {
        'original_config': original_config_path,
        'source_snapshot': snapshot_path,
        'snapshot_time': estimated_time if estimated_time >= 0 else 'unknown',
        'created_date': datetime.now().strftime('%Y-%m-%d'),
        'reason': 'Continue to extended end time'
    }
    
    # KEY CHANGES for continuation
    config['resumeFromSnapshot'] = snapshot_path
    config['resetTimeOnResume'] = False  # CRITICAL: continue from snapshot time
    config['endTime'] = new_end_time
    
    # Update output directory if specified
    if output_directory:
        config['outputDirectory'] = output_directory
    
    # Ensure checkpoint is enabled
    if 'checkpoint' not in config:
        config['checkpoint'] = {}
    config['checkpoint']['enabled'] = True
    config['checkpoint']['saveOnExit'] = True
    
    # Print summary
    print("=" * 60)
    print("CONTINUATION CONFIG GENERATOR")
    print("=" * 60)
    print(f"Original config: {original_config_path}")
    print(f"Source snapshot: {snapshot_path}")
    print(f"Snapshot number: {snap_num}")
    print(f"Estimated time:  {estimated_time:.2f}" if estimated_time >= 0 else "Estimated time:  unknown")
    print(f"New end time:    {new_end_time}")
    print(f"Output path:     {output_path}")
    print("-" * 60)
    print("KEY SETTINGS:")
    print(f"  resumeFromSnapshot: {snapshot_path}")
    print("  resetTimeOnResume:  False  (continues from snapshot time)")
    print(f"  endTime:            {new_end_time}")
    print("=" * 60)
    
    if dry_run:
        print("\n[DRY RUN] Would create config:")
        print(json.dumps(config, indent=2))
        return config
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write config
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✓ Config written to: {output_path}")
    print("\nTo run continuation:")
    print(f"  ./build/sph {output_path}")
    
    return config


def main():
    parser = argparse.ArgumentParser(
        description='Generate continuation config for SPH simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - extend simulation to t=4.0
  python generate_continue_config.py \\
      -c config/presets/simulation_10k_b3pc_nocool.json \\
      -s results/imbh_relaxed_2k_b3pc_nocool/snapshot_0100.csv \\
      -e 4.0

  # Custom output path
  python generate_continue_config.py \\
      -c config/presets/simulation_10k_b3pc_nocool.json \\
      -s results/imbh_relaxed_2k_b3pc_nocool/snapshot_0100.csv \\
      -e 4.0 \\
      -o config/continue/my_custom_name.json

  # Different output directory
  python generate_continue_config.py \\
      -c config/presets/simulation_10k_b3pc_nocool.json \\
      -s results/imbh_relaxed_2k_b3pc_nocool/snapshot_0100.csv \\
      -e 4.0 \\
      --output-dir results/imbh_extended_run/

  # Dry run (preview only)
  python generate_continue_config.py \\
      -c config/presets/simulation_10k_b3pc_nocool.json \\
      -s results/imbh_relaxed_2k_b3pc_nocool/snapshot_0100.csv \\
      -e 4.0 \\
      --dry-run
"""
    )
    
    parser.add_argument(
        '-c', '--original-config',
        required=True,
        help='Path to original simulation config JSON'
    )
    parser.add_argument(
        '-s', '--snapshot',
        required=True,
        help='Path to snapshot file to continue from'
    )
    parser.add_argument(
        '-e', '--end-time',
        type=float,
        required=True,
        help='New simulation end time'
    )
    parser.add_argument(
        '-o', '--output',
        default=None,
        help='Output path for continuation config (auto-generated if not specified)'
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Override output directory for simulation results'
    )
    parser.add_argument(
        '-d', '--description',
        default=None,
        help='Description for the continuation run'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Preview config without writing file'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.original_config):
        print(f"ERROR: Original config not found: {args.original_config}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.snapshot):
        print(f"ERROR: Snapshot not found: {args.snapshot}", file=sys.stderr)
        sys.exit(1)
    
    # Generate config
    generate_continue_config(
        original_config_path=args.original_config,
        snapshot_path=args.snapshot,
        new_end_time=args.end_time,
        output_path=args.output,
        output_directory=args.output_dir,
        description=args.description,
        dry_run=args.dry_run
    )


if __name__ == '__main__':
    main()
