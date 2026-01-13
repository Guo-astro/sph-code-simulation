#!/usr/bin/env python3
"""
Convert SPH snapshot CSV files to compact binary format for web hosting.
Reduces file size by ~5x (11MB CSV → ~2MB binary).

Binary format per snapshot:
- Header: 4 bytes (uint32 particle count)
- Per particle: 40 bytes (10 × float32)
  - x, y, z (position)
  - vx, vy, vz (velocity)
  - dens (density)
  - temp (temperature)
  - mass
  - sound (sound speed)
"""

import numpy as np
import pandas as pd
import struct
import os
import glob
import argparse
from pathlib import Path


def convert_csv_to_binary(csv_path: str, bin_path: str):
    """Convert a single CSV snapshot to binary format."""
    df = pd.read_csv(csv_path, comment='#')

    # Required columns (matching current visualization)
    columns = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'dens', 'temp', 'mass', 'sound']

    # Filter out ghost particles if present
    if 'is_ghost' in df.columns:
        df = df[df['is_ghost'] != 1]

    n_particles = len(df)

    # Write binary file
    with open(bin_path, 'wb') as f:
        # Write particle count (uint32)
        f.write(struct.pack('<I', n_particles))

        # Write particle data as float32 array
        for col in columns:
            if col in df.columns:
                data = df[col].values.astype(np.float32)
            else:
                # Default values if column missing
                data = np.zeros(n_particles, dtype=np.float32)
            f.write(data.tobytes())

    return n_particles, os.path.getsize(bin_path)


def convert_dataset(input_dir: str, output_dir: str, subsample: int = 1):
    """Convert all CSV files in a dataset directory."""
    os.makedirs(output_dir, exist_ok=True)

    # Find all CSV files
    csv_files = sorted(glob.glob(os.path.join(input_dir, '*.csv')))

    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return

    print(f"Found {len(csv_files)} snapshots in {input_dir}")

    # Subsample if requested
    if subsample > 1:
        csv_files = csv_files[::subsample]
        print(f"Subsampling: using every {subsample}th snapshot ({len(csv_files)} total)")

    total_csv_size = 0
    total_bin_size = 0

    for i, csv_path in enumerate(csv_files):
        csv_name = os.path.basename(csv_path)
        bin_name = csv_name.replace('.csv', '.bin')
        bin_path = os.path.join(output_dir, bin_name)

        n_particles, bin_size = convert_csv_to_binary(csv_path, bin_path)
        csv_size = os.path.getsize(csv_path)

        total_csv_size += csv_size
        total_bin_size += bin_size

        ratio = csv_size / bin_size if bin_size > 0 else 0
        print(f"  [{i+1}/{len(csv_files)}] {csv_name} → {bin_name}: "
              f"{n_particles} particles, {csv_size/1024/1024:.1f}MB → {bin_size/1024/1024:.1f}MB "
              f"({ratio:.1f}x compression)")

    print(f"\nTotal: {total_csv_size/1024/1024:.1f}MB → {total_bin_size/1024/1024:.1f}MB "
          f"({total_csv_size/total_bin_size:.1f}x compression)")

    return len(csv_files)


def create_manifest(output_base: str, datasets: dict):
    """Create a manifest JSON file listing all binary snapshots."""
    import json

    manifest = {
        "format": "binary",
        "version": 1,
        "columns": ["x", "y", "z", "vx", "vy", "vz", "dens", "temp", "mass", "sound"],
        "datasets": []
    }

    for dataset_id, info in datasets.items():
        bin_dir = os.path.join(output_base, dataset_id)
        if not os.path.exists(bin_dir):
            continue

        bin_files = sorted(glob.glob(os.path.join(bin_dir, '*.bin')))
        snapshots = [os.path.basename(f) for f in bin_files]

        manifest["datasets"].append({
            "id": dataset_id,
            "name": info.get("name", dataset_id),
            "path": dataset_id,
            "snapshots": snapshots,
            "config": info.get("config", {})
        })

    manifest_path = os.path.join(output_base, "manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"Created manifest: {manifest_path}")
    return manifest_path


def main():
    parser = argparse.ArgumentParser(description='Convert SPH snapshots to binary format')
    parser.add_argument('--subsample', type=int, default=1,
                        help='Use every Nth snapshot (default: 1 = all)')
    parser.add_argument('--output', type=str, default='viz/data',
                        help='Output directory for binary files')
    args = parser.parse_args()

    # Dataset configurations
    datasets = {
        "parabolic_rp0.13pc": {
            "name": "Parabolic r_p = 0.13 pc (Close)",
            "input": "parabolic_rp0.13pc/adiabatic/results",
            "config": {
                "r_peri": 0.13,
                "cloud_pos0": [-15.0, 2.8049, 0.0],
                "cloud_vel0": [7.2211, -2.0553, 0.0]
            }
        },
        "parabolic_rp0.4pc": {
            "name": "Parabolic r_p = 0.4 pc (Wide)",
            "input": "parabolic_rp0.4pc/adiabatic/results",
            "config": {
                "r_peri": 0.4,
                "cloud_pos0": [-14.4675, 4.8773, 0.0],
                "cloud_vel0": [7.4071, -1.2149, 0.0]
            }
        }
    }

    print(f"Converting to binary format (subsample={args.subsample})")
    print(f"Output directory: {args.output}\n")

    os.makedirs(args.output, exist_ok=True)

    for dataset_id, info in datasets.items():
        input_dir = info["input"]
        output_dir = os.path.join(args.output, dataset_id)

        if os.path.exists(input_dir):
            print(f"\n=== {info['name']} ===")
            convert_dataset(input_dir, output_dir, args.subsample)
        else:
            print(f"Skipping {dataset_id}: {input_dir} not found")

    # Create manifest
    create_manifest(args.output, datasets)

    # Calculate total size
    total_size = sum(
        os.path.getsize(os.path.join(dp, f))
        for dp, dn, filenames in os.walk(args.output)
        for f in filenames
    )
    print(f"\n=== Total output size: {total_size/1024/1024:.1f} MB ===")


if __name__ == '__main__':
    main()
