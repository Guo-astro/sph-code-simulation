#!/usr/bin/env python3
"""
SPH Visualization Data Exporter

This script exports SPH simulation output to a format optimized for the
SPH-Viz web visualization tool. It converts CSV particle snapshots to
binary Float32 arrays with accompanying metadata.

Usage:
    python export_viz_data.py <simulation_dir> [options]

Example:
    python export_viz_data.py simulations/benchmarks/sedov/results/gsph_wendland --stride 1
    python export_viz_data.py lane_emden/results/n3_gsph --max-frames 100

The output will be written to <simulation_dir>/viz_data/ with:
- metadata.json: Simulation metadata and field offsets
- frame_XXXXX.bin: Binary Float32 particle data for each frame
- frame_XXXXX.json: Per-frame metadata (time, particle count)
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

# Default field layout (stride = 11 floats per particle)
DEFAULT_FIELDS = [
    'x', 'y', 'z',           # Position
    'vx', 'vy', 'vz',        # Velocity
    'mass',                   # Mass
    'density',               # Density
    'pressure',              # Pressure
    'energy',                # Internal energy
    'smoothing_length',      # Smoothing length
]

# Mapping from various CSV column names to standard names
COLUMN_ALIASES = {
    # Position
    'pos_x': 'x', 'pos_y': 'y', 'pos_z': 'z',
    # Velocity
    'vel_x': 'vx', 'vel_y': 'vy', 'vel_z': 'vz',
    'v_x': 'vx', 'v_y': 'vy', 'v_z': 'vz',
    # Mass
    'm': 'mass',
    # Density
    'rho': 'density', 'dens': 'density',
    # Pressure
    'p': 'pressure', 'P': 'pressure', 'pres': 'pressure',
    # Energy
    'u': 'energy', 'e': 'energy', 'internal_energy': 'energy', 'ene': 'energy',
    # Smoothing length
    'h': 'smoothing_length', 'sml': 'smoothing_length',
    # Sound speed
    'cs': 'sound_speed', 'sound': 'sound_speed',
}


def find_snapshots(sim_dir: Path) -> List[Path]:
    """Find all snapshot files in a simulation directory."""
    snapshots = []
    
    # Look for different snapshot patterns
    patterns = [
        'snapshot_*.csv',
        'output_*.csv',
        'particles_*.csv',
        '*.csv',
    ]
    
    for pattern in patterns:
        found = sorted(sim_dir.glob(pattern))
        if found:
            # Filter out non-snapshot files
            found = [f for f in found if not f.name.startswith('config')]
            if found:
                snapshots = found
                break
    
    return snapshots


def parse_csv_header(filepath: Path) -> Tuple[List[str], str, int]:
    """Parse CSV header and detect delimiter. Returns (columns, delimiter, skip_rows)."""
    skip_rows = 0
    header_line = None
    
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            # Skip comment lines
            if line.startswith('#'):
                skip_rows = i + 1
                continue
            # First non-comment line is the header
            if header_line is None:
                header_line = line
                skip_rows = i + 1
                break
    
    if header_line is None:
        return [], ',', 0
    
    # Detect delimiter
    if '\t' in header_line:
        delimiter = '\t'
    elif ',' in header_line:
        delimiter = ','
    else:
        delimiter = ' '
    
    # Parse header
    columns = [c.strip() for c in header_line.split(delimiter)]
    
    # Normalize column names
    normalized = []
    for col in columns:
        col_lower = col.lower().strip()
        if col_lower in COLUMN_ALIASES:
            normalized.append(COLUMN_ALIASES[col_lower])
        else:
            normalized.append(col_lower)
    
    return normalized, delimiter, skip_rows


def load_snapshot(filepath: Path) -> Tuple[np.ndarray, List[str], float]:
    """Load a snapshot CSV file and return data, columns, and time."""
    columns, delimiter, skip_rows = parse_csv_header(filepath)
    
    # Load data
    try:
        data = np.loadtxt(filepath, delimiter=delimiter, skiprows=skip_rows)
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None, [], 0.0
    
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    # Try to extract time from file header comments
    time = 0.0
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                break
            # Look for time in header
            if 'Time (code):' in line:
                try:
                    time = float(line.split(':')[-1].strip())
                except ValueError:
                    pass
            elif 'Time:' in line and 'code' not in line.lower():
                try:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        time_str = parts[-1].strip().split()[0]
                        time = float(time_str)
                except (ValueError, IndexError):
                    pass
    
    # Fallback: try to extract time from filename
    if time == 0.0:
        stem = filepath.stem
        for prefix in ['snapshot_', 'output_', 'particles_']:
            if stem.startswith(prefix):
                try:
                    time = float(stem[len(prefix):])
                except ValueError:
                    pass
                break
    
    return data, columns, time


def create_binary_frame(
    data: np.ndarray,
    columns: List[str],
    output_fields: List[str]
) -> bytes:
    """Convert snapshot data to binary Float32 format."""
    n_particles = data.shape[0]
    n_fields = len(output_fields)
    
    # Create output array
    output = np.zeros((n_particles, n_fields), dtype=np.float32)
    
    # Map columns to output fields
    for i, field in enumerate(output_fields):
        if field in columns:
            col_idx = columns.index(field)
            output[:, i] = data[:, col_idx].astype(np.float32)
        else:
            # Fill with zeros or defaults
            if field in ['mass']:
                output[:, i] = 1.0  # Default mass
            else:
                output[:, i] = 0.0
    
    return output.tobytes()


def compute_bounding_box(data: np.ndarray, columns: List[str]) -> Dict:
    """Compute bounding box from particle positions."""
    try:
        x_idx = columns.index('x')
        y_idx = columns.index('y')
        z_idx = columns.index('z') if 'z' in columns else None

        min_x, max_x = float(data[:, x_idx].min()), float(data[:, x_idx].max())
        min_y, max_y = float(data[:, y_idx].min()), float(data[:, y_idx].max())

        if z_idx is not None:
            min_z, max_z = float(data[:, z_idx].min()), float(data[:, z_idx].max())
        else:
            min_z, max_z = 0.0, 0.0

        return {
            'min': [min_x, min_y, min_z],
            'max': [max_x, max_y, max_z]
        }
    except (ValueError, IndexError):
        return {
            'min': [-1.0, -1.0, -1.0],
            'max': [1.0, 1.0, 1.0]
        }


def compute_pv_ranges(data: np.ndarray, columns: List[str], v_lsr: float = 0.0) -> Dict:
    """
    Compute position and velocity ranges for P-V diagram.

    Args:
        data: Particle data array
        columns: Column names
        v_lsr: LSR velocity offset (km/s)

    Returns:
        Dictionary with positionRange and velocityRange
    """
    try:
        x_idx = columns.index('x')
        y_idx = columns.index('y')
        vx_idx = columns.index('vx')
        vy_idx = columns.index('vy')
        vz_idx = columns.index('vz') if 'vz' in columns else None

        # Position range (use XY plane extent)
        x = data[:, x_idx]
        y = data[:, y_idx]
        pos_extent = max(np.abs(x).max(), np.abs(y).max())

        # Velocity range (use all velocity components + V_LSR offset)
        vx = data[:, vx_idx]
        vy = data[:, vy_idx]
        if vz_idx is not None:
            vz = data[:, vz_idx]
        else:
            vz = np.zeros_like(vx)

        # Line-of-sight velocity (z-component + V_LSR)
        v_los_min = float(vz.min()) + v_lsr
        v_los_max = float(vz.max()) + v_lsr

        # Add 15% margin
        pos_margin = pos_extent * 0.15
        vel_margin = (v_los_max - v_los_min) * 0.15 or 10.0

        return {
            'positionRange': [float(-pos_extent - pos_margin), float(pos_extent + pos_margin)],
            'velocityRange': [float(v_los_min - vel_margin), float(v_los_max + vel_margin)],
            'velocityStats': {
                'vx': [float(vx.min()), float(vx.max())],
                'vy': [float(vy.min()), float(vy.max())],
                'vz': [float(vz.min()), float(vz.max())],
            }
        }
    except (ValueError, IndexError):
        return {
            'positionRange': [-5.0, 5.0],
            'velocityRange': [-150.0, 50.0],
            'velocityStats': None
        }


def detect_sph_method(sim_dir: Path) -> str:
    """Try to detect SPH method from directory name or config."""
    dir_name = sim_dir.name.lower()
    
    if 'gsph' in dir_name:
        return 'GSPH'
    elif 'ssph' in dir_name:
        return 'SSPH'
    elif 'gdisph' in dir_name:
        return 'GDISPH'
    elif 'disph' in dir_name:
        return 'DISPH'
    elif 'srgsph' in dir_name:
        return 'SRGSPH'
    
    # Check config file
    config_path = sim_dir / 'config.json'
    if config_path.exists():
        try:
            with open(config_path) as f:
                config = json.load(f)
            if 'SPH' in config:
                return config['SPH'].get('method', 'Unknown')
        except Exception:
            pass
    
    return 'Unknown'


def detect_kernel(sim_dir: Path) -> str:
    """Try to detect kernel from directory name or config."""
    dir_name = sim_dir.name.lower()
    
    if 'wendland' in dir_name:
        return 'WendlandC2'
    elif 'cubic' in dir_name:
        return 'CubicSpline'
    
    return 'Unknown'


def export_simulation(
    sim_dir: Path,
    output_dir: Optional[Path] = None,
    max_frames: Optional[int] = None,
    stride: int = 1,
    verbose: bool = True
):
    """Export simulation data to visualization format."""
    if output_dir is None:
        output_dir = sim_dir / 'viz_data'
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find snapshots
    snapshots = find_snapshots(sim_dir)
    if not snapshots:
        print(f"No snapshots found in {sim_dir}")
        return False
    
    if verbose:
        print(f"Found {len(snapshots)} snapshots")
    
    # Apply stride and max frames
    snapshots = snapshots[::stride]
    if max_frames:
        snapshots = snapshots[:max_frames]
    
    if verbose:
        print(f"Exporting {len(snapshots)} frames (stride={stride})")
    
    # Define output fields
    output_fields = DEFAULT_FIELDS
    field_offsets = {field: i for i, field in enumerate(output_fields)}
    
    # Process first frame to get metadata
    first_data, first_columns, first_time = load_snapshot(snapshots[0])
    if first_data is None:
        print("Failed to load first snapshot")
        return False

    n_particles = first_data.shape[0]
    bounding_box = compute_bounding_box(first_data, first_columns)
    pv_ranges = compute_pv_ranges(first_data, first_columns, v_lsr=-120.0)  # Default V_LSR
    
    # Detect simulation properties
    sph_method = detect_sph_method(sim_dir)
    kernel = detect_kernel(sim_dir)
    
    # Process all frames
    times = []
    for i, snapshot_path in enumerate(snapshots):
        if verbose and (i % 10 == 0 or i == len(snapshots) - 1):
            print(f"  Processing frame {i + 1}/{len(snapshots)}: {snapshot_path.name}")
        
        data, columns, time = load_snapshot(snapshot_path)
        if data is None:
            print(f"  Skipping {snapshot_path.name} (load error)")
            continue
        
        # Update bounding box
        frame_bbox = compute_bounding_box(data, columns)
        bounding_box['min'] = [
            min(bounding_box['min'][j], frame_bbox['min'][j])
            for j in range(3)
        ]
        bounding_box['max'] = [
            max(bounding_box['max'][j], frame_bbox['max'][j])
            for j in range(3)
        ]

        # Update P-V ranges (global min/max across all frames)
        frame_pv = compute_pv_ranges(data, columns, v_lsr=-120.0)
        pv_ranges['positionRange'][0] = min(pv_ranges['positionRange'][0], frame_pv['positionRange'][0])
        pv_ranges['positionRange'][1] = max(pv_ranges['positionRange'][1], frame_pv['positionRange'][1])
        pv_ranges['velocityRange'][0] = min(pv_ranges['velocityRange'][0], frame_pv['velocityRange'][0])
        pv_ranges['velocityRange'][1] = max(pv_ranges['velocityRange'][1], frame_pv['velocityRange'][1])
        
        # Export binary frame
        binary_data = create_binary_frame(data, columns, output_fields)
        frame_path = output_dir / f'frame_{i:05d}.bin'
        with open(frame_path, 'wb') as f:
            f.write(binary_data)
        
        # Export frame metadata
        frame_meta = {
            'frameIndex': i,
            'time': time,
            'particleCount': data.shape[0],
            'sourceFile': snapshot_path.name
        }
        frame_meta_path = output_dir / f'frame_{i:05d}.json'
        with open(frame_meta_path, 'w') as f:
            json.dump(frame_meta, f, indent=2)
        
        times.append(time)
    
    # Export main metadata
    metadata = {
        'name': sim_dir.name,
        'description': f'SPH simulation: {sim_dir.name}',
        'method': sph_method,
        'kernel': kernel,
        'dimensions': 3 if 'z' in first_columns else 2,
        'totalFrames': len(snapshots),
        'particleCount': n_particles,
        'timeRange': [times[0] if times else 0, times[-1] if times else 0],
        'boundingBox': bounding_box,
        'pvDiagram': pv_ranges,  # Pre-computed P-V diagram ranges
        'stride': len(output_fields),
        'fieldOffsets': field_offsets,
        'createdAt': str(np.datetime64('now')),
        'sourceDir': str(sim_dir.resolve())
    }
    
    metadata_path = output_dir / 'metadata.json'
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    if verbose:
        print("\n✓ Export complete!")
        print(f"  Output: {output_dir}")
        print(f"  Frames: {len(snapshots)}")
        print(f"  Particles: {n_particles}")
        print(f"  Method: {sph_method}")
        print(f"  Kernel: {kernel}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Export SPH simulation data for web visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        'simulation_dir',
        type=Path,
        help='Path to simulation results directory'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        default=None,
        help='Output directory (default: <sim_dir>/viz_data)'
    )
    parser.add_argument(
        '-m', '--max-frames',
        type=int,
        default=None,
        help='Maximum number of frames to export'
    )
    parser.add_argument(
        '-s', '--stride',
        type=int,
        default=1,
        help='Frame stride (export every Nth frame)'
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress progress output'
    )
    
    args = parser.parse_args()
    
    if not args.simulation_dir.exists():
        print(f"Error: Directory not found: {args.simulation_dir}")
        sys.exit(1)
    
    success = export_simulation(
        args.simulation_dir,
        output_dir=args.output,
        max_frames=args.max_frames,
        stride=args.stride,
        verbose=not args.quiet
    )
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
