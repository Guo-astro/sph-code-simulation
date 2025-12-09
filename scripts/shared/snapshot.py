#!/usr/bin/env python3
"""
Snapshot reading module - Single Source of Truth (SSOT)

This module provides unified snapshot reading functionality for SPH simulation
CSV output files. DO NOT re-implement CSV parsing - use this module instead.

Features:
- Automatic header/comment handling
- Metadata extraction from headers
- Column validation
- Type-safe data access

Author: SSOT refactoring for sphcode project
"""

import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Union
from dataclasses import dataclass


@dataclass
class SnapshotMetadata:
    """Container for snapshot metadata extracted from header."""
    time: Optional[float] = None
    step: Optional[int] = None
    n_particles: Optional[int] = None
    dim: Optional[int] = None
    # Lane-Emden IC parameters (if present)
    rho_c: Optional[float] = None
    R: Optional[float] = None
    M: Optional[float] = None
    K: Optional[float] = None
    alpha: Optional[float] = None
    gamma: Optional[float] = None
    # Additional parameters
    extra: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.extra is None:
            self.extra = {}


def parse_header(filepath: Union[str, Path]) -> SnapshotMetadata:
    """
    Parse snapshot header to extract metadata.
    
    Parameters
    ----------
    filepath : str or Path
        Path to snapshot CSV file
        
    Returns
    -------
    SnapshotMetadata
        Extracted metadata
    """
    meta = SnapshotMetadata()
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                line = line[1:].strip()  # Remove # prefix
                
                # Parse key: value pairs
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip().lower()
                    value = value.strip()
                    
                    try:
                        if 'time' in key and 'code' in key:
                            meta.time = float(value)
                        elif key == 'time':
                            meta.time = float(value)
                        elif 'step' in key:
                            meta.step = int(value)
                        elif 'particles' in key or 'n_particles' in key:
                            meta.n_particles = int(value)
                        elif key == 'dim' or 'dimension' in key:
                            meta.dim = int(value)
                        elif 'rho_c' in key or 'central_density' in key:
                            meta.rho_c = float(value)
                        elif key == 'r' or 'radius' in key:
                            meta.R = float(value)
                        elif key == 'm' or 'mass' in key:
                            meta.M = float(value)
                        elif key == 'k' or 'polytropic_k' in key:
                            meta.K = float(value)
                        elif key == 'alpha':
                            meta.alpha = float(value)
                        elif key == 'gamma':
                            meta.gamma = float(value)
                        else:
                            meta.extra[key] = value
                    except ValueError:
                        meta.extra[key] = value
    except Exception:
        pass
    
    return meta


def read_snapshot(
    filepath: Union[str, Path],
    columns: Optional[List[str]] = None,
    include_metadata: bool = False
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, SnapshotMetadata]]:
    """
    Read a snapshot CSV file with automatic header handling.
    
    This is the RECOMMENDED way to read snapshot files.
    
    Parameters
    ----------
    filepath : str or Path
        Path to snapshot CSV file
    columns : list of str, optional
        Specific columns to read (None = all)
    include_metadata : bool
        If True, also return parsed metadata
        
    Returns
    -------
    df : DataFrame
        Particle data
    metadata : SnapshotMetadata (if include_metadata=True)
        Parsed header metadata
        
    Example
    -------
    >>> df = read_snapshot("snapshot_0010.csv")
    >>> r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    
    >>> df, meta = read_snapshot("snapshot_0010.csv", include_metadata=True)
    >>> print(f"Time: {meta.time}")
    """
    filepath = Path(filepath)
    
    # Parse metadata first
    metadata = parse_header(filepath) if include_metadata else None
    
    # Read CSV, skipping comment lines
    try:
        df = pd.read_csv(filepath, comment='#')
    except Exception as e:
        raise IOError(f"Could not read snapshot {filepath}: {e}")
    
    # Select specific columns if requested
    if columns is not None:
        available = [c for c in columns if c in df.columns]
        if available:
            df = df[available]
    
    if include_metadata:
        return df, metadata
    return df


def get_time_from_snapshot(filepath: Union[str, Path]) -> Optional[float]:
    """
    Quick extraction of time from snapshot header.
    
    Parameters
    ----------
    filepath : str or Path
        Path to snapshot file
        
    Returns
    -------
    time : float or None
        Simulation time, or None if not found
    """
    return parse_header(filepath).time


def find_snapshots(
    directory: Union[str, Path],
    pattern: str = "snapshot_*.csv"
) -> List[Path]:
    """
    Find all snapshot files in a directory.
    
    Parameters
    ----------
    directory : str or Path
        Directory to search
    pattern : str
        Glob pattern for snapshot files
        
    Returns
    -------
    list of Path
        Sorted list of snapshot file paths
    """
    directory = Path(directory)
    files = sorted(directory.glob(pattern))
    return files


def load_snapshot_series(
    directory: Union[str, Path],
    pattern: str = "snapshot_*.csv",
    columns: Optional[List[str]] = None,
    max_snapshots: Optional[int] = None
) -> List[Tuple[Path, pd.DataFrame, SnapshotMetadata]]:
    """
    Load a series of snapshots from a directory.
    
    Parameters
    ----------
    directory : str or Path
        Directory containing snapshots
    pattern : str
        Glob pattern for snapshot files
    columns : list of str, optional
        Specific columns to load
    max_snapshots : int, optional
        Maximum number of snapshots to load
        
    Returns
    -------
    list of (filepath, dataframe, metadata) tuples
    """
    files = find_snapshots(directory, pattern)
    
    if max_snapshots is not None:
        files = files[:max_snapshots]
    
    results = []
    for f in files:
        df, meta = read_snapshot(f, columns=columns, include_metadata=True)
        results.append((f, df, meta))
    
    return results


def compute_radii(df: pd.DataFrame) -> np.ndarray:
    """
    Compute radial distance from origin for all particles.
    
    Parameters
    ----------
    df : DataFrame
        Snapshot data with pos_x, pos_y, [pos_z] columns
        
    Returns
    -------
    r : ndarray
        Radial distances
    """
    r_sq = df['pos_x']**2 + df['pos_y']**2
    if 'pos_z' in df.columns:
        r_sq += df['pos_z']**2
    return np.sqrt(r_sq)


def bin_radial_profile(
    r: np.ndarray,
    values: np.ndarray,
    r_max: float = None,
    n_bins: int = 50
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute binned radial profile of a quantity.
    
    Parameters
    ----------
    r : ndarray
        Radial distances
    values : ndarray
        Values to bin (e.g., density)
    r_max : float, optional
        Maximum radius for bins (default: max(r))
    n_bins : int
        Number of bins
        
    Returns
    -------
    r_centers : ndarray
        Bin centers
    binned_values : ndarray
        Mean value in each bin (NaN for empty bins)
    """
    if r_max is None:
        r_max = r.max()
    
    r_bins = np.linspace(0, r_max, n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    binned = []
    for i in range(n_bins):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if mask.sum() > 0:
            binned.append(np.mean(values[mask]))
        else:
            binned.append(np.nan)
    
    return r_centers, np.array(binned)


if __name__ == "__main__":
    print("=" * 60)
    print("Snapshot Module Test")
    print("=" * 60)
    
    # Test with a sample file if available
    test_paths = [
        "sample/imbh_cloud/results/ic_relax_10k/snapshot_0032.csv",
        "sample/shock_tube/results/snapshot_0000.csv",
    ]
    
    for test_path in test_paths:
        if os.path.exists(test_path):
            print(f"\n✓ Testing with: {test_path}")
            
            # Test metadata parsing
            meta = parse_header(test_path)
            print(f"  - Time: {meta.time}")
            print(f"  - N particles: {meta.n_particles}")
            if meta.rho_c:
                print(f"  - ρ_c: {meta.rho_c}")
            
            # Test full read
            df, meta2 = read_snapshot(test_path, include_metadata=True)
            print(f"  - Columns: {list(df.columns)[:5]}...")
            print(f"  - Shape: {df.shape}")
            
            # Test radii computation
            if 'pos_x' in df.columns:
                r = compute_radii(df)
                print(f"  - r range: [{r.min():.4f}, {r.max():.4f}]")
            
            break
    else:
        print("\n⚠️  No test files found")
    
    print("\n" + "=" * 60)
    print("Module loaded successfully!")
    print("=" * 60)
