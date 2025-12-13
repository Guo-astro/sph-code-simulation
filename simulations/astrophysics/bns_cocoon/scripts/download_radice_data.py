#!/usr/bin/env python3
"""
Download script for Radice et al. (2018) BNS merger ejecta data.

Data source: https://zenodo.org/record/1298474
Paper: https://arxiv.org/abs/1809.11161

This script downloads the ejecta data from numerical relativity simulations
that we use as initial conditions for SR-GSPH shock breakout simulations.
"""

import os
import sys
import argparse
import urllib.request
import hashlib
from pathlib import Path


# Zenodo record for Radice+2018
ZENODO_RECORD = "1298474"
ZENODO_BASE_URL = f"https://zenodo.org/record/{ZENODO_RECORD}/files"

# Available simulation data files
# Format: (filename, description, sha256_hash or None)
AVAILABLE_FILES = {
    # Main data products
    "BLh_equal_mass_outflow.h5": (
        "BLh EOS equal mass (1.35-1.35 Msun) outflow data",
        None
    ),
    "SFHo_equal_mass_outflow.h5": (
        "SFHo EOS equal mass (1.35-1.35 Msun) outflow data", 
        None
    ),
    "DD2_equal_mass_outflow.h5": (
        "DD2 EOS equal mass (1.35-1.35 Msun) outflow data",
        None
    ),
    "LS220_equal_mass_outflow.h5": (
        "LS220 EOS equal mass (1.35-1.35 Msun) outflow data",
        None
    ),
}

# Default files for BNS cocoon simulations
DEFAULT_FILES = ["SFHo_equal_mass_outflow.h5"]


def get_data_dir() -> Path:
    """Get the data directory for downloaded files."""
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / "data" / "radice2018"
    return data_dir


def download_file(url: str, dest: Path, show_progress: bool = True) -> bool:
    """Download a file from URL to destination with progress bar."""
    try:
        if show_progress:
            print(f"Downloading: {url}")
            print(f"Destination: {dest}")
        
        # Create parent directory if needed
        dest.parent.mkdir(parents=True, exist_ok=True)
        
        # Download with progress
        def reporthook(count, block_size, total_size):
            if total_size > 0:
                percent = int(count * block_size * 100 / total_size)
                downloaded = count * block_size / (1024 * 1024)
                total = total_size / (1024 * 1024)
                sys.stdout.write(f"\r  {percent:3d}% ({downloaded:.1f}/{total:.1f} MB)")
                sys.stdout.flush()
        
        if show_progress:
            urllib.request.urlretrieve(url, dest, reporthook)
            print()  # Newline after progress
        else:
            urllib.request.urlretrieve(url, dest)
        
        return True
        
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False


def verify_file(filepath: Path, expected_hash: str | None) -> bool:
    """Verify file integrity using SHA256 hash."""
    if expected_hash is None:
        return True  # Skip verification if no hash provided
    
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256_hash.update(chunk)
    
    actual_hash = sha256_hash.hexdigest()
    if actual_hash != expected_hash:
        print(f"Hash mismatch for {filepath.name}:")
        print(f"  Expected: {expected_hash}")
        print(f"  Got:      {actual_hash}")
        return False
    
    return True


def list_available_files():
    """Print list of available data files."""
    print("Available Radice+2018 data files:")
    print("-" * 60)
    for filename, (description, _) in AVAILABLE_FILES.items():
        default = " [default]" if filename in DEFAULT_FILES else ""
        print(f"  {filename}{default}")
        print(f"    {description}")
    print("-" * 60)
    print(f"\nDownload URL base: {ZENODO_BASE_URL}")


def download_radice_data(
    files: list[str] | None = None,
    data_dir: Path | None = None,
    force: bool = False,
    verify: bool = True
) -> list[Path]:
    """
    Download Radice+2018 ejecta data files.
    
    Parameters
    ----------
    files : list of str, optional
        List of filenames to download. If None, downloads default files.
    data_dir : Path, optional
        Directory to save files. If None, uses default data directory.
    force : bool
        If True, re-download even if file exists.
    verify : bool
        If True, verify file hashes after download.
    
    Returns
    -------
    list of Path
        Paths to downloaded files.
    """
    if files is None:
        files = DEFAULT_FILES
    
    if data_dir is None:
        data_dir = get_data_dir()
    
    data_dir.mkdir(parents=True, exist_ok=True)
    
    downloaded = []
    
    for filename in files:
        if filename not in AVAILABLE_FILES:
            print(f"Warning: Unknown file '{filename}', skipping")
            continue
        
        description, expected_hash = AVAILABLE_FILES[filename]
        dest = data_dir / filename
        
        # Check if already exists
        if dest.exists() and not force:
            print(f"File already exists: {dest}")
            if verify and not verify_file(dest, expected_hash):
                print(f"  Hash verification failed, re-downloading...")
            else:
                downloaded.append(dest)
                continue
        
        # Download
        url = f"{ZENODO_BASE_URL}/{filename}"
        if download_file(url, dest):
            if verify and expected_hash:
                if verify_file(dest, expected_hash):
                    print(f"  ✓ Hash verified")
                else:
                    print(f"  ✗ Hash verification failed")
                    continue
            downloaded.append(dest)
    
    return downloaded


def create_sample_data():
    """Create sample/test data for development without downloading large files."""
    import numpy as np
    
    data_dir = get_data_dir()
    sample_file = data_dir / "sample_ejecta.npz"
    
    print(f"Creating sample data: {sample_file}")
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create synthetic ejecta profile
    n_shells = 100
    v_min = 0.05  # 5% c
    v_max = 0.60  # 60% c
    
    velocity = np.linspace(v_min, v_max, n_shells)
    
    # Power-law mass distribution: dM/dv ~ v^-2
    mass_cumulative = 1e-3 * (v_min / velocity) ** 2  # Total ~1e-3 Msun
    mass_per_shell = np.diff(np.insert(mass_cumulative, 0, 0))
    
    # Angular distribution (polar angle from jet axis)
    n_theta = 20
    theta = np.linspace(0, np.pi, n_theta)
    
    # Ejecta concentrated near equator
    angular_dist = np.sin(theta) ** 2
    angular_dist /= angular_dist.sum()
    
    # Save
    np.savez(
        sample_file,
        velocity=velocity,
        mass_cumulative=mass_cumulative,
        mass_per_shell=mass_per_shell,
        theta=theta,
        angular_distribution=angular_dist,
        total_mass=mass_cumulative[-1],
        description="Synthetic ejecta profile for testing"
    )
    
    print(f"  Created sample data with {n_shells} velocity shells")
    print(f"  Velocity range: {v_min:.2f}c - {v_max:.2f}c")
    print(f"  Total mass: {mass_cumulative[-1]:.2e} Msun")
    
    return sample_file


def main():
    parser = argparse.ArgumentParser(
        description="Download Radice+2018 BNS merger ejecta data"
    )
    
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="List available data files"
    )
    parser.add_argument(
        "--files", "-f", nargs="+",
        help="Specific files to download (default: SFHo equal mass)"
    )
    parser.add_argument(
        "--all", "-a", action="store_true",
        help="Download all available files"
    )
    parser.add_argument(
        "--dir", "-d", type=Path,
        help="Download directory"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Force re-download even if file exists"
    )
    parser.add_argument(
        "--no-verify", action="store_true",
        help="Skip hash verification"
    )
    parser.add_argument(
        "--sample", action="store_true",
        help="Create sample data for testing (no download)"
    )
    
    args = parser.parse_args()
    
    if args.list:
        list_available_files()
        return
    
    if args.sample:
        create_sample_data()
        return
    
    # Determine which files to download
    if args.all:
        files = list(AVAILABLE_FILES.keys())
    elif args.files:
        files = args.files
    else:
        files = DEFAULT_FILES
    
    print("=" * 60)
    print("Radice+2018 BNS Ejecta Data Downloader")
    print("=" * 60)
    print(f"Source: Zenodo record {ZENODO_RECORD}")
    print(f"Paper: arXiv:1809.11161")
    print("=" * 60)
    
    downloaded = download_radice_data(
        files=files,
        data_dir=args.dir,
        force=args.force,
        verify=not args.no_verify
    )
    
    print("\n" + "=" * 60)
    print(f"Downloaded {len(downloaded)} file(s)")
    for p in downloaded:
        print(f"  {p}")
    print("=" * 60)


if __name__ == "__main__":
    main()
