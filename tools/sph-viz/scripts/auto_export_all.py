#!/usr/bin/env python3
"""
Auto-export all SPH simulation results to viz_data format.

This script scans the sphcode repository for simulation results and 
automatically exports them for the visualization tool. It uses tqdm
for progress bars and skips simulations that already have up-to-date
viz_data.

Usage:
    python auto_export_all.py [--force] [--quiet] [--list]

Examples:
    python auto_export_all.py              # Export only new/updated simulations
    python auto_export_all.py --force      # Re-export all simulations
    python auto_export_all.py --list       # List simulations without exporting
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

# Try to import tqdm, fall back to simple progress if not available
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# Get the sphcode root directory (two levels up from this script)
SCRIPT_DIR = Path(__file__).parent.resolve()
SPH_VIZ_DIR = SCRIPT_DIR.parent
SPHCODE_ROOT = SPH_VIZ_DIR.parent.parent

# Directories to scan for simulation results
# Only scan IMBH cloud category results (CAT1, CAT2, CAT3, CAT_OKA)
SCAN_DIRS = [
    'simulations/astrophysics/imbh_cloud/results/CAT1',
    'simulations/astrophysics/imbh_cloud/results/CAT2',
    'simulations/astrophysics/imbh_cloud/results/CAT3',
    'simulations/astrophysics/imbh_cloud/results/CAT_OKA',
]


def find_simulation_dirs(root: Path, quiet: bool = False) -> List[Tuple[Path, bool, int]]:
    """
    Find all directories containing simulation snapshots.
    Returns list of (path, needs_export, snapshot_count) tuples.
    """
    simulations = []
    
    # Use tqdm for scanning progress if available
    scan_iter = SCAN_DIRS
    if TQDM_AVAILABLE and not quiet:
        scan_iter = tqdm(SCAN_DIRS, desc="🔍 Scanning directories", 
                        unit="dir", leave=False, colour="cyan")
    
    for scan_dir in scan_iter:
        base_path = root / scan_dir
        if not base_path.exists():
            continue
        
        # Recursively find directories with snapshot files
        for dirpath in base_path.rglob('*'):
            if not dirpath.is_dir():
                continue
            
            # Check if this directory has snapshots
            snapshots = list(dirpath.glob('snapshot_*.csv'))
            if not snapshots:
                snapshots = list(dirpath.glob('output_*.csv'))
            if not snapshots:
                snapshots = list(dirpath.glob('particles_*.csv'))
            
            if snapshots:
                snapshot_count = len(snapshots)
                
                # Check if viz_data already exists and is up-to-date
                viz_data_dir = dirpath / 'viz_data'
                metadata_file = viz_data_dir / 'metadata.json'
                
                needs_export = True
                if metadata_file.exists():
                    try:
                        import json
                        with open(metadata_file, 'r') as f:
                            metadata = json.load(f)
                        
                        # Check 1: Frame count matches
                        exported_frames = metadata.get('totalFrames', 0)
                        if exported_frames != snapshot_count:
                            needs_export = True  # Frame count mismatch
                        else:
                            # Check 2: No snapshot is newer than metadata
                            metadata_mtime = metadata_file.stat().st_mtime
                            newest_snapshot = max(s.stat().st_mtime for s in snapshots)
                            if newest_snapshot <= metadata_mtime:
                                needs_export = False
                            # else: needs_export stays True
                    except (json.JSONDecodeError, KeyError, IOError):
                        needs_export = True  # Metadata corrupted, re-export
                
                simulations.append((dirpath, needs_export, snapshot_count))
    
    return simulations


def export_simulation(sim_dir: Path, quiet: bool = False) -> Tuple[bool, str]:
    """
    Export a single simulation directory.
    Returns (success, message) tuple.
    """
    export_script = SCRIPT_DIR / 'export_viz_data.py'
    
    cmd = [
        sys.executable,
        str(export_script),
        str(sim_dir),
        '--stride', '1',
        '--quiet',  # Always quiet to avoid interfering with tqdm
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=SPHCODE_ROOT
        )
        if result.returncode == 0:
            return True, "✓"
        else:
            error_msg = result.stderr.strip() if result.stderr else "Unknown error"
            return False, f"✗ {error_msg[:50]}"
    except Exception as e:
        return False, f"✗ {str(e)[:50]}"


def simple_progress(iterable, desc: str, total: int):
    """Simple fallback progress indicator when tqdm is not available."""
    for i, item in enumerate(iterable, 1):
        print(f"\r{desc}: {i}/{total}", end="", flush=True)
        yield item
    print()  # New line after completion


def main():
    parser = argparse.ArgumentParser(
        description='Auto-export all SPH simulation results for visualization'
    )
    parser.add_argument(
        '--force', '-f',
        action='store_true',
        help='Force re-export even if viz_data is up-to-date'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress output (still shows summary)'
    )
    parser.add_argument(
        '--list', '-l',
        action='store_true',
        help='List simulations without exporting'
    )
    
    args = parser.parse_args()
    
    if not TQDM_AVAILABLE and not args.quiet:
        print("💡 Tip: Install tqdm for better progress bars: pip install tqdm\n")
    
    if not args.quiet:
        print("=" * 60)
        print("  🚀 SPH Visualization Data Auto-Export")
        print("=" * 60)
        print(f"\n📂 Scanning: {SPHCODE_ROOT}\n")
    
    # Find all simulations with progress
    simulations = find_simulation_dirs(SPHCODE_ROOT, args.quiet)
    
    if not simulations:
        if not args.quiet:
            print("\n⚠️  No simulation results found.")
            print("   Run some simulations first, e.g.:")
            print("     make -f simulations/benchmarks/sedov/Makefile.sedov sedov_compare_run")
        return 0
    
    # Count stats
    total = len(simulations)
    needs_export_count = sum(1 for _, need, _ in simulations if need or args.force)
    up_to_date_count = total - needs_export_count
    total_snapshots = sum(count for _, _, count in simulations)
    
    if not args.quiet:
        print(f"\n📊 Found {total} simulation(s) with {total_snapshots} total snapshots:")
        print(f"   ├── 🔄 {needs_export_count} need export")
        print(f"   └── ✅ {up_to_date_count} up-to-date (will skip)")
    
    # List mode
    if args.list:
        print("\n📋 Simulations:")
        print("-" * 60)
        for sim_dir, needs, snap_count in simulations:
            rel_path = sim_dir.relative_to(SPHCODE_ROOT)
            if needs:
                status = f"🔄 needs export ({snap_count} snapshots)"
            elif args.force:
                status = f"🔄 will re-export ({snap_count} snapshots)"
            else:
                status = f"✅ up-to-date ({snap_count} snapshots)"
            print(f"   {rel_path}")
            print(f"      {status}")
        return 0
    
    # Check if anything needs export
    if needs_export_count == 0 and not args.force:
        if not args.quiet:
            print("\n✅ All simulations are up-to-date! Nothing to export.")
        return 0
    
    # Filter simulations to export
    to_export = [(sim_dir, snap_count) for sim_dir, needs, snap_count in simulations 
                 if needs or args.force]
    
    if not args.quiet:
        print(f"\n🔄 Exporting {len(to_export)} simulation(s)...\n")
    
    # Export with progress bar
    success_count = 0
    fail_count = 0
    skipped_count = up_to_date_count
    failed_sims = []
    
    if TQDM_AVAILABLE and not args.quiet:
        # Use tqdm with detailed progress
        pbar = tqdm(
            to_export,
            desc="📦 Exporting",
            unit="sim",
            colour="green",
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
        )
        for sim_dir, snap_count in pbar:
            rel_path = sim_dir.relative_to(SPHCODE_ROOT)
            # Update description with current simulation
            short_name = str(rel_path).split('/')[-1][:20]
            pbar.set_postfix_str(f"{short_name} ({snap_count} frames)")
            
            success, msg = export_simulation(sim_dir, args.quiet)
            if success:
                success_count += 1
            else:
                fail_count += 1
                failed_sims.append((rel_path, msg))
        pbar.close()
    else:
        # Fallback progress
        progress_iter = simple_progress(to_export, "Exporting", len(to_export)) if not args.quiet else to_export
        for sim_dir, snap_count in progress_iter:
            success, msg = export_simulation(sim_dir, args.quiet)
            if success:
                success_count += 1
            else:
                fail_count += 1
                failed_sims.append((sim_dir.relative_to(SPHCODE_ROOT), msg))
    
    # Summary
    if not args.quiet:
        print("\n" + "=" * 60)
        print("  📊 Export Summary")
        print("=" * 60)
        print(f"   ✅ Exported:  {success_count}")
        print(f"   ⏭️  Skipped:   {skipped_count} (already up-to-date)")
        if fail_count > 0:
            print(f"   ❌ Failed:    {fail_count}")
            print("\n   Failed simulations:")
            for rel_path, msg in failed_sims:
                print(f"      • {rel_path}: {msg}")
        print("=" * 60)
    
    return 0 if fail_count == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
