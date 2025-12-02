#!/usr/bin/env python3
"""
IMBH Tidal Disruption Animation - Concatenated Multi-Segment Visualization
==========================================================================

Creates animations by concatenating snapshots from multiple run segments
(e.g., original run + continuation runs) into a single unified animation.

Key Features:
-------------
- Concatenate multiple snapshot directories in chronological order
- Automatically re-index frames across segments
- Support for both 'single' (3x3) and 'tidal' (3x4) visualization modes
- Time continuity verification and reporting

Usage:
------
    # Concatenate original run with continuation
    python animate_concatenated.py \
        --input-dirs results/original results/continuation/t2to4 \
        --output-dir results/combined_animation \
        --bh-mass 100000

    # List mode (comma-separated)
    python animate_concatenated.py \
        --input-dirs "results/run1,results/run2,results/run3" \
        --output-dir results/combined

Examples:
---------
    # Combine nocool original + continuation
    python animate_concatenated.py \
        --input-dirs sample/imbh_cloud/results/imbh_relaxed_10k_b3pc_nocool \
                     sample/imbh_cloud/results/continuation/10k_b3pc_nocool_t2to4 \
        --output-dir sample/imbh_cloud/results/combined/nocool_full \
        --bh-mass 100000 --mode tidal

Units: IMBH_ENCOUNTER ([L]=pc, [M]=1000 M_☉, [V]=km/s, [T]=0.978 Myr)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm, Normalize
from matplotlib.gridspec import GridSpec
import argparse
from pathlib import Path
import glob
import re
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# PHYSICAL CONSTANTS AND UNIT CONVERSIONS
# =============================================================================

TIME_UNIT = 0.978  # Myr per code unit (1 pc / 1 km/s = 0.978 Myr)
VELOCITY_UNIT = 1.0  # km/s per code unit
LENGTH_UNIT = 1.0  # pc per code unit
MASS_UNIT = 1000.0  # M_sun per code unit mass

GAMMA = 5.0 / 3.0  # Adiabatic index
MU = 2.3  # Mean molecular weight (molecular cloud)
M_H = 1.67e-24  # Hydrogen mass [g]
K_B = 1.38e-16  # Boltzmann constant [erg/K]
VELOCITY_TO_CGS = 1e5  # km/s to cm/s
TEMP_CONVERSION = MU * M_H / K_B * VELOCITY_TO_CGS**2


# =============================================================================
# DATA LOADING UTILITIES
# =============================================================================

def load_snapshot(filepath):
    """Load a single CSV snapshot with robust header detection."""
    try:
        with open(filepath, 'r') as f:
            skip_lines = 0
            for line in f:
                if line.startswith('#'):
                    skip_lines += 1
                else:
                    break
        
        data = np.genfromtxt(filepath, delimiter=',', skip_header=skip_lines, names=True)
        
        snapshot = {
            'pos_x': data['pos_x'],
            'pos_y': data['pos_y'],
            'pos_z': data['pos_z'],
            'vel_x': data['vel_x'],
            'vel_y': data['vel_y'],
            'vel_z': data['vel_z'],
            'dens': data['dens'],
            'pres': data['pres'],
            'mass': data['mass'],
            'sml': data['sml'],
            'sound': data['sound'],
            'filepath': filepath,
        }
        
        if 'ene' in data.dtype.names:
            snapshot['ene'] = data['ene']
        
        return snapshot
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def extract_time_from_snapshot(filepath):
    """Extract time from snapshot header comment."""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    # Look for "Time (code):" in header
                    match = re.search(r'Time\s*\(code\)\s*:\s*([\d.eE+-]+)', line, re.IGNORECASE)
                    if match:
                        return float(match.group(1))
                    # Also try generic "time=" format
                    match = re.search(r'time\s*=\s*([\d.eE+-]+)', line, re.IGNORECASE)
                    if match:
                        return float(match.group(1))
                else:
                    break
    except Exception:
        pass
    return None


def load_energy_file(filepath):
    """Load energy conservation data from energy.dat."""
    try:
        data = np.loadtxt(filepath, skiprows=1)
        return {
            'time': data[:, 0],
            'E_kin': data[:, 1],
            'E_int': data[:, 2],
            'E_pot': data[:, 3],
            'E_tot': data[:, 4],
        }
    except Exception as e:
        print(f"Could not load energy file: {e}")
        return None


def discover_snapshots(input_dirs):
    """
    Discover and concatenate snapshots from multiple directories.
    
    Returns:
        list of (filepath, estimated_time, segment_index)
    """
    all_snapshots = []
    
    for seg_idx, input_dir in enumerate(input_dirs):
        input_path = Path(input_dir)
        if not input_path.exists():
            print(f"⚠ Warning: Directory not found: {input_dir}")
            continue
            
        snapshot_files = sorted(glob.glob(str(input_path / "snapshot_*.csv")))
        
        if len(snapshot_files) == 0:
            print(f"⚠ Warning: No snapshots found in {input_dir}")
            continue
        
        print(f"📁 Segment {seg_idx + 1}: {input_dir}")
        print(f"   Found {len(snapshot_files)} snapshots")
        
        # Extract times from first and last snapshots
        t_first = extract_time_from_snapshot(snapshot_files[0])
        t_last = extract_time_from_snapshot(snapshot_files[-1])
        
        if t_first is not None and t_last is not None:
            print(f"   Time range: {t_first:.3f} → {t_last:.3f} (code units)")
            print(f"              ({t_first * TIME_UNIT:.2f} → {t_last * TIME_UNIT:.2f} Myr)")
        
        for filepath in snapshot_files:
            t = extract_time_from_snapshot(filepath)
            all_snapshots.append({
                'filepath': filepath,
                'time': t if t is not None else float('inf'),
                'segment': seg_idx,
                'segment_dir': input_dir,
            })
    
    # Sort by time
    all_snapshots.sort(key=lambda x: x['time'])
    
    return all_snapshots


def merge_snapshots_by_time(snapshots, overlap_tolerance=0.001):
    """
    Merge snapshots from multiple segments, handling overlaps.
    
    When continuation starts from a snapshot, there may be duplicate times.
    This function removes duplicates, preferring later segments.
    
    Args:
        snapshots: list of snapshot info dicts
        overlap_tolerance: time difference below which snapshots are considered duplicates
    
    Returns:
        deduplicated list of snapshots
    """
    if not snapshots:
        return []
    
    merged = []
    prev_time = -float('inf')
    
    for snap in snapshots:
        t = snap['time']
        
        # Check for duplicate/overlap
        if abs(t - prev_time) < overlap_tolerance:
            # Replace previous with this one (prefer later segment)
            if snap['segment'] > merged[-1]['segment']:
                merged[-1] = snap
        else:
            merged.append(snap)
            prev_time = t
    
    return merged


# =============================================================================
# PHYSICS CALCULATIONS
# =============================================================================

def compute_temperature(snapshot):
    """Compute temperature from internal energy or pressure."""
    if 'ene' in snapshot:
        u = snapshot['ene']
    else:
        u = snapshot['pres'] / (snapshot['dens'] * (GAMMA - 1))
    
    T = u * TEMP_CONVERSION * (GAMMA - 1)
    return np.clip(T, 10, 1e8)


def compute_mach_number(snapshot):
    """Compute local Mach number."""
    v = np.sqrt(snapshot['vel_x']**2 + snapshot['vel_y']**2 + snapshot['vel_z']**2)
    cs = snapshot['sound']
    return v / np.maximum(cs, 1e-10)


def compute_center_of_mass(snapshot):
    """Compute center of mass position."""
    m = snapshot['mass']
    M_tot = np.sum(m)
    com = np.array([
        np.sum(m * snapshot['pos_x']) / M_tot,
        np.sum(m * snapshot['pos_y']) / M_tot,
        np.sum(m * snapshot['pos_z']) / M_tot
    ])
    return com


# =============================================================================
# MAIN CONCATENATED ANIMATION
# =============================================================================

def create_concatenated_animation(input_dirs, output_dir, bh_mass, fps=10, 
                                   mode='single', dpi=100):
    """Create animation from concatenated multi-segment snapshots."""
    
    print("=" * 70)
    print("  Discovering and concatenating snapshots from multiple segments")
    print("=" * 70)
    
    # Discover all snapshots
    all_snapshots = discover_snapshots(input_dirs)
    
    if len(all_snapshots) == 0:
        print("❌ No snapshots found in any input directory!")
        return False
    
    # Merge and deduplicate
    merged_snapshots = merge_snapshots_by_time(all_snapshots)
    
    print()
    print("=" * 70)
    print(f"  Total frames after merge: {len(merged_snapshots)}")
    print("=" * 70)
    
    # Report segment breakdown
    segment_counts = {}
    for snap in merged_snapshots:
        seg = snap['segment']
        segment_counts[seg] = segment_counts.get(seg, 0) + 1
    
    for seg, count in sorted(segment_counts.items()):
        print(f"   Segment {seg + 1}: {count} frames")
    
    # Time range
    t_start = merged_snapshots[0]['time']
    t_end = merged_snapshots[-1]['time']
    print(f"   Time range: {t_start:.3f} → {t_end:.3f} (code units)")
    print(f"              ({t_start * TIME_UNIT:.2f} → {t_end * TIME_UNIT:.2f} Myr)")
    print()
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Extract file paths for animation
    snapshot_files = [snap['filepath'] for snap in merged_snapshots]
    snapshot_times = [snap['time'] for snap in merged_snapshots]
    snapshot_segments = [snap['segment'] for snap in merged_snapshots]
    
    # Scan for plot extents
    print("Scanning snapshots for plot extents...")
    all_x, all_y, all_z = [], [], []
    sample_indices = list(range(0, len(snapshot_files), max(1, len(snapshot_files) // 10)))
    
    for i in sample_indices:
        snap = load_snapshot(snapshot_files[i])
        if snap:
            all_x.extend(snap['pos_x'])
            all_y.extend(snap['pos_y'])
            all_z.extend(snap['pos_z'])
    
    margin = 1.2
    x_range = max(abs(min(all_x)), abs(max(all_x))) * margin
    y_range = max(abs(min(all_y)), abs(max(all_y))) * margin
    z_range = max(abs(min(all_z)), abs(max(all_z))) * margin
    
    xlim = [-x_range, x_range]
    ylim = [-y_range, y_range]
    zlim = [-z_range, z_range]
    
    print(f"   X: [{xlim[0]:.1f}, {xlim[1]:.1f}] pc")
    print(f"   Y: [{ylim[0]:.1f}, {ylim[1]:.1f}] pc")
    print(f"   Z: [{zlim[0]:.1f}, {zlim[1]:.1f}] pc")
    
    # Load first snapshot for normalization
    snap0 = load_snapshot(snapshot_files[0])
    rho_0 = np.median(snap0['dens'])
    P_0 = np.median(snap0['pres'])
    
    # =========================================================================
    # CREATE ANIMATION BASED ON MODE
    # =========================================================================
    
    if mode == 'tidal':
        fig = plt.figure(figsize=(22, 16))
        gs = GridSpec(3, 4, figure=fig, hspace=0.3, wspace=0.25)
        
        ax_xy_dens = fig.add_subplot(gs[0, 0])
        ax_xz_dens = fig.add_subplot(gs[0, 1])
        ax_yz_dens = fig.add_subplot(gs[0, 2])
        ax_3d = fig.add_subplot(gs[0, 3], projection='3d')
        
        ax_xy_temp = fig.add_subplot(gs[1, 0])
        ax_xz_temp = fig.add_subplot(gs[1, 1])
        ax_yz_temp = fig.add_subplot(gs[1, 2])
        ax_mach_hist = fig.add_subplot(gs[1, 3])
        
        ax_xy_shock = fig.add_subplot(gs[2, 0])
        ax_xz_shock = fig.add_subplot(gs[2, 1])
        ax_yz_shock = fig.add_subplot(gs[2, 2])
        ax_info = fig.add_subplot(gs[2, 3])
        
    else:  # single mode
        fig = plt.figure(figsize=(18, 16))
        gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.25)
        
        ax_3d = fig.add_subplot(gs[0, 0], projection='3d')
        ax_xy = fig.add_subplot(gs[0, 1])
        ax_xz = fig.add_subplot(gs[0, 2])
        
        ax_dens_hist = fig.add_subplot(gs[1, 0])
        ax_temp_hist = fig.add_subplot(gs[1, 1])
        ax_mach_hist = fig.add_subplot(gs[1, 2])
        
        ax_energy = fig.add_subplot(gs[2, 0])
        ax_com = fig.add_subplot(gs[2, 1])
        ax_info = fig.add_subplot(gs[2, 2])
    
    # COM tracking
    com_history = []
    
    def update(frame):
        """Update animation frame."""
        if frame >= len(snapshot_files):
            return
        
        snap = load_snapshot(snapshot_files[frame])
        if snap is None:
            return
        
        time = snapshot_times[frame]
        segment = snapshot_segments[frame]
        time_myr = time * TIME_UNIT
        
        x, y, z = snap['pos_x'], snap['pos_y'], snap['pos_z']
        rho = snap['dens']
        T = compute_temperature(snap)
        mach = compute_mach_number(snap)
        
        com = compute_center_of_mass(snap)
        com_history.append(com)
        
        # Shock detection
        is_shocked = (rho / rho_0 > 4) | (mach > 1.5)
        
        # Clear all axes
        for ax in fig.axes:
            ax.clear()
        
        if mode == 'tidal':
            # === Row 0: Cross-section DENSITY ===
            slice_thickness = 0.5
            
            # XY slice (z ~ 0)
            z_mask = np.abs(z - com[2]) < slice_thickness
            if np.sum(z_mask) > 10:
                sc = ax_xy_dens.scatter(x[z_mask], y[z_mask], c=rho[z_mask]/rho_0, 
                                        s=2, cmap='viridis', norm=LogNorm(0.1, 100))
                ax_xy_dens.plot(0, 0, 'r*', ms=15, label='IMBH')
            ax_xy_dens.set_xlim(xlim)
            ax_xy_dens.set_ylim(ylim)
            ax_xy_dens.set_xlabel('X [pc]')
            ax_xy_dens.set_ylabel('Y [pc]')
            ax_xy_dens.set_title(f'XY Density (z≈{com[2]:.1f})')
            ax_xy_dens.set_aspect('equal')
            ax_xy_dens.legend(loc='upper right')
            
            # XZ slice (y ~ COM_y)
            y_mask = np.abs(y - com[1]) < slice_thickness
            if np.sum(y_mask) > 10:
                ax_xz_dens.scatter(x[y_mask], z[y_mask], c=rho[y_mask]/rho_0, 
                                   s=2, cmap='viridis', norm=LogNorm(0.1, 100))
                ax_xz_dens.plot(0, 0, 'r*', ms=15)
            ax_xz_dens.set_xlim(xlim)
            ax_xz_dens.set_ylim(zlim)
            ax_xz_dens.set_xlabel('X [pc]')
            ax_xz_dens.set_ylabel('Z [pc]')
            ax_xz_dens.set_title(f'XZ Density (y≈{com[1]:.1f})')
            ax_xz_dens.set_aspect('equal')
            
            # YZ slice (x ~ COM_x)
            x_mask = np.abs(x - com[0]) < slice_thickness
            if np.sum(x_mask) > 10:
                ax_yz_dens.scatter(y[x_mask], z[x_mask], c=rho[x_mask]/rho_0, 
                                   s=2, cmap='viridis', norm=LogNorm(0.1, 100))
            ax_yz_dens.set_xlim(ylim)
            ax_yz_dens.set_ylim(zlim)
            ax_yz_dens.set_xlabel('Y [pc]')
            ax_yz_dens.set_ylabel('Z [pc]')
            ax_yz_dens.set_title(f'YZ Density (x≈{com[0]:.1f})')
            ax_yz_dens.set_aspect('equal')
            
            # 3D view
            colors = np.log10(rho / rho_0)
            ax_3d.scatter(x, y, z, c=colors, s=1, cmap='viridis', alpha=0.5)
            ax_3d.scatter(0, 0, 0, c='red', s=200, marker='*', label='IMBH')
            ax_3d.set_xlabel('X [pc]')
            ax_3d.set_ylabel('Y [pc]')
            ax_3d.set_zlabel('Z [pc]')
            ax_3d.set_title('3D View (log ρ/ρ₀)')
            
            # === Row 1: Cross-section TEMPERATURE ===
            if np.sum(z_mask) > 10:
                ax_xy_temp.scatter(x[z_mask], y[z_mask], c=T[z_mask], 
                                   s=2, cmap='hot', norm=LogNorm(10, 1e6))
                ax_xy_temp.plot(0, 0, 'c*', ms=15)
            ax_xy_temp.set_xlim(xlim)
            ax_xy_temp.set_ylim(ylim)
            ax_xy_temp.set_xlabel('X [pc]')
            ax_xy_temp.set_ylabel('Y [pc]')
            ax_xy_temp.set_title('XY Temperature [K]')
            ax_xy_temp.set_aspect('equal')
            
            if np.sum(y_mask) > 10:
                ax_xz_temp.scatter(x[y_mask], z[y_mask], c=T[y_mask], 
                                   s=2, cmap='hot', norm=LogNorm(10, 1e6))
                ax_xz_temp.plot(0, 0, 'c*', ms=15)
            ax_xz_temp.set_xlim(xlim)
            ax_xz_temp.set_ylim(zlim)
            ax_xz_temp.set_xlabel('X [pc]')
            ax_xz_temp.set_ylabel('Z [pc]')
            ax_xz_temp.set_title('XZ Temperature [K]')
            ax_xz_temp.set_aspect('equal')
            
            if np.sum(x_mask) > 10:
                ax_yz_temp.scatter(y[x_mask], z[x_mask], c=T[x_mask], 
                                   s=2, cmap='hot', norm=LogNorm(10, 1e6))
            ax_yz_temp.set_xlim(ylim)
            ax_yz_temp.set_ylim(zlim)
            ax_yz_temp.set_xlabel('Y [pc]')
            ax_yz_temp.set_ylabel('Z [pc]')
            ax_yz_temp.set_title('YZ Temperature [K]')
            ax_yz_temp.set_aspect('equal')
            
            # Mach histogram
            ax_mach_hist.hist(mach, bins=50, range=(0, 5), alpha=0.7, color='purple')
            ax_mach_hist.axvline(1.0, color='red', ls='--', lw=2, label='M=1 (sonic)')
            ax_mach_hist.set_xlabel('Mach Number')
            ax_mach_hist.set_ylabel('Count')
            ax_mach_hist.set_title(f'Mach Distribution (max={mach.max():.1f})')
            ax_mach_hist.legend()
            
            # === Row 2: SHOCK views ===
            shock_colors = np.where(is_shocked, 'red', 'blue')
            
            if np.sum(z_mask) > 10:
                ax_xy_shock.scatter(x[z_mask], y[z_mask], c=shock_colors[z_mask], s=2, alpha=0.6)
                ax_xy_shock.plot(0, 0, 'g*', ms=15)
            ax_xy_shock.set_xlim(xlim)
            ax_xy_shock.set_ylim(ylim)
            ax_xy_shock.set_xlabel('X [pc]')
            ax_xy_shock.set_ylabel('Y [pc]')
            ax_xy_shock.set_title(f'XY Shocks ({np.sum(is_shocked & z_mask)}/{np.sum(z_mask)} shocked)')
            ax_xy_shock.set_aspect('equal')
            
            if np.sum(y_mask) > 10:
                ax_xz_shock.scatter(x[y_mask], z[y_mask], c=shock_colors[y_mask], s=2, alpha=0.6)
                ax_xz_shock.plot(0, 0, 'g*', ms=15)
            ax_xz_shock.set_xlim(xlim)
            ax_xz_shock.set_ylim(zlim)
            ax_xz_shock.set_xlabel('X [pc]')
            ax_xz_shock.set_ylabel('Z [pc]')
            ax_xz_shock.set_title('XZ Shocks')
            ax_xz_shock.set_aspect('equal')
            
            if np.sum(x_mask) > 10:
                ax_yz_shock.scatter(y[x_mask], z[x_mask], c=shock_colors[x_mask], s=2, alpha=0.6)
            ax_yz_shock.set_xlim(ylim)
            ax_yz_shock.set_ylim(zlim)
            ax_yz_shock.set_xlabel('Y [pc]')
            ax_yz_shock.set_ylabel('Z [pc]')
            ax_yz_shock.set_title('YZ Shocks')
            ax_yz_shock.set_aspect('equal')
            
            # Info panel
            ax_info.axis('off')
            info_text = f"""IMBH Tidal Disruption (Concatenated)
---------------------------------------
Frame: {frame + 1}/{len(snapshot_files)}
Segment: {segment + 1} of {len(input_dirs)}
Time: {time:.3f} ({time_myr:.2f} Myr)

M_BH = {bh_mass:.0e} Msun

Particles: {len(x):,}
Shocked: {np.sum(is_shocked):,} ({100*np.sum(is_shocked)/len(x):.1f}%)

rho/rho0: [{rho.min()/rho_0:.2f}, {rho.max()/rho_0:.1f}]
T: [{T.min():.0f}, {T.max():.0e}] K
Mach: [0, {mach.max():.1f}]

COM: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f}) pc
---------------------------------------
Legend:
  [RED]  Shocked (rho/rho0>4 or M>1.5)
  [BLUE] Unshocked
  [*]    IMBH location"""
            ax_info.text(0.05, 0.95, info_text, transform=ax_info.transAxes, fontsize=10,
                        verticalalignment='top', family='monospace',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
        else:  # single mode
            # === 3D view ===
            colors = np.log10(rho / rho_0)
            ax_3d.scatter(x, y, z, c=colors, s=1, cmap='viridis', alpha=0.5)
            ax_3d.scatter(0, 0, 0, c='red', s=200, marker='*', label='IMBH')
            ax_3d.set_xlabel('X [pc]')
            ax_3d.set_ylabel('Y [pc]')
            ax_3d.set_zlabel('Z [pc]')
            ax_3d.set_title(f'3D View (log ρ/ρ₀)\nt = {time_myr:.2f} Myr')
            
            # === XY projection ===
            ax_xy.scatter(x, y, c=np.log10(rho/rho_0), s=1, cmap='viridis', alpha=0.5)
            ax_xy.plot(0, 0, 'r*', ms=15, label='IMBH')
            ax_xy.set_xlim(xlim)
            ax_xy.set_ylim(ylim)
            ax_xy.set_xlabel('X [pc]')
            ax_xy.set_ylabel('Y [pc]')
            ax_xy.set_title('XY Projection')
            ax_xy.set_aspect('equal')
            ax_xy.legend()
            
            # === XZ projection ===
            ax_xz.scatter(x, z, c=np.log10(rho/rho_0), s=1, cmap='viridis', alpha=0.5)
            ax_xz.plot(0, 0, 'r*', ms=15)
            ax_xz.set_xlim(xlim)
            ax_xz.set_ylim(zlim)
            ax_xz.set_xlabel('X [pc]')
            ax_xz.set_ylabel('Z [pc]')
            ax_xz.set_title('XZ Projection')
            ax_xz.set_aspect('equal')
            
            # === Histograms ===
            ax_dens_hist.hist(np.log10(rho/rho_0), bins=50, alpha=0.7, color='blue')
            ax_dens_hist.set_xlabel('log(ρ/ρ₀)')
            ax_dens_hist.set_ylabel('Count')
            ax_dens_hist.set_title('Density Distribution')
            
            ax_temp_hist.hist(np.log10(T), bins=50, alpha=0.7, color='orange')
            ax_temp_hist.set_xlabel('log(T) [K]')
            ax_temp_hist.set_ylabel('Count')
            ax_temp_hist.set_title('Temperature Distribution')
            
            ax_mach_hist.hist(mach, bins=50, range=(0, 5), alpha=0.7, color='purple')
            ax_mach_hist.axvline(1.0, color='red', ls='--', lw=2, label='M=1')
            ax_mach_hist.set_xlabel('Mach Number')
            ax_mach_hist.set_ylabel('Count')
            ax_mach_hist.set_title('Mach Distribution')
            ax_mach_hist.legend()
            
            # === Energy (placeholder) ===
            ax_energy.axis('off')
            ax_energy.text(0.5, 0.5, 'Energy plot\n(if available)', ha='center', va='center',
                          transform=ax_energy.transAxes, fontsize=12)
            
            # === COM trajectory ===
            if len(com_history) > 1:
                com_arr = np.array(com_history)
                ax_com.plot(com_arr[:, 0], com_arr[:, 1], 'b-', lw=1.5, alpha=0.7)
                ax_com.plot(com_arr[-1, 0], com_arr[-1, 1], 'ro', ms=8)
                ax_com.plot(0, 0, 'r*', ms=15, label='IMBH')
            ax_com.set_xlabel('X [pc]')
            ax_com.set_ylabel('Y [pc]')
            ax_com.set_title('COM Trajectory')
            ax_com.set_aspect('equal')
            ax_com.legend()
            
            # === Info panel ===
            ax_info.axis('off')
            info_text = f"""IMBH Tidal Disruption
━━━━━━━━━━━━━━━━━━━━━━
Frame: {frame + 1}/{len(snapshot_files)}
Segment: {segment + 1}/{len(input_dirs)}
Time: {time:.3f} ({time_myr:.2f} Myr)

M_BH = {bh_mass:.0e} M☉
Particles: {len(x):,}

ρ/ρ₀: [{rho.min()/rho_0:.2f}, {rho.max()/rho_0:.1f}]
T: [{T.min():.0f}, {T.max():.0e}] K

COM: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})"""
            ax_info.text(0.05, 0.95, info_text, transform=ax_info.transAxes, fontsize=11,
                        verticalalignment='top', family='monospace',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        if (frame + 1) % 10 == 0:
            print(f"   Frame {frame + 1}/{len(snapshot_files)}")
    
    # Create animation
    print(f"\nCreating animation with {len(snapshot_files)} frames...")
    anim = FuncAnimation(fig, update, frames=len(snapshot_files), interval=1000 // fps, repeat=True)
    
    # Determine output filename
    output_name = f"concatenated_{mode}.gif"
    output_file = output_path / output_name
    
    print(f"Saving to {output_file}...")
    writer = PillowWriter(fps=fps)
    anim.save(str(output_file), writer=writer, dpi=dpi)
    
    plt.close(fig)
    
    print(f"✅ Animation saved: {output_file}")
    return True


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='IMBH Tidal Disruption Animation - Concatenated Multi-Segment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Concatenates snapshots from multiple simulation segments into a single animation.
Useful for combining original runs with continuation runs.

Examples:
  # Combine original + continuation
  python animate_concatenated.py \\
      --input-dirs results/original results/continuation/t2to4 \\
      --output-dir results/combined \\
      --bh-mass 100000

  # With tidal mode
  python animate_concatenated.py \\
      --input-dirs results/nocool results/nocool_cont \\
      --output-dir results/nocool_full \\
      --mode tidal --bh-mass 100000
        """
    )
    
    parser.add_argument('--input-dirs', nargs='+', required=True,
                        help='Input directories with snapshots (in chronological order)')
    parser.add_argument('--output-dir', required=True, 
                        help='Output directory for animation')
    parser.add_argument('--bh-mass', type=float, default=1e5, 
                        help='IMBH mass in M_sun (default: 1e5)')
    parser.add_argument('--fps', type=int, default=10, 
                        help='Frames per second (default: 10)')
    parser.add_argument('--mode', choices=['single', 'tidal'], default='single', 
                        help='Visualization mode')
    parser.add_argument('--dpi', type=int, default=100, 
                        help='Output DPI (default: 100)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("  IMBH Tidal Disruption - Concatenated Animation")
    print("=" * 70)
    print(f"Input directories:")
    for i, d in enumerate(args.input_dirs):
        print(f"  {i + 1}. {d}")
    print(f"Output: {args.output_dir}")
    print(f"Mode:   {args.mode}")
    print(f"IMBH:   {args.bh_mass:.0e} M_sun")
    print(f"FPS:    {args.fps}")
    print(f"DPI:    {args.dpi}")
    print()
    
    success = create_concatenated_animation(
        args.input_dirs, args.output_dir, args.bh_mass, 
        args.fps, args.mode, args.dpi
    )
    
    if success:
        print("\n✅ Concatenated animation complete!")
        return 0
    else:
        print("\n❌ Animation failed!")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
