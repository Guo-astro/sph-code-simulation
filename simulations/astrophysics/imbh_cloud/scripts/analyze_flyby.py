#!/usr/bin/env python3
"""
Analyze IMBH flyby simulation results.
Tracks BH position, cloud response, and tidal effects.
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

def load_snapshot(filepath):
    """Load CSV snapshot data, skipping comment lines."""
    # Read with pandas, skip comment lines
    df = pd.read_csv(filepath, comment='#')
    return df

def analyze_cloud(df):
    """Calculate cloud properties from particle data."""
    x, y, z = df['pos_x'].values, df['pos_y'].values, df['pos_z'].values
    vx, vy, vz = df['vel_x'].values, df['vel_y'].values, df['vel_z'].values
    
    # Center of mass
    cm = np.array([np.mean(x), np.mean(y), np.mean(z)])
    
    # Bulk velocity
    v_bulk = np.array([np.mean(vx), np.mean(vy), np.mean(vz)])
    
    # Velocity dispersion (internal motion only)
    v_internal_x = vx - v_bulk[0]
    v_internal_y = vy - v_bulk[1]
    v_internal_z = vz - v_bulk[2]
    sigma_v = np.sqrt(np.mean(v_internal_x**2 + v_internal_y**2 + v_internal_z**2))
    
    # Cloud radii from CoM
    r = np.sqrt((x-cm[0])**2 + (y-cm[1])**2 + (z-cm[2])**2)
    r_max = np.max(r)
    r_mean = np.mean(r)
    
    return {
        'cm': cm,
        'v_bulk': v_bulk,
        'sigma_v': sigma_v,
        'r_max': r_max,
        'r_mean': r_mean
    }

def main():
    results_dir = Path("/Users/guo-opt-p148/sph-code-simulation/simulations/astrophysics/imbh_cloud/results/phase3_flyby")
    
    # Find all snapshots
    snapshots = sorted(results_dir.glob("snapshot_*.csv"))
    print(f"Found {len(snapshots)} snapshots")
    
    if len(snapshots) == 0:
        print("No snapshots found!")
        return
    
    # BH orbital parameters
    x0_bh, y0_bh = 15.0, 4.07  # Initial BH position
    v_bh = -80.0  # BH velocity in x-direction (km/s = code units)
    
    print("\n" + "="*80)
    print("IMBH FLYBY ANALYSIS")
    print("="*80)
    
    # Analyze key snapshots
    reference = None
    print(f"\n{'Snap':>5} {'Time':>8} {'BH_x':>8} {'BH_r':>8} {'Cloud_x':>8} {'v_bulk_x':>8} {'sigma_v':>8} {'R_max':>8}")
    print("-"*80)
    
    for snap_path in snapshots[::max(1, len(snapshots)//10)]:  # Sample every ~10%
        snap_num = int(snap_path.stem.split('_')[1])
        
        # Read time from header
        with open(snap_path) as f:
            for line in f:
                if 'Time (code):' in line:
                    t = float(line.split(':')[1].strip())
                    break
        
        # Load and analyze cloud
        data = load_snapshot(snap_path)
        props = analyze_cloud(data)
        
        # Calculate BH position at this time
        x_bh = x0_bh + v_bh * t
        y_bh = y0_bh
        r_bh = np.sqrt(x_bh**2 + y_bh**2)  # Distance from origin (cloud)
        
        if reference is None:
            reference = props
        
        print(f"{snap_num:5d} {t:8.4f} {x_bh:8.2f} {r_bh:8.2f} {props['cm'][0]:8.4f} {props['v_bulk'][0]:8.3f} {props['sigma_v']:8.3f} {props['r_max']:8.3f}")
    
    # Also analyze latest snapshot in detail
    print("\n" + "="*80)
    print("DETAILED ANALYSIS OF LATEST SNAPSHOT")
    print("="*80)
    
    latest = snapshots[-1]
    snap_num = int(latest.stem.split('_')[1])
    
    # Read time
    with open(latest) as f:
        for line in f:
            if 'Time (code):' in line:
                t = float(line.split(':')[1].strip())
                break
    
    data = load_snapshot(latest)
    props = analyze_cloud(data)
    
    x_bh = x0_bh + v_bh * t
    y_bh = y0_bh
    r_bh = np.sqrt(x_bh**2 + y_bh**2)
    
    print(f"\nSnapshot: {snap_num} at t = {t:.4f} (code time)")
    print(f"  Physical time: {t * 0.978:.4f} Myr")
    
    print(f"\nBH Position:")
    print(f"  x_BH = {x_bh:.3f} pc")
    print(f"  y_BH = {y_bh:.3f} pc")
    print(f"  Distance from cloud center: {r_bh:.3f} pc")
    
    # Pericenter info
    t_peri = -x0_bh / v_bh  # Time when x_BH = 0
    r_peri = y_bh  # Pericenter distance = impact parameter
    print(f"\nPericenter passage:")
    print(f"  Pericenter time: t = {t_peri:.4f}")
    if t < t_peri:
        print(f"  Status: APPROACHING (t < t_peri)")
    else:
        print(f"  Status: RECEDING (t > t_peri)")
    print(f"  Pericenter distance: {r_peri:.3f} pc")
    
    print(f"\nCloud Properties:")
    print(f"  Center of mass: ({props['cm'][0]:.4f}, {props['cm'][1]:.4f}, {props['cm'][2]:.4f}) pc")
    print(f"  CoM displacement from origin: {np.sqrt(np.sum(props['cm']**2)):.4f} pc")
    print(f"  Bulk velocity: ({props['v_bulk'][0]:.3f}, {props['v_bulk'][1]:.3f}, {props['v_bulk'][2]:.3f}) km/s")
    print(f"  Speed toward BH: {props['v_bulk'][0]:.3f} km/s")
    print(f"  Internal velocity dispersion: {props['sigma_v']:.3f} km/s")
    print(f"  Maximum radius: {props['r_max']:.3f} pc")
    
    # Compare to initial
    if reference:
        print(f"\nChanges from initial:")
        print(f"  CoM displacement: {np.sqrt(np.sum((props['cm']-reference['cm'])**2)):.4f} pc")
        print(f"  Bulk velocity gained: {np.sqrt(np.sum((props['v_bulk']-reference['v_bulk'])**2)):.3f} km/s")
        print(f"  Dispersion change: {props['sigma_v'] - reference['sigma_v']:.3f} km/s ({(props['sigma_v']/reference['sigma_v']-1)*100:.1f}%)")
        print(f"  Radius change: {(props['r_max']/reference['r_max']-1)*100:.1f}%")
    
    # Tidal effects assessment
    print("\n" + "="*80)
    print("TIDAL EFFECTS ASSESSMENT")
    print("="*80)
    
    M_BH = 100000.0  # Solar masses
    M_cloud = 127.5  # Solar masses
    R_cloud = props['r_max']
    r_tidal = R_cloud * (M_BH / M_cloud) ** (1/3)
    
    print(f"\nTidal radius (cloud would be disrupted if closer): {r_tidal:.2f} pc")
    print(f"Current BH distance: {r_bh:.2f} pc")
    print(f"Pericenter distance: {r_peri:.2f} pc")
    
    if r_bh < r_tidal:
        print(f"⚠️  Cloud is INSIDE tidal radius - expect significant disruption")
    else:
        print(f"✓  Cloud is outside tidal radius")
    
    if r_peri < r_tidal:
        print(f"⚠️  Pericenter is INSIDE tidal radius - disruption during close passage")
    
    # Expected final outcome
    print("\n" + "="*80)
    print("EXPECTED OUTCOME")
    print("="*80)
    print(f"Cloud will gain velocity toward +x direction as BH passes")
    print(f"Target v_los ~ 80 km/s for Oka et al. 2017 HVCC simulation")
    print(f"Current v_x of cloud: {props['v_bulk'][0]:.3f} km/s")
    
if __name__ == '__main__':
    main()
