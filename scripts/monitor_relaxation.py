#!/usr/bin/env python3
"""
Monitor K&I 2000 Bonnor-Ebert relaxation convergence.
Compares SPH density to analytical profile as relaxation progresses.
"""

import numpy as np
import pandas as pd
import glob
import os
import sys
import time

# Paths
RESULTS_DIR = "simulations/astrophysics/imbh_cloud/results/hvcc_10k/relaxation"
PROFILE_FILE = os.path.join(RESULTS_DIR, "ki2000_profile.dat")

def load_analytical_profile(filename):
    """Load the K&I 2000 analytical profile."""
    r_table = []
    rho_table = []

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                r_table.append(float(parts[0]))
                rho_table.append(float(parts[1]))

    return np.array(r_table), np.array(rho_table)

def interpolate_analytical_density(r, r_table, rho_table):
    """Interpolate analytical density at radius r."""
    return np.interp(r, r_table, rho_table, left=rho_table[0], right=rho_table[-1])

def analyze_snapshot(csv_file, r_table, rho_table, R_cloud):
    """Analyze a single snapshot and return density error statistics."""
    try:
        df = pd.read_csv(csv_file, comment='#')
    except Exception as e:
        return None

    # Filter cloud particles (not ghost)
    if 'is_ghost' in df.columns:
        cloud = df[df['is_ghost'] == 0].copy()
    else:
        cloud = df.copy()

    # Compute radius
    cloud['r'] = np.sqrt(cloud['pos_x']**2 + cloud['pos_y']**2 + cloud['pos_z']**2)

    # Filter particles within R_cloud
    cloud = cloud[cloud['r'] <= R_cloud * 1.01]

    if len(cloud) == 0:
        return None

    # Get analytical density at each particle position
    rho_analytical = interpolate_analytical_density(cloud['r'].values, r_table, rho_table)
    rho_sph = cloud['dens'].values

    # Compute relative error
    rel_error = np.abs(rho_sph - rho_analytical) / rho_analytical

    # Statistics
    mean_error = np.mean(rel_error) * 100
    max_error = np.max(rel_error) * 100
    std_error = np.std(rel_error) * 100

    # RMS error
    rms_error = np.sqrt(np.mean((rho_sph - rho_analytical)**2 / rho_analytical**2)) * 100

    return {
        'n_particles': len(cloud),
        'mean_error': mean_error,
        'max_error': max_error,
        'std_error': std_error,
        'rms_error': rms_error
    }

def main():
    print("=" * 70)
    print("K&I 2000 Bonnor-Ebert Relaxation Monitor")
    print("=" * 70)

    # Load analytical profile
    if not os.path.exists(PROFILE_FILE):
        print(f"ERROR: Profile file not found: {PROFILE_FILE}")
        sys.exit(1)

    r_table, rho_table = load_analytical_profile(PROFILE_FILE)
    R_cloud = r_table[-1]

    print(f"Analytical profile: {len(r_table)} points, R_cloud = {R_cloud:.6f}")
    print(f"rho_center = {rho_table[0]:.4f}, rho_edge = {rho_table[-1]:.4f}")
    print()

    # Find all snapshots
    snapshot_pattern = os.path.join(RESULTS_DIR, "snapshot_*.csv")

    print(f"{'Snapshot':<20} {'Step':<10} {'Particles':<12} {'Mean Err%':<12} {'Max Err%':<12} {'RMS Err%':<12}")
    print("-" * 78)

    processed = set()

    while True:
        snapshots = sorted(glob.glob(snapshot_pattern))

        for snap in snapshots:
            if snap in processed:
                continue

            basename = os.path.basename(snap)
            step_str = basename.replace('snapshot_', '').replace('.csv', '')
            try:
                step = int(step_str) * 1000  # Assuming output every 1000 steps
            except:
                step = 0

            result = analyze_snapshot(snap, r_table, rho_table, R_cloud)

            if result:
                print(f"{basename:<20} {step:<10} {result['n_particles']:<12} "
                      f"{result['mean_error']:<12.2f} {result['max_error']:<12.2f} "
                      f"{result['rms_error']:<12.2f}")
                sys.stdout.flush()

            processed.add(snap)

        # Check if simulation is done (look for many snapshots or specific indicator)
        if len(snapshots) >= 100:  # 100 snapshots = 100000 steps
            print("\nRelaxation complete!")
            break

        # Wait and check again
        time.sleep(10)

if __name__ == "__main__":
    main()
