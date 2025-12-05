#!/usr/bin/env python3
"""
Compare Self-Gravitating Hydrostatic Equilibrium: With vs Without Grad-H

This script analyzes self-gravitating cloud simulations to show how the
grad-h correction prevents artificial core collapse.

Usage:
    python3 compare_gradh_selfgrav.py --with-gradh <dir1> --no-gradh <dir2> --output <outdir>
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Compare self-gravitating grad-h simulations')
parser.add_argument('--with-gradh', required=True, help='Directory with grad-h results')
parser.add_argument('--no-gradh', required=True, help='Directory without grad-h results')
parser.add_argument('--output', required=True, help='Output directory for figures')
args = parser.parse_args()

# Configuration
with_gradh_dir = args.with_gradh
no_gradh_dir = args.no_gradh
output_dir = args.output

cases = ['with_gradh', 'no_gradh']
case_dirs = [with_gradh_dir, no_gradh_dir]
case_labels = ['With Grad-H (Stable)', 'Without Grad-H (Collapse)']
case_colors = ['#0173B2', '#D55E00']

print('=' * 70)
print('Self-Gravitating Grad-H Study: Hydrostatic Equilibrium Comparison')
print('=' * 70)

os.makedirs(output_dir, exist_ok=True)


def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    with open(filename, 'r') as f:
        metadata = {}
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        f.seek(0)
        lines = [ln for ln in f.readlines() if not ln.startswith('#')]
        header = lines[0].strip().split(',')
        
        for col_name in header:
            data[col_name] = []
        
        for line in lines[1:]:
            values = line.strip().split(',')
            for i, col_name in enumerate(header):
                try:
                    data[col_name].append(float(values[i]))
                except (ValueError, IndexError):
                    pass
    
    for key in data:
        data[key] = np.array(data[key])
    
    return data, metadata


def compute_radial_profile(data, n_bins=50):
    """Compute radial density profile from center of mass"""
    x = data.get('pos_x', np.array([]))
    y = data.get('pos_y', np.array([]))
    z = data.get('pos_z', np.array([]))
    dens = data.get('dens', np.array([]))
    mass = data.get('mass', np.ones_like(dens))
    
    if len(x) == 0:
        return None, None, None
    
    # Center of mass
    total_mass = np.sum(mass)
    x_com = np.sum(mass * x) / total_mass
    y_com = np.sum(mass * y) / total_mass
    z_com = np.sum(mass * z) / total_mass if len(z) > 0 else 0
    
    # Radial distance from COM
    if len(z) > 0:
        r = np.sqrt((x - x_com)**2 + (y - y_com)**2 + (z - z_com)**2)
    else:
        r = np.sqrt((x - x_com)**2 + (y - y_com)**2)
    
    r_max = np.percentile(r, 95)
    r_bins = np.linspace(0, r_max, n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    dens_profile = []
    dens_std = []
    
    for i in range(n_bins):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask) > 0:
            dens_profile.append(np.mean(dens[mask]))
            dens_std.append(np.std(dens[mask]))
        else:
            dens_profile.append(np.nan)
            dens_std.append(0)
    
    return r_centers, np.array(dens_profile), np.array(dens_std)


def compute_central_density(data, r_frac=0.1):
    """Compute central density (mean within r_frac of cloud radius)"""
    x = data.get('pos_x', np.array([]))
    y = data.get('pos_y', np.array([]))
    z = data.get('pos_z', np.array([]))
    dens = data.get('dens', np.array([]))
    mass = data.get('mass', np.ones_like(dens))
    
    if len(x) == 0:
        return None
    
    total_mass = np.sum(mass)
    x_com = np.sum(mass * x) / total_mass
    y_com = np.sum(mass * y) / total_mass
    z_com = np.sum(mass * z) / total_mass if len(z) > 0 else 0
    
    if len(z) > 0:
        r = np.sqrt((x - x_com)**2 + (y - y_com)**2 + (z - z_com)**2)
    else:
        r = np.sqrt((x - x_com)**2 + (y - y_com)**2)
    
    r_cloud = np.percentile(r, 90)
    central_mask = r < r_frac * r_cloud
    
    if np.sum(central_mask) > 0:
        return np.mean(dens[central_mask])
    return None


def compute_cloud_radius(data, frac=0.9):
    """Compute radius containing frac of mass"""
    x = data.get('pos_x', np.array([]))
    y = data.get('pos_y', np.array([]))
    z = data.get('pos_z', np.array([]))
    mass = data.get('mass', np.ones(len(x)))
    
    if len(x) == 0:
        return None
    
    total_mass = np.sum(mass)
    x_com = np.sum(mass * x) / total_mass
    y_com = np.sum(mass * y) / total_mass
    z_com = np.sum(mass * z) / total_mass if len(z) > 0 else 0
    
    if len(z) > 0:
        r = np.sqrt((x - x_com)**2 + (y - y_com)**2 + (z - z_com)**2)
    else:
        r = np.sqrt((x - x_com)**2 + (y - y_com)**2)
    
    # Sort by radius
    sort_idx = np.argsort(r)
    r_sorted = r[sort_idx]
    mass_sorted = mass[sort_idx]
    cumulative_mass = np.cumsum(mass_sorted) / total_mass
    
    idx = np.searchsorted(cumulative_mass, frac)
    if idx < len(r_sorted):
        return r_sorted[idx]
    return r_sorted[-1]


# Scan for snapshots
print("Scanning for snapshot files...")
case_snapshots = {}

for case, case_dir in zip(cases, case_dirs):
    if not os.path.exists(case_dir):
        print(f"ERROR: Directory not found: {case_dir}")
        sys.exit(1)
    
    snapshots = sorted(glob.glob(os.path.join(case_dir, "snapshot_*.csv")))
    snapshot_numbers = []
    for snap in snapshots:
        num_str = os.path.basename(snap).replace('snapshot_', '').replace('.csv', '')
        try:
            snapshot_numbers.append(int(num_str))
        except ValueError:
            pass
    
    case_snapshots[case] = {'files': snapshots, 'numbers': snapshot_numbers}
    print(f"  {case}: {len(snapshot_numbers)} snapshots")

common_numbers = sorted(
    set(case_snapshots['with_gradh']['numbers']) & set(case_snapshots['no_gradh']['numbers'])
)
print(f"Found {len(common_numbers)} common snapshots")

# ================================================================================
# Plot 1: Central density evolution
# ================================================================================
print("\nGenerating central density evolution...")

fig, ax = plt.subplots(figsize=(10, 6))

for case, case_dir, label, color in zip(cases, case_dirs, case_labels, case_colors):
    times = []
    central_densities = []
    
    for snap_num in common_numbers:
        snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, meta = load_snapshot(snap_file)
            t = float(meta.get('time', 0))
            rho_c = compute_central_density(data)
            
            if rho_c is not None:
                times.append(t)
                central_densities.append(rho_c)
    
    if times:
        times = np.array(times)
        central_densities = np.array(central_densities)
        if central_densities[0] > 0:
            central_densities /= central_densities[0]
        ax.plot(times, central_densities, '-', label=label, color=color, linewidth=2)

ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Initial')
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('Central Density / Initial', fontsize=12)
ax.set_title('Self-Gravitating Cloud: Central Density Evolution', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selfgrav_gradh_central_density.png'), dpi=150)
plt.close()
print("  Saved: selfgrav_gradh_central_density.png")

# ================================================================================
# Plot 2: Cloud radius evolution
# ================================================================================
print("Generating cloud radius evolution...")

fig, ax = plt.subplots(figsize=(10, 6))

for case, case_dir, label, color in zip(cases, case_dirs, case_labels, case_colors):
    times = []
    radii = []
    
    for snap_num in common_numbers:
        snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, meta = load_snapshot(snap_file)
            t = float(meta.get('time', 0))
            r = compute_cloud_radius(data)
            
            if r is not None:
                times.append(t)
                radii.append(r)
    
    if times:
        times = np.array(times)
        radii = np.array(radii)
        if radii[0] > 0:
            radii /= radii[0]
        ax.plot(times, radii, '-', label=label, color=color, linewidth=2)

ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Initial')
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('R(90% mass) / Initial', fontsize=12)
ax.set_title('Self-Gravitating Cloud: Radius Evolution', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selfgrav_gradh_radius.png'), dpi=150)
plt.close()
print("  Saved: selfgrav_gradh_radius.png")

# ================================================================================
# Plot 3: Radial density profile at multiple times
# ================================================================================
print("Generating radial profile evolution...")

if len(common_numbers) >= 4:
    plot_indices = [0, len(common_numbers)//3, 2*len(common_numbers)//3, -1]
else:
    plot_indices = list(range(len(common_numbers)))

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.ravel()

for ax_idx, snap_idx in enumerate(plot_indices[:4]):
    ax = axes[ax_idx]
    snap_num = common_numbers[snap_idx]
    
    for case, case_dir, label, color in zip(cases, case_dirs, case_labels, case_colors):
        snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, meta = load_snapshot(snap_file)
            r, dens_profile, dens_std = compute_radial_profile(data)
            
            if r is not None:
                ax.plot(r, dens_profile, '-', color=color, linewidth=2, label=label)
                ax.fill_between(r, dens_profile - dens_std, dens_profile + dens_std, 
                               color=color, alpha=0.2)
            
            time = float(meta.get('time', 0))
    
    ax.set_xlabel('Radius', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f't = {time:.2f}', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

plt.suptitle('Radial Density Profile Evolution', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'selfgrav_gradh_profile_evolution.png'), dpi=150)
plt.close()
print("  Saved: selfgrav_gradh_profile_evolution.png")

print()
print('=' * 70)
print('All self-gravitating comparison plots generated!')
print(f'Output: {output_dir}')
print('=' * 70)
