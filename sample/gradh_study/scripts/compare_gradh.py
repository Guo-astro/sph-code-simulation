#!/usr/bin/env python3
"""
Compare Hydrostatic Equilibrium: With vs Without Grad-H Correction

This script compares two hydrostatic simulations to demonstrate the importance
of the grad-h correction term (Ω_i) in preventing artificial core collapse.

Theory:
    The grad-h correction factor Ω_i = 1 / (1 + (h/Dρ) * dρ/dh) accounts for
    variable smoothing length in the kernel gradient. Without this correction,
    the density summation becomes inconsistent with pressure gradients.

Usage:
    python3 compare_gradh.py --with-gradh <dir1> --no-gradh <dir2> --output <outdir>
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import argparse
from pathlib import Path

# Parse arguments
parser = argparse.ArgumentParser(description='Compare grad-h vs no-grad-h simulations')
parser.add_argument('--with-gradh', required=True, help='Directory with grad-h results')
parser.add_argument('--no-gradh', required=True, help='Directory without grad-h results')
parser.add_argument('--output', required=True, help='Output directory for figures')
args = parser.parse_args()

# Configuration
with_gradh_dir = args.with_gradh
no_gradh_dir = args.no_gradh
output_dir = args.output

# Labels and colors
cases = ['with_gradh', 'no_gradh']
case_dirs = [with_gradh_dir, no_gradh_dir]
case_labels = ['With Grad-H (Ω ≠ 1)', 'Without Grad-H (Ω = 1)']
case_colors = ['#0173B2', '#D55E00']  # Blue, Red-orange

print('=' * 70)
print('Grad-H Correction Study: Hydrostatic Equilibrium Comparison')
print('=' * 70)
print(f'With grad-h:    {with_gradh_dir}')
print(f'Without grad-h: {no_gradh_dir}')
print(f'Output:         {output_dir}')
print()

# Create output directory
os.makedirs(output_dir, exist_ok=True)


def load_snapshot(filename):
    """Load a single snapshot CSV file"""
    data = {}
    with open(filename, 'r') as f:
        # Read metadata lines (start with #)
        metadata = {}
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)
                    metadata[key.strip()] = value.strip()
            else:
                break
        
        # Read header and data
        f.seek(0)
        lines = [l for l in f.readlines() if not l.startswith('#')]
        header = lines[0].strip().split(',')
        
        # Parse data
        for col_name in header:
            data[col_name] = []
        
        for line in lines[1:]:
            values = line.strip().split(',')
            for i, col_name in enumerate(header):
                try:
                    data[col_name].append(float(values[i]))
                except (ValueError, IndexError):
                    pass
    
    # Convert to numpy arrays
    for key in data:
        data[key] = np.array(data[key])
    
    return data, metadata


def load_energy_file(directory):
    """Load energy conservation file"""
    energy_file = os.path.join(directory, 'energy.csv')
    if not os.path.exists(energy_file):
        return None
    
    data = {'time': [], 'kinetic': [], 'thermal': [], 'total': []}
    with open(energy_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split(',')
            if len(parts) >= 4:
                try:
                    data['time'].append(float(parts[0]))
                    data['kinetic'].append(float(parts[1]))
                    data['thermal'].append(float(parts[2]))
                    data['total'].append(float(parts[3]))
                except ValueError:
                    pass
    
    for key in data:
        data[key] = np.array(data[key])
    
    return data


def compute_central_density(data):
    """Compute central density (particles near center)"""
    if 'pos_x' not in data or 'pos_y' not in data:
        return None
    
    x = data['pos_x']
    y = data['pos_y']
    dens = data['dens']
    
    # Find particles within 10% of domain
    r = np.sqrt(x**2 + y**2)
    central_mask = r < 0.1
    
    if np.sum(central_mask) > 0:
        return np.mean(dens[central_mask])
    return None


# Scan for snapshots
print("Scanning for snapshot files...")
case_snapshots = {}

for case, case_dir in zip(cases, case_dirs):
    if not os.path.exists(case_dir):
        print(f"  ERROR: Directory not found: {case_dir}")
        sys.exit(1)
    
    snapshots = sorted(glob.glob(os.path.join(case_dir, "snapshot_*.csv")))
    if not snapshots:
        print(f"  ERROR: No snapshots found in {case_dir}")
        sys.exit(1)
    
    snapshot_numbers = []
    for snap in snapshots:
        basename = os.path.basename(snap)
        num_str = basename.replace('snapshot_', '').replace('.csv', '')
        try:
            snapshot_numbers.append(int(num_str))
        except ValueError:
            pass
    
    case_snapshots[case] = {
        'files': snapshots,
        'numbers': snapshot_numbers
    }
    print(f"  {case:15s}: {len(snapshot_numbers)} snapshots")

# Find common snapshots
common_numbers = set(case_snapshots['with_gradh']['numbers']) & set(case_snapshots['no_gradh']['numbers'])
common_numbers = sorted(common_numbers)
print(f"\nFound {len(common_numbers)} common snapshots")

# ================================================================================
# Plot 1: Side-by-side particle distribution comparison
# ================================================================================
print("\nGenerating particle distribution comparison...")

# Select snapshots to show
if len(common_numbers) >= 4:
    plot_indices = [0, len(common_numbers)//3, 2*len(common_numbers)//3, -1]
    plot_nums = [common_numbers[i] for i in plot_indices]
else:
    plot_nums = common_numbers[:4]

fig, axes = plt.subplots(2, len(plot_nums), figsize=(4*len(plot_nums), 8))
fig.suptitle('Grad-H Correction Effect on Hydrostatic Equilibrium', fontsize=14, fontweight='bold')

for col, snap_num in enumerate(plot_nums):
    for row, (case, case_dir, label, color) in enumerate(zip(cases, case_dirs, case_labels, case_colors)):
        ax = axes[row, col] if len(plot_nums) > 1 else axes[row]
        
        snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
        if os.path.exists(snap_file):
            data, meta = load_snapshot(snap_file)
            
            x = data.get('pos_x', [])
            y = data.get('pos_y', [])
            dens = data.get('dens', np.ones_like(x))
            
            scatter = ax.scatter(x, y, c=dens, s=1, cmap='viridis', 
                               vmin=np.percentile(dens, 5), vmax=np.percentile(dens, 95))
            
            time = float(meta.get('time', 0))
            if row == 0:
                ax.set_title(f't = {time:.2f}', fontsize=11)
            
            if col == 0:
                ax.set_ylabel(f'{label}\ny', fontsize=10)
            
            ax.set_xlim(-0.5, 0.5)
            ax.set_ylim(-0.5, 0.5)
            ax.set_aspect('equal')
            
            if row == 1:
                ax.set_xlabel('x', fontsize=10)

plt.tight_layout()
fig.colorbar(scatter, ax=axes.ravel().tolist(), label='Density', shrink=0.6)
plt.savefig(os.path.join(output_dir, 'hydrostatic_gradh_particles.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: hydrostatic_gradh_particles.png")

# ================================================================================
# Plot 2: Central density evolution
# ================================================================================
print("Generating central density evolution plot...")

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
        # Normalize by initial
        if central_densities[0] > 0:
            central_densities /= central_densities[0]
        
        ax.plot(times, central_densities, '-', label=label, color=color, linewidth=2)

ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Initial')
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('Central Density / Initial', fontsize=12)
ax.set_title('Central Density Evolution: Grad-H vs No Grad-H', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, None)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hydrostatic_gradh_density_evolution.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: hydrostatic_gradh_density_evolution.png")

# ================================================================================
# Plot 3: Density profile comparison at final time
# ================================================================================
print("Generating density profile comparison...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

snap_num = common_numbers[-1]

for ax_idx, (case, case_dir, label, color) in enumerate(zip(cases, case_dirs, case_labels, case_colors)):
    ax = axes[ax_idx]
    
    snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
    if os.path.exists(snap_file):
        data, meta = load_snapshot(snap_file)
        
        x = data.get('pos_x', [])
        y = data.get('pos_y', [])
        dens = data.get('dens', [])
        
        r = np.sqrt(np.array(x)**2 + np.array(y)**2)
        
        # Bin the density profile
        r_bins = np.linspace(0, 0.5, 30)
        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        dens_profile = []
        dens_std = []
        
        for i in range(len(r_bins) - 1):
            mask = (r >= r_bins[i]) & (r < r_bins[i+1])
            if np.sum(mask) > 0:
                dens_profile.append(np.mean(np.array(dens)[mask]))
                dens_std.append(np.std(np.array(dens)[mask]))
            else:
                dens_profile.append(np.nan)
                dens_std.append(0)
        
        dens_profile = np.array(dens_profile)
        dens_std = np.array(dens_std)
        
        ax.plot(r_centers, dens_profile, '-', color=color, linewidth=2, label='Mean')
        ax.fill_between(r_centers, dens_profile - dens_std, dens_profile + dens_std, 
                       color=color, alpha=0.3, label='±1σ')
        
        time = float(meta.get('time', 0))
        ax.set_title(f'{label}\nt = {time:.2f}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Radius', fontsize=11)
        ax.set_ylabel('Density', fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 0.5)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hydrostatic_gradh_density_profile.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: hydrostatic_gradh_density_profile.png")

# ================================================================================
# Plot 4: Energy conservation comparison
# ================================================================================
print("Generating energy conservation comparison...")

fig, ax = plt.subplots(figsize=(10, 6))

for case, case_dir, label, color in zip(cases, case_dirs, case_labels, case_colors):
    energy_data = load_energy_file(case_dir)
    
    if energy_data is not None and len(energy_data['time']) > 0:
        t = energy_data['time']
        E_total = energy_data['total']
        
        # Relative energy change
        if E_total[0] != 0:
            dE = (E_total - E_total[0]) / np.abs(E_total[0])
        else:
            dE = E_total - E_total[0]
        
        ax.plot(t, dE, '-', label=label, color=color, linewidth=2)

ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('(E - E₀) / |E₀|', fontsize=12)
ax.set_title('Energy Conservation: Grad-H vs No Grad-H', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hydrostatic_gradh_energy.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: hydrostatic_gradh_energy.png")

# ================================================================================
# Plot 5: Grad-h factor distribution
# ================================================================================
print("Generating grad-h factor distribution...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

snap_num = common_numbers[-1]

for ax_idx, (case, case_dir, label, color) in enumerate(zip(cases, case_dirs, case_labels, case_colors)):
    ax = axes[ax_idx]
    
    snap_file = os.path.join(case_dir, f"snapshot_{snap_num:04d}.csv")
    if os.path.exists(snap_file):
        data, meta = load_snapshot(snap_file)
        
        gradh = data.get('gradh', [])
        
        if len(gradh) > 0:
            ax.hist(gradh, bins=50, color=color, alpha=0.7, edgecolor='black')
            ax.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='Ω = 1 (no correction)')
            ax.axvline(x=np.mean(gradh), color='green', linestyle='-', linewidth=2, label=f'Mean = {np.mean(gradh):.4f}')
            
            ax.set_xlabel('Grad-H Factor (Ω)', fontsize=11)
            ax.set_ylabel('Count', fontsize=11)
            ax.set_title(f'{label}', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hydrostatic_gradh_factor_dist.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: hydrostatic_gradh_factor_dist.png")

print()
print('=' * 70)
print('✓ All comparison plots generated!')
print(f'📊 Output directory: {output_dir}')
print('=' * 70)
