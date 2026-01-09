#!/usr/bin/env python3
"""Quick visualization of HVCC isothermal sphere relaxation."""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Find snapshot files
results_dir = "../results/hvcc_10k/relaxation"
snapshot_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))

if not snapshot_files:
    print("No snapshot files found!")
    exit(1)

print(f"Found {len(snapshot_files)} snapshots")

# Load first and last snapshots
def load_snapshot(filename):
    # Find header line (starts with 'id,')
    header_line = 0
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('id,'):
                header_line = i
                break

    # Load using pandas-like approach with numpy
    import csv
    with open(filename, 'r') as f:
        for _ in range(header_line):
            next(f)
        reader = csv.DictReader(f)
        rows = list(reader)

    # Convert to structured array
    n = len(rows)
    dtype = [(key, 'f8') for key in rows[0].keys()]
    arr = np.zeros(n, dtype=dtype)
    for i, row in enumerate(rows):
        for key, val in row.items():
            try:
                arr[i][key] = float(val)
            except:
                arr[i][key] = 0.0
    return arr

# Load snapshots
snap0 = load_snapshot(snapshot_files[0])
snap_final = load_snapshot(snapshot_files[-1])

print(f"Initial: {len(snap0)} particles")
print(f"Final: {len(snap_final)} particles")

# Create figure
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Separate cloud and envelope particles
cloud0 = snap0[snap0['is_ghost'] == 0]
ghost0 = snap0[snap0['is_ghost'] == 1]
cloudf = snap_final[snap_final['is_ghost'] == 0]
ghostf = snap_final[snap_final['is_ghost'] == 1]

print(f"Initial: {len(cloud0)} cloud + {len(ghost0)} envelope particles")
print(f"Final: {len(cloudf)} cloud + {len(ghostf)} envelope particles")

# Plot 1: Initial XY scatter with envelope
ax = axes[0, 0]
r0 = np.sqrt(snap0['pos_x']**2 + snap0['pos_y']**2 + snap0['pos_z']**2)
r0_cloud = np.sqrt(cloud0['pos_x']**2 + cloud0['pos_y']**2 + cloud0['pos_z']**2)
r0_ghost = np.sqrt(ghost0['pos_x']**2 + ghost0['pos_y']**2 + ghost0['pos_z']**2) if len(ghost0) > 0 else np.array([])

# Plot cloud particles
sc = ax.scatter(cloud0['pos_x'], cloud0['pos_y'], c=cloud0['dens'], s=1, cmap='viridis', label='Cloud')
# Plot envelope particles in red
if len(ghost0) > 0:
    ax.scatter(ghost0['pos_x'], ghost0['pos_y'], c='red', s=5, marker='x', alpha=0.7, label='Envelope')
ax.set_xlabel('X [pc]')
ax.set_ylabel('Y [pc]')
ax.set_title('Initial: Cloud (color) + Envelope (red x)')
ax.set_aspect('equal')
ax.legend(loc='upper right', fontsize=7)
plt.colorbar(sc, ax=ax, label='Density [M_sun/pc^3]')

# Plot 2: Final XY scatter with envelope
ax = axes[0, 1]
rf = np.sqrt(snap_final['pos_x']**2 + snap_final['pos_y']**2 + snap_final['pos_z']**2)
rf_cloud = np.sqrt(cloudf['pos_x']**2 + cloudf['pos_y']**2 + cloudf['pos_z']**2)
rf_ghost = np.sqrt(ghostf['pos_x']**2 + ghostf['pos_y']**2 + ghostf['pos_z']**2) if len(ghostf) > 0 else np.array([])

sc = ax.scatter(cloudf['pos_x'], cloudf['pos_y'], c=cloudf['dens'], s=1, cmap='viridis', label='Cloud')
if len(ghostf) > 0:
    ax.scatter(ghostf['pos_x'], ghostf['pos_y'], c='red', s=5, marker='x', alpha=0.7, label='Envelope')
ax.set_xlabel('X [pc]')
ax.set_ylabel('Y [pc]')
ax.set_title('Final: Cloud (color) + Envelope (red x)')
ax.set_aspect('equal')
ax.legend(loc='upper right', fontsize=7)
plt.colorbar(sc, ax=ax, label='Density [M_sun/pc^3]')

# Plot 3: Radial density profile comparison with analytical
ax = axes[0, 2]

# Analytical modified isothermal sphere parameters (K&I 2000 equilibrium)
# From IC generation: rho(r) = rho_c / (1 + (r/r_c)^2)
# With n_edge = P_ext/T = 710/10 = 71 cm^-3 for equilibrium
rho_c = 57.5846  # M_sun/pc^3 (n_center = 1000 cm^-3)
n_center = 1000.0  # cm^-3
n_ext = 71.0       # cm^-3 (K&I external at T=10K)
density_contrast = n_center / n_ext  # = 14.08
R_over_rc = np.sqrt(density_contrast - 1.0)  # = 3.62
# Values from actual simulation (M=40 M_sun, K&I equilibrium constraint)
R_cloud = 1.04   # pc (from simulation output)
r_c = 0.288      # pc (from simulation output)

r_analytic = np.linspace(0.01, 2.0, 200)
rho_analytic = rho_c / (1.0 + (r_analytic / r_c)**2)

# Bin simulation data by radius
r_bins = np.linspace(0, 2.0, 80)
r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

# Cloud particles only - Initial
hist0_cloud, _ = np.histogram(r0_cloud, bins=r_bins, weights=cloud0['dens'])
count0_cloud, _ = np.histogram(r0_cloud, bins=r_bins)
rho0_cloud_avg = np.where(count0_cloud > 0, hist0_cloud / count0_cloud, np.nan)

# Envelope particles - Initial
if len(ghost0) > 0:
    hist0_ghost, _ = np.histogram(r0_ghost, bins=r_bins, weights=ghost0['dens'])
    count0_ghost, _ = np.histogram(r0_ghost, bins=r_bins)
    rho0_ghost_avg = np.where(count0_ghost > 0, hist0_ghost / count0_ghost, np.nan)
else:
    rho0_ghost_avg = np.full_like(r_centers, np.nan)

# Cloud particles only - Final
hist_cloud, _ = np.histogram(rf_cloud, bins=r_bins, weights=cloudf['dens'])
count_cloud, _ = np.histogram(rf_cloud, bins=r_bins)
rhof_cloud_avg = np.where(count_cloud > 0, hist_cloud / count_cloud, np.nan)

# Envelope particles - Final
if len(ghostf) > 0:
    histf_ghost, _ = np.histogram(rf_ghost, bins=r_bins, weights=ghostf['dens'])
    countf_ghost, _ = np.histogram(rf_ghost, bins=r_bins)
    rhof_ghost_avg = np.where(countf_ghost > 0, histf_ghost / countf_ghost, np.nan)
else:
    rhof_ghost_avg = np.full_like(r_centers, np.nan)

# K&I external density (for equilibrium)
rho_ext = rho_c / density_contrast  # = rho_c at r=R_cloud

# Plot
ax.semilogy(r_analytic, rho_analytic, 'k-', label='Analytical', linewidth=2.5)
ax.semilogy(r_centers, rho0_cloud_avg, 'b-', label='SPH Initial', linewidth=2, marker='o', markersize=3)
ax.semilogy(r_centers, rhof_cloud_avg, 'r--', label='SPH Final', linewidth=2, marker='s', markersize=3)
ax.axhline(rho_ext, color='orange', linestyle='--', linewidth=1.5, label=f'ρ_ext (n={n_ext:.0f}/cc)')
ax.axvline(R_cloud, color='gray', linestyle=':', linewidth=1.5)
ax.axvline(r_c, color='green', linestyle=':', alpha=0.7)
ax.text(R_cloud+0.02, 50, f'R={R_cloud:.2f}', fontsize=8, color='gray')
ax.text(r_c+0.02, 50, f'r_c={r_c:.2f}', fontsize=8, color='green', alpha=0.7)
ax.set_xlabel('Radius [pc]')
ax.set_ylabel('Density [M_sun/pc^3]')
ax.set_title('Radial Density Profile\n(Modified Isothermal Sphere)')
ax.legend(fontsize=8, loc='upper right')
ax.set_xlim(0, 1.5)
ax.set_ylim(1, 100)
ax.grid(True, alpha=0.3)

# Plot 4: Velocity distribution
ax = axes[1, 0]
v0 = np.sqrt(snap0['vel_x']**2 + snap0['vel_y']**2 + snap0['vel_z']**2)
vf = np.sqrt(snap_final['vel_x']**2 + snap_final['vel_y']**2 + snap_final['vel_z']**2)
ax.hist(v0, bins=50, alpha=0.5, label='Initial', density=True)
ax.hist(vf, bins=50, alpha=0.5, label='Final', density=True)
ax.set_xlabel('Speed [km/s]')
ax.set_ylabel('Probability density')
ax.set_title('Velocity Distribution')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 5: Energy file
ax = axes[1, 1]
try:
    energy_file = f"{results_dir}/energy.dat"
    energy = np.loadtxt(energy_file, skiprows=1)
    if len(energy.shape) == 1:
        energy = energy.reshape(1, -1)
    ax.plot(energy[:, 0], energy[:, 1], 'b-', label='Kinetic', linewidth=2)
    ax.plot(energy[:, 0], energy[:, 2], 'r-', label='Internal', linewidth=2)
    ax.plot(energy[:, 0], energy[:, 3], 'g-', label='Total', linewidth=2)
    ax.set_xlabel('Time [code units]')
    ax.set_ylabel('Energy [code units]')
    ax.set_title('Energy Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
except Exception as e:
    ax.text(0.5, 0.5, f'Energy file error:\n{e}', ha='center', va='center', transform=ax.transAxes)

# Plot 6: Radial velocity vs radius
ax = axes[1, 2]
# Radial velocity: v_r = (r . v) / |r|
vr0 = (snap0['pos_x']*snap0['vel_x'] + snap0['pos_y']*snap0['vel_y'] + snap0['pos_z']*snap0['vel_z']) / np.maximum(r0, 1e-10)
vrf = (snap_final['pos_x']*snap_final['vel_x'] + snap_final['pos_y']*snap_final['vel_y'] + snap_final['pos_z']*snap_final['vel_z']) / np.maximum(rf, 1e-10)

ax.scatter(r0, vr0, s=1, alpha=0.3, label='Initial')
ax.scatter(rf, vrf, s=1, alpha=0.3, label='Final')
ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
ax.set_xlabel('Radius [pc]')
ax.set_ylabel('Radial Velocity [km/s]')
ax.set_title('Radial Velocity vs Radius')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{results_dir}/relaxation_analysis.png', dpi=150)
print(f"Saved: {results_dir}/relaxation_analysis.png")
plt.close()

# Print statistics
print("\n=== Statistics ===")
print(f"Initial: R_max = {r0.max():.3f} pc, <rho> = {snap0['dens'].mean():.2f} M_sun/pc^3")
print(f"Final:   R_max = {rf.max():.3f} pc, <rho> = {snap_final['dens'].mean():.2f} M_sun/pc^3")
print(f"Initial: <v> = {v0.mean():.4f} km/s, v_max = {v0.max():.4f} km/s")
print(f"Final:   <v> = {vf.mean():.4f} km/s, v_max = {vf.max():.4f} km/s")
