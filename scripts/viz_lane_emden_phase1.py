#!/usr/bin/env python3
"""Lane-Emden Phase 1 relaxation GIF with analytic density overlay"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
from scipy.integrate import odeint
import h5py

# Lane-Emden solution for n=1.5
def lane_emden_ode(y, xi):
    theta, dtheta = y
    if xi < 1e-10:
        return [dtheta, -1/3]
    if theta < 0:
        return [0, 0]
    return [dtheta, -theta**1.5 - 2*dtheta/xi]

xi_span = np.linspace(1e-10, 5, 1000)
y0 = [1.0, 0.0]
sol = odeint(lane_emden_ode, y0, xi_span)
theta_sol = np.maximum(sol[:, 0], 0)
dtheta_sol = sol[:, 1]

idx_zero = np.where(theta_sol <= 0)[0]
if len(idx_zero) > 0:
    xi_1 = xi_span[idx_zero[0]]
    dtheta_1 = dtheta_sol[idx_zero[0]]
else:
    xi_1 = 3.65375
    dtheta_1 = -0.20330

# Physical parameters (code units: L=pc, M=M_sun, V=km/s)
R_cloud = 1.0  # pc
M_cloud = 40.0  # M_sun (from config)
n_center_target = 1000.0  # cm^-3

# Unit conversion: n [cm^-3] = rho [M_sun/pc^3] * density_to_n
# density_to_n = M_sun / (mu * m_H) / pc^3
#              = 1.989e33 / (1.27 * 1.67e-24) / (3.086e18)^3
#              = 20.3 cm^-3 per (M_sun/pc^3)
density_to_n = 20.3  # for mu=1.27 (atomic gas)

rho_center = n_center_target / density_to_n  # M_sun/pc^3 (~49.3)
alpha = R_cloud / xi_1

# Verify M_cloud from Lane-Emden solution
M_cloud_calc = 4 * np.pi * alpha**3 * rho_center * xi_1**2 * abs(dtheta_1)
print(f"Lane-Emden: M={M_cloud:.2f} M_sun (config), M_calc={M_cloud_calc:.2f} M_sun")
print(f"Lane-Emden: rho_c={rho_center:.1f} M_sun/pc^3 = n_c={n_center_target:.0f} cm^-3")

r_analytic = xi_span * alpha
rho_analytic = rho_center * theta_sol**1.5
mask = r_analytic <= R_cloud
r_analytic = r_analytic[mask]
rho_analytic = rho_analytic[mask]

# Load phase1 snapshots
data_dir = Path("simulations/astrophysics/imbh_cloud/results/lane_emden/phase1_relaxation")
snapshots = sorted(data_dir.glob("snapshot_*.h5"))
if not snapshots:
    snapshots = sorted(data_dir.glob("snapshot_*.csv"))
use_hdf5 = len(snapshots) > 0 and snapshots[0].suffix == '.h5'
print(f"Found {len(snapshots)} snapshots ({'HDF5' if use_hdf5 else 'CSV'})")

# Pre-scan to get axis limits (excluding outliers for better visualization)
print("Scanning snapshots for axis limits...")
all_r = []
all_rho = []
for snap in snapshots:
    if use_hdf5:
        with h5py.File(snap, 'r') as f:
            x = np.array(f['particles/pos_x'])
            y = np.array(f['particles/pos_y'])
            z = np.array(f['particles/pos_z'])
            rho_code = np.array(f['particles/dens'])
    else:
        df = pd.read_csv(snap, comment='#')
        x, y, z = df['pos_x'].values, df['pos_y'].values, df['pos_z'].values
        rho_code = df['dens'].values
    r = np.sqrt(x**2 + y**2 + z**2)
    # Density is now in physical units (M_sun/pc^3) directly from simulation
    # No mass_scale multiplication needed after lane_emden.cpp fix
    rho = rho_code
    all_r.extend(r)
    all_rho.extend(rho)

all_r = np.array(all_r)
all_rho = np.array(all_rho)

# Filter out NaN values
valid_mask = ~np.isnan(all_r) & ~np.isnan(all_rho)
all_r = all_r[valid_mask]
all_rho = all_rho[valid_mask]

# Use 99th percentile to exclude outliers
global_r_max = np.nanpercentile(all_r, 99) * 1.2
global_rho_max = np.nanpercentile(all_rho, 99) * 1.1
print(f"Axis limits (99th percentile): r_max={global_r_max:.3f} pc, rho_max={global_rho_max:.1f} M☉/pc³")
print(f"True max: r_max={all_r.max():.3f} pc ({(all_r > global_r_max).sum()} outliers)")

fig, ax = plt.subplots(figsize=(10, 7))

def update(frame):
    ax.clear()
    if use_hdf5:
        with h5py.File(snapshots[frame], 'r') as f:
            x = np.array(f['particles/pos_x'])
            y = np.array(f['particles/pos_y'])
            z = np.array(f['particles/pos_z'])
            rho_code = np.array(f['particles/dens'])
    else:
        df = pd.read_csv(snapshots[frame], comment='#')
        x, y, z = df['pos_x'].values, df['pos_y'].values, df['pos_z'].values
        rho_code = df['dens'].values

    # Density is now in physical units (M_sun/pc^3) directly from simulation
    rho = rho_code
    r = np.sqrt(x**2 + y**2 + z**2)

    ax.scatter(r, rho, s=1, alpha=0.3, c='blue', label='SPH particles')
    ax.plot(r_analytic, rho_analytic, 'r-', lw=2, label=f'Lane-Emden (M={M_cloud:.1f}M☉)')

    # Relaxation step
    step = frame * 10  # relaxationOutputFreq = 10
    ax.set_xlabel('Radius [pc]', fontsize=12)
    ax.set_ylabel('Density [M☉/pc³]', fontsize=12)
    ax.set_title(f'Lane-Emden γ=5/3 Relaxation - Step {step}', fontsize=14)
    ax.set_xlim(0, global_r_max)
    ax.set_ylim(0, global_rho_max)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Secondary y-axis for number density (cm^-3)
    ax2 = ax.secondary_yaxis('right', functions=(lambda x: x * density_to_n, lambda x: x / density_to_n))
    ax2.set_ylabel('Number density [cm⁻³]', fontsize=12)

    r_max = r.max()
    rho_max = rho.max()
    n_max = rho_max * density_to_n
    n_center_actual = rho.max() * density_to_n if len(rho) > 0 else 0
    n_outliers = (r > global_r_max).sum()
    ax.text(0.02, 0.98, f'R_max = {r_max:.3f} pc\nρ_max = {rho_max:.1f} M☉/pc³\nn_max = {n_max:.0f} cm⁻³\nTarget: {n_center_target:.0f} cm⁻³',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    return ax,

ani = FuncAnimation(fig, update, frames=len(snapshots), interval=200, blit=False)
output_path = data_dir / "phase1_evolution_density.gif"
ani.save(output_path, writer=PillowWriter(fps=5))
print(f"Saved: {output_path}")
plt.close()
