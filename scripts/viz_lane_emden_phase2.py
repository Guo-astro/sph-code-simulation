#!/usr/bin/env python3
"""Lane-Emden Phase 2 evolution GIF with analytic density overlay"""
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

# =============================================
# PHYSICAL PARAMETERS - Derived from Lane-Emden
# =============================================
R_cloud = 1.0  # pc
n_center_target = 1000.0  # cm^-3
density_to_n = 17.37  # cm^-3 per (M_sun/pc^3)
rho_center = n_center_target / density_to_n  # M_sun/pc^3
alpha = R_cloud / xi_1

# Mass from Lane-Emden
M_cloud = 4 * np.pi * alpha**3 * rho_center * xi_1**2 * abs(dtheta_1)
mass_scale = M_cloud  # Code uses M=1, R=1

print(f"Lane-Emden: M={M_cloud:.2f} M_sun, rho_c={rho_center:.1f} M_sun/pc^3")

r_analytic = xi_span * alpha
rho_analytic = rho_center * theta_sol**1.5
mask = r_analytic <= R_cloud
r_analytic = r_analytic[mask]
rho_analytic = rho_analytic[mask]

# Load snapshots (prefer HDF5, fall back to CSV)
data_dir = Path("simulations/astrophysics/imbh_cloud/results/lane_emden/phase2_hydrostatic")
snapshots = sorted(data_dir.glob("snapshot_*.h5"))
if not snapshots:
    snapshots = sorted(data_dir.glob("snapshot_*.csv"))
use_hdf5 = len(snapshots) > 0 and snapshots[0].suffix == '.h5'
print(f"Found {len(snapshots)} snapshots ({'HDF5' if use_hdf5 else 'CSV'})")

# Pre-scan all snapshots to get global axis limits
print("Scanning snapshots for axis limits...")
global_r_max = 0
global_rho_max = 0
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
    rho = rho_code * mass_scale
    global_r_max = max(global_r_max, r.max())
    global_rho_max = max(global_rho_max, rho.max())

# Add 10% padding
global_r_max *= 1.1
global_rho_max *= 1.1
print(f"Axis limits: r_max={global_r_max:.3f} pc, rho_max={global_rho_max:.1f} M☉/pc³")

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

    # Convert code units to physical
    rho = rho_code * mass_scale
    r = np.sqrt(x**2 + y**2 + z**2)

    ax.scatter(r, rho, s=1, alpha=0.3, c='blue', label='SPH particles')
    ax.plot(r_analytic, rho_analytic, 'r-', lw=2, label=f'Lane-Emden (M={M_cloud:.1f}M☉)')

    # Time from snapshot number (dt ~ 0.5 code units)
    time = frame * 0.5
    ax.set_xlabel('Radius [pc]', fontsize=12)
    ax.set_ylabel('Density [M☉/pc³] / [cm⁻³]', fontsize=12)
    ax.set_title(f'Lane-Emden γ=5/3 with Self-Gravity - t = {time:.1f} [code] ~ {time*0.978:.1f} Myr', fontsize=14)
    ax.set_xlim(0, global_r_max)
    ax.set_ylim(0, global_rho_max)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    r_max = r.max()
    rho_max = rho.max()
    n_max = rho_max * density_to_n
    ax.text(0.02, 0.98, f'R_max = {r_max:.3f} pc\nρ_max = {rho_max:.1f} M☉/pc³\nn_max = {n_max:.0f} cm⁻³',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    return ax,

ani = FuncAnimation(fig, update, frames=len(snapshots), interval=200, blit=False)
output_path = data_dir / "phase2_evolution_density.gif"
ani.save(output_path, writer=PillowWriter(fps=5))
print(f"Saved: {output_path}")
plt.close()
