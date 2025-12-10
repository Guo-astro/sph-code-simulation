#!/usr/bin/env python3
"""
Analyze the initial condition snapshot to extract cloud properties.
"""
import pandas as pd
import numpy as np

snapshot_path = "sample/imbh_cloud/results/hydrostatic/61k/GSPH_kernel_gravity/snapshot_0000.csv"

# Read the snapshot, skip comment lines
df = pd.read_csv(snapshot_path, comment='#')

print("=" * 60)
print("CLOUD PROPERTIES FROM HYDROSTATIC SNAPSHOT")
print("=" * 60)

print(f"\n=== Particle Statistics ===")
print(f"Number of particles: {len(df)}")
print(f"Particle mass (each): {df['mass'].iloc[0]:.6e}")
print(f"Total cloud mass (M_cloud): {df['mass'].sum():.6f} code units")

# Find central density (max density)
rho_max = df['dens'].max()
rho_mean = df['dens'].mean()
rho_min = df['dens'].min()
print(f"\n=== Density ===")
print(f"Central density (rho_c): {rho_max:.6f} code units")
print(f"Mean density: {rho_mean:.6f} code units")
print(f"Min density (edge): {rho_min:.6f} code units")
print(f"Density contrast (rho_c/rho_edge): {rho_max/rho_min:.2f}")

# Calculate cloud radius
r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
R_cloud = r.max()
R_90 = np.percentile(r, 90)
R_50 = np.percentile(r, 50)
print(f"\n=== Cloud Radius ===")
print(f"Maximum radius (R_cloud): {R_cloud:.6f} code units")
print(f"90th percentile radius: {R_90:.6f} code units")
print(f"Half-mass radius (~): {R_50:.6f} code units")

# Sound speed analysis
idx_max_dens = df['dens'].idxmax()
c_s_mean = df['sound'].mean()
c_s_center = df.loc[idx_max_dens, 'sound']
print(f"\n=== Sound Speed ===")
print(f"Central sound speed (c_s,c): {c_s_center:.6f} code units")
print(f"Mean sound speed: {c_s_mean:.6f} code units")

# Pressure analysis
P_center = df.loc[idx_max_dens, 'pres']
P_mean = df['pres'].mean()
print(f"\n=== Pressure ===")
print(f"Central pressure (P_c): {P_center:.6f} code units")
print(f"Mean pressure: {P_mean:.6f} code units")

# Internal energy
e_center = df.loc[idx_max_dens, 'ene']
print(f"\n=== Internal Energy ===")
print(f"Central specific energy: {e_center:.6f} code units")

# Verify polytropic relation: P = K * rho^gamma
gamma = 5.0/3.0
K_measured = P_center / (rho_max ** gamma)
print(f"\n=== Polytrope Verification ===")
print(f"gamma = {gamma:.6f}")
print(f"K = P_c / rho_c^gamma = {K_measured:.6f}")

# Timescales
G = 1.0
t_ff = np.sqrt(3 * np.pi / (32 * G * rho_max))
t_sound = R_cloud / c_s_center
t_dyn = np.sqrt(R_cloud**3 / (G * df['mass'].sum()))
print(f"\n=== Timescales (G=1) ===")
print(f"Free-fall time (central): t_ff = {t_ff:.6f} code units")
print(f"Sound crossing time: t_sound = {t_sound:.6f} code units")
print(f"Dynamical time: t_dyn = {t_dyn:.6f} code units")

# Lane-Emden n=3/2 (gamma=5/3) theoretical values
# For n=3/2: xi_1 = 3.65375, -xi^2 * theta'(xi_1) = 2.71406
xi_1 = 3.65375
omega_n = 2.71406
print(f"\n=== Lane-Emden n=3/2 Theory ===")
print(f"xi_1 (first zero): {xi_1}")
print(f"omega_n = -xi^2 * theta'(xi): {omega_n}")

# Calculate alpha (length scale)
# R = alpha * xi_1, so alpha = R / xi_1
alpha = R_cloud / xi_1
print(f"alpha = R/xi_1 = {alpha:.6f}")

# Verify: alpha^2 = K*(n+1)*rho_c^(1/n - 1) / (4*pi*G)
# For n=3/2: alpha^2 = K*2.5*rho_c^(-1/3) / (4*pi*G)
n = 1.5
alpha_sq_theory = K_measured * (n + 1) * (rho_max ** (1.0/n - 1)) / (4 * np.pi * G)
print(f"alpha^2 (theory) = {alpha_sq_theory:.6f}")
print(f"alpha (theory) = {np.sqrt(alpha_sq_theory):.6f}")
print(f"alpha (measured) = {alpha:.6f}")

# Total mass from Lane-Emden
# M = 4*pi*alpha^3*rho_c*omega_n
M_theory = 4 * np.pi * alpha**3 * rho_max * omega_n
print(f"\n=== Mass Verification ===")
print(f"M (from particles): {df['mass'].sum():.6f}")
print(f"M (Lane-Emden theory): {M_theory:.6f}")

print("\n" + "=" * 60)
print("SUMMARY FOR SIMULATION SETUP")
print("=" * 60)
print(f"""
Cloud Parameters (code units, G=1):
  - Total mass:       M_cloud = {df['mass'].sum():.4f}
  - Radius:           R_cloud = {R_cloud:.4f}
  - Central density:  rho_c   = {rho_max:.4f}
  - Central pressure: P_c     = {P_center:.4f}
  - Central c_s:      c_s     = {c_s_center:.4f}
  - Polytropic K:     K       = {K_measured:.4f}
  - gamma:            γ       = {gamma:.4f}
  
Timescales:
  - Free-fall:        t_ff    = {t_ff:.4f}
  - Sound crossing:   t_sound = {t_sound:.4f}
  - Dynamical:        t_dyn   = {t_dyn:.4f}
""")
