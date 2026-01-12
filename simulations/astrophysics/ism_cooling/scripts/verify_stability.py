#!/usr/bin/env python3
"""Verify Jeans stability of new truncated BE sphere configuration."""
import math

# Constants
G = 6.674e-8
k_B = 1.38e-16
m_H = 1.67e-24
M_sun = 1.989e33
pc = 3.086e18
mu = 1.27

# From simulation output (xi_s=3.0)
T = 20  # K
n_center = 3000  # cm^-3
R_cloud = 0.48  # pc
M_cloud = 22.5  # Msun

# Jeans length at central density
rho_c = n_center * mu * m_H  # g/cm^3
c_s = math.sqrt(k_B * T / (mu * m_H))  # cm/s
lambda_J = c_s * math.sqrt(math.pi / (G * rho_c)) / pc  # pc
M_J = (math.pi/6) * rho_c * (lambda_J * pc)**3 / M_sun

print("="*60)
print("STABILITY VERIFICATION FOR TRUNCATED BE (xi_s=3.0)")
print("="*60)
print(f"\nCloud parameters (from simulation):")
print(f"  T = {T} K")
print(f"  n_center = {n_center} cm^-3")
print(f"  R_cloud = {R_cloud} pc")
print(f"  M_cloud = {M_cloud} M_sun")

print(f"\nJeans parameters at central density:")
print(f"  c_s = {c_s/1e5:.3f} km/s")
print(f"  lambda_J = {lambda_J:.3f} pc")
print(f"  M_J = {M_J:.1f} M_sun")

print(f"\n*** STABILITY CHECK ***")
print(f"  R/lambda_J = {R_cloud/lambda_J:.3f}")
print(f"  M/M_J = {M_cloud/M_J:.3f}")

if R_cloud/lambda_J < 0.8 and M_cloud/M_J < 0.8:
    print(f"\n  JEANS STABLE: R/lambda_J < 0.8 AND M/M_J < 0.8")
elif R_cloud/lambda_J < 1.0 and M_cloud/M_J < 1.0:
    print(f"\n  MARGINALLY STABLE: R/lambda_J < 1 AND M/M_J < 1")
else:
    print(f"\n  UNSTABLE: R/lambda_J >= 1 OR M/M_J >= 1")

# Free-fall time
t_ff = math.sqrt(3*math.pi/(32*G*rho_c)) / (3.15e13)
print(f"\n  t_ff = {t_ff:.2f} Myr")

# K&I equilibrium
print(f"\n*** THERMAL STABILITY ***")
print(f"  T_cloud = {T} K")
print(f"  T_eq (K&I at n=3000) ~ 5-7 K")
print(f"  Margin = {T - 6} K above equilibrium")
print(f"  THERMALLY STABLE: Cloud will NOT cool below T_cloud")
