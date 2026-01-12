#!/usr/bin/env python3
"""Check Jeans mass stability during cooling"""
import math

G = 6.674e-8
k_B = 1.38e-16
m_H = 1.67e-24
pc = 3.086e18
M_sun = 1.989e33
mu = 1.27

n_center = 3000.0
T_init = 20.0
T_eq = 6.0
M_cloud = 22.5  # Msun

rho_c = n_center * mu * m_H

def jeans_mass(T, rho):
    c_s = math.sqrt(k_B * T / (mu * m_H))
    lJ = c_s * math.sqrt(math.pi / (G * rho))
    M_J = (math.pi/6) * rho * lJ**3 / M_sun
    return M_J, lJ/pc

MJ_init, lJ_init = jeans_mass(T_init, rho_c)
MJ_eq, lJ_eq = jeans_mass(T_eq, rho_c)

print("=== INITIAL STATE (T=20K, n=3000) ===")
print(f"M_cloud = {M_cloud} Msun")
print(f"M_Jeans = {MJ_init:.1f} Msun")
print(f"M/M_J = {M_cloud/MJ_init:.2f}")
print(f"lambda_J = {lJ_init:.3f} pc")

print("\n=== AT T=6K (same density) ===")
print(f"M_Jeans = {MJ_eq:.1f} Msun")
print(f"M/M_J = {M_cloud/MJ_eq:.2f}")
print(f"lambda_J = {lJ_eq:.3f} pc")

if M_cloud/MJ_eq > 1:
    print("\n*** MASS EXCEEDS JEANS MASS! ***")
    print("The cloud WILL collapse when cooled to equilibrium!")
    print("This is the TRUE ROOT CAUSE of the simulation failure!")
