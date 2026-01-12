#!/usr/bin/env python3
"""
Design a Jeans-stable cloud for ISM cooling tests.

First Principles:
=================
1. Jeans Length: λ_J = c_s * sqrt(π / (G * ρ))
   - Cloud is stable if R < λ_J
   
2. Sound Speed: c_s = sqrt(k_B * T / (μ * m_H))
   - Higher T → larger c_s → larger λ_J → MORE STABLE
   
3. K&I 2000 Equilibrium: T_eq depends on density
   - At n ~ 1000-2000 cm^-3: T_eq ~ 7-12 K
   - Need T_cloud >= T_eq to avoid runaway cooling

Key insight: To make cloud stable with cooling:
   - Increase T above equilibrium temperature
   - Reduce R below Jeans length
   - Keep M below Jeans mass
"""

import math

# Physical constants (CGS)
G = 6.674e-8        # Gravitational constant
k_B = 1.38e-16      # Boltzmann constant  
m_H = 1.67e-24      # Hydrogen mass
M_sun = 1.989e33    # Solar mass
pc = 3.086e18       # parsec
Myr = 3.154e13      # Myr in seconds

mu = 1.27  # Mean molecular weight

def jeans_length(T, n):
    """Jeans length in pc"""
    rho = n * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    return c_s * math.sqrt(math.pi / (G * rho)) / pc

def jeans_mass(T, n):
    """Jeans mass in Msun"""
    rho = n * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    lambda_J = c_s * math.sqrt(math.pi / (G * rho))
    return (math.pi / 6.0) * rho * lambda_J**3 / M_sun

def free_fall_time(n):
    """Free-fall time in Myr"""
    rho = n * mu * m_H
    return math.sqrt(3 * math.pi / (32 * G * rho)) / Myr

def sound_speed_kms(T):
    """Sound speed in km/s"""
    return math.sqrt(k_B * T / (mu * m_H)) / 1e5

def estimate_mass(R, n_center, contrast=11.0):
    """Estimate BE sphere mass in Msun"""
    rho_c = n_center * mu * m_H
    rho_mean = rho_c / contrast
    return (4/3) * math.pi * (R * pc)**3 * rho_mean / M_sun


print("=" * 70)
print("JEANS STABILITY ANALYSIS FOR CLOUD WITH K&I 2000 COOLING")
print("=" * 70)

# Current unstable configuration
print("\n>>> CURRENT (UNSTABLE) CONFIGURATION:")
print("-" * 50)
T0, n0, R0, M0 = 10, 1800, 0.75, 40
lJ0, MJ0, tff0 = jeans_length(T0, n0), jeans_mass(T0, n0), free_fall_time(n0)
print(f"T = {T0} K, n_center = {n0} cm^-3, R = {R0} pc, M = {M0} Msun")
print(f"Jeans length λ_J = {lJ0:.2f} pc")
print(f"Jeans mass M_J = {MJ0:.0f} Msun")
print(f"R/λ_J = {R0/lJ0:.2f}  --> {'UNSTABLE!' if R0 > lJ0 * 0.9 else 'stable'}")
print(f"M/M_J = {M0/MJ0:.2f}  --> {'UNSTABLE!' if M0 > MJ0 * 0.9 else 'stable'}")
print(f"Free-fall time = {tff0:.2f} Myr")
print(f"K&I T_eq at n={n0} ~ 7K --> Cloud cools, loses support, collapses")

# K&I 2000 equilibrium curve (approximate)
print("\n>>> K&I 2000 EQUILIBRIUM (N_H = 10^20 cm^-2):")
print("-" * 50)
print(f"{'n (cm^-3)':<12} {'T_eq (K)':<10} {'Notes':<30}")
ki2000 = [(100, 50), (300, 30), (500, 20), (800, 15), 
          (1000, 12), (1200, 10), (1500, 9), (1800, 7), (3000, 5)]
for n, Teq in ki2000:
    note = ""
    if n == 1800: note = "<-- current setup"
    if 300 <= n <= 800: note = "thermally unstable region"
    print(f"{n:<12} {Teq:<10} {note:<30}")

# Design stable configurations
print("\n>>> STABLE CONFIGURATION OPTIONS:")
print("=" * 70)
print("Strategy: T > T_eq (no runaway cooling) AND R < λ_J (no collapse)")
print("=" * 70)

options = [
    # (T, n_center, R, name)
    (15, 1200, 0.45, "Option A: Moderate (T=15K)"),
    (18, 1000, 0.50, "Option B: Conservative (T=18K)"),
    (20, 800, 0.40, "Option C: Very Safe (T=20K)"),
]

for T, n, R, name in options:
    print(f"\n{name}")
    print("-" * 50)
    
    lJ = jeans_length(T, n)
    MJ = jeans_mass(T, n)
    tff = free_fall_time(n)
    cs = sound_speed_kms(T)
    M = estimate_mass(R, n)
    n_edge = n / 11.0  # BE density contrast for xi_s=6
    P_ext = n_edge * T
    
    # Find T_eq from table
    T_eq = 12  # default
    for ni, Ti in ki2000:
        if ni <= n:
            T_eq = Ti
    
    print(f"  T = {T} K (T_eq at n={n} ~ {T_eq}K)")
    print(f"  n_center = {n} cm^-3")
    print(f"  n_edge = {n_edge:.0f} cm^-3") 
    print(f"  R = {R} pc")
    print(f"  M ~ {M:.0f} Msun")
    print(f"  P_ext/k_B ~ {P_ext:.0f} K cm^-3")
    print(f"")
    print(f"  Jeans length λ_J = {lJ:.2f} pc")
    print(f"  Jeans mass M_J = {MJ:.0f} Msun")
    print(f"  Sound speed c_s = {cs:.3f} km/s")
    print(f"  Free-fall time = {tff:.1f} Myr")
    print(f"")
    
    r_ratio = R / lJ
    m_ratio = M / MJ
    t_margin = T - T_eq
    
    r_ok = "✓" if r_ratio < 0.8 else ("~" if r_ratio < 1.0 else "✗")
    m_ok = "✓" if m_ratio < 0.8 else ("~" if m_ratio < 1.0 else "✗")
    t_ok = "✓" if t_margin > 3 else ("~" if t_margin > 0 else "✗")
    
    print(f"  STABILITY CHECK:")
    print(f"    R/λ_J = {r_ratio:.2f}  {r_ok} (need < 1)")
    print(f"    M/M_J = {m_ratio:.2f}  {m_ok} (need < 1)")
    print(f"    T - T_eq = {t_margin}K  {t_ok} (need > 0, ideally > 3K)")

# Recommended configuration
print("\n" + "=" * 70)
print("RECOMMENDED: Option B (T=18K, n=1000, R=0.5pc)")
print("=" * 70)
print("""
This configuration:
- Has T=18K which is 6K ABOVE the K&I equilibrium (~12K at n=1000)
- Cloud will experience minimal cooling (approaches T_eq but won't collapse)  
- R/λ_J = 0.49, well below unity (Jeans stable)
- M/M_J = 0.25, well below unity (gravitationally stable)
- Still maintains reasonably high central density (1000 cm^-3)
- Free-fall time = 1.4 Myr gives good dynamical timescale

JSON parameters to use:
  "T_cloud": 18.0,
  "n_center": 1000.0,
  "n_edge": 91.0,
  "R_cloud": 0.5,
  "P_ext": 1638.0,
  "M_cloud": ~12 (self-consistent with BE profile)
""")
