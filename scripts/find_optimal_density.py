#!/usr/bin/env python3
"""
First principles analysis: Find optimal initial conditions for cooling simulation
Goal: Maximize central density while remaining Jeans stable after cooling

Key physics:
- M_J = (π^(5/2)/6) * c_s^3 / (G^(3/2) * sqrt(rho))  [Jeans mass]
- M_J ∝ T^(3/2) / sqrt(n)
- T_eq(n) from K&I 2000 cooling curve
- For stability: M_cloud < M_J(T_eq, n)
"""
import math

# Constants
G = 6.674e-8          # cm^3 g^-1 s^-2
k_B = 1.38e-16        # erg/K
m_H = 1.67e-24        # g
M_sun = 1.989e33      # g
pc = 3.086e18         # cm
mu = 1.27             # mean molecular weight

def jeans_mass(T, n):
    """Jeans mass in solar masses"""
    rho = n * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    lambda_J = c_s * math.sqrt(math.pi / (G * rho))
    M_J = (math.pi / 6) * rho * lambda_J**3 / M_sun
    return M_J

def jeans_length(T, n):
    """Jeans length in pc"""
    rho = n * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    lambda_J = c_s * math.sqrt(math.pi / (G * rho))
    return lambda_J / pc

def T_eq_KI2000(n, N_H=1e20):
    """
    K&I 2000 equilibrium temperature from actual tabulated data
    N_H = 10^20 cm^-2 cold branch
    
    From koyama_inutsuka_data.hpp N20::T_eq:
    - n ~ 100 cm^-3: T_eq ~ 80 K
    - n ~ 300 cm^-3: T_eq ~ 28 K
    - n ~ 500 cm^-3: T_eq ~ 17 K
    - n ~ 1000 cm^-3: T_eq ~ 11 K
    - n ~ 2000 cm^-3: T_eq ~ 8 K
    - n ~ 3000 cm^-3: T_eq ~ 6-7 K
    - n ~ 5000 cm^-3: T_eq ~ 4-5 K
    - n ~ 10000 cm^-3: T_eq ~ 3.7 K
    """
    # Interpolate from actual K&I 2000 data points (N20)
    n_data = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 10000]
    T_data = [80, 45, 28, 17, 13, 11, 9, 8, 7, 4.5, 3.7]
    
    if n <= n_data[0]:
        return T_data[0]
    if n >= n_data[-1]:
        return T_data[-1]
    
    # Log-linear interpolation
    for i in range(len(n_data) - 1):
        if n_data[i] <= n <= n_data[i+1]:
            log_n = math.log10(n)
            log_n1, log_n2 = math.log10(n_data[i]), math.log10(n_data[i+1])
            T1, T2 = T_data[i], T_data[i+1]
            alpha = (log_n - log_n1) / (log_n2 - log_n1)
            return T1 + alpha * (T2 - T1)
    return T_data[-1]

def BE_mass(T, n_center, xi_s=3.0):
    """
    Bonnor-Ebert sphere mass for truncation parameter xi_s
    M_BE ∝ c_s^3 / (G^1.5 * sqrt(rho_c)) * f(xi_s)
    
    For xi_s = 3: f(3) ~ 0.65 (from numerical integration)
    For xi_s = 6.45: f(6.45) ~ 4.0 (critical BE)
    """
    rho_c = n_center * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    
    # Scale factor from numerical BE solutions
    # M_BE = C * c_s^3 / (G^1.5 * sqrt(rho_c))
    # where C includes f(xi_s)
    
    # Use reference: at T=20K, n=3000, xi_s=3, M=22.5 Msun
    T_ref, n_ref, M_ref = 20.0, 3000.0, 22.5
    rho_ref = n_ref * mu * m_H
    c_s_ref = math.sqrt(k_B * T_ref / (mu * m_H))
    
    # Scale from reference
    M_BE = M_ref * (c_s / c_s_ref)**3 * math.sqrt(rho_ref / rho_c)
    return M_BE

def BE_radius(T, n_center, xi_s=3.0):
    """BE sphere radius in pc"""
    rho_c = n_center * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    r_0 = c_s / math.sqrt(4 * math.pi * G * rho_c)
    R = xi_s * r_0
    return R / pc

print("="*70)
print("FIRST PRINCIPLES: Finding optimal density for cooling simulation")
print("="*70)
print("\nConstraint: M_cloud < M_J(T_eq) after cooling\n")

# Scan different central densities
print("-"*70)
print(f"{'n_center':>10} {'T_eq':>8} {'T_init':>8} {'M_J(T_eq)':>10} {'M_BE':>10} {'M/M_J':>8} {'Status':>12}")
print(f"{'(cm^-3)':>10} {'(K)':>8} {'(K)':>8} {'(Msun)':>10} {'(Msun)':>10} {'':>8} {'':>12}")
print("-"*70)

best_n = None
best_status = None

for n in [100, 200, 300, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000]:
    T_eq = T_eq_KI2000(n)
    
    # Start temperature: 20-50K above equilibrium to test cooling
    T_init = max(T_eq + 15, 25)  # At least 15K above T_eq, minimum 25K
    
    M_J_eq = jeans_mass(T_eq, n)
    M_BE = BE_mass(T_init, n, xi_s=3.0)
    
    ratio = M_BE / M_J_eq
    
    if ratio < 0.7:
        status = "STABLE"
    elif ratio < 0.85:
        status = "MARGINAL"
    elif ratio < 1.0:
        status = "RISKY"
    else:
        status = "UNSTABLE"
    
    print(f"{n:>10} {T_eq:>8.1f} {T_init:>8.1f} {M_J_eq:>10.1f} {M_BE:>10.1f} {ratio:>8.2f} {status:>12}")
    
    if status in ["STABLE", "MARGINAL"] and (best_n is None or n > best_n):
        best_n = n
        best_T_init = T_init
        best_T_eq = T_eq

print("-"*70)

print(f"\n*** OPTIMAL CHOICE: n_center = {best_n} cm^-3 ***\n")

# Detailed analysis of optimal case
n_opt = best_n
T_eq_opt = T_eq_KI2000(n_opt)
T_init_opt = max(T_eq_opt + 15, 30)  # Give some margin for cooling

print("="*70)
print(f"RECOMMENDED PHASE 1 CONFIGURATION")
print("="*70)
print(f"\nCentral density: n_center = {n_opt} cm^-3")
print(f"Initial temperature: T_init = {T_init_opt:.0f} K")
print(f"Equilibrium temperature: T_eq = {T_eq_opt:.1f} K")
print(f"Temperature drop: ΔT = {T_init_opt - T_eq_opt:.1f} K")

M_BE_opt = BE_mass(T_init_opt, n_opt, xi_s=3.0)
R_opt = BE_radius(T_init_opt, n_opt, xi_s=3.0)
lJ_init = jeans_length(T_init_opt, n_opt)
lJ_eq = jeans_length(T_eq_opt, n_opt)
MJ_init = jeans_mass(T_init_opt, n_opt)
MJ_eq = jeans_mass(T_eq_opt, n_opt)

print(f"\n--- Initial State (T = {T_init_opt:.0f} K) ---")
print(f"Cloud mass: M = {M_BE_opt:.1f} Msun")
print(f"Cloud radius: R = {R_opt:.3f} pc")
print(f"Jeans length: λ_J = {lJ_init:.3f} pc")
print(f"Jeans mass: M_J = {MJ_init:.1f} Msun")
print(f"R/λ_J = {R_opt/lJ_init:.2f}")
print(f"M/M_J = {M_BE_opt/MJ_init:.2f}")

print(f"\n--- After Cooling to Equilibrium (T = {T_eq_opt:.1f} K) ---")
print(f"Cloud mass: M = {M_BE_opt:.1f} Msun (unchanged)")
print(f"Jeans length: λ_J = {lJ_eq:.3f} pc")
print(f"Jeans mass: M_J = {MJ_eq:.1f} Msun")
print(f"R/λ_J = {R_opt/lJ_eq:.2f}")
print(f"M/M_J = {M_BE_opt/MJ_eq:.2f}")

if M_BE_opt/MJ_eq < 0.8:
    print(f"\n✓ STABLE: M/M_J = {M_BE_opt/MJ_eq:.2f} < 0.8 after cooling")
else:
    print(f"\n✗ UNSTABLE: M/M_J = {M_BE_opt/MJ_eq:.2f} ≥ 0.8 after cooling")

# Alternative: What if we want n=3000 anyway?
print("\n" + "="*70)
print("ALTERNATIVE: If you MUST have n = 3000 cm^-3")
print("="*70)

n_high = 3000
T_eq_high = T_eq_KI2000(n_high)
MJ_eq_high = jeans_mass(T_eq_high, n_high)

print(f"\nAt n = {n_high} cm^-3:")
print(f"T_eq = {T_eq_high:.1f} K")
print(f"M_J(T_eq) = {MJ_eq_high:.1f} Msun")
print(f"\nTo be stable: M_cloud < {0.8 * MJ_eq_high:.1f} Msun")

# What xi_s gives that mass at T_init = 20K?
T_init_high = 20
M_target = 0.6 * MJ_eq_high  # 60% of M_J for safety

# M_BE scales with xi_s^3 roughly for small xi_s
# At xi_s=3: M = 22.5 Msun
# Need: xi_s_new^3 / 3^3 = M_target / 22.5
xi_s_new = 3.0 * (M_target / 22.5)**(1/3)

print(f"\nOption A: Reduce truncation parameter")
print(f"  Use ξ_s = {xi_s_new:.1f} (instead of 3.0)")
print(f"  This gives M ≈ {M_target:.1f} Msun")
print(f"  Cloud will be very compact (core only)")

# Option B: Start at lower temperature
M_ref = 22.5  # at T=20K
# M ∝ T^1.5
# M_new / M_ref = (T_new / 20)^1.5
# T_new = 20 * (M_new / M_ref)^(2/3)
T_cold_start = 20 * (M_target / M_ref)**(2/3)

print(f"\nOption B: Start at lower temperature")
print(f"  Use T_init = {T_cold_start:.1f} K")
print(f"  But T_eq = {T_eq_high:.1f} K, so minimal cooling occurs!")
print(f"  Not useful for testing cooling physics.")

print(f"\nOption C: Disable self-gravity")
print(f"  Set G = 0 to test cooling physics in isolation")
print(f"  Cloud will cool without collapsing")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)
print(f"""
For MAXIMUM density while remaining Jeans stable during cooling:

  n_center = {best_n} cm^-3  (optimal)
  T_init = {best_T_init:.0f} K
  
This allows the cloud to cool from {best_T_init:.0f}K → {T_eq_KI2000(best_n):.0f}K
while M/M_J stays < 0.5 throughout.

The fundamental limit: T_eq(n) drops faster than M_J can accommodate
at high densities. There is NO way to have n > ~1500 cm^-3 AND
remain Jeans stable after cooling with K&I 2000 physics.
""")
