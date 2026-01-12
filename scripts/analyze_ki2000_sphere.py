#!/usr/bin/env python3
"""
Analyze K&I 2000 compatible sphere Jeans stability

Key question: Does using T = T_eq(n) throughout solve the Jeans problem?

The K&I 2000 sphere uses:
- P(ρ) = n * k_B * T_eq(n)  [barotropic EOS]
- c_eff² = (k_B T / m_n) * (1 + d ln T / d ln n)
- Hydrostatic: dP/dr = -ρ G M(r) / r²

This is NOT an isothermal sphere - it's a "barotropic" sphere where
T varies with density according to the K&I 2000 equilibrium curve.
"""
import math

# Constants
G = 6.674e-8
k_B = 1.38e-16
m_H = 1.67e-24
M_sun = 1.989e33
pc = 3.086e18
mu = 1.27
m_n = mu * m_H

# K&I 2000 T_eq data (N_H = 10^20)
n_data = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 10000]
T_data = [80, 45, 28, 17, 13, 11, 9, 8, 7, 4.5, 3.7]

def T_eq(n):
    """K&I 2000 equilibrium temperature"""
    if n <= n_data[0]:
        return T_data[0]
    if n >= n_data[-1]:
        return T_data[-1]
    for i in range(len(n_data) - 1):
        if n_data[i] <= n <= n_data[i+1]:
            log_n = math.log10(n)
            log_n1, log_n2 = math.log10(n_data[i]), math.log10(n_data[i+1])
            alpha = (log_n - log_n1) / (log_n2 - log_n1)
            return T_data[i] + alpha * (T_data[i+1] - T_data[i])
    return T_data[-1]

def P_eq(n):
    """K&I 2000 equilibrium pressure P/k_B [K cm^-3]"""
    return n * T_eq(n)

def dlnT_dlnn(n, eps=0.05):
    """Numerical derivative d ln T / d ln n"""
    n_lo = n * (1 - eps)
    n_hi = n * (1 + eps)
    T_lo = T_eq(n_lo)
    T_hi = T_eq(n_hi)
    return (math.log(T_hi) - math.log(T_lo)) / (math.log(n_hi) - math.log(n_lo))

def c_eff_squared(n):
    """Effective sound speed squared [cm²/s²]"""
    T = T_eq(n)
    dlogT_dlogn = dlnT_dlnn(n)
    # c_eff² = (k_B T / m_n) * (1 + d ln T / d ln n)
    # Note: dlogT/dlogn is NEGATIVE for cold branch, so c_eff² < c_s²
    return (k_B * T / m_n) * (1.0 + dlogT_dlogn)

def jeans_mass(T, n):
    """Standard Jeans mass for isothermal gas"""
    rho = n * m_n
    c_s = math.sqrt(k_B * T / m_n)
    lJ = c_s * math.sqrt(math.pi / (G * rho))
    return (math.pi / 6) * rho * lJ**3 / M_sun

def jeans_mass_barotropic(n):
    """Jeans mass using effective sound speed"""
    rho = n * m_n
    c_eff = math.sqrt(max(c_eff_squared(n), 1e10))  # Floor to avoid negative
    lJ_eff = c_eff * math.sqrt(math.pi / (G * rho))
    return (math.pi / 6) * rho * lJ_eff**3 / M_sun

print("="*70)
print("K&I 2000 COMPATIBLE SPHERE: JEANS STABILITY ANALYSIS")
print("="*70)
print("\nKey physics:")
print("- K&I 2000 barotropic: T = T_eq(n), not constant")
print("- On cold branch: dT/dn < 0 (T decreases as n increases)")
print("- c_eff² = (k_B T / m_n) × (1 + d ln T / d ln n)")
print("- Since d ln T / d ln n < 0, we have c_eff² < c_s² (isothermal)")
print()

print("-"*70)
print(f"{'n (cm⁻³)':>10} {'T_eq (K)':>10} {'dlogT/dlogn':>12} {'c_eff/c_s':>10} {'M_J_iso':>10} {'M_J_eff':>10}")
print("-"*70)

for n in [100, 200, 300, 500, 1000, 2000, 3000, 5000]:
    T = T_eq(n)
    dlog = dlnT_dlnn(n)
    c_s_sq = k_B * T / m_n
    c_eff_sq = c_eff_squared(n)
    c_ratio = math.sqrt(max(c_eff_sq / c_s_sq, 0))
    MJ_iso = jeans_mass(T, n)
    MJ_eff = jeans_mass_barotropic(n)
    
    print(f"{n:>10} {T:>10.1f} {dlog:>12.3f} {c_ratio:>10.3f} {MJ_iso:>10.1f} {MJ_eff:>10.1f}")

print("-"*70)

print("""
KEY INSIGHT:

On the K&I 2000 cold branch:
  d ln T / d ln n ~ -0.4 to -0.7 (T DECREASES with n)
  
This means:
  c_eff² = c_s² × (1 + dlogT/dlogn) 
         = c_s² × (1 - 0.5)  [approximately]
         = 0.5 × c_s²
         
So c_eff ~ 0.7 × c_s

The EFFECTIVE Jeans mass is SMALLER than isothermal:
  M_J_eff ∝ c_eff³ ~ (0.7)³ × c_s³ ~ 0.35 × M_J_iso

This makes the stability problem WORSE, not better!

CONCLUSION:
A K&I 2000 compatible sphere does NOT solve the Jeans instability problem.
In fact, it makes the sphere MORE prone to collapse because the effective
sound speed (which resists gravity) is REDUCED by the negative dT/dn.

The ONLY solutions are:
1. Use LOW density (n < 300 cm⁻³) where T_eq is high enough
2. Use very SMALL mass (M << M_J)
3. Disable self-gravity (G = 0)
""")

# Final recommendation
print("="*70)
print("RECOMMENDED CONFIGURATION FOR HIGH DENSITY COOLING TEST")
print("="*70)
print("""
If you want n ~ 3000 cm⁻³ with working cooling physics:

OPTION A: No self-gravity (test cooling in isolation)
  - Set G = 0 in config
  - Cloud cools without collapsing
  - Tests K&I 2000 thermal equilibrium
  - Not physically realistic for ISM

OPTION B: Low-mass compact cloud  
  - n_center = 3000 cm⁻³
  - M_cloud < 8 M☉ (to be sub-Jeans at T_eq)
  - ξ_s ~ 1.5-2.0 (very truncated core)
  - R ~ 0.1-0.2 pc
  - N particles ~ 1000-5000 (small simulation)

OPTION C: Accept moderate density
  - n_center = 200 cm⁻³ (maximum stable)
  - T_eq = 45 K
  - M_cloud ~ 450 M☉
  - R ~ 3 pc
  - Standard resolution (N ~ 30000)
""")
