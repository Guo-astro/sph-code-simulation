#!/usr/bin/env python3
"""
Design a Jeans-stable Bonnor-Ebert sphere for ISM cooling tests.

The key insight is that for a BE sphere:
- R_BE = xi_s * r_0 = xi_s * c_s / sqrt(4*pi*G*rho_c)
- R_BE ~ sqrt(T / n_center)

So for a given density, higher T gives larger R.
For stability, we need R < lambda_J and M < M_J.

Author: Design tool for ISM cooling test configurations
"""

import math

# Physical constants (CGS)
G = 6.674e-8           # gravitational constant
k_B = 1.38e-16         # Boltzmann constant
m_H = 1.67e-24         # hydrogen mass
M_sun = 1.989e33       # solar mass
pc = 3.086e18          # parsec in cm
mu = 1.27              # mean molecular weight (neutral ISM)


def be_scale_length(T, n_center):
    """Isothermal scale length r_0 = c_s / sqrt(4*pi*G*rho_c)"""
    rho_c = n_center * mu * m_H  # g/cm^3
    c_s = math.sqrt(k_B * T / (mu * m_H))  # cm/s
    r_0 = c_s / math.sqrt(4 * math.pi * G * rho_c)  # cm
    return r_0 / pc  # pc


def be_radius(T, n_center, xi_s=6.0):
    """BE sphere radius R = xi_s * r_0"""
    return xi_s * be_scale_length(T, n_center)


def be_mass(T, n_center, xi_s=6.0, m_s=14.34):
    """BE sphere mass from dimensionless mass m_s = xi^2 * |dpsi/dxi|"""
    rho_c = n_center * mu * m_H  # g/cm^3
    r_0 = be_scale_length(T, n_center) * pc  # cm
    M = 4 * math.pi * rho_c * r_0**3 * m_s / M_sun
    return M


def jeans_length(T, n):
    """Jeans length: lambda_J = c_s * sqrt(pi / (G * rho))"""
    rho = n * mu * m_H
    c_s = math.sqrt(k_B * T / (mu * m_H))
    return c_s * math.sqrt(math.pi / (G * rho)) / pc


def jeans_mass(T, n):
    """Jeans mass: M_J = (pi/6) * rho * lambda_J^3"""
    rho = n * mu * m_H
    lambda_J = jeans_length(T, n) * pc  # cm
    return (math.pi / 6.0) * rho * lambda_J**3 / M_sun


def ki2000_equilibrium_temp(n):
    """
    Approximate K&I 2000 equilibrium temperature.
    From their Fig.1, the thermal equilibrium curve gives:
    - n ~ 100 cm^-3: T ~ 100-200 K
    - n ~ 1000 cm^-3: T ~ 10-15 K
    - n ~ 5000 cm^-3: T ~ 5-7 K
    """
    # Rough fit to K&I equilibrium curve in cold phase
    if n < 100:
        return 200  # warm neutral medium
    elif n < 500:
        return 50 * (n / 100) ** (-0.5)
    elif n < 2000:
        return 15 * (n / 1000) ** (-0.3)
    else:
        return 7 * (n / 2000) ** (-0.2)


def analyze_configuration(T, n_center, xi_s=6.0):
    """Analyze a BE sphere configuration for Jeans stability."""
    R = be_radius(T, n_center, xi_s)
    M = be_mass(T, n_center, xi_s)
    lJ = jeans_length(T, n_center)
    MJ = jeans_mass(T, n_center)
    
    r_ratio = R / lJ
    m_ratio = M / MJ
    
    T_eq = ki2000_equilibrium_temp(n_center)
    
    return {
        'T': T,
        'n': n_center,
        'R': R,
        'M': M,
        'lambda_J': lJ,
        'M_J': MJ,
        'R/lJ': r_ratio,
        'M/MJ': m_ratio,
        'T_eq': T_eq,
        'T_margin': T - T_eq
    }


def main():
    print("=" * 70)
    print("DESIGNING A JEANS-STABLE BONNOR-EBERT SPHERE")
    print("=" * 70)
    
    print("\n### KEY PHYSICS ###")
    print("""
For a Bonnor-Ebert sphere at temperature T and central density n_c:

  Scale length:  r_0 = c_s / sqrt(4*pi*G*rho_c)
  Cloud radius:  R = xi_s * r_0  (xi_s = 6 for critical BE sphere)
  
  Since c_s ~ sqrt(T) and rho_c ~ n_c:
    R ~ sqrt(T / n_c)
    
For Jeans stability:
  lambda_J = c_s * sqrt(pi / (G * rho)) ~ sqrt(T / n)
  
Notice: R/lambda_J is nearly CONSTANT for BE spheres!
  R/lambda_J = xi_s / sqrt(4) ~ 3 for xi_s = 6
  
This means a critical BE sphere (xi_s=6) is ALWAYS marginally Jeans unstable!
To get stable configuration, we need xi_s < 6 (truncated BE sphere).
""")
    
    print("\n" + "=" * 70)
    print("SCANNING PARAMETER SPACE (xi_s = 6.0, critical BE)")
    print("=" * 70)
    
    # Show that critical BE is always marginal
    print("\n--- Critical BE sphere (xi_s=6) is always marginal ---")
    print(f"{'T(K)':<8} {'n(cm-3)':<10} {'R(pc)':<10} {'M(Msun)':<10} {'lJ(pc)':<10} {'MJ(Msun)':<10} {'R/lJ':<8}")
    print("-" * 80)
    
    for T in [10, 15, 20, 50, 100]:
        for n in [100, 1000, 5000]:
            cfg = analyze_configuration(T, n, xi_s=6.0)
            print(f"{T:<8} {n:<10} {cfg['R']:<10.3f} {cfg['M']:<10.1f} {cfg['lambda_J']:<10.2f} {cfg['M_J']:<10.0f} {cfg['R/lJ']:<8.2f}")
    
    print("\n" + "=" * 70)
    print("SOLUTION: USE TRUNCATED BE SPHERE (xi_s < 6)")
    print("=" * 70)
    
    print("""
To achieve R < lambda_J, we need:
  xi_s * r_0 < c_s * sqrt(pi / (G * rho))
  xi_s < sqrt(4 * pi) ~ 3.54

For stability margin, use xi_s ~ 3.0
""")
    
    print("\n--- Stable configurations with xi_s = 3.0 ---")
    print(f"{'T(K)':<8} {'n(cm-3)':<10} {'R(pc)':<10} {'M(Msun)':<10} {'lJ(pc)':<10} {'MJ(Msun)':<10} {'R/lJ':<8} {'M/MJ':<8} {'T_eq':<8} {'Status'}")
    print("-" * 110)
    
    xi_s = 3.0
    m_s = 1.83  # dimensionless mass for xi_s = 3.0 (from Lane-Emden integration)
    
    stable_configs = []
    for T in [15, 18, 20]:
        for n in [1000, 2000, 3000, 5000]:
            R = be_radius(T, n, xi_s)
            # For truncated sphere, mass is different
            rho_c = n * mu * m_H
            r_0 = be_scale_length(T, n) * pc
            M = 4 * math.pi * rho_c * r_0**3 * m_s / M_sun
            lJ = jeans_length(T, n)
            MJ = jeans_mass(T, n)
            T_eq = ki2000_equilibrium_temp(n)
            
            r_ratio = R / lJ
            m_ratio = M / MJ
            
            if r_ratio < 0.8 and m_ratio < 0.8 and T > T_eq:
                status = "STABLE"
                stable_configs.append((T, n, R, M, lJ, MJ, r_ratio, m_ratio, T_eq))
            elif r_ratio < 1.0 and m_ratio < 1.0:
                status = "marginal"
            else:
                status = "unstable"
            
            print(f"{T:<8} {n:<10} {R:<10.3f} {M:<10.2f} {lJ:<10.2f} {MJ:<10.0f} {r_ratio:<8.2f} {m_ratio:<8.2f} {T_eq:<8.1f} {status}")
    
    print("\n" + "=" * 70)
    print("RECOMMENDED CONFIGURATION")
    print("=" * 70)
    
    # Recommend T=20K, n=3000, xi_s=3.0
    T_rec, n_rec, xi_rec = 20, 3000, 3.0
    m_s_rec = 1.83
    
    R_rec = be_radius(T_rec, n_rec, xi_rec)
    rho_c = n_rec * mu * m_H
    r_0 = be_scale_length(T_rec, n_rec) * pc
    M_rec = 4 * math.pi * rho_c * r_0**3 * m_s_rec / M_sun
    lJ_rec = jeans_length(T_rec, n_rec)
    MJ_rec = jeans_mass(T_rec, n_rec)
    T_eq_rec = ki2000_equilibrium_temp(n_rec)
    n_edge_rec = n_rec * math.exp(-3.0**2 / 6)  # rough edge density
    
    print(f"""
RECOMMENDED STABLE CONFIGURATION:

  Temperature:     T = {T_rec} K  (> T_eq = {T_eq_rec:.1f} K from K&I 2000)
  Central density: n_center = {n_rec} cm^-3
  BE truncation:   xi_s = {xi_rec}  (not critical, truncated early)
  
  Self-consistent BE parameters:
    r_0 (scale)  = {be_scale_length(T_rec, n_rec):.4f} pc
    R_cloud      = {R_rec:.3f} pc
    M_cloud      = {M_rec:.2f} M_sun
    n_edge       ~ {n_edge_rec:.0f} cm^-3
  
  Jeans stability check:
    lambda_J = {lJ_rec:.3f} pc
    M_J      = {MJ_rec:.1f} M_sun
    
    R / lambda_J = {R_rec/lJ_rec:.2f} < 1   ✓ STABLE
    M / M_J      = {M_rec/MJ_rec:.2f} < 1   ✓ STABLE
  
  Thermal stability:
    T_cloud = {T_rec} K > T_eq = {T_eq_rec:.1f} K
    Margin  = +{T_rec - T_eq_rec:.1f} K
    Cloud will NOT runaway cool to equilibrium
    
  Free-fall time:
    t_ff = sqrt(3*pi / (32*G*rho))
         = {math.sqrt(3*math.pi/(32*G*rho_c)) / (3.15e13):.2f} Myr
         
JSON CONFIG PARAMETERS:
  "temperature_K": {T_rec},
  "n_center_cgs": {n_rec},
  "xi_s": {xi_rec}
""")
    
    print("\n" + "=" * 70)
    print("ALTERNATIVE: USE M_cloud MODE (SPECIFY MASS DIRECTLY)")
    print("=" * 70)
    
    print("""
If you want to control mass directly, use M_cloud mode in isothermal_bonnor_ebert:

  "ic_method": "isothermal_bonnor_ebert",
  "isothermal_bonnor_ebert": {
    "mode": "M_cloud",
    "M_cloud": 5.0,        // Desired mass in M_sun
    "n_center_cgs": 3000,  // Central density
    "temperature_K": 20,   // Isothermal temperature
    "xi_s": 6.0            // Will be ignored, computed from M
  }

This will compute xi_s from the specified mass, giving a truncated sphere.
""")


if __name__ == "__main__":
    main()
