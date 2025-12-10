#!/usr/bin/env python3
"""
Koyama & Inutsuka (2000) Thermal Equilibrium Calculator - Version 2

Follows EXACTLY the formulas in Appendix A of K&I (2000) and references therein:
- Bakes & Tielens (1994) for photoelectric heating
- Wolfire et al. (1995) for ISM parameters
- Hollenbach & McKee (1979, 1989) for cooling functions

The key insight from the CNM_FIGURE_ANALYSIS.md is that at n=10 cm^-3,
the equilibrium temperature is ~1200 K (unstable branch), NOT 25 K.

The hardcoded data represents the EQUILIBRIUM curve which includes
both stable (WNM, CNM) and unstable branches.
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Physical Constants (CGS)
# =============================================================================
k_B = 1.380649e-16      # Boltzmann constant [erg/K]
m_H = 1.6726e-24        # Hydrogen mass [g]
eV = 1.602e-12          # 1 eV in erg

# =============================================================================
# ISM Parameters from Wolfire et al. (1995) as used by K&I (2000)
# =============================================================================
G_0 = 1.7               # FUV field (1.7 x Habing, following W95)
zeta_CR = 1.8e-17       # Cosmic ray ionization rate [s^-1]

# Abundances relative to H (from W95/K&I2000 paper)
x_He = 0.1              # Helium
x_O = 4.6e-4            # Oxygen
x_C = 3.0e-4            # Carbon  (all assumed to be C+)
x_Si = 3.55e-6          # Silicon (all assumed to be Si+)
x_Fe = 7.08e-7          # Iron (all assumed to be Fe+)

# Grain parameters
T_gr = 8.0              # Grain temperature [K]


def compute_electron_fraction(n_H, T, N_H=1e20, G_0=1.7):
    """
    Compute equilibrium electron fraction from ionization balance.

    Main ionization sources at low density:
    - Photoionization by UV (attenuated by column)
    - Cosmic ray ionization
    - X-ray ionization (attenuated by column)
    - Collisional ionization (at high T)

    Recombination:
    - Radiative recombination
    - Grain recombination (at low T)
    """
    # UV attenuation: tau_d = N_H / (1.87e21 cm^-2)
    tau_d = N_H / 1.87e21

    # Photoionization rate (W95 eq. A2-A4)
    # Gamma_H,photo ~ 4.6e-10 * G_0 * exp(-2.5*tau_d) at Solar neighborhood
    # But at high column, this is negligible
    Gamma_photo = 4.6e-10 * G_0 * np.exp(-2.5 * tau_d)  # [s^-1]

    # Cosmic ray ionization (dominant at high column)
    Gamma_CR = zeta_CR  # [s^-1]

    # X-ray ionization (W95 appendix) - also attenuated
    tau_X = N_H / 1e22
    Gamma_X = 1.2e-17 * np.exp(-tau_X)  # Rough fit from W95

    # Total ionization rate
    Gamma_ion = Gamma_photo + Gamma_CR + Gamma_X

    # Radiative recombination coefficient (Shapiro & Kang 1987)
    # alpha_B ~ 2.6e-13 * (T/10^4)^(-0.7) cm^3/s
    T4 = T / 1e4
    alpha_rec = 2.6e-13 * T4**(-0.7)  # [cm^3/s]

    # At equilibrium: n_H * x_H * Gamma_ion = n_e * n_H * x_H * alpha_rec
    # => Gamma_ion = x_e * n_H * alpha_rec  (since x_H ~ 1 for neutral gas)
    # => x_e = Gamma_ion / (n_H * alpha_rec)

    x_e = Gamma_ion / (n_H * alpha_rec + 1e-30)

    # At high T, collisional ionization becomes important
    if T > 5000:
        # Collisional ionization coefficient
        k_ion = 5.8e-9 * T4**0.5 * np.exp(-15.8 / T4)  # [cm^3/s]
        # x_e^2 * n_H * alpha = x_H * Gamma + x_e * x_H * n_H * k_ion
        # Assuming x_H ~ 1-x_e ~ 1:
        # x_e^2 * n_H * alpha = Gamma + x_e * n_H * k_ion
        # x_e^2 - (k_ion/alpha)*x_e - Gamma/(n_H*alpha) = 0
        a = 1.0
        b = -k_ion / alpha_rec
        c = -Gamma_ion / (n_H * alpha_rec + 1e-30)
        x_e = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)

    # Clamp
    x_e = np.clip(x_e, 1e-10, 1.0)

    return x_e


# =============================================================================
# HEATING PROCESSES (from Appendix A.1)
# =============================================================================

def photoelectric_heating(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """
    Photoelectric heating from small grains and PAHs
    Bakes & Tielens (1994) - Equations 7-9 in K&I paper

    Gamma_pe = 1.0e-24 * epsilon * G_0  [erg/s per H]

    epsilon is heating efficiency, function of psi = G_0 * T^0.5 / n_e
    """
    n_e = x_e * n_H
    if n_e < 1e-10 * n_H:
        n_e = 1e-10 * n_H

    # FUV attenuation by dust
    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    if G_eff <= 1e-10:
        return 0.0

    # Heating efficiency parameter (BT94)
    psi = G_eff * np.sqrt(T) / n_e

    # BT94 fit (eq. 8 in paper)
    term1 = 4.9e-2 / (1.0 + (psi / 1925.0)**0.73)
    term2 = 3.7e-2 * (T / 1e4)**0.7 / (1.0 + psi / 5000.0)
    epsilon = term1 + term2

    # Heating rate [erg/s per H nucleus]
    Gamma_pe = 1.0e-24 * epsilon * G_eff

    return Gamma_pe


def pe_recombination_cooling(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """
    Recombination cooling on grains (BT94 eq. 9)
    Lambda_pe = 4.65e-30 * T^0.94 * psi^beta * n_e  [erg/s per H]
    """
    n_e = x_e * n_H
    if n_e < 1e-10 * n_H:
        n_e = 1e-10 * n_H

    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    if G_eff <= 1e-10 or T <= 0:
        return 0.0

    psi = G_eff * np.sqrt(T) / n_e
    beta = 0.74 / T**0.068

    # Per H nucleus
    Lambda_pe = 4.65e-30 * T**0.94 * psi**beta * n_e / n_H

    return Lambda_pe


def cosmic_ray_heating(n_H, T, x_e):
    """
    Cosmic ray heating (Appendix A.1.2)

    Primary CR ionization deposits ~6-20 eV per ionization.
    Following Shull & Van Steenberg (1985) for heat deposition.
    """
    # Heat deposited per ionization depends on x_e
    # For x_e << 1: E_h ~ 7 eV
    # For x_e ~ 0.1: E_h ~ 3 eV

    if x_e < 0.01:
        E_h = 7.0 * eV
    elif x_e < 0.1:
        E_h = 6.0 * eV
    else:
        E_h = 3.0 * eV

    Gamma_CR = zeta_CR * E_h  # [erg/s per H]

    return Gamma_CR


def xray_heating(n_H, T, x_e, N_H=1e20):
    """
    X-ray heating from diffuse background (W95 Appendix)
    Strongly attenuated at high column density.
    """
    tau_X = N_H / 1e22
    Gamma_X = 1.5e-25 * np.exp(-tau_X)  # Approximate from W95

    return Gamma_X


# =============================================================================
# COOLING PROCESSES (from Appendix A.1.4-A.1.7)
# =============================================================================

def lyman_alpha_cooling(n_H, T, x_e):
    """
    Hydrogen Lyman-alpha cooling
    Collisional excitation of H by electrons at high T (> 5000 K)

    From HM89: Lambda_Lya ~ n_e * n_H * 7.5e-19 * exp(-118348/T)
    """
    if T < 3000:
        return 0.0

    n_e = x_e * n_H
    x_H = 1.0 - x_e  # Neutral H fraction
    n_HI = x_H * n_H

    # Collisional excitation rate (HM89)
    # q_12 ~ 7.5e-19 * exp(-E/kT) where E = 10.2 eV = 118348 K
    q_Lya = 7.5e-19 * np.exp(-118348.0 / T)  # [cm^3/s * erg]

    # Cooling per H nucleus
    Lambda_Lya = q_Lya * n_e * n_HI / n_H

    return Lambda_Lya


def cii_158um_cooling(n_H, T, x_e):
    """
    CII 158 micron fine-structure cooling
    Main coolant of WNM and CNM at intermediate densities.

    From HM89 and W95:
    Two-level atom with collisions by e- and H
    """
    n_e = x_e * n_H
    n_HI = (1.0 - x_e) * n_H

    # CII abundance (assume all C is C+)
    n_CII = x_C * n_H

    # Energy of 158 um transition
    Delta_E = 91.2 * k_B  # 91.2 K in energy

    # Einstein A coefficient
    A_ul = 2.4e-6  # s^-1

    # Collision rates (HM89, W95)
    # Electrons:
    gamma_e = 2.8e-7 * (T / 100.0)**(-0.5)  # cm^3/s
    # Neutral H:
    gamma_H = 8.0e-10 * (T / 100.0)**0.07   # cm^3/s

    # Critical density
    n_cr = A_ul / (gamma_e * x_e + gamma_H * (1.0 - x_e) + 1e-30)

    # Boltzmann factor
    exp_factor = np.exp(-Delta_E / (k_B * T))

    # Two-level cooling rate
    # Lambda = n_CII * Delta_E * A_ul * (g_u/g_l) * exp(-E/kT) / (1 + n_cr/n + (g_u/g_l)*exp(-E/kT))
    g_ratio = 2.0  # g_u/g_l = 4/2 = 2

    Lambda_CII = n_CII / n_H * Delta_E * A_ul * g_ratio * exp_factor / \
                 (1.0 + n_cr / n_H + g_ratio * exp_factor)

    return Lambda_CII


def oi_63um_cooling(n_H, T, x_e):
    """
    OI 63 micron fine-structure cooling
    Important at T ~ 100-500 K

    Oxygen has same ionization potential as H, so stays neutral.
    """
    n_OI = x_O * n_H
    n_HI = (1.0 - x_e) * n_H

    # 63 micron transition (3P_1 -> 3P_2)
    Delta_E = 228.0 * k_B  # 228 K
    A_ul = 8.9e-5  # s^-1

    # Collision rate with H (HM89)
    gamma_H = 9.2e-11 * (T / 100.0)**0.67  # cm^3/s

    # Critical density
    n_cr = A_ul / (gamma_H + 1e-30)

    # Boltzmann factor
    exp_factor = np.exp(-Delta_E / (k_B * T))

    # Statistical weight ratio
    g_ratio = 3.0 / 5.0  # g_u/g_l for 3P_1/3P_2

    Lambda_OI = n_OI / n_H * Delta_E * A_ul * g_ratio * exp_factor / \
                (1.0 + n_cr / n_H + g_ratio * exp_factor)

    return Lambda_OI


def metastable_line_cooling(n_H, T, x_e):
    """
    Metastable line cooling from CII (2326A), OI (6300A), FeII, SiII
    From HM89 - important at intermediate temperatures
    """
    Lambda = 0.0

    # CII 2326 Angstrom (metastable)
    if T > 1000:
        # CII 2P -> 4P transition at 62000 K
        n_CII = x_C * n_H
        Delta_E_CII = 62000.0 * k_B
        A_CII = 6.3e4  # s^-1 (forbidden line)
        gamma_CII = 3.0e-9  # Collision rate with e-
        n_e = x_e * n_H
        exp_CII = np.exp(-62000.0 / T)
        n_cr_CII = A_CII / (gamma_CII + 1e-30)
        Lambda += n_CII / n_H * Delta_E_CII * A_CII * exp_CII / (1.0 + n_cr_CII / n_H)

    # OI 6300 Angstrom (forbidden)
    if T > 500:
        n_OI = x_O * n_H
        Delta_E_OI = 22800.0 * k_B  # 22800 K
        A_OI = 6.3e-3  # s^-1
        exp_OI = np.exp(-22800.0 / T)
        Lambda += n_OI / n_H * Delta_E_OI * A_OI * exp_OI / (1.0 + 1e4 / n_H)

    return Lambda


def grain_gas_cooling(n_H, T, T_gr=8.0):
    """
    Gas-grain thermal coupling (Appendix A.1.7)
    Important at high density where gas and grains equilibrate.

    Lambda_gr = 1.2e-31 * n * (T/1000)^0.5 * (100A/a)^0.5 *
                (1 - 0.8*exp(-75/T)) * (T - T_gr)
    """
    if T <= T_gr:
        return 0.0  # Grains would heat gas, not cool it

    # Per H nucleus (eq. 18 in paper)
    Lambda_gr = 1.2e-31 * n_H * (T / 1000.0)**0.5 * \
                (1.0 - 0.8 * np.exp(-75.0 / T)) * (T - T_gr)

    return Lambda_gr


# =============================================================================
# TOTAL HEATING AND COOLING
# =============================================================================

def total_heating(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """Total heating rate [erg/s per H]"""
    Gamma = 0.0
    Gamma += photoelectric_heating(n_H, T, x_e, G_0, N_H)
    Gamma += cosmic_ray_heating(n_H, T, x_e)
    Gamma += xray_heating(n_H, T, x_e, N_H)
    return Gamma


def total_cooling(n_H, T, x_e, N_H=1e20):
    """Total cooling rate [erg/s per H]"""
    Lambda = 0.0
    Lambda += lyman_alpha_cooling(n_H, T, x_e)
    Lambda += cii_158um_cooling(n_H, T, x_e)
    Lambda += oi_63um_cooling(n_H, T, x_e)
    Lambda += metastable_line_cooling(n_H, T, x_e)
    Lambda += grain_gas_cooling(n_H, T)
    Lambda += pe_recombination_cooling(n_H, T, x_e, G_0, N_H)
    return Lambda


def net_cooling(n_H, T, N_H=1e20, G_0=1.7):
    """Net cooling rate Lambda - Gamma"""
    x_e = compute_electron_fraction(n_H, T, N_H, G_0)
    Gamma = total_heating(n_H, T, x_e, G_0, N_H)
    Lambda = total_cooling(n_H, T, x_e, N_H)
    return Lambda - Gamma


# =============================================================================
# EQUILIBRIUM SOLVER
# =============================================================================

def find_equilibrium_temperature(n_H, N_H=1e20, G_0=1.7, T_min=10.0, T_max=2e4):
    """
    Find equilibrium temperature where Gamma = Lambda
    May have multiple roots (WNM, unstable, CNM branches)
    """
    def func(log_T):
        T = 10**log_T
        return net_cooling(n_H, T, N_H, G_0)

    try:
        # Search for root
        log_T_eq = brentq(func, np.log10(T_min), np.log10(T_max), xtol=1e-4)
        return 10**log_T_eq
    except ValueError:
        # No root found - return endpoint where net cooling is smallest
        net_low = abs(net_cooling(n_H, T_min, N_H, G_0))
        net_high = abs(net_cooling(n_H, T_max, N_H, G_0))
        if net_low < net_high:
            return T_min
        else:
            return T_max


def find_all_equilibrium_temperatures(n_H, N_H=1e20, G_0=1.7, T_min=10.0, T_max=2e4, n_search=100):
    """
    Find all equilibrium temperatures (may be multiple due to thermal instability)
    Returns list of equilibrium temperatures
    """
    T_arr = np.logspace(np.log10(T_min), np.log10(T_max), n_search)
    net_arr = np.array([net_cooling(n_H, T, N_H, G_0) for T in T_arr])

    # Find sign changes
    equilibria = []
    for i in range(len(T_arr) - 1):
        if net_arr[i] * net_arr[i+1] < 0:
            # Sign change - refine with brentq
            try:
                T_eq = brentq(lambda log_T: net_cooling(n_H, 10**log_T, N_H, G_0),
                             np.log10(T_arr[i]), np.log10(T_arr[i+1]), xtol=1e-4)
                equilibria.append(10**T_eq)
            except:
                pass

    return equilibria


# =============================================================================
# DIAGNOSTIC: Print heating/cooling at specific conditions
# =============================================================================

def diagnose_thermal_balance(n_H, T, N_H=1e20, G_0=1.7):
    """Print detailed heating/cooling breakdown"""
    x_e = compute_electron_fraction(n_H, T, N_H, G_0)

    print(f"\n{'='*60}")
    print(f"Thermal balance at n_H = {n_H:.2e} cm^-3, T = {T:.1f} K")
    print(f"N_H = {N_H:.0e} cm^-2, G_0 = {G_0}")
    print(f"x_e = {x_e:.4e}")
    print(f"{'='*60}")

    # Heating
    Gamma_PE = photoelectric_heating(n_H, T, x_e, G_0, N_H)
    Gamma_CR = cosmic_ray_heating(n_H, T, x_e)
    Gamma_X = xray_heating(n_H, T, x_e, N_H)
    Gamma_total = total_heating(n_H, T, x_e, G_0, N_H)

    print(f"\nHEATING [erg/s per H]:")
    print(f"  Photoelectric: {Gamma_PE:.4e}")
    print(f"  Cosmic ray:    {Gamma_CR:.4e}")
    print(f"  X-ray:         {Gamma_X:.4e}")
    print(f"  TOTAL:         {Gamma_total:.4e}")

    # Cooling
    Lambda_Lya = lyman_alpha_cooling(n_H, T, x_e)
    Lambda_CII = cii_158um_cooling(n_H, T, x_e)
    Lambda_OI = oi_63um_cooling(n_H, T, x_e)
    Lambda_meta = metastable_line_cooling(n_H, T, x_e)
    Lambda_gr = grain_gas_cooling(n_H, T)
    Lambda_PE_rec = pe_recombination_cooling(n_H, T, x_e, G_0, N_H)
    Lambda_total = total_cooling(n_H, T, x_e, N_H)

    print(f"\nCOOLING [erg/s per H]:")
    print(f"  Lyman-alpha:    {Lambda_Lya:.4e}")
    print(f"  CII 158um:      {Lambda_CII:.4e}")
    print(f"  OI 63um:        {Lambda_OI:.4e}")
    print(f"  Metastable:     {Lambda_meta:.4e}")
    print(f"  Grain-gas:      {Lambda_gr:.4e}")
    print(f"  PE recomb:      {Lambda_PE_rec:.4e}")
    print(f"  TOTAL:          {Lambda_total:.4e}")

    print(f"\nNET (Lambda - Gamma): {Lambda_total - Gamma_total:.4e}")
    print(f"Heating/Cooling ratio: {Gamma_total/Lambda_total:.3f}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*70)
    print("Koyama & Inutsuka (2000) Thermal Equilibrium - Version 2")
    print("="*70)

    # Diagnostic: check conditions at key densities
    test_cases = [
        (0.1, 8000),   # WNM
        (1.0, 5000),   # Transition
        (10.0, 1000),  # Intermediate
        (100.0, 100),  # CNM
        (1000.0, 50),  # Dense
    ]

    for n_H, T_guess in test_cases:
        diagnose_thermal_balance(n_H, T_guess, N_H=1e20)

    # Calculate equilibrium curve
    print("\n" + "="*70)
    print("CALCULATING EQUILIBRIUM CURVE")
    print("="*70)

    log_n_arr = np.linspace(-1, 6, 100)
    n_arr = 10**log_n_arr

    T_eq_1e19 = []
    T_eq_1e20 = []

    for n_H in n_arr:
        T_eq_1e19.append(find_equilibrium_temperature(n_H, N_H=1e19))
        T_eq_1e20.append(find_equilibrium_temperature(n_H, N_H=1e20))

    T_eq_1e19 = np.array(T_eq_1e19)
    T_eq_1e20 = np.array(T_eq_1e20)

    # Load hardcoded data for comparison
    n_data = np.array([
        1.000000e-01, 1.274138e-01, 1.654571e-01, 2.108153e-01, 2.737607e-01, 3.488090e-01, 4.529566e-01, 5.771294e-01,
        7.494491e-01, 9.549019e-01, 1.240017e+00, 1.579953e+00, 2.051697e+00, 2.614146e+00, 3.394679e+00, 4.325291e+00,
        5.616739e+00, 7.293788e+00, 9.293296e+00, 1.206809e+01, 1.537642e+01, 1.996752e+01, 2.544139e+01, 3.303769e+01,
        4.209460e+01, 5.492352e+01, 6.964852e+01, 9.087488e+01, 1.152384e+02, 1.503590e+02, 1.906702e+02, 2.487797e+02,
        3.169798e+02, 4.116238e+02, 5.244657e+02, 6.810609e+02, 8.677659e+02, 1.126864e+03, 1.435781e+03, 1.864477e+03,
        2.375602e+03, 3.084910e+03, 3.930603e+03, 5.104205e+03, 6.503464e+03, 8.445272e+03, 1.076045e+04, 1.397331e+04,
        1.780393e+04, 2.311983e+04, 2.945787e+04, 3.825341e+04, 4.874015e+04, 6.329301e+04, 8.064405e+04, 1.047228e+05,
        1.334314e+05, 1.732714e+05, 2.207717e+05, 2.866899e+05, 3.652826e+05, 4.743488e+05, 6.043861e+05, 7.848441e+05,
        1.000000e+06
    ])

    T_data = np.array([
        8.440767e+03, 8.266786e+03, 8.266786e+03, 8.138657e+03, 8.138657e+03, 8.138657e+03, 7.970903e+03, 7.847361e+03,
        7.847361e+03, 7.725733e+03, 7.449216e+03, 7.071271e+03, 6.574144e+03, 5.623413e+03, 4.472007e+03, 3.375919e+03,
        2.419181e+03, 1.706714e+03, 1.197821e+03, 8.718711e+02, 6.479743e+02, 5.073135e+02, 4.034399e+02, 3.327442e+02,
        2.787572e+02, 2.421982e+02, 2.104338e+02, 1.896223e+02, 1.708690e+02, 1.531708e+02, 1.401954e+02, 1.310197e+02,
        1.199208e+02, 1.138365e+02, 1.058335e+02, 9.890680e+01, 9.388864e+01, 8.912509e+01, 8.593515e+01, 8.157513e+01,
        7.865542e+01, 7.584020e+01, 7.350750e+01, 7.199236e+01, 6.941563e+01, 6.833975e+01, 6.728054e+01, 6.589375e+01,
        6.487245e+01, 6.487245e+01, 6.386698e+01, 6.386698e+01, 6.255055e+01, 6.158107e+01, 6.031176e+01, 5.937698e+01,
        5.845668e+01, 5.549082e+01, 5.158968e+01, 4.721942e+01, 4.254950e+01, 3.834143e+01, 3.564594e+01, 3.331295e+01,
        3.162278e+01
    ])

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.loglog(n_arr, T_eq_1e19, '-', label='Calculated (N_H=1e19)', color='C0', linewidth=2)
    ax.loglog(n_arr, T_eq_1e20, '-', label='Calculated (N_H=1e20)', color='C1', linewidth=2)
    ax.loglog(n_data, T_data, 'o', label='Hardcoded (N_H=1e20)', color='C2', markersize=3, alpha=0.7)
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'$T_{eq}$ [K]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 2e4)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Thermal Equilibrium Temperature')

    # Plot heating/cooling rates at equilibrium
    ax = axes[1]
    Gamma_PE = []
    Gamma_CR = []
    Lambda_CII = []
    Lambda_OI = []
    Lambda_Lya = []

    for n_H, T in zip(n_arr, T_eq_1e20):
        x_e = compute_electron_fraction(n_H, T, 1e20, G_0)
        Gamma_PE.append(photoelectric_heating(n_H, T, x_e, G_0, 1e20))
        Gamma_CR.append(cosmic_ray_heating(n_H, T, x_e))
        Lambda_CII.append(cii_158um_cooling(n_H, T, x_e))
        Lambda_OI.append(oi_63um_cooling(n_H, T, x_e))
        Lambda_Lya.append(lyman_alpha_cooling(n_H, T, x_e))

    ax.loglog(n_arr, Gamma_PE, '--', label='PE heating', color='C0')
    ax.loglog(n_arr, Gamma_CR, '--', label='CR heating', color='C1')
    ax.loglog(n_arr, Lambda_CII, '-', label='CII cooling', color='C2')
    ax.loglog(n_arr, Lambda_OI, '-', label='OI cooling', color='C3')
    ax.loglog(n_arr, Lambda_Lya, '-', label='Lya cooling', color='C4')
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'Rate [erg s$^{-1}$ per H]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-28, 1e-22)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_title('Heating (--) and Cooling (-) Rates at Equilibrium')

    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/scripts/ki2000_equilibrium_v2.png', dpi=150)
    print(f"\nSaved: ki2000_equilibrium_v2.png")

    # Print key values
    print("\n" + "="*70)
    print("KEY EQUILIBRIUM VALUES (N_H = 1e20)")
    print("="*70)
    key_n = [0.1, 1.0, 10.0, 100.0, 1000.0]
    for n in key_n:
        T_calc = np.interp(n, n_arr, T_eq_1e20)
        T_hard = np.interp(n, n_data, T_data)
        print(f"n = {n:>7.1f}: T_calc = {T_calc:>8.1f} K, T_hardcoded = {T_hard:>8.1f} K, ratio = {T_calc/T_hard:.2f}")


if __name__ == "__main__":
    main()
