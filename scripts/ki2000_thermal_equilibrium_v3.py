#!/usr/bin/env python3
"""
Koyama & Inutsuka (2000) Thermal Equilibrium Calculator - Version 3

CRITICAL FIX: The electron fraction must be properly calculated.
At WNM conditions (n~0.1, T~8000K), x_e ~ 0.01-0.1, NOT 1.0!

The key physics:
- At low density, ionization is by UV + cosmic rays
- At high column (N_H=1e20), UV is blocked, CR dominates
- Recombination is radiative + grain recombination
- For WNM: x_e ~ 0.01-0.1 (from Fig 1b of paper)
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
G_0 = 1.7               # FUV field (1.7 x Habing)
zeta_CR = 1.8e-17       # Cosmic ray ionization rate [s^-1]

# Abundances relative to H
x_He = 0.1
x_O = 4.6e-4
x_C = 3.0e-4
x_Si = 3.55e-6
x_Fe = 7.08e-7

T_gr = 8.0              # Grain temperature [K]


def compute_electron_fraction(n_H, T, N_H=1e20, G_0=1.7):
    """
    Equilibrium electron fraction.

    Key point: For neutral ISM, x_e << 1 (typically 1e-4 to 0.1)
    Ionization sources: UV (at low column), CR (always), X-ray (at low column)
    Recombination: radiative + grain-assisted

    From Fig 1b of K&I 2000:
    - At n=0.1 cm^-3: x_e ~ 0.1
    - At n=10 cm^-3: x_e ~ 0.001
    - At n=1000 cm^-3: x_e ~ 0.0003
    """
    # Attenuation factors
    tau_d = N_H / 1.87e21  # Dust optical depth
    tau_X = N_H / 1e22      # X-ray optical depth

    # Ionization rates [s^-1]
    # UV photoionization (strongly attenuated at high column)
    Gamma_UV = 4.6e-10 * G_0 * np.exp(-2.5 * tau_d)

    # Cosmic ray ionization (penetrating, always present)
    Gamma_CR = zeta_CR

    # X-ray ionization (attenuated)
    Gamma_X = 1.2e-17 * np.exp(-tau_X)

    # Total ionization rate
    Gamma_ion = Gamma_UV + Gamma_CR + Gamma_X

    # Recombination coefficient [cm^3/s]
    # Radiative recombination (case B)
    T4 = T / 1e4
    alpha_rad = 2.6e-13 * T4**(-0.7)

    # Grain-assisted recombination (important at low T)
    # From Weingartner & Draine (2001), roughly:
    # alpha_gr ~ 1e-14 at T ~ 100K, ~1e-15 at T ~ 1000K
    alpha_gr = 1e-14 * (100.0 / max(T, 10.0))**0.5

    # Total recombination
    alpha_rec = alpha_rad + alpha_gr

    # Equilibrium: x_H * Gamma_ion = x_e * n_H * alpha_rec
    # For neutral gas x_H ~ 1: x_e = Gamma_ion / (n_H * alpha_rec)
    x_e = Gamma_ion / (n_H * alpha_rec + 1e-30)

    # Clamp to physical range (neutral ISM: x_e < 1)
    x_e = np.clip(x_e, 1e-10, 0.999)

    return x_e


# =============================================================================
# HEATING
# =============================================================================

def photoelectric_heating(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """
    Photoelectric heating from grains (BT94)
    Gamma_pe = 1.0e-24 * epsilon * G_0  [erg/s per H]
    """
    n_e = x_e * n_H
    if n_e < 1e-6:
        n_e = 1e-6

    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    if G_eff < 1e-10:
        return 0.0

    # Heating efficiency (BT94)
    psi = G_eff * np.sqrt(T) / n_e

    term1 = 4.9e-2 / (1.0 + (psi / 1925.0)**0.73)
    term2 = 3.7e-2 * (T / 1e4)**0.7 / (1.0 + psi / 5000.0)
    epsilon = term1 + term2

    Gamma_pe = 1.0e-24 * epsilon * G_eff

    return Gamma_pe


def cosmic_ray_heating(n_H, T, x_e):
    """CR heating: ~7 eV per ionization at low x_e"""
    if x_e < 0.01:
        E_h = 7.0 * eV
    else:
        E_h = 5.0 * eV
    return zeta_CR * E_h


def xray_heating(n_H, T, N_H=1e20):
    """X-ray heating (attenuated)"""
    tau_X = N_H / 1e22
    return 1.5e-25 * np.exp(-tau_X)


# =============================================================================
# COOLING
# =============================================================================

def lyman_alpha_cooling(n_H, T, x_e):
    """Lya cooling at high T"""
    if T < 5000:
        return 0.0

    n_e = x_e * n_H
    x_H = 1.0 - x_e
    n_HI = x_H * n_H

    # HM89 rate
    Lambda_Lya = 7.5e-19 * np.exp(-118348.0 / T) * n_e * n_HI / n_H

    return Lambda_Lya


def cii_158um_cooling(n_H, T, x_e):
    """
    CII 158 micron fine-structure cooling
    This is THE dominant coolant for WNM/CNM
    """
    n_e = x_e * n_H
    n_HI = (1.0 - x_e) * n_H
    n_CII = x_C * n_H

    # Transition parameters
    Delta_E = 91.2 * k_B  # Energy [erg]
    A_ul = 2.4e-6         # Einstein A [s^-1]
    g_u = 4               # Upper level degeneracy
    g_l = 2               # Lower level degeneracy

    # Collision rates [cm^3/s]
    gamma_e = 2.8e-7 * (T / 100.0)**(-0.5)  # electrons
    gamma_H = 8.0e-10 * (T / 100.0)**0.07    # neutral H

    # Total de-excitation rate
    q_tot = gamma_e * n_e + gamma_H * n_HI

    # Critical density
    n_cr = A_ul / (gamma_e * x_e + gamma_H * (1.0 - x_e) + 1e-30)

    # Boltzmann factor
    exp_dE = np.exp(-Delta_E / (k_B * T))

    # Cooling rate per H nucleus [erg/s]
    # Standard two-level formula
    Lambda = n_CII / n_H * Delta_E * A_ul * (g_u / g_l) * exp_dE / \
             (1.0 + n_cr / n_H + (g_u / g_l) * exp_dE)

    return Lambda


def oi_63um_cooling(n_H, T, x_e):
    """OI 63 micron cooling"""
    n_OI = x_O * n_H
    n_HI = (1.0 - x_e) * n_H

    Delta_E = 228.0 * k_B
    A_ul = 8.9e-5
    gamma_H = 9.2e-11 * (T / 100.0)**0.67
    n_cr = A_ul / (gamma_H + 1e-30)
    exp_dE = np.exp(-Delta_E / (k_B * T))

    Lambda = n_OI / n_H * Delta_E * A_ul * 0.6 * exp_dE / (1.0 + n_cr / n_H + 0.6 * exp_dE)

    return Lambda


def grain_cooling(n_H, T):
    """Gas-grain thermal exchange"""
    if T <= T_gr:
        return 0.0
    return 1.2e-31 * n_H * (T / 1000.0)**0.5 * (1.0 - 0.8 * np.exp(-75.0 / T)) * (T - T_gr)


def pe_recombination_cooling(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """PE recombination cooling on grains (BT94)"""
    n_e = x_e * n_H
    if n_e < 1e-6:
        return 0.0

    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)
    if G_eff < 1e-10:
        return 0.0

    psi = G_eff * np.sqrt(T) / n_e
    beta = 0.74 / T**0.068

    Lambda = 4.65e-30 * T**0.94 * psi**beta * n_e / n_H

    return Lambda


# =============================================================================
# NET COOLING FUNCTION
# =============================================================================

def total_heating(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """Total heating rate [erg/s per H]"""
    return (photoelectric_heating(n_H, T, x_e, G_0, N_H) +
            cosmic_ray_heating(n_H, T, x_e) +
            xray_heating(n_H, T, N_H))


def total_cooling(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """Total cooling rate [erg/s per H]"""
    return (lyman_alpha_cooling(n_H, T, x_e) +
            cii_158um_cooling(n_H, T, x_e) +
            oi_63um_cooling(n_H, T, x_e) +
            grain_cooling(n_H, T) +
            pe_recombination_cooling(n_H, T, x_e, G_0, N_H))


def net_cooling(n_H, T, N_H=1e20, G_0=1.7):
    """Lambda - Gamma"""
    x_e = compute_electron_fraction(n_H, T, N_H, G_0)
    return total_cooling(n_H, T, x_e, G_0, N_H) - total_heating(n_H, T, x_e, G_0, N_H)


def find_equilibrium_T(n_H, N_H=1e20, G_0=1.7, T_min=10.0, T_max=2e4):
    """Find equilibrium temperature"""
    def func(log_T):
        return net_cooling(n_H, 10**log_T, N_H, G_0)

    try:
        log_T_eq = brentq(func, np.log10(T_min), np.log10(T_max), xtol=1e-4)
        return 10**log_T_eq
    except:
        return np.nan


# =============================================================================
# DIAGNOSTIC
# =============================================================================

def diagnose(n_H, T, N_H=1e20, G_0=1.7):
    """Print thermal balance at given conditions"""
    x_e = compute_electron_fraction(n_H, T, N_H, G_0)

    print(f"\nn_H = {n_H:.2e}, T = {T:.0f} K, x_e = {x_e:.4e}")

    Gamma_PE = photoelectric_heating(n_H, T, x_e, G_0, N_H)
    Gamma_CR = cosmic_ray_heating(n_H, T, x_e)
    Gamma_X = xray_heating(n_H, T, N_H)
    Gamma_tot = total_heating(n_H, T, x_e, G_0, N_H)

    Lambda_Lya = lyman_alpha_cooling(n_H, T, x_e)
    Lambda_CII = cii_158um_cooling(n_H, T, x_e)
    Lambda_OI = oi_63um_cooling(n_H, T, x_e)
    Lambda_GR = grain_cooling(n_H, T)
    Lambda_PE = pe_recombination_cooling(n_H, T, x_e, G_0, N_H)
    Lambda_tot = total_cooling(n_H, T, x_e, G_0, N_H)

    print(f"  Heating: PE={Gamma_PE:.2e}, CR={Gamma_CR:.2e}, X={Gamma_X:.2e}, Total={Gamma_tot:.2e}")
    print(f"  Cooling: Lya={Lambda_Lya:.2e}, CII={Lambda_CII:.2e}, OI={Lambda_OI:.2e}, GR={Lambda_GR:.2e}, PE_rec={Lambda_PE:.2e}, Total={Lambda_tot:.2e}")
    print(f"  Net (Lambda-Gamma) = {Lambda_tot - Gamma_tot:.2e}, ratio = {Gamma_tot/(Lambda_tot+1e-40):.3f}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*70)
    print("K&I 2000 Thermal Equilibrium - V3 (Fixed electron fraction)")
    print("="*70)

    # Test electron fraction at key densities
    print("\nElectron fractions at various (n, T):")
    for n, T in [(0.1, 8000), (1.0, 5000), (10.0, 1000), (100.0, 200), (1000.0, 50)]:
        x_e = compute_electron_fraction(n, T, N_H=1e20)
        print(f"  n={n:.1f}, T={T:.0f} K: x_e = {x_e:.4e}")

    # Diagnostic at key conditions
    print("\n" + "="*70)
    diagnose(0.1, 8000, N_H=1e20)
    diagnose(1.0, 5000, N_H=1e20)
    diagnose(10.0, 1000, N_H=1e20)
    diagnose(100.0, 200, N_H=1e20)

    # Calculate equilibrium curve
    print("\n" + "="*70)
    print("Calculating equilibrium curve...")
    print("="*70)

    log_n = np.linspace(-1, 6, 150)
    n_arr = 10**log_n

    T_eq_1e19 = np.array([find_equilibrium_T(n, N_H=1e19) for n in n_arr])
    T_eq_1e20 = np.array([find_equilibrium_T(n, N_H=1e20) for n in n_arr])

    # Hardcoded data
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
    ax.loglog(n_data, T_data, 'o', label='Hardcoded (from paper)', color='C2', markersize=3, alpha=0.7)
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'$T_{eq}$ [K]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 2e4)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Thermal Equilibrium Temperature')

    # Electron fraction
    ax = axes[1]
    x_e_calc = [compute_electron_fraction(n, np.interp(np.log10(n), np.log10(n_arr), T_eq_1e20), 1e20) for n in n_arr]
    ax.loglog(n_arr, x_e_calc, '-', label='Calculated x_e', color='C0', linewidth=2)

    # Hardcoded x_e
    x_e_data = np.array([
        9.124377e-02, 7.929356e-02, 6.944798e-02, 6.035237e-02, 5.244801e-02, 4.557889e-02, 3.991953e-02, 3.469127e-02,
        3.014775e-02, 2.579382e-02, 2.241560e-02, 1.902936e-02, 1.590465e-02, 1.288476e-02, 1.003907e-02, 7.761109e-03,
        5.770578e-03, 4.358015e-03, 3.240293e-03, 2.505039e-03, 1.951784e-03, 1.544625e-03, 1.261136e-03, 1.045864e-03,
        8.741278e-04, 7.596432e-04, 6.448868e-04, 5.736924e-04, 5.143534e-04, 4.684013e-04, 4.265546e-04, 3.976417e-04,
        3.706887e-04, 3.537428e-04, 3.375715e-04, 3.323471e-04, 3.246616e-04, 3.171539e-04, 3.171539e-04, 3.098198e-04,
        3.098198e-04, 3.098198e-04, 3.098198e-04, 3.098198e-04, 3.098198e-04, 3.026553e-04, 3.026553e-04, 3.026553e-04,
        3.026553e-04, 3.026553e-04, 2.956565e-04, 2.956565e-04, 2.888196e-04, 2.821407e-04, 2.756163e-04, 2.692427e-04,
        2.569344e-04, 2.451887e-04, 2.285692e-04, 2.081490e-04, 1.910371e-04, 1.699469e-04, 1.476890e-04, 1.263598e-04,
        1.047905e-04
    ])
    ax.loglog(n_data, x_e_data, 'o', label='Hardcoded x_e', color='C2', markersize=3)
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'$x_e$')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-5, 1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Electron Fraction')

    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/scripts/ki2000_equilibrium_v3.png', dpi=150)
    print("\nSaved: ki2000_equilibrium_v3.png")

    # Print comparison
    print("\n" + "="*70)
    print("COMPARISON (N_H = 1e20)")
    print("="*70)
    key_n = [0.1, 1.0, 10.0, 100.0, 1000.0, 1e4, 1e5]
    for n in key_n:
        T_calc = np.interp(np.log10(n), np.log10(n_arr), T_eq_1e20)
        T_hard = np.interp(np.log10(n), np.log10(n_data), T_data)
        if np.isnan(T_calc):
            print(f"n = {n:>8.1e}: T_calc = NaN, T_hardcoded = {T_hard:>8.1f} K")
        else:
            ratio = T_calc / T_hard if T_hard > 0 else 0
            status = "OK" if 0.5 < ratio < 2.0 else "MISMATCH"
            print(f"n = {n:>8.1e}: T_calc = {T_calc:>8.1f} K, T_hardcoded = {T_hard:>8.1f} K, ratio = {ratio:.2f} [{status}]")


if __name__ == "__main__":
    main()
