#!/usr/bin/env python3
"""
Koyama & Inutsuka (2000) Thermal Equilibrium with Full Chemistry Network

Implements the COMPLETE thermal and chemical processes from Appendix A:
- Appendix A.1: Heating and Cooling processes
- Appendix A.2: Chemical reactions (ionization, H2, CO)

Solves the chemical equilibrium numerically to get correct x_e, x_H2, x_CO,
then uses these to calculate thermal equilibrium.

References:
- K&I 2000: Koyama & Inutsuka (2000) ApJ 532, 980
- W95: Wolfire et al. (1995) ApJ 443, 152
- BT94: Bakes & Tielens (1994) ApJ 427, 822
- HM79: Hollenbach & McKee (1979) ApJS 41, 555
- HM89: Hollenbach & McKee (1989) ApJ 342, 306
- TH85: Tielens & Hollenbach (1985) ApJ 291, 722
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, fsolve
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Physical Constants (CGS)
# =============================================================================
k_B = 1.380649e-16      # Boltzmann constant [erg/K]
m_H = 1.6726e-24        # Hydrogen mass [g]
eV = 1.602e-12          # 1 eV in erg
yr = 3.1557e7           # seconds per year

# =============================================================================
# ISM Parameters (from W95 as used by K&I 2000)
# =============================================================================
G_0 = 1.7               # FUV field (1.7 x Habing)
zeta_CR = 1.8e-17       # Cosmic ray ionization rate [s^-1] per H

# Abundances relative to H (W95, listed in K&I paper Appendix A.1.4)
x_He = 0.1              # Helium
x_O = 4.6e-4            # Oxygen
x_C = 3.0e-4            # Carbon
x_Si = 3.55e-6          # Silicon
x_Fe = 7.08e-7          # Iron

# Grain temperature (for 100 Angstrom grains)
T_gr = 8.0              # K


# =============================================================================
# CHEMICAL NETWORK (Appendix A.2)
# =============================================================================

def solve_electron_fraction(n_H, T, N_H=1e20, G_0=1.7):
    """
    Solve for equilibrium electron fraction (Appendix A.2.1, W95 Appendix)

    IMPORTANT: In the neutral ISM, electrons come from:
    1. Hydrogen ionization (by UV, CR, X-ray, collisional)
    2. Metal ions (C+, Si+, Fe+) - already ionized by UV

    The key physics:
    - At LOW COLUMN (N_H < 1e20): C, Si, Fe are photoionized -> x_e ~ x_C ~ 3e-4
    - At HIGH COLUMN: UV blocked, CR ionization dominates -> x_e ~ sqrt(zeta/(n*alpha))

    Ionization equilibrium for H:
    (1-x_e) * Gamma = x_e * x_HII * n_H * alpha
    where x_HII ~ x_e for charge neutrality

    So: (1-x_e) * Gamma = x_e^2 * n_H * alpha
    """
    # Dust optical depth for UV attenuation
    tau_d = N_H / 1.87e21

    # =========================================================================
    # HYDROGEN IONIZATION RATES
    # =========================================================================
    # UV photoionization (strongly attenuated at high column)
    # The 2.5*tau_d factor is for dust extinction of ionizing UV
    Gamma_H_UV = 4.6e-10 * G_0 * np.exp(-2.5 * tau_d)  # [s^-1]

    # Cosmic ray ionization
    Gamma_H_CR = zeta_CR  # [s^-1]

    # X-ray ionization
    tau_X = N_H / 1e22
    Gamma_H_X = 1.2e-17 * np.exp(-tau_X)  # [s^-1]

    # Total H ionization rate
    Gamma_H = Gamma_H_UV + Gamma_H_CR + Gamma_H_X

    # H recombination coefficient (Shapiro & Kang 1987)
    T4 = T / 1e4
    alpha_H = 2.6e-13 * T4**(-0.7)  # [cm^3/s]

    # Collisional ionization coefficient (HM79)
    if T > 3000:
        k_coll = 5.8e-9 * T4**0.5 * np.exp(-15.8 / T4)  # [cm^3/s]
    else:
        k_coll = 0.0

    # =========================================================================
    # SOLVE IONIZATION EQUILIBRIUM
    # =========================================================================
    # Full equation including collisional ionization:
    # (1-x_e) * (Gamma + x_e*n_H*k_coll) = x_e^2 * n_H * alpha
    #
    # Rearranging: x_e^2 * n_H * alpha + x_e^2 * n_H * k_coll - Gamma - x_e*n_H*k_coll = 0
    # This is complex. Use iterative approach.

    # First, estimate x_e from simple balance (ignoring collisional)
    # (1-x_e) * Gamma = x_e^2 * n_H * alpha
    # For x_e << 1: x_e ~ sqrt(Gamma / (n_H * alpha))

    x_e_H = np.sqrt(Gamma_H / (n_H * alpha_H + 1e-30))

    # If collisional ionization is important, iterate
    if k_coll > 0:
        for _ in range(5):
            # Effective ionization rate including collisions
            Gamma_eff = Gamma_H + x_e_H * n_H * k_coll
            # Solve: (1-x_e)*Gamma_eff = x_e^2 * n_H * alpha
            # x_e^2 * n_H * alpha + x_e * Gamma_eff - Gamma_eff = 0
            a = n_H * alpha_H
            b = Gamma_eff
            c = -Gamma_eff
            if a > 0:
                x_e_H = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)
            x_e_H = max(1e-10, min(0.999, x_e_H))

    # =========================================================================
    # METAL ION CONTRIBUTION
    # =========================================================================
    # At low column density, C, Si, Fe are photoionized
    # Each ionized metal contributes one electron

    # UV factor for metal ionization
    UV_factor = np.exp(-tau_d)  # UV attenuation

    # C+ contributes x_C electrons if UV is present
    # Recombination of C+ is on grains, with rate ~ n_e * alpha_gr
    # At equilibrium: Gamma_C_UV * x_C0 = n_e * alpha_gr * x_C+
    # For strong UV: almost all C is C+, so x_C+ ~ x_C

    # Grain recombination becomes important at high n_e
    # alpha_gr ~ 4e-12 * (T/100)^0.5 * (n_e/n_H)^(-0.7) (simplified)
    # At low density: n_e dominated by C+, so n_e ~ n_H * x_C

    x_metals = (x_C + x_Si + x_Fe) * UV_factor

    # Total electron fraction
    # At low density: x_e ~ x_metals (metal ions dominate)
    # At high density: x_e ~ x_e_H (H ionization dominates)

    # Use simple sum (valid when x_e << 1)
    x_e = x_e_H + x_metals

    # At very high density or T, H ionization dominates
    if x_e_H > x_metals:
        x_e = x_e_H  # Avoid double-counting

    # Clamp to physical range
    x_e = np.clip(x_e, 1e-10, 0.999)

    return x_e


def recombination_rate(T):
    """
    Radiative recombination coefficient for hydrogen (Shapiro & Kang 1987)
    alpha_rec = 2.6e-13 * (T/10^4)^-0.7 cm^3/s
    """
    T4 = T / 1e4
    alpha_rec = 2.6e-13 * T4**(-0.7)  # cm^3/s
    return alpha_rec


def h2_formation_rate_grains(n_H, T, x_H):
    """
    H2 formation rate on grain surfaces (TH85, K&I eq. 16-17)

    R = 0.5 * v_H * n_H * n_d * sigma_d * S(T)
      ≈ 6e-17 * (T/300)^0.5 * n_H * S(T) [cm^-3 s^-1]

    where S(T) is sticking coefficient
    """
    # Sticking coefficient (TH85 eq. 17)
    S_T = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)

    # Formation rate coefficient [s^-1] per H atom
    R_gr = 6e-17 * np.sqrt(T / 300.0) * S_T * n_H * x_H

    return R_gr  # [s^-1] formation rate per H nucleus


def h2_associative_detachment_rate(n_H, T, x_e):
    """
    H2 formation via H + H- -> H2 + e- (HM79)

    This is a minor channel compared to grain formation,
    but included for completeness.
    """
    # H- formation rate limited by radiative attachment
    # Simplified: R_assoc ~ few x 10^-20 * x_e * n_H at typical ISM conditions
    R_assoc = 3e-20 * x_e * n_H  # [s^-1] rough estimate
    return R_assoc


def h2_photodissociation_rate(T, N_H2, N_H=1e20, G_0=1.7):
    """
    H2 photodissociation rate (TH85, K&I eq. 19-20)

    R_pump = 3.4e-10 * G_0 * beta(tau) [s^-1]

    where beta accounts for dust attenuation and H2 self-shielding
    """
    # Dust attenuation of UV
    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    if G_eff < 1e-15:
        return 0.0

    # H2 self-shielding (TH85)
    # tau_H2 = 1.2e-14 * N_H2 / delta_V_D
    # Assume delta_V_D ~ sound speed ~ 1 km/s
    delta_V_D = max(0.1, np.sqrt(k_B * T / m_H) / 1e5)  # km/s
    tau_H2 = 1.2e-14 * N_H2 / delta_V_D

    # Self-shielding function beta (TH85 eq. 37)
    if tau_H2 < 2:
        beta = np.exp(-tau_H2)
    elif tau_H2 < 1e3:
        beta = tau_H2**(-0.75)
    else:
        beta = 4e-3 * tau_H2**(-0.25)

    R_photo = 3.4e-10 * G_eff * beta  # [s^-1]

    return R_photo


def h2_cosmic_ray_destruction_rate():
    """
    H2 destruction by cosmic rays (HM89)
    Rate = 2.29 * zeta_CR per H2 molecule
    """
    return 2.29 * zeta_CR  # [s^-1]


def h2_collisional_dissociation_rate(n_H, T, x_H):
    """
    Collisional dissociation of H2 (HM79)
    k_D = 3.4e-9 * exp(8000/T) * exp(-5.19e4/T) cm^3/s
    """
    if T < 1000:
        return 0.0

    k_D = 3.4e-9 * np.exp(8000.0 / T) * np.exp(-5.19e4 / T)  # cm^3/s
    R_coll = k_D * n_H * x_H  # [s^-1]

    return R_coll


def solve_h2_fraction(n_H, T, x_e, N_H=1e20, G_0=1.7):
    """
    Solve for equilibrium H2 fraction (Appendix A.2.2)

    Formation: grain surface + associative detachment
    Destruction: photodissociation + cosmic ray + collisional

    At equilibrium: 2 * x_H * R_form = x_H2 * R_dest
    """
    x_H = 1.0 - x_e  # Neutral H fraction (ignoring H2 for now)

    # Formation rates [s^-1]
    R_form_gr = h2_formation_rate_grains(n_H, T, x_H)
    R_form_assoc = h2_associative_detachment_rate(n_H, T, x_e)
    R_form = R_form_gr + R_form_assoc

    # Estimate N_H2 for self-shielding (iterative would be better)
    # Start with small guess
    x_H2_guess = 1e-6
    N_H2_guess = x_H2_guess * N_H / 2.0

    # Iterate a few times to get self-consistent N_H2
    for _ in range(5):
        R_photo = h2_photodissociation_rate(T, N_H2_guess, N_H, G_0)
        R_CR = h2_cosmic_ray_destruction_rate()
        R_coll = h2_collisional_dissociation_rate(n_H, T, x_H)
        R_dest = R_photo + R_CR + R_coll

        if R_dest > 0:
            # Equilibrium: d(x_H2)/dt = 2*x_H*R_form - x_H2*R_dest = 0
            # x_H2 = 2 * x_H * R_form / R_dest
            x_H2 = 2.0 * x_H * R_form / R_dest
        else:
            x_H2 = 1.0

        x_H2 = min(1.0, max(0.0, x_H2))
        N_H2_guess = x_H2 * N_H / 2.0

    return x_H2


def solve_co_fraction(n_H, T, x_H2, N_H=1e20, G_0=1.7):
    """
    Solve for equilibrium CO fraction (Appendix A.2.3)

    Using simplified Nelson & Langer (1997) treatment:
    dn(CO)/dt = k_0 * n(C+) * n * beta - Gamma_CO * n(CO)
    """
    if x_H2 < 1e-10:
        return 0.0

    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    # Parameters from K&I eq. 21-22
    k_0 = 5e-16       # cm^3 s^-1
    Gamma_CO = 1e-10 * G_eff  # CO photodissociation rate [s^-1]
    Gamma_CHx = 5e-10 * G_eff  # CHx photodissociation
    k_1 = 5e-10       # O + CHx -> CO rate

    n_H2 = x_H2 * n_H / 2.0

    # beta factor (eq. 22)
    if n_H2 > 0 and G_eff > 0:
        beta = k_1 * x_O / (k_1 * x_O + G_eff * Gamma_CHx / n_H2)
    else:
        beta = 1.0

    # Equilibrium: k_0 * x_C * n_H * n_H * beta = Gamma_CO * x_CO * x_C * n_H
    # x_CO = k_0 * n_H * beta / Gamma_CO  (fraction of C in CO)
    if Gamma_CO > 0:
        x_CO = k_0 * n_H * beta / Gamma_CO
    else:
        x_CO = 1.0

    x_CO = min(1.0, max(0.0, x_CO))

    return x_CO


def solve_chemistry(n_H, T, N_H=1e20, G_0=1.7):
    """
    Solve complete chemical equilibrium for x_e, x_H2, x_CO
    """
    # Electron fraction
    x_e = solve_electron_fraction(n_H, T, N_H, G_0)

    # H2 fraction
    x_H2 = solve_h2_fraction(n_H, T, x_e, N_H, G_0)

    # Update neutral H fraction
    x_H = max(0.0, 1.0 - x_e - x_H2)

    # CO fraction (relative to C abundance)
    x_CO = solve_co_fraction(n_H, T, x_H2, N_H, G_0)

    return x_e, x_H, x_H2, x_CO


# =============================================================================
# HEATING PROCESSES (Appendix A.1.1-A.1.3)
# =============================================================================

def photoelectric_heating(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """
    Photoelectric heating from grains and PAHs (BT94, K&I eq. 7-8)
    Gamma_pe = 1.0e-24 * epsilon * G_0  [erg/s per H]
    """
    n_e = max(x_e * n_H, 1e-6)

    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    if G_eff < 1e-15:
        return 0.0

    psi = G_eff * np.sqrt(T) / n_e

    term1 = 4.9e-2 / (1.0 + (psi / 1925.0)**0.73)
    term2 = 3.7e-2 * (T / 1e4)**0.7 / (1.0 + psi / 5000.0)
    epsilon = term1 + term2

    Gamma_pe = 1.0e-24 * epsilon * G_eff

    return Gamma_pe


def cosmic_ray_heating(n_H, T, x_e):
    """
    Cosmic ray heating (K&I Appendix A.1.2)
    Gamma_CR = zeta_CR * E_h(E, x_e)
    """
    # Heat deposited per ionization (Shull & Van Steenberg 1985)
    if x_e < 0.01:
        E_h = 7.0 * eV
    elif x_e < 0.1:
        E_h = 6.0 * eV
    else:
        E_h = 4.0 * eV

    return zeta_CR * E_h


def xray_heating(n_H, T, N_H=1e20):
    """X-ray heating (W95 Appendix)"""
    tau_X = N_H / 1e22
    return 1.5e-25 * np.exp(-tau_X)


def h2_formation_heating(n_H, T, x_H, x_H2):
    """
    H2 formation heating (HM79, K&I eq. 12)
    """
    if x_H < 1e-10:
        return 0.0

    # Formation rate
    S_T = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)
    R = 6e-17 * np.sqrt(T / 300.0) * S_T * n_H

    # Critical density
    denom = 1.6 * x_H * np.exp(-(400.0 / max(T, 10))**2) + \
            1.4 * x_H2 * np.exp(-12000.0 / (T + 1200.0))
    if denom > 0:
        n_cr = 1e6 * T**(-0.5) / denom
    else:
        n_cr = 1e30

    # Energy: 0.2 eV escapes, up to 4.2 eV heats
    E_heat = (0.2 + 4.2 / (1.0 + n_cr / n_H)) * eV

    Gamma_H2 = R * x_H * E_heat

    return Gamma_H2


def h2_photodissociation_heating(n_H, T, x_H, x_H2, N_H=1e20, G_0=1.7):
    """
    UV pumping / photodissociation heating (HM79, K&I eq. 10)
    """
    if x_H2 < 1e-10:
        return 0.0

    N_H2 = x_H2 * N_H / 2.0
    R_pump = h2_photodissociation_rate(T, N_H2, N_H, G_0)

    denom = 1.6 * x_H * np.exp(-(400.0 / max(T, 10))**2)
    if denom > 0:
        n_cr = 1e6 * T**(-0.5) / denom
    else:
        n_cr = 1e30

    E_heat = 2.2 * eV / (1.0 + n_cr / n_H)

    Gamma_UV = 9.0 * R_pump * E_heat * x_H2

    return Gamma_UV


# =============================================================================
# COOLING PROCESSES (Appendix A.1.4-A.1.7)
# =============================================================================

def lyman_alpha_cooling(n_H, T, x_e, x_H):
    """Lyman-alpha cooling (HM89)"""
    if T < 5000:
        return 0.0

    n_e = x_e * n_H
    n_HI = x_H * n_H

    Lambda = 7.5e-19 * np.exp(-118348.0 / T) * n_e * n_HI / n_H

    return Lambda


def cii_158um_cooling(n_H, T, x_e):
    """
    CII 158 micron fine-structure cooling (HM89)
    Main coolant at intermediate temperatures
    """
    n_e = x_e * n_H
    n_HI = (1.0 - x_e) * n_H
    n_CII = x_C * n_H

    Delta_E = 91.2 * k_B  # 91.2 K
    A_ul = 2.4e-6
    g_ratio = 2.0  # g_u/g_l = 4/2

    gamma_e = 2.8e-7 * (T / 100.0)**(-0.5)
    gamma_H = 8.0e-10 * (T / 100.0)**0.07

    n_cr = A_ul / (gamma_e * x_e + gamma_H * (1.0 - x_e) + 1e-30)
    exp_dE = np.exp(-Delta_E / (k_B * T))

    Lambda = n_CII / n_H * Delta_E * A_ul * g_ratio * exp_dE / \
             (1.0 + n_cr / n_H + g_ratio * exp_dE)

    return Lambda


def oi_63um_cooling(n_H, T, x_e):
    """OI 63 micron fine-structure cooling"""
    n_OI = x_O * n_H

    Delta_E = 228.0 * k_B
    A_ul = 8.9e-5
    gamma_H = 9.2e-11 * (T / 100.0)**0.67

    n_cr = A_ul / (gamma_H + 1e-30)
    exp_dE = np.exp(-Delta_E / (k_B * T))

    Lambda = n_OI / n_H * Delta_E * A_ul * 0.6 * exp_dE / (1.0 + n_cr / n_H)

    return Lambda


def h2_cooling(n_H, T, x_H, x_H2):
    """
    H2 rovibrational cooling (HM79, K&I Appendix A.1.5)
    """
    if x_H2 < 1e-10 or T < 100:
        return 0.0

    n_H2 = x_H2 * n_H / 2.0

    # LTE cooling function (HM79)
    logT = np.log10(T)
    if T < 100:
        L_H_LTE = 0.0
    elif T < 1000:
        L_H_LTE = 10**(-24.216 + 4.595 * (logT - 3) - 3.458 * (logT - 3)**2)
    else:
        L_H_LTE = 10**(-24.216 + 4.595 * (logT - 3) - 3.458 * (logT - 3)**2)

    # Critical density
    n_cr_H = 10**(3.0 + 4.0 * (logT - 3) - 0.8 * (logT - 3)**2) if T > 100 else 1e10

    Lambda = n_H2 * x_H * L_H_LTE / (1.0 + n_cr_H / n_H) / n_H

    return Lambda


def co_cooling(n_H, T, x_CO):
    """CO rotational cooling (McKee et al. 1982, K&I eq. 15)"""
    if x_CO < 1e-10:
        return 0.0

    n_CO = x_CO * x_C * n_H

    A_0 = 9.7e-8
    E_0 = 2.76 * k_B
    T3 = T / 1000.0
    n_cr = 3.3e6 * T3**0.75

    Lambda_rot = 4.0 * (k_B * T)**2 * A_0 / \
                 (E_0 * (1.0 + n_cr / n_H + 1.5 * np.sqrt(n_cr / n_H)))

    return Lambda_rot * n_CO / n_H


def grain_cooling(n_H, T):
    """Gas-grain thermal exchange (HM89, K&I eq. 18)"""
    if T <= T_gr:
        return 0.0

    Lambda = 1.2e-31 * n_H * (T / 1000.0)**0.5 * \
             (1.0 - 0.8 * np.exp(-75.0 / T)) * (T - T_gr)

    return Lambda


def pe_recombination_cooling(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """PE recombination cooling (BT94 eq. 9)"""
    n_e = x_e * n_H
    if n_e < 1e-6:
        return 0.0

    tau_d = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_d)

    if G_eff < 1e-15:
        return 0.0

    psi = G_eff * np.sqrt(T) / n_e
    beta = 0.74 / T**0.068

    Lambda = 4.65e-30 * T**0.94 * psi**beta * n_e / n_H

    return Lambda


# =============================================================================
# TOTAL RATES AND EQUILIBRIUM
# =============================================================================

def total_heating(n_H, T, x_e, x_H, x_H2, G_0=1.7, N_H=1e20):
    """Total heating rate [erg/s per H]"""
    Gamma = 0.0
    Gamma += photoelectric_heating(n_H, T, x_e, G_0, N_H)
    Gamma += cosmic_ray_heating(n_H, T, x_e)
    Gamma += xray_heating(n_H, T, N_H)
    Gamma += h2_formation_heating(n_H, T, x_H, x_H2)
    Gamma += h2_photodissociation_heating(n_H, T, x_H, x_H2, N_H, G_0)
    return Gamma


def total_cooling(n_H, T, x_e, x_H, x_H2, x_CO, G_0=1.7, N_H=1e20):
    """Total cooling rate [erg/s per H]"""
    Lambda = 0.0
    Lambda += lyman_alpha_cooling(n_H, T, x_e, x_H)
    Lambda += cii_158um_cooling(n_H, T, x_e)
    Lambda += oi_63um_cooling(n_H, T, x_e)
    Lambda += h2_cooling(n_H, T, x_H, x_H2)
    Lambda += co_cooling(n_H, T, x_CO)
    Lambda += grain_cooling(n_H, T)
    Lambda += pe_recombination_cooling(n_H, T, x_e, G_0, N_H)
    return Lambda


def net_cooling(n_H, T, N_H=1e20, G_0=1.7):
    """Net cooling rate with self-consistent chemistry"""
    x_e, x_H, x_H2, x_CO = solve_chemistry(n_H, T, N_H, G_0)

    Gamma = total_heating(n_H, T, x_e, x_H, x_H2, G_0, N_H)
    Lambda = total_cooling(n_H, T, x_e, x_H, x_H2, x_CO, G_0, N_H)

    return Lambda - Gamma


def find_equilibrium_temperature(n_H, N_H=1e20, G_0=1.7, T_min=10.0, T_max=1.5e4):
    """Find equilibrium temperature"""
    def func(log_T):
        return net_cooling(n_H, 10**log_T, N_H, G_0)

    try:
        log_T_eq = brentq(func, np.log10(T_min), np.log10(T_max), xtol=1e-4)
        return 10**log_T_eq
    except:
        # Check endpoints
        net_min = net_cooling(n_H, T_min, N_H, G_0)
        net_max = net_cooling(n_H, T_max, N_H, G_0)
        if abs(net_min) < abs(net_max):
            return T_min
        else:
            return T_max


# =============================================================================
# TIMESCALES
# =============================================================================

def cooling_timescale(n_H, T, x_e, x_H, x_H2, x_CO, N_H=1e20, G_0=1.7):
    """Cooling timescale t_cool = (3/2)*n*k*T / (n^2 * Lambda)"""
    Lambda = total_cooling(n_H, T, x_e, x_H, x_H2, x_CO, G_0, N_H)
    if Lambda > 0:
        t_cool = 1.5 * k_B * T / (n_H * Lambda)
        return t_cool / yr  # years
    return 1e15


def recombination_timescale(n_H, T, x_e):
    """Recombination timescale t_rec = 1 / (x_e * n_H * alpha)"""
    alpha = recombination_rate(T)
    if x_e * n_H * alpha > 0:
        t_rec = 1.0 / (x_e * n_H * alpha)
        return t_rec / yr
    return 1e15


def h2_formation_timescale(n_H, T, x_H):
    """H2 formation timescale t_H2 = 1 / (R * n_H * x_H)"""
    S_T = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)
    R = 6e-17 * np.sqrt(T / 300.0) * S_T
    if R * n_H * x_H > 0:
        t_H2 = 1.0 / (R * n_H * x_H)
        return t_H2 / yr
    return 1e15


def freefall_timescale(n_H):
    """Free-fall timescale t_ff = sqrt(3*pi / (32*G*rho))"""
    G = 6.674e-8
    rho = 1.4 * m_H * n_H
    t_ff = np.sqrt(3.0 * np.pi / (32.0 * G * rho))
    return t_ff / yr


# =============================================================================
# MAIN: Generate Figure 1
# =============================================================================

def main():
    print("="*70)
    print("Koyama & Inutsuka (2000) - Full Thermal Equilibrium Calculation")
    print("="*70)

    # Density grid
    log_n = np.linspace(-1, 6, 150)
    n_arr = 10**log_n

    # Results storage
    results = {}

    for N_H_label, N_H in [('1e19', 1e19), ('1e20', 1e20)]:
        print(f"\nCalculating for N_H = {N_H_label} cm^-2...")

        T_eq = []
        P_eq = []
        x_e_eq = []
        x_H2_eq = []
        x_CO_eq = []

        Gamma_PE = []
        Gamma_CR = []
        Gamma_H2 = []
        Lambda_CII = []
        Lambda_OI = []
        Lambda_Lya = []
        Lambda_CO_arr = []
        Lambda_GR = []

        t_cool_arr = []
        t_rec_arr = []
        t_H2_arr = []
        t_ff_arr = []

        for n_H in n_arr:
            # Find equilibrium T
            T = find_equilibrium_temperature(n_H, N_H, G_0)
            T_eq.append(T)

            # Get chemistry at equilibrium
            x_e, x_H, x_H2, x_CO = solve_chemistry(n_H, T, N_H, G_0)

            # Pressure: P/k = (1.1 + x_e - x_H2/2) * n * T
            mu_eff = 1.1 + x_e - x_H2 / 2.0
            P = mu_eff * n_H * T

            P_eq.append(P)
            x_e_eq.append(x_e)
            x_H2_eq.append(x_H2)
            x_CO_eq.append(x_CO)

            # Heating rates
            Gamma_PE.append(photoelectric_heating(n_H, T, x_e, G_0, N_H))
            Gamma_CR.append(cosmic_ray_heating(n_H, T, x_e))
            Gamma_H2.append(h2_formation_heating(n_H, T, x_H, x_H2))

            # Cooling rates
            Lambda_CII.append(cii_158um_cooling(n_H, T, x_e))
            Lambda_OI.append(oi_63um_cooling(n_H, T, x_e))
            Lambda_Lya.append(lyman_alpha_cooling(n_H, T, x_e, x_H))
            Lambda_CO_arr.append(co_cooling(n_H, T, x_CO))
            Lambda_GR.append(grain_cooling(n_H, T))

            # Timescales
            t_cool_arr.append(cooling_timescale(n_H, T, x_e, x_H, x_H2, x_CO, N_H, G_0))
            t_rec_arr.append(recombination_timescale(n_H, T, x_e))
            t_H2_arr.append(h2_formation_timescale(n_H, T, x_H))
            t_ff_arr.append(freefall_timescale(n_H))

        results[N_H_label] = {
            'n': n_arr,
            'T': np.array(T_eq),
            'P': np.array(P_eq),
            'x_e': np.array(x_e_eq),
            'x_H2': np.array(x_H2_eq),
            'x_CO': np.array(x_CO_eq),
            'Gamma_PE': np.array(Gamma_PE),
            'Gamma_CR': np.array(Gamma_CR),
            'Gamma_H2': np.array(Gamma_H2),
            'Lambda_CII': np.array(Lambda_CII),
            'Lambda_OI': np.array(Lambda_OI),
            'Lambda_Lya': np.array(Lambda_Lya),
            'Lambda_CO': np.array(Lambda_CO_arr),
            'Lambda_GR': np.array(Lambda_GR),
            't_cool': np.array(t_cool_arr),
            't_rec': np.array(t_rec_arr),
            't_H2': np.array(t_H2_arr),
            't_ff': np.array(t_ff_arr),
        }

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
    x_H2_data = np.array([
        1.000000e-08, 1.098104e-08, 1.442737e-08, 1.910371e-08, 2.451887e-08, 3.171539e-08, 4.070548e-08, 5.024592e-08,
        5.872729e-08, 6.154060e-08, 6.154060e-08, 6.154060e-08, 6.917770e-08, 8.741278e-08, 1.095965e-07, 1.439927e-07,
        1.862559e-07, 2.447111e-07, 3.165362e-07, 4.158790e-07, 5.506778e-07, 7.235044e-07, 9.580137e-07, 1.258680e-06,
        1.692855e-06, 2.294623e-06, 3.159196e-06, 4.557889e-06, 6.944798e-06, 1.171078e-05, 6.023482e-05, 2.001887e-04,
        4.575697e-04, 9.304020e-04, 1.777413e-03, 3.316997e-03, 6.047015e-03, 1.093831e-02, 1.902936e-02, 3.310537e-02,
        5.626155e-02, 9.340370e-02, 1.479772e-01, 2.185459e-01, 3.080130e-01, 4.078492e-01, 4.995288e-01, 5.884190e-01,
        6.770987e-01, 7.435249e-01, 7.913912e-01, 8.489338e-01, 8.690299e-01, 9.106605e-01, 9.322177e-01, 9.322177e-01,
        9.542853e-01, 9.768753e-01, 9.768753e-01, 9.768753e-01, 9.768753e-01, 9.768753e-01, 9.768753e-01, 9.768753e-01,
        9.768753e-01
    ])

    # =========================================================================
    # Create Figure 1 (4 panels like K&I 2000)
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel (a): Temperature and Pressure vs density
    ax = axes[0, 0]
    r = results['1e20']
    ax.loglog(r['n'], r['T'], '-', label='T (calculated)', color='C0', linewidth=2)
    ax.loglog(n_data, T_data, 'o', label='T (K&I data)', color='C0', markersize=3, alpha=0.5)
    ax.loglog(r['n'], r['P'], '-', label='P/k (calculated)', color='C1', linewidth=2)
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'T [K], P/k [K cm$^{-3}$]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 1e7)
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(a) Equilibrium T and P (N_H=1e20)')
    ax.grid(True, alpha=0.3)

    # Panel (b): Chemical fractions
    ax = axes[0, 1]
    ax.loglog(r['n'], r['x_e'], '-', label=r'$x_e$ (calc)', color='C0', linewidth=2)
    ax.loglog(n_data, x_e_data, 'o', label=r'$x_e$ (K&I)', color='C0', markersize=3, alpha=0.5)
    ax.loglog(r['n'], r['x_H2'], '-', label=r'$x_{H_2}$ (calc)', color='C1', linewidth=2)
    ax.loglog(n_data, x_H2_data, 'o', label=r'$x_{H_2}$ (K&I)', color='C1', markersize=3, alpha=0.5)
    ax.loglog(r['n'], np.maximum(r['x_CO'], 1e-12), '-', label=r'$x_{CO}$ (calc)', color='C2', linewidth=2)
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel('Fraction')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-10, 1)
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(b) Chemical fractions')
    ax.grid(True, alpha=0.3)

    # Panel (c): Heating and Cooling rates
    ax = axes[1, 0]
    ax.loglog(r['n'], r['Gamma_PE'], '--', label='PE (heat)', color='C0')
    ax.loglog(r['n'], r['Gamma_CR'], '--', label='CR (heat)', color='C1')
    ax.loglog(r['n'], np.maximum(r['Gamma_H2'], 1e-30), '--', label='H2 (heat)', color='C2')
    ax.loglog(r['n'], r['Lambda_CII'], '-', label='CII', color='C3')
    ax.loglog(r['n'], r['Lambda_OI'], '-', label='OI', color='C4')
    ax.loglog(r['n'], np.maximum(r['Lambda_Lya'], 1e-30), '-', label='Lya', color='C5')
    ax.loglog(r['n'], np.maximum(r['Lambda_CO'], 1e-30), '-', label='CO', color='C6')
    ax.loglog(r['n'], np.maximum(r['Lambda_GR'], 1e-30), '-', label='GR', color='C7')
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'Rate [erg s$^{-1}$ per H]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-28, 1e-22)
    ax.legend(loc='best', fontsize=7, ncol=2)
    ax.set_title('(c) Heating (--) and Cooling (-) rates')
    ax.grid(True, alpha=0.3)

    # Panel (d): Timescales
    ax = axes[1, 1]
    ax.loglog(r['n'], r['t_cool'], '-', label=r'$t_{cool}$', color='C0')
    ax.loglog(r['n'], r['t_rec'], '--', label=r'$t_{rec}$', color='C1')
    ax.loglog(r['n'], r['t_H2'], '-.', label=r'$t_{H_2}$', color='C2')
    ax.loglog(r['n'], r['t_ff'], ':', label=r'$t_{ff}$', color='C3')
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel('Timescale [yr]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1, 1e12)
    ax.legend(loc='best')
    ax.set_title('(d) Timescales')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/scripts/ki2000_figure1_full.png', dpi=150)
    print("\nSaved: ki2000_figure1_full.png")

    # Print comparison table
    print("\n" + "="*70)
    print("COMPARISON WITH K&I DATA (N_H = 1e20)")
    print("="*70)
    key_n = [0.1, 1.0, 10.0, 100.0, 1000.0, 1e4, 1e5]
    print(f"{'n [cm^-3]':>12s} {'T_calc':>10s} {'T_KI':>10s} {'x_e_calc':>12s} {'x_e_KI':>12s}")
    print("-"*60)
    for n in key_n:
        T_c = np.interp(np.log10(n), np.log10(r['n']), r['T'])
        T_k = np.interp(np.log10(n), np.log10(n_data), T_data)
        x_e_c = np.interp(np.log10(n), np.log10(r['n']), r['x_e'])
        x_e_k = np.interp(np.log10(n), np.log10(n_data), x_e_data)
        print(f"{n:>12.1e} {T_c:>10.1f} {T_k:>10.1f} {x_e_c:>12.4e} {x_e_k:>12.4e}")


if __name__ == "__main__":
    main()
