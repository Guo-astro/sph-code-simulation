#!/usr/bin/env python3
"""
Koyama & Inutsuka (2000) Thermal Equilibrium Calculator

This script implements ALL thermal processes from Appendix A of K&I (2000)
"Molecular Cloud Formation in Shock-Compressed Layers"
ApJ, 532, 980 (2000)

Reproduces Figure 1: Equilibrium temperature, pressure, chemical fractions,
heating/cooling rates, and timescales vs density.
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, fsolve
from scipy.integrate import quad

# =============================================================================
# Physical Constants (CGS)
# =============================================================================
k_B = 1.380649e-16      # Boltzmann constant [erg/K]
m_H = 1.6726e-24        # Hydrogen mass [g]
m_e = 9.109e-28         # Electron mass [g]
e = 4.803e-10           # Electron charge [esu]
c = 2.998e10            # Speed of light [cm/s]
h = 6.626e-27           # Planck constant [erg s]
eV = 1.602e-12          # 1 eV in erg

# =============================================================================
# ISM Parameters (Following Wolfire et al. 1995 as used by K&I 2000)
# =============================================================================
G_0 = 1.7               # FUV field (1.7 x Habing)
zeta_CR = 1.8e-17       # Cosmic ray ionization rate [s^-1] (eq. in paper)

# Abundances relative to H (from W95, listed in paper)
x_He = 0.1              # Helium
x_O = 4.6e-4            # Oxygen
x_C = 3.0e-4            # Carbon
x_Si = 3.55e-6          # Silicon
x_Fe = 7.08e-7          # Iron

# Grain parameters
T_gr = 8.0              # Grain temperature [K] (100 Angstrom grains)
a_grain = 100e-8        # Grain radius [cm] = 100 Angstrom

# =============================================================================
# Heating Processes (Appendix A.1)
# =============================================================================

def photoelectric_heating(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """
    Photoelectric heating from small grains and PAHs
    Bakes & Tielens (1994) - Eq. 7-9 in paper

    Gamma_pe = 1.0e-24 * epsilon * G_0  [erg/s per H]

    epsilon depends on G_0 * T^0.5 / n_e
    """
    n_e = x_e * n_H
    if n_e <= 0:
        n_e = 1e-10 * n_H  # Minimum electron density

    # FUV attenuation by dust (simple exponential)
    # tau_UV ~ N_H / (1.87e21)
    tau_UV = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_UV)

    if G_eff <= 0:
        return 0.0

    # Heating efficiency (BT94 fit, eq. 8)
    psi = G_eff * np.sqrt(T) / n_e

    # Two-term fit from BT94
    term1 = 4.9e-2 / (1.0 + (psi / 1925.0)**0.73)
    term2 = 3.7e-2 * (T / 1e4)**0.7 / (1.0 + psi / 5000.0)
    epsilon = term1 + term2

    Gamma_pe = 1.0e-24 * epsilon * G_eff  # [erg/s per H]
    return Gamma_pe


def photoelectric_cooling(n_H, T, x_e, G_0=1.7, N_H=1e20):
    """
    Recombination cooling on grains (BT94 eq. 9)
    Lambda_pe = 4.65e-30 * T^0.94 * (G_0*T^0.5/n_e)^beta * n_e
    """
    n_e = x_e * n_H
    if n_e <= 0:
        n_e = 1e-10 * n_H

    tau_UV = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_UV)

    if G_eff <= 0 or T <= 0:
        return 0.0

    psi = G_eff * np.sqrt(T) / n_e
    beta = 0.74 / T**0.068

    Lambda_pe = 4.65e-30 * T**0.94 * psi**beta * n_e / n_H  # per H
    return Lambda_pe


def cosmic_ray_heating(n_H, T, x_e, x_H2=0):
    """
    Cosmic ray heating (Appendix A.1.2)
    Primary + secondary ionization heating

    Gamma_CR = zeta_CR * E_h(E, x_e)

    E_h from Shull & Van Steenberg (1985)
    """
    # For primary CR ionization of H, typical energy ~ 20-40 eV
    # Heat deposited per ionization: E_h ~ 6-10 eV for x_e ~ 0.01-0.1
    # SVS85 approximation: E_h ~ 6.5 eV * (1 - x_e^0.4) for x_e < 0.1

    if x_e < 0.01:
        E_h = 7.0 * eV  # ~7 eV per ionization for low ionization
    elif x_e < 0.1:
        E_h = 6.5 * eV * (1.0 - x_e**0.4)
    else:
        E_h = 2.0 * eV  # Most energy goes to ionization at high x_e

    Gamma_CR = zeta_CR * E_h  # [erg/s per H]
    return Gamma_CR


def xray_heating(n_H, T, x_e, N_H=1e20):
    """
    X-ray heating (Appendix A.1.2)
    Following Wolfire et al. 1995, using their analytic fits
    """
    # X-ray flux is attenuated by column density
    # Simple fit from W95: heating rate ~ 1e-25 * exp(-N_H/10^22)
    # Strongly attenuated for N_H > 10^20

    tau_X = N_H / 1e22
    Gamma_X = 1e-25 * np.exp(-tau_X)  # [erg/s per H]

    return Gamma_X


def h2_formation_heating(n_H, T, x_e, x_H, x_H2, N_H=1e20):
    """
    H2 formation heating on grains (Appendix A.1.3)
    From Hollenbach & McKee (1979), eq. 12

    Gamma_gr = R * x_H * {0.2 + 4.2 * [1 + n_cr(H2)/n]^-1}

    where R is the H2 formation rate on grains
    """
    # H2 formation rate on grains (TH85, eq. 16-17)
    S_T = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)
    R = 6e-17 * np.sqrt(T / 300.0) * n_H * S_T  # Formation rate coefficient

    # Critical density for H2 (HM79)
    n_cr_H2 = 1e6 * T**(-0.5) / (1.6 * x_H * np.exp(-(400.0 / T)**2) +
                                  1.4 * x_H2 * np.exp(-12000.0 / (T + 1200.0)))
    if n_cr_H2 <= 0:
        n_cr_H2 = 1e30

    # Heating rate per H nucleus (eq. 12)
    # Energy: 0.2 eV escapes, up to 4.2 eV heats gas
    energy_factor = 0.2 + 4.2 / (1.0 + n_cr_H2 / n_H)
    Gamma_H2_form = R * x_H * energy_factor * eV

    return Gamma_H2_form


def h2_photodissociation_heating(n_H, T, x_H2, N_H2=1e15, G_0=1.7, N_H=1e20):
    """
    UV pumping and photodissociation heating (Appendix A.1.3, eq. 10-11)
    From Hollenbach & McKee (1979)
    """
    # FUV attenuation
    tau_UV = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_UV)

    if G_eff <= 0 or x_H2 <= 0:
        return 0.0

    # H2 self-shielding function beta(tau) from TH85
    # tau = 1.2e-14 * N_H2 / delta_V_D
    # Assume delta_V_D ~ 1 km/s = 1e5 cm/s
    delta_V_D = 1.0  # km/s
    tau_H2 = 1.2e-14 * N_H2 / delta_V_D

    # beta(tau) approximation from TH85
    if tau_H2 < 1:
        beta = 1.0
    elif tau_H2 < 1e4:
        beta = tau_H2**(-0.75)
    else:
        beta = 1e-4

    # Photodissociation rate (eq. 19)
    R_pump = 3.4e-10 * G_eff * beta  # [s^-1]

    # Critical density for UV pumping
    n_cr_H2 = 1e6 * T**(-0.5) / max(0.01, 1.6 * (1 - x_H2) * np.exp(-(400.0 / T)**2))

    # Heating per pumping event ~ 2.2 eV (eq. 10)
    energy_per_pump = 2.2 * eV

    Gamma_UV = 9.0 * R_pump * energy_per_pump / (1.0 + n_cr_H2 / n_H) * x_H2

    return Gamma_UV


# =============================================================================
# Cooling Processes (Appendix A.1.4-A.1.7)
# =============================================================================

def lyman_alpha_cooling(n_H, T, x_e, x_H):
    """
    Hydrogen Lyman-alpha cooling (Appendix A.1.4)
    Collisional excitation of H by electrons
    """
    if T < 3000:
        return 0.0  # Negligible below ~3000 K

    n_e = x_e * n_H
    n_HI = x_H * n_H

    # Lya excitation energy = 10.2 eV
    E_Lya = 10.2 * eV

    # Collisional excitation rate coefficient (HM89)
    # k_Lya ~ 1.9e-8 * (1/sqrt(T)) * exp(-E/kT)
    T4 = T / 1e4
    k_Lya = 7.5e-19 * np.exp(-118348.0 / T)  # [cm^3/s]

    # Cooling rate per H
    Lambda_Lya = k_Lya * n_e * n_HI * E_Lya / n_H

    return Lambda_Lya


def cii_fine_structure_cooling(n_H, T, x_e, x_C=x_C):
    """
    CII 158 micron fine structure cooling (Appendix A.1.4)
    Main coolant for WNM/CNM
    """
    n_e = x_e * n_H
    n_HI = (1.0 - x_e) * n_H  # Neutral H
    n_CII = x_C * n_H  # Assume all C is CII (ionization potential < 13.6 eV)

    # Energy of 158 um transition
    E_CII = 92.0 * k_B  # ~92 K
    Delta_E = E_CII

    # Statistical weights
    g_l = 2  # Lower level (2P_1/2)
    g_u = 4  # Upper level (2P_3/2)

    # Collision rates with electrons and H (from HM89)
    # Electron collisions (dominant at higher T)
    gamma_e = 2.8e-7 * (T / 100.0)**(-0.5)  # [cm^3/s]

    # H collisions
    gamma_H = 8.0e-10 * (T / 100.0)**0.07  # [cm^3/s]

    # Critical density
    A_ul = 2.3e-6  # [s^-1] spontaneous emission
    n_cr = A_ul / (gamma_e * x_e + gamma_H * (1 - x_e))

    # Cooling rate (2-level atom)
    exp_factor = np.exp(-Delta_E / (k_B * T))

    # Effective rate
    q_eff = (gamma_e * n_e + gamma_H * n_HI) / n_H

    # Cooling function per H nucleus
    Lambda_CII = n_CII / n_H * Delta_E * A_ul * exp_factor / (1.0 + n_cr / n_H + (g_u / g_l) * exp_factor)

    return Lambda_CII


def oi_fine_structure_cooling(n_H, T, x_e, x_O=x_O):
    """
    OI 63 micron and 145 micron fine structure cooling (Appendix A.1.4)
    Important at T ~ 100-500 K
    """
    n_HI = (1.0 - x_e) * n_H
    n_OI = x_O * n_H  # O stays neutral (same ionization potential as H)

    # 63 micron transition (dominant)
    E_63 = 228.0 * k_B  # ~228 K
    A_63 = 8.9e-5  # [s^-1]

    # 145 micron transition
    E_145 = 99.0 * k_B  # ~99 K
    A_145 = 1.7e-5  # [s^-1]

    # Collision rate with H (HM89)
    gamma_H_63 = 9.2e-11 * (T / 100.0)**0.67  # [cm^3/s]
    gamma_H_145 = 4.3e-11 * (T / 100.0)**0.80

    # 63 micron cooling
    exp_63 = np.exp(-E_63 / (k_B * T))
    n_cr_63 = A_63 / gamma_H_63
    Lambda_63 = n_OI / n_H * E_63 * A_63 * exp_63 / (1.0 + n_cr_63 / n_H)

    # 145 micron cooling (usually smaller)
    exp_145 = np.exp(-E_145 / (k_B * T))
    n_cr_145 = A_145 / gamma_H_145
    Lambda_145 = n_OI / n_H * E_145 * A_145 * exp_145 / (1.0 + n_cr_145 / n_H)

    return Lambda_63 + Lambda_145


def h2_cooling(n_H, T, x_H, x_H2):
    """
    H2 rovibrational cooling (Appendix A.1.5)
    From Hollenbach & McKee (1979)
    """
    if x_H2 <= 0 or T < 100:
        return 0.0

    n_H2 = x_H2 * n_H / 2.0  # H2 number density
    n_HI = x_H * n_H

    # LTE cooling function from HM79
    # L_vr(LTE) for collisions with H and H2

    # Rotation-vibration cooling with H collisions
    # Low-T fit (T < 1000 K)
    if T < 100:
        L_H_LTE = 0.0
    elif T < 1000:
        L_H_LTE = 10**(-24.216 + 4.595 * np.log10(T/1000) - 3.458 * (np.log10(T/1000))**2)
    else:
        L_H_LTE = 10**(-24.216 + 4.595 * np.log10(T/1000) - 3.458 * (np.log10(T/1000))**2)

    # Critical density for H collisions (HM79)
    n_cr_H = 10**(3.0 + 4.0 * np.log10(T/1000.0) - 0.8 * (np.log10(T/1000.0))**2)

    # H2-H2 collisions (usually subdominant)
    L_H2_LTE = L_H_LTE * 0.1  # Approximate
    n_cr_H2 = n_cr_H * 10.0

    # Total cooling
    Lambda_H2 = n_H2 * (x_H * L_H_LTE / (1.0 + n_cr_H / n_H) +
                        x_H2 * L_H2_LTE / (1.0 + n_cr_H2 / n_H)) / n_H

    return Lambda_H2


def co_cooling(n_H, T, x_CO, x_C=x_C):
    """
    CO rotational and vibrational cooling (Appendix A.1.6)
    From McKee et al. (1982)
    """
    if x_CO <= 0:
        return 0.0

    n_CO = x_CO * x_C * n_H  # CO number density

    # Rotational cooling (eq. 15)
    A_0 = 9.7e-8  # [s^-1]
    E_0 = 2.76 * k_B  # ~2.76 K
    T3 = T / 1000.0
    n_cr_rot = 3.3e6 * T3**0.75  # [cm^-3]

    Lambda_rot = 4.0 * (k_B * T)**2 * A_0 / (E_0 * (1.0 + n_cr_rot / n_H + 1.5 * np.sqrt(n_cr_rot / n_H)))

    # Vibrational cooling (eq. 16, usually subdominant at low T)
    Delta_E_vib = 3080.0 * k_B  # ~3080 K

    # Vibrational rate coefficients from HM89
    gamma_H_vib = 1e-12 * np.exp(-3080.0 / T)
    gamma_H2_vib = 1e-12 * np.exp(-3080.0 / T) * 10.0

    Lambda_vib = n_CO / n_H * Delta_E_vib * (gamma_H_vib * (1.0 - x_CO) + gamma_H2_vib * x_CO * 0.5)

    return (Lambda_rot + Lambda_vib) * n_CO / n_H


def grain_cooling(n_H, T, T_gr=8.0):
    """
    Gas-grain energy exchange (Appendix A.1.7)
    From Hollenbach & McKee (1989)
    """
    if T < T_gr:
        return 0.0  # Grains heat gas

    # Cooling rate per H (eq. 18)
    a_100 = a_grain / 100e-8  # Normalized to 100 Angstrom

    Lambda_gr = 1.2e-31 * n_H * (T / 1000.0)**0.5 * a_100**(-0.5) * \
                (1.0 - 0.8 * np.exp(-75.0 / T)) * (T - T_gr)

    return Lambda_gr


# =============================================================================
# Chemical Equilibrium (Appendix A.2)
# =============================================================================

def electron_equilibrium(n_H, T, N_H=1e20, G_0=1.7):
    """
    Equilibrium electron fraction from ionization balance

    Ionization: photoionization + cosmic ray + X-ray + collisional
    Recombination: radiative + grain recombination
    """
    # Photoionization rate (attenuated by dust)
    tau_UV = N_H / 1.87e21
    Gamma_photo = 4.6e-10 * G_0 * np.exp(-2.5 * tau_UV)  # [s^-1]

    # Cosmic ray ionization
    Gamma_CR = zeta_CR  # [s^-1]

    # X-ray ionization (from W95)
    tau_X = N_H / 1e22
    Gamma_X = 1e-15 * np.exp(-tau_X)  # [s^-1]

    # Collisional ionization (HM79)
    T4 = T / 1e4
    k_coll = 5.8e-9 * T4**0.5 * np.exp(-15.8 / T4) if T > 3000 else 0.0  # [cm^3/s]

    # Radiative recombination (Shapiro & Kang 1987)
    alpha_rec = 2.6e-13 * (T / 1e4)**(-0.7)  # [cm^3/s]

    # Total ionization rate
    Gamma_ion = Gamma_photo + Gamma_CR + Gamma_X

    # Equilibrium: n_H * Gamma + n_e * n_H * k_coll = n_e * n_H * alpha_rec
    # x_e * (Gamma + x_e * n_H * k_coll) = x_e^2 * n_H * alpha_rec
    # Gamma + x_e * n_H * k_coll = x_e * n_H * alpha_rec
    # x_e = Gamma / (n_H * alpha_rec - n_H * k_coll)

    # Solve quadratic for x_e
    a = n_H * (alpha_rec - k_coll)
    if a > 0:
        x_e = Gamma_ion / a
    else:
        # High T limit: collisional ionization dominates
        x_e = np.sqrt(Gamma_ion / (n_H * alpha_rec)) if alpha_rec > 0 else 1.0

    # Clamp to physical range
    x_e = max(1e-10, min(1.0, x_e))

    return x_e


def h2_equilibrium(n_H, T, x_e, N_H=1e20, N_H2=1e15, G_0=1.7):
    """
    Equilibrium H2 fraction (Appendix A.2.2)

    Formation: grain surface + associative detachment
    Destruction: photodissociation + cosmic ray + collisional
    """
    x_H = 1.0 - x_e  # Neutral fraction (ignoring H2)

    # Formation rate on grains (TH85, eq. 16-17)
    S_T = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)
    R_gr = 6e-17 * np.sqrt(T / 300.0) * S_T * n_H  # [s^-1] per H

    # Associative detachment: H + H- -> H2 + e-
    # H- formation/destruction equilibrium (simplified)
    R_assoc = 1e-20 * x_e * n_H  # [s^-1] approximate

    # Total formation rate
    R_form = R_gr + R_assoc

    # Photodissociation (TH85, eq. 19-20)
    tau_UV = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_UV)

    # Self-shielding
    delta_V_D = 1.0  # km/s
    tau_H2 = 1.2e-14 * N_H2 / delta_V_D
    if tau_H2 < 1:
        beta = 1.0
    elif tau_H2 < 1e4:
        beta = tau_H2**(-0.75)
    else:
        beta = 1e-4

    R_photo = 3.4e-10 * G_eff * beta  # [s^-1]

    # Cosmic ray destruction
    R_CR = 2.29 * zeta_CR  # [s^-1]

    # Collisional dissociation (HM79)
    k_diss = 3.4e-9 * np.exp(8000.0 / T) * np.exp(-5.19e4 / T) if T > 2000 else 0.0
    R_coll = k_diss * n_H * x_H  # [s^-1]

    # Total destruction rate
    R_dest = R_photo + R_CR + R_coll

    # Equilibrium: 2 * x_H * R_form = x_H2 * R_dest
    # x_H2 = 2 * x_H * R_form / R_dest
    if R_dest > 0:
        x_H2 = 2.0 * x_H * R_form / (R_dest + 1e-30)
    else:
        x_H2 = 1.0

    # Clamp
    x_H2 = max(0.0, min(1.0, x_H2))

    return x_H2


def co_equilibrium(n_H, x_H2, N_H=1e20, G_0=1.7):
    """
    Equilibrium CO fraction (Appendix A.2.3)
    Simplified Nelson & Langer (1997) treatment
    """
    if x_H2 < 1e-8:
        return 0.0

    tau_UV = N_H / 1.87e21
    G_eff = G_0 * np.exp(-tau_UV)

    # Formation: C+ + O -> CO (through intermediates)
    k_0 = 5e-16  # [cm^3/s]
    Gamma_CO = 1e-10 * G_eff  # [s^-1] photodissociation
    Gamma_CHx = 5e-10 * G_eff  # CHx photodissociation

    k_1 = 5e-10  # O + CHx -> CO reaction rate

    n_H2 = x_H2 * n_H / 2.0

    # beta factor (eq. 22)
    beta = k_1 * x_O / (k_1 * x_O + G_eff * Gamma_CHx / max(n_H2, 1.0))

    # Equilibrium CO fraction
    # Formation: k_0 * n(C+) * n * beta
    # Destruction: Gamma_CO * n(CO)
    if Gamma_CO > 0:
        x_CO = k_0 * n_H * beta / (Gamma_CO + 1e-30)
    else:
        x_CO = 1.0

    # Clamp to carbon abundance
    x_CO = max(0.0, min(1.0, x_CO))

    return x_CO


# =============================================================================
# Total Heating and Cooling
# =============================================================================

def total_heating(n_H, T, x_e, x_H, x_H2, N_H=1e20, G_0=1.7):
    """Total heating rate per H nucleus [erg/s]"""
    Gamma = 0.0

    Gamma += photoelectric_heating(n_H, T, x_e, G_0, N_H)
    Gamma += cosmic_ray_heating(n_H, T, x_e, x_H2)
    Gamma += xray_heating(n_H, T, x_e, N_H)
    Gamma += h2_formation_heating(n_H, T, x_e, x_H, x_H2, N_H)

    # H2 photodissociation heating (small)
    N_H2 = x_H2 * N_H / 2.0
    Gamma += h2_photodissociation_heating(n_H, T, x_H2, N_H2, G_0, N_H)

    return Gamma


def total_cooling(n_H, T, x_e, x_H, x_H2, x_CO, N_H=1e20):
    """Total cooling rate per H nucleus [erg/s]"""
    Lambda = 0.0

    Lambda += lyman_alpha_cooling(n_H, T, x_e, x_H)
    Lambda += cii_fine_structure_cooling(n_H, T, x_e)
    Lambda += oi_fine_structure_cooling(n_H, T, x_e)
    Lambda += h2_cooling(n_H, T, x_H, x_H2)
    Lambda += co_cooling(n_H, T, x_CO)
    Lambda += grain_cooling(n_H, T)

    # PE recombination cooling
    Lambda += photoelectric_cooling(n_H, T, x_e, G_0, N_H)

    return Lambda


def net_cooling(n_H, T, N_H=1e20, G_0=1.7):
    """
    Net cooling rate Lambda - Gamma, with self-consistent chemistry
    """
    # Chemical equilibrium
    x_e = electron_equilibrium(n_H, T, N_H, G_0)
    x_H = 1.0 - x_e
    x_H2 = h2_equilibrium(n_H, T, x_e, N_H, N_H * x_H * 1e-5, G_0)

    # Update neutral H fraction accounting for H2
    x_H = 1.0 - x_e - x_H2
    if x_H < 0:
        x_H = 0.0
        x_H2 = 1.0 - x_e

    x_CO = co_equilibrium(n_H, x_H2, N_H, G_0)

    Gamma = total_heating(n_H, T, x_e, x_H, x_H2, N_H, G_0)
    Lambda = total_cooling(n_H, T, x_e, x_H, x_H2, x_CO, N_H)

    return Lambda - Gamma


# =============================================================================
# Thermal Equilibrium Solver
# =============================================================================

def find_equilibrium_temperature(n_H, N_H=1e20, G_0=1.7, T_min=10.0, T_max=1e5):
    """
    Find equilibrium temperature where Gamma = Lambda
    """
    def func(log_T):
        T = 10**log_T
        return net_cooling(n_H, T, N_H, G_0)

    try:
        log_T_eq = brentq(func, np.log10(T_min), np.log10(T_max))
        return 10**log_T_eq
    except ValueError:
        # No crossing found - check endpoints
        net_low = net_cooling(n_H, T_min, N_H, G_0)
        net_high = net_cooling(n_H, T_max, N_H, G_0)

        if net_low < 0:
            return T_min
        elif net_high > 0:
            return T_max
        else:
            return np.nan


# =============================================================================
# Main: Generate Figure 1
# =============================================================================

def generate_figure1():
    """Generate all panels of K&I 2000 Figure 1"""

    # Density range
    log_n_arr = np.linspace(-1, 6, 200)
    n_arr = 10**log_n_arr

    # Storage for both column densities
    results = {}

    for N_H_label, N_H in [('1e19', 1e19), ('1e20', 1e20)]:
        print(f"\nCalculating equilibrium for N_H = {N_H_label} cm^-2...")

        T_eq = np.zeros_like(n_arr)
        P_eq = np.zeros_like(n_arr)
        x_e_eq = np.zeros_like(n_arr)
        x_H2_eq = np.zeros_like(n_arr)
        x_CO_eq = np.zeros_like(n_arr)

        # Heating/cooling rates at equilibrium
        Gamma_PE = np.zeros_like(n_arr)
        Gamma_CR = np.zeros_like(n_arr)
        Gamma_H2 = np.zeros_like(n_arr)
        Lambda_CII = np.zeros_like(n_arr)
        Lambda_OI = np.zeros_like(n_arr)
        Lambda_Lya = np.zeros_like(n_arr)
        Lambda_CO = np.zeros_like(n_arr)
        Lambda_GR = np.zeros_like(n_arr)

        # Timescales
        t_cool = np.zeros_like(n_arr)
        t_rec = np.zeros_like(n_arr)
        t_H2 = np.zeros_like(n_arr)
        t_ff = np.zeros_like(n_arr)

        for i, n_H in enumerate(n_arr):
            # Find equilibrium temperature
            T = find_equilibrium_temperature(n_H, N_H, G_0)
            if np.isnan(T):
                T = 100.0  # Default
            T_eq[i] = T

            # Pressure P/k = n * T (ignoring He contribution for now)
            # Paper uses P = (1.1 + x_e - x_H2/2) * n * k * T
            x_e = electron_equilibrium(n_H, T, N_H, G_0)
            x_H = 1.0 - x_e
            x_H2 = h2_equilibrium(n_H, T, x_e, N_H)
            x_H = max(0.0, 1.0 - x_e - x_H2)
            x_CO = co_equilibrium(n_H, x_H2, N_H, G_0)

            mu_eff = 1.1 + x_e - x_H2 / 2.0
            P_eq[i] = mu_eff * n_H * T  # P/k [K cm^-3]

            x_e_eq[i] = x_e
            x_H2_eq[i] = x_H2
            x_CO_eq[i] = x_CO

            # Individual rates at equilibrium
            Gamma_PE[i] = photoelectric_heating(n_H, T, x_e, G_0, N_H)
            Gamma_CR[i] = cosmic_ray_heating(n_H, T, x_e, x_H2)
            Gamma_H2[i] = h2_formation_heating(n_H, T, x_e, x_H, x_H2, N_H)

            Lambda_CII[i] = cii_fine_structure_cooling(n_H, T, x_e)
            Lambda_OI[i] = oi_fine_structure_cooling(n_H, T, x_e)
            Lambda_Lya[i] = lyman_alpha_cooling(n_H, T, x_e, x_H)
            Lambda_CO[i] = co_cooling(n_H, T, x_CO)
            Lambda_GR[i] = grain_cooling(n_H, T)

            # Timescales
            # Cooling time: tau_cool = (3/2) * n * k * T / (n^2 * Lambda)
            Lambda_total = Lambda_CII[i] + Lambda_OI[i] + Lambda_Lya[i] + Lambda_CO[i] + Lambda_GR[i]
            if Lambda_total > 0 and n_H > 0:
                t_cool[i] = 1.5 * k_B * T / (n_H * Lambda_total) / (365.25 * 24 * 3600)  # years
            else:
                t_cool[i] = 1e12

            # Recombination time
            alpha_rec = 2.6e-13 * (T / 1e4)**(-0.7)
            if x_e * n_H * alpha_rec > 0:
                t_rec[i] = 1.0 / (x_e * n_H * alpha_rec) / (365.25 * 24 * 3600)
            else:
                t_rec[i] = 1e12

            # H2 formation time
            S_T = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)
            R_H2 = 6e-17 * np.sqrt(T / 300.0) * S_T * n_H
            if R_H2 > 0:
                t_H2[i] = 1.0 / R_H2 / (365.25 * 24 * 3600)
            else:
                t_H2[i] = 1e12

            # Free-fall time
            G = 6.674e-8
            rho = 1.4 * m_H * n_H
            t_ff[i] = np.sqrt(3.0 * np.pi / (32.0 * G * rho)) / (365.25 * 24 * 3600)

        results[N_H_label] = {
            'n': n_arr,
            'T': T_eq,
            'P': P_eq,
            'x_e': x_e_eq,
            'x_H2': x_H2_eq,
            'x_CO': x_CO_eq,
            'Gamma_PE': Gamma_PE,
            'Gamma_CR': Gamma_CR,
            'Gamma_H2': Gamma_H2,
            'Lambda_CII': Lambda_CII,
            'Lambda_OI': Lambda_OI,
            'Lambda_Lya': Lambda_Lya,
            'Lambda_CO': Lambda_CO,
            'Lambda_GR': Lambda_GR,
            't_cool': t_cool,
            't_rec': t_rec,
            't_H2': t_H2,
            't_ff': t_ff,
        }

    return results


def plot_figure1(results, save_path=None):
    """Plot Figure 1 from K&I 2000"""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel (a): Temperature and Pressure
    ax = axes[0, 0]
    for N_H_label, style in [('1e19', '-'), ('1e20', '--')]:
        r = results[N_H_label]
        ax.loglog(r['n'], r['T'], style, label=f'T (N_H={N_H_label})', color='C0')
        ax.loglog(r['n'], r['P'], style, label=f'P/k (N_H={N_H_label})', color='C1')
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'T [K], P/k [K cm$^{-3}$]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 1e7)
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(a) Equilibrium T and P')
    ax.grid(True, alpha=0.3)

    # Panel (b): Chemical fractions
    ax = axes[0, 1]
    r = results['1e20']
    ax.loglog(r['n'], r['x_e'], '-', label=r'$x_e$', color='C0')
    ax.loglog(r['n'], r['x_H2'], '-', label=r'$x_{H_2}$', color='C1')
    ax.loglog(r['n'], r['x_CO'], '-', label=r'$x_{CO}$', color='C2')
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel('Fraction')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-10, 1)
    ax.legend(loc='best')
    ax.set_title('(b) Chemical fractions')
    ax.grid(True, alpha=0.3)

    # Panel (c): Heating and Cooling rates
    ax = axes[1, 0]
    r = results['1e20']
    # Heating (dashed)
    ax.loglog(r['n'], r['Gamma_PE'], '--', label='PE', color='C0')
    ax.loglog(r['n'], r['Gamma_CR'], '--', label='CR', color='C1')
    ax.loglog(r['n'], r['Gamma_H2'], '--', label=r'$H_2$', color='C2')
    # Cooling (solid)
    ax.loglog(r['n'], r['Lambda_CII'], '-', label='CII', color='C3')
    ax.loglog(r['n'], r['Lambda_OI'], '-', label='OI', color='C4')
    ax.loglog(r['n'], r['Lambda_Lya'], '-', label=r'Ly$\alpha$', color='C5')
    ax.loglog(r['n'], r['Lambda_CO'], '-', label='CO', color='C6')
    ax.loglog(r['n'], r['Lambda_GR'], '-', label='GR', color='C7')
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'$\Gamma$, $\Lambda$ [erg s$^{-1}$ per H]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-28, 1e-22)
    ax.legend(loc='best', fontsize=8, ncol=2)
    ax.set_title('(c) Heating (--) and Cooling (-)  rates')
    ax.grid(True, alpha=0.3)

    # Panel (d): Timescales
    ax = axes[1, 1]
    r = results['1e20']
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

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved figure to {save_path}")

    plt.show()
    return fig


def compare_with_hardcoded_data(results):
    """Compare calculated results with hardcoded data from koyama_inutsuka_data.hpp"""

    print("\n" + "="*70)
    print("COMPARISON WITH HARDCODED DATA")
    print("="*70)

    # Hardcoded temperature data for N_H = 1e20 (from the hpp file)
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

    # Interpolate calculated results to comparison points
    r = results['1e20']
    T_calc = np.interp(n_data, r['n'], r['T'])

    # Print comparison table
    print("\nTemperature comparison (N_H = 1e20):")
    print("-" * 60)
    print(f"{'n_H [cm^-3]':>15s} {'T_hardcoded [K]':>18s} {'T_calculated [K]':>18s} {'Ratio':>10s}")
    print("-" * 60)

    for i in range(0, len(n_data), 5):  # Sample every 5th point
        ratio = T_calc[i] / T_data[i] if T_data[i] > 0 else np.nan
        print(f"{n_data[i]:>15.3e} {T_data[i]:>18.2f} {T_calc[i]:>18.2f} {ratio:>10.3f}")

    # Key densities check
    print("\n" + "="*70)
    print("KEY DENSITY CHECKS:")
    print("="*70)

    key_densities = [0.1, 1.0, 10.0, 100.0, 1000.0, 1e4, 1e5, 1e6]

    for n in key_densities:
        T_h = np.interp(n, n_data, T_data)
        T_c = np.interp(n, r['n'], r['T'])
        ratio = T_c / T_h if T_h > 0 else np.nan
        status = "OK" if 0.5 < ratio < 2.0 else "MISMATCH"
        print(f"n = {n:>8.1e} cm^-3: T_hardcoded = {T_h:>8.1f} K, T_calculated = {T_c:>8.1f} K, ratio = {ratio:.2f} [{status}]")

    return n_data, T_data, T_calc


if __name__ == "__main__":
    print("="*70)
    print("Koyama & Inutsuka (2000) Thermal Equilibrium Calculator")
    print("="*70)

    # Generate equilibrium curves
    results = generate_figure1()

    # Compare with hardcoded data
    n_data, T_data, T_calc = compare_with_hardcoded_data(results)

    # Plot
    fig = plot_figure1(results, save_path='/Users/guo/Downloads/sphcode/scripts/ki2000_figure1_calculated.png')

    # Also plot comparison
    fig2, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(n_data, T_data, 'o-', label='Hardcoded (from paper)', markersize=3)
    ax.loglog(results['1e20']['n'], results['1e20']['T'], '-', label='Calculated', linewidth=2)
    ax.set_xlabel(r'$n_H$ [cm$^{-3}$]')
    ax.set_ylabel(r'$T_{eq}$ [K]')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 1e5)
    ax.legend()
    ax.set_title('Comparison: Hardcoded vs Calculated Equilibrium Temperature')
    ax.grid(True, alpha=0.3)
    plt.savefig('/Users/guo/Downloads/sphcode/scripts/ki2000_comparison.png', dpi=150)
    plt.show()

    print("\nDone!")
