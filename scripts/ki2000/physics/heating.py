#!/usr/bin/env python3
"""
First-principles heating rates for ISM thermal equilibrium.

Implements heating processes from Koyama & Inutsuka (2000, ApJ 532, 980),
based on formulas from:
- Bakes & Tielens (1994) - Photoelectric heating
- Wolfire et al. (1995, 2003) - Updated photoelectric rates
- Shull & Van Steenberg (1985) - Cosmic ray heating
- Hollenbach & McKee (1979) - H2 formation heating

Usage:
    from ki2000.physics.heating import (
        photoelectric_heating,
        cosmic_ray_heating,
        xray_heating,
        h2_formation_heating
    )
    
    # Heating rate per hydrogen nucleus [erg/s]
    gamma = photoelectric_heating(n, T, x_e, G0=1.7, N_H=1e19)
"""

import numpy as np
from typing import Union

# Physical constants
K_B = 1.380649e-16  # Boltzmann constant [erg/K]


def photoelectric_efficiency(
    T: Union[float, np.ndarray],
    n_e: Union[float, np.ndarray],
    G0: float = 1.7,
    phi_PAH: float = 0.5
) -> Union[float, np.ndarray]:
    """
    Calculate photoelectric heating efficiency.
    
    From Bakes & Tielens (1994), also used in K&I 2000 and Wolfire et al. (2003).
    
    Parameters
    ----------
    T : float or array
        Gas temperature [K]
    n_e : float or array
        Electron number density [cm^-3]
    G0 : float
        FUV field strength in Habing units (local ISM = 1.0)
        K&I 2000 uses G0 = 1.7
    phi_PAH : float
        PAH abundance scaling factor (default 0.5 from Wolfire 2003)
    
    Returns
    -------
    float or array
        Heating efficiency epsilon (dimensionless, 0 < epsilon < 0.05)
    
    Notes
    -----
    The charging parameter is psi = G0 * T^0.5 / (n_e * phi_PAH)
    For high psi (strong UV, low density), efficiency is low.
    For low psi (weak UV, high density), efficiency approaches ~5%.
    """
    # Charging parameter
    psi = G0 * np.sqrt(T) / (n_e * phi_PAH + 1e-30)
    
    # Bakes & Tielens (1994) formula - K&I 2000 uses slightly different constants
    # K&I 2000: 1925 and 5000 instead of Wolfire's 4e-3/1e-24 normalization
    epsilon = (
        4.9e-2 / (1.0 + (psi / 1925.0)**0.73)
        + 3.7e-2 * (T / 1e4)**0.7 / (1.0 + psi / 5000.0)
    )
    
    return epsilon


def photoelectric_heating(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_e: Union[float, np.ndarray],
    G0: float = 1.7,
    N_H: float = 1e19,
    phi_PAH: float = 0.5
) -> Union[float, np.ndarray]:
    """
    Photoelectric heating rate from small grains and PAHs.
    
    From Bakes & Tielens (1994), as used in K&I 2000 eq. A1-A2.
    This is the dominant heating process in the WNM and CNM.
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    T : float or array
        Gas temperature [K]
    x_e : float or array
        Electron fraction n_e/n
    G0 : float
        FUV field strength in Habing units (local ISM = 1.0)
        K&I 2000 uses G0 = 1.7
    N_H : float
        Shielding column density [cm^-2]
        Affects FUV attenuation
    phi_PAH : float
        PAH abundance scaling (Wolfire 2003 uses 0.5)
    
    Returns
    -------
    float or array
        Heating rate per H nucleus [erg/s]
    
    Notes
    -----
    Gamma_pe = 1.0e-24 * epsilon * G0 * exp(-tau_FUV)
    where tau_FUV is the dust optical depth in FUV.
    """
    n_e = n * x_e
    
    epsilon = photoelectric_efficiency(T, n_e, G0, phi_PAH)
    
    # FUV attenuation by dust (approximate)
    # tau_FUV ~ sigma_d * N_H with sigma_d ~ 1e-21 cm^2 per H
    tau_FUV = 1.0e-21 * N_H
    atten = np.exp(-tau_FUV)
    
    # Heating rate per H nucleus [erg/s]
    # K&I 2000 eq. A1: Gamma_pe = 1.0e-24 * epsilon * G0
    gamma_pe = 1.0e-24 * epsilon * G0 * atten
    
    return gamma_pe


def photoelectric_cooling(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_e: Union[float, np.ndarray],
    G0: float = 1.7,
    phi_PAH: float = 0.5
) -> Union[float, np.ndarray]:
    """
    Recombination cooling on grains and PAHs.
    
    From Bakes & Tielens (1994), as used in K&I 2000.
    This is a cooling term that partially offsets photoelectric heating.
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    T : float or array
        Gas temperature [K]
    x_e : float or array
        Electron fraction n_e/n
    G0 : float
        FUV field strength in Habing units
    phi_PAH : float
        PAH abundance scaling
    
    Returns
    -------
    float or array
        Cooling rate per H nucleus [erg/s]
    """
    n_e = n * x_e
    
    # Charging parameter
    psi = G0 * np.sqrt(T) / (n_e * phi_PAH + 1e-30)
    
    # K&I 2000 eq. A3
    beta = 0.74 / T**0.068
    lambda_pe = 4.65e-30 * T**0.94 * psi**beta * n_e
    
    return lambda_pe


def cosmic_ray_heating(
    n: Union[float, np.ndarray],
    x_e: Union[float, np.ndarray],
    zeta_CR: float = 1.8e-17
) -> Union[float, np.ndarray]:
    """
    Cosmic ray ionization heating.
    
    From Shull & Van Steenberg (1985), as used in K&I 2000.
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    x_e : float or array
        Electron fraction n_e/n
    zeta_CR : float
        Primary cosmic ray ionization rate [s^-1]
        K&I 2000 uses zeta_CR = 1.8e-17 s^-1
    
    Returns
    -------
    float or array
        Heating rate per H nucleus [erg/s]
    
    Notes
    -----
    The heat deposited depends on electron fraction because secondary
    ionizations become more efficient at low ionization.
    E_h ~ 6.5 eV for high x_e, up to ~20 eV for low x_e.
    """
    # Heat deposited per ionization from Shull & Van Steenberg (1985)
    # E_h depends on x_e - approximate fit
    x_e = np.maximum(x_e, 1e-10)  # Avoid division by zero
    
    # Approximate E_h(x_e) from SVS85 - interpolate between limits
    # Low ionization: more energy goes to heating
    # High ionization: more energy goes to ionization
    E_h_eV = 6.5 + 13.5 * (1.0 - x_e) / (1.0 + x_e)  # Approximate
    E_h = E_h_eV * 1.602e-12  # Convert eV to erg
    
    # Heating rate per H nucleus
    # Note: zeta_CR is per H atom, need to account for H+ and H2
    gamma_CR = zeta_CR * E_h
    
    return gamma_CR


def xray_heating(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_e: Union[float, np.ndarray],
    N_H: float = 1e19
) -> Union[float, np.ndarray]:
    """
    Soft X-ray heating.
    
    From Wolfire et al. (1995), as used in K&I 2000.
    X-rays ionize and heat via photoelectron thermalization.
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    T : float or array
        Gas temperature [K]
    x_e : float or array
        Electron fraction
    N_H : float
        Absorbing column density [cm^-2]
    
    Returns
    -------
    float or array
        Heating rate per H nucleus [erg/s]
    
    Notes
    -----
    X-ray heating is subdominant to photoelectric heating in most
    diffuse ISM conditions but becomes important at high column
    densities where FUV is attenuated.
    """
    # Wolfire et al. (1995) X-ray heating fit
    # Attenuated by neutral hydrogen column
    # sigma_X ~ 1.5e-22 cm^2 at 0.5 keV
    tau_X = 1.5e-22 * N_H
    atten = np.exp(-tau_X)
    
    # Local diffuse X-ray intensity heating rate
    # From W95 Appendix - approximately 2e-26 erg/s per H at low column
    gamma_X = 2.0e-26 * atten
    
    return gamma_X


def h2_critical_density(T: Union[float, np.ndarray], x_H: Union[float, np.ndarray],
                        x_H2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Critical density for H2 rovibrational cooling.
    
    From Hollenbach & McKee (1979).
    
    Parameters
    ----------
    T : float or array
        Temperature [K]
    x_H : float or array
        Atomic hydrogen fraction n_H/n
    x_H2 : float or array
        Molecular hydrogen fraction n_H2/n
    
    Returns
    -------
    float or array
        Critical density [cm^-3]
    """
    # HM79 critical density formula (K&I 2000 eq. A4)
    T_safe = np.maximum(T, 10.0)
    
    term1 = 1.6 * x_H * np.exp(-(400.0 / T_safe)**2)
    term2 = 1.4 * x_H2 * np.exp(-12000.0 / (T_safe + 1200.0))
    
    n_cr = 1e6 * T_safe**(-0.5) / (term1 + term2 + 1e-30)
    
    return n_cr


def h2_formation_heating(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_H: Union[float, np.ndarray],
    x_H2: Union[float, np.ndarray],
    T_gr: float = 8.0
) -> Union[float, np.ndarray]:
    """
    H2 formation heating on dust grains.
    
    From Hollenbach & McKee (1979), as used in K&I 2000.
    When H2 forms on grains, ~4.5 eV is released, part of which heats the gas.
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    T : float or array
        Gas temperature [K]
    x_H : float or array
        Atomic hydrogen fraction n_H/n (typically 1 - 2*x_H2 for neutral gas)
    x_H2 : float or array
        Molecular hydrogen fraction n_H2/n
    T_gr : float
        Grain temperature [K]
    
    Returns
    -------
    float or array
        Heating rate per H nucleus [erg/s]
    
    Notes
    -----
    The heating efficiency depends on whether excited H2 can radiate
    before collisional de-excitation (depends on density vs n_cr).
    """
    T_safe = np.maximum(T, 10.0)
    
    # H2 formation rate coefficient on grains (TH85)
    # R = 6e-17 * (T/300)^0.5 * S(T) cm^3/s
    sticking = h2_sticking_coefficient(T_safe, T_gr)
    R = 6.0e-17 * (T_safe / 300.0)**0.5 * sticking
    
    # Critical density
    n_cr = h2_critical_density(T_safe, x_H, x_H2)
    
    # Heating rate (HM79, K&I 2000 eq. A6)
    # Energy available: 0.2 eV direct + 4.2 eV if collisionally de-excited
    E_direct = 0.2 * 1.602e-12  # erg
    E_collisional = 4.2 * 1.602e-12  # erg
    
    # Fraction going to heating vs radiation depends on n/n_cr
    heat_per_formation = E_direct + E_collisional / (1.0 + n_cr / (n + 1e-30))
    
    # Formation rate per H nucleus = R * n_H
    gamma_gr = R * x_H * n * heat_per_formation
    
    return gamma_gr


def h2_sticking_coefficient(T: Union[float, np.ndarray],
                            T_gr: float = 8.0) -> Union[float, np.ndarray]:
    """
    Sticking coefficient for H atoms on grains.
    
    From Tielens & Hollenbach (1985).
    
    Parameters
    ----------
    T : float or array
        Gas temperature [K]
    T_gr : float
        Grain temperature [K]
    
    Returns
    -------
    float or array
        Sticking coefficient S(T) (0 < S < 1)
    """
    # TH85 sticking coefficient (K&I 2000 eq. after A5)
    T_eff = T + T_gr
    S = 1.0 / (1.0 + 0.04 * np.sqrt(T_eff) + 2e-3 * T + 8e-6 * T**2)
    return S


def h2_photodissociation_heating(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_H: Union[float, np.ndarray],
    x_H2: Union[float, np.ndarray],
    G0: float = 1.7,
    N_H2: float = 1e14
) -> Union[float, np.ndarray]:
    """
    H2 photodissociation heating.
    
    From Hollenbach & McKee (1979), as used in K&I 2000.
    FUV photons pump H2 to excited states, some of which dissociate.
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    T : float or array
        Gas temperature [K]
    x_H : float or array
        Atomic hydrogen fraction
    x_H2 : float or array
        Molecular hydrogen fraction
    G0 : float
        FUV field strength in Habing units
    N_H2 : float
        H2 column density for self-shielding [cm^-2]
    
    Returns
    -------
    float or array
        Heating rate per H nucleus [erg/s]
    
    Notes
    -----
    ~2.2 eV is deposited as kinetic energy per photodissociation event.
    The rate depends on self-shielding of H2.
    """
    T_safe = np.maximum(T, 10.0)
    
    # H2 photodissociation rate (self-shielded)
    # R_pump = 3.3e-11 * G0 * f_shield [s^-1]
    f_shield = h2_self_shielding(N_H2)
    R_pump = 3.3e-11 * G0 * f_shield
    
    # Critical density for thermalization
    n_cr = h2_critical_density(T_safe, x_H, x_H2)
    
    # Energy deposited: 2.2 eV thermalized fraction depends on density
    E_heat = 2.2 * 1.602e-12  # erg
    heat_fraction = 1.0 / (1.0 + n_cr / (n + 1e-30))
    
    # Heating rate per H nucleus (K&I 2000 eq. A4)
    # Factor of 9 accounts for multiple pumping events per dissociation
    gamma_UV = 9.0 * R_pump * x_H2 * E_heat * heat_fraction
    
    return gamma_UV


def h2_self_shielding(N_H2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    H2 self-shielding factor.
    
    From Draine & Bertoldi (1996).
    
    Parameters
    ----------
    N_H2 : float or array
        H2 column density [cm^-2]
    
    Returns
    -------
    float or array
        Self-shielding factor (0 < f < 1)
    """
    # Draine & Bertoldi (1996) approximate fit
    x = N_H2 / 5e14
    f_shield = np.where(
        N_H2 < 1e14,
        1.0,
        0.965 / (1.0 + x)**2 + 0.035 * np.exp(-8.5e-4 * np.sqrt(1.0 + x))
    )
    return f_shield


def total_heating_rate(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_e: Union[float, np.ndarray],
    x_H: Union[float, np.ndarray] = None,
    x_H2: Union[float, np.ndarray] = 1e-6,
    G0: float = 1.7,
    N_H: float = 1e19,
    zeta_CR: float = 1.8e-17
) -> Union[float, np.ndarray]:
    """
    Total heating rate per hydrogen nucleus.
    
    Sum of all heating processes:
    - Photoelectric (PE)
    - Cosmic ray (CR)
    - X-ray (XR)
    - H2 formation/photodissociation (H2)
    
    Parameters
    ----------
    n : float or array
        Hydrogen nucleus density [cm^-3]
    T : float or array
        Gas temperature [K]
    x_e : float or array
        Electron fraction n_e/n
    x_H : float or array, optional
        Atomic hydrogen fraction (default: 1 - 2*x_H2)
    x_H2 : float or array
        Molecular hydrogen fraction (default: 1e-6)
    G0 : float
        FUV field strength in Habing units
    N_H : float
        Shielding column density [cm^-2]
    zeta_CR : float
        Cosmic ray ionization rate [s^-1]
    
    Returns
    -------
    float or array
        Total heating rate per H nucleus [erg/s]
    """
    if x_H is None:
        x_H = 1.0 - 2.0 * x_H2
    
    # Sum all heating contributions
    gamma_pe = photoelectric_heating(n, T, x_e, G0, N_H)
    gamma_cr = cosmic_ray_heating(n, x_e, zeta_CR)
    gamma_xr = xray_heating(n, T, x_e, N_H)
    gamma_h2 = h2_formation_heating(n, T, x_H, x_H2)
    gamma_uv = h2_photodissociation_heating(n, T, x_H, x_H2, G0)
    
    # Subtract grain recombination cooling from PE heating
    lambda_pe = photoelectric_cooling(n, T, x_e, G0)
    
    total = gamma_pe - lambda_pe + gamma_cr + gamma_xr + gamma_h2 + gamma_uv
    
    return total


def heating_breakdown(
    n: Union[float, np.ndarray],
    T: Union[float, np.ndarray],
    x_e: Union[float, np.ndarray],
    x_H: Union[float, np.ndarray] = None,
    x_H2: Union[float, np.ndarray] = 1e-6,
    G0: float = 1.7,
    N_H: float = 1e19,
    zeta_CR: float = 1.8e-17
) -> dict:
    """
    Get individual heating rate contributions.
    
    Returns a dictionary with all heating processes for diagnostic purposes.
    
    Parameters
    ----------
    Same as total_heating_rate
    
    Returns
    -------
    dict
        Dictionary with keys 'PE', 'CR', 'XR', 'H2_form', 'H2_photo', 'PE_cool', 'total'
    """
    if x_H is None:
        x_H = 1.0 - 2.0 * x_H2
    
    gamma_pe = photoelectric_heating(n, T, x_e, G0, N_H)
    gamma_cr = cosmic_ray_heating(n, x_e, zeta_CR)
    gamma_xr = xray_heating(n, T, x_e, N_H)
    gamma_h2_form = h2_formation_heating(n, T, x_H, x_H2)
    gamma_h2_photo = h2_photodissociation_heating(n, T, x_H, x_H2, G0)
    lambda_pe = photoelectric_cooling(n, T, x_e, G0)
    
    total = gamma_pe - lambda_pe + gamma_cr + gamma_xr + gamma_h2_form + gamma_h2_photo
    
    return {
        'PE': gamma_pe,
        'PE_cool': lambda_pe,
        'CR': gamma_cr,
        'XR': gamma_xr,
        'H2_form': gamma_h2_form,
        'H2_photo': gamma_h2_photo,
        'total': total
    }


if __name__ == '__main__':
    # Quick test
    import matplotlib.pyplot as plt
    
    log_n = np.linspace(-2, 6, 100)
    n = 10**log_n
    
    # Approximate equilibrium electron fraction (simple fit)
    x_e = 1e-4 * (n / 1.0)**(-0.5)
    x_e = np.clip(x_e, 1e-7, 0.1)
    
    # Approximate equilibrium temperature (very rough)
    T = 8000 * (n / 0.1)**(-0.4)
    T = np.clip(T, 10, 1e4)
    
    for N_H in [1e19, 1e20]:
        rates = heating_breakdown(n, T, x_e, N_H=N_H)
        
        plt.figure(figsize=(10, 6))
        for key, label in [('PE', 'Photoelectric'), ('CR', 'Cosmic Ray'),
                          ('XR', 'X-ray'), ('H2_form', 'H2 formation')]:
            rate = rates[key]
            plt.loglog(n, rate, label=label)
        
        plt.loglog(n, rates['total'], 'k-', lw=2, label='Total')
        plt.xlabel('n [cm$^{-3}$]')
        plt.ylabel('Heating rate [erg s$^{-1}$]')
        plt.title(f'Heating rates (N_H = {N_H:.0e} cm$^{{-2}}$)')
        plt.legend()
        plt.xlim(1e-2, 1e6)
        plt.ylim(1e-30, 1e-22)
        plt.grid(True, alpha=0.3)
        plt.savefig(f'/Users/guo/Downloads/sphcode/scripts/ki2000/physics/heating_test_N{int(np.log10(N_H))}.png', dpi=150)
        plt.close()
    
    print("Heating module test complete - plots saved.")
