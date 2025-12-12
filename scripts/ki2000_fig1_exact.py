#!/usr/bin/env python3
"""
Generate exact reproduction of K&I 2000 Figure 1.

This script plots the extracted data from K&I 2000 PostScript figures
to create a faithful reproduction of their Figure 1.

Figure 1 contains 4 panels:
  (a) Temperature T(n)
  (b) Chemical fractions x_e, x_H2, x_CO vs n
  (c) Heating/cooling rates (computed from first principles)
  (d) Thermal timescales

Each panel shows two cases:
  - Solid lines: N_H = 10^19 cm^-2
  - Dashed lines: N_H = 10^20 cm^-2

References:
  K&I: Koyama & Inutsuka 2000, ApJ 532, 980
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d


# =============================================================================
# Physical Constants
# =============================================================================
k_B = 1.380649e-16  # Boltzmann constant [erg/K]


# =============================================================================
# K&I 2000 Heating and Cooling Functions (Appendix A.1)
# =============================================================================

def heating_photoelectric(n: float, T: float, G_0: float = 1.7) -> float:
    """Photoelectric heating from dust grains (K&I eq. 1).

    Gamma_PE = 1.0e-24 * epsilon * G_0 [erg s^-1 H^-1]
    epsilon = 4.9e-2 / (1 + (psi/1925)^0.73) + 3.7e-2 * (T/1e4)^0.7 / (1 + psi/5000)
    psi = G_0 * sqrt(T) / n_e
    """
    # Need x_e to compute psi - use approximate value
    x_e = 0.01 if n > 10 else 0.1
    n_e = x_e * n

    if n_e < 1e-10:
        n_e = 1e-10

    psi = G_0 * np.sqrt(T) / n_e
    epsilon = 4.9e-2 / (1 + (psi / 1925)**0.73) + 3.7e-2 * (T / 1e4)**0.7 / (1 + psi / 5000)

    return 1.0e-24 * epsilon * G_0


def heating_xray(N_H: float = 1e19) -> float:
    """X-ray heating (K&I eq. 3).

    Gamma_XR = 8.9e-26 * exp(-N_19 / 60) [erg s^-1 H^-1]
    """
    N_19 = N_H / 1e19
    return 8.9e-26 * np.exp(-N_19 / 60)


def heating_cosmic_ray(x_e: float, zeta_CR: float = 1.8e-17) -> float:
    """Cosmic ray heating (K&I eq. 5).

    Gamma_CR = zeta_CR * (5.5e-12 + 8.2e-14/x_e) [erg s^-1 H^-1]

    Note: This includes both primary and secondary heating.
    """
    x_e = max(x_e, 1e-6)
    return zeta_CR * (5.5e-12 + 8.2e-14 / x_e)


def heating_h2_formation(n: float, T: float, x_H2: float) -> float:
    """H2 formation heating on grains (K&I eq. 6).

    Gamma_H2 = 2.4e-12 * R * n_HI [erg s^-1 H^-1]
    where R is the H2 formation rate coefficient.
    """
    # H2 formation rate (K&I eq. 16-17)
    T_gr = 8.0
    S = 1.0 / (1.0 + 0.04 * np.sqrt(T + T_gr) + 2e-3 * T + 8e-6 * T**2)
    R = 6e-17 * np.sqrt(T / 300.0) * S

    x_HI = max(1.0 - 2 * x_H2, 0.01)

    return 2.4e-12 * R * n * x_HI


def cooling_lyman_alpha(n: float, T: float, x_e: float) -> float:
    """Lyman-alpha cooling (K&I eq. 7).

    Lambda_Lya = 7.3e-19 * x_e * exp(-118400/T) [erg s^-1 H^-1]
    """
    if T < 3000:
        return 0.0
    return 7.3e-19 * x_e * np.exp(-118400 / T)


def cooling_OI(n: float, T: float, x_e: float) -> float:
    """OI 63 micron fine structure cooling (K&I eq. 8).

    Lambda_OI = 3.8e-33 * sqrt(T) * (1 - 2.4e-4 * T) * exp(-228/T) * n [erg s^-1 H^-2]
    """
    if T < 50:
        return 0.0
    return 3.8e-33 * np.sqrt(T) * (1 - 2.4e-4 * T) * np.exp(-228 / T) * n


def cooling_CII(n: float, T: float, x_e: float, N_H: float = 1e19) -> float:
    """CII 158 micron fine structure cooling (K&I eq. 9-10).

    Lambda_CII = 6.2e4 * A_C * exp(-92/T) *
                 [2.3e-8 * (psi/T)^0.5 + 1e-12] * exp(-tau_d) [erg s^-1 H^-1]
    """
    A_C = 3.0e-4  # Carbon abundance
    tau_d = N_H / 1.87e21

    # Need psi for rate
    x_e = max(x_e, 1e-6)
    n_e = x_e * n
    G_0 = 1.7
    psi = G_0 * np.sqrt(T) / max(n_e, 1e-10)

    rate = 2.3e-8 * (psi / T)**0.5 + 1e-12

    return 6.2e4 * A_C * np.exp(-92 / T) * rate * np.exp(-tau_d)


def cooling_CO(n: float, T: float, x_CO: float, x_H2: float) -> float:
    """CO rotational cooling (K&I eq. 11-13).

    Lambda_CO is complex, involving multiple terms.
    Simplified version for equilibrium.
    """
    if x_CO < 1e-10 or T < 5:
        return 0.0

    # Simplified CO cooling (order of magnitude)
    # Full expression in K&I eq. 11-13
    n_CO = x_CO * n
    n_H2 = x_H2 * n

    # CO-H2 collisional rate
    L_0 = 1e-24 * (T / 10)**2  # Approximate

    return L_0 * n_CO * n_H2 / n


def cooling_gas_grain(n: float, T: float) -> float:
    """Gas-grain cooling (K&I eq. 14).

    Lambda_gr = 3.8e-33 * sqrt(T) * (T - T_gr) * n [erg s^-1 H^-2]
    """
    T_gr = 8.0
    if T <= T_gr:
        return 0.0
    return 3.8e-33 * np.sqrt(T) * (T - T_gr) * n


def compute_heating_cooling(n: float, T: float, x_e: float, x_H2: float,
                            x_CO: float, N_H: float = 1e19, G_0: float = 1.7):
    """Compute total heating and cooling rates.

    Returns:
        heating: Total heating rate [erg s^-1 H^-1] (= n*Gamma)
        cooling: Total cooling rate [erg s^-1 H^-1] (= n^2*Lambda for collisional)
        components: Dict with individual contributions
    """
    # Heating rates [erg s^-1 H^-1]
    Gamma_PE = heating_photoelectric(n, T, G_0)
    Gamma_XR = heating_xray(N_H)
    Gamma_CR = heating_cosmic_ray(x_e)
    Gamma_H2 = heating_h2_formation(n, T, x_H2)

    # Cooling rates
    Lambda_Lya = cooling_lyman_alpha(n, T, x_e)
    Lambda_OI = cooling_OI(n, T, x_e)
    Lambda_CII = cooling_CII(n, T, x_e, N_H)
    Lambda_CO = cooling_CO(n, T, x_CO, x_H2)
    Lambda_gr = cooling_gas_grain(n, T)

    # Total rates
    total_heating = Gamma_PE + Gamma_XR + Gamma_CR + Gamma_H2
    total_cooling = Lambda_Lya + Lambda_OI + Lambda_CII + Lambda_CO + Lambda_gr

    components = {
        'PE': Gamma_PE,
        'XR': Gamma_XR,
        'CR': Gamma_CR,
        'H2_form': Gamma_H2,
        'Lya': Lambda_Lya,
        'OI': Lambda_OI,
        'CII': Lambda_CII,
        'CO': Lambda_CO,
        'grain': Lambda_gr,
    }

    return total_heating, total_cooling, components


def load_ki_data(data_dir: Path) -> dict:
    """Load all extracted K&I data."""
    data = {}

    def load_file(filename):
        filepath = data_dir / filename
        if filepath.exists():
            raw = np.loadtxt(filepath, comments='#')
            idx = np.argsort(raw[:, 0])
            return raw[idx, 0], raw[idx, 1]
        return None, None

    # Temperature
    data['T_N19'] = load_file('f1a_temperature_N19.txt')
    data['T_N20'] = load_file('f1a_temperature_N20.txt')

    # Pressure
    data['P_N19'] = load_file('f1a_pressure_N19.txt')
    data['P_N20'] = load_file('f1a_pressure_N20.txt')

    # Chemical fractions
    data['xe_N19'] = load_file('f1b_x_e_N19.txt')
    data['xe_N20'] = load_file('f1b_x_e_N20.txt')
    data['xH2_N19'] = load_file('f1b_x_H2_N19.txt')
    data['xH2_N20'] = load_file('f1b_x_H2_N20.txt')
    data['xCO_N19'] = load_file('f1b_x_CO_N19.txt')
    data['xCO_N20'] = load_file('f1b_x_CO_N20.txt')

    # Timescale curves
    for i in range(4):
        lt_map = {0: 1, 1: 2, 2: 4, 3: 6}
        data[f'ts_curve{i}'] = load_file(f'f1d_curve{i}_LT{lt_map[i]}.txt')

    # Panel (c) - Heating/Cooling curves (extracted from f1c.ps)
    # Load all f1c curves for both N19 and N20 cases
    data['f1c_N19'] = []
    data['f1c_N20'] = []
    for i in range(10):  # Up to 10 curves per case
        n19_curve = load_file(f'f1c_N19_curve{i}.txt')
        n20_curve = load_file(f'f1c_N20_curve{i}.txt')
        if n19_curve[0] is not None:
            data['f1c_N19'].append(n19_curve)
        if n20_curve[0] is not None:
            data['f1c_N20'].append(n20_curve)

    return data


def load_extracted_panel_c(data_dir: Path) -> dict:
    """Load the extracted heating/cooling curves from f1c.ps.
    
    Returns dict with curves organized by process and column density case.
    Based on curve analysis and label positions in f1c.ps.
    """
    def load_file(filename):
        filepath = data_dir / filename
        if filepath.exists():
            raw = np.loadtxt(filepath, comments='#')
            idx = np.argsort(raw[:, 0])
            return raw[idx, 0], raw[idx, 1]
        return None, None

    result = {
        'N19': {'curves': []},
        'N20': {'curves': []},
    }
    
    # Load N19 curves (LT1 = solid lines)
    for i in range(10):
        curve = load_file(f'f1c_N19_curve{i}.txt')
        if curve[0] is not None:
            result['N19']['curves'].append(curve)
    
    # Load N20 curves (LT2 = dashed lines)
    for i in range(10):
        curve = load_file(f'f1c_N20_curve{i}.txt')
        if curve[0] is not None:
            result['N20']['curves'].append(curve)
    
    return result


def compute_panel_c_data(data: dict) -> dict:
    """Compute heating/cooling rates from first principles using panel (a) and (b) data.

    Uses the extracted T(n) and x_i(n) data to compute heating and cooling rates
    at thermal equilibrium for each density.

    Returns dict with heating/cooling curves for N19 and N20 cases.
    """
    result = {}

    for case in ['N19', 'N20']:
        N_H = 1e19 if case == 'N19' else 1e20

        # Get data
        n_T, T = data[f'T_{case}']
        n_xe, xe = data[f'xe_{case}']
        n_h2, xh2 = data[f'xH2_{case}']
        n_co, xco = data[f'xCO_{case}']

        # Create interpolators
        T_interp = interp1d(np.log10(n_T), np.log10(T),
                            kind='linear', fill_value='extrapolate')
        xe_interp = interp1d(np.log10(n_xe), np.log10(xe),
                             kind='linear', fill_value='extrapolate')
        xh2_interp = interp1d(np.log10(n_h2), np.log10(xh2),
                              kind='linear', fill_value='extrapolate')
        xco_interp = interp1d(np.log10(n_co), np.log10(xco),
                              kind='linear', fill_value='extrapolate')

        # Compute on grid
        n_grid = np.logspace(-2, 6, 200)
        heating = np.zeros_like(n_grid)
        cooling = np.zeros_like(n_grid)
        components = {k: np.zeros_like(n_grid) for k in
                      ['PE', 'XR', 'CR', 'H2_form', 'Lya', 'OI', 'CII', 'CO', 'grain']}

        for i, n in enumerate(n_grid):
            log_n = np.log10(n)

            # Get interpolated values
            T_val = 10**T_interp(log_n)
            xe_val = 10**xe_interp(log_n)
            xh2_val = 10**xh2_interp(log_n) if log_n >= np.log10(n_h2.min()) else 1e-10
            xco_val = 10**xco_interp(log_n) if log_n >= np.log10(n_co.min()) else 1e-15

            # Clip to physical range
            T_val = np.clip(T_val, 3, 1e5)
            xe_val = np.clip(xe_val, 1e-6, 1)
            xh2_val = np.clip(xh2_val, 0, 0.5)
            xco_val = np.clip(xco_val, 0, 3e-4)

            h, c, comp = compute_heating_cooling(n, T_val, xe_val, xh2_val, xco_val, N_H)
            heating[i] = h
            cooling[i] = c
            for k in components:
                components[k][i] = comp[k]

        result[case] = {
            'n': n_grid,
            'heating': heating,
            'cooling': cooling,
            'components': components,
        }

    return result


def plot_ki_fig1(data: dict, output_file: Path, hc_data: dict = None):
    """Create K&I 2000 Figure 1 reproduction."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Colors and styles
    style_N19 = {'color': 'black', 'linestyle': '-', 'linewidth': 1.5}
    style_N20 = {'color': 'black', 'linestyle': '--', 'linewidth': 1.5}

    # ==========================================================================
    # Panel (a) - Temperature
    # ==========================================================================
    ax = axes[0, 0]

    n_T19, T_19 = data['T_N19']
    n_T20, T_20 = data['T_N20']
    ax.loglog(n_T19, T_19, **style_N19, label=r'$N_H = 10^{19}$ cm$^{-2}$')
    ax.loglog(n_T20, T_20, **style_N20, label=r'$N_H = 10^{20}$ cm$^{-2}$')

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$T$ [K]')
    ax.set_xlim(1e-2, 1e4)
    ax.set_ylim(1, 2e4)
    ax.legend(loc='upper right', fontsize=9)
    ax.text(0.05, 0.95, '(a)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    # ==========================================================================
    # Panel (b) - Chemical fractions
    # ==========================================================================
    ax = axes[0, 1]

    # Electron fraction
    n_xe19, xe_19 = data['xe_N19']
    n_xe20, xe_20 = data['xe_N20']
    ax.loglog(n_xe19, xe_19, **style_N19)
    ax.loglog(n_xe20, xe_20, **style_N20)

    # H2 fraction
    n_h2_19, xh2_19 = data['xH2_N19']
    n_h2_20, xh2_20 = data['xH2_N20']
    ax.loglog(n_h2_19, xh2_19, **style_N19)
    ax.loglog(n_h2_20, xh2_20, **style_N20)

    # CO fraction
    n_co_19, xco_19 = data['xCO_N19']
    n_co_20, xco_20 = data['xCO_N20']
    ax.loglog(n_co_19, xco_19, **style_N19)
    ax.loglog(n_co_20, xco_20, **style_N20)

    # Labels (positioned like K&I)
    ax.text(0.3, 6e-2, 'electron', fontsize=10, style='italic')
    ax.text(3e5, 0.5, r'H$_2$', fontsize=10, style='italic')
    ax.text(3e5, 5e-5, 'CO', fontsize=10, style='italic')

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$x_i$')
    ax.set_xlim(1e-2, 1e6)
    ax.set_ylim(1e-8, 1)
    ax.text(0.05, 0.95, '(b)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    # ==========================================================================
    # Panel (c) - Heating/Cooling rates (extracted from K&I 2000 f1c.ps)
    # ==========================================================================
    ax = axes[1, 0]

    # Use extracted data from f1c.ps instead of computed values
    f1c_curves = data.get('f1c_N19', [])
    f1c_curves_N20 = data.get('f1c_N20', [])
    
    if f1c_curves or f1c_curves_N20:
        # Plot N19 curves (solid lines)
        for i, (n, rate) in enumerate(f1c_curves):
            ax.loglog(n, rate, 'k-', lw=1.5, alpha=0.8)
        
        # Plot N20 curves (dashed lines)
        for i, (n, rate) in enumerate(f1c_curves_N20):
            ax.loglog(n, rate, 'k--', lw=1.5, alpha=0.8)
        
        # Add labels matching K&I figure positions
        # (positions derived from f1c.ps label coordinates)
        ax.text(0.3, 2e-26, 'PE', fontsize=9, style='italic')
        ax.text(1e5, 1e-27, 'XR', fontsize=9, style='italic')
        ax.text(1.5e4, 3e-28, 'CR', fontsize=9, style='italic')
        ax.text(1.5e4, 3e-26, r'H$_2$', fontsize=9, style='italic')
        ax.text(0.02, 8e-27, r'Ly$\alpha$', fontsize=9, style='italic')
        ax.text(10, 4e-27, 'OI', fontsize=9, style='italic')
        ax.text(10, 3e-26, 'CII', fontsize=9, style='italic')
        ax.text(2e4, 3e-27, 'GR', fontsize=9, style='italic')
        ax.text(6e5, 1.5e-25, 'CO', fontsize=9, style='italic')
    
    elif hc_data is not None:
        # Fall back to computed data if extracted data not available
        for case, style in [('N19', '-'), ('N20', '--')]:
            hc = hc_data[case]
            n = hc['n']
            ax.loglog(n, hc['heating'], f'k{style}', lw=2.0, alpha=0.8)
            ax.loglog(n, hc['cooling'], f'k{style}', lw=2.0, alpha=0.8)

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$\Gamma$, $\Lambda$ [erg s$^{-1}$ H$^{-1}$]')
    ax.set_xlim(1e-2, 1e6)
    ax.set_ylim(1e-28, 1e-23)
    ax.text(0.05, 0.95, '(c)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    # ==========================================================================
    # Panel (d) - Timescales
    # ==========================================================================
    ax = axes[1, 1]

    # Timescale curve labels (from K&I 2000)
    ts_labels = [r'$t_{\rm cool}$', r'$t_{\rm rec}$', r'$t_{\rm ff}$', r'$t_{\rm H_2}$']
    ts_styles = ['-', '--', '-.', ':']

    for i in range(4):
        ts = data.get(f'ts_curve{i}')
        if ts[0] is not None:
            n, t = ts
            ax.loglog(n, t, f'k{ts_styles[i]}', lw=1.5, label=ts_labels[i])

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$t$ [yr]')
    ax.set_xlim(1e-2, 1e4)
    ax.set_ylim(1e2, 1e10)
    ax.legend(loc='upper right', fontsize=9)
    ax.text(0.05, 0.95, '(d)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    plt.suptitle('K&I 2000 Figure 1', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.savefig(output_file.with_suffix('.pdf'), bbox_inches='tight')
    print(f'Saved: {output_file}')


def plot_ki_fig1b_detailed(data: dict, output_file: Path):
    """Create detailed Figure 1b showing chemical fractions."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Colors matching K&I style
    colors = {'N19': 'blue', 'N20': 'red'}
    styles = {'N19': '-', 'N20': '--'}

    # Electron fraction
    n_xe19, xe_19 = data['xe_N19']
    n_xe20, xe_20 = data['xe_N20']
    ax.loglog(n_xe19, xe_19, color=colors['N19'], linestyle=styles['N19'],
              lw=2, label=r'$x_e$ ($N_H=10^{19}$)')
    ax.loglog(n_xe20, xe_20, color=colors['N20'], linestyle=styles['N20'],
              lw=2, label=r'$x_e$ ($N_H=10^{20}$)')

    # H2 fraction
    n_h2_19, xh2_19 = data['xH2_N19']
    n_h2_20, xh2_20 = data['xH2_N20']
    ax.loglog(n_h2_19, xh2_19, color=colors['N19'], linestyle=styles['N19'], lw=2)
    ax.loglog(n_h2_20, xh2_20, color=colors['N20'], linestyle=styles['N20'], lw=2)

    # CO fraction
    n_co_19, xco_19 = data['xCO_N19']
    n_co_20, xco_20 = data['xCO_N20']
    ax.loglog(n_co_19, xco_19, color=colors['N19'], linestyle=styles['N19'], lw=2)
    ax.loglog(n_co_20, xco_20, color=colors['N20'], linestyle=styles['N20'], lw=2)

    # Labels
    ax.text(5, 2e-2, 'electron', fontsize=12, fontweight='bold')
    ax.text(2e5, 0.4, r'H$_2$', fontsize=12, fontweight='bold')
    ax.text(2e5, 8e-5, 'CO', fontsize=12, fontweight='bold')

    ax.set_xlabel(r'$\log\, n$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel(r'$\log\, x_i$', fontsize=12)
    ax.set_xlim(1e-2, 1e6)
    ax.set_ylim(1e-8, 1)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('K&I 2000 Figure 1b: Chemical Fractions', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved: {output_file}')


def plot_ki_fig1_with_model(data: dict, output_file: Path, hc_data: dict = None):
    """Create K&I 2000 Figure 1 with model comparison."""
    # Try to import the chemistry model
    try:
        from ki2000_chemistry_ode import ChemParams, solve_chemistry_steady_state
        have_model = True
    except ImportError:
        have_model = False
        print("Warning: Could not import chemistry model")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Colors and styles
    style_N19 = {'color': 'black', 'linestyle': '-', 'linewidth': 1.5}
    style_N20 = {'color': 'black', 'linestyle': '--', 'linewidth': 1.5}
    style_model = {'color': 'red', 'linestyle': '-', 'linewidth': 1.2, 'alpha': 0.8}

    # ==========================================================================
    # Panel (a) - Temperature
    # ==========================================================================
    ax = axes[0, 0]

    n_T19, T_19 = data['T_N19']
    n_T20, T_20 = data['T_N20']
    ax.loglog(n_T19, T_19, **style_N19, label=r'$N_H = 10^{19}$ cm$^{-2}$')
    ax.loglog(n_T20, T_20, **style_N20, label=r'$N_H = 10^{20}$ cm$^{-2}$')

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$T$ [K]')
    ax.set_xlim(1e-2, 1e4)
    ax.set_ylim(1, 2e4)
    ax.legend(loc='upper right', fontsize=9)
    ax.text(0.05, 0.95, '(a)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    # ==========================================================================
    # Panel (b) - Chemical fractions with model comparison
    # ==========================================================================
    ax = axes[0, 1]

    # K&I data
    n_xe19, xe_19 = data['xe_N19']
    n_xe20, xe_20 = data['xe_N20']
    n_h2_19, xh2_19 = data['xH2_N19']
    n_h2_20, xh2_20 = data['xH2_N20']
    n_co_19, xco_19 = data['xCO_N19']
    n_co_20, xco_20 = data['xCO_N20']

    # Plot K&I data
    ax.loglog(n_xe19, xe_19, **style_N19, label='K&I 2000')
    ax.loglog(n_xe20, xe_20, **style_N20)
    ax.loglog(n_h2_19, xh2_19, **style_N19)
    ax.loglog(n_h2_20, xh2_20, **style_N20)
    ax.loglog(n_co_19, xco_19, **style_N19)
    ax.loglog(n_co_20, xco_20, **style_N20)

    # Compute model if available
    if have_model:
        # Create temperature interpolator from K&I data
        T_interp_N19 = interp1d(np.log10(n_T19), np.log10(T_19),
                                 kind='linear', fill_value='extrapolate')

        # Compute chemistry at K&I equilibrium temperatures
        n_model = np.logspace(-2, 6, 100)
        params_N19 = ChemParams(N_H=1e19)

        x_e_model = np.zeros_like(n_model)
        x_H2_model = np.zeros_like(n_model)
        x_CO_model = np.zeros_like(n_model)

        for i, n in enumerate(n_model):
            # Get temperature from K&I thermal equilibrium
            if np.log10(n) <= np.log10(n_T19.max()):
                T = 10**T_interp_N19(np.log10(n))
            else:
                T = 8.0  # CNM temperature for high density

            result = solve_chemistry_steady_state(n, T, params_N19)
            x_e_model[i] = result['x_e']
            x_H2_model[i] = result['x_H2']
            x_CO_model[i] = result['x_CO']

        # Plot model
        ax.loglog(n_model, x_e_model, **style_model, label='Model')
        ax.loglog(n_model, 2 * x_H2_model + 1e-10, **style_model)  # x_2 = 2*x_H2
        ax.loglog(n_model, x_CO_model + 1e-12, **style_model)

    # Labels
    ax.text(0.3, 6e-2, 'electron', fontsize=10, style='italic')
    ax.text(3e5, 0.5, r'H$_2$', fontsize=10, style='italic')
    ax.text(3e5, 5e-5, 'CO', fontsize=10, style='italic')

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$x_i$')
    ax.set_xlim(1e-2, 1e6)
    ax.set_ylim(1e-8, 1)
    ax.legend(loc='lower left', fontsize=9)
    ax.text(0.05, 0.95, '(b)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    # ==========================================================================
    # Panel (c) - Heating/Cooling rates (computed from first principles)
    # ==========================================================================
    ax = axes[1, 0]

    if hc_data is not None:
        # Plot total heating and cooling for both cases
        for case, style in [('N19', '-'), ('N20', '--')]:
            hc = hc_data[case]
            n = hc['n']

            # Total heating (thick line)
            ax.loglog(n, hc['heating'], f'k{style}', lw=2.0, alpha=0.8)

            # Total cooling (thick line)
            ax.loglog(n, hc['cooling'], f'k{style}', lw=2.0, alpha=0.8)

            # Individual components (thinner lines)
            comp = hc['components']

            # Heating components
            for name, color in [('PE', 'orange'), ('XR', 'purple'), ('CR', 'green')]:
                vals = comp[name]
                mask = vals > 1e-30
                if mask.any():
                    ax.loglog(n[mask], vals[mask], color=color, linestyle=style,
                              lw=1.0, alpha=0.6)

            # Cooling components
            for name, color in [('Lya', 'blue'), ('CII', 'red'), ('OI', 'cyan'),
                                ('CO', 'brown'), ('grain', 'gray')]:
                vals = comp[name]
                mask = vals > 1e-30
                if mask.any():
                    ax.loglog(n[mask], vals[mask], color=color, linestyle=style,
                              lw=1.0, alpha=0.6)

        # Add labels
        ax.text(2e-2, 1.5e-25, 'PE', fontsize=9, style='italic', color='orange')
        ax.text(1e4, 1e-26, 'XR', fontsize=9, style='italic', color='purple')
        ax.text(1e3, 3e-26, 'CR', fontsize=9, style='italic', color='green')
        ax.text(3e-1, 3e-26, r'Ly$\alpha$', fontsize=9, style='italic', color='blue')
        ax.text(1e2, 2e-24, 'CII', fontsize=9, style='italic', color='red')
        ax.text(1e5, 1e-24, 'CO', fontsize=9, style='italic', color='brown')

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$\Gamma$, $\Lambda$ [erg s$^{-1}$ H$^{-1}$]')
    ax.set_xlim(1e-2, 1e6)
    ax.set_ylim(1e-28, 1e-23)
    ax.text(0.05, 0.95, '(c)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    # ==========================================================================
    # Panel (d) - Timescales
    # ==========================================================================
    ax = axes[1, 1]

    ts_labels = [r'$t_{\rm cool}$', r'$t_{\rm rec}$', r'$t_{\rm ff}$', r'$t_{\rm H_2}$']
    ts_styles = ['-', '--', '-.', ':']

    for i in range(4):
        ts = data.get(f'ts_curve{i}')
        if ts[0] is not None:
            n, t = ts
            ax.loglog(n, t, f'k{ts_styles[i]}', lw=1.5, label=ts_labels[i])

    ax.set_xlabel(r'$n$ [cm$^{-3}$]')
    ax.set_ylabel(r'$t$ [yr]')
    ax.set_xlim(1e-2, 1e4)
    ax.set_ylim(1e2, 1e10)
    ax.legend(loc='upper right', fontsize=9)
    ax.text(0.05, 0.95, '(d)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.grid(True, alpha=0.3, which='both')

    plt.suptitle('K&I 2000 Figure 1 with First-Principles Model Comparison',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved: {output_file}')


def print_table(data: dict):
    """Print table of values at key densities."""
    print("\n" + "=" * 80)
    print("K&I 2000 Chemical Equilibrium Values (N_H = 10^19 cm^-2)")
    print("=" * 80)
    print(f"{'n [cm^-3]':>12} {'T [K]':>10} {'x_e':>12} {'x_H2':>12} {'x_CO':>12}")
    print("-" * 80)

    n_T, T = data['T_N19']
    n_xe, xe = data['xe_N19']
    n_h2, xh2 = data['xH2_N19']
    n_co, xco = data['xCO_N19']

    target_densities = [0.01, 0.1, 1, 10, 100, 1000, 1e4]

    for target_n in target_densities:
        # Find closest values
        idx_T = np.argmin(np.abs(n_T - target_n))
        idx_xe = np.argmin(np.abs(n_xe - target_n))
        idx_h2 = np.argmin(np.abs(n_h2 - target_n)) if n_h2 is not None else -1
        idx_co = np.argmin(np.abs(n_co - target_n)) if n_co is not None else -1

        T_val = T[idx_T] if abs(np.log10(n_T[idx_T]/target_n)) < 0.5 else np.nan
        xe_val = xe[idx_xe] if abs(np.log10(n_xe[idx_xe]/target_n)) < 0.5 else np.nan
        xh2_val = xh2[idx_h2] if idx_h2 >= 0 and abs(np.log10(n_h2[idx_h2]/target_n)) < 0.5 else np.nan
        xco_val = xco[idx_co] if idx_co >= 0 and abs(np.log10(n_co[idx_co]/target_n)) < 0.5 else np.nan

        print(f"{target_n:>12.2e} {T_val:>10.1f} {xe_val:>12.4e} {xh2_val:>12.4e} {xco_val:>12.4e}")


def main():
    """Main function."""
    import sys
    data_dir = Path(__file__).parent.parent / 'data' / 'ki2000_extracted'

    # Load data
    print("Loading extracted K&I data...")
    data = load_ki_data(data_dir)

    # Compute panel (c) from first principles
    print("Computing heating/cooling rates from first principles...")
    hc_data = compute_panel_c_data(data)

    # Print table
    print_table(data)

    # Create plots
    print("\nGenerating figures...")
    plot_ki_fig1(data, data_dir / 'ki2000_fig1_exact.png', hc_data)
    plot_ki_fig1b_detailed(data, data_dir / 'ki2000_fig1b_detailed.png')

    # Create model comparison
    print("Generating model comparison (this may take a moment)...")
    plot_ki_fig1_with_model(data, data_dir / 'ki2000_fig1_with_model.png', hc_data)

    if '--show' in sys.argv:
        plt.show()


if __name__ == '__main__':
    main()
