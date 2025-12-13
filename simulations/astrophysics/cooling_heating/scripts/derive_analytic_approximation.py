#!/usr/bin/env python3
"""
Derive physically-motivated analytic approximations for T(n) and P(n).
Goal: Simple formulas capturing thermal bistability without solving chemistry.

Physical basis:
- WNM branch: Nearly isothermal at T~8000K (photoelectric heating balances cooling)
- Thermal transition: Rapid drop over narrow density range (thermal instability)
- CNM branch: T decreases slowly with density (dominated by grain cooling)
- Molecular regime: Further cooling as H2 forms

Strategy: Piecewise smooth transitions using hyperbolic tangent functions
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def load_exact_data():
    """Load the exact digitized data."""
    T_1e19 = np.loadtxt('../results/f1a_curve_0.txt')
    T_1e20 = np.loadtxt('../results/f1a_curve_1.txt')
    P_1e19 = np.loadtxt('../results/f1a_curve_2.txt')
    P_1e20 = np.loadtxt('../results/f1a_curve_3.txt')
    
    return T_1e19, T_1e20, P_1e19, P_1e20

def temperature_approximation(n, T_WNM, T_CNM, T_mol, n_trans, delta_n, n_mol, delta_mol):
    """
    Three-phase temperature model with smooth transitions.
    
    Parameters:
    -----------
    n : array
        Density [cm^-3]
    T_WNM : float
        WNM temperature [K] (warm branch)
    T_CNM : float
        CNM temperature [K] (cold branch)
    T_mol : float
        Molecular cloud temperature [K]
    n_trans : float
        Transition density from WNM to CNM [cm^-3]
    delta_n : float
        Transition width [cm^-3]
    n_mol : float
        Molecular transition density [cm^-3]
    delta_mol : float
        Molecular transition width [cm^-3]
    
    Returns:
    --------
    T : array
        Temperature [K]
    """
    
    log_n = np.log10(n)
    log_n_trans = np.log10(n_trans)
    log_delta_n = np.log10(delta_n)
    log_n_mol = np.log10(n_mol)
    log_delta_mol = np.log10(delta_mol)
    
    # WNM to CNM transition (S-curve)
    # Use tanh for smooth transition
    transition_1 = 0.5 * (1 + np.tanh((log_n - log_n_trans) / log_delta_n))
    
    # CNM to molecular transition
    transition_2 = 0.5 * (1 + np.tanh((log_n - log_n_mol) / log_delta_mol))
    
    # Three-phase temperature
    T = T_WNM * (1 - transition_1) + \
        T_CNM * transition_1 * (1 - transition_2) + \
        T_mol * transition_2
    
    return T

def pressure_approximation(n, P_min, n_min, alpha_low, alpha_high):
    """
    Pressure model with minimum at transition.
    
    Physical basis:
    - At low n: P ~ n*T_WNM (isothermal WNM)
    - Minimum at thermal transition
    - At high n: P increases due to compression and self-gravity
    
    Parameters:
    -----------
    n : array
        Density [cm^-3]
    P_min : float
        Minimum pressure at transition [K/cm^3]
    n_min : float
        Density at pressure minimum [cm^-3]
    alpha_low : float
        Power-law index for n < n_min
    alpha_high : float
        Power-law index for n > n_min
    """
    
    P = np.zeros_like(n)
    
    # Low density: P ~ n^alpha_low
    mask_low = n < n_min
    P[mask_low] = P_min * (n[mask_low] / n_min)**alpha_low
    
    # High density: P ~ n^alpha_high
    mask_high = n >= n_min
    P[mask_high] = P_min * (n[mask_high] / n_min)**alpha_high
    
    return P

def fit_temperature_model(n_data, T_data, initial_guess):
    """Fit the temperature approximation to exact data."""
    
    # Use log-space fitting for better behavior
    def fit_func(n, T_WNM, T_CNM, T_mol, n_trans, delta_n, n_mol, delta_mol):
        return np.log10(temperature_approximation(n, T_WNM, T_CNM, T_mol, 
                                                   n_trans, delta_n, n_mol, delta_mol))
    
    popt, pcov = curve_fit(fit_func, n_data, np.log10(T_data), 
                           p0=initial_guess, maxfev=10000)
    
    return popt

def simplified_temperature(n, N_H_col=1e19):
    """
    Simple analytic formula for temperature based on physical phases.
    
    This is the formula that can be used in SPH codes!
    
    Model: Smooth power-law cooling T ~ n^alpha with density-dependent slope
    """
    
    n_arr = np.atleast_1d(n)
    log_n = np.log10(n_arr)
    
    if N_H_col >= 5e19:  # Higher column density (more shielding)
        # Reference point
        T_ref = 7600.0     # Temperature at n=1
        n_ref = 1.0
        
        # Cooling parameters - transitions from slow to fast cooling
        alpha_low = -0.05   # Nearly isothermal at low n
        alpha_high = -0.85  # Steep cooling in transition
        alpha_cnm = -0.15   # Slow cooling in CNM
        
        # Transition points
        n_start = 1.5       # Start significant cooling
        n_steep = 15.0      # Steepest cooling
        n_cnm = 80.0        # CNM regime
        
        width_start = 0.4
        width_steep = 0.35
        
    else:  # Lower column density (less shielding)
        # Reference point
        T_ref = 8000.0
        n_ref = 1.0
        
        # Cooling parameters
        alpha_low = -0.03
        alpha_high = -0.80
        alpha_cnm = -0.12
        
        # Transition points
        n_start = 2.0
        n_steep = 20.0
        n_cnm = 120.0
        
        width_start = 0.45
        width_steep = 0.4
    
    # Smooth transitions for cooling index
    f_start = 0.5 * (1.0 + np.tanh((log_n - np.log10(n_start)) / width_start))
    f_steep = 0.5 * (1.0 + np.tanh((log_n - np.log10(n_steep)) / width_steep))
    
    # Blended cooling index
    alpha = alpha_low * (1.0 - f_start) + alpha_high * f_start * (1.0 - f_steep) + alpha_cnm * f_steep
    
    # Apply power law with variable index
    T = T_ref * (n_arr / n_ref) ** alpha
    
    return np.squeeze(T) if np.isscalar(n) else T

def simplified_pressure(n, N_H_col=1e19):
    """
    Simple analytic formula for pressure P/k_B [K cm^-3].
    
    This is the formula that can be used in SPH codes!
    
    Model: The thermal pressure follows the S-curve with a pressure minimum.
    At the thermal transition, pressure drops due to rapid cooling.
    """
    
    n_arr = np.atleast_1d(n)
    log_n = np.log10(n_arr)
    
    if N_H_col >= 5e19:  # Higher column density
        # Pressure parameters for high shielding
        P_WNM = 2500.0      # WNM pressure level [K cm^-3]
        P_min_factor = 0.25 # Pressure drops to 25% at transition
        P_CNM_slope = 0.7   # CNM pressure slope (P ~ n^0.7)
        n_trans = 2.5       # Transition density
        width = 0.25        # Transition width
    else:  # Lower column density
        # Pressure parameters for low shielding
        P_WNM = 3000.0      # Slightly higher WNM pressure
        P_min_factor = 0.20 # Deeper pressure minimum
        P_CNM_slope = 0.75  # Slightly steeper CNM slope
        n_trans = 4.0       # Later transition
        width = 0.3
    
    # WNM branch: nearly constant pressure
    P_warm = P_WNM * np.ones_like(n_arr)
    
    # Pressure minimum at transition (thermal instability)
    P_min = P_WNM * P_min_factor
    
    # CNM branch: pressure increases with density
    P_cold = P_min * (n_arr / n_trans) ** P_CNM_slope
    
    # Smooth transition using tanh
    f = 0.5 * (1.0 + np.tanh((log_n - np.log10(n_trans)) / width))
    
    # Blend WNM and CNM branches
    P = P_warm * (1.0 - f) + P_cold * f
    
    # Ensure monotonic increase at high densities
    P = np.maximum(P, P_min * (n_arr / n_trans) ** 0.5)
    
    return np.squeeze(P) if np.isscalar(n) else P

def main():
    print("="*70)
    print("DERIVING ANALYTIC APPROXIMATIONS")
    print("="*70)
    
    # Load exact data
    T_1e19, T_1e20, P_1e19, P_1e20 = load_exact_data()
    
    print("\n✓ Loaded exact digitized data")
    
    # Create test points
    n_test = np.logspace(-1, 6, 500)
    
    # Calculate approximations
    T_approx_1e19 = simplified_temperature(n_test, N_H_col=1e19)
    T_approx_1e20 = simplified_temperature(n_test, N_H_col=1e20)
    P_approx_1e19 = simplified_pressure(n_test, N_H_col=1e19)
    P_approx_1e20 = simplified_pressure(n_test, N_H_col=1e20)
    
    print("✓ Calculated analytic approximations")
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # Temperature comparison - N_H = 1e19
    ax1 = axes[0, 0]
    ax1.loglog(T_1e19[:,0], T_1e19[:,1], 'b-', linewidth=3, 
              label='Exact (digitized)', alpha=0.7)
    ax1.loglog(n_test, T_approx_1e19, 'r--', linewidth=2.5,
              label='Analytic approximation')
    ax1.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax1.set_ylabel(r'$T$ [K]', fontsize=12, fontweight='bold')
    ax1.set_xlim(0.1, 1e6)
    ax1.set_ylim(10, 1e5)
    ax1.legend(fontsize=11, loc='best')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_title(r'Temperature: $N_H=10^{19}$ cm$^{-2}$', 
                  fontsize=13, fontweight='bold')
    
    # Temperature comparison - N_H = 1e20
    ax2 = axes[0, 1]
    ax2.loglog(T_1e20[:,0], T_1e20[:,1], 'b-', linewidth=3,
              label='Exact (digitized)', alpha=0.7)
    ax2.loglog(n_test, T_approx_1e20, 'r--', linewidth=2.5,
              label='Analytic approximation')
    ax2.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax2.set_ylabel(r'$T$ [K]', fontsize=12, fontweight='bold')
    ax2.set_xlim(0.1, 1e6)
    ax2.set_ylim(10, 1e5)
    ax2.legend(fontsize=11, loc='best')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_title(r'Temperature: $N_H=10^{20}$ cm$^{-2}$',
                  fontsize=13, fontweight='bold')
    
    # Pressure comparison - N_H = 1e19
    ax3 = axes[1, 0]
    ax3.loglog(P_1e19[:,0], P_1e19[:,1], 'b-', linewidth=3,
              label='Exact (digitized)', alpha=0.7)
    ax3.loglog(n_test, P_approx_1e19, 'r--', linewidth=2.5,
              label='Analytic approximation')
    ax3.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax3.set_ylabel(r'$P/k_B$ [K cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax3.set_xlim(0.1, 1e6)
    ax3.set_ylim(100, 1e8)
    ax3.legend(fontsize=11, loc='best')
    ax3.grid(True, alpha=0.3, which='both')
    ax3.set_title(r'Pressure: $N_H=10^{19}$ cm$^{-2}$',
                  fontsize=13, fontweight='bold')
    
    # Pressure comparison - N_H = 1e20
    ax4 = axes[1, 1]
    ax4.loglog(P_1e20[:,0], P_1e20[:,1], 'b-', linewidth=3,
              label='Exact (digitized)', alpha=0.7)
    ax4.loglog(n_test, P_approx_1e20, 'r--', linewidth=2.5,
              label='Analytic approximation')
    ax4.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax4.set_ylabel(r'$P/k_B$ [K cm$^{-3}$]', fontsize=12, fontweight='bold')
    ax4.set_xlim(0.1, 1e6)
    ax4.set_ylim(100, 1e8)
    ax4.legend(fontsize=11, loc='best')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_title(r'Pressure: $N_H=10^{20}$ cm$^{-2}$',
                  fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('../results/analytic_approximation_comparison.png', 
                dpi=300, bbox_inches='tight')
    print("✓ Saved: ../results/analytic_approximation_comparison.png")
    
    # Calculate RMS errors
    print("\n" + "="*70)
    print("APPROXIMATION QUALITY")
    print("="*70)
    
    # Interpolate approximations to exact data points
    T_approx_at_exact_1e19 = simplified_temperature(T_1e19[:,0], N_H_col=1e19)
    T_approx_at_exact_1e20 = simplified_temperature(T_1e20[:,0], N_H_col=1e20)
    P_approx_at_exact_1e19 = simplified_pressure(P_1e19[:,0], N_H_col=1e19)
    P_approx_at_exact_1e20 = simplified_pressure(P_1e20[:,0], N_H_col=1e20)
    
    # RMS errors in log space
    T_rms_1e19 = np.sqrt(np.mean((np.log10(T_approx_at_exact_1e19) - 
                                  np.log10(T_1e19[:,1]))**2))
    T_rms_1e20 = np.sqrt(np.mean((np.log10(T_approx_at_exact_1e20) - 
                                  np.log10(T_1e20[:,1]))**2))
    P_rms_1e19 = np.sqrt(np.mean((np.log10(P_approx_at_exact_1e19) - 
                                  np.log10(P_1e19[:,1]))**2))
    P_rms_1e20 = np.sqrt(np.mean((np.log10(P_approx_at_exact_1e20) - 
                                  np.log10(P_1e20[:,1]))**2))
    
    print("\nRMS error in log10(T):")
    print(f"  N_H = 1e19: {T_rms_1e19:.3f} dex")
    print(f"  N_H = 1e20: {T_rms_1e20:.3f} dex")
    
    print("\nRMS error in log10(P):")
    print(f"  N_H = 1e19: {P_rms_1e19:.3f} dex")
    print(f"  N_H = 1e20: {P_rms_1e20:.3f} dex")
    
    # Check key physical values
    print("\n" + "="*70)
    print("KEY VALUES COMPARISON")
    print("="*70)
    
    for n_check in [0.1, 1.0, 10.0, 100.0, 1000.0]:
        T_exact = np.interp(n_check, T_1e19[::-1,0], T_1e19[::-1,1])
        T_approx = simplified_temperature(n_check, N_H_col=1e19)
        error = abs(T_exact - T_approx) / T_exact * 100
        print(f"n = {n_check:6.1f} cm^-3: T_exact = {T_exact:7.1f} K, "
              f"T_approx = {T_approx:7.1f} K, error = {error:5.1f}%")
    
    print("\n" + "="*70)
    print("✓ ANALYTIC APPROXIMATION DERIVED SUCCESSFULLY!")
    print("="*70)

if __name__ == '__main__':
    main()
