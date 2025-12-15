#!/usr/bin/env python3
"""
Plot Inoue & Inutsuka (2008) ISM Cooling Function

Visualizes the cooling function properties:
- Cooling/heating coefficients vs temperature
- Net cooling rate vs temperature for different densities
- Thermal equilibrium curve (temperature vs density)
- Thermal stability diagram
- Cooling timescale maps

Reference: Inoue, T. & Inutsuka, S. (2008), ApJ
"Two-Fluid MHD Simulations of Converging HI Flows in the Interstellar Medium"

Usage:
    python scripts/plot_inoue_inutsuka_cooling.py [--output-dir DIR] [--dpi DPI]
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse
from pathlib import Path

# Physical constants (CGS)
k_B = 1.380649e-16      # erg K^-1
m_proton = 1.6726219e-24  # g
m_n = 1.27 * m_proton    # Mean neutral particle mass (91% H + 9% He)
Gamma_0 = 2.0e-26        # erg s^-1 (heating rate constant)
gamma = 5.0 / 3.0        # Adiabatic index

# Temperature bounds
T_MIN = 10.0
T_MAX = 1.0e8


def cooling_coefficient_ratio(T):
    """
    Cooling coefficient ratio Λ/Γ (Eq. 9, corrected)
    
    Args:
        T: Temperature [K]
        
    Returns:
        Λ/Γ [cm^3]
    """
    T = np.maximum(T, T_MIN)
    term1 = 1.0e7 * np.exp(-114800.0 / (T + 1000.0))
    term2 = 1.4e-2 * np.sqrt(T) * np.exp(-92.0 / T)
    return term1 + term2


def cooling_coefficient(T):
    """
    Cooling coefficient Λ(T)
    
    Args:
        T: Temperature [K]
        
    Returns:
        Λ [erg cm^3 s^-1]
    """
    return Gamma_0 * cooling_coefficient_ratio(T)


def net_cooling_rate(n_H, T):
    """
    Net cooling rate per H nucleus
    
    Args:
        n_H: Number density [cm^-3]
        T: Temperature [K]
        
    Returns:
        L [erg s^-1] (positive = cooling, negative = heating)
    """
    Lambda = cooling_coefficient(T)
    return -Gamma_0 + n_H * Lambda


def volumetric_cooling_rate(n_H, T):
    """
    Volumetric cooling rate
    
    Args:
        n_H: Number density [cm^-3]
        T: Temperature [K]
        
    Returns:
        ρ_n L [erg cm^-3 s^-1]
    """
    return n_H * net_cooling_rate(n_H, T)


def equilibrium_temperature(n_H, T_guess=None, tol=1e-8, max_iter=100):
    """
    Solve for thermal equilibrium temperature at given density
    
    Solves: Γ = n Λ(T)
    
    Args:
        n_H: Number density [cm^-3]
        T_guess: Initial temperature guess [K]
        tol: Relative tolerance
        max_iter: Maximum iterations
        
    Returns:
        T_eq [K]
    """
    if n_H <= 0.0:
        return T_MAX
    
    # Initial guess based on typical ISM phases
    if T_guess is None:
        if n_H < 0.1:
            T_guess = 8000.0  # WNM
        elif n_H < 10.0:
            T_guess = 1000.0  # Transition
        else:
            T_guess = 100.0   # CNM
    
    T = T_guess
    T_lo = T_MIN
    T_hi = T_MAX
    
    for iteration in range(max_iter):
        ratio = cooling_coefficient_ratio(T)
        f = n_H * ratio - 1.0
        
        # Check convergence
        if abs(f) < tol:
            return T
        
        # Numerical derivative
        dT = T * 1e-6
        ratio_plus = cooling_coefficient_ratio(T + dT)
        df = n_H * (ratio_plus - ratio) / dT
        
        # Newton step with bounds
        T_new = T - f / df
        
        # Bisection fallback
        if T_new < T_lo or T_new > T_hi or not np.isfinite(T_new):
            T_new = np.sqrt(T_lo * T_hi)
        
        # Update bounds
        if f > 0:
            T_hi = T
        else:
            T_lo = T
        
        T = T_new
    
    return T


def is_thermally_unstable(n_H, T):
    """
    Check thermal stability using Balbus criterion
    
    Args:
        n_H: Number density [cm^-3]
        T: Temperature [K]
        
    Returns:
        bool: True if thermally unstable (∂(L/T)/∂T|_P < 0)
    """
    Lambda = cooling_coefficient(T)
    
    # Numerical derivative dΛ/dT
    dT = T * 1e-6
    Lambda_plus = cooling_coefficient(T + dT)
    dLambda_dT = (Lambda_plus - Lambda) / dT
    
    T2 = T * T
    deriv = Gamma_0 / T2 - 2.0 * n_H * Lambda / T2 + (n_H / T) * dLambda_dT
    
    return deriv < 0.0


def cooling_timescale(n_H, T):
    """
    Cooling/heating timescale
    
    Args:
        n_H: Number density [cm^-3]
        T: Temperature [K]
        
    Returns:
        t_cool [s]
    """
    L = np.abs(net_cooling_rate(n_H, T))
    L = np.maximum(L, 1e-40)
    return k_B * T / (m_n * L * (gamma - 1.0))


def plot_cooling_coefficients(output_dir, dpi=150):
    """Plot cooling/heating coefficients vs temperature"""
    T = np.logspace(np.log10(T_MIN), np.log10(T_MAX), 1000)
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Panel 1: Λ/Γ ratio
    ax = axes[0]
    ratio = cooling_coefficient_ratio(T)
    ax.loglog(T, ratio, 'b-', linewidth=2)
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel(r'$\Lambda / \Gamma$ [cm$^3$]')
    ax.set_title('Inoue & Inutsuka (2008) Cooling Coefficient Ratio')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label=r'$\Lambda/\Gamma = 1$')
    ax.legend()
    
    # Panel 2: Λ(T)
    ax = axes[1]
    Lambda = cooling_coefficient(T)
    ax.loglog(T, Lambda, 'r-', linewidth=2)
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    ax.set_title('Cooling Coefficient')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=Gamma_0, color='gray', linestyle='--', alpha=0.5, label=r'$\Gamma_0 = 2 \times 10^{-26}$ erg s$^{-1}$')
    ax.legend()
    
    plt.tight_layout()
    output_path = output_dir / 'cooling_coefficients.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()


def plot_net_cooling_rate(output_dir, dpi=150):
    """Plot net cooling rate vs temperature for different densities"""
    T = np.logspace(np.log10(T_MIN), 4.5, 1000)
    densities = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]  # cm^-3
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(densities)))
    
    for n_H, color in zip(densities, colors):
        L = net_cooling_rate(n_H, T)
        label = f'n = {n_H:.1f} cm$^{{-3}}$'
        ax.plot(T, L, linewidth=2, color=color, label=label)
        
        # Mark equilibrium point
        T_eq = equilibrium_temperature(n_H)
        ax.plot(T_eq, 0, 'o', color=color, markersize=8, markeredgecolor='white', markeredgewidth=1.5)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('Temperature [K]', fontsize=12)
    ax.set_ylabel(r'Net Cooling Rate $L$ [erg s$^{-1}$]', fontsize=12)
    ax.set_title('Net Cooling Rate per H Nucleus (Inoue & Inutsuka 2008)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=10)
    ax.text(0.98, 0.95, 'Cooling (L > 0)', transform=ax.transAxes, 
            ha='right', va='top', fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.text(0.98, 0.05, 'Heating (L < 0)', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=11, bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.tight_layout()
    output_path = output_dir / 'net_cooling_rate.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()


def plot_equilibrium_curve(output_dir, dpi=150):
    """Plot thermal equilibrium curve and pressure"""
    n_H = np.logspace(-2, 3, 200)  # cm^-3
    
    T_eq = np.array([equilibrium_temperature(n) for n in n_H])
    P_eq = n_H * T_eq  # P/k_B in K cm^-3
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: Temperature vs density
    ax = axes[0]
    ax.loglog(n_H, T_eq, 'b-', linewidth=3)
    ax.set_xlabel('Number Density [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Equilibrium Temperature [K]', fontsize=12)
    ax.set_title('Thermal Equilibrium Curve', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Mark typical ISM phases
    ax.axhline(y=6000, color='orange', linestyle='--', alpha=0.5, label='WNM (~6000 K)')
    ax.axhline(y=100, color='cyan', linestyle='--', alpha=0.5, label='CNM (~100 K)')
    ax.legend(fontsize=10)
    
    # Panel 2: Pressure vs density
    ax = axes[1]
    ax.loglog(n_H, P_eq, 'r-', linewidth=3)
    ax.set_xlabel('Number Density [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel(r'Equilibrium Pressure $P/k_B$ [K cm$^{-3}$]', fontsize=12)
    ax.set_title('Thermal Pressure Equilibrium', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Mark typical ISM phases
    ax.axhline(y=3500, color='orange', linestyle='--', alpha=0.5, label='WNM (P/k ~ 3500)')
    ax.axhline(y=3000, color='cyan', linestyle='--', alpha=0.5, label='CNM (P/k ~ 3000)')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    output_path = output_dir / 'equilibrium_curve.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()


def plot_phase_diagram(output_dir, dpi=150):
    """Plot density-temperature phase diagram with stability"""
    n_grid = np.logspace(-2, 3, 150)
    T_grid = np.logspace(1, 5, 150)
    
    N, T = np.meshgrid(n_grid, T_grid)
    
    # Calculate net cooling rate
    L = net_cooling_rate(N, T)
    
    # Calculate stability
    stability = np.zeros_like(N)
    for i in range(len(T_grid)):
        for j in range(len(n_grid)):
            stability[i, j] = 1.0 if is_thermally_unstable(N[i, j], T[i, j]) else 0.0
    
    # Compute equilibrium curve
    T_eq = np.array([equilibrium_temperature(n) for n in n_grid])
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Panel 1: Net cooling rate
    ax = axes[0]
    levels = np.logspace(-30, -22, 20)
    levels = np.concatenate([-levels[::-1], [0], levels])
    
    cf = ax.contourf(N, T, L, levels=levels, cmap='RdBu_r', norm=plt.Normalize(vmin=-1e-23, vmax=1e-23), extend='both')
    ax.contour(N, T, L, levels=[0], colors='black', linewidths=2)
    ax.plot(n_grid, T_eq, 'k--', linewidth=3, label='Thermal Equilibrium')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number Density [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Temperature [K]', fontsize=12)
    ax.set_title('Net Cooling Rate per H Nucleus', fontsize=14)
    ax.legend(fontsize=10)
    
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label(r'$L$ [erg s$^{-1}$]', fontsize=11)
    
    # Panel 2: Thermal stability
    ax = axes[1]
    im = ax.contourf(N, T, stability, levels=[0, 0.5, 1.0], colors=['lightblue', 'salmon'], alpha=0.7)
    ax.plot(n_grid, T_eq, 'k-', linewidth=3, label='Thermal Equilibrium')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number Density [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Temperature [K]', fontsize=12)
    ax.set_title('Thermal Stability (Balbus Criterion)', fontsize=14)
    
    # Custom legend for stability
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', alpha=0.7, label='Stable'),
        Patch(facecolor='salmon', alpha=0.7, label='Unstable'),
        plt.Line2D([0], [0], color='k', linewidth=3, label='Equilibrium')
    ]
    ax.legend(handles=legend_elements, fontsize=10)
    
    plt.tight_layout()
    output_path = output_dir / 'phase_diagram.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()


def plot_cooling_timescale(output_dir, dpi=150):
    """Plot cooling timescale as function of density and temperature"""
    n_grid = np.logspace(-2, 3, 150)
    T_grid = np.logspace(1.5, 4.5, 150)
    
    N, T = np.meshgrid(n_grid, T_grid)
    
    # Calculate cooling timescale
    t_cool = cooling_timescale(N, T)
    t_cool_yr = t_cool / (365.25 * 24 * 3600)  # Convert to years
    
    # Compute equilibrium curve
    T_eq = np.array([equilibrium_temperature(n) for n in n_grid])
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot cooling timescale
    levels = np.logspace(-2, 8, 50)
    cf = ax.contourf(N, T, t_cool_yr, levels=levels, cmap='viridis', norm=LogNorm(), extend='both')
    
    # Contour lines
    contour_levels = [1e-1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6]
    cs = ax.contour(N, T, t_cool_yr, levels=contour_levels, colors='white', linewidths=0.8, alpha=0.5)
    ax.clabel(cs, inline=True, fontsize=9, fmt='%g yr')
    
    # Equilibrium curve
    ax.plot(n_grid, T_eq, 'r--', linewidth=3, label='Thermal Equilibrium')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number Density [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Temperature [K]', fontsize=12)
    ax.set_title('Cooling Timescale (Inoue & Inutsuka 2008)', fontsize=14)
    ax.legend(fontsize=11)
    
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label(r'$t_{\rm cool}$ [yr]', fontsize=12)
    
    plt.tight_layout()
    output_path = output_dir / 'cooling_timescale.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()


def plot_summary_panel(output_dir, dpi=150):
    """Create a comprehensive 4-panel summary figure"""
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Panel 1: Cooling coefficient ratio
    ax1 = fig.add_subplot(gs[0, 0])
    T = np.logspace(np.log10(T_MIN), 4.5, 1000)
    ratio = cooling_coefficient_ratio(T)
    ax1.loglog(T, ratio, 'b-', linewidth=2)
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Temperature [K]', fontsize=11)
    ax1.set_ylabel(r'$\Lambda / \Gamma$ [cm$^3$]', fontsize=11)
    ax1.set_title('(a) Cooling Coefficient Ratio', fontsize=12)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Net cooling rate for selected densities
    ax2 = fig.add_subplot(gs[0, 1])
    densities = [0.3, 1.0, 3.0, 10.0, 30.0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(densities)))
    
    for n_H, color in zip(densities, colors):
        L = net_cooling_rate(n_H, T)
        ax2.plot(T, L, linewidth=2, color=color, label=f'n = {n_H:.1f}')
        T_eq = equilibrium_temperature(n_H)
        ax2.plot(T_eq, 0, 'o', color=color, markersize=6, markeredgecolor='white', markeredgewidth=1)
    
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_xlabel('Temperature [K]', fontsize=11)
    ax2.set_ylabel(r'$L$ [erg s$^{-1}$]', fontsize=11)
    ax2.set_title('(b) Net Cooling Rate', fontsize=12)
    ax2.legend(fontsize=9, ncol=2)
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Equilibrium curve
    ax3 = fig.add_subplot(gs[1, 0])
    n_H = np.logspace(-2, 3, 200)
    T_eq = np.array([equilibrium_temperature(n) for n in n_H])
    ax3.loglog(n_H, T_eq, 'r-', linewidth=3)
    ax3.axhline(y=6000, color='orange', linestyle='--', alpha=0.5, linewidth=1.5, label='WNM')
    ax3.axhline(y=100, color='cyan', linestyle='--', alpha=0.5, linewidth=1.5, label='CNM')
    ax3.set_xlabel('Number Density [cm$^{-3}$]', fontsize=11)
    ax3.set_ylabel('Equilibrium Temperature [K]', fontsize=11)
    ax3.set_title('(c) Thermal Equilibrium', fontsize=12)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Phase diagram with stability
    ax4 = fig.add_subplot(gs[1, 1])
    n_grid_phase = np.logspace(-2, 3, 100)
    T_grid = np.logspace(1, 5, 100)
    N, T_mesh = np.meshgrid(n_grid_phase, T_grid)
    
    stability = np.zeros_like(N)
    for i in range(len(T_grid)):
        for j in range(len(n_grid_phase)):
            stability[i, j] = 1.0 if is_thermally_unstable(N[i, j], T_mesh[i, j]) else 0.0
    
    # Compute equilibrium curve for this grid
    T_eq_phase = np.array([equilibrium_temperature(n) for n in n_grid_phase])
    
    ax4.contourf(N, T_mesh, stability, levels=[0, 0.5, 1.0], colors=['lightblue', 'salmon'], alpha=0.7)
    ax4.plot(n_grid_phase, T_eq_phase, 'k-', linewidth=3, label='Equilibrium')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Number Density [cm$^{-3}$]', fontsize=11)
    ax4.set_ylabel('Temperature [K]', fontsize=11)
    ax4.set_title('(d) Thermal Stability', fontsize=12)
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', alpha=0.7, label='Stable'),
        Patch(facecolor='salmon', alpha=0.7, label='Unstable'),
        plt.Line2D([0], [0], color='k', linewidth=3, label='Equilibrium')
    ]
    ax4.legend(handles=legend_elements, fontsize=9)
    
    fig.suptitle('Inoue & Inutsuka (2008) ISM Cooling Function', fontsize=16, fontweight='bold', y=0.995)
    
    output_path = output_dir / 'cooling_summary.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Plot Inoue & Inutsuka (2008) ISM cooling function properties',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Generate all plots in default directory
    python scripts/plot_inoue_inutsuka_cooling.py
    
    # Specify output directory
    python scripts/plot_inoue_inutsuka_cooling.py --output-dir results/cooling_plots
    
    # High resolution for publication
    python scripts/plot_inoue_inutsuka_cooling.py --dpi 300
        """
    )
    
    parser.add_argument('--output-dir', type=Path, default=Path('results/cooling_analysis'),
                        help='Output directory for plots (default: results/cooling_analysis)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='DPI for output images (default: 150)')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Inoue & Inutsuka (2008) ISM Cooling Function Visualization")
    print("=" * 70)
    print(f"Output directory: {args.output_dir}")
    print(f"DPI: {args.dpi}")
    print()
    
    # Generate plots
    print("Generating plots...")
    print()
    
    plot_cooling_coefficients(args.output_dir, args.dpi)
    plot_net_cooling_rate(args.output_dir, args.dpi)
    plot_equilibrium_curve(args.output_dir, args.dpi)
    plot_phase_diagram(args.output_dir, args.dpi)
    plot_cooling_timescale(args.output_dir, args.dpi)
    plot_summary_panel(args.output_dir, args.dpi)
    
    print()
    print("=" * 70)
    print("✓ All plots generated successfully!")
    print("=" * 70)
    print(f"\nPlots saved to: {args.output_dir.absolute()}")
    print("\nGenerated files:")
    print("  - cooling_coefficients.png    : Λ/Γ and Λ(T)")
    print("  - net_cooling_rate.png        : L vs T for various densities")
    print("  - equilibrium_curve.png       : T_eq and P_eq vs density")
    print("  - phase_diagram.png           : Cooling rate and stability map")
    print("  - cooling_timescale.png       : t_cool in (n,T) space")
    print("  - cooling_summary.png         : Comprehensive 4-panel summary")


if __name__ == '__main__':
    main()
