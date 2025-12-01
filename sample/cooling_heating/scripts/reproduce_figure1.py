#!/usr/bin/env python3
"""
Reproduce Figure 1 from Koyama & Inutsuka (2000)

Figure 1 shows equilibrium thermal and chemical state of ISM:
(a) Temperature and pressure vs density
(b) Chemical fractions vs density
(c) Heating and cooling rates vs density  
(d) Timescales vs density

Two cases: N_H_col = 1e19 and 1e20 cm^-2 (shown as solid and dashed lines)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from thermal_equilibrium import ThermalEquilibriumSolver
from chemistry_network import ChemistryNetwork, k_B

def main():
    """Generate Figure 1."""
    
    print("=" * 70)
    print("Reproducing Figure 1 from Koyama & Inutsuka (2000)")
    print("Full chemistry network implementation")
    print("=" * 70)
    
    # Density range: 0.1 to 1e6 cm^-3
    n_array = np.logspace(-1, 6, 100)
    
    # Two column density cases
    N_col_1 = 1e19  # cm^-2 (solid lines)
    N_col_2 = 1e20  # cm^-2 (dashed lines)
    
    print(f"\nComputing equilibrium for N_H = {N_col_1:.0e} cm^-2 (solid lines)...")
    solver_1 = ThermalEquilibriumSolver(G0=1.7, N_H_col=N_col_1)
    
    # For simplicity, assume N_HI and N_H2 scale with local conditions
    # In reality, these should be integrated column densities
    # For first pass, use simple estimates
    N_HI_1 = N_col_1 * np.ones_like(n_array)
    N_H2_1 = np.zeros_like(n_array)
    
    results_1 = solver_1.compute_equilibrium_curve(n_array, N_HI_1, N_H2_1)
    
    print(f"\nComputing equilibrium for N_H = {N_col_2:.0e} cm^-2 (dashed lines)...")
    solver_2 = ThermalEquilibriumSolver(G0=1.7, N_H_col=N_col_2)
    
    N_HI_2 = N_col_2 * np.ones_like(n_array)
    N_H2_2 = np.zeros_like(n_array)
    
    results_2 = solver_2.compute_equilibrium_curve(n_array, N_HI_2, N_H2_2)
    
    # Compute heating/cooling rates for panel (c)
    print("\nComputing heating/cooling rates...")
    heating_1, cooling_1 = compute_rates(solver_1, results_1)
    heating_2, cooling_2 = compute_rates(solver_2, results_2)
    
    # Compute timescales for panel (d)
    print("\nComputing timescales...")
    timescales_1 = compute_timescales_array(solver_1, results_1)
    
    # Create the four-panel figure
    print("\nGenerating plots...")
    fig = plt.figure(figsize=(12, 10))
    
    # Panel (a): Temperature and Pressure
    ax1 = fig.add_subplot(2, 2, 1)
    plot_panel_a(ax1, results_1, results_2)
    
    # Panel (b): Chemical fractions
    ax2 = fig.add_subplot(2, 2, 2)
    plot_panel_b(ax2, results_1, results_2)
    
    # Panel (c): Heating and cooling rates
    ax3 = fig.add_subplot(2, 2, 3)
    plot_panel_c(ax3, heating_1, cooling_1, heating_2, cooling_2, n_array)
    
    # Panel (d): Timescales
    ax4 = fig.add_subplot(2, 2, 4)
    plot_panel_d(ax4, timescales_1, n_array)
    
    plt.tight_layout()
    
    # Save figure
    output_file = '../results/figure1_reproduction.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved to: {output_file}")
    
    # Also save individual panels as PostScript for comparison
    save_individual_panels(results_1, results_2, heating_1, cooling_1, 
                          heating_2, cooling_2, timescales_1, n_array)
    
    print("\n" + "=" * 70)
    print("Figure 1 reproduction complete!")
    print("=" * 70)

def compute_rates(solver, results):
    """Compute heating and cooling rates at equilibrium."""
    n_array = results['n']
    n_points = len(n_array)
    
    heating = {
        'PE': np.zeros(n_points),
        'CR': np.zeros(n_points),
        'XR': np.zeros(n_points),
        'H2': np.zeros(n_points),
        'total': np.zeros(n_points)
    }
    
    cooling = {
        'CII': np.zeros(n_points),
        'OI': np.zeros(n_points),
        'Lya': np.zeros(n_points),
        'H2': np.zeros(n_points),
        'CO': np.zeros(n_points),
        'GR': np.zeros(n_points),
        'total': np.zeros(n_points)
    }
    
    for i in range(n_points):
        n = results['n'][i]
        T = results['T'][i]
        x_e = results['x_e'][i]
        x_2 = results['x_2'][i]
        x_CO = results['x_CO'][i]
        N_HI = results['N_HI'][i]
        N_H2 = results['N_H2'][i]
        
        _, h_dict, c_dict = solver.chem.net_heating_cooling(
            n, T, x_e, x_2, x_CO, N_HI, N_H2
        )
        
        heating['PE'][i] = h_dict['PE']
        heating['CR'][i] = h_dict['CR']
        heating['XR'][i] = h_dict['XR']
        heating['H2'][i] = h_dict['H2_form'] + h_dict['H2_diss']
        heating['total'][i] = h_dict['total']
        
        # Cooling (combine components for plotting)
        cooling['CII'][i] = c_dict['CII']
        cooling['OI'][i] = c_dict['OI']
        cooling['Lya'][i] = c_dict['Lya']
        cooling['H2'][i] = c_dict['H2']
        cooling['CO'][i] = c_dict['CO_rot'] + c_dict['CO_vib']
        cooling['GR'][i] = c_dict['GR']
        cooling['total'][i] = c_dict['total']
    
    return heating, cooling

def compute_timescales_array(solver, results):
    """Compute timescales for all points."""
    n_array = results['n']
    n_points = len(n_array)
    
    timescales = {
        't_cool': np.zeros(n_points),
        't_rec': np.zeros(n_points),
        't_ff': np.zeros(n_points),
        't_h2': np.zeros(n_points)
    }
    
    for i in range(n_points):
        n = results['n'][i]
        T = results['T'][i]
        x_e = results['x_e'][i]
        x_2 = results['x_2'][i]
        x_CO = results['x_CO'][i]
        
        t = solver.compute_timescales(n, T, x_e, x_2, x_CO)
        
        timescales['t_cool'][i] = t['t_cool']
        timescales['t_rec'][i] = t['t_rec']
        timescales['t_ff'][i] = t['t_ff']
        timescales['t_h2'][i] = t['t_h2']
    
    return timescales

def plot_panel_a(ax, res1, res2):
    """Panel (a): Temperature and Pressure."""
    
    # Temperature
    ax.loglog(res1['n'], res1['T'], 'b-', linewidth=2, label=r'$T$ ($N_H=10^{19}$)')
    ax.loglog(res2['n'], res2['T'], 'b--', linewidth=2, label=r'$T$ ($N_H=10^{20}$)')
    
    # Pressure (scaled for visibility)
    P_scale = 1e14  # Scale to match temperature range roughly
    ax.loglog(res1['n'], res1['P'] / P_scale, 'r-', linewidth=2, 
              label=r'$P/10^{14}$ ($N_H=10^{19}$)')
    ax.loglog(res2['n'], res2['P'] / P_scale, 'r--', linewidth=2,
              label=r'$P/10^{14}$ ($N_H=10^{20}$)')
    
    ax.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel(r'$T$ [K], $P/10^{14}$ [dyne cm$^{-2}$]', fontsize=12)
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 1e5)
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_title('(a) Equilibrium Temperature and Pressure', fontsize=12, fontweight='bold')

def plot_panel_b(ax, res1, res2):
    """Panel (b): Chemical fractions (using res1, similar for both column densities)."""
    
    ax.loglog(res1['n'], res1['x_e'], 'r-', linewidth=2, label=r'$x_e$')
    ax.loglog(res1['n'], res1['x_2'], 'b-', linewidth=2, label=r'$x_2$ (H$_2$)')
    ax.loglog(res1['n'], res1['x_CO'], 'g-', linewidth=2, label=r'$x_{CO}$')
    ax.loglog(res1['n'], res1['x_HI'], 'k--', linewidth=2, label=r'$x_{HI}$')
    
    ax.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Chemical Fraction', fontsize=12)
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-10, 2)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_title('(b) Chemical Fractions', fontsize=12, fontweight='bold')

def plot_panel_c(ax, h1, c1, h2, c2, n_array):
    """Panel (c): Heating and Cooling rates."""
    
    # Solid lines: N_H = 1e19
    # Heating (dashed)
    ax.loglog(n_array, h1['PE'], 'r:', linewidth=1.5, label='PE (heat)')
    ax.loglog(n_array, h1['CR'], 'b:', linewidth=1.5, label='CR (heat)')
    ax.loglog(n_array, h1['XR'], 'g:', linewidth=1.5, label='XR (heat)')
    ax.loglog(n_array, h1['H2'], 'm:', linewidth=1.5, label='H$_2$ (heat)')
    
    # Cooling (solid)
    ax.loglog(n_array, c1['CII'], 'r-', linewidth=2, label='CII (cool)')
    ax.loglog(n_array, c1['OI'], 'b-', linewidth=2, label='OI (cool)')
    ax.loglog(n_array, c1['Lya'], 'g-', linewidth=2, label=r'Ly-$\alpha$ (cool)')
    ax.loglog(n_array, c1['H2'], 'm-', linewidth=2, label='H$_2$ (cool)')
    ax.loglog(n_array, c1['CO'], 'c-', linewidth=2, label='CO (cool)')
    ax.loglog(n_array, c1['GR'], 'k-', linewidth=2, label='GR (cool)')
    
    ax.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel(r'Rate [erg s$^{-1}$ H$^{-1}$]', fontsize=12)
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-30, 1e-20)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_title('(c) Heating and Cooling Rates', fontsize=12, fontweight='bold')

def plot_panel_d(ax, timescales, n_array):
    """Panel (d): Timescales."""
    
    ax.loglog(n_array, timescales['t_cool'], 'r-', linewidth=2, label=r'$t_{cool}$')
    ax.loglog(n_array, timescales['t_rec'], 'b--', linewidth=2, label=r'$t_{rec}$')
    ax.loglog(n_array, timescales['t_ff'], 'g-.', linewidth=2, label=r'$t_{ff}$')
    ax.loglog(n_array, timescales['t_h2'], 'm:', linewidth=2, label=r'$t_{H_2}$')
    
    ax.set_xlabel(r'$n$ [cm$^{-3}$]', fontsize=12)
    ax.set_ylabel('Timescale [years]', fontsize=12)
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e0, 1e12)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_title('(d) Physical Timescales', fontsize=12, fontweight='bold')

def save_individual_panels(res1, res2, h1, c1, h2, c2, t1, n_array):
    """Save individual panels as separate files for detailed comparison."""
    
    # Panel a
    fig_a, ax_a = plt.subplots(figsize=(6, 5))
    plot_panel_a(ax_a, res1, res2)
    fig_a.savefig('../results/f1a_reproduction.png', dpi=300, bbox_inches='tight')
    plt.close(fig_a)
    
    # Panel b
    fig_b, ax_b = plt.subplots(figsize=(6, 5))
    plot_panel_b(ax_b, res1, res2)
    fig_b.savefig('../results/f1b_reproduction.png', dpi=300, bbox_inches='tight')
    plt.close(fig_b)
    
    # Panel c
    fig_c, ax_c = plt.subplots(figsize=(6, 5))
    plot_panel_c(ax_c, h1, c1, h2, c2, n_array)
    fig_c.savefig('../results/f1c_reproduction.png', dpi=300, bbox_inches='tight')
    plt.close(fig_c)
    
    # Panel d
    fig_d, ax_d = plt.subplots(figsize=(6, 5))
    plot_panel_d(ax_d, t1, n_array)
    fig_d.savefig('../results/f1d_reproduction.png', dpi=300, bbox_inches='tight')
    plt.close(fig_d)
    
    print("✓ Individual panels saved as PNG files")

if __name__ == '__main__':
    main()
