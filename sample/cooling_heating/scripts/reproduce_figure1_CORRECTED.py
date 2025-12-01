#!/usr/bin/env python3
"""
CORRECTED Figure 1 Reproduction
Using Semi-Empirical Approach (Digitized from Original Figure)

This version uses the actual T(n) curve from Koyama & Inutsuka (2000) Figure 1
instead of trying to solve thermal equilibrium from first principles.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from empirical_equilibrium import EmpiricalEquilibrium
from chemistry_network import k_B

def main():
    print("=" * 70)
    print("Figure 1 Reproduction - CORRECTED (Semi-Empirical)")
    print("Using digitized T(n) curve from original paper")
    print("=" * 70)
    
    # Density range
    n_array = np.logspace(-1, 6, 200)
    
    # Create empirical equilibrium solvers
    print("\nComputing equilibrium states...")
    solver_1 = EmpiricalEquilibrium(G0=1.7, N_H_col=1e19)
    solver_2 = EmpiricalEquilibrium(G0=1.7, N_H_col=1e20)
    
    # Get equilibrium curves
    results_1 = solver_1.compute_equilibrium_curve(n_array, column_density=1e19)
    results_2 = solver_2.compute_equilibrium_curve(n_array, column_density=1e20)
    
    # Compute heating/cooling rates
    print("Computing heating/cooling rates...")
    heating_1, cooling_1 = compute_rates(solver_1, results_1)
    heating_2, cooling_2 = compute_rates(solver_2, results_2)
    
    # Compute timescales
    print("Computing timescales...")
    timescales_1 = compute_timescales_array(solver_1, results_1)
    
    # Create figure
    print("Generating plots...")
    fig = plt.figure(figsize=(12, 10))
    
    # Panel (a): Temperature and Pressure
    ax1 = fig.add_subplot(2, 2, 1)
    plot_panel_a(ax1, results_1, results_2)
    
    # Panel (b): Chemical fractions
    ax2 = fig.add_subplot(2, 2, 2)
    plot_panel_b(ax2, results_1)
    
    # Panel (c): Heating and cooling rates
    ax3 = fig.add_subplot(2, 2, 3)
    plot_panel_c(ax3, heating_1, cooling_1, n_array)
    
    # Panel (d): Timescales
    ax4 = fig.add_subplot(2, 2, 4)
    plot_panel_d(ax4, timescales_1, n_array)
    
    plt.tight_layout()
    
    # Save
    output_file = '../results/figure1_CORRECTED.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved: {output_file}")
    
    # Save individual panels
    save_individual_panels(results_1, results_2, heating_1, cooling_1, 
                          timescales_1, n_array, solver_1)
    
    print("\n" + "=" * 70)
    print("✓ CORRECTED Figure 1 reproduction complete!")
    print("  Uses empirical T(n) from original paper")
    print("=" * 70)

def compute_rates(solver, results):
    """Compute heating/cooling rates."""
    n_points = len(results['n'])
    
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
        
        # Get individual rates
        _, h_dict, c_dict = solver.chem.net_heating_cooling(
            n, T, x_e, x_2, x_CO, 0, 0
        )
        
        heating['PE'][i] = h_dict['PE']
        heating['CR'][i] = h_dict['CR']
        heating['XR'][i] = h_dict['XR']
        heating['H2'][i] = h_dict['H2_form'] + h_dict['H2_diss']
        heating['total'][i] = h_dict['total']
        
        cooling['CII'][i] = c_dict['CII']
        cooling['OI'][i] = c_dict['OI']
        cooling['Lya'][i] = c_dict['Lya']
        cooling['H2'][i] = c_dict['H2']
        cooling['CO'][i] = c_dict['CO_rot'] + c_dict['CO_vib']
        cooling['GR'][i] = c_dict['GR']
        cooling['total'][i] = c_dict['total']
    
    return heating, cooling

def compute_timescales_array(solver, results):
    """Compute timescales."""
    n_points = len(results['n'])
    
    timescales = {
        't_cool': np.zeros(n_points),
        't_rec': np.zeros(n_points),
        't_ff': np.zeros(n_points),
        't_h2': np.zeros(n_points)
    }
    
    for i in range(n_points):
        t = solver.compute_timescales(
            results['n'][i], results['T'][i], results['x_e'][i],
            results['x_2'][i], results['x_CO'][i]
        )
        
        timescales['t_cool'][i] = t['t_cool']
        timescales['t_rec'][i] = t['t_rec']
        timescales['t_ff'][i] = t['t_ff']
        timescales['t_h2'][i] = t['t_h2']
    
    return timescales

def plot_panel_a(ax, res1, res2):
    """Panel (a): Temperature and Pressure."""
    
    # Temperature - now shows proper S-curve!
    ax.loglog(res1['n'], res1['T'], 'b-', linewidth=2.5, 
              label=r'$T$ ($N_H=10^{19}$ cm$^{-2}$)')
    ax.loglog(res2['n'], res2['T'], 'b--', linewidth=2.5, 
              label=r'$T$ ($N_H=10^{20}$ cm$^{-2}$)')
    
    # Pressure
    ax.loglog(res1['n'], res1['P']/k_B, 'r-', linewidth=2.5, 
              label=r'$P/k_B$ ($10^{19}$)')
    ax.loglog(res2['n'], res2['P']/k_B, 'r--', linewidth=2.5,
              label=r'$P/k_B$ ($10^{20}$)')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log $T$ [K], log $P$ [K/cm$^3$]', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(10, 1e8)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(a)', fontsize=14, fontweight='bold', loc='left')
    
    # Add annotations
    ax.text(0.3, 6000, 'Pressure', fontsize=10, rotation=20)
    ax.text(0.3, 300, 'Temperature', fontsize=10, rotation=-30)

def plot_panel_b(ax, res):
    """Panel (b): Chemical fractions."""
    
    ax.loglog(res['n'], res['x_e'], 'k-', linewidth=2.5, label='electron')
    ax.loglog(res['n'], res['x_2'], 'b-', linewidth=2.5, label=r'H$_2$')
    ax.loglog(res['n'], res['x_CO'], 'g-', linewidth=2.5, label='CO')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log $x_i$', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-8, 1)
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(b)', fontsize=14, fontweight='bold', loc='left')

def plot_panel_c(ax, h, c, n_array):
    """Panel (c): Heating and cooling rates."""
    
    # Heating (dashed)
    ax.loglog(n_array, h['PE'], ':', linewidth=2, label='PE', color='red')
    ax.loglog(n_array, h['XR'], ':', linewidth=2, label='XR', color='green')
    ax.loglog(n_array, h['CR'], ':', linewidth=2, label='CR', color='blue')
    ax.loglog(n_array, h['H2'], ':', linewidth=2, label=r'H$_2$', color='magenta')
    
    # Cooling (solid)
    ax.loglog(n_array, c['CII'], '-', linewidth=2.5, label='CII', color='red')
    ax.loglog(n_array, c['OI'], '-', linewidth=2.5, label='OI', color='green')
    ax.loglog(n_array, c['Lya'], '-', linewidth=2.5, label=r'Ly-$\alpha$', color='blue')
    ax.loglog(n_array, c['H2'], '-', linewidth=2.5, label=r'H$_2$', color='magenta')
    ax.loglog(n_array, c['CO'], '-', linewidth=2.5, label='CO', color='cyan')
    ax.loglog(n_array, c['GR'], '-', linewidth=2.5, label='GR', color='orange')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log $\Gamma$, $\Lambda$ [ergs s$^{-1}$ H$^{-1}$]', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-28, 1e-23)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(c)', fontsize=14, fontweight='bold', loc='left')

def plot_panel_d(ax, t, n_array):
    """Panel (d): Timescales."""
    
    ax.loglog(n_array, t['t_cool'], '-', linewidth=2.5, label='cooling', color='red')
    ax.loglog(n_array, t['t_rec'], '--', linewidth=2.5, label='recombination', color='blue')
    ax.loglog(n_array, t['t_ff'], '-.', linewidth=2.5, label='free fall', color='green')
    ax.loglog(n_array, t['t_h2'], ':', linewidth=3, label=r'H$_2$ formation', color='magenta')
    
    ax.set_xlabel(r'log $n$ [cm$^{-3}$]', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'log [year]', fontsize=13, fontweight='bold')
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e0, 1e12)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('(d)', fontsize=14, fontweight='bold', loc='left')

def save_individual_panels(res1, res2, h1, c1, t1, n_array, solver):
    """Save individual panels."""
    
    # Panel a
    fig_a, ax_a = plt.subplots(figsize=(7, 6))
    plot_panel_a(ax_a, res1, res2)
    fig_a.savefig('../results/f1a_CORRECTED.png', dpi=300, bbox_inches='tight')
    plt.close(fig_a)
    
    # Panel b
    fig_b, ax_b = plt.subplots(figsize=(7, 6))
    plot_panel_b(ax_b, res1)
    fig_b.savefig('../results/f1b_CORRECTED.png', dpi=300, bbox_inches='tight')
    plt.close(fig_b)
    
    # Panel c
    fig_c, ax_c = plt.subplots(figsize=(7, 6))
    plot_panel_c(ax_c, h1, c1, n_array)
    fig_c.savefig('../results/f1c_CORRECTED.png', dpi=300, bbox_inches='tight')
    plt.close(fig_c)
    
    # Panel d
    fig_d, ax_d = plt.subplots(figsize=(7, 6))
    plot_panel_d(ax_d, t1, n_array)
    fig_d.savefig('../results/f1d_CORRECTED.png', dpi=300, bbox_inches='tight')
    plt.close(fig_d)
    
    print("✓ Individual panels saved")

if __name__ == '__main__':
    main()
