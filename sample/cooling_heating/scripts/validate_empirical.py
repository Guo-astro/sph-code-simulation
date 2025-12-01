#!/usr/bin/env python3
"""
Validate Empirical T(n) Curve
Checks that the digitized data shows proper thermal bistability
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from empirical_equilibrium import EmpiricalEquilibrium

def main():
    print("=" * 70)
    print("Validating Empirical T(n) Curve")
    print("=" * 70)
    
    solver = EmpiricalEquilibrium(G0=1.7, N_H_col=1e19)
    
    # Test specific densities across all regimes
    test_densities = [0.1, 0.5, 1.0, 5.0, 10, 50, 100, 1000, 10000, 1e5, 1e6]
    
    print("\nTemperature vs Density:")
    print(f"{'n (cm^-3)':<12} {'T (K)':<10} {'P/k_B (K/cm^3)':<18} Regime")
    print("-" * 70)
    
    for n in test_densities:
        T = solver.get_equilibrium_T(n, column_density=1e19)
        P = solver.get_equilibrium_P(n, column_density=1e19)
        
        # Determine regime
        if T > 5000:
            regime = "WNM (Warm)"
        elif T > 1000:
            regime = "Transition"
        elif T > 100:
            regime = "CNM (Cold)"
        else:
            regime = "Molecular"
        
        print(f"{n:<12.2e} {T:<10.1f} {P/1.38e-16:<18.2e} {regime}")
    
    # Check thermal bistability
    print("\n" + "=" * 70)
    print("Checking Thermal Bistability:")
    print("=" * 70)
    
    # At n ~ 1-10 cm^-3, we should see S-curve structure
    # WNM: T ~ 7000-8000 K
    # CNM: T ~ 100 K
    
    T_at_1 = solver.get_equilibrium_T(1.0, column_density=1e19)
    T_at_10 = solver.get_equilibrium_T(10.0, column_density=1e19)
    
    print(f"T(n=1 cm^-3)  = {T_at_1:.0f} K")
    print(f"T(n=10 cm^-3) = {T_at_10:.0f} K")
    
    if 6000 < T_at_1 < 9000:
        print("✓ WNM branch present at n=1")
    else:
        print(f"✗ WNM temperature wrong: expected 6000-9000K, got {T_at_1:.0f}K")
    
    if 80 < T_at_10 < 200:
        print("✓ CNM branch present at n=10")
    else:
        print(f"✗ CNM temperature wrong: expected 80-200K, got {T_at_10:.0f}K")
    
    # Molecular regime
    T_at_1000 = solver.get_equilibrium_T(1000, column_density=1e19)
    if 30 < T_at_1000 < 70:
        print(f"✓ Molecular regime at n=1000: T={T_at_1000:.0f}K")
    else:
        print(f"✗ Molecular T wrong: expected 30-70K, got {T_at_1000:.0f}K")
    
    # Plot to visualize
    print("\nGenerating validation plot...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Full curve
    n_array = np.logspace(-1, 6, 200)
    results = solver.compute_equilibrium_curve(n_array, column_density=1e19)
    
    ax1.loglog(results['n'], results['T'], 'b-', linewidth=3, label='Empirical T(n)')
    ax1.scatter(test_densities, 
               [solver.get_equilibrium_T(n, 1e19) for n in test_densities],
               s=100, c='red', zorder=5, label='Test points')
    ax1.axhline(8000, color='green', linestyle='--', alpha=0.5, label='WNM ~8000K')
    ax1.axhline(100, color='orange', linestyle='--', alpha=0.5, label='CNM ~100K')
    ax1.axhline(50, color='purple', linestyle='--', alpha=0.5, label='Molecular ~50K')
    ax1.set_xlabel('n [cm$^{-3}$]', fontsize=13)
    ax1.set_ylabel('T [K]', fontsize=13)
    ax1.set_title('Temperature vs Density', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0.1, 1e6)
    ax1.set_ylim(10, 1e5)
    
    # Pressure
    ax2.loglog(results['n'], results['P']/1.38e-16, 'r-', linewidth=3)
    ax2.set_xlabel('n [cm$^{-3}$]', fontsize=13)
    ax2.set_ylabel('P/k$_B$ [K cm$^{-3}$]', fontsize=13)
    ax2.set_title('Pressure vs Density', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.1, 1e6)
    
    plt.tight_layout()
    plt.savefig('../results/empirical_validation.png', dpi=200, bbox_inches='tight')
    print("✓ Validation plot saved: ../results/empirical_validation.png")
    
    print("\n" + "=" * 70)
    print("✓ Empirical curve validation complete!")
    print("=" * 70)

if __name__ == '__main__':
    main()
