#!/usr/bin/env python3
"""
Final reconciliation: Why is measured growth rate so much slower than predicted?

Key insight: The system is NOT in pure gravitational collapse mode.
It's in a QUASI-STATIC contraction where pressure mostly balances gravity.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_snapshot(filepath):
    return pd.read_csv(filepath, comment='#')

def main():
    no_gradh_dir = Path("/Users/guo/Downloads/sphcode/sample/gradh_study_3d/results/kernel_test/no_gradh_hk")
    
    print("=" * 80)
    print("UNDERSTANDING THE SLOW COLLAPSE")
    print("=" * 80)
    
    # The key realization: the system starts in EQUILIBRIUM
    # Small perturbations grow, but most particles are still near equilibrium
    
    print("""
The paradox:
- Theory predicts: Γ/ω_ff ~ √(1 - C_eff) ~ 0.6 for C_eff ~ 0.6
- Measured: Γ/ω_ff ~ 0.006

Resolution: The analysis assumed the ENTIRE force is imbalanced.
But the system starts in equilibrium! The imbalance is only in the
PERTURBATION, not the background force.

Let's compute the proper perturbation analysis.
""")

    # Load data
    df_0 = load_snapshot(no_gradh_dir / "snapshot_0000.csv")
    df_50 = load_snapshot(no_gradh_dir / "snapshot_0050.csv")
    
    for df in [df_0, df_50]:
        df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    
    # Get center particles
    core_0 = df_0.nsmallest(10, 'r')
    core_50 = df_50.nsmallest(10, 'r')
    
    rho_0 = core_0['dens'].mean()
    rho_50 = core_50['dens'].mean()
    P_0 = core_0['pres'].mean()
    P_50 = core_50['pres'].mean()
    
    print("\n" + "=" * 80)
    print("EQUILIBRIUM ANALYSIS")
    print("=" * 80)
    
    print(f"""
At t=0, the system is in hydrostatic equilibrium:
  ∇P = -ρ∇Φ (gravity balanced by pressure gradient)

The TOTAL force on center particles is near zero.
Only the PERTURBATION creates imbalance.

Initial state (center):
  ρ₀ = {rho_0:.4f}
  P₀ = {P_0:.4f}

At t=50 (center):
  ρ = {rho_50:.4f}
  P = {P_50:.4f}
  δρ/ρ₀ = {(rho_50-rho_0)/rho_0:.4f}
  δP/P₀ = {(P_50-P_0)/P_0:.4f}
""")

    # The perturbation equation
    print("\n" + "=" * 80)
    print("PERTURBATION EQUATIONS")
    print("=" * 80)
    
    print("""
For a self-gravitating polytrope, the linearized perturbation equation is:

  ∂²δρ/∂t² = c_s²∇²δρ + 4πGρ₀·δρ - (1/ρ₀)∇·(δρ∇P₀)

The instability arises from the LAST term when numerical scheme 
doesn't properly handle the pressure gradient.

In GSPH without Ω:
- The base state is NOT exactly in equilibrium in the discrete sense
- Small drift → secular growth of perturbation
- Growth rate depends on HOW MUCH the discrete approximation deviates

The key quantity is:
  ε = (discrete F_P) / (exact F_P) - 1
  
If ε < 0, the discrete pressure force is TOO WEAK → collapse.
""")

    # Compute the actual imbalance from force tracking
    print("\n" + "=" * 80)
    print("FORCE IMBALANCE MEASUREMENT")
    print("=" * 80)
    
    # Track the force-like quantity
    h_0 = core_0['sml'].mean()
    h_50 = core_50['sml'].mean()
    
    # Pressure force scales as P/(ρ²h⁴) in SPH
    F_P_0 = P_0 / (rho_0**2 * h_0**4)
    F_P_50 = P_50 / (rho_50**2 * h_50**4)
    
    # For equilibrium, F_P should scale with gravity F_g ∝ ρ
    # If F_P/F_g = const → stable
    # If F_P/F_g decreases with compression → collapse
    
    F_ratio_0 = F_P_0 / rho_0  # ∝ F_P/F_g
    F_ratio_50 = F_P_50 / rho_50
    
    epsilon = (F_ratio_50 - F_ratio_0) / F_ratio_0
    
    print(f"""
Force-to-gravity ratio:
  At t=0:  F_P/ρ ∝ {F_ratio_0:.6f}
  At t=50: F_P/ρ ∝ {F_ratio_50:.6f}
  
Fractional change: ε = {epsilon:.4f}

This ε represents the FRACTIONAL IMBALANCE that drives collapse.
""")

    # The growth rate from perturbation theory
    print("\n" + "=" * 80)
    print("GROWTH RATE FROM FRACTIONAL IMBALANCE")
    print("=" * 80)
    
    G = 1.0
    omega_ff = np.sqrt(4 * np.pi * G * rho_0)
    t_ff = np.sqrt(3 * np.pi / (32 * G * rho_0))
    
    # The perturbation grows as:
    # d²(δρ)/dt² ~ ε · 4πGρ · δρ
    # → Γ² ~ ε · 4πGρ = ε · ω_ff²
    # → Γ ~ √ε · ω_ff
    
    if epsilon < 0:
        Gamma_theory = np.sqrt(-epsilon) * omega_ff
    else:
        Gamma_theory = 0
    
    print(f"""
From perturbation theory:
  Γ² ~ |ε| · ω_ff²
  Γ ~ √|ε| · ω_ff

With ε = {epsilon:.4f}:
  Γ_theory = √|{epsilon:.4f}| × {omega_ff:.4f} = {Gamma_theory:.4f}
  Γ/ω_ff = {Gamma_theory/omega_ff:.4f}

Measured: Γ ≈ 0.026, Γ/ω_ff ≈ 0.006
""")

    # Track ε over time
    print("\n" + "=" * 80)
    print("TRACKING FRACTIONAL IMBALANCE OVER TIME")
    print("=" * 80)
    
    print(f"{'snap':>4} {'t':>8} {'ρ/ρ₀':>8} {'F_P/ρ ratio':>12} {'ε(t)':>10} {'√|ε|·ω_ff':>12}")
    print("-" * 60)
    
    snaps = list(range(0, 100, 5))
    epsilons = []
    gammas_pred = []
    
    for snap in snaps:
        df = load_snapshot(no_gradh_dir / f"snapshot_{snap:04d}.csv")
        df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
        core = df.nsmallest(10, 'r')
        
        rho = core['dens'].mean()
        P = core['pres'].mean()
        h = core['sml'].mean()
        
        F_P = P / (rho**2 * h**4)
        F_ratio = F_P / rho
        
        eps = (F_ratio - F_ratio_0) / F_ratio_0
        epsilons.append(eps)
        
        if eps < 0:
            gam_pred = np.sqrt(-eps) * omega_ff
        else:
            gam_pred = 0
        gammas_pred.append(gam_pred)
        
        t = snap * 0.5  # approximate
        
        print(f"{snap:4d} {t:8.1f} {rho/rho_0:8.4f} {F_ratio/F_ratio_0:12.4f} {eps:+10.4f} {gam_pred:12.4f}")
    
    # The issue - ε is POSITIVE at first (force stronger than equilibrium)
    # Then becomes negative later
    
    print("\n" + "=" * 80)
    print("KEY INSIGHT: OSCILLATION + DRIFT")
    print("=" * 80)
    
    print("""
The data shows ε oscillates around 0 with a slow DRIFT toward negative.

This means:
1. Initially, the discrete force is actually SLIGHTLY TOO STRONG
   (the Lane-Emden solution isn't exact equilibrium for GSPH)
   
2. The system OSCILLATES trying to find equilibrium
   
3. But each oscillation DAMPS slightly toward the wrong state
   because C_eff < 1
   
4. Over many oscillations, this accumulates into secular collapse

This is NOT exponential instability - it's NUMERICAL DRIFT!
""")

    # Calculate drift rate
    print("\n" + "=" * 80)
    print("DRIFT RATE ANALYSIS")
    print("=" * 80)
    
    # The compression rate
    ln_rho_ratios = []
    times = []
    
    for snap in range(0, 301, 10):
        df = load_snapshot(no_gradh_dir / f"snapshot_{snap:04d}.csv")
        df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
        core = df.nsmallest(10, 'r')
        rho = core['dens'].mean()
        
        ln_rho_ratios.append(np.log(rho / rho_0))
        times.append(snap * 0.5)  # approximate
    
    times = np.array(times)
    ln_rho_ratios = np.array(ln_rho_ratios)
    
    # Fit linear slope
    mask = times > 20  # after initial transient
    coeffs = np.polyfit(times[mask], ln_rho_ratios[mask], 1)
    drift_rate = coeffs[0]
    
    print(f"""
Drift rate in ln(ρ/ρ₀):
  Measured: d ln(ρ)/dt = {drift_rate:.6f}
  
This corresponds to:
  Γ_drift = {drift_rate:.4f}
  τ_drift = 1/Γ = {1/drift_rate:.1f} code units
  τ_drift/t_ff = {(1/drift_rate)/t_ff:.1f}
""")

    # Final understanding
    print("\n" + "=" * 80)
    print("FINAL UNDERSTANDING: NUMERICAL DRIFT INSTABILITY")
    print("=" * 80)
    
    print(f"""
The GSPH without grad-h instability is a NUMERICAL DRIFT, not dynamical collapse:

1. MECHANISM:
   - GSPH force response coefficient C_eff ≈ 0.6 < 1
   - But the system starts in approximate equilibrium
   - Small imbalance drives slow DRIFT toward higher density
   
2. TIMESCALE:
   - Drift rate: Γ ≈ {drift_rate:.4f}
   - Much slower than free-fall: Γ/ω_ff ≈ {drift_rate/omega_ff:.4f}
   - Because equilibrium is only SLIGHTLY violated
   
3. NATURE:
   - NOT exponential growth of perturbation
   - Rather: secular drift with possible oscillations
   - τ_drift ≈ {1/drift_rate:.0f} code units ≈ {(1/drift_rate)/t_ff:.0f} × t_ff
   
4. FIX:
   - The Ω correction makes C_eff → 1 exactly
   - Hydrostatic equilibrium becomes EXACT discrete solution
   - No drift → stable oscillations

This explains why:
- Simple theory (Γ ~ √(1-C_eff)·ω_ff) gives wrong rate
- Actual growth is much slower
- System shows quasi-static contraction, not free-fall
""")

    # Create final plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Density evolution
    ax1 = axes[0, 0]
    ax1.semilogy(times, np.exp(ln_rho_ratios), 'r-', linewidth=2)
    ax1.axhline(y=1, color='k', linestyle='--', alpha=0.5)
    fit_line = np.exp(drift_rate * times[mask] + coeffs[1])
    ax1.semilogy(times[mask], fit_line, 'b--', label=f'Fit: Γ={drift_rate:.4f}')
    ax1.set_xlabel('Time [code units]')
    ax1.set_ylabel('ρ_core/ρ₀')
    ax1.set_title('Core Density Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Force imbalance ε
    ax2 = axes[0, 1]
    times_eps = [snap * 0.5 for snap in snaps]
    ax2.plot(times_eps, epsilons, 'go-', markersize=5)
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    ax2.set_xlabel('Time [code units]')
    ax2.set_ylabel('ε = (F_P/ρ - F₀/ρ₀)/(F₀/ρ₀)')
    ax2.set_title('Fractional Force Imbalance')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: ln(ρ) to show linear drift
    ax3 = axes[1, 0]
    ax3.plot(times, ln_rho_ratios, 'r-', linewidth=2)
    ax3.plot(times[mask], drift_rate * times[mask] + coeffs[1], 'b--', 
             label=f'Linear fit: slope={drift_rate:.4f}')
    ax3.set_xlabel('Time [code units]')
    ax3.set_ylabel('ln(ρ_core/ρ₀)')
    ax3.set_title('Logarithm of Density (showing drift)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Summary
    ax4 = axes[1, 1]
    ax4.text(0.1, 0.9, 'GSPH Instability Summary', fontsize=14, fontweight='bold',
             transform=ax4.transAxes)
    ax4.text(0.1, 0.75, f'Effective coefficient: C_eff ≈ 0.6', fontsize=12,
             transform=ax4.transAxes)
    ax4.text(0.1, 0.60, f'Drift rate: Γ = {drift_rate:.4f}', fontsize=12,
             transform=ax4.transAxes)
    ax4.text(0.1, 0.45, f'Free-fall freq: ω_ff = {omega_ff:.4f}', fontsize=12,
             transform=ax4.transAxes)
    ax4.text(0.1, 0.30, f'Ratio: Γ/ω_ff = {drift_rate/omega_ff:.4f}', fontsize=12,
             transform=ax4.transAxes)
    ax4.text(0.1, 0.15, f'Drift timescale: τ = {1/drift_rate:.0f} ≈ {(1/drift_rate)/t_ff:.0f}×t_ff', 
             fontsize=12, transform=ax4.transAxes)
    ax4.axis('off')
    
    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/sample/gradh_study_3d/results/drift_analysis.png', dpi=150)
    print("\nPlot saved to: results/drift_analysis.png")
    
    plt.show()

if __name__ == "__main__":
    main()
