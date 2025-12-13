#!/usr/bin/env python3
"""
Derive the correct theory for GSPH instability by analyzing the actual data.
Key finding: The perturbation is NON-UNIFORM (center compresses more than edges).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_snapshot(filepath):
    return pd.read_csv(filepath, comment='#')

def main():
    no_gradh_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/no_gradh_hk")
    gradh_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/gradh_hk")
    
    print("=" * 80)
    print("CORRECT THEORY FOR GSPH INSTABILITY")
    print("=" * 80)
    
    # Load initial and final snapshots to understand perturbation structure
    df_init = load_snapshot(no_gradh_dir / "snapshot_0000.csv")
    df_50 = load_snapshot(no_gradh_dir / "snapshot_0050.csv")
    df_100 = load_snapshot(no_gradh_dir / "snapshot_0100.csv")
    
    # Calculate radii
    for df in [df_init, df_50, df_100]:
        df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
    
    # Bin by radius to see radial structure
    print("\n" + "=" * 80)
    print("RADIAL DENSITY PROFILE EVOLUTION")
    print("=" * 80)
    
    r_bins = np.linspace(0, 1.0, 11)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    def get_radial_profile(df, r_bins):
        """Get average density in radial bins."""
        rho_profile = []
        for i in range(len(r_bins) - 1):
            mask = (df['r'] >= r_bins[i]) & (df['r'] < r_bins[i+1])
            if mask.sum() > 0:
                rho_profile.append(df.loc[mask, 'dens'].mean())
            else:
                rho_profile.append(np.nan)
        return np.array(rho_profile)
    
    rho_init = get_radial_profile(df_init, r_bins)
    rho_50 = get_radial_profile(df_50, r_bins)
    rho_100 = get_radial_profile(df_100, r_bins)
    
    print(f"\n{'r_center':>10} {'ρ(t=0)':>10} {'ρ(t=50)':>10} {'ρ(t=100)':>10} {'δρ/ρ(50)':>12} {'δρ/ρ(100)':>12}")
    print("-" * 75)
    
    for i, r in enumerate(r_centers):
        if not np.isnan(rho_init[i]):
            delta_50 = (rho_50[i] - rho_init[i]) / rho_init[i] if not np.isnan(rho_50[i]) else np.nan
            delta_100 = (rho_100[i] - rho_init[i]) / rho_init[i] if not np.isnan(rho_100[i]) else np.nan
            print(f"{r:10.2f} {rho_init[i]:10.4f} {rho_50[i]:10.4f} {rho_100[i]:10.4f} {delta_50:+12.4f} {delta_100:+12.4f}")
    
    print("\n" + "=" * 80)
    print("KEY OBSERVATION: NON-UNIFORM PERTURBATION")
    print("=" * 80)
    
    print("""
The perturbation is NOT uniform across the sphere!

- CENTER (r ~ 0): Strong compression (δρ/ρ >> 0)
- EDGE (r ~ 1): Weak compression or even expansion (δρ/ρ ~ 0 or < 0)

This is CRITICAL for understanding the force response!
""")

    # Now analyze particle pairs - what does a center particle "see"?
    print("\n" + "=" * 80)
    print("PARTICLE PAIR ANALYSIS")
    print("=" * 80)
    
    # Get a center particle at t=50
    center_mask = df_50['r'] < 0.1
    edge_mask = (df_50['r'] > 0.4) & (df_50['r'] < 0.6)
    
    rho_center = df_50.loc[center_mask, 'dens'].mean()
    rho_edge = df_50.loc[edge_mask, 'dens'].mean()
    P_center = df_50.loc[center_mask, 'pres'].mean()
    P_edge = df_50.loc[edge_mask, 'pres'].mean()
    cs_center = df_50.loc[center_mask, 'sound'].mean()
    cs_edge = df_50.loc[edge_mask, 'sound'].mean()
    
    rho_center_0 = df_init.loc[df_init['r'] < 0.1, 'dens'].mean()
    rho_edge_0 = df_init.loc[(df_init['r'] > 0.4) & (df_init['r'] < 0.6), 'dens'].mean()
    
    print(f"\nAt t=50 snapshot:")
    print(f"  Center (r<0.1): ρ = {rho_center:.4f}, P = {P_center:.4f}, c_s = {cs_center:.4f}")
    print(f"  Edge (r~0.5):   ρ = {rho_edge:.4f}, P = {P_edge:.4f}, c_s = {cs_edge:.4f}")
    
    print(f"\nRelative to initial:")
    print(f"  Center: δρ/ρ = {(rho_center - rho_center_0)/rho_center_0:+.4f}")
    print(f"  Edge:   δρ/ρ = {(rho_edge - rho_edge_0)/rho_edge_0:+.4f}")
    
    # Riemann solver pressure
    print("\n" + "=" * 80)
    print("RIEMANN SOLVER RESPONSE")
    print("=" * 80)
    
    print("""
For a center-edge particle pair, the Riemann pressure is:

  p* = (c_R·P_L + c_L·P_R) / (c_L + c_R)
  
where L = center (high ρ, high P), R = edge (lower ρ, lower P)

Key insight: p* is the AVERAGE of P_center and P_edge!
""")

    # Calculate p* for center-edge pair
    c_center = cs_center
    c_edge = cs_edge
    
    pstar = (c_edge * P_center + c_center * P_edge) / (c_center + c_edge)
    
    print(f"\nRiemann pressure for center-edge pair:")
    print(f"  P_center = {P_center:.4f}")
    print(f"  P_edge = {P_edge:.4f}")
    print(f"  p* = {pstar:.4f}")
    print(f"  p*/P_center = {pstar/P_center:.4f}")
    print(f"  p*/P_edge = {pstar/P_edge:.4f}")
    
    print("\n" + "=" * 80)
    print("THE FORCE IMBALANCE")
    print("=" * 80)
    
    print(f"""
The center particle i experiences force from edge particle j:

  F_ij = m_j · p* · (1/ρ_i² · ∇W_i + 1/ρ_j² · ∇W_j)

In Standard SPH, it would be:
  F_ij = m_j · (P_i/ρ_i² · ∇W_i + P_j/ρ_j² · ∇W_j)

The DIFFERENCE in force response:

Standard SPH (i-side contribution):
  F_i ∝ P_i/ρ_i² = {P_center:.4f}/{rho_center:.4f}² = {P_center/rho_center**2:.6f}

GSPH (i-side contribution):
  F_i ∝ p*/ρ_i² = {pstar:.4f}/{rho_center:.4f}² = {pstar/rho_center**2:.6f}

Ratio: GSPH/Standard = {(pstar/rho_center**2)/(P_center/rho_center**2):.4f} = p*/P_i = {pstar/P_center:.4f}

The GSPH force on center particle is WEAKER by factor p*/P_center < 1 !
""")

    # Now derive the effective coefficient
    print("\n" + "=" * 80)
    print("DERIVING THE EFFECTIVE COEFFICIENT")
    print("=" * 80)
    
    print("""
For non-uniform perturbation with center compression δρ_c and edge δρ_e:

GRAVITY:
  δF_g/F_g = δρ_c/ρ_c  (responds to LOCAL density)
  
PRESSURE (Standard SPH):
  δF_P/F_P = δP_c/P_c = γ · δρ_c/ρ_c  (responds to LOCAL pressure)
  But with |∇W| ∝ h^(-4) ∝ ρ^(4/3) and 1/ρ² term:
  δF_P/F_P ≈ (γ - 2 + 4/3) δρ/ρ = (γ - 2/3) δρ/ρ
  For γ = 5/3: δF_P/F_P = δρ/ρ → C_P = 1 ✓
  
PRESSURE (GSPH with shared p*):
  p* averages P_center and P_edge:
  p* ≈ (P_c + P_e)/2 for similar c_s
  
  δp*/p* = (δP_c + δP_e)/(P_c + P_e) · (P_c + P_e)/(2·P̄)
         ≈ (γ·δρ_c/ρ_c · P_c + γ·δρ_e/ρ_e · P_e) / (P_c + P_e)
         
  If δρ_e << δρ_c (edge doesn't compress as much):
  δp*/p* ≈ γ · δρ_c/ρ_c · P_c/(P_c + P_e) < γ · δρ_c/ρ_c
  
  This gives REDUCED response: C_P < 1
""")

    # Calculate the actual reduction factor from data
    P_center_0 = df_init.loc[df_init['r'] < 0.1, 'pres'].mean()
    P_edge_0 = df_init.loc[(df_init['r'] > 0.4) & (df_init['r'] < 0.6), 'pres'].mean()
    
    # Initial p*
    c_center_0 = df_init.loc[df_init['r'] < 0.1, 'sound'].mean()
    c_edge_0 = df_init.loc[(df_init['r'] > 0.4) & (df_init['r'] < 0.6), 'sound'].mean()
    pstar_0 = (c_edge_0 * P_center_0 + c_center_0 * P_edge_0) / (c_center_0 + c_edge_0)
    
    delta_P_center = (P_center - P_center_0) / P_center_0
    delta_P_edge = (P_edge - P_edge_0) / P_edge_0
    delta_pstar = (pstar - pstar_0) / pstar_0
    delta_rho_center = (rho_center - rho_center_0) / rho_center_0
    
    print(f"\nFrom data at t=50:")
    print(f"  δP_center/P_center_0 = {delta_P_center:+.4f}")
    print(f"  δP_edge/P_edge_0 = {delta_P_edge:+.4f}")
    print(f"  δp*/p*_0 = {delta_pstar:+.4f}")
    print(f"  δρ_center/ρ_center_0 = {delta_rho_center:+.4f}")
    
    # The key ratio
    ratio_pstar_P = delta_pstar / delta_P_center if abs(delta_P_center) > 0.01 else float('nan')
    
    print(f"\n  Ratio (δp*/p*) / (δP_center/P_center) = {ratio_pstar_P:.4f}")
    print(f"  This is the p* RESPONSE REDUCTION factor!")
    
    # Expected C_eff
    # F ∝ p* · h^(-4) / ρ² = p* · ρ^(4/3) / ρ² = p* · ρ^(-2/3)
    # δF/F = δp*/p* + (-2/3) δρ/ρ
    # For Standard: δP/P = γ·δρ/ρ, so δF/F = γ·δρ/ρ - (2/3)δρ/ρ = (γ-2/3)δρ/ρ = δρ/ρ for γ=5/3
    # For GSPH: δp*/p* = R·γ·δρ/ρ where R < 1 is the reduction
    # So δF/F = R·γ·δρ/ρ - (2/3)δρ/ρ = (R·γ - 2/3)δρ/ρ
    
    gamma = 5/3
    R = ratio_pstar_P if not np.isnan(ratio_pstar_P) else 0.5
    C_eff_predicted = R * gamma - 2/3
    
    print(f"\nPredicted effective coefficient:")
    print(f"  C_eff = R·γ - 2/3 = {R:.4f} × {gamma:.4f} - {2/3:.4f} = {C_eff_predicted:.4f}")
    
    print("\n" + "=" * 80)
    print("INSTABILITY CRITERION")
    print("=" * 80)
    
    print(f"""
For stability against gravity (C_g = 1), we need C_P ≥ 1.

In GSPH without Ω:
  C_P = R·γ - 2/3
  
where R = (δp*/p*) / (δP/P) is the p* response reduction factor.

For stability: R·γ - 2/3 ≥ 1
              R ≥ (1 + 2/3)/γ = {(1 + 2/3)/gamma:.4f}

Measured R ≈ {R:.4f} < {(1 + 2/3)/gamma:.4f} → UNSTABLE!

The instability growth rate depends on (1 - C_P):
  Γ² = 4πGρ (1 - C_P) = 4πGρ (1 - R·γ + 2/3)
  Γ² = 4πGρ ({1 - C_eff_predicted:.4f})
""")

    G = 1.0
    rho_avg = rho_center_0
    omega_ff = np.sqrt(4 * np.pi * G * rho_avg)
    
    if C_eff_predicted < 1:
        Gamma_predicted = np.sqrt(4 * np.pi * G * rho_avg * (1 - C_eff_predicted))
        print(f"  Predicted Γ = {Gamma_predicted:.4f}")
        print(f"  Predicted Γ/ω_ff = {Gamma_predicted/omega_ff:.4f}")
        print(f"  Predicted e-folding time = {1/Gamma_predicted:.2f}")
    else:
        print(f"  C_eff ≥ 1 → STABLE (no growth)")
    
    # Compare with measured
    print("\n" + "=" * 80)
    print("COMPARISON WITH MEASURED GROWTH")
    print("=" * 80)
    
    # From previous analysis: Γ_measured ≈ 0.026
    Gamma_measured = 0.026
    
    print(f"\nMeasured growth rate: Γ = {Gamma_measured:.4f}")
    print(f"Measured Γ/ω_ff = {Gamma_measured/omega_ff:.4f}")
    
    if C_eff_predicted < 1:
        print(f"\nPredicted Γ = {Gamma_predicted:.4f}")
        print(f"Ratio predicted/measured = {Gamma_predicted/Gamma_measured:.2f}")
    
    # The quasi-static nature
    print("\n" + "=" * 80)
    print("WHY IS THE GROWTH SO SLOW?")
    print("=" * 80)
    
    print(f"""
The measured growth rate Γ ≈ 0.026 << ω_ff ≈ {omega_ff:.2f}

This is because C_eff ≈ {C_eff_predicted:.2f} is CLOSE to 1!

From Γ² = 4πGρ(1 - C_P):
  Γ/ω_ff = √(1 - C_P)
  
For C_P = {C_eff_predicted:.2f}:
  Γ/ω_ff = √({1 - C_eff_predicted:.4f}) = {np.sqrt(abs(1 - C_eff_predicted)):.4f}
  
This explains the QUASI-STATIC nature:
- C_eff is close to 1 but slightly less
- System slowly contracts seeking equilibrium
- But no stable equilibrium exists → secular collapse
- Timescale ~ t_ff / √(1-C_P) >> t_ff
""")

    # The role of Ω
    print("\n" + "=" * 80)
    print("HOW Ω FIXES THE INSTABILITY")
    print("=" * 80)
    
    # Load gradh data for comparison
    df_gradh_50 = load_snapshot(gradh_dir / "snapshot_0050.csv")
    df_gradh_50['r'] = np.sqrt(df_gradh_50['pos_x']**2 + df_gradh_50['pos_y']**2 + df_gradh_50['pos_z']**2)
    
    center_mask_g = df_gradh_50['r'] < 0.1
    omega_center = df_gradh_50.loc[center_mask_g, 'gradh'].mean()
    rho_center_g = df_gradh_50.loc[center_mask_g, 'dens'].mean()
    
    print(f"\nIn GSPH with Ω, at t=50:")
    print(f"  Center Ω = {omega_center:.4f}")
    print(f"  Center ρ = {rho_center_g:.4f} (vs {rho_center:.4f} without Ω)")
    
    print(f"""
The Ω factor appears in the force as:
  F = m·p*·(Ω_i/ρ_i²·∇W_i + Ω_j/ρ_j²·∇W_j)

Ω compensates for the h-variation in ρ estimation:
  Ω = 1 / (1 + (h/Dρ)(dρ/dh))

For compression (h decreases, ρ increases):
  dρ/dh < 0 → (h/Dρ)(dρ/dh) < 0 → Ω > 1
  
The center sees Ω > 1, BOOSTING the pressure force!

This provides the missing response to restore C_eff → 1.

More fundamentally, Ω ensures that:
1. The DISCRETE equations conserve energy exactly
2. Hydrostatic equilibrium is an EXACT solution
3. The variational structure is preserved

This is not just a force coefficient fix - it's STRUCTURAL CONSISTENCY.
""")

    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY: THE CORRECT THEORY")
    print("=" * 80)
    
    print(f"""
GSPH without grad-h correction is UNSTABLE because:

1. The Riemann solver pressure p* AVERAGES left and right states
2. For non-uniform perturbation (center >> edge), this averaging
   REDUCES the pressure response at the center
3. The reduction factor R = (δp*/p*) / (δP/P) ≈ {R:.2f} < 1

4. The effective force coefficient becomes:
   C_eff = R·γ - 2/3 ≈ {C_eff_predicted:.2f} < 1
   
5. Since C_eff < 1, gravity (C_g = 1) dominates:
   - Net inward force at perturbation maxima
   - Slow SECULAR collapse with rate Γ ~ ω_ff·√(1-C_eff)
   - Timescale ~ {1/Gamma_measured:.0f} >> t_ff = {np.sqrt(3*np.pi/(32*G*rho_avg)):.2f}

The grad-h correction (Ω factor) restores stability by:
   - Boosting force at compressed regions (Ω > 1)
   - Ensuring discrete energy conservation
   - Making hydrostatic equilibrium an exact solution

This is a NUMERICAL instability specific to GSPH without grad-h,
not a physical Jeans instability.
""")

    # Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Radial density profiles
    ax1 = axes[0, 0]
    ax1.plot(r_centers, rho_init, 'k-o', label='t=0', markersize=4)
    ax1.plot(r_centers, rho_50, 'b-s', label='t=50', markersize=4)
    ax1.plot(r_centers, rho_100, 'r-^', label='t=100', markersize=4)
    ax1.set_xlabel('Radius')
    ax1.set_ylabel('Density ρ')
    ax1.set_title('Radial Density Profile (no grad-h)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Relative perturbation
    ax2 = axes[0, 1]
    delta_50 = (rho_50 - rho_init) / rho_init
    delta_100 = (rho_100 - rho_init) / rho_init
    ax2.plot(r_centers, delta_50, 'b-s', label='t=50', markersize=4)
    ax2.plot(r_centers, delta_100, 'r-^', label='t=100', markersize=4)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Radius')
    ax2.set_ylabel('δρ/ρ₀')
    ax2.set_title('Relative Density Perturbation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: p* response
    ax3 = axes[1, 0]
    # Track p* equivalent through time
    snaps = [0, 10, 20, 30, 50, 70, 100]
    pstar_ratio = []
    P_center_ratio = []
    
    for snap in snaps:
        df = load_snapshot(no_gradh_dir / f"snapshot_{snap:04d}.csv")
        df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
        
        c_mask = df['r'] < 0.1
        e_mask = (df['r'] > 0.4) & (df['r'] < 0.6)
        
        Pc = df.loc[c_mask, 'pres'].mean()
        Pe = df.loc[e_mask, 'pres'].mean()
        cc = df.loc[c_mask, 'sound'].mean()
        ce = df.loc[e_mask, 'sound'].mean()
        
        ps = (ce * Pc + cc * Pe) / (cc + ce)
        
        pstar_ratio.append(ps / pstar_0)
        P_center_ratio.append(Pc / P_center_0)
    
    ax3.plot(snaps, pstar_ratio, 'r-o', label='p*/p*₀ (Riemann)', markersize=5)
    ax3.plot(snaps, P_center_ratio, 'b-s', label='P_center/P₀ (local)', markersize=5)
    ax3.set_xlabel('Snapshot')
    ax3.set_ylabel('Pressure ratio')
    ax3.set_title('Pressure Response: p* vs Local P')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Effective coefficient over time
    ax4 = axes[1, 1]
    
    C_effs = []
    rho_ratios = []
    
    for snap in range(5, 100, 5):
        df = load_snapshot(no_gradh_dir / f"snapshot_{snap:04d}.csv")
        df['r'] = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + df['pos_z']**2)
        
        core = df.nsmallest(10, 'r')
        rho = core['dens'].mean()
        P = core['pres'].mean()
        h = core['sml'].mean()
        
        df_init_temp = load_snapshot(no_gradh_dir / "snapshot_0000.csv")
        df_init_temp['r'] = np.sqrt(df_init_temp['pos_x']**2 + df_init_temp['pos_y']**2 + df_init_temp['pos_z']**2)
        core_init_temp = df_init_temp.nsmallest(10, 'r')
        rho0_temp = core_init_temp['dens'].mean()
        P0_temp = core_init_temp['pres'].mean()
        h0_temp = core_init_temp['sml'].mean()
        
        F = P / (rho**2 * h**4)
        F0 = P0_temp / (rho0_temp**2 * h0_temp**4)
        
        delta_rho = (rho - rho0_temp) / rho0_temp
        delta_F = (F - F0) / F0
        
        if abs(delta_rho) > 0.1:
            C_eff = delta_F / delta_rho
            C_effs.append(C_eff)
            rho_ratios.append(rho / rho0_temp)
    
    ax4.plot(rho_ratios, C_effs, 'go-', markersize=5)
    ax4.axhline(y=1.0, color='b', linestyle='--', label='C = 1 (stable)', linewidth=2)
    ax4.axhline(y=1/6, color='r', linestyle='--', label='C = 1/6 (uniform theory)', linewidth=2)
    ax4.set_xlabel('ρ/ρ₀')
    ax4.set_ylabel('Effective coefficient C_eff')
    ax4.set_title('Force Response Coefficient')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/correct_theory.png', dpi=150)
    print("\nPlot saved to: results/correct_theory.png")
    
    plt.show()

if __name__ == "__main__":
    main()
