#!/usr/bin/env python3
"""
Deep analysis of grad-h correction effect
Verify theoretical predictions vs actual simulation data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_snapshot(filepath):
    """Load snapshot CSV, skipping comment lines"""
    return pd.read_csv(filepath, comment='#')

def analyze_all_particles(data):
    """Get all particle data with radius"""
    r = np.sqrt(data['pos_x']**2 + data['pos_y']**2 + data['pos_z']**2)
    data = data.copy()
    data['r'] = r
    return data

def main():
    # Paths
    gradh_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/gradh_hk")
    no_gradh_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/no_gradh_hk")
    
    print("=" * 80)
    print("DETAILED GRAD-H ANALYSIS")
    print("=" * 80)
    
    # Load initial state
    gradh_0 = load_snapshot(gradh_dir / "snapshot_0000.csv")
    gradh_0 = analyze_all_particles(gradh_0)
    
    no_gradh_0 = load_snapshot(no_gradh_dir / "snapshot_0000.csv")
    no_gradh_0 = analyze_all_particles(no_gradh_0)
    
    # Focus on innermost particle
    idx_center_gradh = gradh_0['r'].idxmin()
    idx_center_no_gradh = no_gradh_0['r'].idxmin()
    
    print("\n" + "=" * 80)
    print("INNERMOST PARTICLE AT t=0")
    print("=" * 80)
    
    p_gradh = gradh_0.loc[idx_center_gradh]
    p_no_gradh = no_gradh_0.loc[idx_center_no_gradh]
    
    print(f"\nGrad-h case:")
    print(f"  id = {p_gradh['id']}")
    print(f"  r = {p_gradh['r']:.6f}")
    print(f"  ρ = {p_gradh['dens']:.6f}")
    print(f"  P = {p_gradh['pres']:.6f}")
    print(f"  h = {p_gradh['sml']:.6f}")
    print(f"  Ω = {p_gradh['gradh']:.6f}")
    print(f"  neighbors = {p_gradh['neighbor']}")
    
    print(f"\nNo grad-h case:")
    print(f"  id = {p_no_gradh['id']}")
    print(f"  r = {p_no_gradh['r']:.6f}")
    print(f"  ρ = {p_no_gradh['dens']:.6f}")
    print(f"  P = {p_no_gradh['pres']:.6f}")
    print(f"  h = {p_no_gradh['sml']:.6f}")
    print(f"  Ω = {p_no_gradh['gradh']:.6f}")
    
    # Track the same particle through time
    particle_id = int(p_gradh['id'])
    print(f"\n" + "=" * 80)
    print(f"TRACKING PARTICLE ID={particle_id} THROUGH TIME")
    print("=" * 80)
    
    snapshots = [0, 10, 20, 50, 100, 150, 200, 250, 300]
    
    data_gradh = []
    data_no_gradh = []
    
    for snap in snapshots:
        gradh_file = gradh_dir / f"snapshot_{snap:04d}.csv"
        no_gradh_file = no_gradh_dir / f"snapshot_{snap:04d}.csv"
        
        if gradh_file.exists():
            df = load_snapshot(gradh_file)
            df = analyze_all_particles(df)
            p = df[df['id'] == particle_id]
            if len(p) > 0:
                p = p.iloc[0]
                data_gradh.append({
                    'snap': snap,
                    'r': p['r'],
                    'rho': p['dens'],
                    'P': p['pres'],
                    'h': p['sml'],
                    'omega': p['gradh'],
                    'neighbors': p['neighbor']
                })
        
        if no_gradh_file.exists():
            df = load_snapshot(no_gradh_file)
            df = analyze_all_particles(df)
            p = df[df['id'] == particle_id]
            if len(p) > 0:
                p = p.iloc[0]
                data_no_gradh.append({
                    'snap': snap,
                    'r': p['r'],
                    'rho': p['dens'],
                    'P': p['pres'],
                    'h': p['sml'],
                    'omega': p['gradh'],
                    'neighbors': p['neighbor']
                })
    
    df_gradh = pd.DataFrame(data_gradh)
    df_no_gradh = pd.DataFrame(data_no_gradh)
    
    print("\nGrad-h case (same particle):")
    print(df_gradh.to_string(index=False))
    
    print("\nNo grad-h case (same particle):")
    print(df_no_gradh.to_string(index=False))
    
    # Calculate theoretical predictions
    print("\n" + "=" * 80)
    print("THEORETICAL CALCULATIONS")
    print("=" * 80)
    
    gamma = 5/3
    
    # Initial values
    rho_0 = df_gradh['rho'].iloc[0]
    h_0 = df_gradh['h'].iloc[0]
    omega_0 = df_gradh['omega'].iloc[0]
    P_0 = df_gradh['P'].iloc[0]
    
    print(f"\nInitial state:")
    print(f"  ρ₀ = {rho_0:.6f}")
    print(f"  h₀ = {h_0:.6f}")
    print(f"  Ω₀ = {omega_0:.6f}")
    print(f"  P₀ = {P_0:.6f}")
    
    print("\n" + "-" * 60)
    print("Analysis for each snapshot:")
    print("-" * 60)
    
    for i in range(len(df_gradh)):
        row = df_gradh.iloc[i]
        snap = row['snap']
        rho = row['rho']
        h = row['h']
        omega = row['omega']
        P = row['P']
        
        # Relative changes
        delta_rho_rel = (rho - rho_0) / rho_0
        delta_h_rel = (h - h_0) / h_0
        delta_omega_rel = (omega - omega_0) / omega_0
        delta_P_rel = (P - P_0) / P_0
        
        # Theoretical: h ∝ ρ^(-1/3)
        h_theory = h_0 * (rho / rho_0)**(-1/3)
        
        # Theoretical: |∇W| ∝ h^(-4)
        grad_W_rel = (h_0 / h)**4 - 1  # relative change
        
        # Theoretical: P = K ρ^γ (isentropic)
        P_theory = P_0 * (rho / rho_0)**gamma
        
        # For p* (Riemann solver): δp*/p* ≈ (γ/2) δρ/ρ
        pstar_rel_theory = (gamma / 2) * delta_rho_rel
        
        print(f"\nSnapshot {snap}:")
        print(f"  ρ/ρ₀ = {rho/rho_0:.4f}, δρ/ρ₀ = {delta_rho_rel:+.4f}")
        print(f"  h/h₀ = {h/h_0:.4f}, δh/h₀ = {delta_h_rel:+.4f}")
        print(f"    Theory h/h₀ = (ρ/ρ₀)^(-1/3) = {(rho/rho_0)**(-1/3):.4f}")
        print(f"  Ω/Ω₀ = {omega/omega_0:.4f}, δΩ/Ω₀ = {delta_omega_rel:+.4f}")
        print(f"  P/P₀ = {P/P_0:.4f}, δP/P₀ = {delta_P_rel:+.4f}")
        print(f"    Theory P/P₀ = (ρ/ρ₀)^γ = {(rho/rho_0)**gamma:.4f}")
        
        # Force term analysis
        # F ∝ p* · Ω · |∇W| / ρ²
        # Relative to initial:
        # F/F₀ = (p*/p*₀) · (Ω/Ω₀) · (|∇W|/|∇W|₀) · (ρ₀/ρ)²
        
        grad_W_ratio = (h_0 / h)**4
        force_ratio = (P / P_0) * (omega / omega_0) * grad_W_ratio * (rho_0 / rho)**2
        
        # Using p* scaling instead of P
        # If p* = (c_L P_L + c_R P_R)/(c_L + c_R), for symmetric case p* ≈ P
        # But for perturbation δp*/p* ≈ (γ/2) δρ/ρ
        
        print(f"\n  Force analysis (F ∝ p* · Ω · |∇W| / ρ²):")
        print(f"    |∇W|/|∇W|₀ = (h₀/h)⁴ = {grad_W_ratio:.4f}")
        print(f"    (ρ₀/ρ)² = {(rho_0/rho)**2:.4f}")
        print(f"    Ω · |∇W| / |∇W|₀ = {omega/omega_0 * grad_W_ratio:.4f}")
        print(f"    Full F/F₀ (using P) = {force_ratio:.4f}")
    
    # Compare collapse in no_gradh case
    print("\n" + "=" * 80)
    print("COLLAPSE ANALYSIS (no grad-h)")
    print("=" * 80)
    
    rho_0_no = df_no_gradh['rho'].iloc[0]
    h_0_no = df_no_gradh['h'].iloc[0]
    
    for i in range(len(df_no_gradh)):
        row = df_no_gradh.iloc[i]
        snap = row['snap']
        rho = row['rho']
        h = row['h']
        
        delta_rho_rel = (rho - rho_0_no) / rho_0_no
        
        print(f"\nSnapshot {snap}:")
        print(f"  ρ/ρ₀ = {rho/rho_0_no:.2f}, δρ/ρ₀ = {delta_rho_rel:+.2f}")
        print(f"  h/h₀ = {h/h_0_no:.4f}")
        
        # Force ratio without Ω correction
        grad_W_ratio = (h_0_no / h)**4
        force_ratio_no = (rho / rho_0_no)**gamma * grad_W_ratio * (rho_0_no / rho)**2
        # = ρ^γ · h^(-4) · ρ^(-2) 
        # With h ∝ ρ^(-1/3): = ρ^γ · ρ^(4/3) · ρ^(-2) = ρ^(γ + 4/3 - 2) = ρ^(γ - 2/3)
        # For γ = 5/3: = ρ^1
        
        theory_force = (rho / rho_0_no)**(gamma - 2/3)
        print(f"  F/F₀ (computed) = {force_ratio_no:.4f}")
        print(f"  F/F₀ (theory ρ^(γ-2/3)=ρ^1) = {theory_force:.4f}")
    
    # KEY INSIGHT
    print("\n" + "=" * 80)
    print("KEY INSIGHT: WHY Ω STABILIZES")
    print("=" * 80)
    
    print("""
The simple analysis says:
  Without Ω: F ∝ ρ^(γ-2/3) = ρ^1 → coefficient = 1 (marginally stable)
  With Ω: If Ω·|∇W| = const, then F ∝ ρ^(γ/2-2) = ρ^(-7/6) (unstable!)

But the data shows:
  - Without Ω: COLLAPSE occurs (coefficient < 1 in practice)
  - With Ω: STABLE

The resolution: 

1. The "coefficient = 1" for Standard SPH is MARGINAL stability
   - Any numerical error tips it toward instability
   
2. For GSPH without Ω, the coefficient is 1/6 (from p* response)
   - This is clearly < 1 → unstable
   
3. For GSPH with Ω, the LAGRANGIAN consistency matters:
   - Ω ensures the DISCRETE equations conserve energy exactly
   - Ω ensures hydrostatic equilibrium is a DISCRETE steady state
   - This is NOT captured by the simple force coefficient analysis

The Ω correction doesn't just change the force coefficient - it makes the 
NUMERICAL SCHEME variationally consistent, which provides stability through
ENERGY CONSERVATION rather than through force balance response.
""")

    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Density ratio
    ax1 = axes[0, 0]
    ax1.semilogy(df_gradh['snap'], df_gradh['rho']/rho_0, 'b-o', label='with Ω', markersize=4)
    ax1.semilogy(df_no_gradh['snap'], df_no_gradh['rho']/rho_0_no, 'r-s', label='without Ω', markersize=4)
    ax1.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Snapshot')
    ax1.set_ylabel('ρ/ρ₀')
    ax1.set_title('Density Evolution (same particle)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Omega evolution
    ax2 = axes[0, 1]
    ax2.plot(df_gradh['snap'], df_gradh['omega'], 'b-o', markersize=4)
    ax2.axhline(y=1.0, color='r', linestyle='--', label='no correction')
    ax2.set_xlabel('Snapshot')
    ax2.set_ylabel('Ω')
    ax2.set_title('Ω Factor Evolution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: h evolution
    ax3 = axes[1, 0]
    ax3.plot(df_gradh['snap'], df_gradh['h']/h_0, 'b-o', label='with Ω', markersize=4)
    ax3.plot(df_no_gradh['snap'], df_no_gradh['h']/h_0_no, 'r-s', label='without Ω', markersize=4)
    ax3.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Snapshot')
    ax3.set_ylabel('h/h₀')
    ax3.set_title('Smoothing Length Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Ω · |∇W| / |∇W|₀
    ax4 = axes[1, 1]
    omega_gradW_gradh = df_gradh['omega'] * (h_0 / df_gradh['h'])**4
    omega_gradW_no = 1.0 * (h_0_no / df_no_gradh['h'])**4
    
    ax4.plot(df_gradh['snap'], omega_gradW_gradh / omega_gradW_gradh.iloc[0], 'b-o', label='Ω·|∇W| (with Ω)', markersize=4)
    ax4.plot(df_no_gradh['snap'], omega_gradW_no / omega_gradW_no.iloc[0], 'r-s', label='|∇W| (without Ω)', markersize=4)
    ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Snapshot')
    ax4.set_ylabel('Normalized Ω·|∇W|')
    ax4.set_title('Kernel Gradient Scaling')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/gradh_detailed_analysis.png', dpi=150)
    print("\nPlot saved to: results/gradh_detailed_analysis.png")
    
    plt.show()

if __name__ == "__main__":
    main()
