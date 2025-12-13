#!/usr/bin/env python3
"""
Analyze the effect of grad-h correction on center particles
Compare gradh vs no_gradh simulations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_snapshot(filepath):
    """Load snapshot CSV, skipping comment lines"""
    return pd.read_csv(filepath, comment='#')

def analyze_center_particles(data, r_max=0.05):
    """Select particles near center (r < r_max)"""
    r = np.sqrt(data['pos_x']**2 + data['pos_y']**2 + data['pos_z']**2)
    center_mask = r < r_max
    return data[center_mask].copy(), r[center_mask]

def main():
    # Paths
    gradh_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/gradh_hk")
    no_gradh_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/kernel_test/no_gradh_hk")
    
    # Load initial and some later snapshots
    snapshots = [0, 10, 50, 100, 200, 300]
    
    print("=" * 80)
    print("ANALYSIS OF GRAD-H EFFECT ON CENTER PARTICLES")
    print("=" * 80)
    
    # First, analyze initial state to understand Ω values
    print("\n" + "=" * 80)
    print("INITIAL STATE (t=0) - CENTER PARTICLES (r < 0.05)")
    print("=" * 80)
    
    gradh_0 = load_snapshot(gradh_dir / "snapshot_0000.csv")
    no_gradh_0 = load_snapshot(no_gradh_dir / "snapshot_0000.csv")
    
    center_gradh, r_gradh = analyze_center_particles(gradh_0)
    center_no_gradh, r_no_gradh = analyze_center_particles(no_gradh_0)
    
    print(f"\nNumber of center particles (gradh): {len(center_gradh)}")
    print(f"Number of center particles (no_gradh): {len(center_no_gradh)}")
    
    print("\n--- GRAD-H CASE (Ω enabled) ---")
    print(f"  Ω (gradh) range: {center_gradh['gradh'].min():.6f} - {center_gradh['gradh'].max():.6f}")
    print(f"  Ω (gradh) mean:  {center_gradh['gradh'].mean():.6f}")
    print(f"  ρ range: {center_gradh['dens'].min():.6f} - {center_gradh['dens'].max():.6f}")
    print(f"  ρ mean:  {center_gradh['dens'].mean():.6f}")
    print(f"  h range: {center_gradh['sml'].min():.6f} - {center_gradh['sml'].max():.6f}")
    print(f"  h mean:  {center_gradh['sml'].mean():.6f}")
    print(f"  P range: {center_gradh['pres'].min():.6f} - {center_gradh['pres'].max():.6f}")
    
    print("\n--- NO GRAD-H CASE (Ω = 1) ---")
    print(f"  Ω (gradh) range: {center_no_gradh['gradh'].min():.6f} - {center_no_gradh['gradh'].max():.6f}")
    print(f"  ρ range: {center_no_gradh['dens'].min():.6f} - {center_no_gradh['dens'].max():.6f}")
    print(f"  ρ mean:  {center_no_gradh['dens'].mean():.6f}")
    print(f"  h range: {center_no_gradh['sml'].min():.6f} - {center_no_gradh['sml'].max():.6f}")
    print(f"  P range: {center_no_gradh['pres'].min():.6f} - {center_no_gradh['pres'].max():.6f}")
    
    # Track evolution of center density
    print("\n" + "=" * 80)
    print("TIME EVOLUTION OF CENTER DENSITY")
    print("=" * 80)
    
    times_gradh = []
    times_no_gradh = []
    rho_center_gradh = []
    rho_center_no_gradh = []
    omega_center_gradh = []
    h_center_gradh = []
    h_center_no_gradh = []
    
    for snap in snapshots:
        gradh_file = gradh_dir / f"snapshot_{snap:04d}.csv"
        no_gradh_file = no_gradh_dir / f"snapshot_{snap:04d}.csv"
        
        if gradh_file.exists():
            data = load_snapshot(gradh_file)
            center, _ = analyze_center_particles(data)
            rho_center_gradh.append(center['dens'].mean())
            omega_center_gradh.append(center['gradh'].mean())
            h_center_gradh.append(center['sml'].mean())
            times_gradh.append(snap)
        
        if no_gradh_file.exists():
            data = load_snapshot(no_gradh_file)
            center, _ = analyze_center_particles(data)
            rho_center_no_gradh.append(center['dens'].mean())
            h_center_no_gradh.append(center['sml'].mean())
            times_no_gradh.append(snap)
    
    print(f"\n{'Snapshot':<10} {'ρ(gradh)':<12} {'ρ(no_gradh)':<12} {'Ω(gradh)':<12} {'h(gradh)':<12} {'h(no_gradh)':<12}")
    print("-" * 70)
    for i, snap in enumerate(times_gradh):
        if i < len(times_no_gradh):
            print(f"{snap:<10} {rho_center_gradh[i]:<12.6f} {rho_center_no_gradh[i]:<12.6f} {omega_center_gradh[i]:<12.6f} {h_center_gradh[i]:<12.6f} {h_center_no_gradh[i]:<12.6f}")
    
    # Analyze the relationship between Ω and density/h
    print("\n" + "=" * 80)
    print("THEORETICAL ANALYSIS: Ω vs ρ RELATIONSHIP")
    print("=" * 80)
    
    # From the code: Ω = 1 / (1 + h/(D*ρ) * dρ/dh)
    # For h ∝ ρ^(-1/3), we have dh/dρ = -h/(3ρ)
    # And dρ/dh = -3ρ/h
    # So: h/(D*ρ) * dρ/dh = h/(3ρ) * (-3ρ/h) = -1
    # Therefore Ω = 1/(1-1) = undefined? No, let's check the actual formula
    
    # Actually dρ/dh comes from the SPH sum: ρ = Σ m_j W(r_ij, h)
    # dρ/dh = Σ m_j ∂W/∂h
    
    print("\nFrom the code:")
    print("  Ω = 1 / (1 + h/(D*ρ) * dρ/dh)")
    print("  where dρ/dh = Σ m_j ∂W(r_ij, h)/∂h")
    
    print("\nObserved Ω values at center:")
    for i, snap in enumerate(times_gradh[:3]):
        print(f"  Snapshot {snap}: Ω = {omega_center_gradh[i]:.6f}, ρ = {rho_center_gradh[i]:.6f}, h = {h_center_gradh[i]:.6f}")
    
    # Check if Ω * |∇W| ≈ const
    # |∇W| ∝ h^(-4) ∝ ρ^(4/3)
    print("\n" + "=" * 80)
    print("CHECK: Does Ω * |∇W| remain constant?")
    print("=" * 80)
    print("If Ω cancels the h^(-4) scaling, then Ω * h^(-4) should be constant")
    
    print(f"\n{'Snapshot':<10} {'Ω':<12} {'h':<12} {'h^(-4)':<15} {'Ω * h^(-4)':<15}")
    print("-" * 65)
    for i, snap in enumerate(times_gradh[:6]):
        h = h_center_gradh[i]
        omega = omega_center_gradh[i]
        h_inv4 = h**(-4)
        product = omega * h_inv4
        print(f"{snap:<10} {omega:<12.6f} {h:<12.6f} {h_inv4:<15.2f} {product:<15.2f}")
    
    # Normalize to initial value
    if len(times_gradh) > 0:
        h0 = h_center_gradh[0]
        omega0 = omega_center_gradh[0]
        product0 = omega0 * h0**(-4)
        
        print(f"\nNormalized to initial (Ω₀*h₀⁻⁴ = {product0:.2f}):")
        print(f"{'Snapshot':<10} {'(Ω*h⁻⁴)/(Ω₀*h₀⁻⁴)':<20}")
        print("-" * 30)
        for i, snap in enumerate(times_gradh[:6]):
            h = h_center_gradh[i]
            omega = omega_center_gradh[i]
            ratio = (omega * h**(-4)) / product0
            print(f"{snap:<10} {ratio:<20.6f}")
    
    # Calculate the effective force response
    print("\n" + "=" * 80)
    print("FORCE RESPONSE ANALYSIS")
    print("=" * 80)
    
    # GSPH force ∝ p* * Ω * |∇W| / ρ²
    # With Ω*|∇W| ≈ const, force ∝ p* / ρ²
    # Since p* ∝ ρ^(γ/2) for small perturbations (from Riemann solver)
    # Force ∝ ρ^(γ/2) / ρ² = ρ^(γ/2 - 2)
    
    # For γ = 5/3: γ/2 - 2 = 5/6 - 2 = -7/6
    # So force DECREASES with increasing density!
    
    gamma = 5/3
    print(f"\nWith Ω correction (assuming Ω*|∇W| = const):")
    print(f"  F ∝ p* / ρ² ∝ ρ^(γ/2 - 2) = ρ^({gamma/2 - 2:.4f})")
    print(f"  For γ = 5/3: exponent = {gamma/2 - 2:.4f} = -7/6")
    print(f"  → Force DECREASES as density increases!")
    
    print(f"\nWithout Ω correction (Ω = 1):")
    print(f"  F ∝ p* * |∇W| / ρ² ∝ ρ^(γ/2) * ρ^(4/3) / ρ² = ρ^(γ/2 + 4/3 - 2)")
    print(f"  For γ = 5/3: exponent = 5/6 + 4/3 - 2 = {5/6 + 4/3 - 2:.4f} = 1/6")
    print(f"  → Force increases weakly with density (coefficient 1/6)")
    
    print("\n" + "=" * 80)
    print("STABILITY CRITERION")
    print("=" * 80)
    print("For stability against gravity, need δF_P/F_P ≥ δρ/ρ")
    print("i.e., force response coefficient C ≥ 1")
    print(f"\n  Without Ω: C = 1/6 < 1 → UNSTABLE")
    print(f"  With Ω: C = -7/6 < 1 → Also appears unstable by this criterion!")
    
    print("\n" + "=" * 80)
    print("BUT SIMULATIONS SHOW GSPH+Ω IS STABLE!")
    print("=" * 80)
    print("This suggests the simple single-particle analysis is incomplete.")
    print("The stability must come from the FULL Lagrangian formulation.")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Density evolution
    ax1 = axes[0, 0]
    ax1.plot(times_gradh, rho_center_gradh, 'b-o', label='gradh (Ω enabled)', markersize=4)
    ax1.plot(times_no_gradh, rho_center_no_gradh, 'r-s', label='no_gradh (Ω=1)', markersize=4)
    ax1.set_xlabel('Snapshot')
    ax1.set_ylabel('Center Density ρ')
    ax1.set_title('Density Evolution at Center')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Omega evolution
    ax2 = axes[0, 1]
    ax2.plot(times_gradh, omega_center_gradh, 'b-o', markersize=4)
    ax2.axhline(y=1.0, color='r', linestyle='--', label='Ω=1 (no correction)')
    ax2.set_xlabel('Snapshot')
    ax2.set_ylabel('Ω (grad-h factor)')
    ax2.set_title('Ω Evolution at Center')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: h evolution
    ax3 = axes[1, 0]
    ax3.plot(times_gradh, h_center_gradh, 'b-o', label='gradh', markersize=4)
    ax3.plot(times_no_gradh, h_center_no_gradh, 'r-s', label='no_gradh', markersize=4)
    ax3.set_xlabel('Snapshot')
    ax3.set_ylabel('Smoothing Length h')
    ax3.set_title('Smoothing Length Evolution at Center')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Ω * h^(-4) (should be ~constant if theory correct)
    ax4 = axes[1, 1]
    omega_h4_gradh = [omega_center_gradh[i] * h_center_gradh[i]**(-4) for i in range(len(times_gradh))]
    h4_no_gradh = [1.0 * h_center_no_gradh[i]**(-4) for i in range(len(times_no_gradh))]
    
    # Normalize
    if omega_h4_gradh:
        omega_h4_gradh_norm = [x / omega_h4_gradh[0] for x in omega_h4_gradh]
        h4_no_gradh_norm = [x / h4_no_gradh[0] for x in h4_no_gradh]
        
        ax4.plot(times_gradh, omega_h4_gradh_norm, 'b-o', label='Ω*h⁻⁴ (gradh)', markersize=4)
        ax4.plot(times_no_gradh, h4_no_gradh_norm, 'r-s', label='h⁻⁴ (no_gradh)', markersize=4)
        ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
        ax4.set_xlabel('Snapshot')
        ax4.set_ylabel('Normalized Value')
        ax4.set_title('Ω*|∇W| Scaling Check (normalized to t=0)')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/Users/guo/Downloads/sphcode/simulations/stability/gradh_study_3d/results/gradh_analysis.png', dpi=150)
    print(f"\nPlot saved to: results/gradh_analysis.png")
    
    plt.show()

if __name__ == "__main__":
    main()
