#!/usr/bin/env python3
"""
analyze_polytrope_stability.py

Compare grad-h vs no-grad-h for 1D polytropic slab hydrostatic equilibrium.
This tests whether grad-h correction prevents artificial core collapse.

Key metric: Central density evolution over time
- Stable: ρ_center ≈ ρ_center(t=0)
- Collapse: ρ_center increases over time
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob

def read_csv(filepath):
    """Read SPH CSV output (skip comment lines)"""
    return pd.read_csv(filepath, comment='#')

def get_central_density(df):
    """Get density at x=0 (center of slab)"""
    x = df['pos_x'].values
    rho = df['dens'].values
    
    # Find particle closest to center
    idx = np.argmin(np.abs(x))
    return rho[idx]

def get_density_profile(df):
    """Get sorted density profile"""
    x = df['pos_x'].values
    rho = df['dens'].values
    
    # Sort by position
    order = np.argsort(x)
    return x[order], rho[order]

def analyze_stability(results_dir, label, output_time=0.1):
    """Analyze density evolution for stability"""
    csv_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"), 
                       key=lambda x: int(Path(x).stem.split('_')[-1]))
    
    if not csv_files:
        print(f"No CSV files found in {results_dir}")
        return None, None, None
    
    times = []
    central_densities = []
    max_densities = []
    
    for csv_file in csv_files:
        # Extract time from filename
        idx = int(Path(csv_file).stem.split('_')[-1])
        t = idx * output_time
        
        df = read_csv(csv_file)
        rho_center = get_central_density(df)
        rho_max = df['dens'].max()
        
        times.append(t)
        central_densities.append(rho_center)
        max_densities.append(rho_max)
    
    return np.array(times), np.array(central_densities), np.array(max_densities)

def main():
    base_dir = Path("results/gradh_test")
    
    # Analyze both cases
    t_gradh, rho_c_gradh, rho_max_gradh = analyze_stability(
        base_dir / "polytrope_gradh", "grad-h", output_time=0.1)
    
    t_nogradh, rho_c_nogradh, rho_max_nogradh = analyze_stability(
        base_dir / "polytrope_nograd", "no-grad-h", output_time=0.1)
    
    if t_gradh is None or t_nogradh is None:
        print("Failed to load data")
        return
    
    # Initial central density (should be 1.0)
    rho_c_init = rho_c_gradh[0]
    
    print("=" * 60)
    print("Polytropic Slab Stability Analysis")
    print("=" * 60)
    print(f"Initial central density: ρ_c(t=0) = {rho_c_init:.4f}")
    print()
    
    # Calculate collapse metrics
    # Relative change in central density
    delta_rho_gradh = (rho_c_gradh[-1] - rho_c_gradh[0]) / rho_c_gradh[0] * 100
    delta_rho_nogradh = (rho_c_nogradh[-1] - rho_c_nogradh[0]) / rho_c_nogradh[0] * 100
    
    print("Central Density Change (t=0 → t=5):")
    print(f"  With grad-h:    Δρ_c/ρ_c = {delta_rho_gradh:+.1f}%")
    print(f"  Without grad-h: Δρ_c/ρ_c = {delta_rho_nogradh:+.1f}%")
    print()
    
    # Maximum density change
    delta_max_gradh = (rho_max_gradh[-1] - rho_max_gradh[0]) / rho_max_gradh[0] * 100
    delta_max_nogradh = (rho_max_nogradh[-1] - rho_max_nogradh[0]) / rho_max_nogradh[0] * 100
    
    print("Maximum Density Change (t=0 → t=5):")
    print(f"  With grad-h:    Δρ_max/ρ_max = {delta_max_gradh:+.1f}%")
    print(f"  Without grad-h: Δρ_max/ρ_max = {delta_max_nogradh:+.1f}%")
    print()
    
    # Interpret results
    print("=" * 60)
    print("INTERPRETATION:")
    print("=" * 60)
    if abs(delta_rho_gradh) < 10 and abs(delta_rho_nogradh) > 20:
        print("✓ Grad-h PREVENTS collapse: Central density stable")
        print("✗ No-grad-h shows COLLAPSE: Central density increases")
    elif abs(delta_rho_gradh) < 10 and abs(delta_rho_nogradh) < 10:
        print("Both cases stable (1D may be too stable to show effect)")
    elif delta_rho_gradh > 0 and delta_rho_nogradh > 0:
        if delta_rho_nogradh > delta_rho_gradh * 1.5:
            print("Grad-h REDUCES collapse rate compared to no-grad-h")
        else:
            print("Both cases show some collapse (similar rates)")
    else:
        print("Results inconclusive - check density profiles")
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Central density vs time
    ax1 = axes[0, 0]
    ax1.plot(t_gradh, rho_c_gradh, 'b-o', label='grad-h', markersize=4)
    ax1.plot(t_nogradh, rho_c_nogradh, 'r-s', label='no-grad-h', markersize=4)
    ax1.axhline(rho_c_init, color='k', linestyle='--', alpha=0.5, label='Initial')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Central Density ρ_c')
    ax1.set_title('Central Density Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Normalized central density
    ax2 = axes[0, 1]
    ax2.plot(t_gradh, rho_c_gradh / rho_c_init, 'b-o', label='grad-h', markersize=4)
    ax2.plot(t_nogradh, rho_c_nogradh / rho_c_nogradh[0], 'r-s', label='no-grad-h', markersize=4)
    ax2.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('ρ_c(t) / ρ_c(0)')
    ax2.set_title('Normalized Central Density')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Max density vs time
    ax3 = axes[1, 0]
    ax3.plot(t_gradh, rho_max_gradh, 'b-o', label='grad-h', markersize=4)
    ax3.plot(t_nogradh, rho_max_nogradh, 'r-s', label='no-grad-h', markersize=4)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Maximum Density ρ_max')
    ax3.set_title('Maximum Density Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Final density profiles
    ax4 = axes[1, 1]
    
    # Read final snapshots
    gradh_final = sorted(glob.glob("results/gradh_test/polytrope_gradh/snapshot_*.csv"))[-1]
    nogradh_final = sorted(glob.glob("results/gradh_test/polytrope_nograd/snapshot_*.csv"))[-1]
    gradh_init = sorted(glob.glob("results/gradh_test/polytrope_gradh/snapshot_*.csv"))[0]
    
    df_gradh_final = read_csv(gradh_final)
    df_nogradh_final = read_csv(nogradh_final)
    df_init = read_csv(gradh_init)
    
    x_init, rho_init = get_density_profile(df_init)
    x_gradh, rho_gradh = get_density_profile(df_gradh_final)
    x_nogradh, rho_nogradh = get_density_profile(df_nogradh_final)
    
    ax4.plot(x_init, rho_init, 'k--', label='Initial (t=0)', linewidth=2)
    ax4.plot(x_gradh, rho_gradh, 'b-', label='grad-h (t=5)', linewidth=1.5)
    ax4.plot(x_nogradh, rho_nogradh, 'r-', label='no-grad-h (t=5)', linewidth=1.5)
    ax4.set_xlabel('Position x')
    ax4.set_ylabel('Density ρ')
    ax4.set_title('Final Density Profiles')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/gradh_test/polytrope_stability_comparison.png', dpi=150)
    print(f"\nPlot saved to: results/gradh_test/polytrope_stability_comparison.png")
    plt.show()

if __name__ == "__main__":
    main()
