#!/usr/bin/env python3
"""
Analyze sinusoidal perturbation decay for diffusion coefficient measurement.

Theory:
    ρ(x,t) = ρ_mean * (1 + A(t)*sin(2πx/λ))
    A(t) = A₀ * exp(-Γt)
    Γ = D_eff * k²  where k = 2π/λ
    D_eff ≈ ε * c_s * h  (for GSPH without grad-h)

By measuring A(t), we can extract:
    Γ = -d(ln A)/dt
    D_eff = Γ / k²
    ε = D_eff / (c_s * h)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

def extract_amplitude(df, wavelength=1.0):
    """Extract sinusoidal amplitude from density profile."""
    x = df['pos_x'].values
    rho = df['dens'].values
    
    rho_mean = np.mean(rho)
    
    # Fit sine wave: ρ = ρ_mean * (1 + A*sin(2πx/λ + φ))
    # Using Fourier analysis
    k = 2 * np.pi / wavelength
    
    # Compute Fourier coefficients
    cos_coeff = np.mean((rho / rho_mean - 1) * np.cos(k * x)) * 2
    sin_coeff = np.mean((rho / rho_mean - 1) * np.sin(k * x)) * 2
    
    # Amplitude is sqrt(a² + b²)
    amplitude = np.sqrt(cos_coeff**2 + sin_coeff**2)
    
    return amplitude, rho_mean

def main():
    if len(sys.argv) < 2:
        results_dir = Path("results/sinusoidal_test")
    else:
        results_dir = Path(sys.argv[1])
    
    # Find all snapshot files
    snapshots = sorted(glob.glob(str(results_dir / "snapshot_*.csv")))
    
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return
    
    print(f"Found {len(snapshots)} snapshots")
    
    # Parameters (should match config)
    wavelength = 1.0
    gamma = 1.4
    
    # Extract amplitude from each snapshot
    times = []
    amplitudes = []
    
    for snap in snapshots:
        # Read time from header (# Time (code): value)
        t = None
        with open(snap, 'r') as f:
            for line in f:
                if line.startswith('# Time (code):'):
                    t = float(line.split(':')[1].strip())
                    break
        
        if t is None:
            # Fallback to filename-based time
            snap_num = int(Path(snap).stem.split('_')[1])
            t = snap_num * 0.1  # Default
        
        df = pd.read_csv(snap, comment='#')
        
        A, rho_mean = extract_amplitude(df, wavelength)
        
        times.append(t)
        amplitudes.append(A)
        
        if snap == snapshots[0] or snap == snapshots[-1]:
            print(f"  {Path(snap).name}: t={t:.3f}, A={A:.6f}, rho_mean={rho_mean:.4f}")
    
    times = np.array(times)
    amplitudes = np.array(amplitudes)
    
    # Fit exponential decay: A = A0 * exp(-Γ*t)
    # Linear fit to ln(A) vs t
    valid = amplitudes > 0.001  # Filter near-zero amplitudes
    if np.sum(valid) < 3:
        print("Not enough valid data points for fitting")
        return
    
    ln_A = np.log(amplitudes[valid])
    t_fit = times[valid]
    
    # Linear fit: ln(A) = ln(A0) - Γ*t
    coeffs = np.polyfit(t_fit, ln_A, 1)
    Gamma_measured = -coeffs[0]
    A0_fit = np.exp(coeffs[1])
    
    # Compute theory prediction
    # Need to estimate h from data
    df0 = pd.read_csv(snapshots[0], comment='#')
    h_avg = df0['sml'].mean()
    rho_mean = df0['dens'].mean()
    P_avg = df0['pres'].mean()
    c_s = np.sqrt(gamma * P_avg / rho_mean)
    
    k = 2 * np.pi / wavelength
    
    # Theory: D_eff = ε * c_s * h, Γ = D * k²
    D_measured = Gamma_measured / (k**2)
    epsilon_measured = D_measured / (c_s * h_avg)
    
    print("\n=== Diffusion Analysis Results ===")
    print(f"Sound speed c_s = {c_s:.4f}")
    print(f"Average smoothing length h = {h_avg:.4f}")
    print(f"Wavenumber k = {k:.4f}")
    print(f"\nMeasured decay rate Γ = {Gamma_measured:.4f}")
    print(f"Measured diffusivity D = Γ/k² = {D_measured:.4f}")
    print(f"Measured ε = D/(c_s·h) = {epsilon_measured:.4f}")
    print(f"\nExpected ε ≈ 0.35-0.40 for GSPH without grad-h")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot amplitude vs time
    ax1 = axes[0]
    ax1.semilogy(times, amplitudes, 'bo-', label='Measured')
    t_theory = np.linspace(0, times[-1], 100)
    ax1.semilogy(t_theory, A0_fit * np.exp(-Gamma_measured * t_theory), 'r--', 
                 label=f'Fit: Γ={Gamma_measured:.3f}')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Amplitude A(t)')
    ax1.set_title('Amplitude Decay')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot ln(A) vs time (should be linear)
    ax2 = axes[1]
    ax2.plot(times[valid], ln_A, 'bo-', label='Measured')
    ax2.plot(t_fit, coeffs[0] * t_fit + coeffs[1], 'r--', 
             label=f'Linear fit: slope=-{Gamma_measured:.3f}')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('ln(A)')
    ax2.set_title('Log Amplitude (should be linear)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = results_dir / "diffusion_analysis.png"
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to {output_file}")
    
    # Save results
    results_file = results_dir / "diffusion_results.txt"
    with open(results_file, 'w') as f:
        f.write("=== Diffusion Analysis Results ===\n")
        f.write(f"Sound speed c_s = {c_s:.6f}\n")
        f.write(f"Average smoothing length h = {h_avg:.6f}\n")
        f.write(f"Wavenumber k = {k:.6f}\n")
        f.write(f"\nMeasured decay rate Γ = {Gamma_measured:.6f}\n")
        f.write(f"Measured diffusivity D = Γ/k² = {D_measured:.6f}\n")
        f.write(f"Measured ε = D/(c_s·h) = {epsilon_measured:.6f}\n")
        f.write(f"\nExpected ε ≈ 0.35-0.40 for GSPH without grad-h\n")
        f.write(f"\nTime,Amplitude\n")
        for t, A in zip(times, amplitudes):
            f.write(f"{t:.6f},{A:.6f}\n")
    
    print(f"Results saved to {results_file}")
    
    plt.show()

if __name__ == "__main__":
    main()
