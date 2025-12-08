#!/usr/bin/env python3
"""
analyze_jeans_instability.py

Analyze Jeans instability simulation results to compare grad-h vs no-grad-h.

Theory:
- Jeans instability: ω² = c_s²k² - 4πGρ₀
- For λ > λ_J, ω² < 0 and perturbations grow as exp(Γt)
- Growth rate Γ = sqrt(-ω²)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob

def read_csv(filepath):
    """Read SPH CSV output (skip comment lines starting with #)"""
    return pd.read_csv(filepath, comment='#')

def compute_amplitude(df, wavelength):
    """
    Compute the amplitude of the sinusoidal perturbation
    by fitting to ρ(x) = ρ₀(1 + A·sin(kx + φ))
    """
    x = df['pos_x'].values
    rho = df['dens'].values
    
    rho_mean = np.mean(rho)
    rho_normalized = (rho - rho_mean) / rho_mean
    
    # Fit: A*sin(kx + φ) = A*sin(φ)*cos(kx) + A*cos(φ)*sin(kx)
    k = 2 * np.pi / wavelength
    sin_kx = np.sin(k * x)
    cos_kx = np.cos(k * x)
    
    # Least squares fit
    A_sin = np.mean(rho_normalized * sin_kx) * 2  # A*cos(φ)
    A_cos = np.mean(rho_normalized * cos_kx) * 2  # A*sin(φ)
    
    A = np.sqrt(A_sin**2 + A_cos**2)
    return A, rho_mean

def compute_density_contrast(df, wavelength):
    """
    Compute max-min density contrast
    """
    rho = df['dens'].values
    return (np.max(rho) - np.min(rho)) / (2 * np.mean(rho))

def analyze_jeans(results_dir, label, wavelength, ax_amp, ax_growth, color):
    """Analyze a single simulation"""
    csv_files = sorted(glob.glob(f"{results_dir}/snapshot_*.csv"))
    
    if not csv_files:
        print(f"No CSV files found in {results_dir}")
        return None, None
    
    times = []
    amplitudes = []
    contrasts = []
    
    for csv_file in csv_files:
        # Extract time from filename
        idx = int(Path(csv_file).stem.split('_')[-1])
        t = idx * 0.02  # output_time from config
        
        df = read_csv(csv_file)
        A, rho_mean = compute_amplitude(df, wavelength)
        contrast = compute_density_contrast(df, wavelength)
        
        times.append(t)
        amplitudes.append(A)
        contrasts.append(contrast)
    
    times = np.array(times)
    amplitudes = np.array(amplitudes)
    contrasts = np.array(contrasts)
    
    # Plot amplitude vs time
    ax_amp.semilogy(times, amplitudes, 'o-', color=color, label=label, markersize=4)
    
    # Fit exponential growth in early regime
    # A(t) = A₀ * exp(Γt)
    # log(A) = log(A₀) + Γ*t
    
    # Use only early times where amplitude < 0.5 (linear regime)
    mask = (amplitudes > 0.001) & (amplitudes < 0.5) & (times > 0.1)
    if np.sum(mask) > 3:
        log_A = np.log(amplitudes[mask])
        t_fit = times[mask]
        
        # Linear fit
        coeffs = np.polyfit(t_fit, log_A, 1)
        Gamma_measured = coeffs[0]
        A0_measured = np.exp(coeffs[1])
        
        # Plot fit
        t_theory = np.linspace(times[0], times[-1], 100)
        A_fit = A0_measured * np.exp(Gamma_measured * t_theory)
        ax_amp.semilogy(t_theory, A_fit, '--', color=color, alpha=0.7, 
                        label=f'{label} fit: Γ={Gamma_measured:.3f}')
        
        print(f"{label}:")
        print(f"  Measured growth rate Γ = {Gamma_measured:.4f}")
        print(f"  Initial amplitude A₀ = {A0_measured:.4f}")
        
        return Gamma_measured, A0_measured
    
    return None, None

def main():
    # Parameters from config
    wavelength = 4.0
    rho_0 = 1.0
    c_s = 1.0
    G = 1.0
    gamma = 1.4
    
    # Theoretical values
    k = 2 * np.pi / wavelength
    lambda_J = c_s * np.sqrt(np.pi / (G * rho_0))
    omega_sq = c_s**2 * k**2 - 4 * np.pi * G * rho_0
    
    print("=" * 60)
    print("Jeans Instability Analysis")
    print("=" * 60)
    print(f"Wavelength λ = {wavelength}")
    print(f"Jeans length λ_J = {lambda_J:.4f}")
    print(f"λ/λ_J = {wavelength/lambda_J:.4f}")
    print(f"ω² = {omega_sq:.4f}")
    
    if omega_sq < 0:
        Gamma_theory = np.sqrt(-omega_sq)
        print(f"UNSTABLE: Theoretical growth rate Γ = {Gamma_theory:.4f}")
        print(f"e-folding time τ = {1/Gamma_theory:.4f}")
    else:
        print(f"STABLE: ω = {np.sqrt(omega_sq):.4f}")
        Gamma_theory = None
    
    print("=" * 60)
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_amp = axes[0]
    ax_growth = axes[1]
    
    # Analyze both cases
    results = {}
    base_dir = Path("sample/jeans_instability/results")
    
    Gamma_gradh, A0_gradh = analyze_jeans(
        base_dir / "gradh", "grad-h", wavelength, ax_amp, ax_growth, 'blue')
    
    Gamma_nogradh, A0_nogradh = analyze_jeans(
        base_dir / "nogradh", "no-grad-h", wavelength, ax_amp, ax_growth, 'red')
    
    # Add theoretical line
    if Gamma_theory:
        t_theory = np.linspace(0, 3, 100)
        A_theory = 0.01 * np.exp(Gamma_theory * t_theory)
        ax_amp.semilogy(t_theory, A_theory, 'k--', linewidth=2, alpha=0.5,
                        label=f'Theory: Γ={Gamma_theory:.3f}')
    
    ax_amp.set_xlabel('Time')
    ax_amp.set_ylabel('Perturbation Amplitude A')
    ax_amp.set_title('Jeans Instability: Amplitude Growth')
    ax_amp.legend()
    ax_amp.grid(True, alpha=0.3)
    ax_amp.set_ylim(1e-3, 10)
    
    # Bar chart comparison
    methods = []
    gammas = []
    colors = []
    
    if Gamma_theory:
        methods.append('Theory')
        gammas.append(Gamma_theory)
        colors.append('black')
    
    if Gamma_gradh:
        methods.append('grad-h')
        gammas.append(Gamma_gradh)
        colors.append('blue')
    
    if Gamma_nogradh:
        methods.append('no-grad-h')
        gammas.append(Gamma_nogradh)
        colors.append('red')
    
    if gammas:
        x = np.arange(len(methods))
        bars = ax_growth.bar(x, gammas, color=colors, alpha=0.7)
        ax_growth.set_xticks(x)
        ax_growth.set_xticklabels(methods)
        ax_growth.set_ylabel('Growth Rate Γ')
        ax_growth.set_title('Growth Rate Comparison')
        ax_growth.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, g in zip(bars, gammas):
            ax_growth.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                          f'{g:.3f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('sample/jeans_instability/jeans_analysis.png', dpi=150)
    plt.show()
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    if Gamma_theory:
        print(f"Theoretical growth rate: Γ = {Gamma_theory:.4f}")
    if Gamma_gradh:
        print(f"With grad-h:    Γ = {Gamma_gradh:.4f}", end="")
        if Gamma_theory:
            print(f" (ratio to theory: {Gamma_gradh/Gamma_theory:.3f})")
        else:
            print()
    if Gamma_nogradh:
        print(f"Without grad-h: Γ = {Gamma_nogradh:.4f}", end="")
        if Gamma_theory:
            print(f" (ratio to theory: {Gamma_nogradh/Gamma_theory:.3f})")
        else:
            print()
    
    if Gamma_gradh and Gamma_nogradh:
        print(f"\nDifference: {abs(Gamma_gradh - Gamma_nogradh):.4f}")
        print(f"Ratio (gradh/nogradh): {Gamma_gradh/Gamma_nogradh:.4f}")

if __name__ == "__main__":
    main()
