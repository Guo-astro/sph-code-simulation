#!/usr/bin/env python3
"""
Jeans Instability Wavelength Comparison Analysis

This script analyzes and visualizes the results of Jeans instability simulations
with different wavelengths (long λ=4.0 vs short λ=1.0) and grad-h settings.

Theory:
- Jeans length: λ_J = c_s √(π/Gρ₀) ≈ 1.77 for our parameters
- Long wavelength (λ=4.0 > λ_J): gravitationally unstable → exponential growth
- Short wavelength (λ=1.0 < λ_J): pressure-supported → stable oscillations

Expected Results:
- Long wave with grad-h: matches Jeans theory exactly
- Long wave without grad-h: faster/different growth rate (extra error from missing Ω)
- Short wave: stable oscillations regardless of grad-h (pressure dominates)

Usage:
    python analyze_wavelength_comparison.py

Output:
    - jeans_wavelength_comparison.png: Multi-panel comparison plot
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Physics parameters (matching config files)
RHO_0 = 1.0
C_S = 1.0
G = 1.0

# Jeans length and critical wavenumber
LAMBDA_J = C_S * np.sqrt(np.pi / (G * RHO_0))
K_J = 2 * np.pi / LAMBDA_J

print(f"Jeans Length: λ_J = {LAMBDA_J:.3f}")
print(f"Jeans Wavenumber: k_J = {K_J:.3f}")
print()


def load_snapshot(filepath):
    """Load a CSV snapshot and return position and density."""
    try:
        # Find header line (starts with 'id,')
        with open(filepath) as f:
            header_line = 0
            for i, line in enumerate(f):
                if line.startswith('id,'):
                    header_line = i
                    break
        
        import pandas as pd
        df = pd.read_csv(filepath, skiprows=header_line)
        x = df['pos_x'].values
        rho = df['dens'].values
        return x, rho
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None, None


def get_density_amplitude(results_dir, wavelength):
    """
    Extract density perturbation amplitude vs time for a simulation.
    
    Returns:
        times: array of simulation times
        amplitudes: array of density perturbation amplitudes (δρ/ρ₀)
    """
    results_path = Path(results_dir)
    snapshots = sorted(results_path.glob("snapshot_*.csv"))
    
    if not snapshots:
        print(f"No snapshots found in {results_dir}")
        return np.array([]), np.array([])
    
    times = []
    amplitudes = []
    k = 2 * np.pi / wavelength
    
    for snap in snapshots:
        # Extract time from filename (snapshot_XXXX.csv)
        try:
            idx = int(snap.stem.split('_')[1])
            time = idx * 0.02  # outputTime from config
        except:
            continue
        
        x, rho = load_snapshot(snap)
        if x is None:
            continue
        
        # Compute Fourier amplitude at the perturbation wavenumber
        # Using least squares fit to sinusoidal perturbation
        delta_rho = rho - RHO_0
        
        # Fit: δρ = A·sin(kx) + B·cos(kx)
        sin_kx = np.sin(k * x)
        cos_kx = np.cos(k * x)
        
        # Least squares
        A = 2 * np.mean(delta_rho * sin_kx)
        B = 2 * np.mean(delta_rho * cos_kx)
        
        amplitude = np.sqrt(A**2 + B**2) / RHO_0
        
        times.append(time)
        amplitudes.append(amplitude)
    
    return np.array(times), np.array(amplitudes)


def compute_jeans_theory(times, wavelength, initial_amplitude=0.01):
    """
    Compute theoretical Jeans instability growth.
    
    For wavelength > λ_J: exponential growth with rate ω = √(4πGρ₀ - c_s²k²)
    For wavelength < λ_J: oscillation with frequency ω = √(c_s²k² - 4πGρ₀)
    """
    k = 2 * np.pi / wavelength
    omega_sq = C_S**2 * k**2 - 4 * np.pi * G * RHO_0
    
    if omega_sq < 0:
        # Unstable: exponential growth
        gamma = np.sqrt(-omega_sq)
        return initial_amplitude * np.cosh(gamma * times)
    else:
        # Stable: oscillation
        omega = np.sqrt(omega_sq)
        return initial_amplitude * np.abs(np.cos(omega * times))


def main():
    base_dir = Path("/Users/guo/Downloads/sphcode/simulations/stability/jeans_instability/results")
    
    # Define cases
    cases = {
        'long_gradh': {'dir': 'long_wave_gradh', 'wavelength': 4.0, 'label': 'λ=4.0 (with grad-h)'},
        'long_nogradh': {'dir': 'long_wave_nogradh', 'wavelength': 4.0, 'label': 'λ=4.0 (no grad-h)'},
        'short_gradh': {'dir': 'short_wave_gradh', 'wavelength': 1.0, 'label': 'λ=1.0 (with grad-h)'},
        'short_nogradh': {'dir': 'short_wave_nogradh', 'wavelength': 1.0, 'label': 'λ=1.0 (no grad-h)'},
    }
    
    # Load all data
    results = {}
    for key, case in cases.items():
        results_dir = base_dir / case['dir']
        if results_dir.exists():
            times, amplitudes = get_density_amplitude(results_dir, case['wavelength'])
            results[key] = {'times': times, 'amplitudes': amplitudes, **case}
            print(f"Loaded {key}: {len(times)} snapshots")
        else:
            print(f"Warning: {results_dir} not found")
    
    # Create figure
    fig = plt.figure(figsize=(14, 10))
    
    # Panel 1: Long wavelength comparison (log scale)
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_title(f'Long Wavelength (λ=4.0 > λ_J={LAMBDA_J:.2f})\nGravitationally Unstable', fontsize=11)
    
    for key in ['long_gradh', 'long_nogradh']:
        if key in results:
            r = results[key]
            style = '-' if 'gradh' in key and 'no' not in key else '--'
            ax1.semilogy(r['times'], r['amplitudes'], style, linewidth=2, label=r['label'])
    
    # Theoretical prediction
    if 'long_gradh' in results:
        t_theory = np.linspace(0, results['long_gradh']['times'].max(), 200)
        amp_theory = compute_jeans_theory(t_theory, 4.0)
        ax1.semilogy(t_theory, amp_theory, 'k:', linewidth=1.5, alpha=0.7, label='Jeans theory')
    
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Density Amplitude δρ/ρ₀')
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0.005, 10])
    
    # Panel 2: Short wavelength comparison (linear scale)
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.set_title(f'Short Wavelength (λ=1.0 < λ_J={LAMBDA_J:.2f})\nPressure-Supported (Stable)', fontsize=11)
    
    for key in ['short_gradh', 'short_nogradh']:
        if key in results:
            r = results[key]
            style = '-' if 'gradh' in key and 'no' not in key else '--'
            ax2.plot(r['times'], r['amplitudes'], style, linewidth=2, label=r['label'])
    
    # Theoretical prediction
    if 'short_gradh' in results:
        t_theory = np.linspace(0, results['short_gradh']['times'].max(), 200)
        amp_theory = compute_jeans_theory(t_theory, 1.0)
        ax2.plot(t_theory, amp_theory, 'k:', linewidth=1.5, alpha=0.7, label='Theory (oscillation)')
    
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Density Amplitude δρ/ρ₀')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Grad-h effect on long wavelength (difference)
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_title('Grad-h Effect on Unstable Mode (λ=4.0)', fontsize=11)
    
    if 'long_gradh' in results and 'long_nogradh' in results:
        r1 = results['long_gradh']
        r2 = results['long_nogradh']
        
        # Find common time range
        min_len = min(len(r1['times']), len(r2['times']))
        t = r1['times'][:min_len]
        
        ax3.semilogy(t, r1['amplitudes'][:min_len], 'b-', linewidth=2, label='With grad-h')
        ax3.semilogy(t, r2['amplitudes'][:min_len], 'r--', linewidth=2, label='Without grad-h')
        
        # Compute ratio
        ratio = r2['amplitudes'][:min_len] / r1['amplitudes'][:min_len]
        ax3_twin = ax3.twinx()
        ax3_twin.plot(t, ratio, 'g-', linewidth=1.5, alpha=0.6)
        ax3_twin.set_ylabel('Ratio (no-gradh / gradh)', color='g')
        ax3_twin.tick_params(axis='y', labelcolor='g')
    
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Density Amplitude δρ/ρ₀')
    ax3.legend(loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Density profiles at late time
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.set_title('Density Profiles at Late Time', fontsize=11)
    
    # Load final snapshots
    for key, color in [('long_gradh', 'blue'), ('long_nogradh', 'red'), 
                       ('short_gradh', 'green'), ('short_nogradh', 'orange')]:
        if key in results:
            results_dir = base_dir / cases[key]['dir']
            snapshots = sorted(results_dir.glob("snapshot_*.csv"))
            if snapshots:
                # Use a mid-time snapshot to avoid extreme values
                snap_idx = len(snapshots) // 2
                x, rho = load_snapshot(snapshots[snap_idx])
                if x is not None:
                    # Sort by x for proper plotting
                    sort_idx = np.argsort(x)
                    style = '-' if 'gradh' in key and 'no' not in key else '--'
                    ax4.plot(x[sort_idx], rho[sort_idx], style, color=color, 
                            linewidth=1.5, alpha=0.8, label=results[key]['label'])
    
    ax4.axhline(y=RHO_0, color='k', linestyle=':', alpha=0.5, label='ρ₀')
    ax4.set_xlabel('Position x')
    ax4.set_ylabel('Density ρ')
    ax4.legend(loc='upper right', fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_path = base_dir.parent / 'jeans_wavelength_comparison.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_path}")
    
    # Print analysis summary
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    
    print(f"\nPhysics parameters:")
    print(f"  ρ₀ = {RHO_0}, c_s = {C_S}, G = {G}")
    print(f"  λ_J = c_s√(π/Gρ₀) = {LAMBDA_J:.3f}")
    
    print(f"\nTest cases:")
    print(f"  Long wavelength (λ=4.0 > λ_J): Gravitationally UNSTABLE")
    print(f"  Short wavelength (λ=1.0 < λ_J): Pressure-supported STABLE")
    
    if 'long_gradh' in results and 'long_nogradh' in results:
        r1 = results['long_gradh']
        r2 = results['long_nogradh']
        
        # Estimate growth rates by fitting exponential
        if len(r1['amplitudes']) > 10:
            # Use early-time data where perturbation is still small
            mask1 = r1['amplitudes'] < 0.5
            if np.sum(mask1) > 5:
                t1 = r1['times'][mask1]
                a1 = r1['amplitudes'][mask1]
                log_a1 = np.log(a1)
                gamma1, _ = np.polyfit(t1, log_a1, 1)
                
                t2 = r2['times'][mask1[:len(r2['times'])]]
                a2 = r2['amplitudes'][:len(t2)]
                log_a2 = np.log(a2)
                gamma2, _ = np.polyfit(t2, log_a2, 1)
                
                # Theoretical growth rate
                k_long = 2 * np.pi / 4.0
                gamma_theory = np.sqrt(4 * np.pi * G * RHO_0 - C_S**2 * k_long**2)
                
                print(f"\nGrowth rates for long wavelength (λ=4.0):")
                print(f"  Theory:       γ = {gamma_theory:.4f}")
                print(f"  With grad-h:  γ = {gamma1:.4f} (error: {abs(gamma1-gamma_theory)/gamma_theory*100:.1f}%)")
                print(f"  No grad-h:    γ = {gamma2:.4f} (error: {abs(gamma2-gamma_theory)/gamma_theory*100:.1f}%)")
    
    if 'short_gradh' in results:
        print(f"\nShort wavelength (λ=1.0) behavior:")
        print(f"  Shows oscillation (pressure-supported), grad-h has minimal effect")
        print(f"  Theoretical frequency: ω = {np.sqrt(C_S**2 * (2*np.pi)**2 - 4*np.pi*G*RHO_0):.4f}")
    
    print("\n" + "="*60)


if __name__ == "__main__":
    main()
