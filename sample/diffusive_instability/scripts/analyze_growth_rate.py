#!/usr/bin/env python3
"""
Analyze diffusive instability growth rate from simulation data.

Measures the growth rate Γ from density perturbation evolution and compares
with theoretical prediction:
    Γ_theory = ε·c_s·h / R²

where ε ≈ 0.4 without grad-h correction, ε ≈ 0 with grad-h.

Usage:
    python analyze_growth_rate.py results/<run_name> --plot
"""

import argparse
import glob
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import curve_fit

# Add parent scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "scripts"))


def load_snapshot(filepath: str) -> dict:
    """Load particle data from HDF5 snapshot."""
    with h5py.File(filepath, 'r') as f:
        data = {
            'time': f['Header'].attrs['Time'][0] if 'Time' in f['Header'].attrs else 0.0,
            'x': f['PartType0/Coordinates'][:, 0],
            'rho': f['PartType0/Density'][:],
            'm': f['PartType0/Mass'][:],
        }
        # Optional fields
        if 'SmoothingLength' in f['PartType0']:
            data['h'] = f['PartType0/SmoothingLength'][:]
        if 'Velocities' in f['PartType0']:
            data['v'] = f['PartType0/Velocities'][:, 0]
    return data


def load_metadata(results_dir: str) -> dict:
    """Load simulation metadata from JSON output."""
    metadata_file = os.path.join(results_dir, 'output_metadata.json')
    if os.path.exists(metadata_file):
        with open(metadata_file, 'r') as f:
            return json.load(f)
    return {}


def compute_density_perturbation(data: dict, rho_analytic_func=None) -> float:
    """
    Compute density perturbation amplitude δρ/ρ.
    
    For isothermal slab: ρ(x) = ρ₀·sech²(x/H)
    We measure δρ = ρ_sim - ρ_analytic and return RMS(δρ)/mean(ρ).
    """
    x = data['x']
    rho = data['rho']
    
    if rho_analytic_func is not None:
        rho_analytic = rho_analytic_func(x)
        delta_rho = rho - rho_analytic
        return np.sqrt(np.mean(delta_rho**2)) / np.mean(rho)
    else:
        # Without analytic, use deviation from local mean
        # Sort by x and compute local mean
        idx = np.argsort(x)
        x_sorted = x[idx]
        rho_sorted = rho[idx]
        
        # Use kernel smoothing or simple windowed mean
        window = max(len(rho) // 20, 5)
        rho_smoothed = np.convolve(rho_sorted, np.ones(window)/window, mode='same')
        
        delta_rho = rho_sorted - rho_smoothed
        return np.sqrt(np.mean(delta_rho**2)) / np.mean(rho)


def compute_power_spectrum(data: dict, k_bins: int = 50) -> tuple:
    """
    Compute 1D power spectrum of density field.
    
    Returns (k, P(k)) where k is wavenumber and P(k) is power.
    """
    x = data['x']
    rho = data['rho']
    m = data['m']
    
    # Sort by position
    idx = np.argsort(x)
    x = x[idx]
    rho = rho[idx]
    m = m[idx]
    
    # Grid the density onto uniform mesh
    L = x.max() - x.min()
    N_grid = min(len(x) * 2, 1024)  # Grid resolution
    
    x_grid = np.linspace(x.min(), x.max(), N_grid)
    dx = x_grid[1] - x_grid[0]
    
    # Assign to grid using CIC (cloud-in-cell)
    rho_grid = np.zeros(N_grid)
    for i in range(len(x)):
        xi = (x[i] - x.min()) / dx
        j = int(xi)
        w = xi - j
        if j < N_grid - 1:
            rho_grid[j] += (1 - w) * m[i]
            rho_grid[j + 1] += w * m[i]
        else:
            rho_grid[j] += m[i]
    
    rho_grid /= dx  # Convert to density
    
    # FFT
    rho_fft = np.fft.rfft(rho_grid - np.mean(rho_grid))
    k = np.fft.rfftfreq(N_grid, dx) * 2 * np.pi
    Pk = np.abs(rho_fft)**2 / N_grid
    
    return k[1:], Pk[1:]  # Skip k=0


def fit_exponential_growth(times: np.ndarray, perturbations: np.ndarray,
                           fit_range: tuple = None) -> tuple:
    """
    Fit δρ/ρ = A·exp(Γ·t) to extract growth rate Γ.
    
    Returns (Gamma, A, fit_quality).
    """
    if fit_range is not None:
        mask = (times >= fit_range[0]) & (times <= fit_range[1])
        times = times[mask]
        perturbations = perturbations[mask]
    
    if len(times) < 3:
        return np.nan, np.nan, np.nan
    
    # Take log for linear fit
    log_pert = np.log(perturbations)
    
    # Filter out any invalid values
    valid = np.isfinite(log_pert)
    if np.sum(valid) < 3:
        return np.nan, np.nan, np.nan
    
    times = times[valid]
    log_pert = log_pert[valid]
    
    # Linear regression: log(δρ/ρ) = log(A) + Γ·t
    coeffs = np.polyfit(times, log_pert, 1)
    Gamma = coeffs[0]
    A = np.exp(coeffs[1])
    
    # Fit quality (R²)
    residuals = log_pert - (coeffs[0] * times + coeffs[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((log_pert - np.mean(log_pert))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    return Gamma, A, r_squared


def theoretical_growth_rate(c_s: float, h: float, R: float, 
                            epsilon: float = 0.4) -> float:
    """
    Compute theoretical diffusive instability growth rate.
    
    Γ = ε·c_s·h / R²
    
    Parameters:
        c_s: Sound speed
        h: Smoothing length
        R: Perturbation wavelength / 2π
        epsilon: Force error coefficient (≈0.4 without grad-h, ≈0 with)
    
    Returns:
        Γ: Growth rate [1/time]
    """
    return epsilon * c_s * h / R**2


def analyze_run(results_dir: str, plot: bool = False, 
                output_file: str = None) -> dict:
    """
    Analyze a single simulation run.
    
    Returns dictionary with growth rate and comparison to theory.
    """
    # Find snapshots
    snapshot_pattern = os.path.join(results_dir, 'snapshot_*.hdf5')
    snapshots = sorted(glob.glob(snapshot_pattern))
    
    if len(snapshots) < 2:
        print(f"Error: Need at least 2 snapshots, found {len(snapshots)}")
        return None
    
    print(f"Found {len(snapshots)} snapshots in {results_dir}")
    
    # Load metadata
    metadata = load_metadata(results_dir)
    
    # Load all snapshots and compute perturbation amplitude
    times = []
    perturbations = []
    mean_h = []
    mean_rho = []
    
    for snap_file in snapshots:
        data = load_snapshot(snap_file)
        times.append(data['time'])
        perturbations.append(compute_density_perturbation(data))
        if 'h' in data:
            mean_h.append(np.mean(data['h']))
        mean_rho.append(np.mean(data['rho']))
    
    times = np.array(times)
    perturbations = np.array(perturbations)
    
    # Fit growth rate
    # Use latter half of simulation for better exponential regime
    t_mid = 0.5 * (times.min() + times.max())
    Gamma_measured, A_fit, r_squared = fit_exponential_growth(
        times, perturbations, fit_range=(t_mid, times.max())
    )
    
    # Extract parameters from config if available
    config_file = os.path.join(results_dir, '..', 'config.json')
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
    else:
        config = {}
    
    c_s = config.get('sound_speed', metadata.get('sound_speed', 1.0))
    H = config.get('scale_height', metadata.get('scale_height', 0.5))
    use_grad_h = config.get('use_grad_h_correction', 
                            metadata.get('use_grad_h_correction', False))
    
    h_mean = np.mean(mean_h) if mean_h else H / 10  # Estimate
    
    # Theoretical prediction
    # Use scale height H as characteristic wavelength R
    epsilon = 0.0 if use_grad_h else 0.4
    Gamma_theory = theoretical_growth_rate(c_s, h_mean, H, epsilon)
    
    results = {
        'Gamma_measured': Gamma_measured,
        'Gamma_theory': Gamma_theory,
        'epsilon_effective': Gamma_measured * H**2 / (c_s * h_mean) if h_mean > 0 else np.nan,
        'epsilon_expected': epsilon,
        'fit_amplitude': A_fit,
        'fit_r_squared': r_squared,
        'sound_speed': c_s,
        'scale_height': H,
        'mean_h': h_mean,
        'use_grad_h': use_grad_h,
        'n_snapshots': len(snapshots),
        'time_range': [times.min(), times.max()],
    }
    
    print("\n" + "="*60)
    print("DIFFUSIVE INSTABILITY ANALYSIS RESULTS")
    print("="*60)
    print(f"Measured growth rate Γ = {Gamma_measured:.4f}")
    print(f"Theoretical prediction Γ = {Gamma_theory:.4f} (ε = {epsilon})")
    print(f"Effective ε from measurement = {results['epsilon_effective']:.4f}")
    print(f"Fit R² = {r_squared:.4f}")
    print(f"Grad-h correction: {'ON' if use_grad_h else 'OFF'}")
    print("="*60)
    
    if plot:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Perturbation growth
        ax1 = axes[0, 0]
        ax1.semilogy(times, perturbations, 'ko-', label='Measured')
        
        # Plot fit
        t_fit = np.linspace(times.min(), times.max(), 100)
        ax1.semilogy(t_fit, A_fit * np.exp(Gamma_measured * t_fit), 
                     'r--', label=f'Fit: Γ = {Gamma_measured:.4f}')
        
        # Plot theory
        if Gamma_theory > 0:
            pert_init = perturbations[0]
            ax1.semilogy(t_fit, pert_init * np.exp(Gamma_theory * t_fit), 
                         'b:', label=f'Theory: Γ = {Gamma_theory:.4f}')
        
        ax1.set_xlabel('Time')
        ax1.set_ylabel('δρ/ρ')
        ax1.set_title('Density Perturbation Growth')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Final density profile
        ax2 = axes[0, 1]
        final_data = load_snapshot(snapshots[-1])
        idx = np.argsort(final_data['x'])
        ax2.plot(final_data['x'][idx], final_data['rho'][idx], 'b-', alpha=0.7)
        ax2.set_xlabel('x')
        ax2.set_ylabel('ρ')
        ax2.set_title(f'Final Density Profile (t = {final_data["time"]:.2f})')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Power spectrum evolution
        ax3 = axes[1, 0]
        colors = plt.cm.viridis(np.linspace(0, 1, min(5, len(snapshots))))
        snap_indices = np.linspace(0, len(snapshots)-1, min(5, len(snapshots)), dtype=int)
        
        for i, idx in enumerate(snap_indices):
            data = load_snapshot(snapshots[idx])
            k, Pk = compute_power_spectrum(data)
            ax3.loglog(k, Pk, color=colors[i], alpha=0.7,
                      label=f't = {data["time"]:.2f}')
        
        ax3.set_xlabel('k')
        ax3.set_ylabel('P(k)')
        ax3.set_title('Power Spectrum Evolution')
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Smoothing length evolution
        ax4 = axes[1, 1]
        if mean_h:
            ax4.plot(times, mean_h, 'g-o')
            ax4.set_xlabel('Time')
            ax4.set_ylabel('Mean h')
            ax4.set_title('Mean Smoothing Length')
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(0.5, 0.5, 'No smoothing length data', 
                    ha='center', va='center', transform=ax4.transAxes)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {output_file}")
        else:
            plt.show()
    
    return results


def compare_runs(run_dirs: list, labels: list = None, 
                 output_file: str = None) -> dict:
    """
    Compare growth rates across multiple runs.
    """
    import matplotlib.pyplot as plt
    
    if labels is None:
        labels = [os.path.basename(d) for d in run_dirs]
    
    all_results = []
    for run_dir, label in zip(run_dirs, labels):
        print(f"\nAnalyzing: {label}")
        results = analyze_run(run_dir, plot=False)
        if results:
            results['label'] = label
            all_results.append(results)
    
    if not all_results:
        print("No valid results to compare")
        return None
    
    # Create comparison plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Bar chart of growth rates
    ax1 = axes[0]
    x = np.arange(len(all_results))
    width = 0.35
    
    Gamma_measured = [r['Gamma_measured'] for r in all_results]
    Gamma_theory = [r['Gamma_theory'] for r in all_results]
    
    bars1 = ax1.bar(x - width/2, Gamma_measured, width, label='Measured', color='steelblue')
    bars2 = ax1.bar(x + width/2, Gamma_theory, width, label='Theory', color='coral')
    
    ax1.set_ylabel('Growth Rate Γ')
    ax1.set_title('Measured vs Theoretical Growth Rates')
    ax1.set_xticks(x)
    ax1.set_xticklabels([r['label'] for r in all_results], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Effective epsilon
    ax2 = axes[1]
    epsilon_eff = [r['epsilon_effective'] for r in all_results]
    colors = ['green' if r['use_grad_h'] else 'red' for r in all_results]
    
    bars = ax2.bar(x, epsilon_eff, color=colors, alpha=0.7)
    ax2.axhline(y=0.4, color='r', linestyle='--', label='ε = 0.4 (no grad-h)')
    ax2.axhline(y=0.0, color='g', linestyle='--', label='ε = 0 (with grad-h)')
    
    ax2.set_ylabel('Effective ε')
    ax2.set_title('Effective Force Error Coefficient')
    ax2.set_xticks(x)
    ax2.set_xticklabels([r['label'] for r in all_results], rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Comparison plot saved to: {output_file}")
    else:
        plt.show()
    
    return {'runs': all_results}


def main():
    parser = argparse.ArgumentParser(
        description='Analyze diffusive instability growth rate from SPH simulations'
    )
    parser.add_argument('results_dir', nargs='+',
                        help='Results directory (or multiple for comparison)')
    parser.add_argument('--plot', '-p', action='store_true',
                        help='Generate plots')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output file for plot (PNG/PDF)')
    parser.add_argument('--compare', '-c', action='store_true',
                        help='Compare multiple runs')
    parser.add_argument('--labels', '-l', nargs='+', default=None,
                        help='Labels for comparison plot')
    
    args = parser.parse_args()
    
    if args.compare or len(args.results_dir) > 1:
        compare_runs(args.results_dir, args.labels, args.output)
    else:
        analyze_run(args.results_dir[0], args.plot, args.output)


if __name__ == '__main__':
    main()
