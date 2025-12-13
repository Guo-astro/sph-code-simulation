#!/usr/bin/env python3
"""
Animate 1D isothermal slab evolution and compare theory vs simulation.

Creates side-by-side animations showing:
1. Density profile evolution
2. Central density vs time (with theoretical prediction)
3. Phase space (x vs v)
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import glob
import os
import sys
from pathlib import Path

# Parameters
plt.rcParams['figure.dpi'] = 100
plt.rcParams['animation.html'] = 'jshtml'


def load_csv_snapshot(filepath):
    """Load particle data from CSV snapshot."""
    df = pd.read_csv(filepath, comment='#')
    return {
        'x': df['pos_x'].values,
        'rho': df['dens'].values,
        'vel': df['vel_x'].values,
        'mass': df['mass'].values,
        'h': df['sml'].values,
        'pres': df['pres'].values,
        'grav_acc': df['grav_acc_x'].values if 'grav_acc_x' in df.columns else np.zeros(len(df)),
        'phi': df['phi'].values if 'phi' in df.columns else np.zeros(len(df)),
    }


def load_energy_file(filepath):
    """Load energy evolution data."""
    data = np.loadtxt(filepath, comments='#')
    return {
        'time': data[:, 0],
        'kinetic': data[:, 1],
        'thermal': data[:, 2],
        'potential': data[:, 3],
        'total': data[:, 4],
    }


def analytic_sech2_profile(x, rho_center, H):
    """Analytic density profile: ρ(x) = ρ₀ sech²(x/H)"""
    return rho_center / np.cosh(x / H)**2


def theoretical_growth_rate(c_s, h_mean, L, epsilon=0.35):
    """
    Theoretical diffusive instability growth rate.
    
    Γ = ε·c_s·h / L²
    
    For cubic spline: ε ≈ 0.35
    For Wendland C4: ε ≈ 0.4
    
    L is the characteristic wavelength of the perturbation.
    """
    return epsilon * c_s * h_mean / L**2


def create_animation(gradh_dir, no_gradh_dir, output_file='slab_evolution.mp4',
                     fps=10, params=None):
    """
    Create side-by-side animation comparing grad-h vs no-grad-h runs.
    """
    if params is None:
        params = {
            'rho_center': 1.0,
            'H': 0.5,  # Scale height
            'c_s': 1.0,
            'G': 1.0,
            'gamma': 1.4,
        }
    
    # Load all snapshots
    gradh_snaps = sorted(glob.glob(os.path.join(gradh_dir, 'snapshot_*.csv')))
    no_gradh_snaps = sorted(glob.glob(os.path.join(no_gradh_dir, 'snapshot_*.csv')))
    
    n_frames = min(len(gradh_snaps), len(no_gradh_snaps))
    if n_frames < 2:
        print(f"Error: Need at least 2 snapshots. Found {len(gradh_snaps)} grad-h, {len(no_gradh_snaps)} no-grad-h")
        return
    
    print(f"Creating animation with {n_frames} frames...")
    
    # Load all data
    gradh_data = [load_csv_snapshot(f) for f in gradh_snaps[:n_frames]]
    no_gradh_data = [load_csv_snapshot(f) for f in no_gradh_snaps[:n_frames]]
    
    # Load energy files
    gradh_energy = load_energy_file(os.path.join(gradh_dir, 'energy.dat'))
    no_gradh_energy = load_energy_file(os.path.join(no_gradh_dir, 'energy.dat'))
    
    # Compute central density evolution
    gradh_central = []
    no_gradh_central = []
    times = []
    
    for i, (gd, ngd) in enumerate(zip(gradh_data, no_gradh_data)):
        times.append(i)  # Snapshot index as time proxy
        
        # Central density (mean of particles near x=0)
        g_center_mask = np.abs(gd['x']) < 0.1
        ng_center_mask = np.abs(ngd['x']) < 0.1
        
        if np.sum(g_center_mask) > 0:
            gradh_central.append(np.mean(gd['rho'][g_center_mask]))
        else:
            gradh_central.append(np.max(gd['rho']))
            
        if np.sum(ng_center_mask) > 0:
            no_gradh_central.append(np.mean(ngd['rho'][ng_center_mask]))
        else:
            no_gradh_central.append(np.max(ngd['rho']))
    
    gradh_central = np.array(gradh_central)
    no_gradh_central = np.array(no_gradh_central)
    times = np.array(times)
    
    # Compute mean smoothing length for theoretical prediction
    h_mean = np.mean([np.mean(d['h']) for d in gradh_data])
    
    # Create figure
    fig = plt.figure(figsize=(16, 10))
    
    # Layout: 2 rows, 3 columns
    # Row 1: Density profiles (grad-h, no-grad-h, comparison)
    # Row 2: Phase space, Central density, Energy
    ax1 = fig.add_subplot(2, 3, 1)  # Grad-h density
    ax2 = fig.add_subplot(2, 3, 2)  # No-grad-h density
    ax3 = fig.add_subplot(2, 3, 3)  # Overlay comparison
    ax4 = fig.add_subplot(2, 3, 4)  # Phase space grad-h
    ax5 = fig.add_subplot(2, 3, 5)  # Phase space no-grad-h
    ax6 = fig.add_subplot(2, 3, 6)  # Central density evolution
    
    # X range for analytic profile
    x_analytic = np.linspace(-2, 2, 200)
    rho_analytic = analytic_sech2_profile(x_analytic, params['rho_center'], params['H'])
    
    # Initialize plots
    scatter1, = ax1.plot([], [], 'b.', markersize=3, alpha=0.7)
    line1_analytic, = ax1.plot(x_analytic, rho_analytic, 'k--', linewidth=1, alpha=0.5, label='Analytic')
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(0, 1.5)
    ax1.set_xlabel('x')
    ax1.set_ylabel('ρ')
    ax1.set_title('With grad-h correction')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    scatter2, = ax2.plot([], [], 'r.', markersize=3, alpha=0.7)
    line2_analytic, = ax2.plot(x_analytic, rho_analytic, 'k--', linewidth=1, alpha=0.5, label='Analytic')
    ax2.set_xlim(-2, 2)
    ax2.set_ylim(0, 1.5)
    ax2.set_xlabel('x')
    ax2.set_ylabel('ρ')
    ax2.set_title('Without grad-h correction')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    scatter3_gradh, = ax3.plot([], [], 'b.', markersize=3, alpha=0.5, label='With grad-h')
    scatter3_no_gradh, = ax3.plot([], [], 'r.', markersize=3, alpha=0.5, label='Without grad-h')
    line3_analytic, = ax3.plot(x_analytic, rho_analytic, 'k--', linewidth=1, alpha=0.5, label='Analytic')
    ax3.set_xlim(-2, 2)
    ax3.set_ylim(0, 1.5)
    ax3.set_xlabel('x')
    ax3.set_ylabel('ρ')
    ax3.set_title('Comparison')
    ax3.legend(loc='upper right', fontsize=8)
    ax3.grid(True, alpha=0.3)
    
    # Phase space
    scatter4, = ax4.plot([], [], 'b.', markersize=3, alpha=0.7)
    ax4.set_xlim(-2, 2)
    ax4.set_ylim(-0.5, 0.5)
    ax4.set_xlabel('x')
    ax4.set_ylabel('v')
    ax4.set_title('Phase space (grad-h)')
    ax4.grid(True, alpha=0.3)
    
    scatter5, = ax5.plot([], [], 'r.', markersize=3, alpha=0.7)
    ax5.set_xlim(-2, 2)
    ax5.set_ylim(-0.5, 0.5)
    ax5.set_xlabel('x')
    ax5.set_ylabel('v')
    ax5.set_title('Phase space (no grad-h)')
    ax5.grid(True, alpha=0.3)
    
    # Central density evolution
    line6_gradh, = ax6.plot([], [], 'b-', linewidth=2, label='With grad-h')
    line6_no_gradh, = ax6.plot([], [], 'r-', linewidth=2, label='Without grad-h')
    marker6, = ax6.plot([], [], 'ko', markersize=8)
    ax6.set_xlim(0, n_frames)
    ax6.set_ylim(0, max(np.max(gradh_central), np.max(no_gradh_central)) * 1.2)
    ax6.set_xlabel('Snapshot (∝ time)')
    ax6.set_ylabel('Central ρ')
    ax6.set_title('Central Density Evolution')
    ax6.legend(loc='upper left', fontsize=8)
    ax6.grid(True, alpha=0.3)
    
    # Time annotation
    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.suptitle('1D Isothermal Slab: Diffusive Instability Study', fontsize=14)
    
    def init():
        scatter1.set_data([], [])
        scatter2.set_data([], [])
        scatter3_gradh.set_data([], [])
        scatter3_no_gradh.set_data([], [])
        scatter4.set_data([], [])
        scatter5.set_data([], [])
        line6_gradh.set_data([], [])
        line6_no_gradh.set_data([], [])
        marker6.set_data([], [])
        time_text.set_text('')
        return (scatter1, scatter2, scatter3_gradh, scatter3_no_gradh, 
                scatter4, scatter5, line6_gradh, line6_no_gradh, marker6, time_text)
    
    def animate(frame):
        gd = gradh_data[frame]
        ngd = no_gradh_data[frame]
        
        # Sort by x for line plots
        g_idx = np.argsort(gd['x'])
        ng_idx = np.argsort(ngd['x'])
        
        scatter1.set_data(gd['x'][g_idx], gd['rho'][g_idx])
        scatter2.set_data(ngd['x'][ng_idx], ngd['rho'][ng_idx])
        scatter3_gradh.set_data(gd['x'][g_idx], gd['rho'][g_idx])
        scatter3_no_gradh.set_data(ngd['x'][ng_idx], ngd['rho'][ng_idx])
        
        scatter4.set_data(gd['x'], gd['vel'])
        scatter5.set_data(ngd['x'], ngd['vel'])
        
        line6_gradh.set_data(times[:frame+1], gradh_central[:frame+1])
        line6_no_gradh.set_data(times[:frame+1], no_gradh_central[:frame+1])
        marker6.set_data([times[frame]], [gradh_central[frame]])
        
        time_text.set_text(f'Snapshot: {frame} / {n_frames-1}')
        
        return (scatter1, scatter2, scatter3_gradh, scatter3_no_gradh,
                scatter4, scatter5, line6_gradh, line6_no_gradh, marker6, time_text)
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, 
                                   frames=n_frames, interval=1000//fps, blit=True)
    
    # Save animation
    if output_file.endswith('.gif'):
        writer = animation.PillowWriter(fps=fps)
    else:
        writer = animation.FFMpegWriter(fps=fps, bitrate=2000)
    
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer=writer)
    print(f"Animation saved!")
    
    plt.close()
    return anim


def analyze_theory_vs_simulation(gradh_dir, no_gradh_dir, output_file=None):
    """
    Detailed comparison of theory vs simulation.
    """
    # Load snapshots
    gradh_snaps = sorted(glob.glob(os.path.join(gradh_dir, 'snapshot_*.csv')))
    no_gradh_snaps = sorted(glob.glob(os.path.join(no_gradh_dir, 'snapshot_*.csv')))
    
    n_frames = min(len(gradh_snaps), len(no_gradh_snaps))
    
    # Extract parameters from first snapshot
    d0 = load_csv_snapshot(gradh_snaps[0])
    N = len(d0['x'])
    h_mean = np.mean(d0['h'])
    
    # Simulation parameters (from config)
    params = {
        'N': N,
        'rho_center': 1.0,
        'H': 0.5,
        'c_s': 1.0,
        'G': 1.0,
        'gamma': 1.4,
        'epsilon_cubic': 0.35,  # Cubic spline
    }
    
    print("="*70)
    print("THEORY VS SIMULATION ANALYSIS")
    print("="*70)
    print(f"\nSimulation Parameters:")
    print(f"  N = {N} particles")
    print(f"  Mean smoothing length h = {h_mean:.4f}")
    print(f"  Scale height H = {params['H']}")
    print(f"  Sound speed c_s = {params['c_s']}")
    print(f"  ε (cubic spline) = {params['epsilon_cubic']}")
    
    # Theoretical predictions
    # For isothermal slab, characteristic wavelength ~ H
    L = params['H']
    Gamma_theory = theoretical_growth_rate(params['c_s'], h_mean, L, params['epsilon_cubic'])
    
    print(f"\nTheoretical Predictions:")
    print(f"  Diffusion coefficient D = ε·c_s·h = {params['epsilon_cubic'] * params['c_s'] * h_mean:.4f}")
    print(f"  Growth rate Γ = D/L² = {Gamma_theory:.6f}")
    print(f"  e-folding time = 1/Γ = {1/Gamma_theory:.2f}")
    
    # Measure from simulation
    gradh_central = []
    no_gradh_central = []
    
    for gf, ngf in zip(gradh_snaps[:n_frames], no_gradh_snaps[:n_frames]):
        gd = load_csv_snapshot(gf)
        ngd = load_csv_snapshot(ngf)
        
        g_mask = np.abs(gd['x']) < 0.1
        ng_mask = np.abs(ngd['x']) < 0.1
        
        gradh_central.append(np.mean(gd['rho'][g_mask]) if np.sum(g_mask) > 0 else np.max(gd['rho']))
        no_gradh_central.append(np.mean(ngd['rho'][ng_mask]) if np.sum(ng_mask) > 0 else np.max(ngd['rho']))
    
    gradh_central = np.array(gradh_central)
    no_gradh_central = np.array(no_gradh_central)
    times = np.arange(n_frames)  # Snapshot index as time
    
    # Fit exponential to no-grad-h data
    # log(ρ) = log(ρ₀) + Γt
    # Use linear fit on log data
    valid_ng = no_gradh_central > 0
    if np.sum(valid_ng) > 5:
        log_rho = np.log(no_gradh_central[valid_ng])
        t_valid = times[valid_ng]
        
        # Linear regression
        coeffs = np.polyfit(t_valid, log_rho, 1)
        Gamma_measured = coeffs[0]  # This is Γ per snapshot
        rho0_fit = np.exp(coeffs[1])
        
        print(f"\nMeasured from Simulation (no grad-h):")
        print(f"  Fitted Γ = {Gamma_measured:.6f} per snapshot")
        print(f"  Fitted ρ₀ = {rho0_fit:.4f}")
        
        # Compare
        print(f"\nComparison:")
        print(f"  Γ_theory / Γ_measured = {Gamma_theory / Gamma_measured:.4f}")
    
    # Analyze grad-h stability
    valid_g = gradh_central > 0
    if np.sum(valid_g) > 5:
        log_rho_g = np.log(gradh_central[valid_g])
        t_valid_g = times[valid_g]
        coeffs_g = np.polyfit(t_valid_g, log_rho_g, 1)
        Gamma_gradh = coeffs_g[0]
        
        print(f"\nMeasured from Simulation (with grad-h):")
        print(f"  Fitted Γ = {Gamma_gradh:.6f} per snapshot")
        print(f"  Expected: Γ ≈ 0 (stable)")
    
    # Create analysis figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: Central density vs time
    ax1 = axes[0, 0]
    ax1.semilogy(times, gradh_central, 'b-o', label='With grad-h', markersize=3)
    ax1.semilogy(times, no_gradh_central, 'r-o', label='Without grad-h', markersize=3)
    
    # Theoretical prediction for no-grad-h
    t_theory = np.linspace(0, n_frames, 100)
    rho_theory = no_gradh_central[0] * np.exp(Gamma_theory * t_theory)
    ax1.semilogy(t_theory, rho_theory, 'r--', alpha=0.5, label=f'Theory: Γ={Gamma_theory:.4f}')
    
    ax1.set_xlabel('Snapshot')
    ax1.set_ylabel('Central ρ (log scale)')
    ax1.set_title('Central Density Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Density profiles at t=0, mid, final
    ax2 = axes[0, 1]
    snap_indices = [0, n_frames//2, n_frames-1]
    colors = ['green', 'orange', 'purple']
    for idx, c in zip(snap_indices, colors):
        d = load_csv_snapshot(no_gradh_snaps[idx])
        sort_idx = np.argsort(d['x'])
        ax2.plot(d['x'][sort_idx], d['rho'][sort_idx], '-', color=c, 
                 alpha=0.7, label=f't={idx}')
    
    x_analytic = np.linspace(-2, 2, 200)
    rho_analytic = analytic_sech2_profile(x_analytic, params['rho_center'], params['H'])
    ax2.plot(x_analytic, rho_analytic, 'k--', label='Analytic IC')
    ax2.set_xlabel('x')
    ax2.set_ylabel('ρ')
    ax2.set_title('Density Profile (no grad-h)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Same for grad-h
    ax3 = axes[0, 2]
    for idx, c in zip(snap_indices, colors):
        d = load_csv_snapshot(gradh_snaps[idx])
        sort_idx = np.argsort(d['x'])
        ax3.plot(d['x'][sort_idx], d['rho'][sort_idx], '-', color=c,
                 alpha=0.7, label=f't={idx}')
    ax3.plot(x_analytic, rho_analytic, 'k--', label='Analytic IC')
    ax3.set_xlabel('x')
    ax3.set_ylabel('ρ')
    ax3.set_title('Density Profile (with grad-h)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: h distribution
    ax4 = axes[1, 0]
    d_final = load_csv_snapshot(no_gradh_snaps[-1])
    ax4.scatter(d_final['x'], d_final['h'], s=5, alpha=0.7)
    ax4.axhline(h_mean, color='r', linestyle='--', label=f'Mean h={h_mean:.3f}')
    ax4.set_xlabel('x')
    ax4.set_ylabel('Smoothing length h')
    ax4.set_title('Smoothing Length Distribution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Gravity check
    ax5 = axes[1, 1]
    d_init = load_csv_snapshot(gradh_snaps[0])
    sort_idx = np.argsort(d_init['x'])
    ax5.plot(d_init['x'][sort_idx], d_init['grav_acc'][sort_idx], 'b-', label='Simulation')
    
    # Analytic gravity for sech² profile: g = -2πG sign(x) Σ(|x|)
    # where Σ(x) = ∫₀ˣ ρ dx = ρ₀ H tanh(x/H)
    x_th = np.linspace(-2, 2, 200)
    Sigma = params['rho_center'] * params['H'] * np.tanh(np.abs(x_th) / params['H'])
    g_th = -2 * np.pi * params['G'] * np.sign(x_th) * Sigma
    ax5.plot(x_th, g_th, 'r--', label='Analytic')
    ax5.set_xlabel('x')
    ax5.set_ylabel('g(x)')
    ax5.set_title('Gravitational Acceleration')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Energy conservation
    ax6 = axes[1, 2]
    gradh_energy = load_energy_file(os.path.join(gradh_dir, 'energy.dat'))
    no_gradh_energy = load_energy_file(os.path.join(no_gradh_dir, 'energy.dat'))
    
    ax6.plot(gradh_energy['time'], gradh_energy['total'], 'b-', label='With grad-h')
    ax6.plot(no_gradh_energy['time'], no_gradh_energy['total'], 'r-', label='Without grad-h')
    ax6.set_xlabel('Time')
    ax6.set_ylabel('Total Energy')
    ax6.set_title('Energy Conservation')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"\nAnalysis figure saved to: {output_file}")
    else:
        plt.show()
    
    plt.close()
    
    # Return analysis results
    return {
        'Gamma_theory': Gamma_theory,
        'Gamma_measured_no_gradh': Gamma_measured if np.sum(valid_ng) > 5 else None,
        'Gamma_measured_gradh': Gamma_gradh if np.sum(valid_g) > 5 else None,
        'h_mean': h_mean,
        'N': N,
    }


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Animate and analyze isothermal slab simulations')
    parser.add_argument('--gradh-dir', default='simulations/stability/diffusive_instability/results/slab_gradh',
                        help='Directory with grad-h results')
    parser.add_argument('--no-gradh-dir', default='simulations/stability/diffusive_instability/results/slab_no_gradh',
                        help='Directory with no-grad-h results')
    parser.add_argument('--output', '-o', default='simulations/stability/diffusive_instability/results/slab_evolution.gif',
                        help='Output animation file')
    parser.add_argument('--analysis', '-a', default=None,
                        help='Output file for analysis figure')
    parser.add_argument('--fps', type=int, default=5, help='Frames per second')
    parser.add_argument('--analyze-only', action='store_true', help='Only run analysis, skip animation')
    
    args = parser.parse_args()
    
    # Run analysis
    results = analyze_theory_vs_simulation(args.gradh_dir, args.no_gradh_dir,
                                           args.analysis or 'simulations/stability/diffusive_instability/results/theory_vs_simulation.png')
    
    # Create animation
    if not args.analyze_only:
        create_animation(args.gradh_dir, args.no_gradh_dir, args.output, args.fps)
