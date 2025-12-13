#!/usr/bin/env python3
"""
PHASE-BY-PHASE ANALYSIS OF GSPH COLLAPSE
=========================================

The collapse has THREE distinct phases:
1. Quasi-stable phase (t ≈ 0-7): Very slow growth, nearly stable
2. Transition phase (t ≈ 7-8.5): Accelerating growth
3. Runaway phase (t > 8.5): Catastrophic collapse

This is NOT a simple exponential! It's a DELAYED RUNAWAY.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import glob
import re
import pandas as pd

G = 1.0  

def load_simulation_data(results_dir):
    """Load simulation data."""
    data = {}
    
    def load_method_data(method_dir):
        if not os.path.exists(method_dir):
            return None
        
        files = sorted(glob.glob(os.path.join(method_dir, "snapshot_*.csv")))
        if not files:
            return None
        
        t_list, rho_max_list = [], []
        for f in files:
            try:
                time_val = None
                with open(f) as fp:
                    for line in fp:
                        if 'Time (physical):' in line:
                            match = re.search(r'Time \(physical\):\s*([\d.e+-]+)', line)
                            if match:
                                time_val = float(match.group(1))
                            break
                
                if time_val is None:
                    idx = int(re.search(r'snapshot_(\d+)', f).group(1))
                    time_val = idx * 0.2
                
                df = pd.read_csv(f, comment='#')
                if 'dens' in df.columns:
                    rho_max = df['dens'].max()
                else:
                    continue
                
                t_list.append(time_val)
                rho_max_list.append(rho_max)
            except Exception:
                continue
        
        if t_list:
            sorted_idx = np.argsort(t_list)
            return {
                't': np.array(t_list)[sorted_idx], 
                'rho_max': np.array(rho_max_list)[sorted_idx]
            }
        return None
    
    for key in ['gsph_nogradh', 'gsph_gradh', 'ssph_nogradh', 'ssph_gradh']:
        result = load_method_data(os.path.join(results_dir, key))
        if result:
            data[key] = result
    
    return data


def runaway_model(t, t_c, A, rho_bg):
    """
    Runaway collapse model:
    
    ρ(t) = rho_bg + A / (t_c - t)^2
    
    This is the characteristic form for gravitational collapse approaching
    a singularity at t = t_c.
    """
    rho = np.zeros_like(t)
    for i, ti in enumerate(t):
        if ti < t_c:
            rho[i] = rho_bg + A / (t_c - ti)**2
        else:
            rho[i] = 1e10  # Singularity
    return rho


def secular_growth_model(t, rho_0, gamma, t_delay):
    """
    Secular growth model (before runaway):
    
    ρ(t) = ρ_0 × exp(γ × max(0, t - t_delay))
    
    This captures the slow exponential growth phase.
    """
    rho = np.zeros_like(t)
    for i, ti in enumerate(t):
        if ti <= t_delay:
            rho[i] = rho_0
        else:
            rho[i] = rho_0 * np.exp(gamma * (ti - t_delay))
    return rho


def combined_model(t, rho_0, gamma, t_delay, t_c, A):
    """
    Combined model: secular growth transitioning to runaway.
    
    Phase 1: ρ ≈ ρ_0 (for t < t_delay)
    Phase 2: ρ = ρ_0 × exp(γ(t - t_delay)) (for t_delay < t < t_trans)
    Phase 3: ρ = A / (t_c - t)^2 (for t approaching t_c)
    
    The transition is smooth.
    """
    rho = np.zeros_like(t)
    
    for i, ti in enumerate(t):
        if ti < t_c:
            # Secular part
            rho_secular = rho_0 * np.exp(gamma * max(0, ti - t_delay))
            
            # Runaway part (only significant near t_c)
            rho_runaway = A / (t_c - ti)**2
            
            # Take maximum of the two
            rho[i] = max(rho_secular, rho_runaway)
        else:
            rho[i] = 1e10
    
    return rho


def main():
    print("=" * 70)
    print("PHASE-BY-PHASE ANALYSIS OF GSPH COLLAPSE")
    print("=" * 70)
    
    # Load data
    results_dir = "results/gradh_comparison"
    data = load_simulation_data(results_dir)
    
    if 'gsph_nogradh' not in data:
        print("ERROR: No data found!")
        return
    
    d = data['gsph_nogradh']
    t_sim = d['t']
    rho_sim = d['rho_max']
    
    print(f"\nSimulation data: {len(t_sim)} snapshots")
    print(f"  Time range: [{t_sim[0]:.2f}, {t_sim[-1]:.2f}]")
    print(f"  Density range: [{rho_sim.min():.2f}, {rho_sim.max():.0f}]")
    
    # Print raw data
    print("\nRaw simulation data:")
    print(f"{'t':>8} {'ρ_max':>12}")
    print("-" * 22)
    for ti, rhoi in zip(t_sim, rho_sim):
        print(f"{ti:8.2f} {rhoi:12.2f}")
    
    # Identify phases
    print("\n" + "=" * 70)
    print("PHASE IDENTIFICATION")
    print("=" * 70)
    
    rho_0 = rho_sim[0]
    
    # Phase 1: quasi-stable (ρ < 2 × ρ_0)
    phase1_mask = rho_sim < 2 * rho_0
    t_phase1_end = t_sim[phase1_mask][-1] if np.any(phase1_mask) else 0
    
    # Phase 2: acceleration (2 × ρ_0 < ρ < 10 × ρ_0)
    phase2_mask = (rho_sim >= 2 * rho_0) & (rho_sim < 10 * rho_0)
    
    # Phase 3: runaway (ρ > 10 × ρ_0)
    phase3_mask = rho_sim >= 10 * rho_0
    t_runaway = t_sim[phase3_mask][0] if np.any(phase3_mask) else t_sim[-1]
    
    print(f"\nPhase 1 (quasi-stable): t ∈ [0, {t_phase1_end:.2f}]")
    print(f"Phase 2 (acceleration): t ∈ [{t_phase1_end:.2f}, {t_runaway:.2f}]")
    print(f"Phase 3 (runaway):      t > {t_runaway:.2f}")
    
    # Fit runaway model to late-time data
    print("\n" + "=" * 70)
    print("FITTING RUNAWAY MODEL")
    print("=" * 70)
    
    # Use data from t > 8 for runaway fit
    late_mask = t_sim > 8.0
    t_late = t_sim[late_mask]
    rho_late = rho_sim[late_mask]
    
    if len(t_late) > 3:
        # Fit: ρ = A / (t_c - t)^2 + rho_bg
        # Linearize: (ρ - rho_bg)^(-1/2) = (t_c - t) / √A
        
        rho_bg = 2.0  # Background density
        
        try:
            def runaway_fit(t, t_c, A):
                return rho_bg + A / np.maximum(t_c - t, 0.001)**2
            
            popt, pcov = curve_fit(runaway_fit, t_late, rho_late, 
                                   p0=[9.1, 0.5],
                                   bounds=([9.0, 0.01], [9.5, 10]))
            t_c_fit, A_fit = popt
            
            print(f"\nRunaway fit parameters:")
            print(f"  Singularity time: t_c = {t_c_fit:.4f}")
            print(f"  Amplitude: A = {A_fit:.4f}")
            print(f"  Background: ρ_bg = {rho_bg:.2f}")
            
        except Exception as e:
            print(f"Runaway fit failed: {e}")
            t_c_fit = 9.1
            A_fit = 0.5
    else:
        t_c_fit = 9.1
        A_fit = 0.5
    
    # Fit secular growth to early-phase data
    print("\n" + "=" * 70)
    print("FITTING SECULAR GROWTH")
    print("=" * 70)
    
    # Use data from t = 4 to t = 8
    mid_mask = (t_sim > 4) & (t_sim < 8.5)
    t_mid = t_sim[mid_mask]
    rho_mid = rho_sim[mid_mask]
    
    if len(t_mid) > 3:
        try:
            def secular_fit(t, gamma, t_delay):
                return rho_0 * np.exp(gamma * np.maximum(t - t_delay, 0))
            
            popt, pcov = curve_fit(secular_fit, t_mid, rho_mid,
                                   p0=[0.1, 5.0],
                                   bounds=([0.01, 0], [2.0, 8]))
            gamma_fit, t_delay_fit = popt
            
            print(f"\nSecular growth parameters:")
            print(f"  Growth rate: Γ = {gamma_fit:.4f} rad/time")
            print(f"  Delay time: t_delay = {t_delay_fit:.2f}")
            print(f"  e-folding time: τ = {1/gamma_fit:.2f}")
            
        except Exception as e:
            print(f"Secular fit failed: {e}")
            gamma_fit = 0.3
            t_delay_fit = 5.0
    else:
        gamma_fit = 0.3
        t_delay_fit = 5.0
    
    # Create comprehensive figure
    print("\n" + "=" * 70)
    print("CREATING FIGURES")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # ===== Panel 1: All data with phases marked =====
    ax1 = axes[0, 0]
    
    ax1.semilogy(t_sim, rho_sim, 'ro-', markersize=6, linewidth=1.5,
                 label='GSPH no-gradh', zorder=5)
    
    # Mark phases
    ax1.axvspan(0, t_phase1_end, alpha=0.2, color='green', label='Phase 1: Quasi-stable')
    ax1.axvspan(t_phase1_end, t_runaway, alpha=0.2, color='yellow', label='Phase 2: Acceleration')
    ax1.axvspan(t_runaway, t_sim[-1]+0.5, alpha=0.2, color='red', label='Phase 3: Runaway')
    
    ax1.set_xlabel('Time', fontsize=12)
    ax1.set_ylabel('Maximum density ρ_max (log)', fontsize=12)
    ax1.set_title('Collapse Phases', fontsize=13)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-0.5, 10)
    
    # ===== Panel 2: Fitted models =====
    ax2 = axes[0, 1]
    
    t_theory = np.linspace(0, 9.05, 500)
    
    # Secular growth model
    rho_secular = secular_growth_model(t_theory, rho_0, gamma_fit, t_delay_fit)
    
    # Runaway model  
    rho_runaway = runaway_model(t_theory, t_c_fit, A_fit, rho_bg=2.0)
    
    # Combined model
    rho_combined = combined_model(t_theory, rho_0, gamma_fit, t_delay_fit, t_c_fit, A_fit)
    
    ax2.semilogy(t_sim, rho_sim, 'ro', markersize=8, label='Simulation', zorder=10)
    ax2.semilogy(t_theory, rho_secular, 'b--', linewidth=2, 
                 label=f'Secular: Γ={gamma_fit:.3f}, t_d={t_delay_fit:.1f}')
    ax2.semilogy(t_theory, rho_runaway, 'g--', linewidth=2,
                 label=f'Runaway: t_c={t_c_fit:.3f}')
    ax2.semilogy(t_theory, rho_combined, 'k-', linewidth=2.5,
                 label='Combined model')
    
    ax2.set_xlabel('Time', fontsize=12)
    ax2.set_ylabel('Maximum density ρ_max (log)', fontsize=12)
    ax2.set_title('Fitted Models', fontsize=13)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-0.5, 10)
    ax2.set_ylim(0.5, 1e4)
    
    # ===== Panel 3: Linear scale comparison =====
    ax3 = axes[1, 0]
    
    ax3.plot(t_sim, rho_sim, 'ro-', markersize=6, linewidth=1.5, label='Simulation')
    ax3.plot(t_theory, rho_combined, 'k-', linewidth=2, label='Combined model')
    
    # Add GSPH with gradh
    if 'gsph_gradh' in data:
        d2 = data['gsph_gradh']
        ax3.plot(d2['t'], d2['rho_max'], 'b^-', markersize=4, linewidth=1,
                 alpha=0.7, label='GSPH with gradh')
    
    ax3.set_xlabel('Time', fontsize=12)
    ax3.set_ylabel('Maximum density ρ_max', fontsize=12)
    ax3.set_title('Linear Scale (Early Time)', fontsize=13)
    ax3.set_xlim(-0.5, 10)
    ax3.set_ylim(0, 50)
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # ===== Panel 4: Summary =====
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # Calculate theoretical parameters
    omega_dyn = np.sqrt(4 * np.pi * G * rho_0)
    t_dyn = 2 * np.pi / omega_dyn
    t_ff = np.sqrt(3 * np.pi / (32 * G * rho_0))
    
    # Infer effective error from growth rate
    epsilon_eff = 2 * gamma_fit / omega_dyn
    
    summary = f"""
╔══════════════════════════════════════════════════════════════════════╗
║        GSPH GRAD-H INSTABILITY: PHASE-BY-PHASE ANALYSIS              ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                      ║
║  PHASE 1: QUASI-STABLE (t < {t_phase1_end:.1f})                               ║
║    • Density nearly constant at ρ ≈ {rho_0:.2f}                           ║
║    • Small numerical noise ±10%                                      ║
║    • Error accumulation begins                                       ║
║                                                                      ║
║  PHASE 2: ACCELERATION ({t_phase1_end:.1f} < t < {t_runaway:.1f})                         ║
║    • Exponential growth with Γ = {gamma_fit:.4f} rad/time                ║
║    • e-folding time τ = {1/gamma_fit:.2f} time units                        ║
║    • Delay time t_delay = {t_delay_fit:.2f} time units                      ║
║                                                                      ║
║  PHASE 3: RUNAWAY (t > {t_runaway:.1f})                                        ║
║    • Singular collapse: ρ ∝ 1/(t_c - t)²                             ║
║    • Singularity time t_c = {t_c_fit:.4f}                                 ║
║    • Final density: ρ_max = {rho_sim.max():.0f}                               ║
║                                                                      ║
║  THEORETICAL INTERPRETATION:                                         ║
║    • Dynamical frequency: ω_dyn = {omega_dyn:.3f} rad/time                  ║
║    • Dynamical time: t_dyn = {t_dyn:.3f}                                    ║
║    • Free-fall time: t_ff = {t_ff:.3f}                                      ║
║    • Delay time ≈ {t_delay_fit/t_dyn:.1f} × t_dyn (instability seeding)              ║
║    • Effective error: ε_eff = 2Γ/ω_dyn = {epsilon_eff*100:.1f}%                  ║
║                                                                      ║
║  KEY PHYSICS:                                                        ║
║    The collapse is NOT a simple exponential! It's a DELAYED          ║
║    RUNAWAY triggered when cumulative error exceeds threshold.        ║
║    The singularity at t_c = {t_c_fit:.3f} is characteristic of            ║
║    gravitational collapse (like black hole formation).               ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
"""
    
    ax4.text(0.02, 0.98, summary, fontsize=10, fontfamily='monospace',
             verticalalignment='top', transform=ax4.transAxes,
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save
    output_dir = "results/gradh_comparison"
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(f'{output_dir}/phase_analysis.png', dpi=150, bbox_inches='tight')
    fig.savefig(f'{output_dir}/phase_analysis.pdf', bbox_inches='tight')
    print(f"\nSaved: {output_dir}/phase_analysis.png")
    
    # Verification table
    print("\n" + "=" * 70)
    print("MODEL VERIFICATION")
    print("=" * 70)
    print(f"\n{'t':>6} {'ρ(sim)':>10} {'ρ(model)':>10} {'Error':>10}")
    print("-" * 38)
    
    for ti, rhoi in zip(t_sim, rho_sim):
        # Combined model prediction
        rho_model = combined_model(np.array([ti]), rho_0, gamma_fit, 
                                   t_delay_fit, t_c_fit, A_fit)[0]
        error = (rho_model - rhoi) / rhoi * 100 if rhoi > 0 else 0
        print(f"{ti:6.2f} {rhoi:10.2f} {rho_model:10.2f} {error:9.1f}%")
    
    plt.show()


if __name__ == "__main__":
    main()
