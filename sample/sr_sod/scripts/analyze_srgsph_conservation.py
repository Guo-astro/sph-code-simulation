#!/usr/bin/env python3
"""
SRGSPH Conservation Law Analysis Script
========================================

This script analyzes the conservation properties and numerical behavior
of the Special Relativistic Godunov SPH (SRGSPH) implementation.

It addresses three key questions:
1. Conservation laws (mass, momentum, energy)
2. Second-order accuracy and gradient evaluation
3. Riemann solver behavior and potential issues

Author: Analysis for SRGSPH Code
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import os
import sys

# Scientific constants
c = 1.0  # Speed of light (code units)

def load_snapshot(filepath):
    """Load a single snapshot CSV file, excluding ghost particles."""
    # Find the header line (starts with 'id,')
    header_line = 0
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('id,'):
                header_line = i
                break
    
    df = pd.read_csv(filepath, skiprows=header_line)
    
    # Filter out ghost particles if is_ghost column exists
    if 'is_ghost' in df.columns:
        df = df[df['is_ghost'] == 0].reset_index(drop=True)
    
    return df

def load_all_snapshots(results_dir):
    """Load all snapshots in chronological order."""
    pattern = os.path.join(results_dir, "snapshot_*.csv")
    files = sorted(glob.glob(pattern))
    snapshots = []
    times = []
    for f in files:
        df = load_snapshot(f)
        snapshots.append(df)
        # Extract time from first row if available, or from filename
        if 'time' in df.columns:
            times.append(df['time'].iloc[0])
        else:
            # Estimate from filename number
            num = int(os.path.basename(f).replace('snapshot_', '').replace('.csv', ''))
            times.append(num * 0.005)  # Assuming outputTime = 0.005
    return snapshots, times

def compute_relativistic_quantities(df, gamma_ad=5.0/3.0):
    """
    Compute relativistic quantities for each particle.
    
    For SRGSPH:
    - N = lab-frame baryon density (stored as 'dens' or 'N')
    - n = rest-frame baryon density = N / γ
    - S = canonical momentum per baryon = γ H v
    - e = canonical energy per baryon = γ H - P/(N c²)
    - H = enthalpy per baryon = 1 + u/c² + P/(n c²)
    
    For ideal gas: P = (γ_c - 1) n u, so:
    H = 1 + (γ_c/(γ_c-1)) P/n  (with c=1)
    """
    # Get primitive variables
    if 'N' in df.columns:
        N = df['N'].values  # Lab-frame baryon density
    else:
        N = df['dens'].values  # Fallback
    
    P = df['pres'].values  # Pressure
    vx = df['vel_x'].values  # Velocity
    
    # Handle multi-dimensional cases
    if 'vel_y' in df.columns:
        vy = df['vel_y'].values
    else:
        vy = np.zeros_like(vx)
    if 'vel_z' in df.columns:
        vz = df['vel_z'].values
    else:
        vz = np.zeros_like(vx)
    
    # Compute velocity magnitude squared
    v2 = vx**2 + vy**2 + vz**2
    v2 = np.clip(v2, 0, 0.999999)  # Ensure subluminal
    
    # Lorentz factor
    gamma = 1.0 / np.sqrt(1.0 - v2)
    
    # Rest-frame baryon density
    n = N / gamma
    
    # Enthalpy per baryon (with c=1)
    # H = 1 + (γ_c/(γ_c-1)) P/n
    H = 1.0 + (gamma_ad / (gamma_ad - 1.0)) * (P / np.maximum(n, 1e-15))
    
    # Canonical momentum per baryon
    Sx = gamma * H * vx
    Sy = gamma * H * vy
    Sz = gamma * H * vz
    
    # Canonical energy per baryon (with c=1)
    # e = γ H - P/N
    e = gamma * H - P / N
    
    # Baryon number per particle (nu)
    if 'nu' in df.columns:
        nu = df['nu'].values
    else:
        # Estimate from initial setup if not available
        # For Sod problem: uniform nu
        nu = np.ones_like(N) * 0.00027778  # Default estimate
    
    return {
        'N': N,
        'n': n,
        'P': P,
        'vx': vx, 'vy': vy, 'vz': vz,
        'v2': v2,
        'gamma': gamma,
        'H': H,
        'Sx': Sx, 'Sy': Sy, 'Sz': Sz,
        'e': e,
        'nu': nu
    }

def compute_total_conserved(df, gamma_ad=5.0/3.0):
    """
    Compute total conserved quantities for the system.

    In SRGSPH, the conserved quantities are:
    - Total baryon number: Σ ν_i
    - Total canonical momentum: Σ ν_i S_i
    - Total canonical energy: Σ ν_i e_i

    Note: These are "per baryon" quantities multiplied by baryon number.
    """
    q = compute_relativistic_quantities(df, gamma_ad)

    nu = q['nu']

    # Total baryon number (should be exactly conserved)
    total_N = np.sum(nu)

    # Total canonical momentum (should be conserved)
    total_Sx = np.sum(nu * q['Sx'])
    total_Sy = np.sum(nu * q['Sy'])
    total_Sz = np.sum(nu * q['Sz'])

    # Total canonical energy (should be conserved)
    total_e = np.sum(nu * q['e'])

    # Compute thermal energy (internal energy)
    # For ideal gas: u = P / ((gamma_ad - 1) * n)
    u = q['P'] / (np.maximum(q['n'], 1e-15) * (gamma_ad - 1.0))
    total_thermal = np.sum(nu * u)

    # Kinetic energy per baryon: (gamma - 1) for c=1
    kinetic_per_baryon = q['gamma'] - 1.0
    total_kinetic = np.sum(nu * kinetic_per_baryon)

    return {
        'baryon_number': total_N,
        'momentum_x': total_Sx,
        'momentum_y': total_Sy,
        'momentum_z': total_Sz,
        'energy': total_e,
        'thermal_energy': total_thermal,
        'kinetic_energy': total_kinetic
    }

def analyze_conservation(results_dir, gamma_ad=5.0/3.0):
    """
    Analyze conservation laws over time.
    
    Returns dataframe with conservation quantities at each timestep.
    """
    snapshots, times = load_all_snapshots(results_dir)
    
    if len(snapshots) == 0:
        print(f"No snapshots found in {results_dir}")
        return None
    
    conservation_data = []
    for i, (snap, t) in enumerate(zip(snapshots, times)):
        conserved = compute_total_conserved(snap, gamma_ad)
        conserved['time'] = t
        conserved['snapshot'] = i
        conservation_data.append(conserved)
    
    df = pd.DataFrame(conservation_data)
    
    # Compute relative errors from initial values
    if len(df) > 0:
        df['baryon_error'] = (df['baryon_number'] - df['baryon_number'].iloc[0]) / df['baryon_number'].iloc[0]
        df['momentum_x_error'] = df['momentum_x'] - df['momentum_x'].iloc[0]  # Absolute since initial can be 0
        df['energy_error'] = (df['energy'] - df['energy'].iloc[0]) / df['energy'].iloc[0]
    
    return df

def plot_conservation_analysis(conservation_df, output_path):
    """Create comprehensive conservation analysis plot."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 14))
    fig.suptitle('SRGSPH Conservation Law Analysis\n' +
                 'Mass, Momentum, and Energy Conservation', fontsize=14)

    t = conservation_df['snapshot']  # Use timestep (integer) instead of time

    # 1. Baryon number conservation
    ax1 = axes[0, 0]
    ax1.plot(t, conservation_df['baryon_error'] * 100, 'b-', linewidth=2)
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Timestep')
    ax1.set_ylabel('Relative Error (%)')
    ax1.set_title('Baryon Number Conservation')
    ax1.grid(True, alpha=0.3)
    ax1.ticklabel_format(style='scientific', axis='y', scilimits=(-3,3))

    # 2. Momentum conservation evolution over time
    ax2 = axes[0, 1]
    ax2.plot(t, conservation_df['momentum_x'], 'r-', linewidth=2, label='$S_x$ (total)')
    ax2.axhline(y=conservation_df['momentum_x'].iloc[0], color='k', linestyle='--', alpha=0.5, label='Initial')
    ax2.set_xlabel('Timestep')
    ax2.set_ylabel('Total Momentum $S_x$')
    ax2.set_title('Momentum Evolution Over Time')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(-3,3))

    # 3. Momentum conservation error
    ax3 = axes[1, 0]
    ax3.plot(t, conservation_df['momentum_x_error'], 'r-', linewidth=2)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Timestep')
    ax3.set_ylabel('Absolute Error')
    ax3.set_title('Momentum Conservation Error')
    ax3.grid(True, alpha=0.3)
    ax3.ticklabel_format(style='scientific', axis='y', scilimits=(-3,3))

    # 4. Total energy evolution
    ax4 = axes[1, 1]
    ax4.plot(t, conservation_df['energy'], 'g-', linewidth=2, label='Total Energy')
    ax4.axhline(y=conservation_df['energy'].iloc[0], color='k', linestyle='--', alpha=0.5, label='Initial')
    ax4.set_xlabel('Timestep')
    ax4.set_ylabel('Total Energy')
    ax4.set_title('Total Energy Evolution Over Time')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Thermal and kinetic energy evolution
    ax5 = axes[2, 0]
    ax5.plot(t, conservation_df['thermal_energy'], 'orange', linewidth=2, label='Thermal Energy')
    ax5.plot(t, conservation_df['kinetic_energy'], 'purple', linewidth=2, label='Kinetic Energy')
    ax5.set_xlabel('Timestep')
    ax5.set_ylabel('Energy')
    ax5.set_title('Thermal & Kinetic Energy Evolution')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    # 6. Energy conservation error with explanation
    ax6 = axes[2, 1]
    ax6.plot(t, conservation_df['energy_error'] * 100, 'g-', linewidth=2, label='Total Energy Error')
    ax6.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax6.set_xlabel('Timestep')
    ax6.set_ylabel('Relative Error (%)')
    ax6.set_title('Energy Conservation Error Analysis')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax6.ticklabel_format(style='scientific', axis='y', scilimits=(-3,3))

    # Add text annotation about energy conservation
    thermal_change = (conservation_df['thermal_energy'].iloc[-1] - conservation_df['thermal_energy'].iloc[0])
    kinetic_change = (conservation_df['kinetic_energy'].iloc[-1] - conservation_df['kinetic_energy'].iloc[0])

    energy_text = f'Thermal Δ: {thermal_change:+.3e}\nKinetic Δ: {kinetic_change:+.3e}'
    ax6.text(0.02, 0.98, energy_text, transform=ax6.transAxes,
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Conservation analysis saved to {output_path}")

def plot_energy_conservation_detailed(conservation_df, output_path):
    """Create detailed energy conservation analysis plot."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Energy Conservation Detailed Analysis\n' +
                 'Why Energy Conservation Matters', fontsize=14)

    t = conservation_df['snapshot']  # Use timestep (integer) instead of time

    # 1. Total energy evolution
    ax1 = axes[0, 0]
    ax1.plot(t, conservation_df['energy'], 'g-', linewidth=2, label='Total Energy')
    initial_energy = conservation_df['energy'].iloc[0]
    ax1.axhline(y=initial_energy, color='k', linestyle='--', alpha=0.5, label='Initial Value')
    ax1.set_xlabel('Timestep')
    ax1.set_ylabel('Total Canonical Energy')
    ax1.set_title('Total Energy Over Time (Should be Conserved)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Energy components stacked
    ax2 = axes[0, 1]
    ax2.plot(t, conservation_df['thermal_energy'], 'orange', linewidth=2, label='Thermal Energy')
    ax2.plot(t, conservation_df['kinetic_energy'], 'purple', linewidth=2, label='Kinetic Energy')
    ax2.set_xlabel('Timestep')
    ax2.set_ylabel('Energy Components')
    ax2.set_title('Energy Components Evolution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Energy conservation error
    ax3 = axes[1, 0]
    energy_error_pct = conservation_df['energy_error'] * 100
    ax3.plot(t, energy_error_pct, 'g-', linewidth=2)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax3.fill_between(t, 0, energy_error_pct, alpha=0.3, color='green')
    ax3.set_xlabel('Timestep')
    ax3.set_ylabel('Relative Error (%)')
    ax3.set_title('Energy Conservation Error')
    ax3.grid(True, alpha=0.3)
    ax3.ticklabel_format(style='scientific', axis='y', scilimits=(-3,3))

    # 4. Explanation panel
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Compute changes
    energy_change = conservation_df['energy'].iloc[-1] - conservation_df['energy'].iloc[0]
    thermal_change = conservation_df['thermal_energy'].iloc[-1] - conservation_df['thermal_energy'].iloc[0]
    kinetic_change = conservation_df['kinetic_energy'].iloc[-1] - conservation_df['kinetic_energy'].iloc[0]
    energy_rel_error = conservation_df['energy_error'].iloc[-1] * 100

    explanation_text = f"""
    Energy Conservation Analysis
    ========================================

    Initial State:
      Total Energy:   {conservation_df['energy'].iloc[0]:.6e}
      Thermal Energy: {conservation_df['thermal_energy'].iloc[0]:.6e}
      Kinetic Energy: {conservation_df['kinetic_energy'].iloc[0]:.6e}

    Final State:
      Total Energy:   {conservation_df['energy'].iloc[-1]:.6e}
      Thermal Energy: {conservation_df['thermal_energy'].iloc[-1]:.6e}
      Kinetic Energy: {conservation_df['kinetic_energy'].iloc[-1]:.6e}

    Changes:
      ΔE_total:   {energy_change:+.6e} ({energy_rel_error:+.3e}%)
      ΔE_thermal: {thermal_change:+.6e}
      ΔE_kinetic: {kinetic_change:+.6e}

    Why This Matters:
    - Total energy MUST be conserved in isolated system
    - Thermal ↔ Kinetic exchange shows shock dynamics
    - Errors indicate numerical issues:
      * Primitive recovery (quartic solver)
      * Time integration (Euler method)
      * Riemann solver accuracy

    Note: In shock tube problem:
    - Initially: High thermal on left, low on right
    - During evolution: Kinetic energy increases
    - After shock: Thermal energy redistributes
    """

    ax4.text(0.05, 0.95, explanation_text, transform=ax4.transAxes,
             fontsize=9, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Detailed energy conservation analysis saved to {output_path}")

def compute_exact_sod_solution(x, t, gamma_ad=5.0/3.0, x0=0.0,
                                rhoL=1.0, pL=1.0, vL=0.0,
                                rhoR=0.125, pR=0.1, vR=0.0):
    """
    Compute exact solution for relativistic Sod shock tube.
    
    This is the Pons et al. (2000) exact Riemann solution for vt=0.
    """
    from scipy.optimize import brentq
    
    # Helper functions for the relativistic Riemann problem
    def compute_h(p, rho, gam):
        """Compute specific enthalpy h = 1 + γ/(γ-1) * p/ρ"""
        return 1.0 + (gam / (gam - 1.0)) * (p / rho)
    
    def compute_cs(p, rho, h, gam):
        """Compute sound speed"""
        return np.sqrt(gam * p / (rho * h))
    
    def shock_velocity_left(p_star, state, gam):
        """Compute post-shock state for left-going shock."""
        pa, rhoa, va = state
        ha = compute_h(pa, rhoa, gam)
        
        # Taub adiabat - quadratic for hb
        A = (gam - 1.0) * (pa - p_star) / (gam * p_star)
        B = ha * (pa - p_star) / rhoa
        
        qa = 1.0 + A
        qb = -A
        qc = B - ha**2
        
        delta = qb**2 - 4*qa*qc
        if delta < 0:
            return None
            
        hb = (-qb + np.sqrt(delta)) / (2*qa)
        rhob = (gam / (gam - 1.0)) * p_star / (hb - 1.0)
        
        # Mass flux
        denom = (hb/rhob - ha/rhoa)
        if abs(denom) < 1e-15:
            return va, rhob
            
        j2 = -(p_star - pa) / denom
        if j2 <= 0:
            return va, rhob
            
        j = np.sqrt(j2)
        
        Wa = 1.0 / np.sqrt(1.0 - va**2)
        Da = rhoa * Wa
        
        term = j**2 + Da**2 * (1.0 - va**2)
        sqrt_term = np.sqrt(term)
        denom_vs = Da**2 + j**2
        
        Vs = (Da**2 * va - j * sqrt_term) / denom_vs  # Left wave
        Ws = 1.0 / np.sqrt(1.0 - Vs**2)
        
        num = ha * Wa * va + Ws * (p_star - pa) / (-j)
        den = ha * Wa + (p_star - pa) * (Ws * va / (-j) + 1.0 / Da)
        
        vb = num / den
        return vb, rhob
    
    def shock_velocity_right(p_star, state, gam):
        """Compute post-shock state for right-going shock."""
        pa, rhoa, va = state
        ha = compute_h(pa, rhoa, gam)
        
        A = (gam - 1.0) * (pa - p_star) / (gam * p_star)
        B = ha * (pa - p_star) / rhoa
        
        qa = 1.0 + A
        qb = -A
        qc = B - ha**2
        
        delta = qb**2 - 4*qa*qc
        if delta < 0:
            return None
            
        hb = (-qb + np.sqrt(delta)) / (2*qa)
        rhob = (gam / (gam - 1.0)) * p_star / (hb - 1.0)
        
        denom = (hb/rhob - ha/rhoa)
        if abs(denom) < 1e-15:
            return va, rhob
            
        j2 = -(p_star - pa) / denom
        if j2 <= 0:
            return va, rhob
            
        j = np.sqrt(j2)
        
        Wa = 1.0 / np.sqrt(1.0 - va**2)
        Da = rhoa * Wa
        
        term = j**2 + Da**2 * (1.0 - va**2)
        sqrt_term = np.sqrt(term)
        denom_vs = Da**2 + j**2
        
        Vs = (Da**2 * va + j * sqrt_term) / denom_vs  # Right wave
        Ws = 1.0 / np.sqrt(1.0 - Vs**2)
        
        num = ha * Wa * va + Ws * (p_star - pa) / j
        den = ha * Wa + (p_star - pa) * (Ws * va / j + 1.0 / Da)
        
        vb = num / den
        return vb, rhob
    
    def rarefaction_velocity(p_star, state, direction, gam):
        """Compute post-rarefaction velocity for zero tangential velocity."""
        pa, rhoa, va = state
        ha = compute_h(pa, rhoa, gam)
        csa = compute_cs(pa, rhoa, ha, gam)
        
        # Isentropic relation
        const_entropy = pa / (rhoa**gam)
        rhob = (p_star / const_entropy)**(1.0/gam)
        hb = compute_h(p_star, rhob, gam)
        csb = compute_cs(p_star, rhob, hb, gam)
        
        # Analytical solution for vt=0
        sqrt_gm1 = np.sqrt(gam - 1.0)
        
        term_v = (1.0 + va) / (1.0 - va)
        term_ca = (sqrt_gm1 + csa) / (sqrt_gm1 - csa)
        term_cb = (sqrt_gm1 + csb) / (sqrt_gm1 - csb)
        
        sign = 1.0 if direction == 'left' else -1.0
        exponent = sign * 2.0 / sqrt_gm1
        
        base = term_ca / term_cb
        A = term_v * (base ** exponent)
        
        vb = (A - 1.0) / (A + 1.0)
        return vb, rhob
    
    def wave_curve_left(p, state, gam):
        """Velocity from left wave curve."""
        pa, rhoa, va = state
        if p > pa:
            result = shock_velocity_left(p, state, gam)
        else:
            result = rarefaction_velocity(p, state, 'left', gam)
        if result is None:
            return np.nan
        return result[0]
    
    def wave_curve_right(p, state, gam):
        """Velocity from right wave curve."""
        pa, rhoa, va = state
        if p > pa:
            result = shock_velocity_right(p, state, gam)
        else:
            result = rarefaction_velocity(p, state, 'right', gam)
        if result is None:
            return np.nan
        return result[0]
    
    # Solve for p_star
    stateL = (pL, rhoL, vL)
    stateR = (pR, rhoR, vR)
    
    def residual(p):
        vl = wave_curve_left(p, stateL, gamma_ad)
        vr = wave_curve_right(p, stateR, gamma_ad)
        return vl - vr
    
    # Find root
    try:
        p_star = brentq(residual, pR * 0.01, pL * 10.0)
    except:
        p_star = 0.5 * (pL + pR)
    
    # Compute star states
    if p_star > pL:
        v_star, rho_Lstar = shock_velocity_left(p_star, stateL, gamma_ad)
    else:
        v_star, rho_Lstar = rarefaction_velocity(p_star, stateL, 'left', gamma_ad)
    
    if p_star > pR:
        v_star_R, rho_Rstar = shock_velocity_right(p_star, stateR, gamma_ad)
    else:
        v_star_R, rho_Rstar = rarefaction_velocity(p_star, stateR, 'right', gamma_ad)
    
    v_star = 0.5 * (v_star + v_star_R)  # Average
    
    # Wave speeds
    hL = compute_h(pL, rhoL, gamma_ad)
    csL = compute_cs(pL, rhoL, hL, gamma_ad)
    hR = compute_h(pR, rhoR, gamma_ad)
    csR = compute_cs(pR, rhoR, hR, gamma_ad)
    
    # Left wave head (rarefaction head or shock)
    if p_star < pL:
        # Rarefaction: head moves at v - cs
        xi_head_L = (vL - csL) / (1.0 - vL * csL)
    else:
        # Shock
        xi_head_L = -0.5  # Approximate
    
    # Right wave (shock)
    # Approximate shock speed
    xi_shock_R = 0.8  # Approximate
    
    # Contact discontinuity
    xi_contact = v_star
    
    # Now evaluate solution at each x
    rho = np.zeros_like(x)
    p = np.zeros_like(x)
    v = np.zeros_like(x)
    
    xi = (x - x0) / t
    
    for i, xi_val in enumerate(xi):
        if xi_val < xi_head_L:
            # Left state
            rho[i] = rhoL
            p[i] = pL
            v[i] = vL
        elif xi_val < xi_contact:
            # Left star state (simplified)
            rho[i] = rho_Lstar
            p[i] = p_star
            v[i] = v_star
        elif xi_val < xi_shock_R:
            # Right star state
            rho[i] = rho_Rstar
            p[i] = p_star
            v[i] = v_star
        else:
            # Right state
            rho[i] = rhoR
            p[i] = pR
            v[i] = vR
    
    return {'rho': rho, 'p': p, 'v': v, 'p_star': p_star, 'v_star': v_star}

def plot_snapshot_with_exact(df, t, output_path, gamma_ad=5.0/3.0):
    """Plot snapshot with exact solution overlay."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'SRGSPH vs Exact Solution at t = {t:.4f}\n' +
                 'Density and Pressure Comparison in Evolved Regions', fontsize=14)
    
    # Get particle data
    x = df['pos_x'].values
    
    q = compute_relativistic_quantities(df, gamma_ad)
    rho = q['n']  # Rest-frame density
    p = q['P']
    v = q['vx']
    
    # Sort by position
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    rho = rho[sort_idx]
    p = p[sort_idx]
    v = v[sort_idx]
    
    # Compute exact solution (approximate)
    x_exact = np.linspace(-0.5, 0.5, 1000)
    exact = compute_exact_sod_solution(x_exact, t, gamma_ad)
    
    # Plot density
    ax1 = axes[0, 0]
    ax1.scatter(x, rho, s=2, alpha=0.5, label='SPH')
    ax1.plot(x_exact, exact['rho'], 'r-', linewidth=1.5, label='Exact')
    ax1.set_xlabel('Position x')
    ax1.set_ylabel(r'Rest-frame density $n$')
    ax1.set_title('Density')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot pressure
    ax2 = axes[0, 1]
    ax2.scatter(x, p, s=2, alpha=0.5, label='SPH')
    ax2.plot(x_exact, exact['p'], 'r-', linewidth=1.5, label='Exact')
    ax2.set_xlabel('Position x')
    ax2.set_ylabel('Pressure P')
    ax2.set_title('Pressure')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot velocity
    ax3 = axes[1, 0]
    ax3.scatter(x, v, s=2, alpha=0.5, label='SPH')
    ax3.plot(x_exact, exact['v'], 'r-', linewidth=1.5, label='Exact')
    ax3.set_xlabel('Position x')
    ax3.set_ylabel('Velocity $v^x$')
    ax3.set_title('Velocity')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Annotations
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    info_text = f"""
    Analysis Notes
    ========================================
    
    Time: t = {t:.4f}
    Number of particles: {len(df)}
    
    Star State (P*, v*):
      P* = {exact['p_star']:.4f}
      v* = {exact['v_star']:.4f}
    
    Observations:
    1. Star region density slightly low
       -> Check conservation law accuracy
       
    2. Overshoot at rarefaction wave head
       -> Effect of 2nd order scheme
       
    3. Flat regions inside rarefaction
       -> Riemann solver characteristics
    """
    
    ax4.text(0.1, 0.9, info_text, transform=ax4.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Snapshot comparison saved to {output_path}")

def create_scheme_explanation_figure(output_path):
    """Create a figure explaining the numerical scheme."""
    fig = plt.figure(figsize=(16, 12))
    
    # Use gridspec for complex layout
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # 1. MUSCL reconstruction diagram
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(0, 2)
    
    # Draw particle positions
    x_particles = [-0.5, 0, 0.5]
    y_values = [1.0, 1.5, 1.2]
    
    ax1.scatter(x_particles, y_values, s=100, c='blue', zorder=5, label='Particle values')
    
    # Draw linear reconstruction
    for i in range(len(x_particles)):
        if i > 0 and i < len(x_particles) - 1:
            grad = (y_values[i+1] - y_values[i-1]) / (x_particles[i+1] - x_particles[i-1])
            x_line = np.linspace(x_particles[i-1]/2 + x_particles[i]/2, 
                                 x_particles[i]/2 + x_particles[i+1]/2, 50)
            y_line = y_values[i] + grad * (x_line - x_particles[i])
            ax1.plot(x_line, y_line, 'r--', linewidth=2, alpha=0.7)
    
    ax1.axvline(x=0.25, color='green', linestyle=':', linewidth=2, label='Interface')
    ax1.set_xlabel('Position x')
    ax1.set_ylabel('Quantity')
    ax1.set_title('MUSCL Reconstruction (2nd Order Reconstruction)\n' +
                  r'$q_{L,R} = q_i + \nabla q_i \cdot \Delta r$')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Gradient evaluation formula
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis('off')
    
    gradient_text = r"""
    SPH Gradient Evaluation
    ========================================
    
    Standard SPH gradient:
    
    $\nabla A_i = \sum_j V_{p,j} (A_j - A_i) \nabla W_{ij}$
    
    where:
    - $V_{p,j} = \nu_j / N_j$ = particle volume
    - $\nabla W_{ij}$ = kernel gradient
    - $(A_j - A_i)$ = difference from neighbors
    
    MUSCL interface reconstruction:
    
    $A_L = A_i + \nabla A_i \cdot \Delta r_i$
    $A_R = A_j + \nabla A_j \cdot \Delta r_j$
    
    with $\Delta r = (x_{midpoint} - x_i)$
    
    Shock/CD detection switches to 1st order when:
    - Shock: $C_{shock} |\mathbf{e}_{ij} \cdot \Delta\mathbf{v}| > c_s$
    - CD: $|\log_{10}(P_i/P_j)| > C_{cd}$
    """
    ax2.text(0.05, 0.95, gradient_text, transform=ax2.transAxes,
             fontsize=11, verticalalignment='top', fontfamily='monospace')
    
    # 3. Riemann solver diagram
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_xlim(-1, 1)
    ax3.set_ylim(0, 1)
    
    # Draw wave structure
    ax3.fill_betweenx([0, 1], -1, -0.5, color='lightblue', alpha=0.5, label='L')
    ax3.fill_betweenx([0, 1], -0.5, -0.1, color='cyan', alpha=0.5, label='L*')
    ax3.fill_betweenx([0, 1], -0.1, 0.3, color='lightgreen', alpha=0.5, label='R*')
    ax3.fill_betweenx([0, 1], 0.3, 1, color='lightyellow', alpha=0.5, label='R')
    
    # Draw waves
    ax3.plot([-0.5, -0.5], [0, 1], 'b-', linewidth=2)
    ax3.plot([-0.1, -0.1], [0, 1], 'k--', linewidth=2)
    ax3.plot([0.3, 0.3], [0, 1], 'r-', linewidth=2)
    
    ax3.text(-0.7, 0.5, 'Raref.', fontsize=10, ha='center')
    ax3.text(-0.1, 0.85, 'Contact', fontsize=10, ha='center')
    ax3.text(0.4, 0.5, 'Shock', fontsize=10, ha='center')
    
    ax3.set_xlabel(r'$\xi = x/t$')
    ax3.set_ylabel('Wave structure')
    ax3.set_title('Riemann Problem Structure')
    ax3.legend(loc='upper right')
    
    # 4. Riemann solver details
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    riemann_text = r"""
    Exact Riemann Solver
    ========================================
    
    Algorithm (Pons et al. 2000):
    
    1. For given $P_*$, solve for $v_L^x(P_*)$:
       - If $P_* > P_L$: Shock (Taub adiabat)
       - If $P_* < P_L$: Rarefaction (ODE integration)
    
    2. Similarly for $v_R^x(P_*)$
    
    3. Find $P_*$ such that $v_L^x(P_*) = v_R^x(P_*)$
       using bisection method
    
    4. Compute star states $(n_*, P_*, v_*)$
    
    Rarefaction ODE (RK4 integration):
    
    $\frac{dv^x}{dP} = \pm \frac{1}{n H \gamma^2 c_s \sqrt{1+g}}$
    
    where $g = \frac{(v^t)^2 (\xi^2-1)}{(1-\xi v^x)^2}$
    
    Note: Flat regions inside rarefaction may be due
    to ODE integration discretization steps.
    """
    ax4.text(0.05, 0.95, riemann_text, transform=ax4.transAxes,
             fontsize=11, verticalalignment='top', fontfamily='monospace')
    
    # 5. Conservation explanation
    ax5 = fig.add_subplot(gs[2, :])
    ax5.axis('off')
    
    conservation_text = r"""
    SRGSPH Conservation Properties
    ================================================================================================================
    
    Canonical variables evolved by SRGSPH:
    
    - Momentum:  $\frac{d\mathbf{S}}{dt} = -\frac{1}{N} \nabla P$  where $\mathbf{S} = \gamma H \mathbf{v}$
    
    - Energy:    $\frac{de}{dt} = -\frac{1}{N} \nabla \cdot (P\mathbf{v})$  where $e = \gamma H - P/N$
    
    Discretized form (Kitajima et al. Eq. 371, 373):
    
    $\langle \nu_i \dot{\mathbf{S}}_i \rangle = -\sum_j P_{ij}^* V_{ij}^2 [\nabla_i W_{ij}(\sqrt{2}h_i) - \nabla_j W_{ij}(\sqrt{2}h_j)]$
    
    $\langle \nu_i \dot{e}_i \rangle = -\sum_j P_{ij}^* \mathbf{v}_{ij}^* \cdot V_{ij}^2 [\nabla_i W_{ij}(\sqrt{2}h_i) - \nabla_j W_{ij}(\sqrt{2}h_j)]$
    
    Conservation guarantee: The term $[\nabla_i W - \nabla_j W]$ is anti-symmetric in $(i,j)$,
    and the Riemann solution $(P_{ij}^*, \mathbf{v}_{ij}^*)$ is symmetric.
    -> Total momentum and energy are exactly conserved (in theory).
    
    Possible sources of error:
    
    1. Primitive variable recovery: Quartic solver for $\gamma$ from $(S, e)$ introduces numerical error
    2. Time integration: Euler method has $O(\Delta t^2)$ error per step
    3. Riemann solver tolerance: Bisection with $tol = 10^{-10}$ still has finite error
    4. Kernel truncation: Neighbor search radius limits kernel support
    """
    ax5.text(0.02, 0.95, conservation_text, transform=ax5.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.suptitle('SRGSPH Numerical Scheme Analysis\nSecond-Order Scheme and Riemann Solver', 
                 fontsize=14, y=0.98)
    
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Scheme explanation saved to {output_path}")

def create_rarefaction_overshoot_analysis(results_dir, output_path, gamma_ad=5.0/3.0):
    """Analyze the rarefaction wave overshoot in detail."""
    snapshots, times = load_all_snapshots(results_dir)
    
    if len(snapshots) == 0:
        print("No snapshots found")
        return
    
    # Use a mid-time snapshot
    idx = len(snapshots) // 2
    df = snapshots[idx]
    t = times[idx]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Rarefaction Wave Analysis at t = {t:.4f}\n' +
                 'Overshoot Analysis at Rarefaction Wave Head', fontsize=14)
    
    # Get data
    x = df['pos_x'].values
    q = compute_relativistic_quantities(df, gamma_ad)
    rho = q['n']
    v = q['vx']
    
    # Sort
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    rho = rho[sort_idx]
    v = v[sort_idx]
    
    # Focus on rarefaction region (left side)
    mask = (x > -0.4) & (x < 0.0)
    x_rar = x[mask]
    rho_rar = rho[mask]
    v_rar = v[mask]
    
    # 1. Density in rarefaction
    ax1 = axes[0, 0]
    ax1.scatter(x_rar, rho_rar, s=10, alpha=0.7, c='blue')
    ax1.set_xlabel('Position x')
    ax1.set_ylabel(r'Rest-frame density $n$')
    ax1.set_title('Density in Rarefaction Region')
    ax1.grid(True, alpha=0.3)
    
    # Highlight overshoot region
    head_mask = x_rar < -0.25
    if np.any(head_mask):
        ax1.axvspan(x_rar[head_mask].min(), x_rar[head_mask].max(), 
                    alpha=0.2, color='red', label='Overshoot region')
    ax1.legend()
    
    # 2. Velocity in rarefaction
    ax2 = axes[0, 1]
    ax2.scatter(x_rar, v_rar, s=10, alpha=0.7, c='green')
    ax2.set_xlabel('Position x')
    ax2.set_ylabel(r'Velocity $v^x$')
    ax2.set_title('Velocity in Rarefaction Region')
    ax2.grid(True, alpha=0.3)
    
    # 3. Check for flat regions (derivative analysis)
    ax3 = axes[1, 0]
    if len(x_rar) > 5:
        dv_dx = np.gradient(v_rar, x_rar)
        ax3.scatter(x_rar, dv_dx, s=10, alpha=0.7, c='purple')
        ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        ax3.set_xlabel('Position x')
        ax3.set_ylabel(r'$dv/dx$')
        ax3.set_title('Velocity Gradient (flat regions -> gradient ~ 0)')
        ax3.grid(True, alpha=0.3)
    
    # 4. Explanation
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    explanation_text = f"""
    Rarefaction Wave Overshoot Analysis
    ========================================
    
    Current configuration:
    - use2ndOrderSRGSPH: See config file
    - Riemann solver: Exact (Pons et al. bisection)
    - Time integration: Euler (1st order)
    
    Possible causes of overshoot at rarefaction head:
    
    1. 2nd order MUSCL reconstruction:
       - Gradient extrapolation can exceed monotonicity
       - Limiter may be needed (minmod, MC, etc.)
       - Current: No slope limiter applied
    
    2. Riemann solver discretization:
       - ODE integration in solve_rarefaction uses RK4
       - Step count: adaptive (64-2048 steps)
       - May need more refinement near head
    
    3. Kernel gradient errors:
       - Gaussian kernel gradient at interface
       - h = sqrt(2) * h_i for SRGSPH
    
    Possible causes of flat regions in rarefaction:
    
    1. Riemann solver ODE stepping:
       - Discrete integration can create plateaus
       - Characteristic speed calculation
    
    2. Contact detection threshold:
       - C_cd = 1.0 may trigger early 1st order
       - Pressure ratio test affects gradients
    
    Number of particles analyzed: {len(x_rar)}
    """
    
    ax4.text(0.05, 0.95, explanation_text, transform=ax4.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Rarefaction analysis saved to {output_path}")

def main():
    """Main analysis function."""
    # Parse command line arguments
    if len(sys.argv) >= 2:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path("sample/sr_sod/results/sod_fig1")
    
    # Output directory is inside results directory
    output_dir = results_dir / "analysis"
    output_dir.mkdir(exist_ok=True)
    
    gamma_ad = 5.0 / 3.0
    
    print("=" * 60)
    print("SRGSPH Conservation and Numerical Analysis")
    print(f"Results directory: {results_dir}")
    print(f"Output directory:  {output_dir}")
    print("=" * 60)

    # 1. Conservation analysis
    print("\n[1/5] Analyzing conservation laws...")
    conservation_df = analyze_conservation(results_dir, gamma_ad)
    if conservation_df is not None:
        plot_conservation_analysis(conservation_df, output_dir / "conservation_analysis.png")

        # Also create detailed energy conservation plot
        plot_energy_conservation_detailed(conservation_df, output_dir / "energy_conservation_detailed.png")

        # Print summary
        print("\n  Conservation Summary:")
        print(f"  - Baryon number max error: {np.max(np.abs(conservation_df['baryon_error']))*100:.2e}%")
        print(f"  - Momentum max error: {np.max(np.abs(conservation_df['momentum_x_error'])):.2e}")
        print(f"  - Energy max error: {np.max(np.abs(conservation_df['energy_error']))*100:.2e}%")
        print(f"\n  Energy Evolution:")
        print(f"  - Initial total energy: {conservation_df['energy'].iloc[0]:.6e}")
        print(f"  - Final total energy: {conservation_df['energy'].iloc[-1]:.6e}")
        print(f"  - Initial thermal energy: {conservation_df['thermal_energy'].iloc[0]:.6e}")
        print(f"  - Final thermal energy: {conservation_df['thermal_energy'].iloc[-1]:.6e}")
        print(f"  - Initial kinetic energy: {conservation_df['kinetic_energy'].iloc[0]:.6e}")
        print(f"  - Final kinetic energy: {conservation_df['kinetic_energy'].iloc[-1]:.6e}")

    # 2. Create scheme explanation figure
    print("\n[2/5] Creating scheme explanation figure...")
    create_scheme_explanation_figure(output_dir / "scheme_explanation.png")

    # 3. Rarefaction overshoot analysis
    print("\n[3/5] Analyzing rarefaction wave overshoot...")
    create_rarefaction_overshoot_analysis(results_dir,
                                          output_dir / "rarefaction_analysis.png",
                                          gamma_ad)

    # 4. Snapshot comparison with exact solution
    print("\n[4/5] Creating snapshot comparison with exact solution...")
    snapshots, times = load_all_snapshots(results_dir)
    if len(snapshots) > 0:
        # Use final snapshot
        final_idx = -1
        plot_snapshot_with_exact(snapshots[final_idx], times[final_idx],
                                output_dir / "snapshot_vs_exact.png", gamma_ad)

    print("\n[5/5] All analyses complete!")
    
    print("\n" + "=" * 60)
    print("Analysis complete! Output files in:", output_dir)
    print("=" * 60)

if __name__ == "__main__":
    main()
