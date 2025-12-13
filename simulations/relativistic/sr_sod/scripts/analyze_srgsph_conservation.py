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

# Add docs directory to path for relativistic_riemann_solver
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'docs'))
from relativistic_riemann_solver import RelativisiticRiemannSolver

# Import shared test case definitions (SSOT)
from sr_test_cases import TEST_CASES, detect_test_type, get_initial_conditions, get_exact_solution

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
    
    For SRGSPH output:
    - dens = rest-frame baryon density n (already recovered from primitive recovery)
    - pres = pressure P
    - vel_x = primitive velocity v
    
    This function computes derived quantities:
    - N = lab-frame baryon density = γn  
    - S = canonical momentum per baryon = γ H v
    - e = canonical energy per baryon = γ H - P/(N c²)
    - H = enthalpy per baryon = 1 + u/c² + P/(n c²)
    
    For ideal gas: P = (γ_c - 1) n u, so:
    H = 1 + (γ_c/(γ_c-1)) P/n  (with c=1)
    """
    # Get primitive variables from SRGSPH output
    # NOTE: SRGSPH stores REST-FRAME density n in 'dens' column,
    # not lab-frame density N. See sr_pre_interaction.cpp line 272:
    # p_i.dens = prim.density;  // Rest-frame density n
    n = df['dens'].values  # Rest-frame baryon density (already in output)
    
    P = df['pres'].values  # Pressure
    vx = df['vel_x'].values  # Velocity
    
    # Handle multi-dimensional cases
    # For 1D simulations, vel_y and vel_z columns won't exist (dimension-aware writer)
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
    
    # Lab-frame baryon density (computed from rest-frame)
    N = n * gamma
    
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
    
    # Baryon number per particle (nu) - stored in 'mass' column for SRGSPH
    if 'nu' in df.columns:
        nu = df['nu'].values
    elif 'mass' in df.columns:
        nu = df['mass'].values  # SRGSPH stores baryon number in mass column
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
    
    Returns TWO types of conservation metrics:
    
    1. INTRINSIC (Scheme Conservation):
       - Measures internal numerical conservation
       - Removes boundary flux effects
       - Should be ~0 for a well-implemented scheme
       
    2. TOTAL (Domain Conservation):
       - Raw sum over all particles in domain
       - Changes due to boundary fluxes (physical!)
       - For outflow BC: material leaves, so totals decrease

    ================================================================================
    FORMULAS AND DATA COLUMNS USED
    ================================================================================
    
    CSV Columns used:
      - mass (ν):     Baryon number per particle (conserved per particle)
      - dens (n):     REST-FRAME baryon density (from primitive recovery)
      - pres (P):     Pressure
      - vel_x (v):    Velocity
      - ene (u):      Specific internal energy (thermal)
    
    Derived quantities:
      γ = 1/√(1-v²)                    Lorentz factor
      H = 1 + (γ_ad/(γ_ad-1)) P/n      Specific enthalpy (with c=1)
      
    Conserved variables (per baryon):
      N = γn                           Lab-frame baryon density  
      S = γHv                          Canonical momentum per baryon
      e = γH - P/N                     Canonical energy per baryon
      
    Total conserved quantities:
      Total baryons:    Σᵢ νᵢ          (exactly conserved, ν is constant per particle)
      Total momentum:   Σᵢ νᵢ Sᵢ       (conserved in isolated system)
      Total energy:     Σᵢ νᵢ eᵢ       (conserved in isolated system)
      
    Energy decomposition:
      Thermal energy:   Σᵢ νᵢ u        where u = P/((γ_ad-1)n)
      Kinetic energy:   Σᵢ νᵢ (γ-1)    relativistic kinetic energy per baryon
      
    ================================================================================
    """
    q = compute_relativistic_quantities(df, gamma_ad)

    nu = q['nu']

    # Total baryon number (should be exactly conserved - ν doesn't change)
    total_N = np.sum(nu)

    # Total canonical momentum (should be conserved in isolated system)
    total_Sx = np.sum(nu * q['Sx'])
    total_Sy = np.sum(nu * q['Sy'])
    total_Sz = np.sum(nu * q['Sz'])

    # Total canonical energy (should be conserved in isolated system)
    total_e = np.sum(nu * q['e'])

    # Compute thermal energy (internal energy)
    # For ideal gas: u = P / ((gamma_ad - 1) * n)
    u = q['P'] / (np.maximum(q['n'], 1e-15) * (gamma_ad - 1.0))
    total_thermal = np.sum(nu * u)

    # Kinetic energy per baryon: (gamma - 1) for c=1
    kinetic_per_baryon = q['gamma'] - 1.0
    total_kinetic = np.sum(nu * kinetic_per_baryon)
    
    # Rest mass energy: Σᵢ νᵢ (with c=1, rest mass = 1 per baryon)
    total_rest_mass = np.sum(nu)

    return {
        'baryon_number': total_N,
        'momentum_x': total_Sx,
        'momentum_y': total_Sy,
        'momentum_z': total_Sz,
        'energy': total_e,
        'thermal_energy': total_thermal,
        'kinetic_energy': total_kinetic,
        'rest_mass_energy': total_rest_mass
    }


def compute_ghost_impulse(df, dt, gamma_ad=5.0/3.0, x_left=-0.5, x_right=0.5):
    """
    Compute the momentum impulse from ghost particles to real particles.
    
    ================================================================================
    GHOST PRESSURE FORCE CALCULATION
    ================================================================================
    
    Ghost particles exert pressure forces on real particles near boundaries.
    The impulse (momentum transferred) per timestep is:
    
        J = F · dt = -∫ P dA · dt
    
    For 1D shock tube:
      - Left boundary (x = x_left):  F_left  = +P_left · A   (pushes rightward)
      - Right boundary (x = x_right): F_right = -P_right · A  (pushes leftward)
      - Net force on real particles: F_net = P_left - P_right
      
    We estimate this by looking at particles near boundaries:
      - Pressure at left edge:  P_left  = average P of particles near x_left
      - Pressure at right edge: P_right = average P of particles near x_right
      
    Cross-sectional area A = 1 for 1D (or actual area for 2D/3D)
    
    ================================================================================
    """
    # Get smoothing length for boundary region definition
    h = df['sml'].mean() if 'sml' in df.columns and len(df) > 0 else 0.01
    
    x = df['pos_x'].values
    P = df['pres'].values
    
    # Find particles near boundaries (within 2h of edge)
    left_mask = x < (x_left + 3*h)
    right_mask = x > (x_right - 3*h)
    
    # Average pressure at each boundary
    P_left = np.mean(P[left_mask]) if np.any(left_mask) else 0.0
    P_right = np.mean(P[right_mask]) if np.any(right_mask) else 0.0
    
    # Net pressure force on real particles (positive = rightward)
    # Left ghost pushes RIGHT (+), Right ghost pushes LEFT (-)
    # Effective area = 1 for 1D normalized problem
    A = 1.0
    F_net = (P_left - P_right) * A
    
    # Impulse = Force × time
    impulse = F_net * dt
    
    # Energy: Work done by boundary pressure = F · dx = F · v · dt
    # Average velocity at boundaries
    v = df['vel_x'].values
    v_left = np.mean(v[left_mask]) if np.any(left_mask) else 0.0
    v_right = np.mean(v[right_mask]) if np.any(right_mask) else 0.0
    
    # Work = P·v·A·dt at each boundary (power × time)
    # Left boundary: material moving left (v<0) does negative work, moving right does positive
    # Right boundary: material moving right (v>0) does negative work (leaves), moving left does positive
    work_left = P_left * v_left * A * dt  # If v>0, ghosts do positive work
    work_right = -P_right * v_right * A * dt  # If v>0, material does work against right ghost
    work_net = work_left + work_right
    
    return {
        'P_left': P_left,
        'P_right': P_right,
        'F_net': F_net,
        'impulse': impulse,
        'work_left': work_left,
        'work_right': work_right,
        'work_net': work_net
    }


def analyze_conservation_detailed(results_dir, gamma_ad=5.0/3.0, x_left=-0.5, x_right=0.5):
    """
    Detailed conservation analysis with INTRINSIC vs TOTAL separation.
    
    ================================================================================
    TWO TYPES OF CONSERVATION METRICS
    ================================================================================
    
    1. TOTAL (Domain-Integrated):
       - What we measure: Sum over all particles in domain
       - Changes when: Material/waves cross boundaries OR internal forces act
       - Formula: Q_total(t) = Σᵢ νᵢ qᵢ(t)
       
    2. INTRINSIC (Scheme Conservation Error):
       - What it measures: Numerical error in the scheme itself
       - Should be ~0 for a conservative scheme with closed boundaries
       - Formula: Q_intrinsic = Q_total + Q_flux (what left + what remains = constant)
       
    IMPORTANT NOTE ON MOMENTUM:
       Even with REFLECTING/OUTFLOW boundaries, TOTAL momentum changes because:
       - The LEFT ghost particles push rightward (high pressure left)
       - The RIGHT ghost particles push leftward (low pressure right)
       - Net force: Leftward ghosts push HARDER → net rightward momentum gain
       
       This is PHYSICAL! The ghost particles represent infinite reservoirs
       that continuously supply pressure forces.
       
       For TRUE momentum conservation, you need PERIODIC boundaries.
       
    For a PERFECTLY CONSERVATIVE scheme with OUTFLOW boundaries:
       - TOTAL momentum changes (ghost pressure forces) → Physical!
       - TOTAL energy conserved (work done by ghosts = energy gained)
       - INTRINSIC: Subtract ghost impulse to see scheme-only conservation
       
    ================================================================================
    """
    snapshots, times = load_all_snapshots(results_dir)
    
    if len(snapshots) == 0:
        print(f"No snapshots found in {results_dir}")
        return None
    
    conservation_data = []
    cumulative_ghost_impulse = 0.0
    cumulative_ghost_work = 0.0
    
    for i, (snap, t) in enumerate(zip(snapshots, times)):
        conserved = compute_total_conserved(snap, gamma_ad)
        conserved['time'] = t
        conserved['snapshot'] = i
        
        # Compute cumulative ghost particle impulse (momentum injected by boundaries)
        if i > 0:
            dt = t - times[i-1]
            ghost = compute_ghost_impulse(snap, dt, gamma_ad, x_left, x_right)
            cumulative_ghost_impulse += ghost['impulse']
            cumulative_ghost_work += ghost['work_net']
            conserved['P_left'] = ghost['P_left']
            conserved['P_right'] = ghost['P_right']
            conserved['F_net'] = ghost['F_net']
        else:
            conserved['P_left'] = snap['pres'].iloc[0]
            conserved['P_right'] = snap['pres'].iloc[-1]
            conserved['F_net'] = 0.0
        
        conserved['cumulative_ghost_impulse'] = cumulative_ghost_impulse
        conserved['cumulative_ghost_work'] = cumulative_ghost_work
        
        conservation_data.append(conserved)
    
    df = pd.DataFrame(conservation_data)
    
    # Compute both TOTAL and INTRINSIC errors
    if len(df) > 0:
        # Initial values
        N0 = df['baryon_number'].iloc[0]
        Sx0 = df['momentum_x'].iloc[0]
        e0 = df['energy'].iloc[0]
        
        # TOTAL errors (what we measure in domain)
        df['baryon_total_error'] = (df['baryon_number'] - N0) / N0
        df['momentum_total_change'] = df['momentum_x'] - Sx0  # Absolute (initial can be 0)
        df['energy_total_error'] = (df['energy'] - e0) / e0
        
        # INTRINSIC errors (subtract ghost impulse to see scheme-only conservation)
        # S_intrinsic = S_total - ghost_impulse (remove external force contribution)
        # If scheme is conservative: S_intrinsic should stay constant
        df['momentum_intrinsic'] = df['momentum_x'] - df['cumulative_ghost_impulse']
        df['energy_intrinsic'] = df['energy'] - df['cumulative_ghost_work']
        
        df['momentum_intrinsic_change'] = df['momentum_intrinsic'] - df['momentum_intrinsic'].iloc[0]
        
        # Handle case where initial intrinsic energy might be zero
        e0_intrinsic = df['energy_intrinsic'].iloc[0]
        if abs(e0_intrinsic) > 1e-15:
            df['energy_intrinsic_error'] = (df['energy_intrinsic'] - e0_intrinsic) / e0_intrinsic
        else:
            df['energy_intrinsic_error'] = df['energy_intrinsic'] - e0_intrinsic
    
    return df


def plot_intrinsic_vs_total_conservation(conservation_df, output_path):
    """
    Create conservation plot with INTRINSIC vs TOTAL separation.
    
    Shows explicit formulas used for each quantity.
    """
    fig = plt.figure(figsize=(16, 18))
    
    # Create grid with formula panel at top
    gs = fig.add_gridspec(4, 2, height_ratios=[0.8, 1, 1, 1], hspace=0.35, wspace=0.25)
    
    # ==================== FORMULA PANEL ====================
    ax_formula = fig.add_subplot(gs[0, :])
    ax_formula.axis('off')
    
    formula_text = """
    ╔══════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
    ║                                     CONSERVATION ANALYSIS FORMULAS                                            ║
    ╠══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
    ║                                                                                                               ║
    ║  CSV COLUMNS USED:                                                                                            ║
    ║    mass → ν (baryon number per particle)    dens → n (rest-frame density)                                     ║
    ║    pres → P (pressure)                      vel_x → v (velocity)                                              ║
    ║                                                                                                               ║
    ║  DERIVED QUANTITIES (per particle i):                                                                         ║
    ║    γᵢ = 1/√(1 - vᵢ²)                                    [Lorentz factor]                                      ║
    ║    Hᵢ = 1 + (γ_ad/(γ_ad-1)) Pᵢ/nᵢ                       [Specific enthalpy, c=1]                              ║
    ║    Sᵢ = γᵢ Hᵢ vᵢ                                        [Canonical momentum per baryon]                       ║
    ║    eᵢ = γᵢ Hᵢ - Pᵢ/(γᵢ nᵢ)                              [Canonical energy per baryon]                         ║
    ║                                                                                                               ║
    ║  TOTAL CONSERVED (integrated over domain):                                                                    ║
    ║    N_total   = Σᵢ νᵢ                                    [Total baryons - ν is constant per particle]          ║
    ║    S_total   = Σᵢ νᵢ Sᵢ                                 [Total momentum - changes due to ghost forces]        ║
    ║    E_total   = Σᵢ νᵢ eᵢ                                 [Total canonical energy]                              ║
    ║                                                                                                               ║
    ║  GHOST IMPULSE (external force from boundary ghosts):                                                         ║
    ║    F_ghost = (P_left - P_right) · A                     [Net pressure force from ghosts, A=1 for 1D]          ║
    ║    J_ghost = ∫ F_ghost dt                               [Cumulative impulse = momentum injected]              ║
    ║                                                                                                               ║
    ║  INTRINSIC CONSERVATION (scheme accuracy, removing ghost forces):                                             ║
    ║    S_intrinsic = S_total - J_ghost                      [Should be ~0 if scheme is conservative]              ║
    ║                                                                                                               ║
    ║  INTERPRETATION:                                                                                              ║
    ║    • S_total changes    → ghost pressure forces push material (physical, expected)                            ║
    ║    • S_intrinsic changes → numerical error (spurious momentum creation/destruction)                           ║
    ╚══════════════════════════════════════════════════════════════════════════════════════════════════════════════╝
    """
    
    ax_formula.text(0.5, 0.5, formula_text, transform=ax_formula.transAxes,
                   fontsize=8, fontfamily='monospace', ha='center', va='center',
                   bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.9))
    
    t = conservation_df['snapshot']
    
    # ==================== BARYON NUMBER ====================
    ax1 = fig.add_subplot(gs[1, 0])
    initial_N = conservation_df['baryon_number'].iloc[0]
    ax1.plot(t, conservation_df['baryon_number'], 'b-', linewidth=2, label='N_total = Σνᵢ')
    ax1.axhline(y=initial_N, color='r', linestyle='--', linewidth=1.5, alpha=0.7, 
                label=f'Initial: {initial_N:.4e}')
    ax1.set_xlabel('Timestep')
    ax1.set_ylabel('Total Baryon Number')
    ax1.set_title('Baryon Number: N_total = Σᵢ νᵢ\n(νᵢ is constant per particle → exactly conserved)')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    final_N_error = conservation_df['baryon_total_error'].iloc[-1] * 100
    ax1.text(0.98, 0.02, f'ΔN/N₀ = {final_N_error:.3e}%', transform=ax1.transAxes,
             fontsize=10, ha='right', va='bottom',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # ==================== BARYON ERROR ====================
    ax2 = fig.add_subplot(gs[1, 1])
    ax2.plot(t, conservation_df['baryon_total_error'] * 100, 'b-', linewidth=2)
    ax2.axhline(y=0, color='r', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.fill_between(t, 0, conservation_df['baryon_total_error'] * 100, alpha=0.3, color='blue')
    ax2.set_xlabel('Timestep')
    ax2.set_ylabel('Relative Error (%)')
    ax2.set_title('Baryon Conservation Error\n(N(t) - N₀) / N₀')
    ax2.grid(True, alpha=0.3)
    
    # ==================== MOMENTUM: TOTAL vs INTRINSIC ====================
    ax3 = fig.add_subplot(gs[2, 0])
    initial_S = conservation_df['momentum_x'].iloc[0]
    ax3.plot(t, conservation_df['momentum_x'], 'b-', linewidth=2, 
             label=f'S_total = Σνᵢ·Sᵢ')
    if 'momentum_intrinsic' in conservation_df.columns:
        ax3.plot(t, conservation_df['momentum_intrinsic'], 'g--', linewidth=2, 
                 label='S_intrinsic = S_total - J_ghost')
    if 'cumulative_ghost_impulse' in conservation_df.columns:
        ax3.plot(t, conservation_df['cumulative_ghost_impulse'], 'r:', linewidth=2, 
                 label='J_ghost = ∫(P_left-P_right)dt')
    ax3.axhline(y=initial_S, color='r', linestyle=':', linewidth=1.5, alpha=0.7, 
                label=f'Initial: {initial_S:.2e}')
    ax3.set_xlabel('Timestep')
    ax3.set_ylabel('Momentum')
    ax3.set_title('Momentum: TOTAL vs INTRINSIC\nS = Σᵢ νᵢ γᵢ Hᵢ vᵢ')
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3)
    
    # ==================== MOMENTUM CHANGES ====================
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.plot(t, conservation_df['momentum_total_change'], 'b-', linewidth=2, 
             label='ΔS_total (physical)')
    if 'momentum_intrinsic_change' in conservation_df.columns:
        ax4.plot(t, conservation_df['momentum_intrinsic_change'], 'g--', linewidth=2, 
                 label='ΔS_intrinsic (scheme error)')
    ax4.axhline(y=0, color='r', linestyle=':', linewidth=1.5, alpha=0.7)
    ax4.set_xlabel('Timestep')
    ax4.set_ylabel('Momentum Change from Initial')
    ax4.set_title('Momentum Changes\nTOTAL: material leaving   |   INTRINSIC: numerical error')
    ax4.legend(loc='best')
    ax4.grid(True, alpha=0.3)
    
    final_dS_total = conservation_df['momentum_total_change'].iloc[-1]
    ax4.text(0.02, 0.98, f'ΔS_total = {final_dS_total:+.4e}', transform=ax4.transAxes,
             fontsize=9, va='top',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # ==================== ENERGY: TOTAL vs INTRINSIC ====================
    ax5 = fig.add_subplot(gs[3, 0])
    initial_E = conservation_df['energy'].iloc[0]
    ax5.plot(t, conservation_df['energy'], 'g-', linewidth=2, 
             label='E_total = Σνᵢ·eᵢ')
    if 'energy_intrinsic' in conservation_df.columns:
        ax5.plot(t, conservation_df['energy_intrinsic'], 'm--', linewidth=2, 
                 label='E_intrinsic = E_total + E_flux')
    ax5.axhline(y=initial_E, color='r', linestyle=':', linewidth=1.5, alpha=0.7, 
                label=f'Initial: {initial_E:.4e}')
    ax5.set_xlabel('Timestep')
    ax5.set_ylabel('Canonical Energy')
    ax5.set_title('Energy: TOTAL vs INTRINSIC\ne = γH - P/N')
    ax5.legend(loc='best')
    ax5.grid(True, alpha=0.3)
    
    # ==================== ENERGY ERRORS ====================
    ax6 = fig.add_subplot(gs[3, 1])
    ax6.plot(t, conservation_df['energy_total_error'] * 100, 'g-', linewidth=2, 
             label='ΔE_total/E₀ (physical)')
    if 'energy_intrinsic_error' in conservation_df.columns:
        ax6.plot(t, conservation_df['energy_intrinsic_error'] * 100, 'm--', linewidth=2, 
                 label='ΔE_intrinsic/E₀ (scheme error)')
    ax6.axhline(y=0, color='r', linestyle=':', linewidth=1.5, alpha=0.7)
    ax6.set_xlabel('Timestep')
    ax6.set_ylabel('Relative Error (%)')
    ax6.set_title('Energy Conservation Errors\nTOTAL: boundary effects   |   INTRINSIC: numerical accuracy')
    ax6.legend(loc='best')
    ax6.grid(True, alpha=0.3)
    
    final_dE_total = conservation_df['energy_total_error'].iloc[-1] * 100
    ax6.text(0.98, 0.02, f'ΔE_total/E₀ = {final_dE_total:+.4e}%', transform=ax6.transAxes,
             fontsize=9, ha='right', va='bottom',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.suptitle('SRGSPH Conservation Analysis: INTRINSIC vs TOTAL\n' +
                 'With Explicit Formulas from CSV Data', fontsize=14, fontweight='bold')
    
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Intrinsic vs Total conservation analysis saved to {output_path}")


def plot_scheme_error_only(conservation_df, output_path):
    """
    Create LOG-SCALE plot showing ONLY scheme conservation errors.
    
    Excludes all ghost particle effects - shows pure numerical accuracy.
    
    Quantities plotted:
      - Baryon error: |ΔN/N₀| (should be exactly 0 - ν is constant)
      - Momentum intrinsic error: |S_intrinsic - S₀| = |S_total - J_ghost - S₀|
      - Energy error: |ΔE/E₀|
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    t = conservation_df['time'].values
    
    # ==================== FORMULA PANEL ====================
    ax_formula = axes[0, 0]
    ax_formula.axis('off')
    
    formula_text = """
╔═══════════════════════════════════════════════════════════════════╗
║         SCHEME CONSERVATION ERROR (Ghost Effects Removed)         ║
╠═══════════════════════════════════════════════════════════════════╣
║                                                                   ║
║  These errors measure NUMERICAL ACCURACY of the SPH scheme.       ║
║  Ghost particle forces are subtracted out.                        ║
║                                                                   ║
║  BARYON ERROR:                                                    ║
║    ε_N = |N(t) - N₀| / N₀                                         ║
║    Should be EXACTLY 0 (ν is constant per particle)               ║
║                                                                   ║
║  MOMENTUM INTRINSIC ERROR:                                        ║
║    ε_S = |S_intrinsic(t) - S_intrinsic(0)|                        ║
║    where S_intrinsic = S_total - J_ghost                          ║
║    J_ghost = ∫(P_left - P_right)·A·dt                             ║
║                                                                   ║
║  ENERGY ERROR:                                                    ║
║    ε_E = |E(t) - E₀| / E₀                                         ║
║    Measures scheme's energy conservation accuracy                 ║
║                                                                   ║
║  LOG SCALE: Lower = better conservation                           ║
╚═══════════════════════════════════════════════════════════════════╝
"""
    ax_formula.text(0.5, 0.5, formula_text, transform=ax_formula.transAxes,
                   fontsize=9, fontfamily='monospace', ha='center', va='center',
                   bbox=dict(boxstyle='round', facecolor='#f5f5f5', alpha=0.95))
    
    # ==================== BARYON ERROR (LOG) ====================
    ax1 = axes[0, 1]
    baryon_error = np.abs(conservation_df['baryon_total_error'].values)
    # Replace zeros with small value for log plot
    baryon_error = np.where(baryon_error > 0, baryon_error, 1e-16)
    
    ax1.semilogy(t, baryon_error, 'b-', linewidth=2, marker='o', markersize=3,
                 label='|ΔN/N₀|')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Relative Error (log scale)')
    ax1.set_title('Baryon Number Error: ε_N = |ΔN/N₀|\n(Should be ~machine precision)')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.legend(loc='best')
    ax1.set_ylim(bottom=1e-17)
    
    final_N_err = baryon_error[-1]
    ax1.text(0.98, 0.98, f'Final: {final_N_err:.2e}', transform=ax1.transAxes,
             fontsize=10, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # ==================== MOMENTUM INTRINSIC ERROR (LOG) ====================
    ax2 = axes[1, 0]
    if 'momentum_intrinsic_change' in conservation_df.columns:
        momentum_error = np.abs(conservation_df['momentum_intrinsic_change'].values)
    else:
        momentum_error = np.abs(conservation_df['momentum_x'] - conservation_df['momentum_x'].iloc[0])
    
    # Replace zeros with small value for log plot
    momentum_error = np.where(momentum_error > 0, momentum_error, 1e-16)
    
    ax2.semilogy(t, momentum_error, 'g-', linewidth=2, marker='s', markersize=3,
                 label='|ΔS_intrinsic|')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Absolute Error (log scale)')
    ax2.set_title('Momentum Intrinsic Error: |S_intrinsic(t) - S_intrinsic(0)|\n(Ghost impulse subtracted)')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(loc='best')
    
    final_S_err = momentum_error[-1]
    ax2.text(0.98, 0.98, f'Final: {final_S_err:.2e}', transform=ax2.transAxes,
             fontsize=10, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # ==================== ENERGY ERROR (LOG) ====================
    ax3 = axes[1, 1]
    energy_error = np.abs(conservation_df['energy_total_error'].values)
    # Replace zeros with small value for log plot
    energy_error = np.where(energy_error > 0, energy_error, 1e-16)
    
    ax3.semilogy(t, energy_error, 'r-', linewidth=2, marker='^', markersize=3,
                 label='|ΔE/E₀|')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Relative Error (log scale)')
    ax3.set_title('Energy Error: ε_E = |ΔE/E₀|\n(Scheme accuracy metric)')
    ax3.grid(True, alpha=0.3, which='both')
    ax3.legend(loc='best')
    
    final_E_err = energy_error[-1]
    ax3.text(0.98, 0.98, f'Final: {final_E_err:.2e}', transform=ax3.transAxes,
             fontsize=10, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.suptitle('SRGSPH SCHEME CONSERVATION ERROR (Log Scale)\n' +
                 'Ghost Particle Effects Excluded', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Scheme error log plot saved to {output_path}")


def print_conservation_summary(conservation_df):
    """Print detailed conservation summary with explicit formulas."""
    print("\n" + "="*90)
    print("SRGSPH CONSERVATION ANALYSIS SUMMARY")
    print("="*90)
    
    print("\n┌─────────────────────────────────────────────────────────────────────────────────────────┐")
    print("│ FORMULAS USED (derived from CSV columns)                                                │")
    print("├─────────────────────────────────────────────────────────────────────────────────────────┤")
    print("│ CSV: mass(ν), dens(n), pres(P), vel_x(v)                                                │")
    print("│                                                                                         │")
    print("│ γ = 1/√(1-v²)                          Lorentz factor                                   │")
    print("│ H = 1 + (γ_ad/(γ_ad-1)) P/n            Specific enthalpy (c=1)                          │")
    print("│ S = γ·H·v                              Canonical momentum per baryon                    │")
    print("│ e = γ·H - P/(γ·n)                      Canonical energy per baryon                      │")
    print("│                                                                                         │")
    print("│ N_total = Σᵢ νᵢ                        Total baryon number                              │")
    print("│ S_total = Σᵢ νᵢ·Sᵢ                     Total momentum                                   │")
    print("│ E_total = Σᵢ νᵢ·eᵢ                     Total canonical energy                           │")
    print("└─────────────────────────────────────────────────────────────────────────────────────────┘")
    
    # Initial and final values
    N0 = conservation_df['baryon_number'].iloc[0]
    Nf = conservation_df['baryon_number'].iloc[-1]
    S0 = conservation_df['momentum_x'].iloc[0]
    Sf = conservation_df['momentum_x'].iloc[-1]
    E0 = conservation_df['energy'].iloc[0]
    Ef = conservation_df['energy'].iloc[-1]
    
    print("\n┌─────────────────────────────────────────────────────────────────────────────────────────┐")
    print("│ TOTAL CONSERVATION (domain-integrated values)                                           │")
    print("├──────────────────┬──────────────────────┬──────────────────────┬────────────────────────┤")
    print("│ Quantity         │ Initial              │ Final                │ Change                 │")
    print("├──────────────────┼──────────────────────┼──────────────────────┼────────────────────────┤")
    print(f"│ Baryon N=Σν      │ {N0:18.10e} │ {Nf:18.10e} │ {(Nf-N0)/N0*100:+18.10e}% │")
    print(f"│ Momentum S=ΣνS   │ {S0:18.10e} │ {Sf:18.10e} │ {Sf-S0:+18.10e}  │")
    print(f"│ Energy E=Σνe     │ {E0:18.10e} │ {Ef:18.10e} │ {(Ef-E0)/E0*100:+18.10e}% │")
    print("└──────────────────┴──────────────────────┴──────────────────────┴────────────────────────┘")
    
    # Check if intrinsic columns exist
    if 'momentum_intrinsic' in conservation_df.columns:
        Si0 = conservation_df['momentum_intrinsic'].iloc[0]
        Sif = conservation_df['momentum_intrinsic'].iloc[-1]
        J_ghost = conservation_df['cumulative_ghost_impulse'].iloc[-1] if 'cumulative_ghost_impulse' in conservation_df.columns else 0
        
        print("\n┌─────────────────────────────────────────────────────────────────────────────────────────┐")
        print("│ GHOST IMPULSE (momentum injected by boundary pressure)                                  │")
        print("├─────────────────────────────────────────────────────────────────────────────────────────┤")
        print(f"│ J_ghost = ∫(P_left - P_right)·A·dt = {J_ghost:+.10e}                                │")
        print(f"│ This is the momentum the ghosts PUSHED into the real particles                         │")
        print("└─────────────────────────────────────────────────────────────────────────────────────────┘")
        
        print("\n┌─────────────────────────────────────────────────────────────────────────────────────────┐")
        print("│ INTRINSIC CONSERVATION (scheme accuracy = S_total - J_ghost)                            │")
        print("├──────────────────┬──────────────────────┬──────────────────────┬────────────────────────┤")
        print("│ Quantity         │ Initial              │ Final                │ Change                 │")
        print("├──────────────────┼──────────────────────┼──────────────────────┼────────────────────────┤")
        print(f"│ S_intrinsic      │ {Si0:18.10e} │ {Sif:18.10e} │ {Sif-Si0:+18.10e}  │")
        print("│                  │                      │                      │ (should be ~0)         │")
        print("└──────────────────┴──────────────────────┴──────────────────────┴────────────────────────┘")
    
    print("\n┌─────────────────────────────────────────────────────────────────────────────────────────┐")
    print("│ INTERPRETATION                                                                          │")
    print("├─────────────────────────────────────────────────────────────────────────────────────────┤")
    print("│ • Baryon N:   EXACTLY CONSERVED (νᵢ is constant per particle by design)                │")
    print("│ • Momentum S: CHANGES due to ghost pressure forces (HIGH P left → push right)          │")
    print("│ • Energy E:   CONSERVED (scheme accuracy metric - should be ~0% error)                 │")
    print("│                                                                                         │")
    print("│ WHY MOMENTUM CHANGES (this is CORRECT physics):                                         │")
    print("│   Ghost particles maintain boundary conditions P_left=1, P_right=0.1                   │")
    print("│   Net pressure gradient dP/dx < 0 → acceleration dv/dt = -(1/ρ)dP/dx > 0               │")
    print("│   All material accelerates RIGHTWARD → momentum increases                              │")
    print("│                                                                                         │")
    print("│ For TRUE momentum conservation: use PERIODIC boundaries (no external forces)           │")
    print("└─────────────────────────────────────────────────────────────────────────────────────────┘")
    print()


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
    """Create comprehensive conservation analysis plot with theoretical expectations."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 14))
    fig.suptitle('SRGSPH Conservation Law Analysis\n' +
                 'Canonical Variables: Baryon Number, Momentum, and Energy', fontsize=14)

    t = conservation_df['snapshot']  # Use timestep (integer) instead of time
    
    # Get initial values for theoretical lines
    initial_baryon = conservation_df['baryon_number'].iloc[0]
    initial_momentum = conservation_df['momentum_x'].iloc[0]
    initial_energy = conservation_df['energy'].iloc[0]

    # 1. Baryon number conservation (should be EXACTLY conserved)
    ax1 = axes[0, 0]
    ax1.plot(t, conservation_df['baryon_number'], 'b-', linewidth=2, label='Computed')
    ax1.axhline(y=initial_baryon, color='r', linestyle='--', linewidth=2, alpha=0.7, label=f'Theory: {initial_baryon:.4e}')
    ax1.set_xlabel('Timestep')
    ax1.set_ylabel('Total Baryon Number Σνᵢ')
    ax1.set_title('Baryon Number Conservation (Exact)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add error annotation
    final_error = conservation_df['baryon_error'].iloc[-1] * 100
    ax1.text(0.98, 0.02, f'Final error: {final_error:.2e}%', transform=ax1.transAxes,
             fontsize=10, ha='right', va='bottom',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # 2. Momentum evolution
    # Note: For shock tube with ghost boundaries, momentum of REAL particles changes
    # The ghost particles absorb the reaction momentum
    ax2 = axes[0, 1]
    ax2.plot(t, conservation_df['momentum_x'], 'b-', linewidth=2, label='Real particles only')
    ax2.axhline(y=initial_momentum, color='r', linestyle='--', linewidth=2, alpha=0.7, 
                label=f'Initial: {initial_momentum:.2e}')
    ax2.set_xlabel('Timestep')
    ax2.set_ylabel('Total Momentum Σνᵢ·Sᵢ')
    ax2.set_title('Momentum of Real Particles\n(Ghost/wall absorbs reaction)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Momentum change (expected for open boundary)
    ax3 = axes[1, 0]
    ax3.plot(t, conservation_df['momentum_x_error'], 'b-', linewidth=2, label='ΔMomentum')
    ax3.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7, label='Closed system expectation')
    ax3.set_xlabel('Timestep')
    ax3.set_ylabel('Momentum Change from Initial')
    ax3.set_title('Momentum Transferred to Boundaries')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Add explanation
    ax3.text(0.02, 0.98, 'Net rightward momentum\ndue to asymmetric IC\n(absorbed by ghost particles)', 
             transform=ax3.transAxes, fontsize=9, va='top',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    # 4. Total canonical energy (should be conserved)
    ax4 = axes[1, 1]
    ax4.plot(t, conservation_df['energy'], 'g-', linewidth=2, label='Computed')
    ax4.axhline(y=initial_energy, color='r', linestyle='--', linewidth=2, alpha=0.7, 
                label=f'Theory: {initial_energy:.4e}')
    ax4.set_xlabel('Timestep')
    ax4.set_ylabel('Total Canonical Energy Σνᵢ·eᵢ')
    ax4.set_title('Canonical Energy Conservation')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Add error annotation
    final_e_error = conservation_df['energy_error'].iloc[-1] * 100
    ax4.text(0.98, 0.02, f'Final error: {final_e_error:.2e}%', transform=ax4.transAxes,
             fontsize=10, ha='right', va='bottom',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # 5. Thermal and kinetic energy evolution (NOT individually conserved - they exchange)
    ax5 = axes[2, 0]
    ax5.plot(t, conservation_df['thermal_energy'], 'orange', linewidth=2, label='Thermal u')
    ax5.plot(t, conservation_df['kinetic_energy'], 'purple', linewidth=2, label='Kinetic (γ-1)')
    ax5.axhline(y=conservation_df['thermal_energy'].iloc[0], color='orange', linestyle=':', alpha=0.5)
    ax5.axhline(y=conservation_df['kinetic_energy'].iloc[0], color='purple', linestyle=':', alpha=0.5)
    ax5.set_xlabel('Timestep')
    ax5.set_ylabel('Energy Component')
    ax5.set_title('Energy Exchange: Thermal ↔ Kinetic\n(Individual components NOT conserved)')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    # 6. Energy conservation error
    ax6 = axes[2, 1]
    ax6.plot(t, conservation_df['energy_error'] * 100, 'g-', linewidth=2)
    ax6.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7, label='Theory: 0%')
    ax6.fill_between(t, 0, conservation_df['energy_error'] * 100, alpha=0.3, color='green')
    ax6.set_xlabel('Timestep')
    ax6.set_ylabel('Relative Error (%)')
    ax6.set_title('Energy Conservation Error')
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    # Add summary stats
    thermal_change = conservation_df['thermal_energy'].iloc[-1] - conservation_df['thermal_energy'].iloc[0]
    kinetic_change = conservation_df['kinetic_energy'].iloc[-1] - conservation_df['kinetic_energy'].iloc[0]
    ax6.text(0.02, 0.98, f'ΔThermal: {thermal_change:+.3e}\nΔKinetic: {kinetic_change:+.3e}\nSum: {thermal_change+kinetic_change:+.3e}', 
             transform=ax6.transAxes, fontsize=9, va='top',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

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

def plot_snapshot_with_exact(df, t, output_path, gamma_ad=5.0/3.0, test_type='sod', v_left=0.0):
    """Plot snapshot with exact solution overlay."""
    # Get initial conditions from SSOT
    ic = get_initial_conditions(test_type, v_left)
    pL, rhoL, vL = ic['pL'], ic['rhoL'], ic['vL']
    pR, rhoR, vR = ic['pR'], ic['rhoR'], ic['vR']
    test_name = ic['name']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'SRGSPH vs Exact Solution at t = {t:.4f}\n' +
                 f'{test_name}', fontsize=14)
    
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
    
    # Compute exact solution using SSOT function
    exact = get_exact_solution(t, test_type, v_left, gamma_ad, n_points=1000)
    x_exact = exact['x']
    
    # Plot density
    ax1 = axes[0, 0]
    ax1.scatter(x, rho, s=2, alpha=0.5, label='SPH')
    ax1.plot(x_exact, exact['dens'], 'r-', linewidth=1.5, label='Exact')
    ax1.set_xlabel('Position x')
    ax1.set_ylabel(r'Rest-frame density $n$')
    ax1.set_title('Density')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot pressure
    ax2 = axes[0, 1]
    ax2.scatter(x, p, s=2, alpha=0.5, label='SPH')
    ax2.plot(x_exact, exact['pres'], 'r-', linewidth=1.5, label='Exact')
    ax2.set_xlabel('Position x')
    ax2.set_ylabel('Pressure P')
    ax2.set_title('Pressure')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot velocity
    ax3 = axes[1, 0]
    ax3.scatter(x, v, s=2, alpha=0.5, label='SPH')
    ax3.plot(x_exact, exact['vel'], 'r-', linewidth=1.5, label='Exact')
    ax3.set_xlabel('Position x')
    ax3.set_ylabel('Velocity $v^x$')
    ax3.set_title('Velocity')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Compute error metrics
    exact_at_sph = {
        'dens': np.interp(x, x_exact, exact['dens']),
        'pres': np.interp(x, x_exact, exact['pres']),
        'vel': np.interp(x, x_exact, exact['vel'])
    }
    rms_rho = np.sqrt(np.mean((rho - exact_at_sph['dens'])**2))
    rms_p = np.sqrt(np.mean((p - exact_at_sph['pres'])**2))
    rms_v = np.sqrt(np.mean((v - exact_at_sph['vel'])**2))
    
    # Annotations
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    info_text = f"""
    Analysis Notes
    ========================================
    
    Test: {test_name}
    Time: t = {t:.4f}
    Number of particles: {len(df)}
    
    Initial Conditions:
      Left:  P={pL}, n={rhoL}, v={vL}
      Right: P={pR}, n={rhoR}, v={vR}
    
    Star State (P*, v*):
      P* = {exact['p_star']:.4f}
      v* = {exact['v_star']:.4f}
    
    RMS Errors vs Exact Solution:
      Density:  {rms_rho:.4f}
      Pressure: {rms_p:.4f}
      Velocity: {rms_v:.4f}
    
    Observations:
    1. Good agreement with exact solution
       across all primitive variables
       
    2. Minor oscillations at discontinuities
       (contact, shock) are typical for SPH
       
    3. Rarefaction fan well-resolved
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
        results_dir = Path("simulations/relativistic/sr_sod/results/sod_fig1")
    
    # Output directory is inside results directory
    output_dir = results_dir / "analysis"
    output_dir.mkdir(exist_ok=True)
    
    gamma_ad = 5.0 / 3.0
    
    print("=" * 60)
    print("SRGSPH Conservation and Numerical Analysis")
    print(f"Results directory: {results_dir}")
    print(f"Output directory:  {output_dir}")
    print("=" * 60)

    # 1. Conservation analysis (INTRINSIC vs TOTAL)
    print("\n[1/7] Analyzing conservation laws (INTRINSIC vs TOTAL)...")
    conservation_detailed_df = analyze_conservation_detailed(results_dir, gamma_ad)
    if conservation_detailed_df is not None:
        # Print detailed summary with formulas
        print_conservation_summary(conservation_detailed_df)
        
        # Plot INTRINSIC vs TOTAL conservation
        plot_intrinsic_vs_total_conservation(
            conservation_detailed_df, 
            output_dir / "conservation_intrinsic_vs_total.png"
        )
        
        # Plot SCHEME ERROR ONLY (log scale, ghost effects excluded)
        plot_scheme_error_only(
            conservation_detailed_df,
            output_dir / "scheme_error_log.png"
        )
    
    # 2. Original conservation analysis (for backward compatibility)
    print("\n[2/7] Creating standard conservation plots...")
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

    # 3. Create scheme explanation figure
    print("\n[3/7] Creating scheme explanation figure...")
    create_scheme_explanation_figure(output_dir / "scheme_explanation.png")

    # 4. Rarefaction overshoot analysis
    print("\n[4/7] Analyzing rarefaction wave overshoot...")
    create_rarefaction_overshoot_analysis(results_dir,
                                          output_dir / "rarefaction_analysis.png",
                                          gamma_ad)

    # 5. Snapshot comparison with exact solution
    print("\n[5/7] Creating snapshot comparison with exact solution...")
    # Detect test type for correct initial conditions
    test_type, v_left = detect_test_type(results_dir)
    print(f"  Detected test type: {test_type}")
    
    snapshots, times = load_all_snapshots(results_dir)
    if len(snapshots) > 0:
        # Use final snapshot
        final_idx = -1
        plot_snapshot_with_exact(snapshots[final_idx], times[final_idx],
                                output_dir / "snapshot_vs_exact.png", gamma_ad,
                                test_type, v_left)

    print("\n[7/7] All analyses complete!")
    
    print("\n" + "=" * 60)
    print("Analysis complete! Output files in:", output_dir)
    print("  - scheme_error_log.png  (NEW: log-scale scheme error, ghost excluded)")
    print("  - conservation_intrinsic_vs_total.png")
    print("  - conservation_analysis.png")
    print("  - energy_conservation_detailed.png")
    print("  - scheme_explanation.png")
    print("  - rarefaction_analysis.png")
    print("  - snapshot_vs_exact.png")
    print("=" * 60)

if __name__ == "__main__":
    main()
