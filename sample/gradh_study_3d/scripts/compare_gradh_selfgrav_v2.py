#!/usr/bin/env python3
"""
Compare Self-Gravitating Cloud Evolution: Grad-H vs No Grad-H
=============================================================
Visualizes the effect of grad-h correction on preventing artificial core collapse.
Includes: density plots, radial profiles, momentum conservation, energy conservation.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import glob
import os
import base64
import sys

# SSOT: Use shared Lane-Emden module
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.shared.lane_emden import load_lane_emden_solution, get_density_profile


# ============================================================================
# Color Palette (High Contrast)
# ============================================================================
COLOR_WITH_GRADH = '#00FF00'     # Bright green
COLOR_NO_GRADH = '#FF3333'       # Bright red
COLOR_ANALYTIC = '#00FFFF'       # Cyan
COLOR_TITLE = '#FFD700'          # Gold
COLOR_KINETIC = '#FF6B6B'        # Coral
COLOR_POTENTIAL = '#4ECDC4'      # Teal
COLOR_INTERNAL = '#FFE66D'       # Yellow
COLOR_TOTAL = '#FFFFFF'          # White


def setup_vtk_style():
    """Configure matplotlib for VTK-like appearance."""
    plt.rcParams.update({
        'figure.facecolor': '#1a1a1a',
        'axes.facecolor': '#252525',
        'axes.edgecolor': '#888888',
        'axes.labelcolor': '#ffffff',
        'axes.titlecolor': '#ffffff',
        'text.color': '#ffffff',
        'xtick.color': '#ffffff',
        'ytick.color': '#ffffff',
        'grid.color': '#444444',
        'grid.alpha': 0.5,
        'legend.facecolor': '#333333',
        'legend.edgecolor': '#666666',
        'legend.labelcolor': '#ffffff',
        'font.family': 'sans-serif',
        'font.size': 10,
        'text.usetex': False,
        'mathtext.fontset': 'cm',
    })


# ============================================================================
# Physics Functions
# ============================================================================
def load_snapshot(filepath):
    """Load a snapshot CSV file."""
    return pd.read_csv(filepath, comment="#")


def compute_energies(df):
    """Compute energy components: kinetic, potential (gravity), internal."""
    # Kinetic energy: E_k = 0.5 * m * v^2
    v2 = df['vel_x']**2 + df['vel_y']**2
    if 'vel_z' in df.columns:
        v2 += df['vel_z']**2
    E_kinetic = 0.5 * (df['mass'] * v2).sum()
    
    # Gravitational potential energy: E_pot = 0.5 * sum(m * phi)
    # Factor of 0.5 to avoid double counting
    E_potential = 0.5 * (df['mass'] * df['phi']).sum()
    
    # Internal (thermal) energy: E_int = sum(m * u)
    E_internal = (df['mass'] * df['ene']).sum()
    
    # Total energy
    E_total = E_kinetic + E_potential + E_internal
    
    return E_kinetic, E_potential, E_internal, E_total


def compute_total_momentum(df):
    """Compute total momentum vector and magnitude."""
    px = (df['mass'] * df['vel_x']).sum()
    py = (df['mass'] * df['vel_y']).sum()
    pz = (df['mass'] * df['vel_z']).sum() if 'vel_z' in df.columns else 0
    return np.sqrt(px**2 + py**2 + pz**2), px, py, pz


def compute_central_density(df, r_fraction=0.1):
    """Compute central density (mean of inner region)."""
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + (df['pos_z']**2 if 'pos_z' in df.columns else 0))
    r_max = r.max()
    mask = r < r_fraction * r_max
    return df.loc[mask, 'dens'].mean() if mask.sum() > 0 else df['dens'].max()


def compute_radial_profile(df, n_bins=50):
    """Compute radial density profile."""
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 + (df['pos_z']**2 if 'pos_z' in df.columns else 0))
    r_max = r.max()
    bins = np.linspace(0, r_max, n_bins + 1)
    r_centers = 0.5 * (bins[:-1] + bins[1:])
    
    rho_profile = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (r >= bins[i]) & (r < bins[i+1])
        if mask.sum() > 0:
            rho_profile[i] = df.loc[mask, 'dens'].mean()
    
    return r_centers, rho_profile


def lane_emden_profile(rho_c, R_cloud, n_points=500):
    """Get Lane-Emden density profile for n=3/2 polytrope (gamma=5/3).
    
    Uses the shared lane_emden module (SSOT).
    
    For n=1.5, ξ₁ ≈ 3.6538 (first zero of θ).
    """
    solution = load_lane_emden_solution(n=1.5, dim=3)
    return get_density_profile(solution, rho_c, R_cloud, n_points=n_points)


def get_snapshots(results_dir):
    """Get sorted list of snapshot files."""
    return sorted(glob.glob(os.path.join(results_dir, "snapshot_*.csv")))


# ============================================================================
# Main Rendering Function
# ============================================================================
def render_frame(with_file, no_file, output_path, frame_num, time,
                 initial_rho_c, initial_R, rho_max_global,
                 time_history, momentum_with, momentum_no,
                 energy_with, energy_no):
    """
    Render comparison frame with 6 panels:
    - Row 1: Scatter (with), Scatter (no), Radial profile
    - Row 2: Momentum conservation, Energy components, Derivation text
    """
    setup_vtk_style()
    
    fig = plt.figure(figsize=(24, 13))
    fig.patch.set_facecolor('#1a1a1a')
    
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.2, 1.2], height_ratios=[1, 1],
                          hspace=0.25, wspace=0.30,
                          left=0.04, right=0.96, top=0.93, bottom=0.05)
    
    ax1 = fig.add_subplot(gs[0, 0])  # Scatter with
    ax2 = fig.add_subplot(gs[0, 1])  # Scatter no
    ax3 = fig.add_subplot(gs[0, 2])  # Radial profile
    ax4 = fig.add_subplot(gs[1, 0])  # Momentum
    ax5 = fig.add_subplot(gs[1, 1])  # Energy
    ax6 = fig.add_subplot(gs[1, 2])  # Derivation
    
    df_with = load_snapshot(with_file)
    df_no = load_snapshot(no_file)
    
    vmin, vmax = 0.1, rho_max_global
    
    # ========================================================================
    # Panel 1 & 2: Scatter plots
    # ========================================================================
    for ax, df, title, color in [(ax1, df_with, 'WITH Grad-H (Stable)', COLOR_WITH_GRADH),
                                  (ax2, df_no, 'NO Grad-H (Collapsing)', COLOR_NO_GRADH)]:
        sc = ax.scatter(df['pos_x'], df['pos_y'], c=df['dens'], 
                       s=2, cmap='inferno', vmin=vmin, vmax=vmax, alpha=0.9)
        ax.set_title(title, fontsize=12, fontweight='bold', color=color)
        ax.set_xlabel('$x$', fontsize=10)
        ax.set_ylabel('$y$', fontsize=10)
        ax.set_aspect('equal')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.grid(True, alpha=0.3)
        cb = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        cb.set_label(r'$\rho$', fontsize=9)
    
    # ========================================================================
    # Panel 3: Radial density profile (log scale)
    # ========================================================================
    r_with, rho_with = compute_radial_profile(df_with)
    r_no, rho_no = compute_radial_profile(df_no)
    
    # Analytic Lane-Emden
    if initial_rho_c and initial_R:
        r_le, rho_le = lane_emden_profile(initial_rho_c, initial_R)
        ax3.plot(r_le, rho_le, '-', color=COLOR_ANALYTIC, lw=3, label='Lane-Emden', zorder=10)
    
    ax3.plot(r_with, rho_with, '-', color=COLOR_WITH_GRADH, lw=2.5, label='With Grad-H')
    ax3.plot(r_no, rho_no, '--', color=COLOR_NO_GRADH, lw=2.5, label='No Grad-H')
    
    ax3.set_xlabel('Radius $r$', fontsize=10)
    ax3.set_ylabel('Density $\\rho$', fontsize=10)
    ax3.set_title('Radial Density Profile', fontsize=12, fontweight='bold')
    ax3.set_yscale('log')
    ax3.set_xlim(0, 1.2)
    ax3.set_ylim(0.01, rho_max_global * 1.5)
    ax3.legend(fontsize=9, loc='upper right')
    ax3.grid(True, alpha=0.3, which='both')
    
    # ========================================================================
    # Panel 4: Momentum Conservation (both cases)
    # ========================================================================
    if len(time_history) > 1:
        # Normalize to initial momentum
        p0_with = momentum_with[0] if momentum_with[0] > 1e-15 else 1e-15
        p0_no = momentum_no[0] if momentum_no[0] > 1e-15 else 1e-15
        
        ax4.plot(time_history, np.array(momentum_with) / p0_with, '-', 
                    color=COLOR_WITH_GRADH, lw=2.5, label='With Grad-H')
        ax4.plot(time_history, np.array(momentum_no) / p0_no, '--', 
                    color=COLOR_NO_GRADH, lw=2.5, label='No Grad-H')
        
        # Mark current time
        ax4.axvline(time, color='#888888', ls=':', lw=1.5, alpha=0.7)
        
        # Reference line at 1
        ax4.axhline(1.0, color='#666666', ls='-', lw=1, alpha=0.5)
    
    ax4.set_xlabel('Time $t$', fontsize=10)
    ax4.set_ylabel('$|\\vec{p}_{tot}|/|\\vec{p}_0|$', fontsize=10)
    ax4.set_title('Momentum Conservation', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9, loc='best')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_xlim(0, max(time_history) if time_history else 100)
    
    # ========================================================================
    # Panel 5: Energy Components (both cases - side by side style)
    # ========================================================================
    if len(time_history) > 1:
        E_k_with, E_p_with, E_i_with, E_t_with = zip(*energy_with)
        E_k_no, E_p_no, E_i_no, E_t_no = zip(*energy_no)
        
        # Normalize to initial total energy magnitude
        E0_with = abs(E_t_with[0]) if abs(E_t_with[0]) > 1e-15 else 1
        E0_no = abs(E_t_no[0]) if abs(E_t_no[0]) > 1e-15 else 1
        
        # Plot WITH grad-h components (solid lines)
        ax5.plot(time_history, np.array(E_k_with)/E0_with, '-', color='#FF6B6B', lw=1.5, alpha=0.8, label='$E_k$ (with)')
        ax5.plot(time_history, np.array(E_p_with)/E0_with, '-', color='#4ECDC4', lw=1.5, alpha=0.8, label='$E_{grav}$ (with)')
        ax5.plot(time_history, np.array(E_i_with)/E0_with, '-', color='#FFE66D', lw=1.5, alpha=0.8, label='$E_{int}$ (with)')
        ax5.plot(time_history, np.array(E_t_with)/E0_with, '-', color=COLOR_WITH_GRADH, lw=3, label='$E_{tot}$ (with)')
        
        # Plot NO grad-h components (dashed lines)
        ax5.plot(time_history, np.array(E_k_no)/E0_no, '--', color='#FF6B6B', lw=1.5, alpha=0.8, label='$E_k$ (no)')
        ax5.plot(time_history, np.array(E_p_no)/E0_no, '--', color='#4ECDC4', lw=1.5, alpha=0.8, label='$E_{grav}$ (no)')
        ax5.plot(time_history, np.array(E_i_no)/E0_no, '--', color='#FFE66D', lw=1.5, alpha=0.8, label='$E_{int}$ (no)')
        ax5.plot(time_history, np.array(E_t_no)/E0_no, '--', color=COLOR_NO_GRADH, lw=3, label='$E_{tot}$ (no)')
        
        ax5.axvline(time, color='#888888', ls=':', lw=1.5, alpha=0.7)
    
    ax5.set_xlabel('Time $t$', fontsize=10)
    ax5.set_ylabel('Energy / $|E_0|$', fontsize=10)
    ax5.set_title('Energy Components: solid=with, dashed=no grad-h', fontsize=11, fontweight='bold')
    ax5.legend(fontsize=7, loc='center left', bbox_to_anchor=(1.01, 0.5), ncol=1)
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim(0, max(time_history) if time_history else 100)
    
    # ========================================================================
    # Panel 6: Derivations (Ω grad-h + Hernquist-Katz softening)
    # ========================================================================
    ax6.axis('off')
    ax6.set_facecolor('#1a1a1a')
    
    y = 0.99
    dy = 0.038  # Reduced spacing to fit more content
    
    # ---- PART A: Grad-h (Ω) Derivation ----
    ax6.text(0.5, y, r'$\mathbf{A.\ Grad}$-$\mathbf{h\ Correction\ (\Omega)\ Derivation}$', 
             transform=ax6.transAxes, fontsize=10, ha='center', color=COLOR_TITLE, fontweight='bold')
    y -= 0.035
    
    ax6.text(0.02, y, r'SPH Lagrangian: $L = \sum_i m_i \left[\frac{1}{2}|\dot{\mathbf{x}}_i|^2 - u_i(\rho_i)\right]$',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy
    ax6.text(0.02, y, r'Density: $\rho_i = \sum_j m_j W_{ij}(h_i)$ with $h_i = \eta(m_i/\rho_i)^{1/D}$ (implicit!)',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy
    ax6.text(0.02, y, r'Chain rule: $\frac{\partial \rho_k}{\partial \mathbf{x}_i} = \underbrace{\sum_j m_j \nabla_i W_{kj}}_{\text{direct}} + \underbrace{\sum_j m_j \frac{\partial W_{kj}}{\partial h_k}\frac{\partial h_k}{\partial \mathbf{x}_i}}_{\text{via } h(\rho)}$',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy
    ax6.text(0.02, y, r'Solving: $\frac{\partial \rho_k}{\partial \mathbf{x}_i} = \Omega_k \sum_j m_j \nabla_i W_{kj}$',
             transform=ax6.transAxes, fontsize=8, ha='left', color=COLOR_WITH_GRADH)
    y -= dy * 0.9
    
    # Omega result box (smaller)
    rect = plt.Rectangle((0.15, y - 0.025), 0.70, 0.035, facecolor='#252525', 
                          edgecolor=COLOR_TITLE, lw=1.5, transform=ax6.transAxes)
    ax6.add_patch(rect)
    ax6.text(0.5, y - 0.008, r'$\Omega_i = \left[1 + \frac{h_i}{D\rho_i}\sum_j m_j \frac{\partial W_{ij}}{\partial h_i}\right]^{-1}$',
             transform=ax6.transAxes, fontsize=9, ha='center', color='#ffffff')
    y -= 0.055
    
    # ---- Separator ----
    ax6.plot([0.05, 0.95], [y+0.01, y+0.01], color='#555555', lw=1, transform=ax6.transAxes)
    y -= 0.015
    
    # ---- PART B: Hernquist-Katz Derivation ----
    ax6.text(0.5, y, r'$\mathbf{B.\ Hernquist}$-$\mathbf{Katz\ (1989)\ Gravitational\ Softening}$', 
             transform=ax6.transAxes, fontsize=10, ha='center', color=COLOR_ANALYTIC, fontweight='bold')
    y -= 0.035
    
    ax6.text(0.02, y, r'Goal: Soften $1/r$ singularity while conserving energy. Define $\epsilon = h/2$ (softening length).',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy
    ax6.text(0.02, y, r'Potential $\phi(r,\epsilon)$ from spline kernel (cubic B-spline), $u = r/\epsilon$:',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy * 0.9
    
    # H-K potential formulas
    ax6.text(0.04, y, r'$\phi = \frac{1}{\epsilon}\left[-\frac{u^2}{2}\left(\frac{1}{3} - \frac{3u^2}{20} + \frac{u^3}{20}\right) + \frac{7}{5}\right]$, $u < 1$',
             transform=ax6.transAxes, fontsize=7, ha='left', color='#aaaaaa')
    y -= dy * 0.85
    ax6.text(0.04, y, r'$\phi = \frac{1}{\epsilon}\left[-\frac{u^2}{2}\left(\frac{4}{3} - u + \frac{3u^2}{10} - \frac{u^3}{30}\right) + \frac{8}{5}\right] - \frac{1}{15r}$, $1 \le u < 2$',
             transform=ax6.transAxes, fontsize=7, ha='left', color='#aaaaaa')
    y -= dy * 0.85
    ax6.text(0.04, y, r'$\phi = 1/r$, $u \ge 2$ (Newtonian)',
             transform=ax6.transAxes, fontsize=7, ha='left', color='#aaaaaa')
    y -= dy
    
    ax6.text(0.02, y, r'Force $g(r) = -d\phi/dr$ (per unit mass, normalized):',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy * 0.9
    ax6.text(0.04, y, r'$g = \frac{1}{\epsilon^3}\left(\frac{4}{3} - \frac{6u^2}{5} + \frac{u^3}{2}\right)$, $u < 1$',
             transform=ax6.transAxes, fontsize=7, ha='left', color='#aaaaaa')
    y -= dy * 0.85
    ax6.text(0.04, y, r'$g = \frac{1}{r^3}\left(-\frac{1}{15} + \frac{8u^3}{3} - 3u^4 + \frac{6u^5}{5} - \frac{u^6}{6}\right)$, $1 \le u < 2$',
             transform=ax6.transAxes, fontsize=7, ha='left', color='#aaaaaa')
    y -= dy
    
    # H-K result box
    rect2 = plt.Rectangle((0.10, y - 0.025), 0.80, 0.035, facecolor='#252525', 
                           edgecolor=COLOR_ANALYTIC, lw=1.5, transform=ax6.transAxes)
    ax6.add_patch(rect2)
    ax6.text(0.5, y - 0.008, r'$\mathbf{F}_{ij} = -G m_i m_j\, g(r_{ij}, \epsilon_{ij})\, \hat{\mathbf{r}}_{ij}$, $\epsilon_{ij} = (h_i + h_j)/4$',
             transform=ax6.transAxes, fontsize=8, ha='center', color='#ffffff')
    y -= 0.05
    
    # ---- Separator ----
    ax6.plot([0.05, 0.95], [y+0.01, y+0.01], color='#555555', lw=1, transform=ax6.transAxes)
    y -= 0.015
    
    # ---- PART C: Why Ω matters for self-gravity ----
    ax6.text(0.5, y, r'$\mathbf{C.\ Why\ \Omega\ Matters\ for\ Self}$-$\mathbf{Gravity}$',
             transform=ax6.transAxes, fontsize=10, ha='center', color=COLOR_TITLE, fontweight='bold')
    y -= 0.03
    ax6.text(0.02, y, r'• H-K softening: $\epsilon = h/2$ $\Rightarrow$ $\phi_i$ depends on $h_i(\rho_i)$',
             transform=ax6.transAxes, fontsize=8, ha='left', color='#ffffff')
    y -= dy * 0.8
    ax6.text(0.02, y, r'• $\Omega=1$: ignores $\partial\phi/\partial h$ $\Rightarrow$ grav. energy NOT conserved $\Rightarrow$ artificial collapse',
             transform=ax6.transAxes, fontsize=8, ha='left', color=COLOR_NO_GRADH)
    y -= dy * 0.8
    ax6.text(0.02, y, r'• $\Omega\neq1$: accounts for $h(\rho)$ in $\phi$ $\Rightarrow$ total energy conserved $\Rightarrow$ stable',
             transform=ax6.transAxes, fontsize=8, ha='left', color=COLOR_WITH_GRADH)
    
    # Main title
    fig.suptitle(r'Grad-H Correction: Lagrangian Consistency $\Rightarrow$ Energy Conservation', 
                 fontsize=16, fontweight='bold', color=COLOR_TITLE)
    
    plt.savefig(output_path, dpi=120, facecolor='#1a1a1a')
    plt.close(fig)
    return output_path


# ============================================================================
# Main Processing
# ============================================================================
def process_all_frames(with_dir, no_dir, output_dir, create_gif=True, create_html=True):
    """Process all snapshots and create visualizations."""
    
    with_files = get_snapshots(with_dir)
    no_files = get_snapshots(no_dir)
    
    if len(with_files) != len(no_files):
        print(f"Warning: Different number of snapshots ({len(with_files)} vs {len(no_files)})")
    
    n_frames = min(len(with_files), len(no_files))
    print(f"Processing {n_frames} frames...")
    
    # Create output directories
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    
    # First pass: collect all data for time histories
    print("Collecting time history data...")
    time_history = []
    momentum_with = []
    momentum_no = []
    energy_with = []
    energy_no = []
    
    initial_rho_c = None
    initial_R = None
    rho_max_global = 0
    
    for i in range(n_frames):
        df_w = load_snapshot(with_files[i])
        df_n = load_snapshot(no_files[i])
        
        # Extract time from filename (snapshot_XXXX.csv -> XXXX * dt)
        time = i * 2.0  # Assuming outputTime = 2.0
        time_history.append(time)
        
        # Momentum
        p_w, _, _, _ = compute_total_momentum(df_w)
        p_n, _, _, _ = compute_total_momentum(df_n)
        momentum_with.append(p_w)
        momentum_no.append(p_n)
        
        # Energy
        energy_with.append(compute_energies(df_w))
        energy_no.append(compute_energies(df_n))
        
        # Initial conditions
        if i == 0:
            initial_rho_c = compute_central_density(df_w)
            r = np.sqrt(df_w['pos_x']**2 + df_w['pos_y']**2 + 
                       (df_w['pos_z']**2 if 'pos_z' in df_w.columns else 0))
            initial_R = r.max()  # Full cloud radius for Lane-Emden scaling
        
        # Track max density
        rho_max_global = max(rho_max_global, df_n['dens'].max())
    
    rho_max_global *= 1.1  # Add margin
    
    # Second pass: render frames
    print("Rendering frames...")
    frame_paths = []
    for i in range(n_frames):
        output_path = os.path.join(frames_dir, f"frame_{i:04d}.png")
        render_frame(
            with_files[i], no_files[i], output_path, i, time_history[i],
            initial_rho_c, initial_R, rho_max_global,
            time_history[:i+1], momentum_with[:i+1], momentum_no[:i+1],
            energy_with[:i+1], energy_no[:i+1]
        )
        frame_paths.append(output_path)
        print(f"  Frame {i+1}/{n_frames}")
    
    # Create GIF
    if create_gif:
        print("Creating GIF animation...")
        try:
            from PIL import Image
            images = [Image.open(p) for p in frame_paths]
            gif_path = os.path.join(output_dir, "gradh_comparison.gif")
            images[0].save(gif_path, save_all=True, append_images=images[1:], 
                          duration=150, loop=0, optimize=True)
            print(f"✓ GIF saved: {gif_path}")
        except ImportError:
            print("PIL not available, skipping GIF creation")
    
    # Create HTML viewer
    if create_html:
        print("Creating interactive HTML viewer...")
        create_interactive_viewer(frames_dir, os.path.join(output_dir, "interactive_viewer.html"))
    
    print("Done!")


def create_interactive_viewer(frames_dir, output_html):
    """Create interactive HTML viewer with embedded frames."""
    frame_files = sorted(glob.glob(os.path.join(frames_dir, "frame_*.png")))
    
    encoded_frames = []
    for f in frame_files:
        with open(f, 'rb') as fp:
            encoded_frames.append(f'data:image/png;base64,{base64.b64encode(fp.read()).decode()}')
    
    html = f'''<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Grad-H Comparison</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ background: #1a1a1a; color: #fff; font-family: Arial, sans-serif; 
       display: flex; flex-direction: column; align-items: center; padding: 20px; }}
h1 {{ color: #FFD700; margin-bottom: 20px; }}
#viewer {{ max-width: 100%; border: 2px solid #444; }}
.controls {{ margin: 20px 0; display: flex; gap: 15px; align-items: center; }}
button {{ background: #333; color: #fff; border: 1px solid #666; padding: 10px 20px;
         cursor: pointer; font-size: 14px; border-radius: 4px; }}
button:hover {{ background: #444; }}
input[type="range"] {{ width: 400px; }}
#info {{ color: #aaa; font-size: 14px; }}
</style></head>
<body>
<h1>Grad-H Effect on Self-Gravitating Cloud</h1>
<img id="viewer" src="" alt="Frame">
<div class="controls">
  <button onclick="prev()">◀ Prev</button>
  <button onclick="togglePlay()">▶ Play</button>
  <button onclick="next()">Next ▶</button>
  <input type="range" id="slider" min="0" max="{len(encoded_frames)-1}" value="0" oninput="goTo(this.value)">
  <span id="info">Frame 1/{len(encoded_frames)}</span>
</div>
<script>
const frames = {encoded_frames};
let idx = 0, playing = false, timer = null;
const img = document.getElementById('viewer');
const slider = document.getElementById('slider');
const info = document.getElementById('info');

function show() {{
  img.src = frames[idx];
  slider.value = idx;
  info.textContent = `Frame ${{idx+1}}/${{frames.length}} (t=${{(idx*2).toFixed(0)}})`;
}}
function next() {{ idx = (idx + 1) % frames.length; show(); }}
function prev() {{ idx = (idx - 1 + frames.length) % frames.length; show(); }}
function goTo(i) {{ idx = parseInt(i); show(); }}
function togglePlay() {{
  playing = !playing;
  if (playing) {{ timer = setInterval(next, 150); }}
  else {{ clearInterval(timer); }}
}}
show();
</script>
</body></html>'''
    
    with open(output_html, 'w') as f:
        f.write(html)
    print(f"✓ HTML viewer: {output_html}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare grad-h vs no grad-h simulations")
    parser.add_argument("--with-gradh", required=True, help="Directory with grad-h results")
    parser.add_argument("--no-gradh", required=True, help="Directory without grad-h results")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--no-gif", action="store_true", help="Skip GIF creation")
    parser.add_argument("--no-html", action="store_true", help="Skip HTML viewer creation")
    
    args = parser.parse_args()
    
    process_all_frames(
        args.with_gradh, 
        args.no_gradh, 
        args.output,
        create_gif=not args.no_gif,
        create_html=not args.no_html
    )
