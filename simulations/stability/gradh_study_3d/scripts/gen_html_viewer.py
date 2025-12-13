#!/usr/bin/env python3
"""
Interactive HTML Viewer Generator for Grad-H Study
===================================================

Generates an interactive Plotly.js-based visualization comparing
self-gravitating cloud simulations with and without grad-h correction.

Features:
- Animated particle scatter plots
- Radial density profiles with Lane-Emden analytic solution overlay
- Grad-h (Omega) distribution vs radius
- Energy conservation comparison

Physics:
- Lane-Emden n=1.5 polytrope (gamma=5/3 ideal gas)
- Hernquist-Katz gravitational softening: epsilon = h/2
- Grad-h term: Omega = [1 + (h/D*rho) * sum_j m_j dW/dh]^(-1)

Uses SSOT module from scripts.shared.lane_emden for Lane-Emden solutions.

Usage:
    python gen_html_viewer.py [--case1 DIR] [--case2 DIR] [--output FILE]

Default paths:
    case1: results/case1_gradh_hk (with grad-h)
    case2: results/case2_no_gradh_hk (no grad-h)
    output: figures_v2/interactive_viewer.html
"""

import sys
from pathlib import Path

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import json
import os
import glob
import argparse

from scripts.shared.lane_emden import solve_lane_emden_spherical


def load_snapshot(filepath):
    """Load SPH snapshot CSV file."""
    return pd.read_csv(filepath, comment="#")


def load_energy_file(filepath):
    """Load energy.dat file with columns: time kinetic thermal potential total."""
    return pd.read_csv(filepath, comment="#", sep=r'\s+',
                       names=['time', 'kinetic', 'thermal', 'potential', 'total'])


def lane_emden_solution(n, xi_max=8.0, num_points=2000):
    """
    Get the Lane-Emden solution using SSOT solver.
    
    Returns:
        xi: dimensionless radial coordinate
        theta: dimensionless density parameter (rho/rho_c = theta^n)
        xi_1: first zero of theta (surface radius)
    """
    xi, theta = solve_lane_emden_spherical(n, xi_max=xi_max, n_points=num_points)
    theta = np.maximum(theta, 0)
    
    # Find surface (first zero of theta)
    idx = np.where(theta < 1e-8)[0]
    xi_1 = xi[idx[0]] if len(idx) > 0 else xi_max
    
    return xi, theta, xi_1


def lane_emden_density_profile(rho_c, R_cloud, n=1.5, n_points=100):
    """
    Compute physical density profile from Lane-Emden solution.
    
    For n=1.5 (gamma=5/3): xi_1 ≈ 3.6538
    
    Args:
        rho_c: Central density
        R_cloud: Physical cloud radius
        n: Polytropic index (default 1.5 for gamma=5/3)
        n_points: Number of output points
        
    Returns:
        r: Physical radius array
        rho: Density array
    """
    xi, theta, xi_1 = lane_emden_solution(n)
    
    # Scale to physical coordinates
    # alpha = R_cloud / xi_1 (dimensional scaling)
    alpha = R_cloud / xi_1
    r_phys = xi * alpha
    rho = rho_c * theta**n  # rho/rho_c = theta^n
    
    # Keep only interior points (r <= R with small buffer)
    mask = (r_phys <= R_cloud * 1.05) & (rho > 1e-6 * rho_c)
    r_out = r_phys[mask]
    rho_out = rho[mask]
    
    # Resample to uniform n_points
    r_interp = np.linspace(0, r_out.max(), n_points)
    rho_interp = np.interp(r_interp, r_out, rho_out)
    
    return r_interp.tolist(), rho_interp.tolist()


def compute_radial_profile(df):
    """Compute radial distance and return raw density/gradh data for scatter plot."""
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +
                (df['pos_z']**2 if 'pos_z' in df.columns else 0))
    rho = df['dens'].values
    gradh = df['gradh'].values if 'gradh' in df.columns else np.ones(len(df))
    
    return r.values, rho, gradh


def compute_acceleration_profile(df, max_points=500):
    """
    Compute radial acceleration profiles for scatter plot.
    
    Returns raw data: radius, |a_grav|, |a_pressure| for each particle.
    """
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +
                (df['pos_z']**2 if 'pos_z' in df.columns else 0))
    
    # Gravity acceleration magnitude
    if all(c in df.columns for c in ['grav_acc_x', 'grav_acc_y', 'grav_acc_z']):
        a_grav = np.sqrt(df['grav_acc_x']**2 + df['grav_acc_y']**2 + df['grav_acc_z']**2)
        # Fluid (pressure) acceleration = total - gravity
        fluid_acc_x = df['acc_x'] - df['grav_acc_x']
        fluid_acc_y = df['acc_y'] - df['grav_acc_y']
        fluid_acc_z = df['acc_z'] - df['grav_acc_z'] if 'acc_z' in df.columns else 0
        a_pres = np.sqrt(fluid_acc_x**2 + fluid_acc_y**2 + fluid_acc_z**2)
    else:
        a_grav = np.zeros(len(df))
        a_pres = np.sqrt(df['acc_x']**2 + df['acc_y']**2 + 
                         (df['acc_z']**2 if 'acc_z' in df.columns else 0))
    
    # Subsample if too many points
    if len(r) > max_points:
        idx = np.random.choice(len(r), max_points, replace=False)
        return r.iloc[idx].values, a_grav.iloc[idx].values, a_pres.iloc[idx].values
    
    return r.values, a_grav.values, a_pres.values


def compute_statistics(df):
    """
    Compute summary statistics for a snapshot.
    
    Includes:
    - Energetics: kinetic (K), thermal (U), potential (W), virial ratio (2K/|W|)
    - Momentum: linear (P), angular (L)
    - Force decomposition: gravity vs fluid (pressure) forces
    """
    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +
                (df['pos_z']**2 if 'pos_z' in df.columns else 0))
    r_max = r.max()
    
    # Central density (inner 10% of radius)
    central_mask = r < 0.1 * r_max
    rho_central = df.loc[central_mask, 'dens'].mean() if central_mask.sum() > 0 else df['dens'].max()
    
    # Grad-h statistics
    gradh_mean = df['gradh'].mean() if 'gradh' in df.columns else 1.0
    gradh_min = df['gradh'].min() if 'gradh' in df.columns else 1.0
    gradh_max = df['gradh'].max() if 'gradh' in df.columns else 1.0
    
    mass = df['mass']
    
    # === ENERGIES (computed from particle data) ===
    # Kinetic energy: K = (1/2) sum(m * v^2)
    v2 = df['vel_x']**2 + df['vel_y']**2 + (df['vel_z']**2 if 'vel_z' in df.columns else 0)
    E_kinetic = 0.5 * (mass * v2).sum()
    
    # Thermal energy: U = sum(m * u) where u is specific internal energy
    E_thermal = (mass * df['ene']).sum() if 'ene' in df.columns else 0.0
    
    # Gravitational potential: W = (1/2) sum(m * phi)
    # Factor of 1/2 to avoid double counting pairs
    E_potential = 0.5 * (mass * df['phi']).sum() if 'phi' in df.columns else 0.0
    
    # Total energy
    E_total = E_kinetic + E_thermal + E_potential
    
    # Virial ratio: should be 2K/|W| = 1 for equilibrium
    # For polytrope: 2K + W = 0, and 3(gamma-1)U + W = 0
    virial_ratio = 2.0 * E_kinetic / abs(E_potential) if E_potential != 0 else 0.0
    
    # === MOMENTUM ===
    # Linear momentum (P = sum(m * v)) - should be ~0 for symmetric system
    px = (mass * df['vel_x']).sum() if 'vel_x' in df.columns else 0.0
    py = (mass * df['vel_y']).sum() if 'vel_y' in df.columns else 0.0
    pz = (mass * df['vel_z']).sum() if 'vel_z' in df.columns else 0.0
    p_total = np.sqrt(px**2 + py**2 + pz**2)
    
    # Angular momentum (L = sum(r x p))
    if all(c in df.columns for c in ['pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']):
        Lx = (mass * (df['pos_y'] * df['vel_z'] - df['pos_z'] * df['vel_y'])).sum()
        Ly = (mass * (df['pos_z'] * df['vel_x'] - df['pos_x'] * df['vel_z'])).sum()
        Lz = (mass * (df['pos_x'] * df['vel_y'] - df['pos_y'] * df['vel_x'])).sum()
    else:
        Lx = Ly = Lz = 0.0
    L_total = np.sqrt(Lx**2 + Ly**2 + Lz**2)
    
    # === FORCE ANALYSIS ===
    # Total acceleration magnitude
    a_total = np.sqrt(df['acc_x']**2 + df['acc_y']**2 + 
                      (df['acc_z']**2 if 'acc_z' in df.columns else 0))
    
    # Gravity acceleration magnitude
    if all(c in df.columns for c in ['grav_acc_x', 'grav_acc_y', 'grav_acc_z']):
        a_grav = np.sqrt(df['grav_acc_x']**2 + df['grav_acc_y']**2 + df['grav_acc_z']**2)
        # Fluid (pressure) acceleration = total - gravity
        fluid_acc_x = df['acc_x'] - df['grav_acc_x']
        fluid_acc_y = df['acc_y'] - df['grav_acc_y']
        fluid_acc_z = df['acc_z'] - df['grav_acc_z'] if 'acc_z' in df.columns else 0
        a_fluid = np.sqrt(fluid_acc_x**2 + fluid_acc_y**2 + fluid_acc_z**2)
    else:
        a_grav = np.zeros(len(df))
        a_fluid = a_total
    
    # RMS force magnitudes
    F_grav_rms = np.sqrt((mass * a_grav**2).sum())
    F_fluid_rms = np.sqrt((mass * a_fluid**2).sum())
    F_total_rms = np.sqrt((mass * a_total**2).sum())
    
    return {
        'rho_central': float(rho_central),
        'rho_max': float(df['dens'].max()),
        'gradh_mean': float(gradh_mean),
        'gradh_min': float(gradh_min),
        'gradh_max': float(gradh_max),
        'r_max': float(r_max),
        # Energies from particle data
        'E_kinetic': float(E_kinetic),
        'E_thermal': float(E_thermal),
        'E_potential': float(E_potential),
        'E_total': float(E_total),
        'virial_ratio': float(virial_ratio),
        # Momentum
        'px': float(px),
        'py': float(py),
        'pz': float(pz),
        'p_total': float(p_total),
        'Lx': float(Lx),
        'Ly': float(Ly),
        'Lz': float(Lz),
        'L_total': float(L_total),
        # Forces
        'F_grav_rms': float(F_grav_rms),
        'F_fluid_rms': float(F_fluid_rms),
        'F_total_rms': float(F_total_rms),
    }


def process_case(results_dir, sample_particles=500):
    """
    Process all snapshots in a results directory.
    
    Returns dict with frames, energy data, and initial conditions.
    """
    snapshots = sorted(glob.glob(os.path.join(results_dir, "snapshot_*.csv")))
    energy_file = os.path.join(results_dir, "energy.dat")
    
    if not snapshots:
        raise FileNotFoundError(f"No snapshots found in {results_dir}")
    
    energy_data = load_energy_file(energy_file) if os.path.exists(energy_file) else None
    
    frames = []
    initial_rho_c = None
    initial_R = None
    
    for i, snap_file in enumerate(snapshots):
        df = load_snapshot(snap_file)
        stats = compute_statistics(df)
        r_raw, rho_raw, gradh_raw = compute_radial_profile(df)
        r_acc, a_grav, a_pres = compute_acceleration_profile(df)
        
        if i == 0:
            initial_rho_c = stats['rho_central']
            initial_R = stats['r_max']
        
        # Sample particles for visualization
        n_sample = min(sample_particles, len(df))
        indices = np.random.choice(len(df), n_sample, replace=False)
        sampled = df.iloc[indices]
        
        frame = {
            'index': i,
            'time': i * 2.0,  # Assuming dt_snapshot = 2.0
            'stats': stats,
            'radial': {
                'r': r_raw.tolist(),
                'rho': rho_raw.tolist(),
                'gradh': gradh_raw.tolist()
            },
            'accel': {
                'r': r_acc.tolist(),
                'a_grav': a_grav.tolist(),
                'a_pres': a_pres.tolist()
            },
            'particles': {
                'x': sampled['pos_x'].tolist(),
                'y': sampled['pos_y'].tolist(),
                'rho': sampled['dens'].tolist(),
                'gradh': sampled['gradh'].tolist() if 'gradh' in sampled.columns else [1.0] * len(sampled),
            }
        }
        frames.append(frame)
    
    # Energy time series
    if energy_data is not None:
        energy_series = {
            'time': energy_data['time'].tolist(),
            'kinetic': energy_data['kinetic'].tolist(),
            'thermal': energy_data['thermal'].tolist(),
            'potential': energy_data['potential'].tolist(),
            'total': energy_data['total'].tolist(),
        }
    else:
        energy_series = {
            'time': [f['time'] for f in frames],
            'total': [0.0] * len(frames)
        }
    
    return {
        'frames': frames,
        'energy': energy_series,
        'initial_rho_c': initial_rho_c,
        'initial_R': initial_R
    }


def generate_html(data_with, data_no, lane_emden_data, output_path):
    """Generate the interactive HTML viewer with improved plots."""
    
    # Extract momentum time series from frames
    momentum_with = {
        'time': [f['time'] for f in data_with['frames']],
        'px': [f['stats']['px'] for f in data_with['frames']],
        'py': [f['stats']['py'] for f in data_with['frames']],
        'pz': [f['stats']['pz'] for f in data_with['frames']],
        'p_total': [f['stats']['p_total'] for f in data_with['frames']],
        'Lx': [f['stats']['Lx'] for f in data_with['frames']],
        'Ly': [f['stats']['Ly'] for f in data_with['frames']],
        'Lz': [f['stats']['Lz'] for f in data_with['frames']],
        'L_total': [f['stats']['L_total'] for f in data_with['frames']],
    }
    momentum_no = {
        'time': [f['time'] for f in data_no['frames']],
        'px': [f['stats']['px'] for f in data_no['frames']],
        'py': [f['stats']['py'] for f in data_no['frames']],
        'pz': [f['stats']['pz'] for f in data_no['frames']],
        'p_total': [f['stats']['p_total'] for f in data_no['frames']],
        'Lx': [f['stats']['Lx'] for f in data_no['frames']],
        'Ly': [f['stats']['Ly'] for f in data_no['frames']],
        'Lz': [f['stats']['Lz'] for f in data_no['frames']],
        'L_total': [f['stats']['L_total'] for f in data_no['frames']],
    }
    
    # Extract particle-computed energy and virial ratio
    virial_with = {
        'time': [f['time'] for f in data_with['frames']],
        'E_kinetic': [f['stats']['E_kinetic'] for f in data_with['frames']],
        'E_thermal': [f['stats']['E_thermal'] for f in data_with['frames']],
        'E_potential': [f['stats']['E_potential'] for f in data_with['frames']],
        'E_total': [f['stats']['E_total'] for f in data_with['frames']],
        'virial_ratio': [f['stats']['virial_ratio'] for f in data_with['frames']],
        'F_grav_rms': [f['stats']['F_grav_rms'] for f in data_with['frames']],
        'F_fluid_rms': [f['stats']['F_fluid_rms'] for f in data_with['frames']],
    }
    virial_no = {
        'time': [f['time'] for f in data_no['frames']],
        'E_kinetic': [f['stats']['E_kinetic'] for f in data_no['frames']],
        'E_thermal': [f['stats']['E_thermal'] for f in data_no['frames']],
        'E_potential': [f['stats']['E_potential'] for f in data_no['frames']],
        'E_total': [f['stats']['E_total'] for f in data_no['frames']],
        'virial_ratio': [f['stats']['virial_ratio'] for f in data_no['frames']],
        'F_grav_rms': [f['stats']['F_grav_rms'] for f in data_no['frames']],
        'F_fluid_rms': [f['stats']['F_fluid_rms'] for f in data_no['frames']],
    }
    
    html_template = '''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Grad-H Effect Analysis with Lane-Emden Solution</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',sans-serif;background:linear-gradient(135deg,#1a1a2e,#16213e,#0f3460);color:#e8e8e8;min-height:100vh;padding:20px}
h1{text-align:center;color:#ffd700;margin-bottom:10px;font-size:2em}
.subtitle{text-align:center;color:#aaa;margin-bottom:20px}
.container{display:grid;grid-template-columns:1fr 1fr;gap:20px;max-width:1800px;margin:0 auto}
.panel{background:rgba(30,30,50,0.9);border-radius:12px;padding:15px;box-shadow:0 4px 20px rgba(0,0,0,0.4);border:1px solid rgba(100,100,150,0.3)}
.panel h3{margin-bottom:10px;font-size:1.1em;border-bottom:1px solid rgba(100,100,150,0.3);padding-bottom:5px}
.with-gradh h3{color:#00ff88}
.no-gradh h3{color:#ff4444}
.theory{grid-column:span 2}
.controls{display:flex;justify-content:center;align-items:center;gap:20px;margin:20px 0;padding:15px;background:rgba(30,30,50,0.9);border-radius:12px;flex-wrap:wrap}
button{background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;padding:10px 20px;border-radius:8px;cursor:pointer;font-size:1em}
button:hover{transform:translateY(-2px);box-shadow:0 4px 15px rgba(102,126,234,0.4)}
#slider{width:400px;height:8px;-webkit-appearance:none;background:linear-gradient(90deg,#667eea,#764ba2);border-radius:4px}
#slider::-webkit-slider-thumb{-webkit-appearance:none;width:20px;height:20px;background:#ffd700;border-radius:50%;cursor:pointer}
.info-box{display:flex;gap:30px;color:#ccc;font-size:1.1em}
.info-box span{color:#ffd700;font-weight:bold}
.stats-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin-top:10px}
.stat-item{background:rgba(50,50,80,0.5);padding:8px;border-radius:6px;text-align:center}
.stat-label{color:#888;font-size:0.85em}
.stat-value{color:#fff;font-size:1.1em;font-weight:bold}
.stat-value.good{color:#00ff88}
.stat-value.bad{color:#ff4444}
.theory-content{display:grid;grid-template-columns:1fr 1fr;gap:20px}
.derivation{padding:15px;background:rgba(40,40,70,0.5);border-radius:8px;line-height:1.8}
.derivation h4{color:#ffd700;margin-bottom:10px}
.equation{background:rgba(60,60,100,0.5);padding:10px;border-radius:6px;margin:8px 0;text-align:center}
.key-insight{background:linear-gradient(135deg,rgba(0,255,136,0.1),rgba(0,200,100,0.1));border-left:3px solid #00ff88;padding:10px 15px;margin:10px 0;border-radius:0 6px 6px 0}
.warning{background:linear-gradient(135deg,rgba(255,68,68,0.1),rgba(200,50,50,0.1));border-left:3px solid #ff4444;padding:10px 15px;margin:10px 0;border-radius:0 6px 6px 0}
.analytic{background:linear-gradient(135deg,rgba(0,255,255,0.1),rgba(0,200,200,0.1));border-left:3px solid #00ffff;padding:10px 15px;margin:10px 0;border-radius:0 6px 6px 0}
.plot-container{height:350px}
.full-width{grid-column:span 2}
</style>
</head>
<body>
<h1>Grad-H Correction Effect on Self-Gravitating Cloud</h1>
<p class="subtitle">First Principles Analysis with Lane-Emden Analytic Solution (n=1.5 polytrope)</p>
<div class="controls">
<button onclick="prevFrame()">&#9664; Prev</button>
<button onclick="togglePlay()" id="playBtn">&#9654; Play</button>
<button onclick="nextFrame()">Next &#9654;</button>
<input type="range" id="slider" min="0" max="10" value="0" oninput="goToFrame(this.value)">
<div class="info-box"><div>Frame: <span id="frameNum">1</span>/<span id="totalFrames">11</span></div><div>Time: <span id="timeVal">0.0</span></div></div>
</div>
<div class="container">
<div class="panel with-gradh"><h3>WITH Ω ≠ 1 (Stable): ρ(x,y) scatter</h3><div id="scatter-with" class="plot-container"></div>
<div class="stats-grid"><div class="stat-item"><div class="stat-label">ρ_central</div><div class="stat-value good" id="rho-central-with">-</div></div>
<div class="stat-item"><div class="stat-label">Ω_mean</div><div class="stat-value good" id="gradh-mean-with">-</div></div>
<div class="stat-item"><div class="stat-label">Ω_range</div><div class="stat-value" id="gradh-range-with">-</div></div></div></div>
<div class="panel no-gradh"><h3>NO Ω (=1, Collapsing): ρ(x,y) scatter</h3><div id="scatter-no" class="plot-container"></div>
<div class="stats-grid"><div class="stat-item"><div class="stat-label">ρ_central</div><div class="stat-value bad" id="rho-central-no">-</div></div>
<div class="stat-item"><div class="stat-label">Ω (forced)</div><div class="stat-value bad" id="gradh-mean-no">1.00</div></div>
<div class="stat-item"><div class="stat-label">ρ increase</div><div class="stat-value bad" id="rho-increase">-</div></div></div></div>
<div class="panel"><h3>Density Profile: ρ(r) vs Lane-Emden θ<sup>n</sup></h3><div id="radial-profile" class="plot-container"></div></div>
<div class="panel"><h3>Grad-H Factor: Ω(r) = [1 + (h/Dρ)Σm∂W/∂h]⁻¹</h3><div id="gradh-distribution" class="plot-container"></div></div>
<div class="panel"><h3>Energy Components (With Ω): K, U, W, E<sub>tot</sub></h3><div id="energy-components-with" class="plot-container"></div></div>
<div class="panel"><h3>Energy Components (No Ω): K, U, W, E<sub>tot</sub></h3><div id="energy-components-no" class="plot-container"></div></div>
<div class="panel"><h3>Energy Error: ΔE/|E₀| = (E(t)−E₀)/|E₀|</h3><div id="energy-total" class="plot-container"></div></div>
<div class="panel"><h3>Virial Ratio: 2K/|W| (equilibrium → small)</h3><div id="virial-plot" class="plot-container"></div></div>
<div class="panel"><h3>Linear Momentum: P = Σmv</h3><div id="momentum-plot" class="plot-container"></div></div>
<div class="panel"><h3>Angular Momentum: L = Σ(r×p)</h3><div id="angular-momentum-plot" class="plot-container"></div></div>
<div class="panel"><h3>Acceleration (With Ω): |a<sub>grav</sub>|(r) vs |a<sub>pres</sub>|(r)</h3><div id="accel-with-plot" class="plot-container"></div></div>
<div class="panel"><h3>Acceleration (No Ω): |a<sub>grav</sub>|(r) vs |a<sub>pres</sub>|(r)</h3><div id="accel-no-plot" class="plot-container"></div></div>
<div class="panel"><h3>Particle Energy: E = ½Σmv² + Σmu + ½Σmφ</h3><div id="particle-energy-plot" class="plot-container"></div></div>
<div class="panel theory"><h3>First Principles: Why Grad-H Matters</h3>
<div class="theory-content">
<div class="derivation"><h4>The Grad-H Correction Term</h4>
<p>In SPH, smoothing length adapts to density: \\( h_i = \\eta (m_i/\\rho_i)^{1/D} \\)</p>
<p>This creates an implicit dependence: \\( \\rho_i \\leftrightarrow h_i \\)</p>
<div class="equation">\\( \\rho_i = \\sum_j m_j W(|\\mathbf{r}_i - \\mathbf{r}_j|, h_i) \\)</div>
<p>Taking the derivative properly (chain rule):</p>
<div class="equation">\\( \\frac{\\partial \\rho_k}{\\partial \\mathbf{x}_i} = \\Omega_k \\sum_j m_j \\nabla_i W_{kj} \\)</div>
<div class="key-insight"><strong>The Omega factor:</strong><br>\\( \\displaystyle \\Omega_i = \\left[1 + \\frac{h_i}{D\\rho_i}\\sum_j m_j \\frac{\\partial W_{ij}}{\\partial h_i}\\right]^{-1} \\)</div>
<div class="analytic"><strong>Lane-Emden Equation (n=1.5):</strong><br>\\( \\frac{1}{\\xi^2}\\frac{d}{d\\xi}\\left(\\xi^2\\frac{d\\theta}{d\\xi}\\right) = -\\theta^{3/2} \\)<br>
Physical density: \\( \\rho(r) = \\rho_c \\theta(\\xi)^{3/2} \\)</div></div>
<div class="derivation"><h4>Effect on Gravitational Potential</h4>
<p>Hernquist-Katz softening uses: \\( \\epsilon = h/2 \\)</p>
<p>So gravitational potential \\( \\phi_i \\) depends on \\( h_i(\\rho_i) \\)!</p>
<div class="warning"><strong>Without Ω correction (Ω = 1):</strong><br>- Ignores \\( \\partial\\phi/\\partial h \\) contribution<br>- Gravitational energy NOT conserved<br>- Artificial numerical collapse occurs<br>- Deviates from Lane-Emden equilibrium</div>
<div class="key-insight"><strong>With Ω correction (Ω ≠ 1):</strong><br>- Properly accounts for h(ρ) in potential<br>- Total energy exactly conserved<br>- Cloud remains in hydrostatic equilibrium<br>- Matches Lane-Emden analytic solution</div>
<p style="margin-top:15px;color:#ffd700"><strong>Physical insight:</strong> The cyan curve shows the exact Lane-Emden n=1.5 polytrope solution. With proper grad-h correction, the SPH simulation maintains this equilibrium.</p></div>
</div></div>
</div>
<script>
var dataWith = ''' + json.dumps(data_with) + ''';
var dataNo = ''' + json.dumps(data_no) + ''';
var laneEmden = ''' + json.dumps(lane_emden_data) + ''';
var momentumWith = ''' + json.dumps(momentum_with) + ''';
var momentumNo = ''' + json.dumps(momentum_no) + ''';
var virialWith = ''' + json.dumps(virial_with) + ''';
var virialNo = ''' + json.dumps(virial_no) + ''';
var currentFrame=0,playing=false,playInterval=null,totalFrames=dataWith.frames.length;
document.getElementById('totalFrames').textContent=totalFrames;
document.getElementById('slider').max=totalFrames-1;
var colorscale=[[0,'#000033'],[0.2,'#003366'],[0.4,'#006699'],[0.6,'#ff9900'],[0.8,'#ff3300'],[1,'#ffff00']];
var rho0_with=dataWith.frames[0].stats.rho_central,rho0_no=dataNo.frames[0].stats.rho_central;
var E0_with=Math.abs(dataWith.energy.total[0]),E0_no=Math.abs(dataNo.energy.total[0]);
var plotLayout={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',font:{color:'#ccc'},margin:{t:10,b:40,l:60,r:20}};

function updatePlots(){
var fw=dataWith.frames[currentFrame],fn=dataNo.frames[currentFrame];
document.getElementById('frameNum').textContent=currentFrame+1;
document.getElementById('timeVal').textContent=fw.time.toFixed(1);
document.getElementById('rho-central-with').textContent=fw.stats.rho_central.toFixed(3);
document.getElementById('gradh-mean-with').textContent=fw.stats.gradh_mean.toFixed(3);
document.getElementById('gradh-range-with').textContent='['+fw.stats.gradh_min.toFixed(2)+','+fw.stats.gradh_max.toFixed(2)+']';
document.getElementById('rho-central-no').textContent=fn.stats.rho_central.toFixed(3);
var rhoInc=((fn.stats.rho_central/rho0_no-1)*100).toFixed(1);
document.getElementById('rho-increase').textContent='+'+rhoInc+'%';
var allRho=fw.particles.rho.concat(fn.particles.rho),rhoMin=Math.min.apply(null,allRho),rhoMax=Math.max.apply(null,allRho);

// Particle scatter plots
Plotly.react('scatter-with',[{x:fw.particles.x,y:fw.particles.y,mode:'markers',type:'scatter',marker:{size:4,color:fw.particles.rho,colorscale:colorscale,cmin:rhoMin,cmax:rhoMax,colorbar:{title:'ρ',thickness:15}}}],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'x',range:[-1.3,1.3],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'y',range:[-1.3,1.3],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)',scaleanchor:'x'},margin:{t:10,b:40,l:50,r:80},font:{color:'#ccc'}},{responsive:true});

Plotly.react('scatter-no',[{x:fn.particles.x,y:fn.particles.y,mode:'markers',type:'scatter',marker:{size:4,color:fn.particles.rho,colorscale:colorscale,cmin:rhoMin,cmax:rhoMax,colorbar:{title:'ρ',thickness:15}}}],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'x',range:[-1.3,1.3],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'y',range:[-1.3,1.3],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)',scaleanchor:'x'},margin:{t:10,b:40,l:50,r:80},font:{color:'#ccc'}},{responsive:true});

// Radial density profile (raw scatter) - FIXED Y RANGE
Plotly.react('radial-profile',[
{x:laneEmden.r,y:laneEmden.rho,mode:'lines',name:'Lane-Emden (n=1.5)',line:{color:'#00ffff',width:3}},
{x:fw.radial.r,y:fw.radial.rho,mode:'markers',name:'With Ω',marker:{size:4,color:'#00ff88',opacity:0.6}},
{x:fn.radial.r,y:fn.radial.rho,mode:'markers',name:'No Ω',marker:{size:4,color:'#ff4444',opacity:0.6,symbol:'x'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Radius r',range:[0,1.2],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Density ρ',type:'log',range:[-2,1],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.6,y:0.95,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'}},{responsive:true});

// Grad-h distribution (raw scatter) - FIXED Y RANGE
Plotly.react('gradh-distribution',[{x:fw.radial.r,y:fw.radial.gradh,mode:'markers',name:'Ω(r)',marker:{size:4,color:'#00ff88',opacity:0.6}},{x:fn.radial.r,y:fn.radial.r.map(function(){return 1}),mode:'lines',name:'Ω=1 (forced)',line:{color:'#ff4444',width:2,dash:'dash'}}],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Radius r',range:[0,1.2],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Ω',range:[0.8,1.5],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.6,y:0.95,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:0,x1:1.2,y0:1,y1:1,line:{color:'#888',width:1,dash:'dot'}}]},{responsive:true});

var ct=fw.time;

// Energy components - WITH grad-h - FIXED Y RANGE
Plotly.react('energy-components-with',[
{x:dataWith.energy.time,y:dataWith.energy.kinetic,mode:'lines',name:'Kinetic',line:{color:'#ff6b6b',width:2}},
{x:dataWith.energy.time,y:dataWith.energy.thermal,mode:'lines',name:'Thermal',line:{color:'#ffd93d',width:2}},
{x:dataWith.energy.time,y:dataWith.energy.potential,mode:'lines',name:'Potential (grav)',line:{color:'#6bcb77',width:2}},
{x:dataWith.energy.time,y:dataWith.energy.total,mode:'lines',name:'Total',line:{color:'#4d96ff',width:3}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Energy',range:[-1.0,0.5],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.02,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:-1,y1:0.5,line:{color:'#ffd700',width:2,dash:'dot'}}]},{responsive:true});

// Energy components - NO grad-h - FIXED Y RANGE
Plotly.react('energy-components-no',[
{x:dataNo.energy.time,y:dataNo.energy.kinetic,mode:'lines',name:'Kinetic',line:{color:'#ff6b6b',width:2}},
{x:dataNo.energy.time,y:dataNo.energy.thermal,mode:'lines',name:'Thermal',line:{color:'#ffd93d',width:2}},
{x:dataNo.energy.time,y:dataNo.energy.potential,mode:'lines',name:'Potential (grav)',line:{color:'#6bcb77',width:2}},
{x:dataNo.energy.time,y:dataNo.energy.total,mode:'lines',name:'Total',line:{color:'#4d96ff',width:3}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Energy',range:[-1.0,0.5],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.02,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:-1,y1:0.5,line:{color:'#ffd700',width:2,dash:'dot'}}]},{responsive:true});

// Total energy comparison (normalized) - FIXED Y RANGE
Plotly.react('energy-total',[
{x:dataWith.energy.time,y:dataWith.energy.total.map(function(e,i){return (e-dataWith.energy.total[0])/E0_with*100}),mode:'lines',name:'With Ω (ΔE/|E₀|%)',line:{color:'#00ff88',width:3}},
{x:dataNo.energy.time,y:dataNo.energy.total.map(function(e,i){return (e-dataNo.energy.total[0])/E0_no*100}),mode:'lines',name:'No Ω (ΔE/|E₀|%)',line:{color:'#ff4444',width:3,dash:'dash'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Energy Error ΔE/|E₀| (%)',range:[-1,1],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.02,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:-1,y1:1,line:{color:'#ffd700',width:2,dash:'dot'}},{type:'line',x0:0,x1:dataWith.energy.time[dataWith.energy.time.length-1],y0:0,y1:0,line:{color:'#888',width:1,dash:'dot'}}]},{responsive:true});

// Momentum conservation (linear) - FIXED Y RANGE
Plotly.react('momentum-plot',[
{x:momentumWith.time,y:momentumWith.px,mode:'lines',name:'Px (with Ω)',line:{color:'#00ff88',width:2}},
{x:momentumWith.time,y:momentumWith.py,mode:'lines',name:'Py (with Ω)',line:{color:'#00cc66',width:2}},
{x:momentumWith.time,y:momentumWith.pz,mode:'lines',name:'Pz (with Ω)',line:{color:'#009944',width:2}},
{x:momentumNo.time,y:momentumNo.px,mode:'lines',name:'Px (no Ω)',line:{color:'#ff4444',width:2,dash:'dash'}},
{x:momentumNo.time,y:momentumNo.py,mode:'lines',name:'Py (no Ω)',line:{color:'#cc3333',width:2,dash:'dash'}},
{x:momentumNo.time,y:momentumNo.pz,mode:'lines',name:'Pz (no Ω)',line:{color:'#992222',width:2,dash:'dash'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Linear Momentum P = Σm·v',range:[-0.01,0.01],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.02,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:-0.01,y1:0.01,line:{color:'#ffd700',width:2,dash:'dot'}},{type:'line',x0:0,x1:momentumWith.time[momentumWith.time.length-1],y0:0,y1:0,line:{color:'#888',width:1,dash:'dot'}}]},{responsive:true});

// Angular momentum (all components) - FIXED Y RANGE
Plotly.react('angular-momentum-plot',[
{x:momentumWith.time,y:momentumWith.Lx,mode:'lines',name:'Lx (with Ω)',line:{color:'#00ff88',width:2}},
{x:momentumWith.time,y:momentumWith.Ly,mode:'lines',name:'Ly (with Ω)',line:{color:'#00cc66',width:2}},
{x:momentumWith.time,y:momentumWith.Lz,mode:'lines',name:'Lz (with Ω)',line:{color:'#009944',width:2}},
{x:momentumNo.time,y:momentumNo.Lx,mode:'lines',name:'Lx (no Ω)',line:{color:'#ff4444',width:2,dash:'dash'}},
{x:momentumNo.time,y:momentumNo.Ly,mode:'lines',name:'Ly (no Ω)',line:{color:'#cc3333',width:2,dash:'dash'}},
{x:momentumNo.time,y:momentumNo.Lz,mode:'lines',name:'Lz (no Ω)',line:{color:'#992222',width:2,dash:'dash'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Angular Momentum L = Σ(r×mv)',range:[-0.005,0.005],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.02,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:-0.005,y1:0.005,line:{color:'#ffd700',width:2,dash:'dot'}},{type:'line',x0:0,x1:momentumWith.time[momentumWith.time.length-1],y0:0,y1:0,line:{color:'#888',width:1,dash:'dot'}}]},{responsive:true});

// Virial ratio - FIXED Y RANGE
Plotly.react('virial-plot',[
{x:virialWith.time,y:virialWith.virial_ratio,mode:'lines',name:'2K/|W| (with Ω)',line:{color:'#00ff88',width:3}},
{x:virialNo.time,y:virialNo.virial_ratio,mode:'lines',name:'2K/|W| (no Ω)',line:{color:'#ff4444',width:3,dash:'dash'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Virial Ratio 2K/|W|',range:[0,0.05],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.7,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:0,y1:0.05,line:{color:'#ffd700',width:2,dash:'dot'}}]},{responsive:true});

// Acceleration scatter - WITH grad-h (raw particle data) - FIXED Y RANGE
Plotly.react('accel-with-plot',[
{x:fw.accel.r,y:fw.accel.a_grav,mode:'markers',name:'|a_grav|',marker:{size:4,color:'#6bcb77',opacity:0.6}},
{x:fw.accel.r,y:fw.accel.a_pres,mode:'markers',name:'|a_pres|',marker:{size:4,color:'#ff6b6b',opacity:0.6,symbol:'x'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Radius r',range:[0,1.2],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Acceleration |a|',type:'log',range:[-2,1],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.7,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'}},{responsive:true});

// Acceleration scatter - NO grad-h (raw particle data) - FIXED Y RANGE
Plotly.react('accel-no-plot',[
{x:fn.accel.r,y:fn.accel.a_grav,mode:'markers',name:'|a_grav|',marker:{size:4,color:'#6bcb77',opacity:0.6}},
{x:fn.accel.r,y:fn.accel.a_pres,mode:'markers',name:'|a_pres|',marker:{size:4,color:'#ff6b6b',opacity:0.6,symbol:'x'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Radius r',range:[0,1.2],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Acceleration |a|',type:'log',range:[-2,1],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.7,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'}},{responsive:true});

// Particle-computed energies comparison - FIXED Y RANGE
var Ep0_with=Math.abs(virialWith.E_total[0]),Ep0_no=Math.abs(virialNo.E_total[0]);
Plotly.react('particle-energy-plot',[
{x:virialWith.time,y:virialWith.E_total.map(function(e){return (e-virialWith.E_total[0])/Ep0_with*100}),mode:'lines',name:'ΔE/|E₀|% (with Ω)',line:{color:'#00ff88',width:3}},
{x:virialNo.time,y:virialNo.E_total.map(function(e){return (e-virialNo.E_total[0])/Ep0_no*100}),mode:'lines',name:'ΔE/|E₀|% (no Ω)',line:{color:'#ff4444',width:3,dash:'dash'}}
],{paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(30,30,50,0.5)',xaxis:{title:'Time',color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},yaxis:{title:'Energy Error (from Σm(½v²+u)+½Σmφ)',range:[-1,1],color:'#ccc',gridcolor:'rgba(100,100,150,0.3)'},legend:{x:0.02,y:0.98,bgcolor:'rgba(30,30,50,0.8)'},margin:{t:10,b:40,l:60,r:20},font:{color:'#ccc'},shapes:[{type:'line',x0:ct,x1:ct,y0:-1,y1:1,line:{color:'#ffd700',width:2,dash:'dot'}},{type:'line',x0:0,x1:virialWith.time[virialWith.time.length-1],y0:0,y1:0,line:{color:'#888',width:1,dash:'dot'}}]},{responsive:true});
}

function goToFrame(i){currentFrame=parseInt(i);document.getElementById('slider').value=currentFrame;updatePlots()}
function nextFrame(){currentFrame=(currentFrame+1)%totalFrames;document.getElementById('slider').value=currentFrame;updatePlots()}
function prevFrame(){currentFrame=(currentFrame-1+totalFrames)%totalFrames;document.getElementById('slider').value=currentFrame;updatePlots()}
function togglePlay(){playing=!playing;var b=document.getElementById('playBtn');if(playing){b.innerHTML='&#9208; Pause';playInterval=setInterval(nextFrame,300)}else{b.innerHTML='&#9654; Play';clearInterval(playInterval)}}
updatePlots();
</script>
</body></html>'''
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(html_template)
    print(f"Generated: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate interactive HTML viewer for Grad-H study',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--case1', default='results/case1_gradh_hk',
                        help='Directory with grad-h ON results (default: results/case1_gradh_hk)')
    parser.add_argument('--case2', default='results/case2_no_gradh_hk',
                        help='Directory with grad-h OFF results (default: results/case2_no_gradh_hk)')
    parser.add_argument('--output', '-o', default='figures_v2/interactive_viewer.html',
                        help='Output HTML file path (default: figures_v2/interactive_viewer.html)')
    parser.add_argument('--particles', type=int, default=600,
                        help='Number of particles to sample for visualization (default: 600)')
    
    args = parser.parse_args()
    
    # Get script directory for relative paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gradh_study_dir = os.path.dirname(script_dir)
    
    # Resolve paths relative to gradh_study directory
    case1_dir = os.path.join(gradh_study_dir, args.case1)
    case2_dir = os.path.join(gradh_study_dir, args.case2)
    output_path = os.path.join(gradh_study_dir, args.output)
    
    print(f"Processing {args.case1}...")
    data_with = process_case(case1_dir, args.particles)
    
    print(f"Processing {args.case2}...")
    data_no = process_case(case2_dir, args.particles)
    
    print(f"Loaded {len(data_with['frames'])} frames")
    
    # Compute Lane-Emden analytic solution
    print("Computing Lane-Emden analytic solution...")
    le_r, le_rho = lane_emden_density_profile(
        data_with['initial_rho_c'],
        data_with['initial_R'],
        n=1.5,
        n_points=80
    )
    lane_emden_data = {'r': le_r, 'rho': le_rho}
    print(f"Lane-Emden: rho_c={data_with['initial_rho_c']:.4f}, R={data_with['initial_R']:.4f}")
    
    # Generate HTML
    generate_html(data_with, data_no, lane_emden_data, output_path)
    
    return 0


if __name__ == '__main__':
    exit(main())
