#!/usr/bin/env python3#!/usr/bin/env python3#!/usr/bin/env python3

"""

Generate Interactive HTML Viewer for Grad-H Study.""""""



Creates a standalone HTML file that visualizes the grad-h term effectGenerate Interactive HTML Viewer for Grad-H StudyInteractive HTML Animation Viewer for Self-Gravitating Cloud Simulations

from first principles using raw simulation data.

===========================================================================================================================

The grad-h correction term accounts for variable smoothing length h(rho):

    Omega_i = [1 + (h_i / D*rho_i) sum_j m_j dW_ij/dh_i]^{-1}Creates a standalone HTML file that visualizes the grad-h (Ω) term effectCreates a self-contained HTML file with interactive controls for viewing



When Omega != 1 (with correction): Energy is conserved -> Cloud stablefrom first principles using raw simulation data.simulation snapshots with density profiles and energy conservation.

When Omega = 1 (no correction): Energy NOT conserved -> Artificial collapse

""""""



import numpy as npThe grad-h correction term Ω accounts for variable smoothing length h(ρ):

import pandas as pd

import json    Ω_i = [1 + (h_i / D·ρ_i) Σ_j m_j ∂W_ij/∂h_i]^{-1}import numpy as np

import os

import sysimport json

import glob

import argparseWhen Ω ≠ 1 (with correction): Energy is conserved → Cloud stablefrom pathlib import Path



When Ω = 1 (no correction): Energy NOT conserved → Artificial collapseimport sys

def load_snapshot(filepath):

    """Load a snapshot CSV file, skipping comment lines.""""""

    return pd.read_csv(filepath, comment="#")

# Add shared module path



def load_energy_file(filepath):import numpy as npscript_dir = Path(__file__).parent

    """Load energy.dat file."""

    return pd.read_csv(filepath, comment="#", sep=r'\s+',import pandas as pdshared_path = script_dir.parent.parent.parent / "scripts" / "shared"

                       names=['time', 'kinetic', 'thermal', 'potential', 'total'])

import jsonif shared_path.exists():



def compute_radial_profile(df, n_bins=30):import os    sys.path.insert(0, str(shared_path))

    """Compute radial density profile with statistics."""

    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +import sys    from lane_emden import load_lane_emden_solution, get_density_profile

                (df['pos_z']**2 if 'pos_z' in df.columns else 0))

    r_max = r.max()from pathlib import Path    HAS_LANE_EMDEN = True

    bins = np.linspace(0, r_max * 1.05, n_bins + 1)

    r_centers = 0.5 * (bins[:-1] + bins[1:])import argparseelse:



    rho_mean = np.zeros(n_bins)    HAS_LANE_EMDEN = False

    rho_std = np.zeros(n_bins)

    gradh_mean = np.zeros(n_bins)



    for i in range(n_bins):def load_snapshot(filepath):

        mask = (r >= bins[i]) & (r < bins[i + 1])

        if mask.sum() > 0:    """Load a snapshot CSV file, skipping comment lines."""def load_snapshot(filepath):

            rho_mean[i] = df.loc[mask, 'dens'].mean()

            rho_std[i] = df.loc[mask, 'dens'].std()    return pd.read_csv(filepath, comment="#")    """Load a single CSV snapshot."""

            if 'gradh' in df.columns:

                gradh_mean[i] = df.loc[mask, 'gradh'].mean()    try:



    return r_centers, rho_mean, rho_std, gradh_mean        with open(filepath, 'r') as f:



def load_energy_file(filepath):            skiprows = 0

def compute_statistics(df):

    """Compute key statistics from a snapshot."""    """Load energy.dat file."""            time_code = 0.0

    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +

                (df['pos_z']**2 if 'pos_z' in df.columns else 0))    return pd.read_csv(filepath, comment="#", delim_whitespace=True,            for line in f:



    # Central density (inner 10%)                       names=['time', 'kinetic', 'thermal', 'potential', 'total'])                if line.startswith('# Time (code):'):

    r_max = r.max()

    central_mask = r < 0.1 * r_max                    time_code = float(line.split(':')[1].strip())

    rho_central = df.loc[central_mask, 'dens'].mean() if central_mask.sum() > 0 else df['dens'].max()

                if line.startswith('id,'):

    # Grad-h statistics

    gradh_mean = df['gradh'].mean() if 'gradh' in df.columns else 1.0def compute_radial_profile(df, n_bins=30):                    break

    gradh_min = df['gradh'].min() if 'gradh' in df.columns else 1.0

    gradh_max = df['gradh'].max() if 'gradh' in df.columns else 1.0    """Compute radial density profile with statistics."""                skiprows += 1



    # Energy from phi column    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +         

    E_kinetic = 0.5 * (df['mass'] * (df['vel_x']**2 + df['vel_y']**2 +

                                      (df['vel_z']**2 if 'vel_z' in df.columns else 0))).sum()                (df['pos_z']**2 if 'pos_z' in df.columns else 0))        data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows+1)

    E_potential = 0.5 * (df['mass'] * df['phi']).sum()

    E_internal = (df['mass'] * df['ene']).sum()    r_max = r.max()        return {

    E_total = E_kinetic + E_potential + E_internal

    bins = np.linspace(0, r_max * 1.05, n_bins + 1)            'time': time_code,

    # Total momentum

    px = (df['mass'] * df['vel_x']).sum()    r_centers = 0.5 * (bins[:-1] + bins[1:])            'pos': data[:, 1:4],

    py = (df['mass'] * df['vel_y']).sum()

    pz = (df['mass'] * df['vel_z']).sum() if 'vel_z' in df.columns else 0                'vel': data[:, 4:7],

    p_total = np.sqrt(px**2 + py**2 + pz**2)

    rho_mean = np.zeros(n_bins)            'mass': data[:, 10],

    return {

        'rho_central': float(rho_central),    rho_std = np.zeros(n_bins)            'rho': data[:, 11],

        'rho_max': float(df['dens'].max()),

        'rho_mean': float(df['dens'].mean()),    gradh_mean = np.zeros(n_bins)            'P': data[:, 12],

        'gradh_mean': float(gradh_mean),

        'gradh_min': float(gradh_min),                'u': data[:, 13],

        'gradh_max': float(gradh_max),

        'E_kinetic': float(E_kinetic),    for i in range(n_bins):            'h': data[:, 14]

        'E_potential': float(E_potential),

        'E_internal': float(E_internal),        mask = (r >= bins[i]) & (r < bins[i+1])        }

        'E_total': float(E_total),

        'p_total': float(p_total),        if mask.sum() > 0:    except Exception as e:

    }

            rho_mean[i] = df.loc[mask, 'dens'].mean()        print(f"Error loading {filepath}: {e}")



def process_case(results_dir, case_name, sample_particles=500):            rho_std[i] = df.loc[mask, 'dens'].std()        return None

    """Process all snapshots from a case directory."""

    snapshots = sorted(glob.glob(os.path.join(results_dir, "snapshot_*.csv")))            if 'gradh' in df.columns:

    if not snapshots:

        print(f"Warning: No snapshots found in {results_dir}")                gradh_mean[i] = df.loc[mask, 'gradh'].mean()

        return None

    def load_energy_file(filepath):

    # Load energy file if available

    energy_file = os.path.join(results_dir, "energy.dat")    return r_centers, rho_mean, rho_std, gradh_mean    """Load energy.dat file."""

    energy_data = None

    if os.path.exists(energy_file):    try:

        try:

            energy_data = load_energy_file(energy_file)        data = []

        except Exception as e:

            print(f"Warning: Could not load energy file: {e}")def compute_statistics(df):        with open(filepath, 'r') as f:



    frames = []    """Compute key statistics from a snapshot."""            for line in f:

    df = None

    for i, snap_file in enumerate(snapshots):    r = np.sqrt(df['pos_x']**2 + df['pos_y']**2 +                 if line.strip() and not line.startswith('#'):

        print(f"  Processing snapshot {i + 1}/{len(snapshots)}: {os.path.basename(snap_file)}")

        df = load_snapshot(snap_file)                (df['pos_z']**2 if 'pos_z' in df.columns else 0))                    values = line.split()



        # Compute statistics                        if len(values) >= 5:

        stats = compute_statistics(df)

    # Central density (inner 10%)                        data.append([float(x) for x in values[:5]])

        # Compute radial profile

        r_centers, rho_mean, rho_std, gradh_mean = compute_radial_profile(df)    r_max = r.max()        if not data:



        # Sample particles for scatter plot    central_mask = r < 0.1 * r_max            return None

        indices = np.random.choice(len(df), min(sample_particles, len(df)), replace=False)

        sampled = df.iloc[indices]    rho_central = df.loc[central_mask, 'dens'].mean() if central_mask.sum() > 0 else df['dens'].max()        data = np.array(data)



        frame = {            return {

            'index': i,

            'time': i * 2.0,  # Assuming outputTime = 2.0    # Grad-h statistics            'time': data[:, 0].tolist(),

            'stats': stats,

            'radial': {    gradh_mean = df['gradh'].mean() if 'gradh' in df.columns else 1.0            'kinetic': data[:, 1].tolist(),

                'r': r_centers.tolist(),

                'rho': rho_mean.tolist(),    gradh_min = df['gradh'].min() if 'gradh' in df.columns else 1.0            'thermal': data[:, 2].tolist(),

                'rho_std': rho_std.tolist(),

                'gradh': gradh_mean.tolist(),    gradh_max = df['gradh'].max() if 'gradh' in df.columns else 1.0            'potential': data[:, 3].tolist(),

            },

            'particles': {                'total': data[:, 4].tolist()

                'x': sampled['pos_x'].tolist(),

                'y': sampled['pos_y'].tolist(),    # Energy from phi column        }

                'rho': sampled['dens'].tolist(),

                'gradh': sampled['gradh'].tolist() if 'gradh' in sampled.columns else [1.0] * len(sampled),    E_kinetic = 0.5 * (df['mass'] * (df['vel_x']**2 + df['vel_y']**2 +     except Exception as e:

            }

        }                       (df['vel_z']**2 if 'vel_z' in df.columns else 0))).sum()        print(f"Error loading energy file: {e}")

        frames.append(frame)

    E_potential = 0.5 * (df['mass'] * df['phi']).sum()        return None

    # Use energy.dat for time series if available

    if energy_data is not None:    E_internal = (df['mass'] * df['ene']).sum()

        energy_series = {

            'time': energy_data['time'].tolist(),    E_total = E_kinetic + E_potential + E_internal

            'kinetic': energy_data['kinetic'].tolist(),

            'thermal': energy_data['thermal'].tolist(),    def compute_radial_profile(pos, rho, n_bins=50):

            'potential': energy_data['potential'].tolist(),

            'total': energy_data['total'].tolist(),    # Total momentum    """Compute radial density profile."""

        }

    else:    px = (df['mass'] * df['vel_x']).sum()    r = np.sqrt(np.sum(pos**2, axis=1))

        # Fall back to computed values

        energy_series = {    py = (df['mass'] * df['vel_y']).sum()    r_max = np.percentile(r, 99)

            'time': [f['time'] for f in frames],

            'kinetic': [f['stats']['E_kinetic'] for f in frames],    pz = (df['mass'] * df['vel_z']).sum() if 'vel_z' in df.columns else 0    bins = np.linspace(0, r_max, n_bins + 1)

            'thermal': [f['stats']['E_internal'] for f in frames],

            'potential': [f['stats']['E_potential'] for f in frames],    p_total = np.sqrt(px**2 + py**2 + pz**2)    r_centers = 0.5 * (bins[:-1] + bins[1:])

            'total': [f['stats']['E_total'] for f in frames],

        }        



    return {    return {    rho_profile = np.zeros(n_bins)

        'name': case_name,

        'frames': frames,        'rho_central': float(rho_central),    for i in range(n_bins):

        'energy': energy_series,

        'n_particles': len(df) if df is not None else 0,        'rho_max': float(df['dens'].max()),        mask = (r >= bins[i]) & (r < bins[i+1])

    }

        'rho_mean': float(df['dens'].mean()),        if mask.sum() > 0:



def generate_html(data_with_gradh, data_no_gradh, output_path):        'gradh_mean': float(gradh_mean),            rho_profile[i] = rho[mask].mean()

    """Generate interactive HTML viewer."""

        'gradh_min': float(gradh_min),    

    html_content = '''<!DOCTYPE html>

<html lang="en">        'gradh_max': float(gradh_max),    return r_centers.tolist(), rho_profile.tolist()

<head>

    <meta charset="UTF-8">        'E_kinetic': float(E_kinetic),

    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <title>Grad-H Effect: First Principles Analysis</title>        'E_potential': float(E_potential),

    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>

    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>        'E_internal': float(E_internal),def get_lane_emden_profile(rho_c, R, n_points=100):

    <style>

        * { margin: 0; padding: 0; box-sizing: border-box; }        'E_total': float(E_total),    """Get theoretical Lane-Emden profile."""

        body {

            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;        'p_total': float(p_total),    if HAS_LANE_EMDEN:

            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);

            color: #e8e8e8;    }        try:

            min-height: 100vh;

            padding: 20px;            solution = load_lane_emden_solution(n=1.5, dim=3)

        }

        h1 {            r, rho = get_density_profile(solution, rho_c, R, n_points=n_points)

            text-align: center;

            color: #ffd700;def process_case(results_dir, case_name, sample_particles=500):            return r.tolist(), rho.tolist()

            margin-bottom: 10px;

            font-size: 2em;    """Process all snapshots from a case directory."""        except Exception:

            text-shadow: 2px 2px 4px rgba(0,0,0,0.5);

        }    import glob            pass

        .subtitle {

            text-align: center;        return None, None

            color: #aaa;

            margin-bottom: 20px;    snapshots = sorted(glob.glob(os.path.join(results_dir, "snapshot_*.csv")))

            font-size: 1.1em;

        }    if not snapshots:

        .container {

            display: grid;        print(f"Warning: No snapshots found in {results_dir}")def generate_html(results_dirs, output_path, title="Self-Gravitating Cloud Simulation"):

            grid-template-columns: 1fr 1fr;

            gap: 20px;        return None    """Generate interactive HTML viewer."""

            max-width: 1800px;

            margin: 0 auto;        

        }

        .panel {    # Load energy file if available    all_data = {}

            background: rgba(30, 30, 50, 0.9);

            border-radius: 12px;    energy_file = os.path.join(results_dir, "energy.dat")    

            padding: 15px;

            box-shadow: 0 4px 20px rgba(0,0,0,0.4);    energy_data = None    for name, results_dir in results_dirs.items():

            border: 1px solid rgba(100,100,150,0.3);

        }    if os.path.exists(energy_file):        results_path = Path(results_dir)

        .panel h3 {

            color: #00ff88;        try:        if not results_path.exists():

            margin-bottom: 10px;

            font-size: 1.1em;            energy_data = load_energy_file(energy_file)            print(f"Warning: {results_dir} does not exist, skipping")

            border-bottom: 1px solid rgba(100,100,150,0.3);

            padding-bottom: 5px;        except:            continue

        }

        .panel.with-gradh h3 { color: #00ff88; }            pass            

        .panel.no-gradh h3 { color: #ff4444; }

        .panel.theory { grid-column: span 2; }            snapshots = sorted(results_path.glob('snapshot_*.csv'))

        .controls {

            display: flex;    frames = []        if not snapshots:

            justify-content: center;

            align-items: center;    for i, snap_file in enumerate(snapshots):            print(f"Warning: No snapshots in {results_dir}, skipping")

            gap: 20px;

            margin: 20px 0;        df = load_snapshot(snap_file)            continue

            padding: 15px;

            background: rgba(30, 30, 50, 0.9);                

            border-radius: 12px;

            flex-wrap: wrap;        # Compute statistics        print(f"Processing {name}: {len(snapshots)} snapshots")

        }

        button {        stats = compute_statistics(df)        

            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);

            color: white;                frames = []

            border: none;

            padding: 10px 20px;        # Compute radial profile        for snap_path in snapshots:

            border-radius: 8px;

            cursor: pointer;        r_centers, rho_mean, rho_std, gradh_mean = compute_radial_profile(df)            snap = load_snapshot(snap_path)

            font-size: 1em;

            transition: transform 0.2s, box-shadow 0.2s;                    if snap is None:

        }

        button:hover {        # Sample particles for scatter plot (stratified by radius)                continue

            transform: translateY(-2px);

            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);        r = np.sqrt(df['pos_x']**2 + df['pos_y']**2)            

        }

        button:active { transform: translateY(0); }        indices = np.random.choice(len(df), min(sample_particles, len(df)), replace=False)            r_profile, rho_profile = compute_radial_profile(snap['pos'], snap['rho'])

        #slider {

            width: 400px;        sampled = df.iloc[indices]            

            height: 8px;

            -webkit-appearance: none;                    # Compute stats

            background: linear-gradient(90deg, #667eea, #764ba2);

            border-radius: 4px;        frame = {            r = np.sqrt(np.sum(snap['pos']**2, axis=1))

            outline: none;

        }            'index': i,            rho_max = snap['rho'].max()

        #slider::-webkit-slider-thumb {

            -webkit-appearance: none;            'time': i * 2.0,  # Assuming outputTime = 2.0            rho_center = snap['rho'][r < 0.1].mean() if (r < 0.1).sum() > 0 else rho_max

            width: 20px;

            height: 20px;            'stats': stats,            

            background: #ffd700;

            border-radius: 50%;            'radial': {            frames.append({

            cursor: pointer;

            box-shadow: 0 2px 10px rgba(255,215,0,0.5);                'r': r_centers.tolist(),                'time': round(snap['time'], 2),

        }

        .info-box {                'rho': rho_mean.tolist(),                'r_profile': r_profile,

            display: flex;

            gap: 30px;                'rho_std': rho_std.tolist(),                'rho_profile': rho_profile,

            justify-content: center;

            color: #ccc;                'gradh': gradh_mean.tolist(),                'rho_max': round(float(rho_max), 4),

            font-size: 1.1em;

        }            },                'rho_center': round(float(rho_center), 4),

        .info-box span { color: #ffd700; font-weight: bold; }

        .stats-grid {            'particles': {                'r_max': round(float(r.max()), 4)

            display: grid;

            grid-template-columns: repeat(3, 1fr);                'x': sampled['pos_x'].tolist(),            })

            gap: 10px;

            margin-top: 10px;                'y': sampled['pos_y'].tolist(),        

        }

        .stat-item {                'rho': sampled['dens'].tolist(),        # Load energy data

            background: rgba(50, 50, 80, 0.5);

            padding: 8px;                'gradh': sampled['gradh'].tolist() if 'gradh' in sampled.columns else [1.0] * len(sampled),        energy_file = results_path / 'energy.dat'

            border-radius: 6px;

            text-align: center;            }        energy = load_energy_file(energy_file) if energy_file.exists() else None

        }

        .stat-label { color: #888; font-size: 0.85em; }        }        

        .stat-value { color: #fff; font-size: 1.1em; font-weight: bold; }

        .stat-value.good { color: #00ff88; }        frames.append(frame)        # Get initial parameters for Lane-Emden

        .stat-value.bad { color: #ff4444; }

        .theory-content {            if frames:

            display: grid;

            grid-template-columns: 1fr 1fr;    # Use energy.dat for time series if available            initial_rho_c = frames[0]['rho_center']

            gap: 20px;

        }    if energy_data is not None:            initial_R = frames[0]['r_max']

        .derivation {

            padding: 15px;        energy_series = {            le_r, le_rho = get_lane_emden_profile(initial_rho_c, initial_R)

            background: rgba(40, 40, 70, 0.5);

            border-radius: 8px;            'time': energy_data['time'].tolist(),        else:

            line-height: 1.8;

        }            'kinetic': energy_data['kinetic'].tolist(),            le_r, le_rho = None, None

        .derivation h4 {

            color: #ffd700;            'thermal': energy_data['thermal'].tolist(),        

            margin-bottom: 10px;

        }            'potential': energy_data['potential'].tolist(),        all_data[name] = {

        .equation {

            background: rgba(60, 60, 100, 0.5);            'total': energy_data['total'].tolist(),            'frames': frames,

            padding: 10px;

            border-radius: 6px;        }            'energy': energy,

            margin: 8px 0;

            text-align: center;    else:            'lane_emden_r': le_r,

            font-size: 1.1em;

        }        # Fall back to computed values            'lane_emden_rho': le_rho

        .key-insight {

            background: linear-gradient(135deg, rgba(0,255,136,0.1) 0%, rgba(0,200,100,0.1) 100%);        energy_series = {        }

            border-left: 3px solid #00ff88;

            padding: 10px 15px;            'time': [f['time'] for f in frames],    

            margin: 10px 0;

            border-radius: 0 6px 6px 0;            'kinetic': [f['stats']['E_kinetic'] for f in frames],    if not all_data:

        }

        .warning {            'thermal': [f['stats']['E_internal'] for f in frames],        print("Error: No data found!")

            background: linear-gradient(135deg, rgba(255,68,68,0.1) 0%, rgba(200,50,50,0.1) 100%);

            border-left: 3px solid #ff4444;            'potential': [f['stats']['E_potential'] for f in frames],        return

            padding: 10px 15px;

            margin: 10px 0;            'total': [f['stats']['E_total'] for f in frames],    

            border-radius: 0 6px 6px 0;

        }        }    # Generate HTML

        .plot-container { height: 350px; }

        .full-width { grid-column: span 2; }        html_content = f'''<!DOCTYPE html>

        @media (max-width: 1200px) {

            .container { grid-template-columns: 1fr; }    return {<html lang="en">

            .panel.theory, .full-width { grid-column: span 1; }

            .theory-content { grid-template-columns: 1fr; }        'name': case_name,<head>

        }

    </style>        'frames': frames,    <meta charset="UTF-8">

</head>

<body>        'energy': energy_series,    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <h1>Grad-H Correction Effect on Self-Gravitating Cloud</h1>

    <p class="subtitle">First Principles Analysis: How Variable Smoothing Length Affects Energy Conservation</p>        'n_particles': len(df),    <title>{title}</title>



    <div class="controls">    }    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>

        <button onclick="prevFrame()">&#9664; Previous</button>

        <button onclick="togglePlay()" id="playBtn">&#9654; Play</button>    <style>

        <button onclick="nextFrame()">Next &#9654;</button>

        <input type="range" id="slider" min="0" max="10" value="0" oninput="goToFrame(this.value)">        * {{ margin: 0; padding: 0; box-sizing: border-box; }}

        <div class="info-box">

            <div>Frame: <span id="frameNum">1</span>/<span id="totalFrames">11</span></div>def generate_html(data_with_gradh, data_no_gradh, output_path):        body {{

            <div>Time: <span id="timeVal">0.0</span></div>

        </div>    """Generate interactive HTML viewer."""            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;

    </div>

                background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);

    <div class="container">

        <div class="panel with-gradh">    html_template = '''<!DOCTYPE html>            color: #e0e0e0;

            <h3>WITH Grad-H (Omega != 1) - Stable</h3>

            <div id="scatter-with" class="plot-container"></div><html lang="en">            min-height: 100vh;

            <div class="stats-grid">

                <div class="stat-item"><head>            padding: 20px;

                    <div class="stat-label">rho_central</div>

                    <div class="stat-value good" id="rho-central-with">-</div>    <meta charset="UTF-8">        }}

                </div>

                <div class="stat-item">    <meta name="viewport" content="width=device-width, initial-scale=1.0">        .container {{ max-width: 1600px; margin: 0 auto; }}

                    <div class="stat-label">Omega_mean</div>

                    <div class="stat-value good" id="gradh-mean-with">-</div>    <title>Grad-H Effect: First Principles Analysis</title>        h1 {{

                </div>

                <div class="stat-item">    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>            text-align: center;

                    <div class="stat-label">Omega_range</div>

                    <div class="stat-value" id="gradh-range-with">-</div>    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>            color: #00d4ff;

                </div>

            </div>    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>            margin-bottom: 20px;

        </div>

    <style>            font-size: 2em;

        <div class="panel no-gradh">

            <h3>NO Grad-H (Omega = 1) - Collapsing</h3>        * { margin: 0; padding: 0; box-sizing: border-box; }            text-shadow: 0 0 20px rgba(0, 212, 255, 0.5);

            <div id="scatter-no" class="plot-container"></div>

            <div class="stats-grid">        body {        }}

                <div class="stat-item">

                    <div class="stat-label">rho_central</div>            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;        .controls {{

                    <div class="stat-value bad" id="rho-central-no">-</div>

                </div>            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);            background: rgba(255,255,255,0.05);

                <div class="stat-item">

                    <div class="stat-label">Omega (forced)</div>            color: #e8e8e8;            padding: 20px;

                    <div class="stat-value bad" id="gradh-mean-no">1.00</div>

                </div>            min-height: 100vh;            border-radius: 15px;

                <div class="stat-item">

                    <div class="stat-label">rho increase</div>            padding: 20px;            margin-bottom: 20px;

                    <div class="stat-value bad" id="rho-increase">-</div>

                </div>        }            display: flex;

            </div>

        </div>        h1 {            flex-wrap: wrap;



        <div class="panel">            text-align: center;            gap: 20px;

            <h3>Radial Density Profile</h3>

            <div id="radial-profile" class="plot-container"></div>            color: #ffd700;            align-items: center;

        </div>

            margin-bottom: 10px;            justify-content: center;

        <div class="panel">

            <h3>Grad-H Term (Omega) Distribution</h3>            font-size: 2em;        }}

            <div id="gradh-distribution" class="plot-container"></div>

        </div>            text-shadow: 2px 2px 4px rgba(0,0,0,0.5);        .control-group {{



        <div class="panel full-width">        }            display: flex;

            <h3>Energy Conservation Comparison</h3>

            <div id="energy-plot" style="height: 300px;"></div>        .subtitle {            align-items: center;

        </div>

            text-align: center;            gap: 10px;

        <div class="panel theory">

            <h3>First Principles: Why Grad-H Matters</h3>            color: #aaa;        }}

            <div class="theory-content">

                <div class="derivation">            margin-bottom: 20px;        label {{ color: #aaa; font-weight: 500; }}

                    <h4>The Grad-H Correction Term</h4>

                    <p>In SPH, smoothing length adapts to density: \\( h_i = \\eta (m_i/\\rho_i)^{1/D} \\)</p>            font-size: 1.1em;        select, button {{

                    <p>This creates an implicit dependence: \\( \\rho_i \\leftrightarrow h_i \\)</p>

        }            padding: 10px 20px;

                    <div class="equation">

                        \\( \\rho_i = \\sum_j m_j W(|\\mathbf{r}_i - \\mathbf{r}_j|, h_i) \\)        .container {            border: none;

                    </div>

            display: grid;            border-radius: 8px;

                    <p>Taking the derivative properly (chain rule):</p>

                    <div class="equation">            grid-template-columns: 1fr 1fr;            font-size: 14px;

                        \\( \\frac{\\partial \\rho_k}{\\partial \\mathbf{x}_i} = \\Omega_k \\sum_j m_j \\nabla_i W_{kj} \\)

                    </div>            gap: 20px;            cursor: pointer;



                    <div class="key-insight">            max-width: 1800px;            transition: all 0.3s;

                        <strong>The Omega factor:</strong><br>

                        \\( \\displaystyle \\Omega_i = \\left[1 + \\frac{h_i}{D\\rho_i}\\sum_j m_j \\frac{\\partial W_{ij}}{\\partial h_i}\\right]^{-1} \\)            margin: 0 auto;        }}

                    </div>

                </div>        }        select {{



                <div class="derivation">        .panel {            background: #2a2a4a;

                    <h4>Effect on Gravitational Potential</h4>

                    <p>Hernquist-Katz softening uses: \\( \\epsilon = h/2 \\)</p>            background: rgba(30, 30, 50, 0.9);            color: #fff;

                    <p>So gravitational potential \\( \\phi_i \\) depends on \\( h_i(\\rho_i) \\)!</p>

            border-radius: 12px;            border: 1px solid #444;

                    <div class="warning">

                        <strong>Without Omega correction (Omega = 1):</strong><br>            padding: 15px;        }}

                        - Ignores \\( \\partial\\phi/\\partial h \\) contribution<br>

                        - Gravitational energy NOT conserved<br>            box-shadow: 0 4px 20px rgba(0,0,0,0.4);        button {{

                        - Artificial numerical collapse occurs

                    </div>            border: 1px solid rgba(100,100,150,0.3);            background: linear-gradient(135deg, #00d4ff, #0099cc);



                    <div class="key-insight">        }            color: #000;

                        <strong>With Omega correction (Omega != 1):</strong><br>

                        - Properly accounts for h(rho) in potential<br>        .panel h3 {            font-weight: bold;

                        - Total energy exactly conserved<br>

                        - Cloud remains in hydrostatic equilibrium            color: #00ff88;        }}

                    </div>

            margin-bottom: 10px;        button:hover {{ transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0, 212, 255, 0.4); }}

                    <p style="margin-top: 15px; color: #ffd700;">

                        <strong>Physical insight:</strong> The correction ensures that the pressure-density            font-size: 1.1em;        button.active {{ background: linear-gradient(135deg, #ff6b6b, #cc5555); }}

                        relation from the Lagrangian is thermodynamically consistent.

                    </p>            border-bottom: 1px solid rgba(100,100,150,0.3);        .slider-container {{

                </div>

            </div>            padding-bottom: 5px;            display: flex;

        </div>

    </div>        }            align-items: center;



    <script>        .panel.with-gradh h3 { color: #00ff88; }            gap: 15px;

    var dataWith = ''' + json.dumps(data_with_gradh) + ''';

    var dataNo = ''' + json.dumps(data_no_gradh) + ''';        .panel.no-gradh h3 { color: #ff4444; }            flex: 1;



    var currentFrame = 0;        .panel.theory { grid-column: span 2; }            min-width: 300px;

    var playing = false;

    var playInterval = null;        .controls {        }}

    var totalFrames = dataWith.frames.length;

            display: flex;        input[type="range"] {{

    document.getElementById('totalFrames').textContent = totalFrames;

    document.getElementById('slider').max = totalFrames - 1;            justify-content: center;            flex: 1;



    var colorscale = [            align-items: center;            height: 8px;

        [0, '#000033'],

        [0.2, '#003366'],            gap: 20px;            border-radius: 4px;

        [0.4, '#006699'],

        [0.6, '#ff9900'],            margin: 20px 0;            background: #333;

        [0.8, '#ff3300'],

        [1, '#ffff00']            padding: 15px;            outline: none;

    ];

            background: rgba(30, 30, 50, 0.9);            -webkit-appearance: none;

    var rho0_with = dataWith.frames[0].stats.rho_central;

    var rho0_no = dataNo.frames[0].stats.rho_central;            border-radius: 12px;        }}

    var E0_with = Math.abs(dataWith.energy.total[0]);

    var E0_no = Math.abs(dataNo.energy.total[0]);        }        input[type="range"]::-webkit-slider-thumb {{



    function updatePlots() {        button {            -webkit-appearance: none;

        var fw = dataWith.frames[currentFrame];

        var fn = dataNo.frames[currentFrame];            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);            width: 20px;



        document.getElementById('frameNum').textContent = currentFrame + 1;            color: white;            height: 20px;

        document.getElementById('timeVal').textContent = fw.time.toFixed(1);

            border: none;            border-radius: 50%;

        document.getElementById('rho-central-with').textContent = fw.stats.rho_central.toFixed(3);

        document.getElementById('gradh-mean-with').textContent = fw.stats.gradh_mean.toFixed(3);            padding: 10px 20px;            background: #00d4ff;

        document.getElementById('gradh-range-with').textContent =

            '[' + fw.stats.gradh_min.toFixed(2) + ', ' + fw.stats.gradh_max.toFixed(2) + ']';            border-radius: 8px;            cursor: pointer;



        document.getElementById('rho-central-no').textContent = fn.stats.rho_central.toFixed(3);            cursor: pointer;            box-shadow: 0 0 10px rgba(0, 212, 255, 0.5);

        var rhoIncrease = ((fn.stats.rho_central / rho0_no - 1) * 100).toFixed(1);

        document.getElementById('rho-increase').textContent = '+' + rhoIncrease + '%';            font-size: 1em;        }}



        var allRho = fw.particles.rho.concat(fn.particles.rho);            transition: transform 0.2s, box-shadow 0.2s;        .time-display {{

        var rhoMin = Math.min.apply(null, allRho);

        var rhoMax = Math.max.apply(null, allRho);        }            background: #2a2a4a;



        Plotly.react('scatter-with', [{        button:hover {            padding: 8px 16px;

            x: fw.particles.x,

            y: fw.particles.y,            transform: translateY(-2px);            border-radius: 8px;

            mode: 'markers',

            type: 'scatter',            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);            font-family: monospace;

            marker: {

                size: 4,        }            font-size: 16px;

                color: fw.particles.rho,

                colorscale: colorscale,        button:active { transform: translateY(0); }            color: #00ff88;

                cmin: rhoMin,

                cmax: rhoMax,        #slider {            min-width: 120px;

                colorbar: { title: 'rho', thickness: 15 }

            },            width: 400px;            text-align: center;

            hovertemplate: 'x: %{x:.3f}<br>y: %{y:.3f}<br>rho: %{marker.color:.3f}<extra></extra>'

        }], {            height: 8px;        }}

            paper_bgcolor: 'rgba(0,0,0,0)',

            plot_bgcolor: 'rgba(30,30,50,0.5)',            -webkit-appearance: none;        .plots {{

            xaxis: { title: 'x', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },

            yaxis: { title: 'y', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)', scaleanchor: 'x' },            background: linear-gradient(90deg, #667eea, #764ba2);            display: grid;

            margin: { t: 10, b: 40, l: 50, r: 80 },

            font: { color: '#ccc' }            border-radius: 4px;            grid-template-columns: 1fr 1fr;

        }, { responsive: true });

            outline: none;            gap: 20px;

        Plotly.react('scatter-no', [{

            x: fn.particles.x,        }        }}

            y: fn.particles.y,

            mode: 'markers',        #slider::-webkit-slider-thumb {        .plot-container {{

            type: 'scatter',

            marker: {            -webkit-appearance: none;            background: rgba(0,0,0,0.3);

                size: 4,

                color: fn.particles.rho,            width: 20px;            border-radius: 15px;

                colorscale: colorscale,

                cmin: rhoMin,            height: 20px;            padding: 10px;

                cmax: rhoMax,

                colorbar: { title: 'rho', thickness: 15 }            background: #ffd700;            box-shadow: 0 10px 30px rgba(0,0,0,0.3);

            },

            hovertemplate: 'x: %{x:.3f}<br>y: %{y:.3f}<br>rho: %{marker.color:.3f}<extra></extra>'            border-radius: 50%;        }}

        }], {

            paper_bgcolor: 'rgba(0,0,0,0)',            cursor: pointer;        .stats {{

            plot_bgcolor: 'rgba(30,30,50,0.5)',

            xaxis: { title: 'x', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },            box-shadow: 0 2px 10px rgba(255,215,0,0.5);            display: grid;

            yaxis: { title: 'y', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)', scaleanchor: 'x' },

            margin: { t: 10, b: 40, l: 50, r: 80 },        }            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));

            font: { color: '#ccc' }

        }, { responsive: true });        .info-box {            gap: 15px;



        Plotly.react('radial-profile', [            display: flex;            margin-top: 20px;

            {

                x: fw.radial.r,            gap: 30px;        }}

                y: fw.radial.rho,

                mode: 'lines',            justify-content: center;        .stat-card {{

                name: 'With Omega',

                line: { color: '#00ff88', width: 3 }            color: #ccc;            background: rgba(255,255,255,0.05);

            },

            {            font-size: 1.1em;            padding: 15px;

                x: fn.radial.r,

                y: fn.radial.rho,        }            border-radius: 10px;

                mode: 'lines',

                name: 'No Omega',        .info-box span { color: #ffd700; font-weight: bold; }            text-align: center;

                line: { color: '#ff4444', width: 3, dash: 'dash' }

            }        .stats-grid {        }}

        ], {

            paper_bgcolor: 'rgba(0,0,0,0)',            display: grid;        .stat-value {{

            plot_bgcolor: 'rgba(30,30,50,0.5)',

            xaxis: { title: 'Radius r', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },            grid-template-columns: repeat(3, 1fr);            font-size: 1.8em;

            yaxis: { title: 'Density rho', type: 'log', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },

            legend: { x: 0.7, y: 0.95, bgcolor: 'rgba(30,30,50,0.8)' },            gap: 10px;            font-weight: bold;

            margin: { t: 10, b: 40, l: 60, r: 20 },

            font: { color: '#ccc' }            margin-top: 10px;            color: #00d4ff;

        }, { responsive: true });

        }        }}

        Plotly.react('gradh-distribution', [

            {        .stat-item {        .stat-label {{ color: #888; font-size: 0.9em; }}

                x: fw.radial.r,

                y: fw.radial.gradh,            background: rgba(50, 50, 80, 0.5);        .info {{

                mode: 'lines+markers',

                name: 'Omega(r) with correction',            padding: 8px;            background: rgba(0, 212, 255, 0.1);

                line: { color: '#00ff88', width: 2 },

                marker: { size: 6 }            border-radius: 6px;            border-left: 4px solid #00d4ff;

            },

            {            text-align: center;            padding: 15px;

                x: fn.radial.r,

                y: fn.radial.r.map(function() { return 1.0; }),        }            margin-top: 20px;

                mode: 'lines',

                name: 'Omega = 1 (forced)',        .stat-label { color: #888; font-size: 0.85em; }            border-radius: 0 10px 10px 0;

                line: { color: '#ff4444', width: 2, dash: 'dash' }

            }        .stat-value { color: #fff; font-size: 1.1em; font-weight: bold; }        }}

        ], {

            paper_bgcolor: 'rgba(0,0,0,0)',        .stat-value.good { color: #00ff88; }        .info h3 {{ color: #00d4ff; margin-bottom: 10px; }}

            plot_bgcolor: 'rgba(30,30,50,0.5)',

            xaxis: { title: 'Radius r', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },        .stat-value.bad { color: #ff4444; }        .info p {{ line-height: 1.6; color: #aaa; }}

            yaxis: { title: 'Omega (grad-h factor)', range: [0.8, 1.5], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },

            legend: { x: 0.6, y: 0.95, bgcolor: 'rgba(30,30,50,0.8)' },        .theory-content {    </style>

            margin: { t: 10, b: 40, l: 60, r: 20 },

            font: { color: '#ccc' },            display: grid;</head>

            shapes: [{

                type: 'line',            grid-template-columns: 1fr 1fr;<body>

                x0: 0, x1: 1.2,

                y0: 1, y1: 1,            gap: 20px;    <div class="container">

                line: { color: '#888', width: 1, dash: 'dot' }

            }]        }        <h1>🌌 {title}</h1>

        }, { responsive: true });

        .derivation {        

        var currentTime = fw.time;

        Plotly.react('energy-plot', [            padding: 15px;        <div class="controls">

            {

                x: dataWith.energy.time,            background: rgba(40, 40, 70, 0.5);            <div class="control-group">

                y: dataWith.energy.total.map(function(e) { return e / E0_with; }),

                mode: 'lines',            border-radius: 8px;                <label>Dataset:</label>

                name: 'E_total (with Omega)',

                line: { color: '#00ff88', width: 3 }            line-height: 1.8;                <select id="dataset-select" onchange="changeDataset()">

            },

            {        }                    {"".join(f'<option value="{name}">{name.replace("_", " ").title()}</option>' for name in all_data.keys())}

                x: dataNo.energy.time,

                y: dataNo.energy.total.map(function(e) { return e / E0_no; }),        .derivation h4 {                </select>

                mode: 'lines',

                name: 'E_total (no Omega)',            color: #ffd700;            </div>

                line: { color: '#ff4444', width: 3, dash: 'dash' }

            },            margin-bottom: 10px;            

            {

                x: dataWith.energy.time,        }            <div class="control-group">

                y: dataWith.energy.kinetic.map(function(e) { return e / E0_with; }),

                mode: 'lines',        .equation {                <button id="play-btn" onclick="togglePlay()">▶ Play</button>

                name: 'E_kin (with)',

                line: { color: '#4ecdc4', width: 1.5 },            background: rgba(60, 60, 100, 0.5);                <button onclick="stepBackward()">⏮ -1</button>

                visible: 'legendonly'

            },            padding: 10px;                <button onclick="stepForward()">⏭ +1</button>

            {

                x: dataWith.energy.time,            border-radius: 6px;            </div>

                y: dataWith.energy.potential.map(function(e) { return e / E0_with; }),

                mode: 'lines',            margin: 8px 0;            

                name: 'E_grav (with)',

                line: { color: '#ff6b6b', width: 1.5 },            text-align: center;            <div class="slider-container">

                visible: 'legendonly'

            }            font-size: 1.1em;                <label>Frame:</label>

        ], {

            paper_bgcolor: 'rgba(0,0,0,0)',        }                <input type="range" id="frame-slider" min="0" value="0" onchange="updateFrame(this.value)">

            plot_bgcolor: 'rgba(30,30,50,0.5)',

            xaxis: { title: 'Time', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },        .key-insight {                <div class="time-display" id="time-display">t = 0.00</div>

            yaxis: { title: 'Energy / |E0|', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },

            legend: { x: 0.7, y: 0.95, bgcolor: 'rgba(30,30,50,0.8)' },            background: linear-gradient(135deg, rgba(0,255,136,0.1) 0%, rgba(0,200,100,0.1) 100%);            </div>

            margin: { t: 10, b: 40, l: 60, r: 20 },

            font: { color: '#ccc' },            border-left: 3px solid #00ff88;            

            shapes: [{

                type: 'line',            padding: 10px 15px;            <div class="control-group">

                x0: currentTime, x1: currentTime,

                y0: -2.5, y1: 0.5,            margin: 10px 0;                <label>Speed:</label>

                line: { color: '#ffd700', width: 2, dash: 'dot' }

            }],            border-radius: 0 6px 6px 0;                <select id="speed-select">

            annotations: [{

                x: currentTime, y: 0.3,        }                    <option value="2000">0.5x</option>

                text: 't=' + currentTime.toFixed(1),

                showarrow: true,        .warning {                    <option value="1000" selected>1x</option>

                arrowhead: 2,

                arrowcolor: '#ffd700',            background: linear-gradient(135deg, rgba(255,68,68,0.1) 0%, rgba(200,50,50,0.1) 100%);                    <option value="500">2x</option>

                font: { color: '#ffd700' }

            }]            border-left: 3px solid #ff4444;                    <option value="200">5x</option>

        }, { responsive: true });

    }            padding: 10px 15px;                </select>



    function goToFrame(idx) {            margin: 10px 0;            </div>

        currentFrame = parseInt(idx);

        document.getElementById('slider').value = currentFrame;            border-radius: 0 6px 6px 0;        </div>

        updatePlots();

    }        }        



    function nextFrame() {        .plot-container { height: 350px; }        <div class="plots">

        currentFrame = (currentFrame + 1) % totalFrames;

        document.getElementById('slider').value = currentFrame;        .full-width { grid-column: span 2; }            <div class="plot-container">

        updatePlots();

    }    </style>                <div id="density-plot" style="height: 450px;"></div>



    function prevFrame() {</head>            </div>

        currentFrame = (currentFrame - 1 + totalFrames) % totalFrames;

        document.getElementById('slider').value = currentFrame;<body>            <div class="plot-container">

        updatePlots();

    }    <h1>🌌 Grad-H Correction (Ω) Effect on Self-Gravitating Cloud</h1>                <div id="energy-plot" style="height: 450px;"></div>



    function togglePlay() {    <p class="subtitle">First Principles Analysis: How Variable Smoothing Length Affects Energy Conservation</p>            </div>

        playing = !playing;

        var btn = document.getElementById('playBtn');            </div>

        if (playing) {

            btn.innerHTML = '&#9208; Pause';    <div class="controls">        

            playInterval = setInterval(nextFrame, 300);

        } else {        <button onclick="prevFrame()">◀ Previous</button>        <div class="stats">

            btn.innerHTML = '&#9654; Play';

            clearInterval(playInterval);        <button onclick="togglePlay()" id="playBtn">▶ Play</button>            <div class="stat-card">

        }

    }        <button onclick="nextFrame()">Next ▶</button>                <div class="stat-value" id="stat-time">0.00</div>



    updatePlots();        <input type="range" id="slider" min="0" max="10" value="0" oninput="goToFrame(this.value)">                <div class="stat-label">Simulation Time</div>

    </script>

</body>        <div class="info-box">            </div>

</html>

'''            <div>Frame: <span id="frameNum">1</span>/<span id="totalFrames">11</span></div>            <div class="stat-card">



    with open(output_path, 'w') as f:            <div>Time: <span id="timeVal">0.0</span></div>                <div class="stat-value" id="stat-rho-max">0.00</div>

        f.write(html_content)

        </div>                <div class="stat-label">Max Density</div>

    print(f"Interactive viewer generated: {output_path}")

    </div>            </div>



def main():                <div class="stat-card">

    parser = argparse.ArgumentParser(description="Generate interactive HTML viewer for grad-h study")

    parser.add_argument("--with-gradh", required=True, help="Directory with grad-h results")    <div class="container">                <div class="stat-value" id="stat-rho-center">0.00</div>

    parser.add_argument("--no-gradh", required=True, help="Directory without grad-h results")

    parser.add_argument("--output", required=True, help="Output HTML file path")        <!-- Scatter plots row -->                <div class="stat-label">Central Density</div>

    parser.add_argument("--sample", type=int, default=500, help="Number of particles to sample per frame")

        <div class="panel with-gradh">            </div>

    args = parser.parse_args()

            <h3>✅ WITH Grad-H (Ω ≠ 1) - Stable</h3>            <div class="stat-card">

    print("Loading data with grad-h correction...")

    data_with = process_case(args.with_gradh, "With Grad-H", args.sample)            <div id="scatter-with" class="plot-container"></div>                <div class="stat-value" id="stat-frame">0 / 0</div>



    print("Loading data without grad-h correction...")            <div class="stats-grid">                <div class="stat-label">Frame</div>

    data_no = process_case(args.no_gradh, "No Grad-H", args.sample)

                <div class="stat-item">            </div>

    if data_with and data_no:

        output_dir = os.path.dirname(args.output)                    <div class="stat-label">ρ_central</div>        </div>

        if output_dir:

            os.makedirs(output_dir, exist_ok=True)                    <div class="stat-value good" id="rho-central-with">-</div>        

        generate_html(data_with, data_no, args.output)

        print(f"\nDone! Open in browser: {args.output}")                </div>        <div class="info">

    else:

        print("Error: Could not load data from both directories")                <div class="stat-item">            <h3>ℹ️ About This Simulation</h3>

        sys.exit(1)

                    <div class="stat-label">Ω_mean</div>            <p>



if __name__ == "__main__":                    <div class="stat-value good" id="gradh-mean-with">-</div>                This visualization shows the evolution of a self-gravitating gas cloud using SPH (Smoothed Particle Hydrodynamics).

    main()

                </div>                The left plot shows the radial density profile compared to the theoretical Lane-Emden solution (cyan dashed line).

                <div class="stat-item">                The right plot shows energy conservation over time.

                    <div class="stat-label">Ω_range</div>                <br><br>

                    <div class="stat-value" id="gradh-range-with">-</div>                <strong>Key hypothesis:</strong> Wendland C4 kernel-convolved gravity should maintain better energy conservation

                </div>                without the grad-h (Ω) correction compared to traditional Hernquist-Katz softening.

            </div>            </p>

        </div>        </div>

            </div>

        <div class="panel no-gradh">    

            <h3>❌ NO Grad-H (Ω = 1) - Collapsing</h3>    <script>

            <div id="scatter-no" class="plot-container"></div>        // Embedded simulation data

            <div class="stats-grid">        const allData = {json.dumps(all_data)};

                <div class="stat-item">        

                    <div class="stat-label">ρ_central</div>        let currentDataset = Object.keys(allData)[0];

                    <div class="stat-value bad" id="rho-central-no">-</div>        let currentFrame = 0;

                </div>        let isPlaying = false;

                <div class="stat-item">        let playInterval = null;

                    <div class="stat-label">Ω (forced)</div>        

                    <div class="stat-value bad" id="gradh-mean-no">1.00</div>        // Initialize

                </div>        document.addEventListener('DOMContentLoaded', function() {{

                <div class="stat-item">            initPlots();

                    <div class="stat-label">ρ increase</div>            updateFrame(0);

                    <div class="stat-value bad" id="rho-increase">-</div>        }});

                </div>        

            </div>        function changeDataset() {{

        </div>            currentDataset = document.getElementById('dataset-select').value;

                    currentFrame = 0;

        <!-- Radial profile and Ω distribution -->            const slider = document.getElementById('frame-slider');

        <div class="panel">            slider.max = allData[currentDataset].frames.length - 1;

            <h3>📊 Radial Density Profile</h3>            slider.value = 0;

            <div id="radial-profile" class="plot-container"></div>            updateFrame(0);

        </div>            updateEnergyPlot();

                }}

        <div class="panel">        

            <h3>🔬 Grad-H Term (Ω) Distribution</h3>        function initPlots() {{

            <div id="gradh-distribution" class="plot-container"></div>            const slider = document.getElementById('frame-slider');

        </div>            slider.max = allData[currentDataset].frames.length - 1;

                    

        <!-- Energy conservation -->            // Density plot

        <div class="panel full-width">            Plotly.newPlot('density-plot', [], {{

            <h3>⚡ Energy Conservation Comparison</h3>                title: {{ text: 'Radial Density Profile', font: {{ color: '#fff' }} }},

            <div id="energy-plot" style="height: 300px;"></div>                paper_bgcolor: 'rgba(0,0,0,0)',

        </div>                plot_bgcolor: 'rgba(0,0,0,0.2)',

                        xaxis: {{ title: 'Radius r', color: '#aaa', gridcolor: '#333' }},

        <!-- Theory panel -->                yaxis: {{ title: 'Density ρ', type: 'log', color: '#aaa', gridcolor: '#333' }},

        <div class="panel theory">                legend: {{ font: {{ color: '#fff' }} }},

            <h3>📐 First Principles: Why Grad-H Matters</h3>                margin: {{ t: 50, r: 20, b: 50, l: 60 }}

            <div class="theory-content">            }});

                <div class="derivation">            

                    <h4>The Grad-H Correction Term</h4>            // Energy plot

                    <p>In SPH, smoothing length adapts to density: \\( h_i = \\eta (m_i/\\rho_i)^{1/D} \\)</p>            updateEnergyPlot();

                    <p>This creates an implicit dependence: \\( \\rho_i \\leftrightarrow h_i \\)</p>        }}

                            

                    <div class="equation">        function updateEnergyPlot() {{

                        \\( \\rho_i = \\sum_j m_j W(|\\mathbf{r}_i - \\mathbf{r}_j|, h_i) \\)            const energy = allData[currentDataset].energy;

                    </div>            if (!energy) {{

                                    Plotly.react('energy-plot', [{{

                    <p>Taking the derivative properly (chain rule):</p>                    x: [0], y: [0], type: 'scatter', name: 'No energy data'

                    <div class="equation">                }}], {{

                        \\( \\frac{\\partial \\rho_k}{\\partial \\mathbf{x}_i} = \\Omega_k \\sum_j m_j \\nabla_i W_{kj} \\)                    title: {{ text: 'Energy Conservation', font: {{ color: '#fff' }} }},

                    </div>                    paper_bgcolor: 'rgba(0,0,0,0)',

                                        plot_bgcolor: 'rgba(0,0,0,0.2)'

                    <div class="key-insight">                }});

                        <strong>The Ω factor:</strong><br>                return;

                        \\( \\displaystyle \\Omega_i = \\left[1 + \\frac{h_i}{D\\rho_i}\\sum_j m_j \\frac{\\partial W_{ij}}{\\partial h_i}\\right]^{-1} \\)            }}

                    </div>            

                </div>            const E0 = energy.total[0];

                            const relError = energy.total.map(e => ((e - E0) / Math.abs(E0)) * 100);

                <div class="derivation">            

                    <h4>Effect on Gravitational Potential</h4>            Plotly.react('energy-plot', [

                    <p>Hernquist-Katz softening uses: \\( \\epsilon = h/2 \\)</p>                {{ x: energy.time, y: energy.kinetic, name: 'Kinetic', line: {{ color: '#ff6b6b' }} }},

                    <p>So gravitational potential \\( \\phi_i \\) depends on \\( h_i(\\rho_i) \\)!</p>                {{ x: energy.time, y: energy.thermal, name: 'Thermal', line: {{ color: '#ffe66d' }} }},

                                    {{ x: energy.time, y: energy.potential, name: 'Potential', line: {{ color: '#4ecdc4' }} }},

                    <div class="warning">                {{ x: energy.time, y: energy.total, name: 'Total', line: {{ color: '#fff', width: 3 }} }}

                        <strong>Without Ω correction (Ω = 1):</strong><br>            ], {{

                        • Ignores \\( \\partial\\phi/\\partial h \\) contribution<br>                title: {{ text: 'Energy Components', font: {{ color: '#fff' }} }},

                        • Gravitational energy NOT conserved<br>                paper_bgcolor: 'rgba(0,0,0,0)',

                        • Artificial numerical collapse occurs                plot_bgcolor: 'rgba(0,0,0,0.2)',

                    </div>                xaxis: {{ title: 'Time', color: '#aaa', gridcolor: '#333' }},

                                    yaxis: {{ title: 'Energy', color: '#aaa', gridcolor: '#333' }},

                    <div class="key-insight">                legend: {{ font: {{ color: '#fff' }}, x: 0.02, y: 0.98 }},

                        <strong>With Ω correction (Ω ≠ 1):</strong><br>                margin: {{ t: 50, r: 20, b: 50, l: 60 }}

                        • Properly accounts for h(ρ) in potential<br>            }});

                        • Total energy exactly conserved<br>        }}

                        • Cloud remains in hydrostatic equilibrium        

                    </div>        function updateFrame(frameIdx) {{

                                frameIdx = parseInt(frameIdx);

                    <p style="margin-top: 15px; color: #ffd700;">            currentFrame = frameIdx;

                        <strong>Physical insight:</strong> The correction ensures that the pressure-density             

                        relation from the Lagrangian is thermodynamically consistent.            const data = allData[currentDataset];

                    </p>            const frame = data.frames[frameIdx];

                </div>            if (!frame) return;

            </div>            

        </div>            // Update slider

    </div>            document.getElementById('frame-slider').value = frameIdx;

                

    <script>            // Update time display

    // Embedded simulation data            document.getElementById('time-display').textContent = `t = ${{frame.time.toFixed(2)}}`;

    const dataWith = DATA_WITH_PLACEHOLDER;            

    const dataNo = DATA_NO_PLACEHOLDER;            // Update stats

                document.getElementById('stat-time').textContent = frame.time.toFixed(2);

    let currentFrame = 0;            document.getElementById('stat-rho-max').textContent = frame.rho_max.toFixed(3);

    let playing = false;            document.getElementById('stat-rho-center').textContent = frame.rho_center.toFixed(3);

    let playInterval = null;            document.getElementById('stat-frame').textContent = `${{frameIdx + 1}} / ${{data.frames.length}}`;

    const totalFrames = dataWith.frames.length;            

                // Update density plot

    document.getElementById('totalFrames').textContent = totalFrames;            const traces = [{{

    document.getElementById('slider').max = totalFrames - 1;                x: frame.r_profile,

                    y: frame.rho_profile,

    // Color scale for density                type: 'scatter',

    const colorscale = [                mode: 'lines+markers',

        [0, '#000033'],                name: 'SPH',

        [0.2, '#003366'],                line: {{ color: '#00ff88', width: 3 }},

        [0.4, '#006699'],                marker: {{ size: 4 }}

        [0.6, '#ff9900'],            }}];

        [0.8, '#ff3300'],            

        [1, '#ffff00']            // Add Lane-Emden reference

    ];            if (data.lane_emden_r && data.lane_emden_rho) {{

                    traces.push({{

    // Initial ρ for normalization                    x: data.lane_emden_r,

    const rho0_with = dataWith.frames[0].stats.rho_central;                    y: data.lane_emden_rho,

    const rho0_no = dataNo.frames[0].stats.rho_central;                    type: 'scatter',

    const E0_with = Math.abs(dataWith.energy.total[0]);                    mode: 'lines',

    const E0_no = Math.abs(dataNo.energy.total[0]);                    name: 'Lane-Emden',

                        line: {{ color: '#00d4ff', width: 2, dash: 'dash' }}

    function updatePlots() {                }});

        const fw = dataWith.frames[currentFrame];            }}

        const fn = dataNo.frames[currentFrame];            

                    Plotly.react('density-plot', traces, {{

        // Update info                title: {{ text: `Radial Density Profile (t = ${{frame.time.toFixed(2)}})`, font: {{ color: '#fff' }} }},

        document.getElementById('frameNum').textContent = currentFrame + 1;                paper_bgcolor: 'rgba(0,0,0,0)',

        document.getElementById('timeVal').textContent = fw.time.toFixed(1);                plot_bgcolor: 'rgba(0,0,0,0.2)',

                        xaxis: {{ title: 'Radius r', color: '#aaa', gridcolor: '#333', range: [0, frame.r_profile[frame.r_profile.length-1] * 1.1] }},

        // Update statistics                yaxis: {{ title: 'Density ρ', type: 'log', color: '#aaa', gridcolor: '#333' }},

        document.getElementById('rho-central-with').textContent = fw.stats.rho_central.toFixed(3);                legend: {{ font: {{ color: '#fff' }}, x: 0.7, y: 0.95 }},

        document.getElementById('gradh-mean-with').textContent = fw.stats.gradh_mean.toFixed(3);                margin: {{ t: 50, r: 20, b: 50, l: 60 }}

        document.getElementById('gradh-range-with').textContent =             }});

            '[' + fw.stats.gradh_min.toFixed(2) + ', ' + fw.stats.gradh_max.toFixed(2) + ']';            

                    // Add time marker on energy plot

        document.getElementById('rho-central-no').textContent = fn.stats.rho_central.toFixed(3);            const energy = data.energy;

        const rhoIncrease = ((fn.stats.rho_central / rho0_no - 1) * 100).toFixed(1);            if (energy) {{

        document.getElementById('rho-increase').textContent = '+' + rhoIncrease + '%';                Plotly.relayout('energy-plot', {{

                            shapes: [{{

        // Get global density range for consistent coloring                        type: 'line',

        const allRho = [...fw.particles.rho, ...fn.particles.rho];                        x0: frame.time, x1: frame.time,

        const rhoMin = Math.min(...allRho);                        y0: 0, y1: 1,

        const rhoMax = Math.max(...allRho);                        yref: 'paper',

                                line: {{ color: '#ff0', width: 2, dash: 'dot' }}

        // Scatter plot - WITH grad-h                    }}]

        Plotly.react('scatter-with', [{                }});

            x: fw.particles.x,            }}

            y: fw.particles.y,        }}

            mode: 'markers',        

            type: 'scatter',        function togglePlay() {{

            marker: {            isPlaying = !isPlaying;

                size: 4,            const btn = document.getElementById('play-btn');

                color: fw.particles.rho,            

                colorscale: colorscale,            if (isPlaying) {{

                cmin: rhoMin,                btn.textContent = '⏸ Pause';

                cmax: rhoMax,                btn.classList.add('active');

                colorbar: { title: 'ρ', thickness: 15 }                const speed = parseInt(document.getElementById('speed-select').value);

            },                playInterval = setInterval(() => {{

            hovertemplate: 'x: %{x:.3f}<br>y: %{y:.3f}<br>ρ: %{marker.color:.3f}<extra></extra>'                    const maxFrame = allData[currentDataset].frames.length - 1;

        }], {                    if (currentFrame >= maxFrame) {{

            paper_bgcolor: 'rgba(0,0,0,0)',                        currentFrame = 0;

            plot_bgcolor: 'rgba(30,30,50,0.5)',                    }} else {{

            xaxis: { title: 'x', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },                        currentFrame++;

            yaxis: { title: 'y', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)', scaleanchor: 'x' },                    }}

            margin: { t: 10, b: 40, l: 50, r: 80 },                    updateFrame(currentFrame);

            font: { color: '#ccc' }                }}, speed);

        }, { responsive: true });            }} else {{

                        btn.textContent = '▶ Play';

        // Scatter plot - NO grad-h                btn.classList.remove('active');

        Plotly.react('scatter-no', [{                clearInterval(playInterval);

            x: fn.particles.x,            }}

            y: fn.particles.y,        }}

            mode: 'markers',        

            type: 'scatter',        function stepForward() {{

            marker: {            const maxFrame = allData[currentDataset].frames.length - 1;

                size: 4,            if (currentFrame < maxFrame) {{

                color: fn.particles.rho,                updateFrame(currentFrame + 1);

                colorscale: colorscale,            }}

                cmin: rhoMin,        }}

                cmax: rhoMax,        

                colorbar: { title: 'ρ', thickness: 15 }        function stepBackward() {{

            },            if (currentFrame > 0) {{

            hovertemplate: 'x: %{x:.3f}<br>y: %{y:.3f}<br>ρ: %{marker.color:.3f}<extra></extra>'                updateFrame(currentFrame - 1);

        }], {            }}

            paper_bgcolor: 'rgba(0,0,0,0)',        }}

            plot_bgcolor: 'rgba(30,30,50,0.5)',    </script>

            xaxis: { title: 'x', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },</body>

            yaxis: { title: 'y', range: [-1.3, 1.3], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)', scaleanchor: 'x' },</html>

            margin: { t: 10, b: 40, l: 50, r: 80 },'''

            font: { color: '#ccc' }    

        }, { responsive: true });    # Write HTML file

            output_path = Path(output_path)

        // Radial density profile    output_path.parent.mkdir(parents=True, exist_ok=True)

        Plotly.react('radial-profile', [    

            {    with open(output_path, 'w') as f:

                x: fw.radial.r,        f.write(html_content)

                y: fw.radial.rho,    

                mode: 'lines',    print(f"\n✓ Generated interactive HTML viewer: {output_path}")

                name: 'With Ω',    print("  Open in browser to view animation")

                line: { color: '#00ff88', width: 3 }    

            },    return output_path

            {

                x: fn.radial.r,

                y: fn.radial.rho,if __name__ == "__main__":

                mode: 'lines',    import argparse

                name: 'No Ω (Ω=1)',    

                line: { color: '#ff4444', width: 3, dash: 'dash' }    parser = argparse.ArgumentParser(description='Generate interactive HTML animation viewer')

            }    parser.add_argument('--output', '-o', default='simulations/stability/gradh_study_3d/results/animation_viewer.html',

        ], {                        help='Output HTML file path')

            paper_bgcolor: 'rgba(0,0,0,0)',    parser.add_argument('--title', '-t', default='Self-Gravitating Cloud: Grad-H Study',

            plot_bgcolor: 'rgba(30,30,50,0.5)',                        help='Title for the viewer')

            xaxis: { title: 'Radius r', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },    

            yaxis: { title: 'Density ρ', type: 'log', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },    args = parser.parse_args()

            legend: { x: 0.7, y: 0.95, bgcolor: 'rgba(30,30,50,0.8)' },    

            margin: { t: 10, b: 40, l: 60, r: 20 },    # Default results directories

            font: { color: '#ccc' }    results_dirs = {

        }, { responsive: true });        'no_gradh_hk': 'simulations/stability/gradh_study_3d/results/long_test/no_gradh_hk',

                'no_gradh_wendland': 'simulations/stability/gradh_study_3d/results/long_test/no_gradh_wendland',

        // Grad-h distribution (Ω vs radius)    }

        Plotly.react('gradh-distribution', [    

            {    generate_html(results_dirs, args.output, args.title)

                x: fw.radial.r,
                y: fw.radial.gradh,
                mode: 'lines+markers',
                name: 'Ω(r) with correction',
                line: { color: '#00ff88', width: 2 },
                marker: { size: 6 }
            },
            {
                x: fn.radial.r,
                y: fn.radial.r.map(function() { return 1.0; }),
                mode: 'lines',
                name: 'Ω = 1 (forced)',
                line: { color: '#ff4444', width: 2, dash: 'dash' }
            }
        ], {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(30,30,50,0.5)',
            xaxis: { title: 'Radius r', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },
            yaxis: { title: 'Ω (grad-h factor)', range: [0.8, 1.5], color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },
            legend: { x: 0.6, y: 0.95, bgcolor: 'rgba(30,30,50,0.8)' },
            margin: { t: 10, b: 40, l: 60, r: 20 },
            font: { color: '#ccc' },
            shapes: [{
                type: 'line',
                x0: 0, x1: 1.2,
                y0: 1, y1: 1,
                line: { color: '#888', width: 1, dash: 'dot' }
            }],
            annotations: [{
                x: 0.3, y: 1.02,
                text: 'Ω = 1 (no correction)',
                showarrow: false,
                font: { color: '#888', size: 10 }
            }]
        }, { responsive: true });
        
        // Energy plot with current time marker
        var currentTime = fw.time;
        Plotly.react('energy-plot', [
            {
                x: dataWith.energy.time,
                y: dataWith.energy.total.map(function(e) { return e / E0_with; }),
                mode: 'lines',
                name: 'E_total (with Ω)',
                line: { color: '#00ff88', width: 3 }
            },
            {
                x: dataNo.energy.time,
                y: dataNo.energy.total.map(function(e) { return e / E0_no; }),
                mode: 'lines',
                name: 'E_total (Ω=1)',
                line: { color: '#ff4444', width: 3, dash: 'dash' }
            },
            {
                x: dataWith.energy.time,
                y: dataWith.energy.kinetic.map(function(e) { return e / E0_with; }),
                mode: 'lines',
                name: 'E_kin (with)',
                line: { color: '#4ecdc4', width: 1.5 },
                visible: 'legendonly'
            },
            {
                x: dataWith.energy.time,
                y: dataWith.energy.potential.map(function(e) { return e / E0_with; }),
                mode: 'lines',
                name: 'E_grav (with)',
                line: { color: '#ff6b6b', width: 1.5 },
                visible: 'legendonly'
            }
        ], {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(30,30,50,0.5)',
            xaxis: { title: 'Time', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },
            yaxis: { title: 'Energy / |E₀|', color: '#ccc', gridcolor: 'rgba(100,100,150,0.3)' },
            legend: { x: 0.7, y: 0.95, bgcolor: 'rgba(30,30,50,0.8)' },
            margin: { t: 10, b: 40, l: 60, r: 20 },
            font: { color: '#ccc' },
            shapes: [{
                type: 'line',
                x0: currentTime, x1: currentTime,
                y0: -2.5, y1: 0.5,
                line: { color: '#ffd700', width: 2, dash: 'dot' }
            }],
            annotations: [{
                x: currentTime, y: 0.3,
                text: 't=' + currentTime.toFixed(1),
                showarrow: true,
                arrowhead: 2,
                arrowcolor: '#ffd700',
                font: { color: '#ffd700' }
            }]
        }, { responsive: true });
    }
    
    function goToFrame(idx) {
        currentFrame = parseInt(idx);
        document.getElementById('slider').value = currentFrame;
        updatePlots();
    }
    
    function nextFrame() {
        currentFrame = (currentFrame + 1) % totalFrames;
        document.getElementById('slider').value = currentFrame;
        updatePlots();
    }
    
    function prevFrame() {
        currentFrame = (currentFrame - 1 + totalFrames) % totalFrames;
        document.getElementById('slider').value = currentFrame;
        updatePlots();
    }
    
    function togglePlay() {
        playing = !playing;
        var btn = document.getElementById('playBtn');
        if (playing) {
            btn.textContent = '⏸ Pause';
            playInterval = setInterval(nextFrame, 300);
        } else {
            btn.textContent = '▶ Play';
            clearInterval(playInterval);
        }
    }
    
    // Initialize
    updatePlots();
    </script>
</body>
</html>
'''
    
    # Replace placeholders with actual data
    html = html_template.replace('DATA_WITH_PLACEHOLDER', json.dumps(data_with_gradh)).replace('DATA_NO_PLACEHOLDER', json.dumps(data_no_gradh))
    
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"✓ Interactive viewer generated: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate interactive HTML viewer for grad-h study")
    parser.add_argument("--with-gradh", required=True, help="Directory with grad-h results")
    parser.add_argument("--no-gradh", required=True, help="Directory without grad-h results")
    parser.add_argument("--output", required=True, help="Output HTML file path")
    parser.add_argument("--sample", type=int, default=500, help="Number of particles to sample per frame")
    
    args = parser.parse_args()
    
    print("Loading data with grad-h correction...")
    data_with = process_case(args.with_gradh, "With Grad-H (Ω ≠ 1)", args.sample)
    
    print("Loading data without grad-h correction...")
    data_no = process_case(args.no_gradh, "No Grad-H (Ω = 1)", args.sample)
    
    if data_with and data_no:
        os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)
        generate_html(data_with, data_no, args.output)
        print(f"\n✓ Done! Open in browser: {args.output}")
    else:
        print("Error: Could not load data from both directories")
        sys.exit(1)


if __name__ == "__main__":
    main()
