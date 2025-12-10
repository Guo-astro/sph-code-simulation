#!/usr/bin/env python3
"""
Generate Standalone Three.js SPH Viewer
========================================

Creates a self-contained HTML file with embedded Three.js visualization
for SPH simulation data. No server required - just open the HTML file.

Usage:
    python generate_threejs_viewer.py <results_dir> [-o output.html]
    python generate_threejs_viewer.py ../results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph

Features:
- 3D particle cloud with orbit controls
- Color by density, velocity, pressure, energy
- Animation playback
- No server required - fully standalone HTML
"""

import numpy as np
import json
import argparse
import base64
import struct
from pathlib import Path
import glob


def load_viz_data(viz_dir: Path):
    """Load all frames from viz_data directory."""
    metadata_path = viz_dir / 'metadata.json'
    if not metadata_path.exists():
        raise FileNotFoundError(f"No metadata.json in {viz_dir}")
    
    with open(metadata_path) as f:
        metadata = json.load(f)
    
    frames = []
    frame_files = sorted(viz_dir.glob('frame_*.bin'))
    
    for bin_path in frame_files:
        json_path = bin_path.with_suffix('.json')
        
        # Read binary data
        with open(bin_path, 'rb') as f:
            data = np.frombuffer(f.read(), dtype=np.float32)
        
        # Read frame metadata
        frame_info = {'time': 0.0}
        if json_path.exists():
            with open(json_path) as f:
                frame_info = json.load(f)
        
        stride = metadata.get('stride', 11)
        n_particles = len(data) // stride
        data = data.reshape(n_particles, stride)
        
        offsets = metadata.get('fieldOffsets', {
            'x': 0, 'y': 1, 'z': 2,
            'vx': 3, 'vy': 4, 'vz': 5,
            'mass': 6, 'density': 7, 'pressure': 8,
            'energy': 9, 'smoothing_length': 10
        })
        
        frame = {
            'time': frame_info.get('time', 0.0),
            'x': data[:, offsets['x']].tolist(),
            'y': data[:, offsets['y']].tolist(),
            'z': data[:, offsets['z']].tolist(),
            'vx': data[:, offsets['vx']].tolist(),
            'vy': data[:, offsets['vy']].tolist(),
            'vz': data[:, offsets['vz']].tolist(),
            'density': data[:, offsets['density']].tolist(),
            'pressure': data[:, offsets['pressure']].tolist(),
            'energy': data[:, offsets['energy']].tolist(),
        }
        frames.append(frame)
    
    return metadata, frames


def generate_html(metadata, frames, output_path: Path):
    """Generate standalone HTML viewer with embedded data."""
    
    # Compute global statistics for color scaling
    all_density = []
    all_velocity = []
    all_pressure = []
    all_energy = []
    
    for f in frames:
        all_density.extend(f['density'])
        vel = np.sqrt(np.array(f['vx'])**2 + np.array(f['vy'])**2 + np.array(f['vz'])**2)
        all_velocity.extend(vel.tolist())
        all_pressure.extend(f['pressure'])
        all_energy.extend(f['energy'])
    
    stats = {
        'density_min': float(np.percentile(all_density, 1)),
        'density_max': float(np.percentile(all_density, 99)),
        'velocity_min': float(np.percentile(all_velocity, 1)),
        'velocity_max': float(np.percentile(all_velocity, 99)),
        'pressure_min': float(np.percentile(all_pressure, 1)),
        'pressure_max': float(np.percentile(all_pressure, 99)),
        'energy_min': float(np.percentile(all_energy, 1)),
        'energy_max': float(np.percentile(all_energy, 99)),
    }
    
    # Subsample if too many particles (for browser performance)
    max_particles = 20000
    n_particles = len(frames[0]['x'])
    stride = max(1, n_particles // max_particles)
    
    if stride > 1:
        print(f"Subsampling from {n_particles} to ~{n_particles//stride} particles for browser performance")
        for f in frames:
            for key in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'pressure', 'energy']:
                f[key] = f[key][::stride]
    
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SPH Visualization - {metadata.get('name', 'Simulation')}</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ 
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    background: #111; 
    color: #eee;
    overflow: hidden;
}}
#container {{ width: 100vw; height: 100vh; }}
#controls {{
    position: fixed;
    bottom: 20px;
    left: 50%;
    transform: translateX(-50%);
    background: rgba(30,30,40,0.95);
    padding: 15px 25px;
    border-radius: 12px;
    display: flex;
    gap: 15px;
    align-items: center;
    box-shadow: 0 4px 20px rgba(0,0,0,0.5);
    z-index: 100;
}}
#settings {{
    position: fixed;
    top: 20px;
    right: 20px;
    background: rgba(30,30,40,0.95);
    padding: 15px;
    border-radius: 12px;
    min-width: 200px;
    box-shadow: 0 4px 20px rgba(0,0,0,0.5);
    z-index: 100;
}}
#info {{
    position: fixed;
    top: 20px;
    left: 20px;
    background: rgba(30,30,40,0.95);
    padding: 15px;
    border-radius: 12px;
    box-shadow: 0 4px 20px rgba(0,0,0,0.5);
    z-index: 100;
}}
button {{
    background: linear-gradient(135deg, #4a90d9, #357abd);
    color: white;
    border: none;
    padding: 8px 16px;
    border-radius: 6px;
    cursor: pointer;
    font-size: 14px;
}}
button:hover {{ background: linear-gradient(135deg, #5aa0e9, #458acd); }}
button:active {{ transform: scale(0.98); }}
input[type="range"] {{
    width: 300px;
    height: 6px;
    -webkit-appearance: none;
    background: #444;
    border-radius: 3px;
}}
input[type="range"]::-webkit-slider-thumb {{
    -webkit-appearance: none;
    width: 16px;
    height: 16px;
    background: #4a90d9;
    border-radius: 50%;
    cursor: pointer;
}}
select {{
    background: #333;
    color: #eee;
    border: 1px solid #555;
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 14px;
}}
label {{ display: block; margin: 8px 0 4px; color: #aaa; font-size: 12px; }}
.stat {{ margin: 5px 0; }}
.stat-label {{ color: #888; font-size: 11px; }}
.stat-value {{ color: #fff; font-weight: bold; }}
#colorbar {{
    position: fixed;
    right: 240px;
    top: 50%;
    transform: translateY(-50%);
    width: 20px;
    height: 200px;
    border-radius: 4px;
    z-index: 100;
}}
#colorbar-labels {{
    position: fixed;
    right: 265px;
    top: 50%;
    transform: translateY(-50%);
    height: 200px;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    font-size: 11px;
    color: #aaa;
    z-index: 100;
}}
</style>
</head>
<body>
<div id="container"></div>

<div id="info">
    <h3 style="margin-bottom:10px;color:#4a90d9">{metadata.get('name', 'SPH Simulation')}</h3>
    <div class="stat">
        <span class="stat-label">Method:</span>
        <span class="stat-value">{metadata.get('method', 'Unknown')}</span>
    </div>
    <div class="stat">
        <span class="stat-label">Particles:</span>
        <span class="stat-value" id="particleCount">{len(frames[0]['x']):,}</span>
    </div>
    <div class="stat">
        <span class="stat-label">Frame:</span>
        <span class="stat-value"><span id="frameNum">1</span> / {len(frames)}</span>
    </div>
    <div class="stat">
        <span class="stat-label">Time:</span>
        <span class="stat-value" id="timeVal">0.000</span>
    </div>
</div>

<div id="settings">
    <label>Color By:</label>
    <select id="colorField" onchange="updateColors()">
        <option value="density">Density</option>
        <option value="velocity">Velocity</option>
        <option value="pressure">Pressure</option>
        <option value="energy">Energy</option>
    </select>
    
    <label>Color Map:</label>
    <select id="colorMap" onchange="updateColors()">
        <option value="viridis">Viridis</option>
        <option value="plasma">Plasma</option>
        <option value="inferno">Inferno</option>
        <option value="turbo">Turbo</option>
    </select>
    
    <label>Point Size:</label>
    <input type="range" id="pointSize" min="0.5" max="5" step="0.1" value="1.5" 
           style="width:100%" onchange="updatePointSize()">
    
    <label>
        <input type="checkbox" id="logScale" onchange="updateColors()"> Log Scale
    </label>
</div>

<div id="colorbar"></div>
<div id="colorbar-labels">
    <span id="colorMax">1.0</span>
    <span id="colorMid">0.5</span>
    <span id="colorMin">0.0</span>
</div>

<div id="controls">
    <button onclick="firstFrame()">⏮</button>
    <button onclick="prevFrame()">◀</button>
    <button onclick="togglePlay()" id="playBtn">▶</button>
    <button onclick="nextFrame()">▶</button>
    <button onclick="lastFrame()">⏭</button>
    <input type="range" id="slider" min="0" max="{len(frames)-1}" value="0" 
           oninput="goToFrame(parseInt(this.value))">
    <select id="speed" style="width:80px">
        <option value="50">20 fps</option>
        <option value="100" selected>10 fps</option>
        <option value="200">5 fps</option>
        <option value="500">2 fps</option>
    </select>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
// Embedded simulation data
const FRAMES = {json.dumps(frames)};
const STATS = {json.dumps(stats)};
const METADATA = {json.dumps(metadata)};

// Color maps
const COLORMAPS = {{
    viridis: [[0.267,0.004,0.329],[0.282,0.140,0.458],[0.253,0.265,0.530],[0.206,0.372,0.553],[0.163,0.471,0.558],[0.128,0.567,0.551],[0.134,0.658,0.517],[0.267,0.749,0.441],[0.478,0.821,0.318],[0.741,0.873,0.150],[0.993,0.906,0.144]],
    plasma: [[0.050,0.030,0.528],[0.294,0.012,0.615],[0.492,0.012,0.658],[0.658,0.134,0.588],[0.798,0.280,0.470],[0.899,0.434,0.329],[0.964,0.596,0.168],[0.988,0.773,0.012],[0.940,0.975,0.131]],
    inferno: [[0.001,0.000,0.014],[0.122,0.047,0.282],[0.341,0.062,0.429],[0.533,0.134,0.416],[0.712,0.213,0.338],[0.863,0.344,0.235],[0.955,0.541,0.145],[0.987,0.756,0.215],[0.988,0.998,0.645]],
    turbo: [[0.190,0.072,0.232],[0.204,0.379,0.831],[0.149,0.658,0.857],[0.304,0.870,0.573],[0.617,0.949,0.318],[0.891,0.876,0.196],[0.987,0.645,0.039],[0.920,0.375,0.009],[0.659,0.122,0.059]]
}};

// Three.js setup
let scene, camera, renderer, points, controls;
let currentFrame = 0;
let isPlaying = false;
let playInterval = null;

function init() {{
    const container = document.getElementById('container');
    
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x111111);
    
    camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 1000);
    
    // Position camera based on bounding box
    const bbox = METADATA.boundingBox || {{min: [-10,-10,-10], max: [10,10,10]}};
    const center = [
        (bbox.min[0] + bbox.max[0]) / 2,
        (bbox.min[1] + bbox.max[1]) / 2,
        (bbox.min[2] + bbox.max[2]) / 2
    ];
    const size = Math.max(
        bbox.max[0] - bbox.min[0],
        bbox.max[1] - bbox.min[1],
        bbox.max[2] - bbox.min[2]
    );
    camera.position.set(center[0] + size, center[1] + size * 0.5, center[2] + size);
    camera.lookAt(center[0], center[1], center[2]);
    
    renderer = new THREE.WebGLRenderer({{ antialias: true }});
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    container.appendChild(renderer.domElement);
    
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.target.set(center[0], center[1], center[2]);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    
    // Create initial points
    createPoints();
    
    // Add axes helper
    const axesHelper = new THREE.AxesHelper(size * 0.3);
    scene.add(axesHelper);
    
    // Add grid
    const gridHelper = new THREE.GridHelper(size, 20, 0x444444, 0x222222);
    gridHelper.position.y = bbox.min[1];
    scene.add(gridHelper);
    
    window.addEventListener('resize', onWindowResize);
    animate();
}}

function createPoints() {{
    const frame = FRAMES[currentFrame];
    const n = frame.x.length;
    
    const geometry = new THREE.BufferGeometry();
    const positions = new Float32Array(n * 3);
    const colors = new Float32Array(n * 3);
    
    for (let i = 0; i < n; i++) {{
        positions[i * 3] = frame.x[i];
        positions[i * 3 + 1] = frame.y[i];
        positions[i * 3 + 2] = frame.z[i];
    }}
    
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    
    const material = new THREE.PointsMaterial({{
        size: parseFloat(document.getElementById('pointSize').value) * 0.1,
        vertexColors: true,
        transparent: true,
        opacity: 0.8,
        sizeAttenuation: true
    }});
    
    if (points) scene.remove(points);
    points = new THREE.Points(geometry, material);
    scene.add(points);
    
    updateColors();
}}

function updatePoints() {{
    const frame = FRAMES[currentFrame];
    const positions = points.geometry.attributes.position.array;
    
    for (let i = 0; i < frame.x.length; i++) {{
        positions[i * 3] = frame.x[i];
        positions[i * 3 + 1] = frame.y[i];
        positions[i * 3 + 2] = frame.z[i];
    }}
    
    points.geometry.attributes.position.needsUpdate = true;
    updateColors();
    
    document.getElementById('frameNum').textContent = currentFrame + 1;
    document.getElementById('timeVal').textContent = frame.time.toFixed(4);
    document.getElementById('slider').value = currentFrame;
}}

function updateColors() {{
    const frame = FRAMES[currentFrame];
    const colors = points.geometry.attributes.color.array;
    const field = document.getElementById('colorField').value;
    const cmapName = document.getElementById('colorMap').value;
    const logScale = document.getElementById('logScale').checked;
    const cmap = COLORMAPS[cmapName];
    
    let values;
    let vmin, vmax;
    
    if (field === 'velocity') {{
        values = frame.vx.map((vx, i) => 
            Math.sqrt(vx*vx + frame.vy[i]*frame.vy[i] + frame.vz[i]*frame.vz[i])
        );
        vmin = STATS.velocity_min;
        vmax = STATS.velocity_max;
    }} else {{
        values = frame[field];
        vmin = STATS[field + '_min'];
        vmax = STATS[field + '_max'];
    }}
    
    if (logScale && vmin > 0) {{
        vmin = Math.log10(vmin);
        vmax = Math.log10(vmax);
    }}
    
    for (let i = 0; i < values.length; i++) {{
        let v = values[i];
        if (logScale && v > 0) v = Math.log10(v);
        
        let t = (v - vmin) / (vmax - vmin);
        t = Math.max(0, Math.min(1, t));
        
        const idx = Math.min(Math.floor(t * (cmap.length - 1)), cmap.length - 2);
        const f = t * (cmap.length - 1) - idx;
        
        colors[i * 3] = cmap[idx][0] * (1 - f) + cmap[idx + 1][0] * f;
        colors[i * 3 + 1] = cmap[idx][1] * (1 - f) + cmap[idx + 1][1] * f;
        colors[i * 3 + 2] = cmap[idx][2] * (1 - f) + cmap[idx + 1][2] * f;
    }}
    
    points.geometry.attributes.color.needsUpdate = true;
    
    // Update colorbar
    const colorbar = document.getElementById('colorbar');
    const gradient = cmap.map((c, i) => 
        `rgb(${{Math.round(c[0]*255)}},${{Math.round(c[1]*255)}},${{Math.round(c[2]*255)}}) ${{100 - i/(cmap.length-1)*100}}%`
    ).join(', ');
    colorbar.style.background = `linear-gradient(${{gradient}})`;
    
    document.getElementById('colorMax').textContent = logScale ? `10^${{vmax.toFixed(1)}}` : vmax.toExponential(2);
    document.getElementById('colorMid').textContent = logScale ? `10^${{((vmin+vmax)/2).toFixed(1)}}` : ((vmin+vmax)/2).toExponential(2);
    document.getElementById('colorMin').textContent = logScale ? `10^${{vmin.toFixed(1)}}` : vmin.toExponential(2);
}}

function updatePointSize() {{
    points.material.size = parseFloat(document.getElementById('pointSize').value) * 0.1;
}}

function animate() {{
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}}

function onWindowResize() {{
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}}

// Playback controls
function goToFrame(n) {{
    currentFrame = Math.max(0, Math.min(FRAMES.length - 1, n));
    updatePoints();
}}

function nextFrame() {{ goToFrame(currentFrame + 1); }}
function prevFrame() {{ goToFrame(currentFrame - 1); }}
function firstFrame() {{ goToFrame(0); }}
function lastFrame() {{ goToFrame(FRAMES.length - 1); }}

function togglePlay() {{
    isPlaying = !isPlaying;
    document.getElementById('playBtn').textContent = isPlaying ? '⏸' : '▶';
    
    if (isPlaying) {{
        const speed = parseInt(document.getElementById('speed').value);
        playInterval = setInterval(() => {{
            if (currentFrame < FRAMES.length - 1) {{
                nextFrame();
            }} else {{
                goToFrame(0);
            }}
        }}, speed);
    }} else {{
        clearInterval(playInterval);
    }}
}}

// Keyboard shortcuts
document.addEventListener('keydown', (e) => {{
    if (e.code === 'Space') {{ togglePlay(); e.preventDefault(); }}
    if (e.code === 'ArrowRight') {{ nextFrame(); }}
    if (e.code === 'ArrowLeft') {{ prevFrame(); }}
    if (e.code === 'Home') {{ firstFrame(); }}
    if (e.code === 'End') {{ lastFrame(); }}
}});

init();
</script>
</body>
</html>'''
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"✓ Generated: {output_path}")
    print(f"  Open in browser to view the simulation")


def main():
    parser = argparse.ArgumentParser(description='Generate standalone Three.js SPH viewer')
    parser.add_argument('results_dir', type=str, help='Directory containing viz_data/')
    parser.add_argument('-o', '--output', type=str, default=None, 
                        help='Output HTML file (default: <results_dir>/viewer.html)')
    parser.add_argument('--max-frames', type=int, default=None,
                        help='Maximum number of frames to include')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    viz_dir = results_dir / 'viz_data'
    
    if not viz_dir.exists():
        # Maybe results_dir is actually the viz_data directory
        if (results_dir / 'metadata.json').exists():
            viz_dir = results_dir
            results_dir = results_dir.parent
        else:
            print(f"Error: No viz_data directory found in {results_dir}")
            print("Run export_viz_data.py first to prepare the data")
            return 1
    
    output_path = Path(args.output) if args.output else results_dir / 'viewer.html'
    
    print(f"Loading data from {viz_dir}...")
    metadata, frames = load_viz_data(viz_dir)
    
    if args.max_frames and len(frames) > args.max_frames:
        stride = len(frames) // args.max_frames
        frames = frames[::stride][:args.max_frames]
        print(f"Subsampled to {len(frames)} frames")
    
    print(f"Generating HTML viewer...")
    print(f"  Frames: {len(frames)}")
    print(f"  Particles per frame: {len(frames[0]['x']):,}")
    
    generate_html(metadata, frames, output_path)
    
    return 0


if __name__ == '__main__':
    exit(main())
