# IMBH-Cloud Tidal Disruption Visualization

Interactive 3D visualization of SPH simulations showing gas cloud tidal disruption by an intermediate-mass black hole (IMBH).

## Features

- 3D particle visualization with Three.js
- Position-Velocity (PV) diagram
- Particle selection and analysis
- Multiple color modes (density, temperature, velocity, Mach number)
- Draggable UI panels
- VTK-style Jet colormap

## GitHub Pages Deployment

### Step 1: Create a GitHub Repository

```bash
# Create a new repo for the visualization (or use existing)
gh repo create imbh-cloud-viz --public --description "IMBH-Cloud Tidal Disruption Visualization"
```

### Step 2: Initialize and Push

```bash
# Navigate to viz folder
cd viz

# Initialize git (if not already)
git init

# Add all files
git add .

# Commit
git commit -m "Initial commit: IMBH-Cloud visualization"

# Add remote
git remote add origin https://github.com/YOUR_USERNAME/imbh-cloud-viz.git

# Push to main branch
git push -u origin main
```

### Step 3: Enable GitHub Pages

1. Go to your repository on GitHub
2. Click **Settings** > **Pages** (in left sidebar)
3. Under "Source", select **Deploy from a branch**
4. Select **main** branch and **/ (root)** folder
5. Click **Save**

Your site will be available at: `https://YOUR_USERNAME.github.io/imbh-cloud-viz/`

## File Structure

```
viz/
├── imbh_cloud_viz.html    # Main HTML file (entry point)
├── datasets.json          # Dataset configuration (for CSV fallback)
├── css/
│   └── style.css          # Styles
├── js/
│   ├── config.js          # Configuration and state
│   ├── visualization.js   # 3D rendering
│   ├── pv-diagram.js      # PV diagram and selection
│   ├── views.js           # Camera views
│   └── main.js            # Entry point and data loading
└── data/
    ├── manifest.json      # Binary data manifest
    ├── parabolic_rp0.13pc/
    │   └── *.bin          # Binary snapshot files
    └── parabolic_rp0.4pc/
        └── *.bin          # Binary snapshot files
```

## Data Format

Binary format (`.bin`) provides ~20x compression over CSV:
- Header: 4 bytes (uint32 particle count)
- Data: 10 × float32 arrays (x, y, z, vx, vy, vz, dens, temp, mass, sound)

Total size: ~39 MB (vs ~740 MB CSV)

## Converting New Data

To add new simulation data:

```bash
# Run conversion script (from parameter_survey directory)
python3 convert_to_binary.py --subsample 2 --output viz/data
```

Options:
- `--subsample N`: Use every Nth snapshot
- `--output DIR`: Output directory

## Local Development

```bash
# Start local server
cd viz
python3 -m http.server 8080

# Open in browser
open http://localhost:8080/imbh_cloud_viz.html
```

## Credits

SPH simulation code and visualization developed for astrophysical research on tidal disruption events.
