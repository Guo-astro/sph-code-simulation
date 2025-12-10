# SPH Visualization Tool

Interactive web-based visualization dashboard for Smoothed Particle Hydrodynamics (SPH) simulation data.

## Features

- **3D Particle Visualization**: Real-time WebGL rendering using React Three Fiber
- **2D Projections**: XY, XZ, and YZ plane projections
- **Time Evolution Charts**: Energy, momentum, and conservation tracking
- **Animation Playback**: Frame-by-frame or continuous playback with adjustable speed
- **Multiple Color Maps**: Viridis, Plasma, Inferno, and custom color maps
- **Color by Field**: Density, pressure, velocity, energy, Mach number
- **Logarithmic Scaling**: Support for large dynamic range data

## Getting Started

### Prerequisites

- Node.js 18+
- Python 3.8+ (for data export)
- NumPy

### Installation

\`\`\`bash
cd tools/sph-viz
npm install
\`\`\`

### Preparing Simulation Data

Before visualizing, you need to export your simulation data to the visualization format:

\`\`\`bash
# Export a single simulation
python scripts/export_viz_data.py /path/to/simulation/results

# Export with options
python scripts/export_viz_data.py sample/sedov/results/gsph_wendland \\
    --stride 2 \\      # Export every 2nd frame
    --max-frames 100  # Limit to 100 frames
\`\`\`

The exporter will create a \`viz_data/\` directory containing:
- \`metadata.json\`: Simulation metadata
- \`frame_XXXXX.bin\`: Binary particle data
- \`frame_XXXXX.json\`: Per-frame metadata

### Running the Dashboard

\`\`\`bash
# Development mode
npm run dev

# Production build
npm run build
npm run start
\`\`\`

Open http://localhost:3000 in your browser.

## Project Structure

\`\`\`
sph-viz/
├── src/
│   ├── components/
│   │   ├── viewer/
│   │   │   ├── ParticleViewer3D.tsx   # WebGL 3D particle cloud
│   │   │   └── Projection2D.tsx       # 2D canvas projections
│   │   ├── charts/
│   │   │   └── Charts.tsx             # Recharts energy/momentum plots
│   │   ├── controls/
│   │   │   ├── PlaybackControls.tsx   # Animation playback
│   │   │   └── VisualizationSettings.tsx
│   │   └── dashboard/
│   │       └── Dashboard.tsx          # Main dashboard layout
│   ├── routes/
│   │   ├── __root.tsx                 # Root layout
│   │   ├── index.tsx                  # Home page
│   │   ├── viz/
│   │   │   └── index.tsx              # Visualization page
│   │   └── api/
│   │       ├── simulations.ts         # List simulations
│   │       ├── simulations.\$simId.ts  # Get simulation metadata
│   │       └── simulations.\$simId.frames.\$frameId.ts  # Get frame data
│   └── types/
│       └── sph.ts                     # TypeScript type definitions
├── scripts/
│   └── export_viz_data.py             # Python data exporter
└── package.json
\`\`\`

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| Space | Play/Pause |
| ← | Previous frame |
| → | Next frame |
| Home | Go to first frame |
| End | Go to last frame |

## Supported SPH Methods

- GSPH (Godunov SPH)
- SSPH (Standard SPH)
- DISPH (Density-Independent SPH)
- GDISPH (Gradient-Corrected DISPH)
- SRGSPH (Special Relativistic GSPH)

## Technology Stack

- [TanStack Start](https://tanstack.com/start) - React full-stack framework
- [React Three Fiber](https://docs.pmnd.rs/react-three-fiber) - React renderer for Three.js
- [drei](https://github.com/pmndrs/drei) - React Three Fiber helpers
- [Recharts](https://recharts.org/) - React charting library
- [Tailwind CSS](https://tailwindcss.com/) - Utility-first CSS
- [TypeScript](https://www.typescriptlang.org/) - Type-safe JavaScript

## License

MIT
