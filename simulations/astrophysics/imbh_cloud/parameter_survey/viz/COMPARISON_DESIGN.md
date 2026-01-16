# Cooling vs Adiabatic Comparison Design

## Overview
Side-by-side frame-by-frame comparison of cooling and adiabatic simulations with synchronized time.

## Data Summary
| Property | Cooling | Adiabatic |
|----------|---------|-----------|
| Snapshots | 204 | 70 |
| Time range | 0 → 6.01 | 0 → 1.38 |
| Output interval | ~0.029 | ~0.020 |
| Physics | Koyama-Inutsuka cooling | No cooling |

**Comparison window**: 0 to 1.38 code time (~1.35 Myr)

## Architecture Design

### 1. UI Layout: Split-Screen View
```
┌────────────────────────────────────────────────────────────────┐
│  [Dataset: Comparison Mode ▼]  [Color: density ▼]  [Exit]      │
├─────────────────────────────┬──────────────────────────────────┤
│                             │                                  │
│       COOLING               │         ADIABATIC                │
│       (3D View)             │         (3D View)                │
│                             │                                  │
│   N=28,542  t=0.50 Myr      │    N=31,063  t=0.50 Myr          │
│                             │                                  │
├─────────────────────────────┴──────────────────────────────────┤
│ ◀ ▶ ▶▶  ═══════════●════════════════  t = 0.50 / 1.35 Myr     │
├─────────────────────────────┬──────────────────────────────────┤
│     PV Diagram (Cooling)    │     PV Diagram (Adiabatic)       │
├─────────────────────────────┼──────────────────────────────────┤
│     Profile Panel           │     Comparison Metrics           │
│     (Column/Plane)          │     ΔN, ΔM, ΔE, etc.             │
└─────────────────────────────┴──────────────────────────────────┘
```

### 2. Time Synchronization Strategy
- **Master timeline** based on simulation time (not frame index)
- **Slider range**: 0 to min(t_cooling_max, t_adiabatic_max) = 1.38 code time
- **Frame lookup**: Find nearest snapshot by time for each simulation
- **Interpolation**: Optional linear interpolation for smoother playback

```javascript
// Time-based frame lookup
function getFrameAtTime(snapshots, targetTime) {
    let bestIdx = 0;
    let bestDiff = Infinity;
    for (let i = 0; i < snapshots.length; i++) {
        const diff = Math.abs(snapshots[i].time - targetTime);
        if (diff < bestDiff) {
            bestDiff = diff;
            bestIdx = i;
        }
    }
    return bestIdx;
}
```

### 3. Camera Synchronization
- **Shared OrbitControls**: Both views share same camera orientation
- **Implementation**: Single camera state, apply to both renderers

```javascript
// Sync cameras
function syncCameras() {
    STATE.cameraRight.position.copy(STATE.cameraLeft.position);
    STATE.cameraRight.quaternion.copy(STATE.cameraLeft.quaternion);
    STATE.controlsRight.target.copy(STATE.controlsLeft.target);
}
```

### 4. Data Loading Strategy
Load both datasets in parallel, build time index:

```javascript
STATE.comparison = {
    enabled: false,
    cooling: {
        snapshots: [],
        timeIndex: []  // [{ time: 0.0, frameIdx: 0 }, ...]
    },
    adiabatic: {
        snapshots: [],
        timeIndex: []
    },
    maxTime: 0,  // min of both max times
    currentTime: 0
};
```

### 5. Comparison Metrics Panel
Show difference metrics at each time step:
- **ΔN**: Particle count difference (accretion effect)
- **ΔM_total**: Total mass difference
- **ΔT_mean**: Mean temperature difference
- **ΔE_thermal**: Thermal energy difference
- **COM distance**: How far have the cloud centers drifted apart

### 6. datasets.json Extension
```json
{
    "datasets": [
        {
            "id": "cooling_results",
            "name": "Parabolic r_p = 0.4 pc + Cooling",
            "path": "data/cooling_results",
            "forceCSV": true,
            "config": { ... }
        },
        {
            "id": "adiabatic_results",
            "name": "Parabolic r_p = 0.4 pc (Adiabatic)",
            "path": "data/adiabatic_results",
            "forceCSV": true,
            "comparisonPair": "cooling_results",
            "config": {
                "r_peri": 0.4,
                "cloud_pos0": [-14.4675, 4.8773, 0.0],
                "cloud_vel0": [7.4071, -1.2149, 0.0],
                "description": "Adiabatic reference run (no cooling)"
            }
        }
    ],
    "comparisons": [
        {
            "id": "cooling_vs_adiabatic",
            "name": "Cooling vs Adiabatic Comparison",
            "left": "cooling_results",
            "right": "adiabatic_results",
            "description": "Effect of radiative cooling on tidal disruption"
        }
    ]
}
```

## Implementation Plan

### Phase 1: Data Preparation
1. Create symlink or copy adiabatic data to viz/data/adiabatic_results/
2. Update datasets.json with adiabatic dataset
3. Parse time from CSV headers during loading

### Phase 2: Core Comparison Infrastructure
1. Add comparison state to config.js
2. Modify data-loader.js to load both datasets
3. Build time-indexed frame lookup

### Phase 3: Split-Screen UI
1. Modify HTML for dual canvas layout
2. Create second Three.js renderer
3. Implement camera synchronization
4. Add comparison timeline slider

### Phase 4: Synchronized Rendering
1. Time-based frame selection
2. Synchronized particle updates
3. Parallel PV diagrams

### Phase 5: Metrics Panel
1. Compute comparison metrics
2. Display difference statistics
3. Time evolution plots (ΔN vs t, etc.)

## Files to Modify/Create

### New Files:
- `js/comparison.js` - Comparison mode logic
- `css/comparison.css` - Split-screen styles

### Modified Files:
- `imbh_cloud_viz.html` - Add dual canvas, comparison UI
- `js/config.js` - Add comparison state
- `js/data-loader.js` - Time parsing, dual loading
- `js/main.js` - Comparison mode toggle
- `datasets.json` - Add adiabatic dataset

## Physical Insights to Highlight

The comparison will reveal:
1. **Cooling-induced fragmentation**: Cooling run shows clumping
2. **Enhanced accretion**: Cooling removes pressure support → more accretion
3. **Temperature evolution**: Cooling maintains low T, adiabatic heats up
4. **Morphology difference**: Stream structure differs significantly
5. **Tidal radius effect**: Cooling changes effective cloud compressibility
