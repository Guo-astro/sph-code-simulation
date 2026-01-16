// ============================================================
// IMBH-Cloud Visualization - Color Utilities
// Color mapping and colorbar functions
// ============================================================

// Current colormap selection (default: plasma for high contrast on dark bg)
let CURRENT_COLORMAP = 'plasma';

// Colormap definitions - all optimized for dark backgrounds
const COLORMAPS = {
    // Plasma: purple -> pink -> orange -> yellow (high contrast, perceptually uniform)
    plasma: [
        { t: 0.00, r: 0.05, g: 0.03, b: 0.53 },  // Deep purple
        { t: 0.25, r: 0.42, g: 0.05, b: 0.68 },  // Purple
        { t: 0.50, r: 0.80, g: 0.20, b: 0.50 },  // Magenta/Pink
        { t: 0.75, r: 0.97, g: 0.55, b: 0.20 },  // Orange
        { t: 1.00, r: 0.94, g: 0.97, b: 0.13 }   // Yellow
    ],
    // Inferno: dark purple -> red -> orange -> yellow (dramatic)
    inferno: [
        { t: 0.00, r: 0.07, g: 0.04, b: 0.15 },  // Very dark purple
        { t: 0.25, r: 0.34, g: 0.06, b: 0.42 },  // Purple
        { t: 0.50, r: 0.73, g: 0.21, b: 0.33 },  // Red-pink
        { t: 0.75, r: 0.97, g: 0.56, b: 0.12 },  // Orange
        { t: 1.00, r: 0.99, g: 1.00, b: 0.64 }   // Light yellow
    ],
    // Viridis: teal -> green -> yellow (perceptually uniform)
    viridis: [
        { t: 0.00, r: 0.27, g: 0.00, b: 0.33 },  // Dark purple
        { t: 0.25, r: 0.28, g: 0.36, b: 0.55 },  // Blue-purple
        { t: 0.50, r: 0.13, g: 0.56, b: 0.55 },  // Teal
        { t: 0.75, r: 0.48, g: 0.77, b: 0.25 },  // Green
        { t: 1.00, r: 0.99, g: 0.91, b: 0.15 }   // Yellow
    ],
    // Turbo: improved rainbow (better than jet)
    turbo: [
        { t: 0.00, r: 0.19, g: 0.07, b: 0.23 },  // Dark purple
        { t: 0.15, r: 0.14, g: 0.35, b: 0.85 },  // Blue
        { t: 0.35, r: 0.10, g: 0.75, b: 0.75 },  // Cyan
        { t: 0.50, r: 0.45, g: 0.90, b: 0.30 },  // Green
        { t: 0.65, r: 0.90, g: 0.85, b: 0.10 },  // Yellow
        { t: 0.80, r: 0.98, g: 0.50, b: 0.08 },  // Orange
        { t: 1.00, r: 0.90, g: 0.12, b: 0.15 }   // Red
    ],
    // Hot: black -> red -> yellow -> white (thermal)
    hot: [
        { t: 0.00, r: 0.10, g: 0.00, b: 0.00 },  // Very dark red
        { t: 0.35, r: 0.80, g: 0.00, b: 0.00 },  // Red
        { t: 0.65, r: 1.00, g: 0.55, b: 0.00 },  // Orange
        { t: 0.85, r: 1.00, g: 0.90, b: 0.30 },  // Yellow
        { t: 1.00, r: 1.00, g: 1.00, b: 0.80 }   // White-yellow
    ],
    // Jet: classic rainbow (kept for compatibility)
    jet: [
        { t: 0.00, r: 0.0,  g: 0.0,  b: 0.5  },  // Dark blue
        { t: 0.15, r: 0.0,  g: 0.0,  b: 1.0  },  // Blue
        { t: 0.35, r: 0.0,  g: 1.0,  b: 1.0  },  // Cyan
        { t: 0.50, r: 0.0,  g: 1.0,  b: 0.0  },  // Green
        { t: 0.65, r: 1.0,  g: 1.0,  b: 0.0  },  // Yellow
        { t: 0.85, r: 1.0,  g: 0.5,  b: 0.0  },  // Orange
        { t: 1.00, r: 1.0,  g: 0.0,  b: 0.0  }   // Red
    ]
};

// CSS gradients for colorbar display
const COLORMAP_CSS = {
    plasma: 'linear-gradient(to bottom, #F0F921, #F89441, #CC4778, #7E03A8, #0D0887)',
    inferno: 'linear-gradient(to bottom, #FCFFA4, #F7D13D, #BC3754, #57106E, #000004)',
    viridis: 'linear-gradient(to bottom, #FDE725, #7AD151, #22A884, #2A788E, #414487, #440154)',
    turbo: 'linear-gradient(to bottom, #E61F1D, #FB8022, #E6D617, #72CA48, #19BFBF, #2458A4, #30123B)',
    hot: 'linear-gradient(to bottom, #FFFFE6, #FFE64D, #FF8C00, #CC0000, #1A0000)',
    jet: 'linear-gradient(to bottom, #FF0000, #FF7F00, #FFFF00, #00FF00, #00FFFF, #0000FF, #00007F)'
};

function valueToColor(t, colormapName) {
    // Clamp t to [0, 1]
    t = Math.max(0, Math.min(1, t));

    // Use specified colormap or current selection
    const mapName = colormapName || CURRENT_COLORMAP;
    const stops = COLORMAPS[mapName] || COLORMAPS.plasma;

    // Find the two stops to interpolate between
    let i = 0;
    while (i < stops.length - 1 && stops[i + 1].t < t) i++;

    const s0 = stops[i];
    const s1 = stops[Math.min(i + 1, stops.length - 1)];

    // Linear interpolation
    const dt = s1.t - s0.t;
    const f = dt > 0 ? (t - s0.t) / dt : 0;

    const r = s0.r + f * (s1.r - s0.r);
    const g = s0.g + f * (s1.g - s0.g);
    const b = s0.b + f * (s1.b - s0.b);

    return new THREE.Color(r, g, b);
}

function setColormap(name) {
    if (COLORMAPS[name]) {
        CURRENT_COLORMAP = name;
        console.log('Colormap set to:', name);
        return true;
    }
    return false;
}

function updateColorbar() {
    const range = STATE.colorRanges[STATE.colorMode];

    let title, minLabel, maxLabel;
    switch (STATE.colorMode) {
        case 'temperature':
            title = 'T (K)';
            // Show as powers of 10
            minLabel = `10^${range.min}`;
            maxLabel = `10^${range.max}`;
            break;
        case 'velocity':
            title = 'σ_v (km/s)';
            minLabel = `${Math.pow(10, range.min).toFixed(0)}`;
            maxLabel = `${Math.pow(10, range.max).toFixed(0)}`;
            break;
        case 'mach':
            title = 'Mach';
            minLabel = range.min.toFixed(1);
            maxLabel = range.max.toFixed(1);
            break;
        case 'vrel':
            title = '|v-v_COM| (km/s)';
            minLabel = range.min.toFixed(1);
            maxLabel = range.max.toFixed(1);
            break;
        case 'density':
        default:
            title = 'n_H2 (cm⁻³)';
            minLabel = `10^${range.min}`;
            maxLabel = `10^${range.max}`;
            break;
    }

    document.getElementById('colorbar-title').textContent = title;
    document.getElementById('colorbar-min').textContent = minLabel;
    document.getElementById('colorbar-max').textContent = maxLabel;

    // Use current colormap's CSS gradient
    const colorbarEl = document.getElementById('colorbar');
    if (colorbarEl) {
        colorbarEl.style.background = COLORMAP_CSS[CURRENT_COLORMAP] || COLORMAP_CSS.plasma;
    }
}
