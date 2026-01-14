// ============================================================
// IMBH-Cloud Visualization - Color Utilities
// Color mapping and colorbar functions
// ============================================================

// VTK-style Jet colormap: dark blue -> blue -> cyan -> green -> yellow -> orange -> red
function valueToColor(t) {
    // Clamp t to [0, 1]
    t = Math.max(0, Math.min(1, t));

    // Jet colormap control points
    const stops = [
        { t: 0.00, r: 0.0,  g: 0.0,  b: 0.5  },  // Dark blue
        { t: 0.15, r: 0.0,  g: 0.0,  b: 1.0  },  // Blue
        { t: 0.35, r: 0.0,  g: 1.0,  b: 1.0  },  // Cyan
        { t: 0.50, r: 0.0,  g: 1.0,  b: 0.0  },  // Green
        { t: 0.65, r: 1.0,  g: 1.0,  b: 0.0  },  // Yellow
        { t: 0.85, r: 1.0,  g: 0.5,  b: 0.0  },  // Orange
        { t: 1.00, r: 1.0,  g: 0.0,  b: 0.0  }   // Red
    ];

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

    // VTK-style Jet colorbar matching valueToColor function
    // Gradient from top (max=1.0) to bottom (min=0.0)
    document.getElementById('colorbar').style.background =
        'linear-gradient(to bottom, ' +
        '#FF0000 0%, ' +    // Red (t=1.00)
        '#FF7F00 15%, ' +   // Orange (t=0.85)
        '#FFFF00 35%, ' +   // Yellow (t=0.65)
        '#00FF00 50%, ' +   // Green (t=0.50)
        '#00FFFF 65%, ' +   // Cyan (t=0.35)
        '#0000FF 85%, ' +   // Blue (t=0.15)
        '#00007F 100%)';
}
