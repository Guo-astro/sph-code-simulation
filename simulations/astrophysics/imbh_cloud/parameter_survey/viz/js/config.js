// ============================================================
// IMBH-Cloud Tidal Disruption - Configuration
// Multi-dataset support with dynamic loading
// ============================================================

// Physical constants (code units: L=pc, M=M_sun, V=km/s)
const CONFIG = {
    G: 0.00430091,
    M_BH: 100000.0,

    // Will be set dynamically from selected dataset
    r_peri: 0.13,
    cloud_pos0: [-15.0, 2.8049, 0.0],
    cloud_vel0: [7.2211, -2.0553, 0.0],

    // Derived orbital parameters (computed after dataset load)
    L_vec: null,
    L: null,
    p: null,

    // Temperature conversion: T = (gamma-1) * mu * m_H * u / k_B
    tempConversion: 186.0,

    // Density conversion: M_sun/pc^3 to n_H2 (cm^-3)
    densityToNH2: 20.3,

    units: {
        density: 'cm^-3',
        temperature: 'K',
        velocity: 'km/s',
        position: 'pc',
        time: 'Myr'
    },

    timeToMyr: 0.978,
    bulkVelocity: [7.2211, -2.0553, 0.0]
};

// State variables
const STATE = {
    scene: null,
    camera: null,
    renderer: null,
    controls: null,

    particleSystem: null,
    bhMesh: null,
    orbitLine: null,
    losArrows: null,
    losSphere: null,
    periMarker: null,

    snapshots: [],
    currentFrame: 0,
    isPlaying: false,
    playInterval: null,

    simType: 'adiabatic',
    colorMode: 'density',

    // Orbital pole direction (controlled by draggable cyan sphere)
    // This represents the orbital angular momentum direction
    losVector: new THREE.Vector3(0, 0, 1),  // Current orbital pole (user-controlled)

    // Original orbital pole from simulation (computed from CONFIG.L_vec)
    // This is the "default" orientation - rotation is computed relative to this
    originalOrbitalPole: new THREE.Vector3(0, 0, 1),  // Will be set from L_vec

    // Group containing all orbital elements (orbit line, particles)
    // This group is rotated when the orbital pole changes
    orbitalPlaneGroup: null,

    // Fixed LOS direction toward Sun (from HVCC)
    // This is the actual observation direction from Earth - FIXED, not user-controllable
    // Direction computed from HVCC's galactic position: (l=-0.40°, b=-0.22°)
    fixedLOStoSun: new THREE.Vector3(0, 0, 1),  // Will be updated when galaxy context loads

    isDraggingLOS: false,
    pvZoom: 1.0,  // PV diagram zoom level

    // PV selection state
    pvSelection: null,       // { x1, y1, x2, y2 } in canvas coords
    pvSelecting: false,      // Currently dragging selection
    selectedParticles: [],   // Array of selected particle indices
    pvRanges: null,          // { minPos, maxPos, minVel, maxVel } for coord conversion

    // Dynamic data
    cloudMass: 0,
    cloudRadius: 0,
    initialCloudParticles: 0,  // Initial cloud particle count (for accretion tracking)
    currentDataset: null,
    datasets: [],
    basePath: '',
    isLoading: false,  // Guard against concurrent loads
    loadId: 0,  // Unique ID for each load operation

    // Global color ranges - use log scale for wide-ranging quantities
    colorRanges: {
        density: { min: 3, max: 6, unit: 'log10(cm^-3)', isLog: true },
        temperature: { min: 1, max: 5, unit: 'log10(K)', isLog: true },
        velocity: { min: 0, max: 1.5, unit: 'log10(km/s)', isLog: true },
        vrel: { min: 0, max: 5, unit: 'km/s', isLog: false },  // Linear |v - v_COM|
        mach: { min: 0, max: 3, unit: 'Mach', isLog: false }   // Linear scale for Mach
    }
};

// Load datasets configuration
async function loadDatasets() {
    try {
        const response = await fetch('datasets.json');
        const data = await response.json();
        STATE.datasets = data.datasets;

        // Apply common config
        if (data.common) {
            Object.assign(CONFIG, data.common);
        }

        // Populate dataset selector
        const selector = document.getElementById('dataset-select');
        if (selector) {
            selector.innerHTML = '';
            STATE.datasets.forEach(ds => {
                const option = document.createElement('option');
                option.value = ds.id;
                option.textContent = ds.name;
                selector.appendChild(option);
            });
        }

        // Check URL parameter first, then global variable, then use first dataset (cooling)
        const urlParams = new URLSearchParams(window.location.search);
        const urlDataset = urlParams.get('dataset');
        let datasetId;
        if (urlDataset && STATE.datasets.find(d => d.id === urlDataset)) {
            datasetId = urlDataset;
        } else if (typeof SELECTED_DATASET !== 'undefined') {
            datasetId = SELECTED_DATASET;
        } else {
            datasetId = STATE.datasets[0].id;  // First dataset is cooling
        }
        console.log('Selecting dataset:', datasetId);
        selectDataset(datasetId);

        // Update selector to show selected dataset
        if (selector) {
            selector.value = datasetId;
        }

        console.log('Loaded dataset:', datasetId);
    } catch (error) {
        console.error('Failed to load datasets:', error);
    }
}

// Select and configure a dataset
function selectDataset(datasetId) {
    const dataset = STATE.datasets.find(ds => ds.id === datasetId);
    if (!dataset) {
        console.error('Dataset not found:', datasetId);
        return;
    }

    STATE.currentDataset = dataset;
    STATE.basePath = dataset.path;

    // Update CONFIG with dataset-specific values
    if (dataset.config) {
        CONFIG.r_peri = dataset.config.r_peri;
        CONFIG.cloud_pos0 = dataset.config.cloud_pos0;
        CONFIG.cloud_vel0 = dataset.config.cloud_vel0;
        CONFIG.bulkVelocity = dataset.config.cloud_vel0;
    }

    // Set simType based on dataset (cooling datasets default to 'cooling')
    if (datasetId.includes('cooling')) {
        STATE.simType = 'cooling';
        const simTypeSelect = document.getElementById('sim-type');
        if (simTypeSelect) simTypeSelect.value = 'cooling';
    }

    // Recompute orbital parameters
    computeOrbitalParams();

    // Update dataset info display
    updateDatasetInfo(dataset);

    console.log('Selected dataset:', datasetId);
}

// Compute derived orbital parameters
function computeOrbitalParams() {
    const pos = CONFIG.cloud_pos0;
    const vel = CONFIG.cloud_vel0;

    CONFIG.L_vec = [
        pos[1] * vel[2] - pos[2] * vel[1],
        pos[2] * vel[0] - pos[0] * vel[2],
        pos[0] * vel[1] - pos[1] * vel[0]
    ];

    CONFIG.L = Math.sqrt(
        CONFIG.L_vec[0]**2 +
        CONFIG.L_vec[1]**2 +
        CONFIG.L_vec[2]**2
    );

    CONFIG.p = CONFIG.L * CONFIG.L / (CONFIG.G * CONFIG.M_BH);
}

// Update dataset info in UI
function updateDatasetInfo(dataset) {
    const periEl = document.getElementById('pericenter-display');
    if (periEl) {
        periEl.textContent = `${CONFIG.r_peri} pc`;
    }
    console.log(`Dataset: ${dataset.name}, r_peri = ${CONFIG.r_peri} pc`);
}

// Update penetration factor display (called after data loads)
function updatePenetrationFactor() {
    if (STATE.cloudMass <= 0 || STATE.cloudRadius <= 0) return;

    // Calculate penetration factor β = r_t / r_p
    // r_t = R_cloud * (M_BH / M_cloud)^(1/3)
    const r_tidal = STATE.cloudRadius * Math.pow(CONFIG.M_BH / STATE.cloudMass, 1/3);
    const beta = r_tidal / CONFIG.r_peri;

    const betaEl = document.getElementById('penetration-display');
    if (betaEl) {
        betaEl.textContent = `${beta.toFixed(2)}`;
    }

    console.log(`Penetration factor: β = ${beta.toFixed(2)} (r_t = ${r_tidal.toFixed(3)} pc, r_p = ${CONFIG.r_peri} pc)`);
}

// Precompute global ranges from all snapshots
function computeGlobalColorRanges() {
    let logDensMin = Infinity, logDensMax = -Infinity;
    let logTempMin = Infinity, logTempMax = -Infinity;
    let logVelDispMin = Infinity, logVelDispMax = -Infinity;
    let vrelMax = 0;
    let machMax = 0;

    for (const snapshot of STATE.snapshots) {
        let comVx = 0, comVy = 0, comVz = 0, totalMass = 0;
        for (const p of snapshot.particles) {
            if (p.is_ghost === 1) continue;
            comVx += p.vx * p.mass;
            comVy += p.vy * p.mass;
            comVz += p.vz * p.mass;
            totalMass += p.mass;
        }
        comVx /= totalMass;
        comVy /= totalMass;
        comVz /= totalMass;

        for (const p of snapshot.particles) {
            if (p.is_ghost === 1) continue;

            // Log density
            const n_H2 = p.dens * CONFIG.densityToNH2;
            const logDens = Math.log10(Math.max(n_H2, 1));
            logDensMin = Math.min(logDensMin, logDens);
            logDensMax = Math.max(logDensMax, logDens);

            // Log temperature
            const logTemp = Math.log10(Math.max(p.temp, 1));
            logTempMin = Math.min(logTempMin, logTemp);
            logTempMax = Math.max(logTempMax, logTemp);

            // Velocity dispersion (relative to COM)
            const dvx = p.vx - comVx;
            const dvy = p.vy - comVy;
            const dvz = p.vz - comVz;
            const velDisp = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

            // Log velocity for log mode
            const logVelDisp = Math.log10(Math.max(velDisp, 0.1));
            logVelDispMin = Math.min(logVelDispMin, logVelDisp);
            logVelDispMax = Math.max(logVelDispMax, logVelDisp);

            // Linear |v - v_COM|
            vrelMax = Math.max(vrelMax, velDisp);

            // Mach number
            if (p.sound > 0) {
                const mach = velDisp / p.sound;
                machMax = Math.max(machMax, mach);
            }
        }
    }

    // Use sensible astrophysical ranges (clamp to reasonable values)
    STATE.colorRanges.density = {
        min: Math.max(2, Math.floor(logDensMin)),
        max: Math.min(7, Math.ceil(logDensMax)),
        unit: 'log10(cm^-3)',
        isLog: true
    };

    STATE.colorRanges.temperature = {
        min: Math.max(0.5, Math.floor(logTempMin * 2) / 2),
        max: Math.min(6, Math.ceil(logTempMax)),
        unit: 'log10(K)',
        isLog: true
    };

    STATE.colorRanges.velocity = {
        min: Math.max(-1, Math.floor(logVelDispMin * 2) / 2),
        max: Math.min(2, Math.ceil(logVelDispMax)),
        unit: 'log10(km/s)',
        isLog: true
    };

    // Linear |v - v_COM| in km/s
    STATE.colorRanges.vrel = {
        min: 0,
        max: Math.max(1, Math.ceil(vrelMax)),
        unit: 'km/s',
        isLog: false
    };

    // Mach number: linear scale, 0 to rounded max
    STATE.colorRanges.mach = {
        min: 0,
        max: Math.max(3, Math.ceil(machMax)),
        unit: 'Mach',
        isLog: false
    };

    console.log('Global color ranges computed:', STATE.colorRanges);
}

// Get normalized color value [0,1] for a particle
function getColorForValue(particle, comVelocity) {
    const range = STATE.colorRanges[STATE.colorMode];
    let value;

    switch (STATE.colorMode) {
        case 'temperature':
            value = Math.log10(Math.max(particle.temp, 1));
            break;
        case 'velocity':
            let velDisp;
            if (comVelocity) {
                const dvx = particle.vx - comVelocity.x;
                const dvy = particle.vy - comVelocity.y;
                const dvz = particle.vz - comVelocity.z;
                velDisp = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            } else {
                velDisp = particle.vel_mag;
            }
            value = Math.log10(Math.max(velDisp, 0.1));
            break;
        case 'vrel':
            // Linear |v - v_COM| in km/s
            if (comVelocity) {
                const dvx = particle.vx - comVelocity.x;
                const dvy = particle.vy - comVelocity.y;
                const dvz = particle.vz - comVelocity.z;
                value = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            } else {
                value = particle.vel_mag;
            }
            break;
        case 'mach':
            // Mach number: relative velocity / sound speed (linear scale)
            if (comVelocity && particle.sound > 0) {
                const dvx = particle.vx - comVelocity.x;
                const dvy = particle.vy - comVelocity.y;
                const dvz = particle.vz - comVelocity.z;
                const vRel = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
                value = vRel / particle.sound;
            } else {
                value = 0;
            }
            break;
        case 'density':
        default:
            const n_H2 = particle.dens * CONFIG.densityToNH2;
            value = Math.log10(Math.max(n_H2, 1));
            break;
    }

    const t = Math.max(0, Math.min(1, (value - range.min) / (range.max - range.min)));
    return t;
}

// Initialize orbital parameters with defaults
computeOrbitalParams();
