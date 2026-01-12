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
    comMarker: null,
    losArrows: null,
    losSphere: null,
    tidalCircle: null,
    hillCirclePeri: null,

    snapshots: [],
    currentFrame: 0,
    isPlaying: false,
    playInterval: null,

    simType: 'adiabatic',
    colorMode: 'density',

    losVector: new THREE.Vector3(0, 0, 1),
    isDraggingLOS: false,

    // Dynamic data
    cloudMass: 0,
    currentDataset: null,
    datasets: [],
    basePath: '',

    // Global color ranges - use log scale for wide-ranging quantities
    colorRanges: {
        density: { min: 3, max: 6, unit: 'log10(cm^-3)', isLog: true },
        temperature: { min: 1, max: 5, unit: 'log10(K)', isLog: true },
        velocity: { min: 0, max: 1.5, unit: 'log10(km/s)', isLog: true }
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

        // Select first dataset by default
        if (STATE.datasets.length > 0) {
            selectDataset(STATE.datasets[0].id);
        }
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
    const infoEl = document.getElementById('dataset-info');
    if (infoEl) {
        infoEl.innerHTML = `
            <div class="info-row"><span class="label">Pericenter:</span><span class="value">r_p = ${CONFIG.r_peri} pc</span></div>
            <div class="info-row"><span class="label">Description:</span><span class="value">${dataset.config.description || ''}</span></div>
        `;
    }
}

// Precompute global ranges from all snapshots
function computeGlobalColorRanges() {
    let logDensMin = Infinity, logDensMax = -Infinity;
    let logTempMin = Infinity, logTempMax = -Infinity;
    let logVelDispMin = Infinity, logVelDispMax = -Infinity;

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

            // Log velocity dispersion
            const dvx = p.vx - comVx;
            const dvy = p.vy - comVy;
            const dvz = p.vz - comVz;
            const velDisp = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            const logVelDisp = Math.log10(Math.max(velDisp, 0.1));
            logVelDispMin = Math.min(logVelDispMin, logVelDisp);
            logVelDispMax = Math.max(logVelDispMax, logVelDisp);
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
