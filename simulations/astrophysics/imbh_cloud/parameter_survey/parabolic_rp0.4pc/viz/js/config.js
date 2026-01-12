// ============================================================
// IMBH-Cloud Tidal Disruption - Configuration
// Physical constants and orbital parameters
// ============================================================

// Physical constants (code units: L=pc, M=M_sun, V=km/s)
const CONFIG = {
    G: 0.00430091,      // Gravitational constant in code units
    M_BH: 100000.0,     // Black hole mass [M_sun]
    r_peri: 0.13,       // Pericenter distance [pc]

    // Initial conditions from config
    cloud_pos0: [-15.0, 2.8049, 0.0],
    cloud_vel0: [7.2211, -2.0553, 0.0],

    // Derived orbital parameters
    L_vec: null,
    L: null,
    p: null,

    // Temperature conversion: T = (gamma-1) * mu * m_H * u / k_B
    // gamma = 5/3, mu = 2.3 (molecular), u in (km/s)^2
    // T = 0.67 * 2.3 * 1.67e-24 * u * 1e10 / 1.38e-16 = 186 * u
    tempConversion: 186.0,  // K per (km/s)^2

    // Density conversion: M_sun/pc^3 to n_H2 (cm^-3)
    // n_H2 = rho / (2 * m_H) = rho [g/cm^3] / 3.34e-24
    // 1 M_sun/pc^3 = 6.77e-23 g/cm^3
    // n_H2 = 6.77e-23 / 3.34e-24 = 20.3 cm^-3 per M_sun/pc^3
    densityToNH2: 20.3,

    // Physical units
    units: {
        density: 'cm^-3',
        temperature: 'K',
        velocity: 'km/s',
        position: 'pc',
        time: 'Myr'
    },

    // Time conversion: 1 code unit = 0.978 Myr
    timeToMyr: 0.978,

    // Bulk velocity of cloud (for velocity dispersion calculation)
    bulkVelocity: [7.2211, -2.0553, 0.0]
};

// Dynamic cloud properties computed from data
STATE.cloudMass = 0;  // Will be computed from first snapshot

// Compute derived orbital parameters
(function computeOrbitalParams() {
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
})();

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

    snapshots: [],
    currentFrame: 0,
    isPlaying: false,
    playInterval: null,

    simType: 'adiabatic',
    colorMode: 'density',

    losVector: new THREE.Vector3(0, 0, 1),
    isDraggingLOS: false,

    // Global color ranges - use log scale for wide-ranging quantities
    // These are sensible astrophysical ranges, actual data may exceed but will be clipped
    colorRanges: {
        density: { min: 3, max: 6, unit: 'log10(cm^-3)', isLog: true },      // 1e3 - 1e6 cm^-3
        temperature: { min: 1, max: 5, unit: 'log10(K)', isLog: true },       // 10 - 100,000 K
        velocity: { min: 0, max: 1.5, unit: 'log10(km/s)', isLog: true }      // 1 - 30 km/s dispersion
    }
};

// Precompute global ranges from all snapshots
function computeGlobalColorRanges() {
    let logDensMin = Infinity, logDensMax = -Infinity;
    let tempMin = Infinity, tempMax = -Infinity;
    let velDispMin = Infinity, velDispMax = -Infinity;

    // Get bulk velocity at each frame (COM velocity)
    for (const snapshot of STATE.snapshots) {
        // Compute COM velocity for this snapshot
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

            // Log density in n_H2 (cm^-3)
            const n_H2 = p.dens * CONFIG.densityToNH2;
            const logDens = Math.log10(Math.max(n_H2, 1));
            logDensMin = Math.min(logDensMin, logDens);
            logDensMax = Math.max(logDensMax, logDens);

            // Temperature (K)
            tempMin = Math.min(tempMin, p.temp);
            tempMax = Math.max(tempMax, p.temp);

            // Velocity dispersion (relative to COM)
            const dvx = p.vx - comVx;
            const dvy = p.vy - comVy;
            const dvz = p.vz - comVz;
            const velDisp = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            velDispMin = Math.min(velDispMin, velDisp);
            velDispMax = Math.max(velDispMax, velDisp);
        }
    }

    // Set ranges with nice round numbers
    STATE.colorRanges.density = {
        min: Math.floor(logDensMin),
        max: Math.ceil(logDensMax),
        unit: 'log10(cm^-3)',
        isLog: true
    };

    STATE.colorRanges.temperature = {
        min: 0,
        max: Math.ceil(tempMax / 10) * 10,
        unit: 'K',
        isLog: false
    };

    STATE.colorRanges.velocity = {
        min: 0,
        max: Math.ceil(velDispMax),
        unit: 'km/s (dispersion)',
        isLog: false
    };

    console.log('Global color ranges computed:', STATE.colorRanges);
}

// Get normalized color value [0,1] for a particle
// All quantities use log scale due to large dynamic range in tidal disruption
function getColorForValue(particle, comVelocity) {
    const range = STATE.colorRanges[STATE.colorMode];
    let value;

    switch (STATE.colorMode) {
        case 'temperature':
            // Log temperature (K)
            value = Math.log10(Math.max(particle.temp, 1));
            break;
        case 'velocity':
            // Log velocity dispersion relative to COM (km/s)
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
            // Log density in n_H2 (cm^-3)
            const n_H2 = particle.dens * CONFIG.densityToNH2;
            value = Math.log10(Math.max(n_H2, 1));
            break;
    }

    // Normalize to [0,1], clipped to range
    const t = Math.max(0, Math.min(1, (value - range.min) / (range.max - range.min)));
    return t;
}
