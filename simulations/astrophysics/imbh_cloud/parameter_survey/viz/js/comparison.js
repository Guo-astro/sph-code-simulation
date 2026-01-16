// ============================================================
// IMBH-Cloud Visualization - Comparison Mode
// Side-by-side comparison of cooling vs adiabatic simulations
// ============================================================

const COMPARISON = {
    enabled: false,

    // Left panel (primary - cooling)
    left: {
        datasetId: null,
        basePath: null,
        snapshotFiles: [],      // List of filenames
        snapshotCache: {},      // { frameIdx: snapshotData } - sparse cache
        timeIndex: [],          // [{ time, frameIdx, filename }, ...]
        currentFrame: 0,
        scene: null,
        camera: null,
        renderer: null,
        controls: null,
        particleSystem: null
    },

    // Right panel (secondary - adiabatic)
    right: {
        datasetId: null,
        basePath: null,
        snapshotFiles: [],
        snapshotCache: {},
        timeIndex: [],
        currentFrame: 0,
        scene: null,
        camera: null,
        renderer: null,
        controls: null,
        particleSystem: null
    },

    // Time synchronization
    currentTime: 0,
    maxTime: 0,
    minTime: 0,

    // Cache settings (keep only a few frames in memory)
    maxCacheSize: 10,  // Max frames to keep per side

    // Comparison metrics at current time
    metrics: {
        leftParticleCount: 0,
        rightParticleCount: 0,
        leftTotalMass: 0,
        rightTotalMass: 0,
        leftMeanTemp: 0,
        rightMeanTemp: 0,
        leftMeanDens: 0,
        rightMeanDens: 0
    },

    // UI state
    isLoading: false,
    loadProgress: { left: 0, right: 0 },

    // Synchronized selection (shared between both views)
    selection: {
        mode: 'column',  // 'column' or 'plane'
        active: false,
        // Column selection
        column: {
            center: null,  // {x, y}
            radius: 0.3
        },
        // Plane selection
        plane: {
            zPosition: 0,
            thickness: 0.2
        },
        // Visualization meshes for both sides
        leftMesh: null,
        rightMesh: null,
        // Lagrangian tracking (track same particles across frames)
        trackingMode: 'lagrangian',  // 'spatial' or 'lagrangian'
        trackedIds: {
            left: new Set(),   // Particle IDs from left (cooling) dataset
            right: new Set()   // Particle IDs from right (adiabatic) dataset
        },
        selectionFrame: { left: 0, right: 0 },  // Frame at which selection was made
        showMesh: true  // Whether to show selection cylinder/plane
    },

    // Color mode and ranges for comparison view
    colorMode: 'density',
    colorRanges: {
        density: { min: 3, max: 6, isLog: true },
        temperature: { min: 1, max: 4, isLog: true },
        velocity: { min: -1, max: 2, isLog: true }
    }
};

// ============================================================
// Initialization
// ============================================================

async function initComparisonMode(comparisonId) {
    console.log('Initializing comparison mode:', comparisonId);

    // Find comparison config
    const response = await fetch('datasets.json');
    const data = await response.json();
    const comparison = data.comparisons?.find(c => c.id === comparisonId);

    if (!comparison) {
        console.error('Comparison not found:', comparisonId);
        return false;
    }

    const leftDataset = data.datasets.find(d => d.id === comparison.left);
    const rightDataset = data.datasets.find(d => d.id === comparison.right);

    if (!leftDataset || !rightDataset) {
        console.error('Dataset not found for comparison');
        return false;
    }

    COMPARISON.left.datasetId = leftDataset.id;
    COMPARISON.left.basePath = leftDataset.path;
    COMPARISON.right.datasetId = rightDataset.id;
    COMPARISON.right.basePath = rightDataset.path;

    // Set up CONFIG with dataset values (both datasets have same orbit params)
    // Use left dataset config as the reference
    if (leftDataset.config) {
        CONFIG.cloud_pos0 = leftDataset.config.cloud_pos0;
        CONFIG.cloud_vel0 = leftDataset.config.cloud_vel0;
        CONFIG.r_peri = leftDataset.config.r_peri;
    }
    // Common physics constants
    if (data.common) {
        CONFIG.G = data.common.G;
        CONFIG.M_BH = data.common.M_BH;
        CONFIG.tempConversion = data.common.tempConversion;
        CONFIG.densityToNH2 = data.common.densityToNH2;
        CONFIG.timeToMyr = data.common.timeToMyr;
    }
    // Compute derived orbital parameters (L_vec, L, p)
    if (typeof computeOrbitalParams === 'function') {
        computeOrbitalParams();
        console.log('Comparison mode: orbital params computed - r_peri:', CONFIG.r_peri, 'p:', CONFIG.p?.toFixed(4));
    }

    // Show comparison UI
    showComparisonUI();

    // Initialize dual renderers
    initDualRenderers();

    // Initialize selection handlers
    initComparisonSelection();

    // Initialize color range sliders
    initComparisonRangeSliders();

    // Initialize profile panel drag/resize
    initComparisonProfilePanel();

    // Scan for files and build time index (lightweight - only reads headers)
    COMPARISON.isLoading = true;
    updateComparisonLoadingUI();

    await Promise.all([
        scanDatasetFiles('left', leftDataset),
        scanDatasetFiles('right', rightDataset)
    ]);

    // Compute overlapping time range
    computeTimeRange();

    // Load first frames
    await Promise.all([
        loadFrameIfNeeded('left', 0),
        loadFrameIfNeeded('right', 0)
    ]);

    COMPARISON.isLoading = false;
    hideComparisonLoading();

    // Initialize to t=0
    COMPARISON.currentTime = COMPARISON.minTime;
    await syncToTime(COMPARISON.currentTime);

    COMPARISON.enabled = true;
    console.log('Comparison mode initialized');

    // Start render loop
    animateComparison();

    return true;
}

// ============================================================
// Data Loading (Lazy - scan headers first, load data on demand)
// ============================================================

async function scanDatasetFiles(side, dataset) {
    const basePath = dataset.path;
    const timeIndex = [];
    const snapshotFiles = [];

    console.log(`Scanning ${side} dataset: ${basePath}`);

    // Scan for snapshot files and read only headers to get times
    let fileIndex = 1;
    let consecutiveFails = 0;

    while (consecutiveFails < 3) {
        const filename = `snapshot_${String(fileIndex).padStart(4, '0')}.csv`;
        const url = `${basePath}/${filename}`;

        try {
            const response = await fetch(url);
            if (!response.ok) {
                consecutiveFails++;
                fileIndex++;
                continue;
            }

            // Only read first 2KB to get header with time info
            const reader = response.body.getReader();
            const { value } = await reader.read();
            reader.cancel();  // Cancel the rest of the download

            const headerText = new TextDecoder().decode(value.slice(0, 2000));
            const time = extractTimeFromHeader(headerText, fileIndex - 1);

            snapshotFiles.push(filename);
            timeIndex.push({
                time: time,
                frameIdx: fileIndex - 1,
                filename: filename
            });

            consecutiveFails = 0;
            COMPARISON.loadProgress[side] = fileIndex;
            updateComparisonLoadingUI();

        } catch (error) {
            consecutiveFails++;
        }

        fileIndex++;
        if (fileIndex > 300) break;  // Safety limit
    }

    COMPARISON[side].snapshotFiles = snapshotFiles;
    COMPARISON[side].timeIndex = timeIndex.sort((a, b) => a.time - b.time);
    COMPARISON[side].snapshotCache = {};

    console.log(`Scanned ${snapshotFiles.length} files for ${side}`);
}

function extractTimeFromHeader(headerText, fallbackIndex) {
    // Extract time from CSV header: "# Time (code): 1.234"
    const match = headerText.match(/# Time \(code\):\s*([\d.e+-]+)/);
    if (match) {
        return parseFloat(match[1]);
    }
    return fallbackIndex * 0.02;  // Fallback
}

async function loadFrameIfNeeded(side, frameIdx) {
    const data = COMPARISON[side];

    // Already cached?
    if (data.snapshotCache[frameIdx]) {
        return data.snapshotCache[frameIdx];
    }

    // Find the filename for this frame
    const timeEntry = data.timeIndex.find(t => t.frameIdx === frameIdx);
    if (!timeEntry) {
        console.warn(`Frame ${frameIdx} not found in ${side} time index`);
        return null;
    }

    const url = `${data.basePath}/${timeEntry.filename}`;
    console.log(`Loading ${side} frame ${frameIdx}: ${url}`);

    try {
        const response = await fetch(url);
        if (!response.ok) {
            console.error(`Failed to load ${url}`);
            return null;
        }

        const csvText = await response.text();
        const snapshot = parseComparisonSnapshot(csvText, frameIdx);

        // Cache management: remove oldest if cache is full
        const cacheKeys = Object.keys(data.snapshotCache).map(Number);
        if (cacheKeys.length >= COMPARISON.maxCacheSize) {
            // Remove the frame furthest from current
            const currentFrame = data.currentFrame;
            cacheKeys.sort((a, b) => Math.abs(b - currentFrame) - Math.abs(a - currentFrame));
            const toRemove = cacheKeys[0];
            delete data.snapshotCache[toRemove];
            console.log(`Cache evicted frame ${toRemove} from ${side}`);
        }

        data.snapshotCache[frameIdx] = snapshot;
        return snapshot;

    } catch (error) {
        console.error(`Error loading frame ${frameIdx}:`, error);
        return null;
    }
}

function getSnapshot(side, frameIdx) {
    return COMPARISON[side].snapshotCache[frameIdx] || null;
}

function parseComparisonSnapshot(csvText, frameIndex) {
    const lines = csvText.split('\n');

    // Extract time from header comments
    let time = frameIndex * 0.02;  // Default fallback
    for (const line of lines) {
        if (line.startsWith('# Time (code):')) {
            time = parseFloat(line.split(':')[1].trim());
            break;
        }
    }

    // Find header line (first non-comment line)
    let headerIdx = 0;
    while (headerIdx < lines.length && lines[headerIdx].startsWith('#')) {
        headerIdx++;
    }

    const headers = lines[headerIdx].split(',').map(h => h.trim());

    // Parse particles
    const particles = [];
    for (let i = headerIdx + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (!line) continue;

        const values = line.split(',').map(v => parseFloat(v.trim()));
        const particle = {};

        headers.forEach((h, idx) => {
            // Map CSV column names to our standard names
            const mapping = {
                'pos_x': 'x', 'pos_y': 'y', 'pos_z': 'z',
                'vel_x': 'vx', 'vel_y': 'vy', 'vel_z': 'vz',
                'dens': 'dens', 'pres': 'pres', 'ene': 'temp',
                'mass': 'mass', 'sound': 'sound', 'id': 'id'
            };
            const key = mapping[h] || h;
            particle[key] = values[idx];
        });

        // Convert internal energy to temperature if needed
        // T = (gamma - 1) * mu * m_H * u / k_B ≈ 186 * u for molecular gas
        if (particle.temp !== undefined && particle.temp > 100) {
            // Already seems like temperature or internal energy
            // If it's internal energy (ene), convert
            particle.temp = particle.temp * CONFIG.tempConversion;
        }

        // Compute velocity magnitude
        if (particle.vx !== undefined) {
            particle.vel_mag = Math.sqrt(
                particle.vx**2 + particle.vy**2 + particle.vz**2
            );
        }

        // Ensure ID
        if (particle.id === undefined) {
            particle.id = i - headerIdx - 1;
        }

        particle.is_ghost = 0;
        particles.push(particle);
    }

    return { time, particles, frameIndex };
}

// ============================================================
// Time Indexing (time index is built during scan, not here)
// ============================================================

function computeTimeRange() {
    const leftTimes = COMPARISON.left.timeIndex;
    const rightTimes = COMPARISON.right.timeIndex;

    if (leftTimes.length === 0 || rightTimes.length === 0) {
        console.error('No snapshots found');
        return;
    }

    // Find overlapping time range
    const leftMin = leftTimes[0].time;
    const leftMax = leftTimes[leftTimes.length - 1].time;
    const rightMin = rightTimes[0].time;
    const rightMax = rightTimes[rightTimes.length - 1].time;

    COMPARISON.minTime = Math.max(leftMin, rightMin);
    COMPARISON.maxTime = Math.min(leftMax, rightMax);

    console.log('Comparison time range:', COMPARISON.minTime.toFixed(3), 'to', COMPARISON.maxTime.toFixed(3));
    console.log('Left:', leftTimes.length, 'frames, Right:', rightTimes.length, 'frames');

    // Update timeline slider
    updateComparisonTimeline();
}

function getFrameAtTime(side, targetTime) {
    const timeIndex = COMPARISON[side].timeIndex;

    let bestIdx = 0;
    let bestDiff = Infinity;

    for (let i = 0; i < timeIndex.length; i++) {
        const diff = Math.abs(timeIndex[i].time - targetTime);
        if (diff < bestDiff) {
            bestDiff = diff;
            bestIdx = timeIndex[i].frameIdx;
        }
    }

    return bestIdx;
}

async function syncToTime(time) {
    COMPARISON.currentTime = Math.max(COMPARISON.minTime, Math.min(COMPARISON.maxTime, time));

    // Find nearest frames for both sides
    const leftFrame = getFrameAtTime('left', COMPARISON.currentTime);
    const rightFrame = getFrameAtTime('right', COMPARISON.currentTime);

    COMPARISON.left.currentFrame = leftFrame;
    COMPARISON.right.currentFrame = rightFrame;

    // Load frames if needed (lazy loading)
    await Promise.all([
        loadFrameIfNeeded('left', leftFrame),
        loadFrameIfNeeded('right', rightFrame)
    ]);

    // Update particle systems
    updateComparisonParticles('left');
    updateComparisonParticles('right');

    // Compute comparison metrics
    computeComparisonMetrics();

    // Update selection visuals, info and profiles if active
    if (COMPARISON.selection.active) {
        updateComparisonSelectionVisuals();  // Move cylinder/plane with tracked particles
        updateComparisonSelectionInfo();
        updateComparisonProfiles();
    }

    // Update UI
    updateComparisonTimeDisplay();
}

// ============================================================
// Dual Renderer Setup
// ============================================================

function initDualRenderers() {
    const leftContainer = document.getElementById('comparison-left-canvas');
    const rightContainer = document.getElementById('comparison-right-canvas');

    if (!leftContainer || !rightContainer) {
        console.error('Comparison containers not found');
        return;
    }

    // Left renderer (cooling)
    COMPARISON.left.scene = new THREE.Scene();
    COMPARISON.left.scene.background = new THREE.Color(0x0a0c14);

    COMPARISON.left.camera = new THREE.PerspectiveCamera(60, 1, 0.1, 1000);
    COMPARISON.left.camera.position.set(0, 0, 25);

    COMPARISON.left.renderer = new THREE.WebGLRenderer({ antialias: true });
    COMPARISON.left.renderer.setPixelRatio(window.devicePixelRatio);
    leftContainer.appendChild(COMPARISON.left.renderer.domElement);

    COMPARISON.left.controls = new THREE.OrbitControls(
        COMPARISON.left.camera,
        COMPARISON.left.renderer.domElement
    );
    COMPARISON.left.controls.enableDamping = true;

    // Right renderer (adiabatic)
    COMPARISON.right.scene = new THREE.Scene();
    COMPARISON.right.scene.background = new THREE.Color(0x0a0c14);

    COMPARISON.right.camera = new THREE.PerspectiveCamera(60, 1, 0.1, 1000);
    COMPARISON.right.camera.position.set(0, 0, 25);

    COMPARISON.right.renderer = new THREE.WebGLRenderer({ antialias: true });
    COMPARISON.right.renderer.setPixelRatio(window.devicePixelRatio);
    rightContainer.appendChild(COMPARISON.right.renderer.domElement);

    COMPARISON.right.controls = new THREE.OrbitControls(
        COMPARISON.right.camera,
        COMPARISON.right.renderer.domElement
    );
    COMPARISON.right.controls.enableDamping = true;

    // Add IMBH markers to both scenes
    addIMBHMarker(COMPARISON.left.scene);
    addIMBHMarker(COMPARISON.right.scene);

    // Add orbit lines to both scenes
    addOrbitLine(COMPARISON.left.scene);
    addOrbitLine(COMPARISON.right.scene);

    // Sync cameras: when left camera moves, update right camera
    COMPARISON.left.controls.addEventListener('change', syncRightCamera);

    // Initial resize
    resizeComparisonRenderers();
    window.addEventListener('resize', resizeComparisonRenderers);
}

function addIMBHMarker(scene) {
    const bhGeometry = new THREE.SphereGeometry(0.15, 32, 32);
    const bhMaterial = new THREE.MeshBasicMaterial({ color: 0xffffff });
    const bhMesh = new THREE.Mesh(bhGeometry, bhMaterial);
    bhMesh.position.set(0, 0, 0);
    scene.add(bhMesh);

    // Accretion ring
    const ringGeometry = new THREE.RingGeometry(0.2, 0.35, 32);
    const ringMaterial = new THREE.MeshBasicMaterial({
        color: 0xffaa00,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: 0.5
    });
    const ring = new THREE.Mesh(ringGeometry, ringMaterial);
    ring.rotation.x = Math.PI / 2;
    scene.add(ring);
}

function addOrbitLine(scene) {
    // Use CONFIG values (should be initialized from datasets.json)
    if (!CONFIG || !CONFIG.cloud_pos0 || !CONFIG.cloud_vel0) {
        console.warn('CONFIG not ready for orbit creation');
        return;
    }

    const points = [];
    const colors = [];

    const r0 = Math.sqrt(
        CONFIG.cloud_pos0[0]**2 +
        CONFIG.cloud_pos0[1]**2 +
        CONFIG.cloud_pos0[2]**2
    );

    const L_norm = [
        CONFIG.L_vec[0]/CONFIG.L,
        CONFIG.L_vec[1]/CONFIG.L,
        CONFIG.L_vec[2]/CONFIG.L
    ];

    const v_cross_L = [
        CONFIG.cloud_vel0[1] * CONFIG.L_vec[2] - CONFIG.cloud_vel0[2] * CONFIG.L_vec[1],
        CONFIG.cloud_vel0[2] * CONFIG.L_vec[0] - CONFIG.cloud_vel0[0] * CONFIG.L_vec[2],
        CONFIG.cloud_vel0[0] * CONFIG.L_vec[1] - CONFIG.cloud_vel0[1] * CONFIG.L_vec[0]
    ];
    const r_hat = [
        CONFIG.cloud_pos0[0]/r0,
        CONFIG.cloud_pos0[1]/r0,
        CONFIG.cloud_pos0[2]/r0
    ];
    const e_vec = [
        v_cross_L[0]/(CONFIG.G*CONFIG.M_BH) - r_hat[0],
        v_cross_L[1]/(CONFIG.G*CONFIG.M_BH) - r_hat[1],
        v_cross_L[2]/(CONFIG.G*CONFIG.M_BH) - r_hat[2]
    ];
    const e_mag = Math.sqrt(e_vec[0]**2 + e_vec[1]**2 + e_vec[2]**2);
    const e_hat = [e_vec[0]/e_mag, e_vec[1]/e_mag, e_vec[2]/e_mag];

    const p_hat = [
        L_norm[1]*e_hat[2] - L_norm[2]*e_hat[1],
        L_norm[2]*e_hat[0] - L_norm[0]*e_hat[2],
        L_norm[0]*e_hat[1] - L_norm[1]*e_hat[0]
    ];

    for (let i = 0; i <= 200; i++) {
        const theta = -Math.PI * 0.95 + i * (Math.PI * 1.9) / 200;
        const r = CONFIG.p / (1 + Math.cos(theta));

        if (r > 50 || r < 0) continue;

        const x_orb = r * Math.cos(theta);
        const y_orb = r * Math.sin(theta);

        const x = x_orb * e_hat[0] + y_orb * p_hat[0];
        const y = x_orb * e_hat[1] + y_orb * p_hat[1];
        const z = x_orb * e_hat[2] + y_orb * p_hat[2];

        points.push(new THREE.Vector3(x, y, z));

        const hue = 0.1 + 0.2 * (1 - Math.abs(theta) / Math.PI);
        colors.push(new THREE.Color().setHSL(hue, 1, 0.5));
    }

    if (points.length === 0) return;

    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const colorArray = new Float32Array(colors.length * 3);
    colors.forEach((c, i) => {
        colorArray[i*3] = c.r;
        colorArray[i*3+1] = c.g;
        colorArray[i*3+2] = c.b;
    });
    geometry.setAttribute('color', new THREE.BufferAttribute(colorArray, 3));

    const material = new THREE.LineBasicMaterial({
        vertexColors: true,
        linewidth: 2
    });

    const orbitLine = new THREE.Line(geometry, material);
    scene.add(orbitLine);

    // Pericenter marker
    const periPos = new THREE.Vector3(
        CONFIG.r_peri * e_hat[0],
        CONFIG.r_peri * e_hat[1],
        CONFIG.r_peri * e_hat[2]
    );
    const periGeom = new THREE.RingGeometry(0.08, 0.12, 32);
    const periMat = new THREE.MeshBasicMaterial({ color: 0x00ff88, side: THREE.DoubleSide });
    const periMarker = new THREE.Mesh(periGeom, periMat);
    periMarker.position.copy(periPos);
    periMarker.lookAt(0, 0, 0);
    scene.add(periMarker);

    console.log('Orbit line added to comparison scene');
}

function syncRightCamera() {
    if (!COMPARISON.enabled) return;

    COMPARISON.right.camera.position.copy(COMPARISON.left.camera.position);
    COMPARISON.right.camera.quaternion.copy(COMPARISON.left.camera.quaternion);
    COMPARISON.right.controls.target.copy(COMPARISON.left.controls.target);
    COMPARISON.right.controls.update();
}

function resizeComparisonRenderers() {
    const leftContainer = document.getElementById('comparison-left-canvas');
    const rightContainer = document.getElementById('comparison-right-canvas');

    if (!leftContainer || !rightContainer) return;

    const leftW = leftContainer.clientWidth;
    const leftH = leftContainer.clientHeight;
    const rightW = rightContainer.clientWidth;
    const rightH = rightContainer.clientHeight;

    if (COMPARISON.left.renderer && leftW > 0 && leftH > 0) {
        COMPARISON.left.renderer.setSize(leftW, leftH);
        COMPARISON.left.camera.aspect = leftW / leftH;
        COMPARISON.left.camera.updateProjectionMatrix();
    }

    if (COMPARISON.right.renderer && rightW > 0 && rightH > 0) {
        COMPARISON.right.renderer.setSize(rightW, rightH);
        COMPARISON.right.camera.aspect = rightW / rightH;
        COMPARISON.right.camera.updateProjectionMatrix();
    }
}

// ============================================================
// Particle Rendering
// ============================================================

function updateComparisonParticles(side) {
    const data = COMPARISON[side];
    const snapshot = data.snapshotCache[data.currentFrame];

    if (!snapshot || !data.scene) {
        console.log(`No snapshot for ${side} frame ${data.currentFrame}`);
        return;
    }

    // Remove old particle system
    if (data.particleSystem) {
        data.scene.remove(data.particleSystem);
        data.particleSystem.geometry.dispose();
        data.particleSystem.material.dispose();
    }

    const particles = snapshot.particles;
    const positions = new Float32Array(particles.length * 3);
    const colors = new Float32Array(particles.length * 3);

    // Compute COM velocity for color mode
    let comVx = 0, comVy = 0, comVz = 0, totalMass = 0;
    for (const p of particles) {
        comVx += p.vx * p.mass;
        comVy += p.vy * p.mass;
        comVz += p.vz * p.mass;
        totalMass += p.mass;
    }
    const comVel = totalMass > 0 ? {
        x: comVx / totalMass,
        y: comVy / totalMass,
        z: comVz / totalMass
    } : { x: 0, y: 0, z: 0 };

    for (let i = 0; i < particles.length; i++) {
        const p = particles[i];

        positions[i * 3] = p.x;
        positions[i * 3 + 1] = p.y;
        positions[i * 3 + 2] = p.z;

        // Check if particle is selected
        const selected = isParticleSelected(p, side);

        // Color by density (same as main viz), or highlight if selected
        const colorValue = getComparisonColorValue(p, comVel);
        const color = valueToColor(colorValue);

        if (selected) {
            // Brighten selected particles
            colors[i * 3] = Math.min(1, color.r + 0.3);
            colors[i * 3 + 1] = Math.min(1, color.g + 0.3);
            colors[i * 3 + 2] = Math.min(1, color.b + 0.3);
        } else if (COMPARISON.selection.active) {
            // Dim non-selected particles when selection is active
            colors[i * 3] = color.r * 0.3;
            colors[i * 3 + 1] = color.g * 0.3;
            colors[i * 3 + 2] = color.b * 0.3;
        } else {
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
        }
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

    const material = new THREE.PointsMaterial({
        size: 0.08,
        vertexColors: true,
        transparent: true,
        opacity: 0.8,
        sizeAttenuation: true
    });

    data.particleSystem = new THREE.Points(geometry, material);
    data.scene.add(data.particleSystem);
}

function getComparisonColorValue(particle, comVel) {
    const colorMode = COMPARISON.colorMode || 'density';
    const range = COMPARISON.colorRanges[colorMode] || { min: 3, max: 6, isLog: true };
    const isLog = range.isLog !== false; // default to log

    let rawValue;
    switch (colorMode) {
        case 'temperature':
            rawValue = particle.temp;
            break;
        case 'velocity':
            const dvx = particle.vx - comVel.x;
            const dvy = particle.vy - comVel.y;
            const dvz = particle.vz - comVel.z;
            rawValue = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            break;
        case 'density':
        default:
            rawValue = particle.dens * (CONFIG?.densityToNH2 || 20.3);
            break;
    }

    // Apply log transform if needed
    const value = isLog ? Math.log10(Math.max(rawValue, 1e-10)) : rawValue;

    return Math.max(0, Math.min(1, (value - range.min) / (range.max - range.min)));
}

// ============================================================
// Comparison Metrics
// ============================================================

function computeComparisonMetrics() {
    const leftSnap = COMPARISON.left.snapshotCache[COMPARISON.left.currentFrame];
    const rightSnap = COMPARISON.right.snapshotCache[COMPARISON.right.currentFrame];

    if (!leftSnap || !rightSnap) return;

    // Left (cooling) metrics
    const leftParticles = leftSnap.particles;
    let leftMass = 0, leftTempSum = 0, leftDensSum = 0;
    for (const p of leftParticles) {
        leftMass += p.mass;
        leftTempSum += p.temp;
        leftDensSum += p.dens * (CONFIG?.densityToNH2 || 20.3);
    }

    // Right (adiabatic) metrics
    const rightParticles = rightSnap.particles;
    let rightMass = 0, rightTempSum = 0, rightDensSum = 0;
    for (const p of rightParticles) {
        rightMass += p.mass;
        rightTempSum += p.temp;
        rightDensSum += p.dens * (CONFIG?.densityToNH2 || 20.3);
    }

    COMPARISON.metrics = {
        leftParticleCount: leftParticles.length,
        rightParticleCount: rightParticles.length,
        leftTotalMass: leftMass,
        rightTotalMass: rightMass,
        leftMeanTemp: leftTempSum / leftParticles.length,
        rightMeanTemp: rightTempSum / rightParticles.length,
        leftMeanDens: leftDensSum / leftParticles.length,
        rightMeanDens: rightDensSum / rightParticles.length
    };

    updateComparisonMetricsUI();
}

// ============================================================
// UI Management
// ============================================================

function showComparisonUI() {
    // Hide main viz, show comparison container
    const mainContainer = document.getElementById('container');
    const compContainer = document.getElementById('comparison-container');

    if (mainContainer) mainContainer.style.display = 'none';
    if (compContainer) compContainer.style.display = 'flex';

    // Hide regular controls
    const controls = document.getElementById('controls');
    if (controls) controls.style.display = 'none';
}

function hideComparisonUI() {
    const mainContainer = document.getElementById('container');
    const compContainer = document.getElementById('comparison-container');

    if (mainContainer) mainContainer.style.display = 'block';
    if (compContainer) compContainer.style.display = 'none';

    const controls = document.getElementById('controls');
    if (controls) controls.style.display = 'block';
}

function updateComparisonLoadingUI() {
    const loadingEl = document.getElementById('comparison-loading');
    if (!loadingEl) return;

    if (COMPARISON.isLoading) {
        loadingEl.style.display = 'flex';
        loadingEl.innerHTML = `
            <div class="spinner"></div>
            <div>Scanning files...</div>
            <div>Cooling: ${COMPARISON.loadProgress.left} files</div>
            <div>Adiabatic: ${COMPARISON.loadProgress.right} files</div>
        `;
    } else {
        loadingEl.style.display = 'none';
    }
}

function hideComparisonLoading() {
    const loadingEl = document.getElementById('comparison-loading');
    if (loadingEl) loadingEl.style.display = 'none';
}

function updateComparisonTimeline() {
    const slider = document.getElementById('comparison-timeline');
    if (!slider) return;

    slider.min = COMPARISON.minTime;
    slider.max = COMPARISON.maxTime;
    slider.step = (COMPARISON.maxTime - COMPARISON.minTime) / 100;
    slider.value = COMPARISON.currentTime;
}

function updateComparisonTimeDisplay() {
    const timeEl = document.getElementById('comparison-time-display');
    if (!timeEl) return;

    const timeMyr = COMPARISON.currentTime * (CONFIG?.timeToMyr || 0.978);
    const maxTimeMyr = COMPARISON.maxTime * (CONFIG?.timeToMyr || 0.978);

    timeEl.textContent = `t = ${timeMyr.toFixed(3)} / ${maxTimeMyr.toFixed(3)} Myr`;

    // Update individual panel times
    const leftTimeEl = document.getElementById('comparison-left-time');
    const rightTimeEl = document.getElementById('comparison-right-time');

    if (leftTimeEl) {
        const leftSnap = COMPARISON.left.snapshotCache[COMPARISON.left.currentFrame];
        leftTimeEl.textContent = `t = ${(leftSnap?.time * 0.978).toFixed(3)} Myr | N = ${leftSnap?.particles.length || 0}`;
    }

    if (rightTimeEl) {
        const rightSnap = COMPARISON.right.snapshotCache[COMPARISON.right.currentFrame];
        rightTimeEl.textContent = `t = ${(rightSnap?.time * 0.978).toFixed(3)} Myr | N = ${rightSnap?.particles.length || 0}`;
    }
}

function updateComparisonMetricsUI() {
    const metricsEl = document.getElementById('comparison-metrics');
    if (!metricsEl) return;

    const m = COMPARISON.metrics;
    const deltaN = m.leftParticleCount - m.rightParticleCount;
    const deltaM = m.leftTotalMass - m.rightTotalMass;
    const deltaT = m.leftMeanTemp - m.rightMeanTemp;

    metricsEl.innerHTML = `
        <div class="metric-row">
            <span class="metric-label">ΔN (particles)</span>
            <span class="metric-value ${deltaN < 0 ? 'negative' : ''}">${deltaN > 0 ? '+' : ''}${deltaN}</span>
        </div>
        <div class="metric-row">
            <span class="metric-label">ΔM (M☉)</span>
            <span class="metric-value ${deltaM < 0 ? 'negative' : ''}">${deltaM.toFixed(2)}</span>
        </div>
        <div class="metric-row">
            <span class="metric-label">ΔT_mean (K)</span>
            <span class="metric-value ${deltaT < 0 ? 'negative' : ''}">${deltaT.toFixed(0)}</span>
        </div>
        <div class="metric-row">
            <span class="metric-label">Cooling particles</span>
            <span class="metric-value">${m.leftParticleCount}</span>
        </div>
        <div class="metric-row">
            <span class="metric-label">Adiabatic particles</span>
            <span class="metric-value">${m.rightParticleCount}</span>
        </div>
    `;
}

// ============================================================
// Timeline Controls
// ============================================================

async function onComparisonTimelineChange(value) {
    await syncToTime(parseFloat(value));
}

let comparisonPlaying = false;
let comparisonPlayInterval = null;

function toggleComparisonPlay() {
    comparisonPlaying = !comparisonPlaying;

    const btn = document.getElementById('comparison-play-btn');
    if (btn) btn.textContent = comparisonPlaying ? '⏸' : '▶';

    if (comparisonPlaying) {
        const timeStep = (COMPARISON.maxTime - COMPARISON.minTime) / 100;

        const playStep = async () => {
            if (!comparisonPlaying) return;

            let newTime = COMPARISON.currentTime + timeStep;
            if (newTime > COMPARISON.maxTime) {
                newTime = COMPARISON.minTime;
            }
            await syncToTime(newTime);

            const slider = document.getElementById('comparison-timeline');
            if (slider) slider.value = COMPARISON.currentTime;

            if (comparisonPlaying) {
                comparisonPlayInterval = setTimeout(playStep, 150);  // Slower for lazy loading
            }
        };

        playStep();
    } else {
        if (comparisonPlayInterval) {
            clearTimeout(comparisonPlayInterval);
            comparisonPlayInterval = null;
        }
    }
}

async function stepComparisonFrame(direction) {
    const timeStep = (COMPARISON.maxTime - COMPARISON.minTime) / 50;
    let newTime = COMPARISON.currentTime + direction * timeStep;
    newTime = Math.max(COMPARISON.minTime, Math.min(COMPARISON.maxTime, newTime));
    await syncToTime(newTime);

    const slider = document.getElementById('comparison-timeline');
    if (slider) slider.value = COMPARISON.currentTime;
}

// ============================================================
// Render Loop
// ============================================================

function animateComparison() {
    if (!COMPARISON.enabled) return;

    requestAnimationFrame(animateComparison);

    if (COMPARISON.left.controls) COMPARISON.left.controls.update();
    if (COMPARISON.right.controls) COMPARISON.right.controls.update();

    if (COMPARISON.left.renderer && COMPARISON.left.scene && COMPARISON.left.camera) {
        COMPARISON.left.renderer.render(COMPARISON.left.scene, COMPARISON.left.camera);
    }

    if (COMPARISON.right.renderer && COMPARISON.right.scene && COMPARISON.right.camera) {
        COMPARISON.right.renderer.render(COMPARISON.right.scene, COMPARISON.right.camera);
    }
}

// ============================================================
// Exit Comparison Mode
// ============================================================

function exitComparisonMode() {
    COMPARISON.enabled = false;

    // Stop playback
    if (comparisonPlayInterval) {
        clearTimeout(comparisonPlayInterval);
        comparisonPlayInterval = null;
    }
    comparisonPlaying = false;

    // Clean up renderers
    if (COMPARISON.left.renderer) {
        COMPARISON.left.renderer.dispose();
        COMPARISON.left.renderer.domElement.remove();
    }
    if (COMPARISON.right.renderer) {
        COMPARISON.right.renderer.dispose();
        COMPARISON.right.renderer.domElement.remove();
    }

    // Clear caches to free memory
    COMPARISON.left.snapshotCache = {};
    COMPARISON.left.snapshotFiles = [];
    COMPARISON.left.timeIndex = [];
    COMPARISON.right.snapshotCache = {};
    COMPARISON.right.snapshotFiles = [];
    COMPARISON.right.timeIndex = [];

    // Reset references
    COMPARISON.left.renderer = null;
    COMPARISON.right.renderer = null;
    COMPARISON.left.scene = null;
    COMPARISON.right.scene = null;

    // Clear selection
    clearComparisonSelection();

    // Restore main UI
    hideComparisonUI();

    console.log('Exited comparison mode');
}

// ============================================================
// Synchronized Selection
// ============================================================

function initComparisonSelection() {
    // Add click handlers to both canvases
    const leftCanvas = document.getElementById('comparison-left-canvas');
    const rightCanvas = document.getElementById('comparison-right-canvas');

    if (leftCanvas) {
        leftCanvas.addEventListener('click', (e) => onComparisonCanvasClick(e, 'left'));
    }
    if (rightCanvas) {
        rightCanvas.addEventListener('click', (e) => onComparisonCanvasClick(e, 'right'));
    }

    console.log('Comparison selection handlers initialized');
}

function onComparisonCanvasClick(event, side) {
    if (!COMPARISON.enabled) return;

    // Skip if clicking on UI elements
    if (event.target.closest('button, select, input')) return;

    const data = COMPARISON[side];
    if (!data.renderer) return;

    const container = event.target.closest('#comparison-left-canvas, #comparison-right-canvas');
    if (!container) return;

    const rect = container.getBoundingClientRect();
    const mouse = new THREE.Vector2(
        ((event.clientX - rect.left) / rect.width) * 2 - 1,
        -((event.clientY - rect.top) / rect.height) * 2 + 1
    );

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(mouse, data.camera);

    // Try to hit particles
    if (data.particleSystem) {
        raycaster.params.Points.threshold = 0.3;
        const intersects = raycaster.intersectObject(data.particleSystem);

        if (intersects.length > 0) {
            const point = intersects[0].point;

            if (COMPARISON.selection.mode === 'column') {
                // Set column center
                COMPARISON.selection.column.center = { x: point.x, y: point.y };
                COMPARISON.selection.active = true;

                console.log(`Column selection at (${point.x.toFixed(2)}, ${point.y.toFixed(2)})`);
            } else {
                // Set plane z position
                COMPARISON.selection.plane.zPosition = point.z;
                COMPARISON.selection.active = true;

                console.log(`Plane selection at z = ${point.z.toFixed(2)}`);
            }

            // Capture particle IDs for Lagrangian tracking
            captureTrackedParticleIds();

            updateComparisonSelectionVisuals();
            updateComparisonParticlesWithSelection();
        }
    }
}

function updateComparisonSelectionVisuals() {
    // Remove old meshes
    removeComparisonSelectionMeshes();

    if (!COMPARISON.selection.active) return;
    if (!COMPARISON.selection.showMesh) return;  // Skip if mesh hidden

    if (COMPARISON.selection.mode === 'column') {
        createComparisonColumnMeshes();
    } else {
        createComparisonPlaneMeshes();
    }
}

function toggleComparisonSelectionMesh() {
    COMPARISON.selection.showMesh = !COMPARISON.selection.showMesh;

    // Update button appearance
    const btn = document.getElementById('comp-show-mesh');
    if (btn) {
        if (COMPARISON.selection.showMesh) {
            btn.textContent = 'Mesh ✓';
            btn.style.background = '#446688';
        } else {
            btn.textContent = 'Mesh ✗';
            btn.style.background = '#334';
        }
    }

    // Update visuals
    updateComparisonSelectionVisuals();
}

function removeComparisonSelectionMeshes() {
    if (COMPARISON.selection.leftMesh) {
        COMPARISON.left.scene?.remove(COMPARISON.selection.leftMesh);
        COMPARISON.selection.leftMesh.geometry?.dispose();
        COMPARISON.selection.leftMesh.material?.dispose();
        COMPARISON.selection.leftMesh = null;
    }
    if (COMPARISON.selection.rightMesh) {
        COMPARISON.right.scene?.remove(COMPARISON.selection.rightMesh);
        COMPARISON.selection.rightMesh.geometry?.dispose();
        COMPARISON.selection.rightMesh.material?.dispose();
        COMPARISON.selection.rightMesh = null;
    }
}

// Compute centroid of tracked particles for a given side
function computeTrackedParticleCentroid(side) {
    const snap = COMPARISON[side].snapshotCache[COMPARISON[side].currentFrame];
    if (!snap) return null;

    const trackedIds = COMPARISON.selection.trackedIds[side];
    if (trackedIds.size === 0) return null;

    let sumX = 0, sumY = 0, sumZ = 0, count = 0;
    for (const p of snap.particles) {
        if (trackedIds.has(p.id)) {
            sumX += p.x;
            sumY += p.y;
            sumZ += p.z;
            count++;
        }
    }

    if (count === 0) return null;
    return { x: sumX / count, y: sumY / count, z: sumZ / count };
}

function createComparisonColumnMeshes() {
    const { center, radius } = COMPARISON.selection.column;
    if (!center) return;

    // Estimate column height from data
    const height = 10;  // Default height

    const geometry = new THREE.CylinderGeometry(radius, radius, height, 32, 1, true);
    const material = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.15,
        side: THREE.DoubleSide,
        depthWrite: false
    });

    // In Lagrangian mode, position cylinder at tracked particles' centroid
    const isLagrangian = COMPARISON.selection.trackingMode === 'lagrangian';

    // Left scene
    if (COMPARISON.left.scene) {
        COMPARISON.selection.leftMesh = new THREE.Mesh(geometry.clone(), material.clone());
        COMPARISON.selection.leftMesh.rotation.x = Math.PI / 2;

        let posX = center.x, posY = center.y;
        if (isLagrangian) {
            const centroid = computeTrackedParticleCentroid('left');
            if (centroid) {
                posX = centroid.x;
                posY = centroid.y;
            }
        }
        COMPARISON.selection.leftMesh.position.set(posX, posY, 0);
        COMPARISON.left.scene.add(COMPARISON.selection.leftMesh);
    }

    // Right scene
    if (COMPARISON.right.scene) {
        COMPARISON.selection.rightMesh = new THREE.Mesh(geometry.clone(), material.clone());
        COMPARISON.selection.rightMesh.rotation.x = Math.PI / 2;

        let posX = center.x, posY = center.y;
        if (isLagrangian) {
            const centroid = computeTrackedParticleCentroid('right');
            if (centroid) {
                posX = centroid.x;
                posY = centroid.y;
            }
        }
        COMPARISON.selection.rightMesh.position.set(posX, posY, 0);
        COMPARISON.right.scene.add(COMPARISON.selection.rightMesh);
    }
}

function createComparisonPlaneMeshes() {
    const { zPosition, thickness } = COMPARISON.selection.plane;

    const geometry = new THREE.PlaneGeometry(20, 20);
    const material = new THREE.MeshBasicMaterial({
        color: 0x00ff88,
        transparent: true,
        opacity: 0.15,
        side: THREE.DoubleSide,
        depthWrite: false
    });

    // In Lagrangian mode, position plane at tracked particles' centroid
    const isLagrangian = COMPARISON.selection.trackingMode === 'lagrangian';

    // Left scene
    if (COMPARISON.left.scene) {
        COMPARISON.selection.leftMesh = new THREE.Mesh(geometry.clone(), material.clone());

        let posZ = zPosition;
        if (isLagrangian) {
            const centroid = computeTrackedParticleCentroid('left');
            if (centroid) {
                posZ = centroid.z;
                // Also move plane x,y to follow the cloud
                COMPARISON.selection.leftMesh.position.set(centroid.x, centroid.y, posZ);
            } else {
                COMPARISON.selection.leftMesh.position.set(0, 0, posZ);
            }
        } else {
            COMPARISON.selection.leftMesh.position.set(0, 0, posZ);
        }
        COMPARISON.left.scene.add(COMPARISON.selection.leftMesh);
    }

    // Right scene
    if (COMPARISON.right.scene) {
        COMPARISON.selection.rightMesh = new THREE.Mesh(geometry.clone(), material.clone());

        let posZ = zPosition;
        if (isLagrangian) {
            const centroid = computeTrackedParticleCentroid('right');
            if (centroid) {
                posZ = centroid.z;
                COMPARISON.selection.rightMesh.position.set(centroid.x, centroid.y, posZ);
            } else {
                COMPARISON.selection.rightMesh.position.set(0, 0, posZ);
            }
        } else {
            COMPARISON.selection.rightMesh.position.set(0, 0, posZ);
        }
        COMPARISON.right.scene.add(COMPARISON.selection.rightMesh);
    }
}

function updateComparisonParticlesWithSelection() {
    // Re-render particles with selection highlighting
    updateComparisonParticles('left');
    updateComparisonParticles('right');

    // Update selection visualization (moves with tracked particles in Lagrangian mode)
    updateComparisonSelectionVisuals();

    // Update selection info
    updateComparisonSelectionInfo();

    // Update profile plots
    updateComparisonProfiles();
}

// Check if particle is within spatial selection (column or plane)
function isParticleInSpatialSelection(particle) {
    if (COMPARISON.selection.mode === 'column') {
        const { center, radius } = COMPARISON.selection.column;
        if (!center) return false;
        const dx = particle.x - center.x;
        const dy = particle.y - center.y;
        return (dx * dx + dy * dy) <= (radius * radius);
    } else {
        const { zPosition, thickness } = COMPARISON.selection.plane;
        const halfThick = thickness / 2;
        return particle.z >= (zPosition - halfThick) && particle.z <= (zPosition + halfThick);
    }
}

// Capture particle IDs at current frame for Lagrangian tracking
function captureTrackedParticleIds() {
    // Clear previous tracked IDs
    COMPARISON.selection.trackedIds.left.clear();
    COMPARISON.selection.trackedIds.right.clear();

    // Record selection frame
    COMPARISON.selection.selectionFrame.left = COMPARISON.left.currentFrame;
    COMPARISON.selection.selectionFrame.right = COMPARISON.right.currentFrame;

    // Get particles from current frame
    const leftSnap = COMPARISON.left.snapshotCache[COMPARISON.left.currentFrame];
    const rightSnap = COMPARISON.right.snapshotCache[COMPARISON.right.currentFrame];

    // Capture IDs of particles within spatial selection
    if (leftSnap) {
        for (const p of leftSnap.particles) {
            if (isParticleInSpatialSelection(p)) {
                COMPARISON.selection.trackedIds.left.add(p.id);
            }
        }
    }

    if (rightSnap) {
        for (const p of rightSnap.particles) {
            if (isParticleInSpatialSelection(p)) {
                COMPARISON.selection.trackedIds.right.add(p.id);
            }
        }
    }

    console.log(`Captured ${COMPARISON.selection.trackedIds.left.size} left, ${COMPARISON.selection.trackedIds.right.size} right particles for tracking`);
}

// Check if particle is selected (uses Lagrangian or spatial mode)
function isParticleSelected(particle, side) {
    if (!COMPARISON.selection.active) return false;

    if (COMPARISON.selection.trackingMode === 'lagrangian') {
        // Use tracked particle IDs
        const trackedSet = side === 'left'
            ? COMPARISON.selection.trackedIds.left
            : COMPARISON.selection.trackedIds.right;
        return trackedSet.has(particle.id);
    } else {
        // Use spatial selection
        return isParticleInSpatialSelection(particle);
    }
}

function updateComparisonSelectionInfo() {
    const infoEl = document.getElementById('comparison-selection-info');
    if (!infoEl) return;

    if (!COMPARISON.selection.active) {
        infoEl.innerHTML = '<span style="color: #666;">Click on particles to select</span>';
        return;
    }

    // Count selected particles in both datasets
    const leftSnap = COMPARISON.left.snapshotCache[COMPARISON.left.currentFrame];
    const rightSnap = COMPARISON.right.snapshotCache[COMPARISON.right.currentFrame];

    let leftCount = 0, rightCount = 0;
    let leftMass = 0, rightMass = 0;

    if (leftSnap) {
        for (const p of leftSnap.particles) {
            if (isParticleSelected(p, 'left')) {
                leftCount++;
                leftMass += p.mass;
            }
        }
    }

    if (rightSnap) {
        for (const p of rightSnap.particles) {
            if (isParticleSelected(p, 'right')) {
                rightCount++;
                rightMass += p.mass;
            }
        }
    }

    const modeLabel = COMPARISON.selection.mode === 'column' ? 'Column' : 'Plane';
    const trackLabel = COMPARISON.selection.trackingMode === 'lagrangian' ? 'Lagrangian' : 'Spatial';
    const deltaN = leftCount - rightCount;
    const deltaM = leftMass - rightMass;

    // For Lagrangian mode, show initial vs current count
    let trackInfo = '';
    if (COMPARISON.selection.trackingMode === 'lagrangian') {
        const initLeft = COMPARISON.selection.trackedIds.left.size;
        const initRight = COMPARISON.selection.trackedIds.right.size;
        trackInfo = `<span style="color: #aaa; font-size: 10px;">(tracking ${initLeft}/${initRight} IDs)</span>`;
    }

    infoEl.innerHTML = `
        <div style="display: flex; gap: 16px; flex-wrap: wrap; align-items: center;">
            <span style="color: #88ddff;">${modeLabel} [${trackLabel}]</span>
            ${trackInfo}
            <span>Cooling: <b style="color: #66ddff;">${leftCount}</b> (${leftMass.toFixed(1)} M☉)</span>
            <span>Adiabatic: <b style="color: #ffcc88;">${rightCount}</b> (${rightMass.toFixed(1)} M☉)</span>
            <span>ΔN: <b style="color: ${deltaN < 0 ? '#faa' : '#afa'};">${deltaN > 0 ? '+' : ''}${deltaN}</b></span>
            <span>ΔM: <b style="color: ${deltaM < 0 ? '#faa' : '#afa'};">${deltaM.toFixed(2)} M☉</b></span>
        </div>
    `;
}

function setComparisonSelectionMode(mode) {
    COMPARISON.selection.mode = mode;
    COMPARISON.selection.active = false;
    removeComparisonSelectionMeshes();
    updateComparisonSelectionInfo();

    // Update mode button styles
    const colBtn = document.getElementById('comp-sel-column');
    const planeBtn = document.getElementById('comp-sel-plane');

    if (colBtn) colBtn.style.background = mode === 'column' ? '#446688' : '#334';
    if (planeBtn) planeBtn.style.background = mode === 'plane' ? '#446688' : '#334';

    console.log('Comparison selection mode:', mode);
}

function setComparisonTrackingMode(mode) {
    COMPARISON.selection.trackingMode = mode;

    // Update button styles
    const lagBtn = document.getElementById('comp-track-lagrangian');
    const spatBtn = document.getElementById('comp-track-spatial');

    if (lagBtn) lagBtn.style.background = mode === 'lagrangian' ? '#446688' : '#334';
    if (spatBtn) spatBtn.style.background = mode === 'spatial' ? '#446688' : '#334';

    // If switching to lagrangian and we have an active selection, recapture IDs
    if (mode === 'lagrangian' && COMPARISON.selection.active) {
        captureTrackedParticleIds();
    }

    // Re-render particles with updated selection mode
    if (COMPARISON.selection.active) {
        updateComparisonParticlesWithSelection();
    }

    console.log('Comparison tracking mode:', mode);
}

function updateComparisonSelectionRadius(value) {
    COMPARISON.selection.column.radius = parseFloat(value);

    const display = document.getElementById('comp-sel-radius-value');
    if (display) display.textContent = value + ' pc';

    if (COMPARISON.selection.active && COMPARISON.selection.mode === 'column') {
        // Recapture particle IDs if in Lagrangian mode
        if (COMPARISON.selection.trackingMode === 'lagrangian') {
            captureTrackedParticleIds();
        }
        updateComparisonSelectionVisuals();
        updateComparisonParticlesWithSelection();
    }
}

function updateComparisonSelectionThickness(value) {
    COMPARISON.selection.plane.thickness = parseFloat(value);

    const display = document.getElementById('comp-sel-thickness-value');
    if (display) display.textContent = value + ' pc';

    if (COMPARISON.selection.active && COMPARISON.selection.mode === 'plane') {
        // Recapture particle IDs if in Lagrangian mode
        if (COMPARISON.selection.trackingMode === 'lagrangian') {
            captureTrackedParticleIds();
        }
        updateComparisonParticlesWithSelection();
    }
}

function clearComparisonSelection() {
    COMPARISON.selection.active = false;
    COMPARISON.selection.column.center = null;
    // Clear tracked particle IDs
    COMPARISON.selection.trackedIds.left.clear();
    COMPARISON.selection.trackedIds.right.clear();
    removeComparisonSelectionMeshes();
    updateComparisonParticles('left');
    updateComparisonParticles('right');
    updateComparisonSelectionInfo();
    updateComparisonProfiles();  // Hide profile panel
}

// ============================================================
// Color Mode and Range Control
// ============================================================

function onComparisonColorModeChange(mode) {
    COMPARISON.colorMode = mode;

    // Update range sliders to match mode's range
    updateComparisonRangeSliders();

    // Refresh visualization
    syncToTime(COMPARISON.currentTime);
}

function onComparisonColormapChange(colormapName) {
    if (setColormap(colormapName)) {
        // Refresh visualization with new colormap
        updateComparisonParticles('left');
        updateComparisonParticles('right');
    }
}

function onComparisonRangeChange() {
    const minSlider = document.getElementById('comp-range-min');
    const maxSlider = document.getElementById('comp-range-max');

    if (!minSlider || !maxSlider) return;

    const minVal = parseFloat(minSlider.value);
    const maxVal = parseFloat(maxSlider.value);

    // Update display
    document.getElementById('comp-range-min-value').textContent = minVal.toFixed(1);
    document.getElementById('comp-range-max-value').textContent = maxVal.toFixed(1);

    // Ensure min < max
    const actualMin = Math.min(minVal, maxVal);
    const actualMax = Math.max(minVal, maxVal);

    // Update current color mode's range
    COMPARISON.colorRanges[COMPARISON.colorMode].min = actualMin;
    COMPARISON.colorRanges[COMPARISON.colorMode].max = actualMax;

    // Refresh visualization
    updateComparisonParticles('left');
    updateComparisonParticles('right');
}

function updateComparisonRangeSliders() {
    const range = COMPARISON.colorRanges[COMPARISON.colorMode];
    const isLog = range.isLog;

    const minSlider = document.getElementById('comp-range-min');
    const maxSlider = document.getElementById('comp-range-max');
    const label = document.getElementById('comp-range-label');

    if (!minSlider || !maxSlider) return;

    // Update slider attributes based on log vs linear
    if (isLog) {
        minSlider.min = -2;
        minSlider.max = 8;
        minSlider.step = 0.5;
        maxSlider.min = -2;
        maxSlider.max = 8;
        maxSlider.step = 0.5;
        if (label) label.textContent = 'Range (log₁₀):';
    } else {
        minSlider.min = 0;
        minSlider.max = 100;
        minSlider.step = 1;
        maxSlider.min = 0;
        maxSlider.max = 100;
        maxSlider.step = 1;
        if (label) label.textContent = 'Range:';
    }

    minSlider.value = range.min;
    maxSlider.value = range.max;
    document.getElementById('comp-range-min-value').textContent = range.min.toFixed(1);
    document.getElementById('comp-range-max-value').textContent = range.max.toFixed(1);

    // Update log toggle button
    const logBtn = document.getElementById('comp-log-toggle');
    if (logBtn) {
        if (isLog) {
            logBtn.textContent = 'Log';
            logBtn.style.background = '#446688';
        } else {
            logBtn.textContent = 'Linear';
            logBtn.style.background = '#664488';
        }
    }
}

function toggleComparisonLogScale() {
    const range = COMPARISON.colorRanges[COMPARISON.colorMode];
    const wasLog = range.isLog;

    // Toggle the scale
    range.isLog = !wasLog;

    // Convert range values between log and linear
    if (wasLog) {
        // Converting from log to linear: 10^x
        range.min = Math.pow(10, range.min);
        range.max = Math.pow(10, range.max);
    } else {
        // Converting from linear to log: log10(x)
        range.min = range.min > 0 ? Math.log10(range.min) : 0;
        range.max = range.max > 0 ? Math.log10(range.max) : 1;
    }

    // Update UI
    updateComparisonRangeSliders();

    // Refresh visualization
    updateComparisonParticles('left');
    updateComparisonParticles('right');
}

function initComparisonRangeSliders() {
    // Set initial values based on default color mode
    updateComparisonRangeSliders();

    // Sync dropdown with state
    const dropdown = document.getElementById('comparison-color-mode');
    if (dropdown) {
        dropdown.value = COMPARISON.colorMode;
    }
}

// ============================================================
// Profile Plots for Selected Particles
// ============================================================

const COMP_PROFILE_CONFIG = {
    density: { label: 'log n_H₂ (cm⁻³)', color: '#ff9966', isLog: true,
               getValue: (p) => Math.log10(Math.max(p.dens * 20.3, 1)) },
    temperature: { label: 'log T (K)', color: '#ff6688', isLog: true,
                   getValue: (p) => Math.log10(Math.max(p.temp, 1)) },
    mach: { label: 'Mach', color: '#ffff66', isLog: false,
            getValue: (p) => p.mach || 0 },
    entropy: { label: 'log s', color: '#66aaff', isLog: true,
               getValue: (p) => {
                   const gamma = 5/3;
                   const s = p.pres / Math.pow(p.dens, gamma);
                   return s > 0 ? Math.log10(s) : -10;
               }},
    pressure: { label: 'log P', color: '#ff66ff', isLog: true,
                getValue: (p) => p.pres > 0 ? Math.log10(p.pres) : -10 }
};

function updateComparisonProfiles() {
    if (!COMPARISON.selection.active) {
        hideComparisonProfilePanel();
        return;
    }

    showComparisonProfilePanel();

    // Get selected particles from both datasets
    const leftSnap = COMPARISON.left.snapshotCache[COMPARISON.left.currentFrame];
    const rightSnap = COMPARISON.right.snapshotCache[COMPARISON.right.currentFrame];

    const leftSelected = leftSnap ? getComparisonSelectedParticles(leftSnap.particles, 'left') : [];
    const rightSelected = rightSnap ? getComparisonSelectedParticles(rightSnap.particles, 'right') : [];

    // Get variable config
    const varKey = document.getElementById('comp-profile-var')?.value || 'density';
    const config = COMP_PROFILE_CONFIG[varKey];

    // Compute shared range for both plots (for fair comparison)
    const allSelected = [...leftSelected, ...rightSelected];
    let yMin = Infinity, yMax = -Infinity;
    for (const p of allSelected) {
        const val = config.getValue(p);
        if (isFinite(val)) {
            if (val < yMin) yMin = val;
            if (val > yMax) yMax = val;
        }
    }

    // Handle empty selection or invalid range
    if (!isFinite(yMin) || !isFinite(yMax)) {
        yMin = 0;
        yMax = 1;
    }

    // Add padding
    const yPad = (yMax - yMin) * 0.1 || 0.5;
    yMin -= yPad;
    yMax += yPad;

    // Draw both profiles
    drawComparisonProfilePlot('comp-profile-left', leftSelected, config, yMin, yMax, '#66ddff');
    drawComparisonProfilePlot('comp-profile-right', rightSelected, config, yMin, yMax, '#ffcc88');
}

function getComparisonSelectedParticles(particles, side) {
    if (!particles) return [];
    return particles.filter(p => isParticleSelected(p, side));
}

function drawComparisonProfilePlot(canvasId, particles, config, fixedYMin, fixedYMax, accentColor) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;

    // Get display size and set canvas buffer size for sharp rendering
    const rect = canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    const W = rect.width;
    const H = rect.height;

    // Update canvas internal size if needed
    if (canvas.width !== W * dpr || canvas.height !== H * dpr) {
        canvas.width = W * dpr;
        canvas.height = H * dpr;
    }

    const ctx = canvas.getContext('2d');
    ctx.setTransform(1, 0, 0, 1, 0, 0);  // Reset transform
    ctx.scale(dpr, dpr);

    // Clear
    ctx.fillStyle = '#0a0c12';
    ctx.fillRect(0, 0, W, H);

    if (particles.length === 0) {
        ctx.fillStyle = '#666';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('No particles selected', W/2, H/2);
        return;
    }

    const margin = { left: 50, right: 15, top: 20, bottom: 30 };
    const plotW = W - margin.left - margin.right;
    const plotH = H - margin.top - margin.bottom;

    // Determine x-axis based on selection mode
    const isColumn = COMPARISON.selection.mode === 'column';
    const xKey = isColumn ? 'z' : 'x';
    const xLabel = isColumn ? 'z (pc)' : 'x (pc)';

    // Get data
    const xVals = particles.map(p => p[xKey]);
    const yVals = particles.map(p => config.getValue(p));

    // X range
    let xMin = Math.min(...xVals.filter(v => isFinite(v)));
    let xMax = Math.max(...xVals.filter(v => isFinite(v)));
    const xPad = (xMax - xMin) * 0.05 || 1;
    xMin -= xPad;
    xMax += xPad;

    const xRange = xMax - xMin || 1;
    const yRange = fixedYMax - fixedYMin || 1;

    // Draw grid
    ctx.strokeStyle = '#222';
    ctx.lineWidth = 0.5;
    for (let i = 0; i <= 4; i++) {
        const gx = margin.left + plotW * i / 4;
        const gy = margin.top + plotH * i / 4;
        ctx.beginPath();
        ctx.moveTo(gx, margin.top);
        ctx.lineTo(gx, margin.top + plotH);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(margin.left, gy);
        ctx.lineTo(margin.left + plotW, gy);
        ctx.stroke();
    }

    // Plot border
    ctx.strokeStyle = '#334';
    ctx.lineWidth = 1;
    ctx.strokeRect(margin.left, margin.top, plotW, plotH);

    // Draw points
    ctx.fillStyle = config.color;
    ctx.globalAlpha = 0.4;
    for (let i = 0; i < particles.length; i++) {
        const px = margin.left + plotW * (xVals[i] - xMin) / xRange;
        const py = margin.top + plotH * (1 - (yVals[i] - fixedYMin) / yRange);
        if (isFinite(px) && isFinite(py) && px >= margin.left && px <= margin.left + plotW) {
            ctx.beginPath();
            ctx.arc(px, Math.max(margin.top, Math.min(margin.top + plotH, py)), 2, 0, Math.PI * 2);
            ctx.fill();
        }
    }
    ctx.globalAlpha = 1;

    // Labels
    ctx.fillStyle = config.color;
    ctx.font = 'bold 10px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(config.label, margin.left + 5, margin.top + 14);

    // Particle count
    ctx.fillStyle = accentColor;
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(`N = ${particles.length}`, W - margin.right - 5, margin.top + 12);

    // X-axis label
    ctx.fillStyle = '#aaa';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(xLabel, margin.left + plotW / 2, H - 8);

    // X-axis ticks
    ctx.fillStyle = '#888';
    ctx.font = '8px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(xMin.toFixed(1), margin.left, H - 18);
    ctx.fillText(xMax.toFixed(1), margin.left + plotW, H - 18);

    // Y-axis ticks
    ctx.textAlign = 'right';
    ctx.fillText(fixedYMax.toFixed(1), margin.left - 4, margin.top + 10);
    ctx.fillText(fixedYMin.toFixed(1), margin.left - 4, margin.top + plotH);
}

function drawComparisonMedian(ctx, x0, y0, w, h, xVals, yVals, xMin, xMax, yMin, yMax, color) {
    const sorted = xVals.map((x, i) => ({ x, y: yVals[i] }))
        .filter(p => isFinite(p.x) && isFinite(p.y))
        .sort((a, b) => a.x - b.x);

    if (sorted.length < 5) return;

    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const numBins = Math.min(20, Math.floor(sorted.length / 3));
    const binSize = sorted.length / numBins;

    ctx.strokeStyle = color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();

    let started = false;
    for (let i = 0; i < numBins; i++) {
        const binData = sorted.slice(Math.floor(i * binSize), Math.floor((i + 1) * binSize));
        if (binData.length === 0) continue;

        const xMed = binData[Math.floor(binData.length / 2)].x;
        const yValues = binData.map(p => p.y).sort((a, b) => a - b);
        const yMed = yValues[Math.floor(yValues.length / 2)];

        const px = x0 + w * (xMed - xMin) / xRange;
        const py = y0 + h * (1 - (yMed - yMin) / yRange);
        const clippedPy = Math.max(y0, Math.min(y0 + h, py));

        if (!started) { ctx.moveTo(px, clippedPy); started = true; }
        else { ctx.lineTo(px, clippedPy); }
    }
    ctx.stroke();
}

// ============================================================
// Profile Panel Drag & Resize
// ============================================================

let profilePanelDrag = { active: false, startX: 0, startY: 0, startLeft: 0, startTop: 0 };
let profilePanelResize = { active: false, startX: 0, startY: 0, startW: 0, startH: 0 };

function initComparisonProfilePanel() {
    const panel = document.getElementById('comparison-profile-panel');
    const header = document.getElementById('comp-profile-header');
    const resizeHandle = document.getElementById('comp-profile-resize');

    if (!panel || !header || !resizeHandle) return;

    // Drag functionality
    header.addEventListener('mousedown', (e) => {
        if (e.target.tagName === 'SELECT' || e.target.tagName === 'BUTTON') return;
        profilePanelDrag.active = true;
        profilePanelDrag.startX = e.clientX;
        profilePanelDrag.startY = e.clientY;
        profilePanelDrag.startLeft = panel.offsetLeft;
        profilePanelDrag.startTop = panel.offsetTop;
        e.preventDefault();
    });

    // Resize functionality
    resizeHandle.addEventListener('mousedown', (e) => {
        profilePanelResize.active = true;
        profilePanelResize.startX = e.clientX;
        profilePanelResize.startY = e.clientY;
        profilePanelResize.startW = panel.offsetWidth;
        profilePanelResize.startH = panel.offsetHeight;
        e.preventDefault();
        e.stopPropagation();
    });

    // Mouse move handler
    document.addEventListener('mousemove', (e) => {
        if (profilePanelDrag.active) {
            const dx = e.clientX - profilePanelDrag.startX;
            const dy = e.clientY - profilePanelDrag.startY;
            panel.style.left = (profilePanelDrag.startLeft + dx) + 'px';
            panel.style.top = (profilePanelDrag.startTop + dy) + 'px';
            panel.style.right = 'auto';
        }
        if (profilePanelResize.active) {
            const dx = e.clientX - profilePanelResize.startX;
            const dy = e.clientY - profilePanelResize.startY;
            const newW = Math.max(400, profilePanelResize.startW + dx);
            const newH = Math.max(200, profilePanelResize.startH + dy);
            panel.style.width = newW + 'px';
            panel.style.height = newH + 'px';
            // Update canvas sizes
            resizeProfileCanvases();
        }
    });

    // Mouse up handler
    document.addEventListener('mouseup', () => {
        if (profilePanelDrag.active || profilePanelResize.active) {
            profilePanelDrag.active = false;
            profilePanelResize.active = false;
            // Redraw after resize
            if (COMPARISON.selection.active) {
                updateComparisonProfiles();
            }
        }
    });

    console.log('Profile panel drag/resize initialized');
}

function resizeProfileCanvases() {
    const leftCanvas = document.getElementById('comp-profile-left');
    const rightCanvas = document.getElementById('comp-profile-right');

    if (leftCanvas) {
        const rect = leftCanvas.getBoundingClientRect();
        leftCanvas.width = rect.width * window.devicePixelRatio;
        leftCanvas.height = rect.height * window.devicePixelRatio;
        const ctx = leftCanvas.getContext('2d');
        ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
    }
    if (rightCanvas) {
        const rect = rightCanvas.getBoundingClientRect();
        rightCanvas.width = rect.width * window.devicePixelRatio;
        rightCanvas.height = rect.height * window.devicePixelRatio;
        const ctx = rightCanvas.getContext('2d');
        ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
    }
}

function showComparisonProfilePanel() {
    const panel = document.getElementById('comparison-profile-panel');
    if (panel) {
        panel.style.display = 'block';
        // Initialize canvas sizes
        resizeProfileCanvases();
    }
}

function hideComparisonProfilePanel() {
    const panel = document.getElementById('comparison-profile-panel');
    if (panel) {
        panel.style.display = 'none';
    }
}
