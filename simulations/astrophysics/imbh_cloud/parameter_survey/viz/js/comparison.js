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
    loadProgress: { left: 0, right: 0 }
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

    // Show comparison UI
    showComparisonUI();

    // Initialize dual renderers
    initDualRenderers();

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

        // Color by density (same as main viz)
        const colorValue = getComparisonColorValue(p, comVel);
        const color = valueToColor(colorValue);

        colors[i * 3] = color.r;
        colors[i * 3 + 1] = color.g;
        colors[i * 3 + 2] = color.b;
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
    const colorMode = STATE?.colorMode || 'density';
    const range = STATE?.colorRanges?.[colorMode] || { min: 3, max: 6 };

    let value;
    switch (colorMode) {
        case 'temperature':
            value = Math.log10(Math.max(particle.temp, 1));
            break;
        case 'velocity':
            const dvx = particle.vx - comVel.x;
            const dvy = particle.vy - comVel.y;
            const dvz = particle.vz - comVel.z;
            value = Math.log10(Math.max(Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz), 0.1));
            break;
        case 'density':
        default:
            const n_H2 = particle.dens * (CONFIG?.densityToNH2 || 20.3);
            value = Math.log10(Math.max(n_H2, 1));
            break;
    }

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

    // Restore main UI
    hideComparisonUI();

    console.log('Exited comparison mode');
}
