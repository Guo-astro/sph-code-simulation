// ============================================================
// IMBH-Cloud Visualization - Main Entry Point
// Multi-dataset support
// ============================================================

// ============================================================
// Draggable Panel Functionality
// ============================================================
function initDraggablePanels() {
    const draggables = document.querySelectorAll('.draggable');
    
    draggables.forEach(panel => {
        const header = panel.querySelector('.drag-header');
        if (!header) return;
        
        let isDragging = false;
        let startX, startY;
        let initialLeft, initialTop;
        
        header.addEventListener('mousedown', (e) => {
            if (e.button !== 0) return; // Only left click
            
            isDragging = true;
            panel.classList.add('dragging');
            
            // Get current computed position
            const rect = panel.getBoundingClientRect();
            const style = window.getComputedStyle(panel);
            
            // Handle transform (e.g., translateX(-50%))
            const transform = style.transform;
            let translateX = 0, translateY = 0;
            if (transform && transform !== 'none') {
                const matrix = new DOMMatrix(transform);
                translateX = matrix.m41;
                translateY = matrix.m42;
            }
            
            // Convert to left/top positioning
            // Remove any centering transforms by using actual bounding rect
            initialLeft = rect.left;
            initialTop = rect.top;
            
            startX = e.clientX;
            startY = e.clientY;
            
            // Clear positioning properties that might interfere
            panel.style.right = 'auto';
            panel.style.bottom = 'auto';
            panel.style.transform = 'none';
            panel.style.left = initialLeft + 'px';
            panel.style.top = initialTop + 'px';
            
            e.preventDefault();
        });
        
        document.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            
            const dx = e.clientX - startX;
            const dy = e.clientY - startY;
            
            let newLeft = initialLeft + dx;
            let newTop = initialTop + dy;
            
            // Constrain to viewport
            const rect = panel.getBoundingClientRect();
            const maxLeft = window.innerWidth - rect.width;
            const maxTop = window.innerHeight - rect.height;
            
            newLeft = Math.max(0, Math.min(newLeft, maxLeft));
            newTop = Math.max(0, Math.min(newTop, maxTop));
            
            panel.style.left = newLeft + 'px';
            panel.style.top = newTop + 'px';
        });
        
        document.addEventListener('mouseup', () => {
            if (isDragging) {
                isDragging = false;
                panel.classList.remove('dragging');
            }
        });
    });
}

// Load binary format snapshots
async function loadBinarySnapshots(manifest, dsInfo, loadId) {
    const basePath = `data/${dsInfo.path}`;
    const columns = manifest.columns; // ["x", "y", "z", "vx", "vy", "vz", "dens", "temp", "mass", "sound"]
    const snapshots = dsInfo.snapshots;
    const total = snapshots.length;

    for (let i = 0; i < total; i++) {
        if (STATE.loadId !== loadId) return;

        const filename = snapshots[i];
        try {
            const response = await fetch(`${basePath}/${filename}`);
            if (!response.ok) continue;

            const buffer = await response.arrayBuffer();
            const snapshot = parseBinarySnapshot(buffer, columns, i);
            STATE.snapshots.push(snapshot);

            document.getElementById('loading').textContent = `Loading... ${i + 1}/${total} snapshots`;
        } catch (error) {
            console.error(`Error loading ${filename}:`, error);
        }
    }
    console.log(`Loaded ${STATE.snapshots.length} binary snapshots`);
}

// Parse binary snapshot buffer
function parseBinarySnapshot(buffer, columns, frameIndex) {
    const dataView = new DataView(buffer);
    const nParticles = dataView.getUint32(0, true);
    let offset = 4;

    // Read column data
    const columnData = {};
    for (const col of columns) {
        columnData[col] = new Float32Array(buffer, offset, nParticles);
        offset += nParticles * 4;
    }

    // Convert to particle objects
    const particles = [];
    for (let i = 0; i < nParticles; i++) {
        const vx = columnData.vx[i];
        const vy = columnData.vy[i];
        const vz = columnData.vz[i];
        const dens = columnData.dens[i];
        const temp = columnData.temp[i];

        // Estimate pressure from ideal gas: P = (gamma-1) * rho * e, where T = (gamma-1)*mu*mH*e/kB
        // For simplicity, use P ~ rho * T (arbitrary units for entropy calculation)
        const pres = dens * temp;

        particles.push({
            x: columnData.x[i],
            y: columnData.y[i],
            z: columnData.z[i],
            vx: vx, vy: vy, vz: vz,
            vel_mag: Math.sqrt(vx*vx + vy*vy + vz*vz),
            dens: dens,
            temp: temp,
            mass: columnData.mass[i],
            sound: columnData.sound[i],
            pres: pres,
            ene: temp / CONFIG.tempConversion,
            is_ghost: 0
        });
    }

    return { time: frameIndex * 0.05, particles };
}

// Load CSV format snapshots
async function loadCSVSnapshots(loadId) {
    const basePath = STATE.basePath + '/' + (STATE.simType === 'adiabatic'
        ? 'adiabatic/results/'
        : 'cooling/results/');

    let snapshotNum = 1;
    let consecutiveFailures = 0;

    while (consecutiveFailures < 3) {
        if (STATE.loadId !== loadId) return;

        const filename = `snapshot_${String(snapshotNum).padStart(4, '0')}.csv`;
        try {
            const response = await fetch(basePath + filename);
            if (!response.ok) {
                consecutiveFailures++;
                snapshotNum++;
                continue;
            }

            const text = await response.text();
            const data = parseCSV(text);
            if (data.particles.length > 0) {
                STATE.snapshots.push(data);
                consecutiveFailures = 0;
                document.getElementById('loading').textContent = `Loading... ${STATE.snapshots.length} snapshots`;
            }
        } catch (e) {
            consecutiveFailures++;
        }
        snapshotNum++;
        if (snapshotNum > 500) break;
    }
    console.log(`Loaded ${STATE.snapshots.length} CSV snapshots`);
}

async function loadSnapshots() {
    // Cancel any previous load in progress
    STATE.loadId++;
    const currentLoadId = STATE.loadId;

    // Wait if another load is in progress
    if (STATE.isLoading) {
        console.log('Cancelling previous load...');
        // Give previous load time to notice cancellation
        await new Promise(resolve => setTimeout(resolve, 100));
    }

    STATE.isLoading = true;
    document.getElementById('loading').classList.remove('hidden');

    if (!STATE.currentDataset) {
        console.error('No dataset selected');
        document.getElementById('loading').classList.add('hidden');
        STATE.isLoading = false;
        return;
    }

    STATE.snapshots = [];

    // Try binary format first (check for manifest.json)
    let usedBinary = false;
    try {
        const manifestResponse = await fetch('data/manifest.json');
        if (manifestResponse.ok) {
            const manifest = await manifestResponse.json();
            const dsInfo = manifest.datasets.find(d => d.id === STATE.currentDataset.id);
            if (dsInfo && dsInfo.snapshots && dsInfo.snapshots.length > 0) {
                console.log('Binary format detected, loading from manifest');
                await loadBinarySnapshots(manifest, dsInfo, currentLoadId);
                usedBinary = true;
            }
        }
    } catch (e) {
        console.log('No binary manifest, falling back to CSV');
    }

    // Fall back to CSV format
    if (!usedBinary) {
        await loadCSVSnapshots(currentLoadId);
    }

    // Check again if cancelled before processing
    if (currentLoadId !== STATE.loadId) {
        console.log('Load cancelled after fetch');
        return;
    }

    // Sort snapshots by time to ensure correct order
    STATE.snapshots.sort((a, b) => a.time - b.time);

    if (STATE.snapshots.length > 0) {
        // Compute cloud mass, radius, and initial particle count from first snapshot
        let cloudMass = 0;
        let cloudParticles = 0;
        let comX = 0, comY = 0, comZ = 0;

        // First pass: compute mass and COM
        for (const p of STATE.snapshots[0].particles) {
            if (p.is_ghost !== 1) {
                cloudMass += p.mass;
                cloudParticles++;
                comX += p.x * p.mass;
                comY += p.y * p.mass;
                comZ += p.z * p.mass;
            }
        }
        comX /= cloudMass;
        comY /= cloudMass;
        comZ /= cloudMass;

        // Second pass: compute cloud radius (90th percentile distance from COM)
        const distances = [];
        for (const p of STATE.snapshots[0].particles) {
            if (p.is_ghost !== 1) {
                const dx = p.x - comX;
                const dy = p.y - comY;
                const dz = p.z - comZ;
                distances.push(Math.sqrt(dx*dx + dy*dy + dz*dz));
            }
        }
        distances.sort((a, b) => a - b);
        const cloudRadius = distances[Math.floor(distances.length * 0.9)];  // 90th percentile

        STATE.cloudMass = cloudMass;
        STATE.cloudRadius = cloudRadius;
        STATE.initialCloudParticles = cloudParticles;
        console.log(`Cloud mass: ${cloudMass.toFixed(2)} M_sun, R_cloud: ${cloudRadius.toFixed(3)} pc, particles: ${cloudParticles}`);

        // Update penetration factor display
        updatePenetrationFactor();

        // Update Hill circle at pericenter
        if (STATE.hillCirclePeri) {
            STATE.scene.remove(STATE.hillCirclePeri);
            STATE.hillCirclePeri.geometry.dispose();

            const r_peri = CONFIG.r_peri;
            const r_hill_peri = r_peri * Math.pow(cloudMass / (3 * CONFIG.M_BH), 1/3);
            const hillCircleGeom = new THREE.RingGeometry(r_hill_peri * 0.95, r_hill_peri * 1.05, 64);
            const hillCircleMat = new THREE.MeshBasicMaterial({
                color: 0x00ff88,
                side: THREE.DoubleSide,
                transparent: true,
                opacity: 0.5
            });
            STATE.hillCirclePeri = new THREE.Mesh(hillCircleGeom, hillCircleMat);
            STATE.hillCirclePeri.position.set(r_peri, 0, 0);
            STATE.scene.add(STATE.hillCirclePeri);
            console.log(`Hill radius at pericenter: ${r_hill_peri.toFixed(4)} pc`);
        }

        // Compute global color ranges
        computeGlobalColorRanges();

        // Update range sliders to match computed ranges
        updateRangeSliders();

        // Update colorbar
        updateColorbar();

        // Update UI
        document.getElementById('timeline').max = STATE.snapshots.length;
        document.getElementById('frame-display').textContent = `1 / ${STATE.snapshots.length}`;

        // Show first frame
        updateVisualization(0);
    } else {
        console.warn('No snapshots found for', STATE.currentDataset.id, STATE.simType);
        document.getElementById('particle-count').textContent = 'No data';
    }

    STATE.isLoading = false;
    document.getElementById('loading').textContent = 'Loading...';
    document.getElementById('loading').classList.add('hidden');
    console.log(`Loaded ${STATE.snapshots.length} snapshots`);
}

function parseCSV(text) {
    const lines = text.split('\n');
    let time = 0;
    let particleCount = 0;
    const particles = [];

    for (const line of lines) {
        if (line.startsWith('# Time (code):')) {
            time = parseFloat(line.split(':')[1].trim());
        } else if (line.startsWith('# Particle Count:')) {
            particleCount = parseInt(line.split(':')[1].trim());
        } else if (!line.startsWith('#') && line.trim() && !line.startsWith('id,')) {
            const parts = line.split(',');
            if (parts.length >= 30) {
                const is_ghost = parseInt(parts[29]);

                const vx = parseFloat(parts[4]);
                const vy = parseFloat(parts[5]);
                const vz = parseFloat(parts[6]);
                const ene = parseFloat(parts[13]);

                particles.push({
                    id: parseInt(parts[0]),
                    x: parseFloat(parts[1]),
                    y: parseFloat(parts[2]),
                    z: parseFloat(parts[3]),
                    vx: vx,
                    vy: vy,
                    vz: vz,
                    vel_mag: Math.sqrt(vx*vx + vy*vy + vz*vz),
                    mass: parseFloat(parts[10]),
                    dens: parseFloat(parts[11]),
                    pres: parseFloat(parts[12]),      // pressure for entropy
                    ene: ene,
                    sound: parseFloat(parts[15]),     // sound speed for Mach
                    temp: ene * CONFIG.tempConversion,
                    is_ghost: is_ghost
                });
            }
        }
    }

    return { time, particleCount, particles };
}

// Playback controls
function togglePlay() {
    STATE.isPlaying = !STATE.isPlaying;
    const btn = document.getElementById('btn-play');

    if (STATE.isPlaying) {
        btn.textContent = 'Pause';
        btn.classList.add('active');
        STATE.playInterval = setInterval(() => {
            if (STATE.currentFrame < STATE.snapshots.length - 1) {
                updateVisualization(STATE.currentFrame + 1);
            } else {
                updateVisualization(0);
            }
        }, 100);
    } else {
        btn.textContent = 'Play';
        btn.classList.remove('active');
        clearInterval(STATE.playInterval);
    }
}

function prevFrame() {
    if (STATE.currentFrame > 0) {
        updateVisualization(STATE.currentFrame - 1);
    }
}

function nextFrame() {
    if (STATE.currentFrame < STATE.snapshots.length - 1) {
        updateVisualization(STATE.currentFrame + 1);
    }
}

function onTimelineChange(e) {
    updateVisualization(parseInt(e.target.value) - 1);
}

function onSimTypeChange(e) {
    STATE.simType = e.target.value;
    // Reset state but don't recreate orbit (same dataset, different sim type)
    if (STATE.isPlaying) {
        togglePlay();
    }
    STATE.currentFrame = 0;
    STATE.snapshots = [];
    STATE.initialCloudParticles = 0;

    // Clear particle system
    if (STATE.particleSystem) {
        STATE.scene.remove(STATE.particleSystem);
        STATE.particleSystem.geometry.dispose();
        STATE.particleSystem.material.dispose();
        STATE.particleSystem = null;
    }

    // Reset UI
    document.getElementById('timeline').value = 1;
    document.getElementById('frame-display').textContent = '1 / --';
    document.getElementById('accreted-count').textContent = '0';

    loadSnapshots();
}

function onColorModeChange(e) {
    STATE.colorMode = e.target.value;
    // Update sliders to match current mode's range
    updateRangeSliders();
    updateColorbar();
    if (STATE.snapshots.length > 0) {
        updateVisualization(STATE.currentFrame);
    }
}

function onRangeChange() {
    const minVal = parseFloat(document.getElementById('range-min').value);
    const maxVal = parseFloat(document.getElementById('range-max').value);

    // Update display
    document.getElementById('range-min-value').textContent = minVal.toFixed(1);
    document.getElementById('range-max-value').textContent = maxVal.toFixed(1);

    // Ensure min < max
    const actualMin = Math.min(minVal, maxVal);
    const actualMax = Math.max(minVal, maxVal);

    // Update current color mode's range
    STATE.colorRanges[STATE.colorMode].min = actualMin;
    STATE.colorRanges[STATE.colorMode].max = actualMax;

    // Refresh visualization
    updateColorbar();
    if (STATE.snapshots.length > 0) {
        updateVisualization(STATE.currentFrame);
    }
}

function updateRangeSliders() {
    const range = STATE.colorRanges[STATE.colorMode];
    const isLog = range.isLog;

    // Update slider attributes based on log vs linear
    const minSlider = document.getElementById('range-min');
    const maxSlider = document.getElementById('range-max');

    if (isLog) {
        minSlider.min = 0;
        minSlider.max = 8;
        minSlider.step = 0.5;
        maxSlider.min = 0;
        maxSlider.max = 8;
        maxSlider.step = 0.5;
        document.getElementById('range-min-label').textContent = 'Min (log10)';
        document.getElementById('range-max-label').textContent = 'Max (log10)';
    } else {
        // Linear mode (e.g., Mach number)
        minSlider.min = 0;
        minSlider.max = 10;
        minSlider.step = 0.5;
        maxSlider.min = 0;
        maxSlider.max = 10;
        maxSlider.step = 0.5;
        document.getElementById('range-min-label').textContent = 'Min';
        document.getElementById('range-max-label').textContent = 'Max';
    }

    document.getElementById('range-min').value = range.min;
    document.getElementById('range-max').value = range.max;
    document.getElementById('range-min-value').textContent = range.min.toFixed(1);
    document.getElementById('range-max-value').textContent = range.max.toFixed(1);
}

function onDatasetChange(e) {
    const datasetId = e.target.value;
    selectDataset(datasetId);
    resetVisualizationState();
    loadSnapshots();
}

function resetVisualizationState() {
    // Stop playback
    if (STATE.isPlaying) {
        togglePlay();
    }

    // Reset frame counter
    STATE.currentFrame = 0;
    STATE.snapshots = [];
    STATE.initialCloudParticles = 0;
    STATE.cloudMass = 0;

    // Clear old orbit line and pericenter marker
    if (STATE.orbitLine) {
        STATE.scene.remove(STATE.orbitLine);
        STATE.orbitLine.geometry.dispose();
        STATE.orbitLine = null;
    }
    if (STATE.periMarker) {
        STATE.scene.remove(STATE.periMarker);
        STATE.periMarker.geometry.dispose();
        STATE.periMarker.material.dispose();
        STATE.periMarker = null;
    }

    // Clear particle system
    if (STATE.particleSystem) {
        STATE.scene.remove(STATE.particleSystem);
        STATE.particleSystem.geometry.dispose();
        STATE.particleSystem.material.dispose();
        STATE.particleSystem = null;
    }

    // Recreate orbit for new dataset
    createAnalyticOrbit();

    // Reset UI
    document.getElementById('timeline').value = 1;
    document.getElementById('frame-display').textContent = '1 / --';
    document.getElementById('accreted-count').textContent = '0';
}


// Initialization
async function init() {
    // Initialize draggable panels
    initDraggablePanels();
    
    // Load datasets configuration first
    await loadDatasets();

    // Initialize Three.js scene
    initScene();

    // Event listeners
    window.addEventListener('resize', onWindowResize);
    document.getElementById('btn-play').addEventListener('click', togglePlay);
    document.getElementById('btn-prev').addEventListener('click', prevFrame);
    document.getElementById('btn-next').addEventListener('click', nextFrame);
    document.getElementById('timeline').addEventListener('input', onTimelineChange);
    document.getElementById('sim-type').addEventListener('change', onSimTypeChange);
    document.getElementById('color-mode').addEventListener('change', onColorModeChange);
    document.getElementById('range-min').addEventListener('input', onRangeChange);
    document.getElementById('range-max').addEventListener('input', onRangeChange);
    document.getElementById('dataset-select').addEventListener('change', onDatasetChange);

    // Initialize range sliders
    updateRangeSliders();

    // Load data
    loadSnapshots();

    // Start animation loop
    animate();
}

// Start when DOM is ready
document.addEventListener('DOMContentLoaded', init);
