// ============================================================
// IMBH-Cloud Visualization - Data Loader
// Supports both CSV and binary formats
// ============================================================

// Detect data format and load appropriately
async function loadSnapshotData() {
    const loadId = ++STATE.loadId;
    STATE.isLoading = true;

    try {
        // Check for manifest.json (binary format)
        const manifestPath = STATE.basePath.replace(/\.\.\//g, '') + '/../data/manifest.json';
        let useBinary = false;
        let manifest = null;

        try {
            const response = await fetch('data/manifest.json');
            if (response.ok) {
                manifest = await response.json();
                // Check if current dataset is in manifest
                const dsInfo = manifest.datasets.find(d => d.id === STATE.currentDataset?.id);
                if (dsInfo && dsInfo.snapshots && dsInfo.snapshots.length > 0) {
                    useBinary = true;
                    console.log('Binary format detected, using manifest');
                }
            }
        } catch (e) {
            console.log('No binary manifest found, using CSV format');
        }

        if (useBinary && manifest) {
            await loadBinarySnapshots(manifest, loadId);
        } else {
            await loadCSVSnapshots(loadId);
        }
    } catch (error) {
        console.error('Error loading data:', error);
    } finally {
        STATE.isLoading = false;
    }
}

// Load snapshots in binary format
async function loadBinarySnapshots(manifest, loadId) {
    const dsInfo = manifest.datasets.find(d => d.id === STATE.currentDataset?.id);
    if (!dsInfo) {
        console.error('Dataset not found in manifest');
        return;
    }

    const basePath = `data/${dsInfo.path}`;
    const snapshots = dsInfo.snapshots;
    const columns = manifest.columns; // ["x", "y", "z", "vx", "vy", "vz", "dens", "temp", "mass", "sound"]

    STATE.snapshots = [];
    const total = snapshots.length;

    for (let i = 0; i < total; i++) {
        if (STATE.loadId !== loadId) {
            console.log('Load cancelled');
            return;
        }

        const filename = snapshots[i];
        const url = `${basePath}/${filename}`;

        try {
            const response = await fetch(url);
            const buffer = await response.arrayBuffer();
            const snapshot = parseBinarySnapshot(buffer, columns, i);

            STATE.snapshots.push(snapshot);
            updateLoadingProgress(i + 1, total);
        } catch (error) {
            console.error(`Error loading ${filename}:`, error);
        }
    }

    console.log(`Loaded ${STATE.snapshots.length} binary snapshots`);
    onDataLoaded();
}

// Parse binary snapshot buffer
function parseBinarySnapshot(buffer, columns, frameIndex) {
    const dataView = new DataView(buffer);

    // Read particle count (uint32, little-endian)
    const nParticles = dataView.getUint32(0, true);
    let offset = 4;

    // Read column data (each column is n_particles × float32)
    const columnData = {};
    for (const col of columns) {
        const floatArray = new Float32Array(buffer, offset, nParticles);
        columnData[col] = floatArray;
        offset += nParticles * 4;
    }

    // Convert to particle objects (matching CSV format)
    const particles = [];
    for (let i = 0; i < nParticles; i++) {
        const vx = columnData.vx[i];
        const vy = columnData.vy[i];
        const vz = columnData.vz[i];

        particles.push({
            id: columnData.id ? columnData.id[i] : i,  // Use file ID or array index
            x: columnData.x[i],
            y: columnData.y[i],
            z: columnData.z[i],
            vx: vx,
            vy: vy,
            vz: vz,
            dens: columnData.dens[i],
            temp: columnData.temp[i],
            mass: columnData.mass[i],
            sound: columnData.sound[i],
            pres: columnData.pres ? columnData.pres[i] : 0,
            vel_mag: Math.sqrt(vx*vx + vy*vy + vz*vz),
            is_ghost: 0
        });
    }

    // Calculate time from frame index (approximate)
    const time = frameIndex * 0.05; // Adjust based on your output frequency

    return { time, particles };
}

// Load snapshots in CSV format (original method)
async function loadCSVSnapshots(loadId) {
    const basePath = STATE.basePath;
    const simType = STATE.simType;

    // Try to get snapshot list from index file or scan directory
    let snapshotFiles = [];

    try {
        const indexUrl = `${basePath}/${simType}/results/index.json`;
        const response = await fetch(indexUrl);
        if (response.ok) {
            const index = await response.json();
            snapshotFiles = index.files;
        }
    } catch (e) {
        // Fall back to numbered snapshots
        for (let i = 1; i <= 100; i++) {
            snapshotFiles.push(`snapshot_${String(i).padStart(4, '0')}.csv`);
        }
    }

    STATE.snapshots = [];

    for (let i = 0; i < snapshotFiles.length; i++) {
        if (STATE.loadId !== loadId) {
            console.log('Load cancelled');
            return;
        }

        const filename = snapshotFiles[i];
        const url = `${basePath}/${simType}/results/${filename}`;

        try {
            const response = await fetch(url);
            if (!response.ok) {
                if (i > 0) break; // Stop at first missing file (after at least one loaded)
                continue;
            }

            const csvText = await response.text();
            const snapshot = parseCSVSnapshot(csvText, i);

            STATE.snapshots.push(snapshot);
            updateLoadingProgress(i + 1, snapshotFiles.length);
        } catch (error) {
            if (i > 0) break;
        }
    }

    console.log(`Loaded ${STATE.snapshots.length} CSV snapshots`);
    onDataLoaded();
}

// Parse CSV snapshot (existing logic)
function parseCSVSnapshot(csvText, frameIndex) {
    const lines = csvText.trim().split('\n');
    const headers = lines[0].split(',').map(h => h.trim());

    const particles = [];
    for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => parseFloat(v.trim()));
        const particle = {};

        headers.forEach((h, idx) => {
            particle[h] = values[idx];
        });

        // Ensure particle has an ID for Lagrangian tracking
        if (particle.id === undefined) {
            particle.id = i - 1;  // Use array index as ID
        }

        // Add derived values
        if (particle.vx !== undefined) {
            particle.vel_mag = Math.sqrt(
                particle.vx**2 + particle.vy**2 + particle.vz**2
            );
        }

        particles.push(particle);
    }

    // Get time from first particle if available
    const time = particles[0]?.time ?? frameIndex * 0.05;

    return { time, particles };
}

// Update loading progress display
function updateLoadingProgress(current, total) {
    const percent = Math.round((current / total) * 100);
    const loadingEl = document.getElementById('loading');
    if (loadingEl) {
        const text = loadingEl.querySelector('span') || loadingEl;
        text.textContent = `Loading snapshots... ${current}/${total} (${percent}%)`;
    }
}

// Called after data is fully loaded
function onDataLoaded() {
    // Extract cloud properties from first snapshot
    if (STATE.snapshots.length > 0) {
        const firstSnap = STATE.snapshots[0];
        const cloudParticles = firstSnap.particles.filter(p => p.is_ghost !== 1);

        STATE.cloudMass = cloudParticles.reduce((sum, p) => sum + (p.mass || 0), 0);
        STATE.initialCloudParticles = cloudParticles.length;

        // Estimate cloud radius (half-mass radius approximation)
        const cx = cloudParticles.reduce((s, p) => s + p.x, 0) / cloudParticles.length;
        const cy = cloudParticles.reduce((s, p) => s + p.y, 0) / cloudParticles.length;
        const cz = cloudParticles.reduce((s, p) => s + p.z, 0) / cloudParticles.length;

        const distances = cloudParticles.map(p =>
            Math.sqrt((p.x-cx)**2 + (p.y-cy)**2 + (p.z-cz)**2)
        ).sort((a, b) => a - b);

        STATE.cloudRadius = distances[Math.floor(distances.length * 0.9)] || 1.0;
    }

    // Compute global color ranges (with error handling)
    try {
        computeGlobalColorRanges();
    } catch (e) {
        console.warn('Error computing color ranges:', e);
    }

    // Compute global profile ranges for fixed plot scales (with error handling)
    try {
        if (typeof computeGlobalProfileRanges === 'function') {
            computeGlobalProfileRanges();
        }
    } catch (e) {
        console.warn('Error computing profile ranges:', e);
    }

    // Update UI
    updatePenetrationFactor();

    // Initialize display
    if (STATE.snapshots.length > 0) {
        STATE.currentFrame = 0;
        updateParticles(STATE.snapshots[0]);
        updatePVDiagram(STATE.snapshots[0]);
        updateInfoPanel(STATE.snapshots[0]);

        // Update timeline
        const timeline = document.getElementById('timeline');
        if (timeline) {
            timeline.max = STATE.snapshots.length - 1;
            timeline.value = 0;
        }
    }

    // Hide loading overlay
    const loadingEl = document.getElementById('loading');
    if (loadingEl) {
        loadingEl.classList.add('hidden');
    }
}
