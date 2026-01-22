// ============================================================
// IMBH-Cloud Visualization - Data Loader
// Supports HDF5 format via h5wasm
// ============================================================

// Global h5wasm instance
let h5wasmModule = null;

// Initialize h5wasm
async function initH5wasm() {
    if (h5wasmModule) return h5wasmModule;

    try {
        // Load h5wasm from CDN
        const { h5wasm } = await import('https://cdn.jsdelivr.net/npm/h5wasm@0.6.8/+esm');
        h5wasmModule = await h5wasm.ready;
        console.log('h5wasm initialized');
        return h5wasmModule;
    } catch (error) {
        console.error('Failed to load h5wasm:', error);
        return null;
    }
}

// Detect data format and load appropriately
async function loadSnapshotData() {
    const loadId = ++STATE.loadId;
    STATE.isLoading = true;

    try {
        // Initialize h5wasm
        await initH5wasm();

        // Load HDF5 snapshots
        await loadHDF5Snapshots(loadId);
    } catch (error) {
        console.error('Error loading data:', error);
    } finally {
        STATE.isLoading = false;
    }
}

// Load snapshots in HDF5 format
async function loadHDF5Snapshots(loadId) {
    const basePath = STATE.basePath;
    const simType = STATE.simType;

    STATE.snapshots = [];
    let fileIndex = 1;
    let consecutiveFails = 0;

    while (consecutiveFails < 3) {
        if (STATE.loadId !== loadId) {
            console.log('Load cancelled');
            return;
        }

        const filename = `snapshot_${String(fileIndex).padStart(4, '0')}.h5`;
        const url = `${basePath}/${simType}/results/${filename}`;

        try {
            const response = await fetch(url);
            if (!response.ok) {
                consecutiveFails++;
                fileIndex++;
                continue;
            }

            const arrayBuffer = await response.arrayBuffer();
            const snapshot = await parseHDF5Snapshot(arrayBuffer, fileIndex - 1);

            if (snapshot) {
                STATE.snapshots.push(snapshot);
                updateLoadingProgress(STATE.snapshots.length, fileIndex);
            }

            consecutiveFails = 0;
        } catch (error) {
            consecutiveFails++;
        }

        fileIndex++;
        if (fileIndex > 500) break;  // Safety limit
    }

    console.log(`Loaded ${STATE.snapshots.length} HDF5 snapshots`);
    onDataLoaded();
}

// Parse HDF5 snapshot
async function parseHDF5Snapshot(arrayBuffer, frameIndex) {
    if (!h5wasmModule) {
        await initH5wasm();
    }

    if (!h5wasmModule) {
        console.error('h5wasm not available');
        return null;
    }

    try {
        // Create a file from the array buffer
        const filename = `temp_${Date.now()}.h5`;
        h5wasmModule.FS.writeFile(filename, new Uint8Array(arrayBuffer));

        const file = new h5wasmModule.File(filename, 'r');

        // Read time from header
        let time = frameIndex * 0.02;
        try {
            if (file.get('header')) {
                const header = file.get('header');
                if (header.attrs && header.attrs['time']) {
                    time = header.attrs['time'].value;
                }
            }
        } catch (e) {
            console.log('No time attribute found, using default');
        }

        // Read particle data
        const particles = [];
        const particlesGroup = file.get('particles');

        if (!particlesGroup) {
            console.error('No particles group in HDF5 file');
            file.close();
            h5wasmModule.FS.unlink(filename);
            return null;
        }

        // Read datasets
        const pos_x = particlesGroup.get('pos_x')?.value || [];
        const pos_y = particlesGroup.get('pos_y')?.value || [];
        const pos_z = particlesGroup.get('pos_z')?.value || [];
        const vel_x = particlesGroup.get('vel_x')?.value || [];
        const vel_y = particlesGroup.get('vel_y')?.value || [];
        const vel_z = particlesGroup.get('vel_z')?.value || [];
        const dens = particlesGroup.get('dens')?.value || [];
        const pres = particlesGroup.get('pres')?.value || [];
        const ene = particlesGroup.get('ene')?.value || [];
        const mass = particlesGroup.get('mass')?.value || [];
        const sml = particlesGroup.get('sml')?.value || [];

        const nParticles = pos_x.length;

        for (let i = 0; i < nParticles; i++) {
            const vx = vel_x[i] || 0;
            const vy = vel_y[i] || 0;
            const vz = vel_z[i] || 0;

            particles.push({
                id: i,
                x: pos_x[i],
                y: pos_y[i],
                z: pos_z[i],
                vx: vx,
                vy: vy,
                vz: vz,
                dens: dens[i] || 0,
                pres: pres[i] || 0,
                temp: (ene[i] || 0) * (CONFIG?.tempConversion || 186),
                mass: mass[i] || 0,
                sound: sml[i] || 0,
                vel_mag: Math.sqrt(vx*vx + vy*vy + vz*vz),
                is_ghost: 0
            });
        }

        file.close();
        h5wasmModule.FS.unlink(filename);

        console.log(`Parsed HDF5 frame ${frameIndex}: ${nParticles} particles, time=${time}`);
        return { time, particles };

    } catch (error) {
        console.error('Error parsing HDF5:', error);
        return null;
    }
}

// Update loading progress display
function updateLoadingProgress(current, total) {
    const percent = Math.round((current / total) * 100);
    const loadingEl = document.getElementById('loading');
    if (loadingEl) {
        const text = loadingEl.querySelector('span') || loadingEl;
        text.textContent = `Loading HDF5 snapshots... ${current} loaded`;
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
