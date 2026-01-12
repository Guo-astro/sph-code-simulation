// ============================================================
// IMBH-Cloud Visualization - Main Entry Point
// ============================================================

async function loadSnapshots() {
    document.getElementById('loading').classList.remove('hidden');

    const basePath = STATE.simType === 'adiabatic'
        ? 'adiabatic_results/'
        : 'cooling_results/';

    STATE.snapshots = [];
    let snapshotNum = 1;
    let consecutiveFailures = 0;

    while (consecutiveFailures < 3) {
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
            }
        } catch (e) {
            consecutiveFailures++;
        }
        snapshotNum++;

        // Safety limit
        if (snapshotNum > 500) break;
    }

    if (STATE.snapshots.length > 0) {
        // Compute cloud mass from first snapshot (sum of non-ghost particles)
        let cloudMass = 0;
        for (const p of STATE.snapshots[0].particles) {
            if (p.is_ghost !== 1) {
                cloudMass += p.mass;
            }
        }
        STATE.cloudMass = cloudMass;
        console.log(`Cloud mass computed from data: ${cloudMass.toFixed(2)} M_sun`);

        // Update Hill circle at pericenter with actual cloud mass
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

        // Compute global color ranges across all snapshots
        computeGlobalColorRanges();

        // Update colorbar with new ranges
        updateColorbar();

        // Update UI
        document.getElementById('timeline').max = STATE.snapshots.length;
        document.getElementById('frame-display').textContent = `1 / ${STATE.snapshots.length}`;

        // Show first frame
        updateVisualization(0);
    }

    document.getElementById('loading').classList.add('hidden');
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
                    ene: ene,
                    // Temperature: T = (gamma-1) * mu * m_H * u / k_B
                    // With gamma=5/3, mu=2.3, u in (km/s)^2: T ≈ 186*u K
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
                updateVisualization(0);  // Loop
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
    loadSnapshots();
}

function onColorModeChange(e) {
    STATE.colorMode = e.target.value;

    // Update colorbar with new mode
    updateColorbar();

    // Refresh visualization with new colors
    if (STATE.snapshots.length > 0) {
        updateVisualization(STATE.currentFrame);
    }
}

// Initialization
function init() {
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

    // Load data
    loadSnapshots();

    // Start animation loop
    animate();
}

// Start when DOM is ready
document.addEventListener('DOMContentLoaded', init);
