// ============================================================
// IMBH-Cloud Visualization - Column & Plane Selection with Profiles
// Drag-to-select with fixed plot ranges across all frames
// ============================================================

const TRACKING = {
    // Current display mode for profiles
    displayMode: 'column',  // 'column' or 'plane'

    // Selection mode (must be enabled to drag-select)
    selectionModeActive: false,

    // Column selection (cylinder along z)
    column: {
        active: false,
        center: null,           // {x, y}
        radius: 0.3,            // pc
        selectedIds: new Set(),
        mesh: null
    },

    // Plane selection (horizontal slice)
    plane: {
        active: false,
        zMin: 0,
        zMax: 0,
        zPosition: 0,           // Current z position of draggable plane
        selectedIds: new Set(),
        mesh: null,
        edges: null,
        // Draggable plane state
        draggableMesh: null,
        draggableEdges: null,
        isDragging: false,
        dragStartY: 0,
        dragStartZ: 0
    },

    // For plane mode: which variable to color 2D scatter
    colorVariable: 'density',

    // Computed profiles
    profiles: null,

    // GLOBAL RANGES for fixed plot scales (computed once from all data)
    globalRanges: null
};

// Physical constants
const PHYS = {
    G: 0.00430091,
    M_BH: 100000.0,
    k_B: 1.381e-16,
    m_H: 1.673e-24,
    pc_to_cm: 3.086e18,
    Msun_to_g: 1.989e33,
    Myr_to_s: 3.156e13
};

// ============================================================
// Compute Global Ranges (optimized - samples frames to avoid memory issues)
// ============================================================

function computeGlobalProfileRanges() {
    if (!STATE.snapshots || STATE.snapshots.length === 0) return;

    const numFrames = STATE.snapshots.length;

    // Sample frames: first, last, and evenly spaced in between (max 10 frames)
    const maxSamples = Math.min(10, numFrames);
    const sampleIndices = new Set();
    sampleIndices.add(0);                    // First frame
    sampleIndices.add(numFrames - 1);        // Last frame

    // Add evenly spaced frames
    for (let i = 1; i < maxSamples - 1; i++) {
        sampleIndices.add(Math.floor(i * numFrames / (maxSamples - 1)));
    }

    const framesToScan = Array.from(sampleIndices).sort((a, b) => a - b);
    console.log(`Computing global ranges from ${framesToScan.length} sampled frames...`);

    const ranges = {
        z: { min: Infinity, max: -Infinity },
        x: { min: Infinity, max: -Infinity },
        y: { min: Infinity, max: -Infinity },
        density: { min: Infinity, max: -Infinity },
        temperature: { min: Infinity, max: -Infinity },
        vz: { min: Infinity, max: -Infinity },
        mach: { min: 0, max: -Infinity },
        a_BH: { min: Infinity, max: -Infinity },
        coolingEff: { min: Infinity, max: -Infinity },
        entropy: { min: Infinity, max: -Infinity }
    };

    // Process each sampled frame
    for (const frameIdx of framesToScan) {
        const snapshot = STATE.snapshots[frameIdx];
        if (!snapshot || !snapshot.particles) continue;

        const particles = snapshot.particles;
        const nParticles = particles.length;

        // Sample particles within frame (max 2000 particles per frame)
        const particleSampleRate = Math.max(1, Math.floor(nParticles / 2000));

        // Quick COM estimate from sampled particles
        let comVx = 0, comVy = 0, comVz = 0, totalMass = 0;
        for (let i = 0; i < nParticles; i += particleSampleRate) {
            const p = particles[i];
            if (p.is_ghost === 1) continue;
            comVx += p.vx * p.mass;
            comVy += p.vy * p.mass;
            comVz += p.vz * p.mass;
            totalMass += p.mass;
        }
        if (totalMass > 0) {
            comVx /= totalMass;
            comVy /= totalMass;
            comVz /= totalMass;
        }

        // Scan sampled particles
        for (let i = 0; i < nParticles; i += particleSampleRate) {
            const p = particles[i];
            if (p.is_ghost === 1) continue;

            // Position ranges
            if (p.z < ranges.z.min) ranges.z.min = p.z;
            if (p.z > ranges.z.max) ranges.z.max = p.z;
            if (p.x < ranges.x.min) ranges.x.min = p.x;
            if (p.x > ranges.x.max) ranges.x.max = p.x;
            if (p.y < ranges.y.min) ranges.y.min = p.y;
            if (p.y > ranges.y.max) ranges.y.max = p.y;

            // Density (log)
            const n_H2 = p.dens * CONFIG.densityToNH2;
            if (n_H2 > 0) {
                const logDens = Math.log10(n_H2);
                if (logDens < ranges.density.min) ranges.density.min = logDens;
                if (logDens > ranges.density.max) ranges.density.max = logDens;
            }

            // Temperature (log)
            if (p.temp > 0) {
                const logTemp = Math.log10(p.temp);
                if (logTemp < ranges.temperature.min) ranges.temperature.min = logTemp;
                if (logTemp > ranges.temperature.max) ranges.temperature.max = logTemp;
            }

            // v_z (linear)
            if (p.vz < ranges.vz.min) ranges.vz.min = p.vz;
            if (p.vz > ranges.vz.max) ranges.vz.max = p.vz;

            // Mach number
            const dvx = p.vx - comVx;
            const dvy = p.vy - comVy;
            const dvz = p.vz - comVz;
            const v_rel = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            const mach = p.sound > 0 ? v_rel / p.sound : 0;
            if (mach > ranges.mach.max) ranges.mach.max = mach;

            // IMBH acceleration (log)
            const r2 = p.x*p.x + p.y*p.y + p.z*p.z;
            if (r2 > 0) {
                const a_BH = PHYS.G * PHYS.M_BH / r2;
                const a_BH_cgs = a_BH * PHYS.pc_to_cm / (PHYS.Myr_to_s * PHYS.Myr_to_s);
                if (a_BH_cgs > 0) {
                    const logA = Math.log10(a_BH_cgs);
                    if (logA < ranges.a_BH.min) ranges.a_BH.min = logA;
                    if (logA > ranges.a_BH.max) ranges.a_BH.max = logA;
                }
            }

            // Cooling efficiency (log)
            if (p.temp > 0 && n_H2 > 0) {
                const Lambda = 1e-27 * Math.pow(p.temp / 10, 2);
                const t_cool = (1.5 * PHYS.k_B * p.temp) / (n_H2 * Lambda);
                const rho_cgs = p.dens * PHYS.Msun_to_g / Math.pow(PHYS.pc_to_cm, 3);
                const t_dyn = 1.0 / Math.sqrt(6.674e-8 * rho_cgs);
                const eta = t_dyn / t_cool;
                if (eta > 0) {
                    const logEta = Math.log10(eta);
                    if (logEta < ranges.coolingEff.min) ranges.coolingEff.min = logEta;
                    if (logEta > ranges.coolingEff.max) ranges.coolingEff.max = logEta;
                }
            }

            // Entropy s = P / rho^gamma (log scale)
            if (p.pres !== undefined && p.dens > 0) {
                const gamma = 5 / 3;
                const entropy = p.pres / Math.pow(p.dens, gamma);
                if (entropy > 0) {
                    const logS = Math.log10(entropy);
                    if (logS < ranges.entropy.min) ranges.entropy.min = logS;
                    if (logS > ranges.entropy.max) ranges.entropy.max = logS;
                }
            }
        }
    }

    // Add padding to ranges (10% to be safe with sampling)
    const pad = (r, fraction = 0.1) => {
        if (!isFinite(r.min) || !isFinite(r.max)) return;
        const range = r.max - r.min;
        if (range > 0) {
            r.min -= range * fraction;
            r.max += range * fraction;
        }
    };

    pad(ranges.z);
    pad(ranges.x);
    pad(ranges.y);
    pad(ranges.density);
    pad(ranges.temperature);
    pad(ranges.vz);
    if (ranges.mach.max > 0) ranges.mach.max *= 1.2;
    pad(ranges.a_BH);
    pad(ranges.coolingEff);
    pad(ranges.entropy);

    TRACKING.globalRanges = ranges;
    console.log('Global ranges computed:', ranges);
}

// ============================================================
// Display Mode Switching
// ============================================================

function setTrackingMode(mode) {
    TRACKING.displayMode = mode;

    const columnBtn = document.getElementById('btn-column-mode');
    const planeBtn = document.getElementById('btn-plane-mode');

    if (columnBtn && planeBtn) {
        columnBtn.classList.toggle('active', mode === 'column');
        planeBtn.classList.toggle('active', mode === 'plane');
    }

    const columnControls = document.getElementById('column-controls');
    const planeControls = document.getElementById('plane-controls');

    if (columnControls) columnControls.style.display = mode === 'column' ? 'block' : 'none';
    if (planeControls) planeControls.style.display = mode === 'plane' ? 'block' : 'none';

    // Update instruction text
    updateInstructionText();

    // Handle draggable plane for plane mode
    if (TRACKING.selectionModeActive) {
        if (mode === 'plane') {
            createDraggablePlane();
            // Disable orbit controls for plane dragging
            if (STATE.controls) {
                STATE.controls.enabled = false;
                console.log('OrbitControls disabled for plane mode');
            }
        } else {
            removeDraggablePlanePreview();
            // Re-enable orbit controls for column mode (click-based)
            if (STATE.controls) {
                STATE.controls.enabled = true;
                console.log('OrbitControls re-enabled for column mode');
            }
        }
    }

    computeProfiles();
    updateProfileCanvas();
}

// ============================================================
// Selection Mode Toggle
// ============================================================

function toggleSelectionMode() {
    TRACKING.selectionModeActive = !TRACKING.selectionModeActive;

    const btn = document.getElementById('btn-select-mode');
    const instruction = document.getElementById('drag-instruction');

    if (TRACKING.selectionModeActive) {
        btn.textContent = 'Exit Selection Mode';
        btn.style.background = '#664433';
        btn.style.borderColor = '#886655';
        btn.style.color = '#ffddaa';

        instruction.style.display = 'block';
        updateInstructionText();

        // For plane mode, create draggable plane and disable orbit controls
        if (TRACKING.displayMode === 'plane') {
            createDraggablePlane();
            // Disable orbit controls in plane mode so we can drag the plane
            if (STATE.controls) {
                STATE.controls.enabled = false;
                console.log('OrbitControls disabled for plane dragging');
            }
        }
    } else {
        btn.textContent = 'Enter Selection Mode';
        btn.style.background = '#334455';
        btn.style.borderColor = '#556677';
        btn.style.color = '#aabbcc';

        instruction.style.display = 'none';
        document.body.style.cursor = 'default';

        // Remove draggable plane preview
        removeDraggablePlanePreview();

        // Re-enable orbit controls
        if (STATE.controls) {
            STATE.controls.enabled = true;
            console.log('OrbitControls re-enabled');
        }
    }

    console.log('Selection mode:', TRACKING.selectionModeActive ? 'ON' : 'OFF');
}

function updateInstructionText() {
    const instruction = document.getElementById('drag-instruction');
    if (!instruction) return;

    if (TRACKING.displayMode === 'column') {
        instruction.textContent = 'Click on cloud to place column center';
    } else {
        instruction.textContent = 'Click and drag up/down to move the plane';
    }
}

// ============================================================
// Mouse Handlers
// ============================================================

function initSelectionHandlers() {
    const container = document.getElementById('container');
    if (!container) return;

    // Use 'click' for column selection (not intercepted by OrbitControls)
    window.addEventListener('click', onGlobalClick, false);

    // Use mousedown/move/up for plane dragging
    container.addEventListener('mousedown', onMouseDown, false);
    window.addEventListener('mousemove', onMouseMove, false);
    window.addEventListener('mouseup', onMouseUp, false);

    // ESC key to exit selection mode
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && TRACKING.selectionModeActive) {
            toggleSelectionMode();
        }
    });

    console.log('Selection handlers initialized');
}

function onGlobalClick(event) {
    console.log('onGlobalClick:', event.target.tagName, 'selectionMode:', TRACKING.selectionModeActive, 'displayMode:', TRACKING.displayMode);

    // Only handle if in column selection mode
    if (!TRACKING.selectionModeActive) return;
    if (TRACKING.displayMode !== 'column') return;

    // Check if click is on the 3D canvas area
    const container = document.getElementById('container');
    if (!container) return;

    // Skip if clicking on UI elements
    if (event.target.closest('.draggable, #controls, #timeline-container, button, input, select')) {
        console.log('Click on UI element, ignoring');
        return;
    }

    // Check if click is within the container bounds
    const rect = container.getBoundingClientRect();
    if (event.clientX < rect.left || event.clientX > rect.right ||
        event.clientY < rect.top || event.clientY > rect.bottom) {
        console.log('Click outside container');
        return;
    }

    console.log('Processing column selection click');
    handleColumnClick(event);
}

function onMouseDown(event) {
    console.log('onMouseDown:', {
        selectionMode: TRACKING.selectionModeActive,
        displayMode: TRACKING.displayMode,
        button: event.button,
        target: event.target.tagName
    });

    if (!TRACKING.selectionModeActive) {
        console.log('Selection mode not active, ignoring');
        return;
    }
    if (event.button !== 0) {
        console.log('Not left button, ignoring');
        return;
    }
    if (event.target.closest('.draggable, #controls, #timeline-container')) {
        console.log('Click on UI element, ignoring');
        return;
    }

    // Plane mode: check if clicking on the draggable plane
    if (TRACKING.displayMode === 'plane') {
        console.log('Calling handlePlaneDragStart');
        handlePlaneDragStart(event);
    }
}

function onMouseMove(event) {
    if (!TRACKING.selectionModeActive) return;

    if (TRACKING.displayMode === 'plane' && TRACKING.plane.isDragging) {
        handlePlaneDrag(event);
    }
}

function onMouseUp(event) {
    if (TRACKING.displayMode === 'plane' && TRACKING.plane.isDragging) {
        handlePlaneDragEnd(event);
    }
}

// ============================================================
// Column Selection - Click to Place
// ============================================================

function handleColumnClick(event) {
    console.log('handleColumnClick called');

    const point = raycastToPoint(event);
    if (!point) {
        console.log('No point from raycast');
        return;
    }

    console.log('Raycast hit:', point);

    TRACKING.column.center = { x: point.x, y: point.y };

    // Use radius from slider
    const radiusSlider = document.getElementById('column-radius');
    if (radiusSlider) {
        TRACKING.column.radius = parseFloat(radiusSlider.value);
    }

    selectParticlesInColumn();
    console.log(`Column placed at (${point.x.toFixed(2)}, ${point.y.toFixed(2)}), r=${TRACKING.column.radius.toFixed(2)}`);
}

function raycastToPoint(event) {
    const container = document.getElementById('container');
    const rect = container.getBoundingClientRect();
    const mouse = new THREE.Vector2(
        ((event.clientX - rect.left) / rect.width) * 2 - 1,
        -((event.clientY - rect.top) / rect.height) * 2 + 1
    );

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(mouse, STATE.camera);

    // Try to hit particles first (use larger threshold for easier selection)
    if (STATE.particleSystem) {
        raycaster.params.Points.threshold = 0.2;
        const intersects = raycaster.intersectObject(STATE.particleSystem);
        if (intersects.length > 0) {
            console.log('Raycast hit particle at:', intersects[0].point);
            return intersects[0].point;
        }
    }

    // Fallback: project ray onto a plane through cloud center
    // Use z=0 or compute cloud center z
    let centerZ = 0;
    if (STATE.snapshots && STATE.snapshots[STATE.currentFrame]) {
        const particles = STATE.snapshots[STATE.currentFrame].particles;
        let sumZ = 0, count = 0;
        for (const p of particles) {
            if (p.is_ghost !== 1) {
                sumZ += p.z;
                count++;
            }
        }
        if (count > 0) centerZ = sumZ / count;
    }

    const plane = new THREE.Plane(new THREE.Vector3(0, 0, 1), -centerZ);
    const intersection = new THREE.Vector3();
    if (raycaster.ray.intersectPlane(plane, intersection)) {
        console.log('Raycast fallback to plane at z=' + centerZ + ':', intersection);
        return intersection;
    }

    console.log('Raycast failed completely');
    return null;
}

// ============================================================
// Plane Selection - Draggable Plane
// ============================================================

function createDraggablePlane() {
    removeDraggablePlanePreview();

    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) return;

    // Get cloud extent (local variables, don't overwrite selection range)
    let cloudXMin = Infinity, cloudXMax = -Infinity;
    let cloudYMin = Infinity, cloudYMax = -Infinity;
    let cloudZMin = Infinity, cloudZMax = -Infinity;

    for (const p of data.particles) {
        if (p.is_ghost === 1) continue;
        cloudXMin = Math.min(cloudXMin, p.x);
        cloudXMax = Math.max(cloudXMax, p.x);
        cloudYMin = Math.min(cloudYMin, p.y);
        cloudYMax = Math.max(cloudYMax, p.y);
        cloudZMin = Math.min(cloudZMin, p.z);
        cloudZMax = Math.max(cloudZMax, p.z);
    }

    const width = (cloudXMax - cloudXMin) * 1.5;
    const height = (cloudYMax - cloudYMin) * 1.5;
    const cx = (cloudXMax + cloudXMin) / 2;
    const cy = (cloudYMax + cloudYMin) / 2;

    // Start at middle z (only set zPosition, not selection range)
    const zMid = (cloudZMax + cloudZMin) / 2;
    TRACKING.plane.zPosition = zMid;
    // Store cloud extent for reference (not selection range)
    TRACKING.plane.cloudZMin = cloudZMin;
    TRACKING.plane.cloudZMax = cloudZMax;

    // Create draggable plane mesh
    const geometry = new THREE.PlaneGeometry(width, height);
    const material = new THREE.MeshBasicMaterial({
        color: 0x00ff88,
        transparent: true,
        opacity: 0.25,
        side: THREE.DoubleSide,
        depthWrite: false
    });

    TRACKING.plane.draggableMesh = new THREE.Mesh(geometry, material);
    TRACKING.plane.draggableMesh.position.set(cx, cy, zMid);
    TRACKING.plane.draggableMesh.userData.isDraggablePlane = true;

    // Add edge outline
    const edgeGeom = new THREE.EdgesGeometry(geometry);
    const edgeMat = new THREE.LineBasicMaterial({ color: 0x00ff88, linewidth: 2 });
    TRACKING.plane.draggableEdges = new THREE.LineSegments(edgeGeom, edgeMat);
    TRACKING.plane.draggableEdges.position.copy(TRACKING.plane.draggableMesh.position);

    const group = STATE.orbitalPlaneGroup || STATE.scene;
    group.add(TRACKING.plane.draggableMesh);
    group.add(TRACKING.plane.draggableEdges);

    // Update UI
    updatePlaneZDisplay();

    console.log('Draggable plane created at z =', zMid.toFixed(2));
}

function removeDraggablePlanePreview() {
    const group = STATE.orbitalPlaneGroup || STATE.scene;

    if (TRACKING.plane.draggableMesh) {
        group.remove(TRACKING.plane.draggableMesh);
        TRACKING.plane.draggableMesh.geometry.dispose();
        TRACKING.plane.draggableMesh.material.dispose();
        TRACKING.plane.draggableMesh = null;
    }
    if (TRACKING.plane.draggableEdges) {
        group.remove(TRACKING.plane.draggableEdges);
        TRACKING.plane.draggableEdges.geometry.dispose();
        TRACKING.plane.draggableEdges.material.dispose();
        TRACKING.plane.draggableEdges = null;
    }
}

function handlePlaneDragStart(event) {
    console.log('handlePlaneDragStart called');
    console.log('draggableMesh exists:', !!TRACKING.plane.draggableMesh);
    console.log('zPosition:', TRACKING.plane.zPosition);

    if (!TRACKING.plane.draggableMesh) {
        console.log('No draggable mesh - aborting');
        return;
    }

    // In plane selection mode, start dragging from anywhere in the canvas
    // (no need to click exactly on the plane mesh)
    TRACKING.plane.isDragging = true;
    TRACKING.plane.dragStartY = event.clientY;
    TRACKING.plane.dragStartZ = TRACKING.plane.zPosition;

    document.body.style.cursor = 'ns-resize';
    console.log('Plane drag started at z =', TRACKING.plane.zPosition);
}

function handlePlaneDrag(event) {
    if (!TRACKING.plane.isDragging) return;

    // Calculate z change based on mouse Y movement
    const deltaY = event.clientY - TRACKING.plane.dragStartY;

    // Convert screen pixels to world units (approximate)
    const sensitivity = 0.01;
    const deltaZ = -deltaY * sensitivity;

    // Update z position
    TRACKING.plane.zPosition = TRACKING.plane.dragStartZ + deltaZ;

    // Update mesh position
    if (TRACKING.plane.draggableMesh) {
        TRACKING.plane.draggableMesh.position.z = TRACKING.plane.zPosition;
    }
    if (TRACKING.plane.draggableEdges) {
        TRACKING.plane.draggableEdges.position.z = TRACKING.plane.zPosition;
    }

    updatePlaneZDisplay();
}

function handlePlaneDragEnd(event) {
    if (!TRACKING.plane.isDragging) return;

    TRACKING.plane.isDragging = false;
    document.body.style.cursor = 'default';

    // Don't re-enable orbit controls here - keep them disabled until exiting selection mode
    // if (STATE.controls) STATE.controls.enabled = true;

    // Select particles at this z with current thickness
    selectParticlesAtPlanePosition();

    console.log('Plane drag ended at z =', TRACKING.plane.zPosition.toFixed(2));
}

function selectParticlesAtPlanePosition() {
    const thickness = parseFloat(document.getElementById('plane-thickness')?.value || 0.2);
    const zCenter = TRACKING.plane.zPosition;

    console.log('selectParticlesAtPlanePosition: zCenter=' + zCenter + ', thickness=' + thickness);

    TRACKING.plane.zMin = zCenter - thickness / 2;
    TRACKING.plane.zMax = zCenter + thickness / 2;

    console.log('Selection range: z=[' + TRACKING.plane.zMin.toFixed(3) + ', ' + TRACKING.plane.zMax.toFixed(3) + ']');

    selectParticlesInPlane();
}

function updatePlaneZDisplay() {
    const zDisplay = document.getElementById('plane-z-value');
    if (zDisplay && TRACKING.plane.zPosition !== undefined) {
        zDisplay.textContent = TRACKING.plane.zPosition.toFixed(2) + ' pc';
    }
}

function updatePlaneThickness(value) {
    const thickness = parseFloat(value);
    const thicknessDisplay = document.getElementById('plane-thickness-value');
    if (thicknessDisplay) {
        thicknessDisplay.textContent = thickness.toFixed(2) + ' pc';
    }

    // Re-select if plane is active
    if (TRACKING.plane.active && TRACKING.plane.zPosition !== undefined) {
        selectParticlesAtPlanePosition();
    }
}

// ============================================================
// Particle Selection
// ============================================================

function selectParticlesInColumn() {
    console.log('selectParticlesInColumn called');
    console.log('Column center:', TRACKING.column.center);
    console.log('Column radius:', TRACKING.column.radius);

    if (!TRACKING.column.center) {
        console.log('No column center set');
        return;
    }

    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) {
        console.log('No snapshot data');
        return;
    }

    const cx = TRACKING.column.center.x;
    const cy = TRACKING.column.center.y;
    const r2 = TRACKING.column.radius ** 2;

    console.log(`Selecting particles: center=(${cx.toFixed(3)}, ${cy.toFixed(3)}), r²=${r2.toFixed(3)}`);
    console.log(`Total particles: ${data.particles.length}`);

    // Check first particle to verify data structure
    if (data.particles.length > 0) {
        const p0 = data.particles[0];
        console.log('First particle:', { id: p0.id, x: p0.x, y: p0.y, z: p0.z });
    }

    TRACKING.column.selectedIds.clear();

    data.particles.forEach(p => {
        if (p.is_ghost === 1) return;
        const dx = p.x - cx;
        const dy = p.y - cy;
        if (dx*dx + dy*dy <= r2) {
            TRACKING.column.selectedIds.add(p.id);
        }
    });

    TRACKING.column.active = TRACKING.column.selectedIds.size > 0;
    console.log(`Column: ${TRACKING.column.selectedIds.size} particles selected`);

    // Update UI
    const radiusEl = document.getElementById('column-radius-value');
    if (radiusEl) radiusEl.textContent = TRACKING.column.radius.toFixed(2) + ' pc';
    const radiusSlider = document.getElementById('column-radius');
    if (radiusSlider) radiusSlider.value = Math.min(1, TRACKING.column.radius);

    updateCombinedSelection();
    updateColumnVisualization();
    computeProfiles();
}

function selectParticlesInPlane() {
    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) return;

    const zMin = TRACKING.plane.zMin;
    const zMax = TRACKING.plane.zMax;

    TRACKING.plane.selectedIds.clear();

    data.particles.forEach(p => {
        if (p.is_ghost === 1) return;
        if (p.z >= zMin && p.z <= zMax) {
            TRACKING.plane.selectedIds.add(p.id);
        }
    });

    TRACKING.plane.active = TRACKING.plane.selectedIds.size > 0;
    console.log(`Plane: ${TRACKING.plane.selectedIds.size} particles selected`);

    // Update UI
    const zDisplay = document.getElementById('plane-z-value');
    if (zDisplay) zDisplay.textContent = `${zMin.toFixed(2)} to ${zMax.toFixed(2)} pc`;

    updateCombinedSelection();
    updatePlaneVisualization();
    computeProfiles();
}

// ============================================================
// Visualization Updates
// ============================================================

function updateColumnVisualization() {
    removeColumnVisualization();

    if (!TRACKING.column.active || TRACKING.column.selectedIds.size === 0) return;

    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) return;

    let zMin = Infinity, zMax = -Infinity;
    let sumX = 0, sumY = 0, count = 0;

    data.particles.forEach(p => {
        if (TRACKING.column.selectedIds.has(p.id)) {
            zMin = Math.min(zMin, p.z);
            zMax = Math.max(zMax, p.z);
            sumX += p.x;
            sumY += p.y;
            count++;
        }
    });

    if (count === 0) return;

    const cx = sumX / count;
    const cy = sumY / count;
    const r = TRACKING.column.radius;
    const height = Math.max(0.5, (zMax - zMin) * 1.3);
    const zCenter = (zMax + zMin) / 2;

    const geometry = new THREE.CylinderGeometry(r, r, height, 32, 1, true);
    const material = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.12,
        side: THREE.DoubleSide,
        depthWrite: false
    });

    TRACKING.column.mesh = new THREE.Mesh(geometry, material);
    TRACKING.column.mesh.rotation.x = Math.PI / 2;
    TRACKING.column.mesh.position.set(cx, cy, zCenter);

    const group = STATE.orbitalPlaneGroup || STATE.scene;
    group.add(TRACKING.column.mesh);
}

function removeColumnVisualization() {
    if (TRACKING.column.mesh) {
        const group = STATE.orbitalPlaneGroup || STATE.scene;
        group.remove(TRACKING.column.mesh);
        TRACKING.column.mesh.geometry.dispose();
        TRACKING.column.mesh.material.dispose();
        TRACKING.column.mesh = null;
    }
}

function updatePlaneVisualization() {
    removePlaneVisualization();

    if (!TRACKING.plane.active || TRACKING.plane.selectedIds.size === 0) return;

    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) return;

    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;
    let sumZ = 0, count = 0;

    data.particles.forEach(p => {
        if (TRACKING.plane.selectedIds.has(p.id)) {
            xMin = Math.min(xMin, p.x);
            xMax = Math.max(xMax, p.x);
            yMin = Math.min(yMin, p.y);
            yMax = Math.max(yMax, p.y);
            sumZ += p.z;
            count++;
        }
    });

    if (count === 0) return;

    const meanZ = sumZ / count;
    const width = Math.max(0.5, (xMax - xMin) * 1.3);
    const height = Math.max(0.5, (yMax - yMin) * 1.3);
    const cx = (xMax + xMin) / 2;
    const cy = (yMax + yMin) / 2;

    const geometry = new THREE.PlaneGeometry(width, height);
    const material = new THREE.MeshBasicMaterial({
        color: 0x00ff88,
        transparent: true,
        opacity: 0.1,
        side: THREE.DoubleSide,
        depthWrite: false
    });

    TRACKING.plane.mesh = new THREE.Mesh(geometry, material);
    TRACKING.plane.mesh.position.set(cx, cy, meanZ);

    const group = STATE.orbitalPlaneGroup || STATE.scene;
    group.add(TRACKING.plane.mesh);

    const edgeGeom = new THREE.EdgesGeometry(geometry);
    const edgeMat = new THREE.LineBasicMaterial({ color: 0x00ff88, opacity: 0.6, transparent: true });
    TRACKING.plane.edges = new THREE.LineSegments(edgeGeom, edgeMat);
    TRACKING.plane.edges.position.set(cx, cy, meanZ);
    group.add(TRACKING.plane.edges);
}

function removePlaneVisualization() {
    const group = STATE.orbitalPlaneGroup || STATE.scene;

    if (TRACKING.plane.mesh) {
        group.remove(TRACKING.plane.mesh);
        TRACKING.plane.mesh.geometry.dispose();
        TRACKING.plane.mesh.material.dispose();
        TRACKING.plane.mesh = null;
    }
    if (TRACKING.plane.edges) {
        group.remove(TRACKING.plane.edges);
        TRACKING.plane.edges.geometry.dispose();
        TRACKING.plane.edges.material.dispose();
        TRACKING.plane.edges = null;
    }
}

// ============================================================
// Combined Selection
// ============================================================

function updateCombinedSelection() {
    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) {
        STATE.selectedParticles = [];
        return;
    }

    STATE.selectedParticles = [];
    data.particles.forEach((p, idx) => {
        if (TRACKING.column.selectedIds.has(p.id) || TRACKING.plane.selectedIds.has(p.id)) {
            STATE.selectedParticles.push(idx);
        }
    });

    updateVisualization(STATE.currentFrame);
}

// ============================================================
// Clear Functions
// ============================================================

function clearProfileSelection() {
    console.log('Clearing profile selections...');

    TRACKING.column.selectedIds.clear();
    TRACKING.column.center = null;
    TRACKING.column.active = false;
    removeColumnVisualization();

    TRACKING.plane.selectedIds.clear();
    TRACKING.plane.active = false;
    TRACKING.plane.zPosition = undefined;
    removePlaneVisualization();

    // Also remove draggable plane preview but recreate if still in selection mode
    removeDraggablePlanePreview();
    if (TRACKING.selectionModeActive && TRACKING.displayMode === 'plane') {
        createDraggablePlane();
    }

    TRACKING.profiles = null;
    STATE.selectedParticles = [];

    if (STATE.snapshots && STATE.snapshots.length > 0) {
        updateVisualization(STATE.currentFrame);
    }
    updateProfileCanvas();

    // Reset z display
    const zDisplay = document.getElementById('plane-z-value');
    if (zDisplay) zDisplay.textContent = '--';

    console.log('Profile selections cleared');
}

// ============================================================
// Parameter Updates (for manual adjustment)
// ============================================================

function updateColumnRadius(value) {
    TRACKING.column.radius = parseFloat(value);
    const el = document.getElementById('column-radius-value');
    if (el) el.textContent = parseFloat(value).toFixed(2) + ' pc';

    if (TRACKING.column.center) {
        selectParticlesInColumn();
    }
}

function setColorVariable(variable) {
    TRACKING.colorVariable = variable;
    updateProfileCanvas();
}

// ============================================================
// Physics Computations
// ============================================================

function computeParticlePhysics(p, comVel) {
    const r = Math.sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    const a_BH = PHYS.G * PHYS.M_BH / (r * r);
    const a_BH_cgs = a_BH * PHYS.pc_to_cm / (PHYS.Myr_to_s * PHYS.Myr_to_s);

    const n_H2 = p.dens * CONFIG.densityToNH2;
    const T = p.temp;

    const dvx = p.vx - comVel.x;
    const dvy = p.vy - comVel.y;
    const dvz = p.vz - comVel.z;
    const v_rel = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

    const mach = p.sound > 0 ? v_rel / p.sound : 0;
    const gamma = 5 / 3;
    const entropy = p.pres / Math.pow(p.dens, gamma);

    let coolingEff = 0;
    if (T > 0 && n_H2 > 0) {
        const Lambda = 1e-27 * Math.pow(T / 10, 2);
        const t_cool = (1.5 * PHYS.k_B * T) / (n_H2 * Lambda);
        const rho_cgs = p.dens * PHYS.Msun_to_g / Math.pow(PHYS.pc_to_cm, 3);
        const t_dyn = 1.0 / Math.sqrt(6.674e-8 * rho_cgs);
        coolingEff = t_dyn / t_cool;
    }

    return {
        r, z: p.z, x: p.x, y: p.y,
        vz: p.vz, v_rel,
        density: n_H2, temperature: T, entropy, mach,
        a_BH: a_BH_cgs, coolingEff, pressure: p.pres
    };
}

// ============================================================
// Profile Computation
// ============================================================

function computeProfiles() {
    const selectedIds = TRACKING.displayMode === 'column'
        ? TRACKING.column.selectedIds
        : TRACKING.plane.selectedIds;

    if (selectedIds.size === 0) {
        TRACKING.profiles = null;
        updateProfileCanvas();
        return;
    }

    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) return;

    const found = [];
    data.particles.forEach(p => {
        if (selectedIds.has(p.id)) found.push(p);
    });

    if (found.length === 0) {
        TRACKING.profiles = { count: 0 };
        updateProfileCanvas();
        return;
    }

    let comVx = 0, comVy = 0, comVz = 0, totalMass = 0;
    for (const p of found) {
        comVx += p.vx * p.mass;
        comVy += p.vy * p.mass;
        comVz += p.vz * p.mass;
        totalMass += p.mass;
    }
    const comVel = {
        x: comVx / totalMass,
        y: comVy / totalMass,
        z: comVz / totalMass
    };

    const particles = found.map(p => computeParticlePhysics(p, comVel));

    TRACKING.profiles = { count: particles.length, particles, comVel };
    updateProfileCanvas();
}

// ============================================================
// Canvas Drawing with FIXED RANGES
// ============================================================

function updateProfileCanvas() {
    const canvas = document.getElementById('profile-canvas');
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    const W = canvas.width;
    const H = canvas.height;

    ctx.fillStyle = '#12141a';
    ctx.fillRect(0, 0, W, H);

    const columnCount = TRACKING.column.selectedIds.size;
    const planeCount = TRACKING.plane.selectedIds.size;

    if (!TRACKING.profiles || TRACKING.profiles.count === 0) {
        ctx.fillStyle = '#aaa';
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'center';
        const msg = TRACKING.displayMode === 'column'
            ? 'Click to place column center'
            : 'Drag the green plane to select';
        ctx.fillText(msg, W / 2, H / 2 - 10);
        ctx.fillStyle = '#888';
        ctx.font = '9px sans-serif';
        ctx.fillText(`Column: ${columnCount} | Plane: ${planeCount}`, W / 2, H / 2 + 10);
        return;
    }

    const p = TRACKING.profiles.particles;

    if (TRACKING.displayMode === 'column') {
        drawColumnProfiles(ctx, W, H, p);
    } else {
        drawPlaneProfiles(ctx, W, H, p);
    }

    ctx.fillStyle = '#bbb';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(`N=${TRACKING.profiles.count}`, 5, H - 5);
    ctx.textAlign = 'right';
    ctx.fillText(`Col:${columnCount} Pln:${planeCount}`, W - 5, H - 5);
}

function drawColumnProfiles(ctx, W, H, particles) {
    const margin = { left: 50, right: 10, top: 8, bottom: 25, gap: 4 };
    const numPlots = 6;
    const plotH = (H - margin.top - margin.bottom - (numPlots - 1) * margin.gap) / numPlots;
    const plotW = W - margin.left - margin.right;

    // Use GLOBAL ranges for z-axis and all y-axes
    const gr = TRACKING.globalRanges;
    const zMin = gr ? gr.z.min : Math.min(...particles.map(pt => pt.z));
    const zMax = gr ? gr.z.max : Math.max(...particles.map(pt => pt.z));

    const variables = [
        { key: 'density', label: 'log n_H2', color: '#ff9966', log: true,
          yMin: gr?.density.min, yMax: gr?.density.max },
        { key: 'temperature', label: 'log T', color: '#ff6688', log: true,
          yMin: gr?.temperature.min, yMax: gr?.temperature.max },
        { key: 'vz', label: 'v_z', color: '#66ff88', log: false,
          yMin: gr?.vz.min, yMax: gr?.vz.max },
        { key: 'mach', label: 'Mach', color: '#ffff66', log: false,
          yMin: gr?.mach.min ?? 0, yMax: gr?.mach.max },
        { key: 'a_BH', label: 'log a_BH', color: '#66aaff', log: true,
          yMin: gr?.a_BH.min, yMax: gr?.a_BH.max },
        { key: 'coolingEff', label: 'log η_cool', color: '#ff66ff', log: true,
          yMin: gr?.coolingEff.min, yMax: gr?.coolingEff.max }
    ];

    variables.forEach((v, i) => {
        const y0 = margin.top + i * (plotH + margin.gap);
        drawScatterPlotFixed(ctx, margin.left, y0, plotW, plotH, particles,
            'z', v.key, v.label, v.color, v.log, zMin, zMax, v.yMin, v.yMax);
    });

    ctx.fillStyle = '#bbb';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('z (pc)', W / 2, H - 8);
}

function drawPlaneProfiles(ctx, W, H, particles) {
    const margin = { left: 45, right: 10, top: 8, bottom: 25, gap: 8 };
    const profileH = 60;
    const profileW = (W - margin.left - margin.right - margin.gap) / 2;
    const scatterH = H - margin.top - profileH - margin.gap - margin.bottom - 15;
    const scatterW = W - margin.left - margin.right;

    const gr = TRACKING.globalRanges;
    const xMin = gr ? gr.x.min : Math.min(...particles.map(pt => pt.x));
    const xMax = gr ? gr.x.max : Math.max(...particles.map(pt => pt.x));
    const yMin = gr ? gr.y.min : Math.min(...particles.map(pt => pt.y));
    const yMax = gr ? gr.y.max : Math.max(...particles.map(pt => pt.y));

    const varConfig = {
        density: { label: 'log n_H2', color: '#ff9966', log: true,
                   yMin: gr?.density.min, yMax: gr?.density.max },
        temperature: { label: 'log T', color: '#ff6688', log: true,
                       yMin: gr?.temperature.min, yMax: gr?.temperature.max },
        mach: { label: 'Mach', color: '#ffff66', log: false,
                yMin: gr?.mach.min ?? 0, yMax: gr?.mach.max },
        entropy: { label: 'log s', color: '#66aaff', log: true,
                   yMin: gr?.entropy.min, yMax: gr?.entropy.max },
        coolingEff: { label: 'log η', color: '#ff66ff', log: true,
                      yMin: gr?.coolingEff.min, yMax: gr?.coolingEff.max }
    };
    const vc = varConfig[TRACKING.colorVariable] || varConfig.density;

    drawScatterPlotFixed(ctx, margin.left, margin.top, profileW, profileH,
        particles, 'x', TRACKING.colorVariable, vc.label, vc.color, vc.log,
        xMin, xMax, vc.yMin, vc.yMax);

    drawScatterPlotFixed(ctx, margin.left + profileW + margin.gap, margin.top, profileW, profileH,
        particles, 'y', TRACKING.colorVariable, vc.label, vc.color, vc.log,
        yMin, yMax, vc.yMin, vc.yMax);

    const scatterY = margin.top + profileH + margin.gap + 10;
    draw2DScatterFixed(ctx, margin.left, scatterY, scatterW, scatterH,
        particles, xMin, xMax, yMin, yMax, TRACKING.colorVariable, vc.log, vc.yMin, vc.yMax);

    ctx.fillStyle = '#bbb';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('x (pc)', margin.left + profileW / 2, margin.top + profileH + 12);
    ctx.fillText('y (pc)', margin.left + profileW + margin.gap + profileW / 2, margin.top + profileH + 12);
    ctx.fillText('x (pc)', W / 2, H - 8);

    ctx.save();
    ctx.translate(12, scatterY + scatterH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('y (pc)', 0, 0);
    ctx.restore();
}

function drawScatterPlotFixed(ctx, x0, y0, w, h, particles, xKey, yKey, label, color, useLog, xMin, xMax, fixedYMin, fixedYMax) {
    ctx.fillStyle = '#0a0c12';
    ctx.fillRect(x0, y0, w, h);

    let yVals = particles.map(pt => pt[yKey]);
    let xVals = particles.map(pt => pt[xKey]);

    if (useLog) {
        yVals = yVals.map(v => v > 0 ? Math.log10(v) : -10);
    }

    // Use FIXED ranges if provided, otherwise compute from data
    let yMin, yMax;
    if (fixedYMin !== undefined && fixedYMax !== undefined && isFinite(fixedYMin) && isFinite(fixedYMax)) {
        yMin = fixedYMin;
        yMax = fixedYMax;
    } else {
        const validY = yVals.filter(v => isFinite(v));
        if (validY.length === 0) return;
        yMin = Math.min(...validY);
        yMax = Math.max(...validY);
    }

    const yRange = yMax - yMin || 1;
    const xRange = xMax - xMin || 1;

    // Grid
    ctx.strokeStyle = '#222';
    ctx.lineWidth = 0.5;
    for (let i = 0; i <= 2; i++) {
        const gy = y0 + h * i / 2;
        ctx.beginPath();
        ctx.moveTo(x0, gy);
        ctx.lineTo(x0 + w, gy);
        ctx.stroke();
    }

    // Zero line for linear scales
    if (!useLog && yMin < 0 && yMax > 0) {
        const zeroY = y0 + h * (1 - (0 - yMin) / yRange);
        ctx.strokeStyle = '#444';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(x0, zeroY);
        ctx.lineTo(x0 + w, zeroY);
        ctx.stroke();
    }

    // Points
    ctx.fillStyle = color;
    ctx.globalAlpha = 0.5;
    for (let i = 0; i < particles.length; i++) {
        const px = x0 + w * (xVals[i] - xMin) / xRange;
        const py = y0 + h * (1 - (yVals[i] - yMin) / yRange);
        if (isFinite(px) && isFinite(py) && px >= x0 && px <= x0+w && py >= y0 && py <= y0+h) {
            ctx.beginPath();
            ctx.arc(px, py, 1.5, 0, Math.PI * 2);
            ctx.fill();
        }
    }
    ctx.globalAlpha = 1;

    // Running median
    drawRunningMedianFixed(ctx, x0, y0, w, h, xVals, yVals, xMin, xMax, yMin, yMax, color);

    // Label
    ctx.fillStyle = color;
    ctx.font = 'bold 8px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(label, x0 + 3, y0 + 9);

    // Y-axis labels
    ctx.fillStyle = '#aaa';
    ctx.font = '7px sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(yMax.toFixed(1), x0 - 2, y0 + 8);
    ctx.fillText(yMin.toFixed(1), x0 - 2, y0 + h - 2);
}

function drawRunningMedianFixed(ctx, x0, y0, w, h, xVals, yVals, xMin, xMax, yMin, yMax, color) {
    const sorted = xVals.map((x, i) => ({ x, y: yVals[i] }))
        .filter(p => isFinite(p.x) && isFinite(p.y))
        .sort((a, b) => a.x - b.x);

    if (sorted.length < 5) return;

    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const numBins = Math.min(15, Math.floor(sorted.length / 3));
    const binSize = sorted.length / numBins;

    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
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

        // Clip to plot area
        const clippedPy = Math.max(y0, Math.min(y0 + h, py));

        if (!started) { ctx.moveTo(px, clippedPy); started = true; }
        else { ctx.lineTo(px, clippedPy); }
    }
    ctx.stroke();
}

function draw2DScatterFixed(ctx, x0, y0, w, h, particles, xMin, xMax, yMin, yMax, colorVar, useLog, cMin, cMax) {
    ctx.fillStyle = '#0a0c12';
    ctx.fillRect(x0, y0, w, h);

    // Grid
    ctx.strokeStyle = '#222';
    ctx.lineWidth = 0.5;
    for (let i = 0; i <= 4; i++) {
        const gx = x0 + w * i / 4;
        const gy = y0 + h * i / 4;
        ctx.beginPath(); ctx.moveTo(gx, y0); ctx.lineTo(gx, y0 + h); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(x0, gy); ctx.lineTo(x0 + w, gy); ctx.stroke();
    }

    let colorVals = particles.map(pt => pt[colorVar]);
    if (useLog) colorVals = colorVals.map(v => v > 0 ? Math.log10(v) : -10);

    // Use fixed color range
    const colorMin = (cMin !== undefined && isFinite(cMin)) ? cMin : Math.min(...colorVals.filter(v => isFinite(v)));
    const colorMax = (cMax !== undefined && isFinite(cMax)) ? cMax : Math.max(...colorVals.filter(v => isFinite(v)));
    const cRange = colorMax - colorMin || 1;

    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;

    for (let i = 0; i < particles.length; i++) {
        const pt = particles[i];
        const px = x0 + w * (pt.x - xMin) / xRange;
        const py = y0 + h * (1 - (pt.y - yMin) / yRange);
        const t = (colorVals[i] - colorMin) / cRange;

        if (!isFinite(px) || !isFinite(py) || !isFinite(t)) continue;
        if (px < x0 || px > x0+w || py < y0 || py > y0+h) continue;

        const color = valueToColor(Math.max(0, Math.min(1, t)));
        ctx.fillStyle = `rgb(${Math.floor(color.r*255)}, ${Math.floor(color.g*255)}, ${Math.floor(color.b*255)})`;
        ctx.beginPath();
        ctx.arc(px, py, 2.5, 0, Math.PI * 2);
        ctx.fill();
    }

    // Colorbar
    const cbW = 12, cbH = h - 20, cbX = x0 + w + 5, cbY = y0 + 10;
    for (let i = 0; i < cbH; i++) {
        const t = 1 - i / cbH;
        const color = valueToColor(t);
        ctx.fillStyle = `rgb(${Math.floor(color.r*255)}, ${Math.floor(color.g*255)}, ${Math.floor(color.b*255)})`;
        ctx.fillRect(cbX, cbY + i, cbW, 1);
    }

    ctx.fillStyle = '#bbb';
    ctx.font = '7px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(colorMax.toFixed(1), cbX + cbW + 2, cbY + 6);
    ctx.fillText(colorMin.toFixed(1), cbX + cbW + 2, cbY + cbH);

    ctx.fillStyle = '#aaa';
    ctx.textAlign = 'center';
    ctx.fillText(xMin.toFixed(1), x0, y0 + h + 10);
    ctx.fillText(xMax.toFixed(1), x0 + w, y0 + h + 10);
    ctx.textAlign = 'right';
    ctx.fillText(yMax.toFixed(1), x0 - 2, y0 + 8);
    ctx.fillText(yMin.toFixed(1), x0 - 2, y0 + h);
}

// ============================================================
// Frame Change Hooks
// ============================================================

function syncSelectedParticlesBeforeRender() {
    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) {
        STATE.selectedParticles = [];
        return;
    }

    STATE.selectedParticles = [];
    data.particles.forEach((p, idx) => {
        if (TRACKING.column.selectedIds.has(p.id) || TRACKING.plane.selectedIds.has(p.id)) {
            STATE.selectedParticles.push(idx);
        }
    });
}

function onFrameChangeUpdateTracking() {
    if (TRACKING.column.active) {
        updateColumnVisualization();
    }
    if (TRACKING.plane.active) {
        updatePlaneVisualization();
    }

    if (TRACKING.column.selectedIds.size > 0 || TRACKING.plane.selectedIds.size > 0) {
        computeProfiles();
    }
}

// ============================================================
// Sync with PV Selection
// ============================================================

function syncTrackingFromPVSelection() {
    const data = STATE.snapshots[STATE.currentFrame];
    if (!data) return;

    const targetIds = TRACKING.displayMode === 'column'
        ? TRACKING.column.selectedIds
        : TRACKING.plane.selectedIds;

    for (const idx of STATE.selectedParticles) {
        const p = data.particles[idx];
        if (p && p.id !== undefined) {
            targetIds.add(p.id);
        }
    }

    if (TRACKING.displayMode === 'column') {
        TRACKING.column.active = targetIds.size > 0;
    } else {
        TRACKING.plane.active = targetIds.size > 0;
    }

    computeProfiles();
}

// ============================================================
// Initialization
// ============================================================

function initTracking() {
    setTrackingMode('column');
    initSelectionHandlers();

    // Compute global ranges after a short delay (data should be loaded)
    setTimeout(() => {
        if (STATE.snapshots && STATE.snapshots.length > 0 && !TRACKING.globalRanges) {
            computeGlobalProfileRanges();
        }
    }, 1000);

    console.log('Tracking initialized');
}

// Hook to compute ranges when data loads
function onDataLoadedComputeRanges() {
    computeGlobalProfileRanges();
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initTracking);
} else {
    setTimeout(initTracking, 100);
}
