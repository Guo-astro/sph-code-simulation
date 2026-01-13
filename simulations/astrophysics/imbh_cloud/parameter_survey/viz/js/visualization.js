// ============================================================
// IMBH-Cloud Visualization - Particle Rendering & Colorbar
// ============================================================

let raycaster, mouse, losDragPlane;

function initScene() {
    // Scene
    STATE.scene = new THREE.Scene();
    STATE.scene.background = new THREE.Color(0x0a0a15);

    // Camera
    STATE.camera = new THREE.PerspectiveCamera(
        60, window.innerWidth / window.innerHeight, 0.01, 1000
    );
    STATE.camera.position.set(0, 0, 25);
    STATE.camera.lookAt(0, 0, 0);

    // Renderer
    STATE.renderer = new THREE.WebGLRenderer({ antialias: true });
    STATE.renderer.setSize(window.innerWidth, window.innerHeight);
    STATE.renderer.setPixelRatio(window.devicePixelRatio);
    document.getElementById('container').appendChild(STATE.renderer.domElement);

    // Controls
    STATE.controls = new THREE.OrbitControls(STATE.camera, STATE.renderer.domElement);
    STATE.controls.enableDamping = true;
    STATE.controls.dampingFactor = 0.05;
    STATE.controls.minDistance = 1;
    STATE.controls.maxDistance = 100;

    // Raycaster for LOS interaction
    raycaster = new THREE.Raycaster();
    mouse = new THREE.Vector2();

    // Lighting
    const ambientLight = new THREE.AmbientLight(0x404040, 0.5);
    STATE.scene.add(ambientLight);

    // Black hole sphere = accretion/sink radius (0.01 pc)
    const sinkRadius = 0.01;
    const bhGeometry = new THREE.SphereGeometry(sinkRadius, 32, 32);
    const bhMaterial = new THREE.MeshBasicMaterial({
        color: 0xff4444,
        transparent: true,
        opacity: 0.7
    });
    STATE.bhMesh = new THREE.Mesh(bhGeometry, bhMaterial);
    STATE.bhMesh.position.set(0, 0, 0);
    STATE.scene.add(STATE.bhMesh);

    // Grid helper (orbital plane)
    const gridHelper = new THREE.GridHelper(40, 20, 0x2a2a4a, 0x1a1a3a);
    gridHelper.rotation.x = Math.PI / 2;
    STATE.scene.add(gridHelper);

    // Axes helper
    const axesHelper = new THREE.AxesHelper(5);
    STATE.scene.add(axesHelper);

    // COM marker removed per user request

    // Tidal and Hill circles removed per user request

    // Create analytic orbit
    createAnalyticOrbit();

    // Create draggable LOS indicator
    createLOSControl();

    // Initialize colorbar
    updateColorbar();

    // Add pointer event listeners for LOS dragging (pointer events work better with Three.js)
    // Use window for move/up to capture even when pointer leaves the canvas
    STATE.renderer.domElement.addEventListener('pointerdown', onPointerDown);
    window.addEventListener('pointermove', onPointerMove);
    window.addEventListener('pointerup', onPointerUp);

    // Add hover effect
    STATE.renderer.domElement.addEventListener('pointermove', onPointerHover);
}

function createAnalyticOrbit() {
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

    STATE.orbitLine = new THREE.Line(geometry, material);
    STATE.scene.add(STATE.orbitLine);

    // Pericenter position
    const periPos = new THREE.Vector3(
        CONFIG.r_peri * e_hat[0],
        CONFIG.r_peri * e_hat[1],
        CONFIG.r_peri * e_hat[2]
    );

    // Pericenter marker
    if (STATE.periMarker) {
        STATE.scene.remove(STATE.periMarker);
        STATE.periMarker.geometry.dispose();
        STATE.periMarker.material.dispose();
    }
    const periGeom = new THREE.RingGeometry(0.08, 0.12, 32);
    const periMat = new THREE.MeshBasicMaterial({ color: 0x00ff88, side: THREE.DoubleSide });
    STATE.periMarker = new THREE.Mesh(periGeom, periMat);
    STATE.periMarker.position.copy(periPos);
    STATE.periMarker.lookAt(0, 0, 0);
    STATE.scene.add(STATE.periMarker);

    console.log(`Orbit created: r_peri=${CONFIG.r_peri} pc, p=${CONFIG.p.toFixed(4)} pc, pericenter at (${periPos.x.toFixed(3)}, ${periPos.y.toFixed(3)}, ${periPos.z.toFixed(3)})`);
}

function createImpactParameterViz() {
    // Initial cloud position and velocity from config
    const pos0 = new THREE.Vector3(
        CONFIG.cloud_pos0[0],
        CONFIG.cloud_pos0[1],
        CONFIG.cloud_pos0[2]
    );
    const vel0 = new THREE.Vector3(
        CONFIG.cloud_vel0[0],
        CONFIG.cloud_vel0[1],
        CONFIG.cloud_vel0[2]
    );
    const vel0_norm = vel0.clone().normalize();

    // Initial velocity arrow (yellow)
    const velArrow = new THREE.ArrowHelper(
        vel0_norm,
        pos0,
        3,  // length
        0xffff00,  // yellow
        0.4,
        0.2
    );
    STATE.scene.add(velArrow);

    // Extend velocity line backward to show approach direction
    const lineGeom = new THREE.BufferGeometry().setFromPoints([
        pos0.clone().sub(vel0_norm.clone().multiplyScalar(8)),  // extended backward
        pos0.clone().add(vel0_norm.clone().multiplyScalar(5))   // extended forward
    ]);
    const lineMat = new THREE.LineDashedMaterial({
        color: 0xffff00,
        dashSize: 0.3,
        gapSize: 0.15,
        transparent: true,
        opacity: 0.5
    });
    const velLine = new THREE.Line(lineGeom, lineMat);
    velLine.computeLineDistances();
    STATE.scene.add(velLine);

    // Impact parameter: perpendicular distance from BH to asymptotic velocity line
    // For a line through pos0 with direction vel0_norm:
    // Distance from origin = |pos0 - (pos0 · vel0_norm) * vel0_norm|
    const projLength = pos0.dot(vel0_norm);
    const projPoint = vel0_norm.clone().multiplyScalar(projLength);
    const perpVector = pos0.clone().sub(projPoint);
    const impactParam = perpVector.length();

    // Draw impact parameter line (cyan dashed)
    const bLineGeom = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(0, 0, 0),
        projPoint.clone().negate()  // closest point on extended line to origin
    ]);
    const bLineMat = new THREE.LineDashedMaterial({
        color: 0x00ffff,
        dashSize: 0.2,
        gapSize: 0.1,
        linewidth: 2
    });
    const bLine = new THREE.Line(bLineGeom, bLineMat);
    bLine.computeLineDistances();
    STATE.scene.add(bLine);

    // Label for impact parameter (using sprite or text)
    console.log(`Impact parameter b = ${impactParam.toFixed(3)} pc (from velocity asymptote)`);

    // Pericenter markers removed per user request
}

function createLOSControl() {
    // Create a group for LOS visualization
    STATE.losArrows = new THREE.Group();

    // Draggable sphere on a unit sphere around origin - make it larger for easier clicking
    const sphereGeom = new THREE.SphereGeometry(0.5, 16, 16);  // Larger radius
    const sphereMat = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.9
    });
    STATE.losSphere = new THREE.Mesh(sphereGeom, sphereMat);
    STATE.losSphere.name = 'losSphere';

    // Position on unit sphere in direction of losVector
    const radius = 5;  // Distance from origin for draggable sphere
    STATE.losSphere.position.copy(STATE.losVector.clone().multiplyScalar(radius));

    // Add sphere directly to scene (not to group) for easier raycasting
    STATE.scene.add(STATE.losSphere);

    // Arrow showing LOS direction
    updateLOSArrows();

    STATE.scene.add(STATE.losArrows);
}

function updateLOSArrows() {
    // Remove old arrows from group
    while (STATE.losArrows.children.length > 0) {
        STATE.losArrows.remove(STATE.losArrows.children[0]);
    }

    // Arrow from origin to sphere
    const arrowLength = 4.5;
    const losArrow = new THREE.ArrowHelper(
        STATE.losVector.clone(),
        new THREE.Vector3(0, 0, 0),
        arrowLength,
        0x00ffff,
        0.5,
        0.25
    );
    STATE.losArrows.add(losArrow);

    // Update sphere position (sphere is in scene, not group)
    if (STATE.losSphere) {
        STATE.losSphere.position.copy(STATE.losVector.clone().multiplyScalar(5));
    }

    // Update LOS info display
    updateLOSInfo();
}

function updateLOSInfo() {
    const v = STATE.losVector;
    const theta = Math.acos(v.z) * 180 / Math.PI;
    const phi = Math.atan2(v.y, v.x) * 180 / Math.PI;

    // Update dropdown to show custom if not matching preset
    const select = document.getElementById('los-direction');
    if (Math.abs(v.x - 1) < 0.01 && Math.abs(v.y) < 0.01 && Math.abs(v.z) < 0.01) {
        select.value = 'x';
    } else if (Math.abs(v.x) < 0.01 && Math.abs(v.y - 1) < 0.01 && Math.abs(v.z) < 0.01) {
        select.value = 'y';
    } else if (Math.abs(v.x) < 0.01 && Math.abs(v.y) < 0.01 && Math.abs(v.z - 1) < 0.01) {
        select.value = 'z';
    }
    // Otherwise leave as is (custom direction)
}

// Pointer handlers for LOS dragging (more reliable than mouse events)
function onPointerDown(event) {
    // Only handle left mouse button
    if (event.button !== 0) return;

    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

    raycaster.setFromCamera(mouse, STATE.camera);

    // Check if losSphere exists
    if (!STATE.losSphere) {
        console.log('LOS sphere not found');
        return;
    }

    // Update world matrix before raycasting
    STATE.losSphere.updateMatrixWorld();

    // Intersect with the sphere
    const intersects = raycaster.intersectObject(STATE.losSphere, false);

    if (intersects.length > 0) {
        event.preventDefault();
        event.stopPropagation();

        STATE.isDraggingLOS = true;
        STATE.controls.enabled = false;
        STATE.losSphere.material.color.setHex(0xff00ff);  // Highlight when dragging
        document.body.style.cursor = 'grabbing';
        console.log('Started dragging LOS');
    }
}

function onPointerMove(event) {
    if (!STATE.isDraggingLOS) return;

    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

    raycaster.setFromCamera(mouse, STATE.camera);

    // Project to imaginary sphere of radius 5 centered at origin
    const sphereRadius = 5;
    const ray = raycaster.ray;

    // Ray-sphere intersection: solve |origin + t*dir|^2 = R^2
    const oc = ray.origin.clone();
    const d = ray.direction.clone();

    const a = d.dot(d);
    const b = 2.0 * oc.dot(d);
    const c = oc.dot(oc) - sphereRadius * sphereRadius;
    const discriminant = b * b - 4 * a * c;

    if (discriminant >= 0) {
        // Use the closer intersection
        const t = (-b - Math.sqrt(discriminant)) / (2 * a);
        if (t > 0) {
            const hitPoint = ray.origin.clone().add(d.multiplyScalar(t));

            // Update LOS vector (normalize to unit vector)
            STATE.losVector.copy(hitPoint.normalize());

            // Update visualization
            updateLOSArrows();

            // Update PV diagram
            if (STATE.snapshots.length > 0) {
                updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
            }
        }
    }
}

function onPointerUp(event) {
    if (STATE.isDraggingLOS) {
        STATE.isDraggingLOS = false;
        STATE.controls.enabled = true;
        STATE.losSphere.material.color.setHex(0x00ffff);  // Reset color
        document.body.style.cursor = 'default';
    }
}

function onPointerHover(event) {
    if (STATE.isDraggingLOS) return;
    if (!STATE.losSphere) return;

    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

    raycaster.setFromCamera(mouse, STATE.camera);
    STATE.losSphere.updateMatrixWorld();

    const intersects = raycaster.intersectObject(STATE.losSphere, false);

    if (intersects.length > 0) {
        document.body.style.cursor = 'grab';
        STATE.losSphere.material.opacity = 1.0;
    } else {
        document.body.style.cursor = 'default';
        STATE.losSphere.material.opacity = 0.9;
    }
}

// Convert normalized value [0,1] to HSL color (blue to red)
function valueToColor(t) {
    const hue = 0.7 * (1 - t);
    return new THREE.Color().setHSL(hue, 1, 0.5);
}

function updateVisualization(frameIndex) {
    if (frameIndex < 0 || frameIndex >= STATE.snapshots.length) return;

    STATE.currentFrame = frameIndex;
    const data = STATE.snapshots[frameIndex];

    // Remove old particle system
    if (STATE.particleSystem) {
        STATE.scene.remove(STATE.particleSystem);
        STATE.particleSystem.geometry.dispose();
        STATE.particleSystem.material.dispose();
    }

    // First pass: compute COM position and velocity (for velocity dispersion)
    let comX = 0, comY = 0, comZ = 0;
    let comVx = 0, comVy = 0, comVz = 0;
    let totalMass = 0;

    for (const p of data.particles) {
        if (p.is_ghost === 1) continue;
        comX += p.x * p.mass;
        comY += p.y * p.mass;
        comZ += p.z * p.mass;
        comVx += p.vx * p.mass;
        comVy += p.vy * p.mass;
        comVz += p.vz * p.mass;
        totalMass += p.mass;
    }

    comX /= totalMass;
    comY /= totalMass;
    comZ /= totalMass;
    comVx /= totalMass;
    comVy /= totalMass;
    comVz /= totalMass;

    const comVelocity = { x: comVx, y: comVy, z: comVz };

    // Prepare arrays
    const positions = new Float32Array(data.particles.length * 3);
    const colors = new Float32Array(data.particles.length * 3);
    const sizes = new Float32Array(data.particles.length);

    // Second pass: set colors using GLOBAL ranges
    data.particles.forEach((p, i) => {
        positions[i*3] = p.x;
        positions[i*3+1] = p.y;
        positions[i*3+2] = p.z;

        let color;

        if (p.is_ghost === 1) {
            color = new THREE.Color(0.4, 0.4, 0.5);
            sizes[i] = 0.02;
        } else {
            const t = getColorForValue(p, comVelocity);
            color = valueToColor(t);
            sizes[i] = 0.03 + 0.02 * t;
        }

        colors[i*3] = color.r;
        colors[i*3+1] = color.g;
        colors[i*3+2] = color.b;
    });

    // Create particle system
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.setAttribute('size', new THREE.BufferAttribute(sizes, 1));

    const material = new THREE.PointsMaterial({
        size: 0.05,
        vertexColors: true,
        transparent: true,
        opacity: 0.8,
        sizeAttenuation: true
    });

    STATE.particleSystem = new THREE.Points(geometry, material);
    STATE.scene.add(STATE.particleSystem);

    // COM distance for display
    const comDist = Math.sqrt(comX*comX + comY*comY + comZ*comZ);

    // Count current cloud particles (non-ghost)
    let currentCloudParticles = 0;
    for (const p of data.particles) {
        if (p.is_ghost !== 1) {
            currentCloudParticles++;
        }
    }
    const accretedCount = STATE.initialCloudParticles - currentCloudParticles;

    // Calculate shock tracers (Godunov SPH - no artificial viscosity)
    let maxMach = 0;
    let supersonicCount = 0;
    let maxEntropy = 0;
    let maxDensity = 0;
    let sumDensity = 0;
    let validParticles = 0;
    const gamma = 1.6667;

    for (const p of data.particles) {
        if (p.is_ghost === 1) continue;
        validParticles++;

        // Mach number: |v - v_COM| / c_s
        const dvx = p.vx - comVelocity.x;
        const dvy = p.vy - comVelocity.y;
        const dvz = p.vz - comVelocity.z;
        const v_rel = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        const mach = v_rel / p.sound;
        if (mach > maxMach) maxMach = mach;
        if (mach > 1) supersonicCount++;

        // Entropy: s = P / rho^gamma
        const entropy = p.pres / Math.pow(p.dens, gamma);
        if (entropy > maxEntropy) maxEntropy = entropy;

        // Density for compression ratio
        if (p.dens > maxDensity) maxDensity = p.dens;
        sumDensity += p.dens;
    }

    const supersonicFrac = (supersonicCount / validParticles * 100).toFixed(1);
    const meanDensity = sumDensity / validParticles;
    const densityRatio = maxDensity / meanDensity;

    // Update shock panel
    document.getElementById('max-mach').textContent = maxMach.toFixed(2);
    document.getElementById('supersonic-frac').textContent = `${supersonicFrac}%`;
    document.getElementById('density-ratio').textContent = densityRatio.toFixed(1);
    document.getElementById('max-entropy').textContent = maxEntropy.toExponential(2);

    // Update UI with physical units
    const timeMyr = data.time * CONFIG.timeToMyr;

    document.getElementById('time-display').textContent = `${timeMyr.toFixed(3)} Myr`;
    document.getElementById('particle-count').textContent = data.particles.length.toLocaleString();
    document.getElementById('accreted-count').textContent = accretedCount.toLocaleString();
    document.getElementById('snapshot-num').textContent = `${frameIndex + 1}`;
    document.getElementById('com-distance').textContent = `${comDist.toFixed(2)} pc`;
    document.getElementById('timeline').value = frameIndex + 1;
    document.getElementById('frame-display').textContent = `${frameIndex + 1} / ${STATE.snapshots.length}`;

    // Update PV diagram
    updatePVDiagram(data);
}

function updateColorbar() {
    const range = STATE.colorRanges[STATE.colorMode];

    let title, minLabel, maxLabel;
    switch (STATE.colorMode) {
        case 'temperature':
            title = 'T (K)';
            // Show as powers of 10
            minLabel = `10^${range.min}`;
            maxLabel = `10^${range.max}`;
            break;
        case 'velocity':
            title = 'σ_v (km/s)';
            minLabel = `${Math.pow(10, range.min).toFixed(0)}`;
            maxLabel = `${Math.pow(10, range.max).toFixed(0)}`;
            break;
        case 'density':
        default:
            title = 'n_H2 (cm⁻³)';
            minLabel = `10^${range.min}`;
            maxLabel = `10^${range.max}`;
            break;
    }

    document.getElementById('colorbar-title').textContent = title;
    document.getElementById('colorbar-min').textContent = minLabel;
    document.getElementById('colorbar-max').textContent = maxLabel;

    // Colorbar gradient: red at top (max) to blue at bottom (min)
    document.getElementById('colorbar').style.background =
        'linear-gradient(to bottom, hsl(0, 100%, 50%), hsl(42, 100%, 50%), hsl(84, 100%, 50%), hsl(126, 100%, 50%), hsl(168, 100%, 50%), hsl(210, 100%, 50%), hsl(252, 100%, 50%))';
}

function animate() {
    requestAnimationFrame(animate);
    STATE.controls.update();
    STATE.renderer.render(STATE.scene, STATE.camera);
}

function onWindowResize() {
    STATE.camera.aspect = window.innerWidth / window.innerHeight;
    STATE.camera.updateProjectionMatrix();
    STATE.renderer.setSize(window.innerWidth, window.innerHeight);
}
