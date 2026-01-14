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

    // Axes helper (basic)
    const axesHelper = new THREE.AxesHelper(5);
    STATE.scene.add(axesHelper);

    // Create labeled coordinate system for LOS understanding
    createLabeledCoordinates();

    // Create galaxy context visualization (hidden by default)
    createGalaxyContext();

    // COM marker removed per user request

    // Tidal and Hill circles removed per user request

    // === ORBITAL PLANE GROUP ===
    // This group contains all elements that rotate with the orbital plane
    // (analytic orbit, particles, COM marker, etc.)
    STATE.orbitalPlaneGroup = new THREE.Group();
    STATE.scene.add(STATE.orbitalPlaneGroup);

    // Set the original orbital pole from simulation's angular momentum
    // This will be updated after CONFIG.L_vec is computed
    updateOriginalOrbitalPole();

    // Create analytic orbit (added to orbitalPlaneGroup)
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

// ============================================================
// Orbital Plane Rotation Functions
// ============================================================

// Set the original orbital pole from simulation's angular momentum vector
function updateOriginalOrbitalPole() {
    if (CONFIG.L_vec && CONFIG.L > 0) {
        // Normalize L_vec to get the orbital pole direction
        STATE.originalOrbitalPole.set(
            CONFIG.L_vec[0] / CONFIG.L,
            CONFIG.L_vec[1] / CONFIG.L,
            CONFIG.L_vec[2] / CONFIG.L
        );

        // Initialize the current orbital pole to the original
        STATE.losVector.copy(STATE.originalOrbitalPole);

        console.log('Original orbital pole set from L_vec:', STATE.originalOrbitalPole);
    } else {
        // Default to +Z if L_vec not available
        STATE.originalOrbitalPole.set(0, 0, 1);
        STATE.losVector.set(0, 0, 1);
        console.log('Original orbital pole defaulted to +Z');
    }
}

// Update the orbital plane group's rotation based on the current orbital pole
// This rotates the entire orbit and particle system to match the new orientation
function updateOrbitalPlaneRotation() {
    if (!STATE.orbitalPlaneGroup) return;

    // Compute the rotation from original orbital pole to current orbital pole
    // This quaternion transforms the original orientation to the new one
    const originalPole = STATE.originalOrbitalPole;
    const currentPole = STATE.losVector;

    // Use quaternion to represent the rotation
    const quaternion = new THREE.Quaternion();
    quaternion.setFromUnitVectors(originalPole, currentPole);

    // Apply the rotation to the orbital plane group
    STATE.orbitalPlaneGroup.quaternion.copy(quaternion);

    console.log('Orbital plane rotated to match pole:', currentPole);
}

// Reset orbital pole to original simulation orientation
function resetOrbitalPoleToOriginal() {
    STATE.losVector.copy(STATE.originalOrbitalPole);
    updateLOSArrows();
    updatePVDiagram();
    console.log('Orbital pole reset to original simulation orientation');
}

// Create text sprite for 3D labels
function createTextSprite(text, color = '#ffffff', fontSize = 48) {
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    canvas.width = 256;
    canvas.height = 128;

    context.fillStyle = 'transparent';
    context.fillRect(0, 0, canvas.width, canvas.height);

    context.font = `Bold ${fontSize}px Arial`;
    context.fillStyle = color;
    context.textAlign = 'center';
    context.textBaseline = 'middle';
    context.fillText(text, canvas.width / 2, canvas.height / 2);

    const texture = new THREE.CanvasTexture(canvas);
    const material = new THREE.SpriteMaterial({ map: texture, transparent: true });
    const sprite = new THREE.Sprite(material);
    sprite.scale.set(2, 1, 1);
    return sprite;
}

// Create labeled coordinate system to help understand LOS convention
function createLabeledCoordinates() {
    const axisLength = 8;
    const labelOffset = 1.2;

    // X-axis label (red) - "East" in sky plane convention
    const xLabel = createTextSprite('X (East)', '#ff6666', 36);
    xLabel.position.set(axisLength * labelOffset, 0, 0);
    STATE.scene.add(xLabel);

    // Y-axis label (green) - "North" in sky plane convention
    const yLabel = createTextSprite('Y (North)', '#66ff66', 36);
    yLabel.position.set(0, axisLength * labelOffset, 0);
    STATE.scene.add(yLabel);

    // Z-axis label (blue) - "LOS" when face-on
    const zLabel = createTextSprite('Z (Face-on)', '#6666ff', 36);
    zLabel.position.set(0, 0, axisLength * labelOffset);
    STATE.scene.add(zLabel);

    // Origin label
    const originLabel = createTextSprite('BH', '#ff4444', 32);
    originLabel.position.set(0.5, 0.5, 0.5);
    STATE.scene.add(originLabel);

    // Add orbital plane label
    const orbitalPlaneLabel = createTextSprite('Orbital Plane (X-Y)', '#888888', 28);
    orbitalPlaneLabel.position.set(6, 6, 0);
    STATE.scene.add(orbitalPlaneLabel);

    // Draw extended axes with arrows
    const arrowSize = 0.3;

    // X-axis arrow (red)
    const xDir = new THREE.Vector3(1, 0, 0);
    const xArrow = new THREE.ArrowHelper(xDir, new THREE.Vector3(0, 0, 0), axisLength, 0xff4444, arrowSize * 2, arrowSize);
    STATE.scene.add(xArrow);

    // Y-axis arrow (green)
    const yDir = new THREE.Vector3(0, 1, 0);
    const yArrow = new THREE.ArrowHelper(yDir, new THREE.Vector3(0, 0, 0), axisLength, 0x44ff44, arrowSize * 2, arrowSize);
    STATE.scene.add(yArrow);

    // Z-axis arrow (blue)
    const zDir = new THREE.Vector3(0, 0, 1);
    const zArrow = new THREE.ArrowHelper(zDir, new THREE.Vector3(0, 0, 0), axisLength, 0x4444ff, arrowSize * 2, arrowSize);
    STATE.scene.add(zArrow);

    // Draw θ (theta) arc - polar angle from Z-axis
    const thetaArcRadius = 3;
    const thetaArcGeometry = new THREE.BufferGeometry();
    const thetaPoints = [];
    for (let i = 0; i <= 20; i++) {
        const angle = (i / 20) * Math.PI / 2;  // 0 to 90°
        thetaPoints.push(new THREE.Vector3(
            thetaArcRadius * Math.sin(angle),
            0,
            thetaArcRadius * Math.cos(angle)
        ));
    }
    thetaArcGeometry.setFromPoints(thetaPoints);
    const thetaArcMaterial = new THREE.LineBasicMaterial({ color: 0xffff00, linewidth: 2 });
    const thetaArc = new THREE.Line(thetaArcGeometry, thetaArcMaterial);
    STATE.scene.add(thetaArc);

    // θ label
    const thetaLabel = createTextSprite('θ', '#ffff00', 40);
    thetaLabel.position.set(1.8, 0, 2.2);
    STATE.scene.add(thetaLabel);

    // θ annotation
    const thetaAnnotation = createTextSprite('(0°=face-on)', '#aaaa00', 24);
    thetaAnnotation.position.set(2.5, 0, 1.5);
    STATE.scene.add(thetaAnnotation);

    // Draw φ (phi) arc - azimuthal angle in X-Y plane
    const phiArcRadius = 4;
    const phiArcGeometry = new THREE.BufferGeometry();
    const phiPoints = [];
    for (let i = 0; i <= 20; i++) {
        const angle = (i / 20) * Math.PI / 2;  // 0 to 90°
        phiPoints.push(new THREE.Vector3(
            phiArcRadius * Math.cos(angle),
            phiArcRadius * Math.sin(angle),
            0
        ));
    }
    phiArcGeometry.setFromPoints(phiPoints);
    const phiArcMaterial = new THREE.LineBasicMaterial({ color: 0xff00ff, linewidth: 2 });
    const phiArc = new THREE.Line(phiArcGeometry, phiArcMaterial);
    STATE.scene.add(phiArc);

    // φ label
    const phiLabel = createTextSprite('φ', '#ff00ff', 40);
    phiLabel.position.set(3, 1.5, 0);
    STATE.scene.add(phiLabel);

    // φ annotation
    const phiAnnotation = createTextSprite('(from X)', '#aa00aa', 24);
    phiAnnotation.position.set(3.8, 2.2, 0);
    STATE.scene.add(phiAnnotation);

    // Edge-on indicator (in X-Y plane)
    const edgeOnLabel = createTextSprite('← Edge-on (i=90°)', '#88ffff', 28);
    edgeOnLabel.position.set(6, 0, 0);
    STATE.scene.add(edgeOnLabel);

    // Face-on indicator (along Z)
    const faceOnLabel = createTextSprite('↑ Face-on (i=0°)', '#88ffff', 28);
    faceOnLabel.position.set(0, 0, 7);
    STATE.scene.add(faceOnLabel);

    // === INCLINATION (i) - same as θ but with astronomy label ===
    // Draw inclination arc (orange) - angle from Z-axis (face-on direction)
    const incArcRadius = 5;
    const incArcGeometry = new THREE.BufferGeometry();
    const incPoints = [];
    for (let i = 0; i <= 30; i++) {
        const angle = (i / 30) * Math.PI / 2;  // 0 to 90°
        incPoints.push(new THREE.Vector3(
            0,
            incArcRadius * Math.sin(angle),
            incArcRadius * Math.cos(angle)
        ));
    }
    incArcGeometry.setFromPoints(incPoints);
    const incArcMaterial = new THREE.LineBasicMaterial({ color: 0xff8800, linewidth: 2 });
    const incArc = new THREE.Line(incArcGeometry, incArcMaterial);
    STATE.scene.add(incArc);

    // Inclination label
    const incLabel = createTextSprite('i (Inclination)', '#ff8800', 32);
    incLabel.position.set(0, 3.5, 3.5);
    STATE.scene.add(incLabel);

    // Inclination annotations
    const inc0Label = createTextSprite('i=0° (face-on)', '#ffaa44', 24);
    inc0Label.position.set(0, 0.5, 5.5);
    STATE.scene.add(inc0Label);

    const inc90Label = createTextSprite('i=90° (edge-on)', '#ffaa44', 24);
    inc90Label.position.set(0, 5.5, 0.5);
    STATE.scene.add(inc90Label);

    // === POSITION ANGLE (PA) - angle in sky plane from North to East ===
    // Draw PA arc (cyan) - in the X-Y plane showing N to E direction
    const paArcRadius = 6;
    const paArcGeometry = new THREE.BufferGeometry();
    const paPoints = [];
    for (let i = 0; i <= 30; i++) {
        const angle = (i / 30) * Math.PI / 2;  // 0 to 90° (N to E)
        paPoints.push(new THREE.Vector3(
            paArcRadius * Math.sin(angle),  // toward East (+X)
            paArcRadius * Math.cos(angle),  // from North (+Y)
            0
        ));
    }
    paArcGeometry.setFromPoints(paPoints);
    const paArcMaterial = new THREE.LineDashedMaterial({
        color: 0x00ffaa,
        linewidth: 2,
        dashSize: 0.3,
        gapSize: 0.15
    });
    const paArc = new THREE.Line(paArcGeometry, paArcMaterial);
    paArc.computeLineDistances();
    STATE.scene.add(paArc);

    // PA label
    const paLabel = createTextSprite('PA (Position Angle)', '#00ffaa', 32);
    paLabel.position.set(4, 5, 0);
    STATE.scene.add(paLabel);

    // PA direction indicators
    const paNorthLabel = createTextSprite('N (PA=0°)', '#00ffaa', 24);
    paNorthLabel.position.set(0, 6.8, 0);
    STATE.scene.add(paNorthLabel);

    const paEastLabel = createTextSprite('E (PA=90°)', '#00ffaa', 24);
    paEastLabel.position.set(6.8, 0, 0);
    STATE.scene.add(paEastLabel);

    // PA explanation
    const paExplainLabel = createTextSprite('PA: N→E in sky', '#008866', 22);
    paExplainLabel.position.set(5, 3, 0);
    STATE.scene.add(paExplainLabel);

    // Store labels in STATE for potential toggling
    STATE.coordinateLabels = {
        xLabel, yLabel, zLabel, originLabel, orbitalPlaneLabel,
        thetaLabel, thetaAnnotation, phiLabel, phiAnnotation,
        edgeOnLabel, faceOnLabel, thetaArc, phiArc,
        xArrow, yArrow, zArrow,
        // Inclination labels
        incArc, incLabel, inc0Label, inc90Label,
        // Position Angle labels
        paArc, paLabel, paNorthLabel, paEastLabel, paExplainLabel
    };

    console.log('Coordinate labels created for LOS convention visualization');
}

// Create a galaxy context visualization showing HVCC location
// Using SIMULATION SCALE (1 unit = 1 pc)
function createGalaxyContext() {
    STATE.galaxyGroup = new THREE.Group();
    STATE.galaxyGroup.visible = false;  // Hidden by default

    // ================================================================
    // SCALE: 1 unit = 1 pc (same as simulation)
    // HVCC is at l=-0.40°, b=-0.22°, ~60 pc from Sgr A*
    // The simulation origin IS the HVCC/IMBH location
    // ================================================================

    // HVCC galactic coordinates
    const l_deg = -0.40;
    const b_deg = -0.22;
    const r_proj = 60;  // pc from Sgr A*

    // Convert to Cartesian (in simulation scale, 1 unit = 1 pc)
    const l_rad = l_deg * Math.PI / 180;
    const b_rad = b_deg * Math.PI / 180;

    // Sgr A* position relative to simulation origin (HVCC)
    // The 60 pc is the total projected distance
    // HVCC is at l=-0.40°, b=-0.22° from Sgr A* (as seen from Earth)
    // Total angular separation ≈ sqrt(0.40² + 0.22²) = 0.456°
    //
    // X offset (East) = 60 pc * (l / total_angle) = 60 * (0.40/0.456) ≈ 52.6 pc
    // Y offset (North) = 60 pc * (b / total_angle) = 60 * (0.22/0.456) ≈ 28.9 pc
    //
    // HVCC is WEST and SOUTH of Sgr A*, so from HVCC, Sgr A* is EAST and NORTH
    const totalAngle = Math.sqrt(l_deg * l_deg + b_deg * b_deg);
    const sgrAX = r_proj * (-l_deg / totalAngle);  // Positive = East (Sgr A* is East of HVCC)
    const sgrAY = r_proj * (-b_deg / totalAngle);  // Positive = North (Sgr A* is North of HVCC)
    const sgrAZ = 0;

    // === GALAXY DISK (centered on Sgr A*, at galactic plane b=0°) ===
    // The galactic plane passes through Sgr A* at Y = sgrAY (since HVCC is below the plane)
    const diskRadius = 80;  // 80 pc visible region
    const diskGeometry = new THREE.RingGeometry(5, diskRadius, 64);
    const diskMaterial = new THREE.MeshBasicMaterial({
        color: 0x997744,  // Higher contrast warm color
        side: THREE.DoubleSide,
        transparent: true,
        opacity: 0.25
    });
    const galaxyDisk = new THREE.Mesh(diskGeometry, diskMaterial);
    galaxyDisk.rotation.x = Math.PI / 2;  // Flat in X-Z plane (galactic plane)
    galaxyDisk.position.set(sgrAX, sgrAY, sgrAZ);  // Centered on Sgr A*, at galactic plane height
    STATE.galaxyGroup.add(galaxyDisk);

    // Galactic plane boundary ring for visibility
    const diskBoundaryGeometry = new THREE.RingGeometry(diskRadius - 1, diskRadius, 64);
    const diskBoundaryMaterial = new THREE.MeshBasicMaterial({
        color: 0xddaa66,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: 0.6
    });
    const diskBoundary = new THREE.Mesh(diskBoundaryGeometry, diskBoundaryMaterial);
    diskBoundary.rotation.x = Math.PI / 2;
    diskBoundary.position.set(sgrAX, sgrAY, sgrAZ);
    STATE.galaxyGroup.add(diskBoundary);

    // Galactic plane label
    const galPlaneLabel = createTextSprite('Galactic Plane (b=0°)', '#ffcc88', 22);
    galPlaneLabel.position.set(sgrAX + 40, sgrAY + 3, 0);
    STATE.galaxyGroup.add(galPlaneLabel);

    // Spiral arms hint (around Sgr A*, in the galactic plane)
    const spiralMaterial = new THREE.LineBasicMaterial({
        color: 0xaa8855,
        transparent: true,
        opacity: 0.5,
        linewidth: 2
    });
    for (let arm = 0; arm < 4; arm++) {
        const spiralPoints = [];
        const armOffset = (arm * Math.PI / 2);
        for (let i = 0; i < 40; i++) {
            const t = i / 40;
            const r = 10 + t * 60;
            const theta = armOffset + t * Math.PI * 1.2;
            spiralPoints.push(new THREE.Vector3(
                sgrAX + r * Math.cos(theta),
                sgrAY,  // At galactic plane height
                sgrAZ + r * Math.sin(theta)
            ));
        }
        const spiralGeometry = new THREE.BufferGeometry().setFromPoints(spiralPoints);
        const spiral = new THREE.Line(spiralGeometry, spiralMaterial);
        STATE.galaxyGroup.add(spiral);
    }

    // Show HVCC offset from galactic plane - HIGH VISIBILITY
    const hvccOffsetTube = new THREE.TubeGeometry(
        new THREE.LineCurve3(
            new THREE.Vector3(0, 0, 0),
            new THREE.Vector3(0, sgrAY, 0)
        ), 1, 0.35, 8, false
    );
    const hvccOffsetMaterial = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.8
    });
    STATE.galaxyGroup.add(new THREE.Mesh(hvccOffsetTube, hvccOffsetMaterial));

    // Dashed line overlay
    const hvccOffsetLineGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(0, sgrAY, 0)
    ]);
    const hvccOffsetLine = new THREE.Line(hvccOffsetLineGeometry, new THREE.LineDashedMaterial({
        color: 0x00ffff,
        dashSize: 2,
        gapSize: 1,
        transparent: true,
        opacity: 0.6
    }));
    hvccOffsetLine.computeLineDistances();
    STATE.galaxyGroup.add(hvccOffsetLine);

    const hvccOffsetLabel = createTextSprite('b=-0.22°', '#00ffff', 22);
    hvccOffsetLabel.position.set(5, sgrAY / 2, 0);
    STATE.galaxyGroup.add(hvccOffsetLabel);

    // === SGR A* (at ~60 pc from simulation) ===
    const sgrAGeometry = new THREE.SphereGeometry(3, 32, 32);
    const sgrAMaterial = new THREE.MeshBasicMaterial({
        color: 0xffcc00,
        transparent: true,
        opacity: 0.95
    });
    const sgrA = new THREE.Mesh(sgrAGeometry, sgrAMaterial);
    sgrA.position.set(sgrAX, sgrAY, sgrAZ);
    STATE.galaxyGroup.add(sgrA);

    // Sgr A* glow
    const sgrAGlowGeometry = new THREE.SphereGeometry(5, 32, 32);
    const sgrAGlowMaterial = new THREE.MeshBasicMaterial({
        color: 0xffaa00,
        transparent: true,
        opacity: 0.3
    });
    const sgrAGlow = new THREE.Mesh(sgrAGlowGeometry, sgrAGlowMaterial);
    sgrAGlow.position.set(sgrAX, sgrAY, sgrAZ);
    STATE.galaxyGroup.add(sgrAGlow);

    const sgrALabel = createTextSprite('Sgr A* (4×10⁶ M☉)', '#ffcc00', 30);
    sgrALabel.position.set(sgrAX, sgrAY + 5, sgrAZ);
    STATE.galaxyGroup.add(sgrALabel);

    // === SIMULATION/HVCC LOCATION (at origin) ===
    // Circle marking the simulation area
    const simCircleGeometry = new THREE.RingGeometry(8, 10, 64);
    const simCircleMaterial = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: 0.3
    });
    const simCircle = new THREE.Mesh(simCircleGeometry, simCircleMaterial);
    STATE.galaxyGroup.add(simCircle);

    const hvccLabel = createTextSprite('HVCC CO-0.40-0.22', '#00ffff', 24);
    hvccLabel.position.set(0, 14, 0);
    STATE.galaxyGroup.add(hvccLabel);

    const hvccCoordLabel = createTextSprite('(l=-0.40°, b=-0.22°)', '#00aaaa', 18);
    hvccCoordLabel.position.set(0, 11, 0);
    STATE.galaxyGroup.add(hvccCoordLabel);

    const simLabel = createTextSprite('IMBH + Cloud Simulation', '#00ffff', 20);
    simLabel.position.set(0, -12, 0);
    STATE.galaxyGroup.add(simLabel);

    // === DISTANCE LINE from Sgr A* to HVCC ===
    const distLineGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(sgrAX, sgrAY, sgrAZ),
        new THREE.Vector3(0, 0, 0)
    ]);
    const distLineMaterial = new THREE.LineDashedMaterial({
        color: 0x00aaaa,
        dashSize: 2,
        gapSize: 1
    });
    const distLine = new THREE.Line(distLineGeometry, distLineMaterial);
    distLine.computeLineDistances();
    STATE.galaxyGroup.add(distLine);

    const distLabel = createTextSprite('~60 pc', '#00aaaa', 22);
    distLabel.position.set(sgrAX / 2, sgrAY / 2 + 3, 0);
    STATE.galaxyGroup.add(distLabel);

    // === GALACTIC COORDINATE AXES (at simulation scale) ===
    const axisLength = 50;  // 50 pc

    // X-axis: Galactic East
    const xArrow = new THREE.ArrowHelper(
        new THREE.Vector3(1, 0, 0), new THREE.Vector3(0, 0, 0),
        axisLength, 0xff6666, 2, 1
    );
    STATE.galaxyGroup.add(xArrow);
    const xLabel = createTextSprite('Gal. East (+l)', '#ff6666', 22);
    xLabel.position.set(axisLength + 5, 0, 0);
    STATE.galaxyGroup.add(xLabel);

    // Y-axis: Galactic North
    const yArrow = new THREE.ArrowHelper(
        new THREE.Vector3(0, 1, 0), new THREE.Vector3(0, 0, 0),
        axisLength, 0x66ff66, 2, 1
    );
    STATE.galaxyGroup.add(yArrow);
    const yLabel = createTextSprite('Gal. North (+b)', '#66ff66', 22);
    yLabel.position.set(0, axisLength + 5, 0);
    STATE.galaxyGroup.add(yLabel);

    // Z-axis: Toward Sun (LOS) - actual direction, not just +Z
    // Sun is at (sgrAX, sgrAY, 70), so direction from HVCC is normalized(sgrAX, sgrAY, 70)
    const sunDir = new THREE.Vector3(sgrAX, sgrAY, 70).normalize();
    const zArrow = new THREE.ArrowHelper(
        sunDir, new THREE.Vector3(0, 0, 0),
        axisLength, 0x6666ff, 2, 1
    );
    STATE.galaxyGroup.add(zArrow);
    const zLabelPos = sunDir.clone().multiplyScalar(axisLength + 5);
    const zLabel = createTextSprite('→ Sun (LOS)', '#6666ff', 22);
    zLabel.position.copy(zLabelPos);
    STATE.galaxyGroup.add(zLabel);

    // === SUN/EARTH INDICATOR ===
    // From HVCC's perspective:
    // - Sun is ~8 kpc away, in the galactic plane (b=0°)
    // - HVCC is at l=-0.40°, b=-0.22° from Sun
    // - So from HVCC, Sun appears at (Δl≈+0.40°, Δb≈+0.22°) - same direction as Sgr A*!
    // - Both Sun and Sgr A* are in the galactic plane at Y = sgrAY
    // - Sun is shown at scaled Z=70, but at correct (X, Y) position
    const sunZ = 70;  // Schematic distance (actual ~8000 pc)
    const sunGeometry = new THREE.SphereGeometry(3, 32, 32);
    const sunMaterial = new THREE.MeshBasicMaterial({
        color: 0xffff00,
        transparent: true,
        opacity: 0.9
    });
    const sunMarker = new THREE.Mesh(sunGeometry, sunMaterial);
    // Sun is in the galactic plane (Y = sgrAY) and nearly same direction as Sgr A*
    sunMarker.position.set(sgrAX, sgrAY, sunZ);
    STATE.galaxyGroup.add(sunMarker);

    // Sun glow effect
    const sunGlowGeometry = new THREE.SphereGeometry(5, 32, 32);
    const sunGlowMaterial = new THREE.MeshBasicMaterial({
        color: 0xffff44,
        transparent: true,
        opacity: 0.2
    });
    const sunGlow = new THREE.Mesh(sunGlowGeometry, sunGlowMaterial);
    sunGlow.position.set(sgrAX, sgrAY, sunZ);
    STATE.galaxyGroup.add(sunGlow);

    const sunLabel = createTextSprite('☉ Sun / Earth (b=0°)', '#ffff00', 24);
    sunLabel.position.set(sgrAX, sgrAY + 8, sunZ);
    STATE.galaxyGroup.add(sunLabel);

    const sunDistLabel = createTextSprite('~8 kpc (schematic Z)', '#aaaa00', 18);
    sunDistLabel.position.set(sgrAX, sgrAY + 5, sunZ);
    STATE.galaxyGroup.add(sunDistLabel);

    // LOS line from HVCC toward Sun (actual direction, not just +Z)
    const losLineGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(sgrAX, sgrAY, sunZ - 5)
    ]);
    const losLineMaterial = new THREE.LineDashedMaterial({
        color: 0xffff44,
        dashSize: 3,
        gapSize: 1.5,
        transparent: true,
        opacity: 0.5
    });
    const losLine = new THREE.Line(losLineGeometry, losLineMaterial);
    losLine.computeLineDistances();
    STATE.galaxyGroup.add(losLine);

    const losLabel = createTextSprite('Line of Sight to Sun', '#ffff44', 18);
    losLabel.position.set(sgrAX / 2, sgrAY / 2 + 5, sunZ / 2);
    STATE.galaxyGroup.add(losLabel);

    // Line from Sgr A* to Sun (both in galactic plane, along LOS direction)
    // Now that Sun is at (sgrAX, sgrAY, sunZ), this is a line along Z
    const gcToSunGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(sgrAX, sgrAY, 0),    // Sgr A*
        new THREE.Vector3(sgrAX, sgrAY, sunZ)  // Sun (same X,Y, different Z)
    ]);
    const gcToSunMaterial = new THREE.LineDashedMaterial({
        color: 0xffaa00,
        dashSize: 4,
        gapSize: 2,
        transparent: true,
        opacity: 0.3
    });
    const gcToSunLine = new THREE.Line(gcToSunGeometry, gcToSunMaterial);
    gcToSunLine.computeLineDistances();
    STATE.galaxyGroup.add(gcToSunLine);

    // Label explaining that Sgr A* and Sun are nearly collinear from HVCC
    const collinearLabel = createTextSprite('Sgr A* → Sun (~8 kpc)', '#ffaa00', 16);
    collinearLabel.position.set(sgrAX + 8, sgrAY, sunZ / 2);
    STATE.galaxyGroup.add(collinearLabel);

    // === SCALE BARS ===
    // 10 pc scale bar
    const scale10Geometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-40, -25, 0),
        new THREE.Vector3(-30, -25, 0)
    ]);
    const scaleMaterial = new THREE.LineBasicMaterial({ color: 0xaaaaaa });
    const scale10Bar = new THREE.Line(scale10Geometry, scaleMaterial);
    STATE.galaxyGroup.add(scale10Bar);

    const scale10Label = createTextSprite('10 pc', '#aaaaaa', 18);
    scale10Label.position.set(-35, -22, 0);
    STATE.galaxyGroup.add(scale10Label);

    // 50 pc scale bar
    const scale50Geometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-40, -30, 0),
        new THREE.Vector3(10, -30, 0)
    ]);
    const scale50Bar = new THREE.Line(scale50Geometry, scaleMaterial);
    STATE.galaxyGroup.add(scale50Bar);

    const scale50Label = createTextSprite('50 pc', '#aaaaaa', 18);
    scale50Label.position.set(-15, -27, 0);
    STATE.galaxyGroup.add(scale50Label);

    // === GALACTIC l, b COORDINATE GRID ===
    // Angular scale at GC distance (~8 kpc): 1° ≈ 140 pc, 0.1° ≈ 14 pc
    const degToPC = sgrAX / Math.abs(l_deg);  // ~132 pc per degree

    // High contrast colors for coordinate grid
    const gridAxisColor = 0xffdd99;  // Bright warm yellow
    const gridTickColor = 0xffcc77;  // Warm orange-yellow
    const gridLabelColor = '#ffeeaa';  // Light yellow for labels

    // l coordinate axis (along galactic plane at b=0°) - THICK LINE
    const lAxisPoints = [];
    for (let x = sgrAX - 70; x <= sgrAX + 70; x += 0.5) {
        lAxisPoints.push(new THREE.Vector3(x, sgrAY, 0));
    }
    const lAxisGeometry = new THREE.BufferGeometry().setFromPoints(lAxisPoints);
    const coordAxisMaterial = new THREE.LineBasicMaterial({
        color: gridAxisColor,
        transparent: true,
        opacity: 0.9,
        linewidth: 3
    });
    const lAxis = new THREE.Line(lAxisGeometry, coordAxisMaterial);
    STATE.galaxyGroup.add(lAxis);

    // Create thick line using tube geometry for l-axis
    const lAxisTubeGeom = new THREE.TubeGeometry(
        new THREE.LineCurve3(
            new THREE.Vector3(sgrAX - 70, sgrAY, 0),
            new THREE.Vector3(sgrAX + 70, sgrAY, 0)
        ), 1, 0.3, 8, false
    );
    const lAxisTube = new THREE.Mesh(lAxisTubeGeom, new THREE.MeshBasicMaterial({
        color: gridAxisColor,
        transparent: true,
        opacity: 0.8
    }));
    STATE.galaxyGroup.add(lAxisTube);

    // l tick marks and labels at b=0° (galactic plane)
    const lValues = [-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3];
    const tickMaterial = new THREE.MeshBasicMaterial({ color: gridTickColor });

    for (const lVal of lValues) {
        // Position: lVal=0 is at Sgr A* (sgrAX), HVCC is at l=-0.40°
        const xPos = sgrAX + lVal * degToPC;
        const tickLen = (lVal === 0 || Math.abs(lVal - l_deg) < 0.01) ? 5 : 2.5;
        const tickWidth = (lVal === 0 || Math.abs(lVal - l_deg) < 0.01) ? 0.4 : 0.25;

        // Tick mark using cylinder for thickness
        const tickGeom = new THREE.CylinderGeometry(tickWidth, tickWidth, tickLen * 2, 8);
        const tick = new THREE.Mesh(tickGeom, tickMaterial);
        tick.position.set(xPos, sgrAY, 0);
        STATE.galaxyGroup.add(tick);

        // Label (only major ticks)
        if (lVal === 0 || lVal === -0.4 || Math.abs(lVal) === 0.2 || Math.abs(lVal) === 0.5) {
            const labelText = lVal === 0 ? 'l=0°' : `l=${lVal > 0 ? '+' : ''}${lVal}°`;
            const labelColor = lVal === 0 ? '#ffcc00' : (Math.abs(lVal - l_deg) < 0.01 ? '#00ffff' : gridLabelColor);
            const fontSize = (lVal === 0 || Math.abs(lVal - l_deg) < 0.01) ? 22 : 18;
            const lLabel = createTextSprite(labelText, labelColor, fontSize);
            lLabel.position.set(xPos, sgrAY + 8, 0);
            STATE.galaxyGroup.add(lLabel);
        }
    }

    // l-axis label
    const lAxisLabel = createTextSprite('← −l (West)          Galactic Longitude          +l (East) →', gridLabelColor, 20);
    lAxisLabel.position.set(sgrAX, sgrAY + 14, 0);
    STATE.galaxyGroup.add(lAxisLabel);

    // b coordinate axis (perpendicular to galactic plane at l=0°, through Sgr A*) - THICK LINE
    const bAxisTubeGeom = new THREE.TubeGeometry(
        new THREE.LineCurve3(
            new THREE.Vector3(sgrAX, sgrAY - 40, 0),
            new THREE.Vector3(sgrAX, sgrAY + 40, 0)
        ), 1, 0.3, 8, false
    );
    const bAxisTube = new THREE.Mesh(bAxisTubeGeom, new THREE.MeshBasicMaterial({
        color: gridAxisColor,
        transparent: true,
        opacity: 0.8
    }));
    STATE.galaxyGroup.add(bAxisTube);

    // b tick marks and labels at l=0° (through Sgr A*)
    const bValues = [-0.3, -0.22, -0.2, -0.1, 0, 0.1, 0.2, 0.3];

    for (const bVal of bValues) {
        // Position: bVal=0 is at galactic plane (sgrAY), HVCC is at b=-0.22°
        const yPos = sgrAY + bVal * degToPC;
        const tickLen = (bVal === 0 || Math.abs(bVal - b_deg) < 0.01) ? 5 : 2.5;
        const tickWidth = (bVal === 0 || Math.abs(bVal - b_deg) < 0.01) ? 0.4 : 0.25;

        // Tick mark using cylinder (rotated for horizontal)
        const tickGeom = new THREE.CylinderGeometry(tickWidth, tickWidth, tickLen * 2, 8);
        const tick = new THREE.Mesh(tickGeom, tickMaterial);
        tick.rotation.z = Math.PI / 2;  // Rotate to horizontal
        tick.position.set(sgrAX, yPos, 0);
        STATE.galaxyGroup.add(tick);

        // Label
        if (bVal === 0 || Math.abs(bVal - b_deg) < 0.01 || Math.abs(bVal) === 0.2 || Math.abs(bVal) === 0.3) {
            const labelText = bVal === 0 ? 'b=0°' : `b=${bVal > 0 ? '+' : ''}${bVal}°`;
            const labelColor = bVal === 0 ? '#ffcc00' : (Math.abs(bVal - b_deg) < 0.01 ? '#00ffff' : gridLabelColor);
            const fontSize = (bVal === 0 || Math.abs(bVal - b_deg) < 0.01) ? 22 : 18;
            const bLabel = createTextSprite(labelText, labelColor, fontSize);
            bLabel.position.set(sgrAX + 10, yPos, 0);
            STATE.galaxyGroup.add(bLabel);
        }
    }

    // b-axis label (vertical text) - high contrast
    const bAxisLabel = createTextSprite('+b (North) ↑', gridLabelColor, 20);
    bAxisLabel.position.set(sgrAX + 14, sgrAY + 35, 0);
    STATE.galaxyGroup.add(bAxisLabel);

    const bAxisLabelNeg = createTextSprite('↓ −b (South)', gridLabelColor, 20);
    bAxisLabelNeg.position.set(sgrAX + 14, sgrAY - 35, 0);
    STATE.galaxyGroup.add(bAxisLabelNeg);

    // Mark HVCC position on coordinate grid - HIGH VISIBILITY
    // Vertical line from HVCC to galactic plane
    const hvccVerticalTube = new THREE.TubeGeometry(
        new THREE.LineCurve3(
            new THREE.Vector3(0, 0, 0),
            new THREE.Vector3(0, sgrAY, 0)
        ), 1, 0.2, 8, false
    );
    const hvccGridMaterial = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.7
    });
    STATE.galaxyGroup.add(new THREE.Mesh(hvccVerticalTube, hvccGridMaterial));

    // Horizontal line from HVCC projection to l=0° axis
    const hvccHorizontalTube = new THREE.TubeGeometry(
        new THREE.LineCurve3(
            new THREE.Vector3(0, sgrAY, 0),
            new THREE.Vector3(sgrAX, sgrAY, 0)
        ), 1, 0.2, 8, false
    );
    STATE.galaxyGroup.add(new THREE.Mesh(hvccHorizontalTube, hvccGridMaterial));

    // Add dashed line version for additional visibility
    const hvccGridMarker = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(0, sgrAY, 0),
        new THREE.Vector3(sgrAX, sgrAY, 0)
    ]);
    const hvccGridLine = new THREE.Line(hvccGridMarker, new THREE.LineDashedMaterial({
        color: 0x00ffff,
        dashSize: 2,
        gapSize: 1,
        transparent: true,
        opacity: 0.5
    }));
    hvccGridLine.computeLineDistances();
    STATE.galaxyGroup.add(hvccGridLine);

    // === INFO TEXT ===
    const infoLabel = createTextSprite('HVCC-centered view | Scale: 1 unit = 1 pc', '#666666', 18);
    infoLabel.position.set(0, -35, 0);
    STATE.galaxyGroup.add(infoLabel);

    const infoLabel2 = createTextSprite('Galactic plane (b=0°) contains both Sgr A* and Sun', '#666666', 16);
    infoLabel2.position.set(0, -39, 0);
    STATE.galaxyGroup.add(infoLabel2);

    const infoLabel3 = createTextSprite('HVCC is at b=-0.22° (below galactic plane)', '#666666', 16);
    infoLabel3.position.set(0, -43, 0);
    STATE.galaxyGroup.add(infoLabel3);

    // Galaxy group centered on simulation origin
    STATE.galaxyGroup.position.set(0, 0, 0);

    STATE.scene.add(STATE.galaxyGroup);

    console.log('Galaxy context visualization created (simulation scale: 1 unit = 1 pc)');
    console.log(`Simulation (HVCC) at origin, Sgr A* at (${sgrAX.toFixed(1)}, ${sgrAY.toFixed(1)}, 0) pc`);
    console.log(`Sun at (${sgrAX.toFixed(1)}, ${sgrAY.toFixed(1)}, 70) - same direction as Sgr A*, ~8 kpc away`);
}

// Toggle galaxy context visibility
function toggleGalaxyContext() {
    if (!STATE.galaxyGroup) return;
    STATE.galaxyGroup.visible = !STATE.galaxyGroup.visible;
    console.log('Galaxy context:', STATE.galaxyGroup.visible ? 'visible' : 'hidden');
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
    // Add to orbital plane group (so it rotates with the orbital pole)
    if (STATE.orbitalPlaneGroup) {
        STATE.orbitalPlaneGroup.add(STATE.orbitLine);
    } else {
        STATE.scene.add(STATE.orbitLine);
    }

    // Pericenter position
    const periPos = new THREE.Vector3(
        CONFIG.r_peri * e_hat[0],
        CONFIG.r_peri * e_hat[1],
        CONFIG.r_peri * e_hat[2]
    );

    // Pericenter marker - add to orbital plane group so it rotates with orbit
    if (STATE.periMarker) {
        if (STATE.orbitalPlaneGroup) {
            STATE.orbitalPlaneGroup.remove(STATE.periMarker);
        } else {
            STATE.scene.remove(STATE.periMarker);
        }
        STATE.periMarker.geometry.dispose();
        STATE.periMarker.material.dispose();
    }
    const periGeom = new THREE.RingGeometry(0.08, 0.12, 32);
    const periMat = new THREE.MeshBasicMaterial({ color: 0x00ff88, side: THREE.DoubleSide });
    STATE.periMarker = new THREE.Mesh(periGeom, periMat);
    STATE.periMarker.position.copy(periPos);
    STATE.periMarker.lookAt(0, 0, 0);
    if (STATE.orbitalPlaneGroup) {
        STATE.orbitalPlaneGroup.add(STATE.periMarker);
    } else {
        STATE.scene.add(STATE.periMarker);
    }

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
    // Create a group for orbital pole and fixed LOS visualization
    STATE.losArrows = new THREE.Group();

    // === DRAGGABLE ORBITAL POLE SPHERE (cyan) ===
    // This represents the direction of orbital angular momentum
    const sphereGeom = new THREE.SphereGeometry(0.5, 16, 16);
    const sphereMat = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.9
    });
    STATE.losSphere = new THREE.Mesh(sphereGeom, sphereMat);
    STATE.losSphere.name = 'losSphere';

    // Position on sphere in direction of orbital pole (losVector)
    const radius = 5;
    STATE.losSphere.position.copy(STATE.losVector.clone().multiplyScalar(radius));
    STATE.scene.add(STATE.losSphere);

    // === FIXED LOS TO SUN INDICATOR (yellow) ===
    // This shows the actual observation direction from Earth - user cannot change this
    STATE.fixedLOSGroup = new THREE.Group();

    // Compute fixed LOS direction (toward Sun from HVCC)
    // Using the same geometry as galaxy context: Sun is at (sgrAX, sgrAY, ~70)
    const l_deg = -0.40;
    const b_deg = -0.22;
    const r_proj = 60;
    const totalAngle = Math.sqrt(l_deg * l_deg + b_deg * b_deg);
    const sgrAX = r_proj * (-l_deg / totalAngle);
    const sgrAY = r_proj * (-b_deg / totalAngle);
    const sunDir = new THREE.Vector3(sgrAX, sgrAY, 70).normalize();

    // Store the fixed LOS direction
    STATE.fixedLOStoSun.copy(sunDir);

    // Yellow arrow showing fixed LOS to Sun
    const fixedLOSArrow = new THREE.ArrowHelper(
        sunDir,
        new THREE.Vector3(0, 0, 0),
        4,
        0xffff00,
        0.6,
        0.3
    );
    STATE.fixedLOSGroup.add(fixedLOSArrow);

    // Small yellow sphere at end of fixed LOS arrow
    const fixedLOSSphereGeom = new THREE.SphereGeometry(0.25, 16, 16);
    const fixedLOSSphereMat = new THREE.MeshBasicMaterial({
        color: 0xffff00,
        transparent: true,
        opacity: 0.7
    });
    const fixedLOSSphere = new THREE.Mesh(fixedLOSSphereGeom, fixedLOSSphereMat);
    fixedLOSSphere.position.copy(sunDir.clone().multiplyScalar(4.5));
    STATE.fixedLOSGroup.add(fixedLOSSphere);

    STATE.scene.add(STATE.fixedLOSGroup);

    // Arrow showing orbital pole direction (cyan)
    updateLOSArrows();

    STATE.scene.add(STATE.losArrows);
}

function updateLOSArrows() {
    // Remove old arrows from group
    while (STATE.losArrows.children.length > 0) {
        STATE.losArrows.remove(STATE.losArrows.children[0]);
    }

    // Arrow showing orbital pole direction (cyan)
    const arrowLength = 4.5;
    const orbitalPoleArrow = new THREE.ArrowHelper(
        STATE.losVector.clone(),  // losVector = orbital pole direction
        new THREE.Vector3(0, 0, 0),
        arrowLength,
        0x00ffff,
        0.5,
        0.25
    );
    STATE.losArrows.add(orbitalPoleArrow);

    // Add a small orbital plane indicator (disk perpendicular to orbital pole)
    const planeGeom = new THREE.RingGeometry(1.5, 2, 32);
    const planeMat = new THREE.MeshBasicMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.2,
        side: THREE.DoubleSide
    });
    const planeMesh = new THREE.Mesh(planeGeom, planeMat);
    // Orient plane perpendicular to orbital pole
    planeMesh.lookAt(STATE.losVector);
    STATE.losArrows.add(planeMesh);

    // Update sphere position (sphere is in scene, not group)
    if (STATE.losSphere) {
        STATE.losSphere.position.copy(STATE.losVector.clone().multiplyScalar(5));
    }

    // Update the orbital plane group rotation to match the new orbital pole
    updateOrbitalPlaneRotation();

    // Update orbital pole info display
    updateLOSInfo();
}

function updateLOSInfo() {
    // v = orbital pole direction (controlled by user via cyan sphere)
    const v = STATE.losVector;
    // n_sun = fixed LOS direction toward Sun (constant)
    const n_sun = STATE.fixedLOStoSun;

    // Spherical coordinates of orbital pole (physics convention)
    // θ (theta): polar angle from +z axis
    // φ (phi): azimuthal angle from +x axis
    const theta = Math.acos(Math.max(-1, Math.min(1, v.z))) * 180 / Math.PI;
    let phi = Math.atan2(v.y, v.x) * 180 / Math.PI;
    if (phi < 0) phi += 360;

    // === INCLINATION (relative to fixed LOS to Sun) ===
    // i = angle between orbital pole and the fixed LOS direction
    // i = 0°: face-on (orbital pole aligned with Earth LOS, we look down the pole)
    // i = 90°: edge-on (orbital pole perpendicular to Earth LOS, we look in the orbital plane)
    const dotProduct = v.x * n_sun.x + v.y * n_sun.y + v.z * n_sun.z;
    const inclination = Math.acos(Math.max(-1, Math.min(1, dotProduct))) * 180 / Math.PI;

    // === POSITION ANGLE (orientation of orbital plane in the sky) ===
    // PA is the angle of the line of nodes (intersection of orbital plane with sky plane)
    // measured from North (galactic +b direction) toward East (galactic +l direction)
    //
    // The line of nodes is perpendicular to both the orbital pole and the LOS
    // lineOfNodes = orbitalPole × LOS (cross product)
    const lineOfNodes = new THREE.Vector3().crossVectors(v, n_sun).normalize();

    // Project onto sky plane (perpendicular to LOS)
    // In the sky plane, define North as the projection of +Y (galactic North)
    // and East as projection of +X (galactic East)
    // PA = angle from projected North to line of nodes, measured toward East

    // For simplicity, compute PA from the azimuthal angle of the line of nodes
    // when viewed along the LOS direction
    let pa = Math.atan2(lineOfNodes.x, lineOfNodes.y) * 180 / Math.PI;
    if (pa < 0) pa += 360;

    // Update dropdown to show preset if matching
    const select = document.getElementById('los-direction');
    if (select) {
        if (Math.abs(v.x - 1) < 0.01 && Math.abs(v.y) < 0.01 && Math.abs(v.z) < 0.01) {
            select.value = 'x';
        } else if (Math.abs(v.x) < 0.01 && Math.abs(v.y - 1) < 0.01 && Math.abs(v.z) < 0.01) {
            select.value = 'y';
        } else if (Math.abs(v.x) < 0.01 && Math.abs(v.y) < 0.01 && Math.abs(v.z - 1) < 0.01) {
            select.value = 'z';
        }
    }

    // Update display
    const vectorEl = document.getElementById('los-vector');
    const thetaEl = document.getElementById('los-theta');
    const phiEl = document.getElementById('los-phi');
    const inclinationEl = document.getElementById('los-inclination');
    const paEl = document.getElementById('los-pa');

    if (vectorEl) vectorEl.textContent = `(${v.x.toFixed(2)}, ${v.y.toFixed(2)}, ${v.z.toFixed(2)})`;
    if (thetaEl) thetaEl.textContent = `${theta.toFixed(1)}°`;
    if (phiEl) phiEl.textContent = `${phi.toFixed(1)}°`;
    if (inclinationEl) inclinationEl.textContent = `${inclination.toFixed(1)}°`;
    if (paEl) paEl.textContent = `${pa.toFixed(1)}°`;
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

// VTK-style Jet colormap: dark blue -> blue -> cyan -> green -> yellow -> orange -> red
function valueToColor(t) {
    // Clamp t to [0, 1]
    t = Math.max(0, Math.min(1, t));

    // Jet colormap control points
    const stops = [
        { t: 0.00, r: 0.0,  g: 0.0,  b: 0.5  },  // Dark blue
        { t: 0.15, r: 0.0,  g: 0.0,  b: 1.0  },  // Blue
        { t: 0.35, r: 0.0,  g: 1.0,  b: 1.0  },  // Cyan
        { t: 0.50, r: 0.0,  g: 1.0,  b: 0.0  },  // Green
        { t: 0.65, r: 1.0,  g: 1.0,  b: 0.0  },  // Yellow
        { t: 0.85, r: 1.0,  g: 0.5,  b: 0.0  },  // Orange
        { t: 1.00, r: 1.0,  g: 0.0,  b: 0.0  }   // Red
    ];

    // Find the two stops to interpolate between
    let i = 0;
    while (i < stops.length - 1 && stops[i + 1].t < t) i++;

    const s0 = stops[i];
    const s1 = stops[Math.min(i + 1, stops.length - 1)];

    // Linear interpolation
    const dt = s1.t - s0.t;
    const f = dt > 0 ? (t - s0.t) / dt : 0;

    const r = s0.r + f * (s1.r - s0.r);
    const g = s0.g + f * (s1.g - s0.g);
    const b = s0.b + f * (s1.b - s0.b);

    return new THREE.Color(r, g, b);
}

function updateVisualization(frameIndex) {
    if (frameIndex < 0 || frameIndex >= STATE.snapshots.length) return;

    STATE.currentFrame = frameIndex;
    const data = STATE.snapshots[frameIndex];

    // Remove old particle system from orbital plane group
    if (STATE.particleSystem) {
        if (STATE.orbitalPlaneGroup) {
            STATE.orbitalPlaneGroup.remove(STATE.particleSystem);
        } else {
            STATE.scene.remove(STATE.particleSystem);
        }
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

    // Create set of selected indices for fast lookup
    const selectedSet = new Set(STATE.selectedParticles || []);

    // Second pass: set colors using GLOBAL ranges
    data.particles.forEach((p, i) => {
        positions[i*3] = p.x;
        positions[i*3+1] = p.y;
        positions[i*3+2] = p.z;

        let color;
        const isSelected = selectedSet.has(i);

        if (p.is_ghost === 1) {
            color = new THREE.Color(0.4, 0.4, 0.5);
            sizes[i] = 0.02;
        } else if (isSelected) {
            // Highlight selected particles in orange
            color = new THREE.Color(1.0, 0.53, 0.4);  // #ff8866
            sizes[i] = 0.06;  // Larger size
        } else {
            const t = getColorForValue(p, comVelocity);
            color = valueToColor(t);
            // Dim non-selected when selection exists
            if (selectedSet.size > 0) {
                color.multiplyScalar(0.4);
                sizes[i] = 0.02;
            } else {
                sizes[i] = 0.03 + 0.02 * t;
            }
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
    // Add to orbital plane group (so particles rotate with the orbital pole)
    if (STATE.orbitalPlaneGroup) {
        STATE.orbitalPlaneGroup.add(STATE.particleSystem);
    } else {
        STATE.scene.add(STATE.particleSystem);
    }

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
        case 'mach':
            title = 'Mach';
            minLabel = range.min.toFixed(1);
            maxLabel = range.max.toFixed(1);
            break;
        case 'vrel':
            title = '|v-v_COM| (km/s)';
            minLabel = range.min.toFixed(1);
            maxLabel = range.max.toFixed(1);
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

    // VTK-style Jet colorbar matching valueToColor function
    // Gradient from top (max=1.0) to bottom (min=0.0)
    document.getElementById('colorbar').style.background =
        'linear-gradient(to bottom, ' +
        '#FF0000 0%, ' +    // Red (t=1.00)
        '#FF7F00 15%, ' +   // Orange (t=0.85)
        '#FFFF00 35%, ' +   // Yellow (t=0.65)
        '#00FF00 50%, ' +   // Green (t=0.50)
        '#00FFFF 65%, ' +   // Cyan (t=0.35)
        '#0000FF 85%, ' +   // Blue (t=0.15)
        '#00007F 100%)';
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
