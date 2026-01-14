// ============================================================
// IMBH-Cloud Visualization - Scene Initialization
// Core THREE.js setup: scene, camera, renderer, controls
// ============================================================

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

    // Add pointer event listeners for LOS dragging
    STATE.renderer.domElement.addEventListener('pointerdown', onPointerDown);
    window.addEventListener('pointermove', onPointerMove);
    window.addEventListener('pointerup', onPointerUp);

    // Add hover effect
    STATE.renderer.domElement.addEventListener('pointermove', onPointerHover);
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
