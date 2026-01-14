// ============================================================
// IMBH-Cloud Visualization - LOS/Orbital Pole Control
// Draggable orbital pole sphere and interaction handlers
// ============================================================

// Module-level variables for raycasting
let raycaster, mouse;

function initLOSRaycaster() {
    raycaster = new THREE.Raycaster();
    mouse = new THREE.Vector2();
}

function createLOSControl() {
    // Initialize raycaster if not done
    if (!raycaster) initLOSRaycaster();

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

    // === SKY PLANE (perpendicular to LOS) ===
    // This is the plane where PA is measured
    const skyPlaneRadius = 6;
    const skyPlaneGeom = new THREE.RingGeometry(2, skyPlaneRadius, 64);
    const skyPlaneMat = new THREE.MeshBasicMaterial({
        color: 0xffff00,
        transparent: true,
        opacity: 0.15,
        side: THREE.DoubleSide
    });
    const skyPlaneMesh = new THREE.Mesh(skyPlaneGeom, skyPlaneMat);
    // Orient perpendicular to LOS (plane normal = sunDir)
    skyPlaneMesh.lookAt(sunDir);
    STATE.fixedLOSGroup.add(skyPlaneMesh);

    // Sky plane border circle
    const skyPlaneBorderGeom = new THREE.BufferGeometry();
    const borderPoints = [];
    for (let i = 0; i <= 64; i++) {
        const angle = (i / 64) * Math.PI * 2;
        borderPoints.push(new THREE.Vector3(
            skyPlaneRadius * Math.cos(angle),
            skyPlaneRadius * Math.sin(angle),
            0
        ));
    }
    skyPlaneBorderGeom.setFromPoints(borderPoints);
    const skyPlaneBorder = new THREE.Line(
        skyPlaneBorderGeom,
        new THREE.LineBasicMaterial({ color: 0xffff00, transparent: true, opacity: 0.5 })
    );
    skyPlaneBorder.lookAt(sunDir);
    STATE.fixedLOSGroup.add(skyPlaneBorder);

    // === NORTH AND EAST IN SKY PLANE ===
    // Project galactic North (+Y) onto sky plane
    const northRef = new THREE.Vector3(0, 1, 0);  // Galactic +b
    const northInSky = northRef.clone().sub(
        sunDir.clone().multiplyScalar(northRef.dot(sunDir))
    ).normalize();

    // East in sky plane (perpendicular to both LOS and North)
    const eastInSky = new THREE.Vector3().crossVectors(sunDir, northInSky).normalize();

    // North arrow in sky plane (green)
    const northArrowLength = skyPlaneRadius * 0.9;
    const northArrow = new THREE.ArrowHelper(
        northInSky,
        new THREE.Vector3(0, 0, 0),
        northArrowLength,
        0x88ff88,  // Light green
        0.4,
        0.2
    );
    STATE.fixedLOSGroup.add(northArrow);

    // East arrow in sky plane (red)
    const eastArrow = new THREE.ArrowHelper(
        eastInSky,
        new THREE.Vector3(0, 0, 0),
        northArrowLength,
        0xff8888,  // Light red
        0.4,
        0.2
    );
    STATE.fixedLOSGroup.add(eastArrow);

    // Labels for N and E at arrow tips
    const northLabel = createTextSprite('N', '#88ff88', 32);
    northLabel.position.copy(northInSky.clone().multiplyScalar(northArrowLength + 0.8));
    STATE.fixedLOSGroup.add(northLabel);

    const eastLabel = createTextSprite('E', '#ff8888', 32);
    eastLabel.position.copy(eastInSky.clone().multiplyScalar(northArrowLength + 0.8));
    STATE.fixedLOSGroup.add(eastLabel);

    // Label for LOS direction
    const losLabel = createTextSprite('LOS→Sun', '#ffff00', 24);
    losLabel.position.copy(sunDir.clone().multiplyScalar(5.5));
    STATE.fixedLOSGroup.add(losLabel);

    // Store directions for reference
    STATE.skyPlaneLabels = { northInSky, eastInSky };

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
    // ================================================================
    // ORBITAL ORIENTATION ANGLES
    // Following standard astronomical conventions from:
    // - Ballone et al. 2018 MNRAS 480, 4684 (CO-0.40-0.22 simulation)
    // - orbitize documentation (orbital element definitions)
    // ================================================================

    // L = orbital pole direction (orbital angular momentum, controlled by user)
    const L = STATE.losVector;
    // n = fixed Line-of-Sight direction toward Sun/Earth (constant)
    const n = STATE.fixedLOStoSun;

    // Spherical coordinates of orbital pole in simulation frame
    // θ (theta): polar angle from +z axis (colatitude)
    // φ (phi): azimuthal angle from +x axis in x-y plane
    const theta = Math.acos(Math.max(-1, Math.min(1, L.z))) * 180 / Math.PI;
    let phi = Math.atan2(L.y, L.x) * 180 / Math.PI;
    if (phi < 0) phi += 360;

    // ================================================================
    // INCLINATION (i) - Standard astronomical definition
    // ================================================================
    // i = angle between orbital angular momentum vector (L) and line of sight (n)
    //
    // Definition: cos(i) = L · n / (|L| |n|)
    //
    // Physical meaning:
    //   i = 0°:   Face-on view (looking down the orbital pole, orbit appears circular)
    //   i = 90°:  Edge-on view (orbital pole ⊥ to LOS, orbit appears as a line)
    //   i = 180°: Face-on view (looking up at orbital pole, opposite sense)
    //
    // Range: 0° to 180° (determines if orbit appears prograde or retrograde)
    //   0° < i < 90°:  Prograde (counterclockwise when viewed from +LOS)
    //   90° < i < 180°: Retrograde (clockwise when viewed from +LOS)
    // ================================================================
    const cosInc = L.x * n.x + L.y * n.y + L.z * n.z;
    const inclination = Math.acos(Math.max(-1, Math.min(1, cosInc))) * 180 / Math.PI;

    // ================================================================
    // POSITION ANGLE OF ASCENDING NODE (Ω or PA)
    // ================================================================
    // The ascending node is where the orbit crosses the sky plane
    // moving AWAY from the observer (receding).
    //
    // PA is measured in the sky plane (perpendicular to LOS):
    //   - From North (galactic +b direction, +Y in our coords)
    //   - Toward East (galactic +l direction, +X in our coords)
    //   - Range: 0° to 360°
    //
    // Calculation:
    // 1. Line of nodes = L × n (intersection of orbital plane and sky plane)
    // 2. Project reference North onto sky plane
    // 3. Measure angle from projected North to line of nodes
    // ================================================================

    // Line of nodes vector (points toward ascending node)
    // This is the intersection of orbital plane with sky plane
    const lineOfNodes = new THREE.Vector3().crossVectors(L, n);

    // Handle edge case: if orbital pole is parallel to LOS (face-on),
    // line of nodes is undefined (PA is degenerate)
    let pa = 0;
    if (lineOfNodes.length() > 0.001) {
        lineOfNodes.normalize();

        // Reference directions in sky plane:
        // North = +Y (galactic +b), East = +X (galactic +l)
        //
        // Project North onto sky plane (remove component along LOS)
        const northRef = new THREE.Vector3(0, 1, 0);  // Galactic North
        const northInSky = northRef.clone().sub(
            n.clone().multiplyScalar(northRef.dot(n))
        );

        // If North is nearly parallel to LOS, use alternative reference
        if (northInSky.length() < 0.001) {
            // Fall back to East direction
            const eastRef = new THREE.Vector3(1, 0, 0);
            northInSky.copy(eastRef.clone().sub(
                n.clone().multiplyScalar(eastRef.dot(n))
            ));
        }
        northInSky.normalize();

        // East in sky plane (perpendicular to both LOS and projected North)
        const eastInSky = new THREE.Vector3().crossVectors(n, northInSky).normalize();

        // PA = angle from North to line of nodes, measured toward East
        // Using atan2 for proper quadrant handling
        const northComponent = lineOfNodes.dot(northInSky);
        const eastComponent = lineOfNodes.dot(eastInSky);
        pa = Math.atan2(eastComponent, northComponent) * 180 / Math.PI;
        if (pa < 0) pa += 360;
    }

    // Update dropdown to show preset if matching
    const select = document.getElementById('los-direction');
    if (select) {
        if (Math.abs(L.x - 1) < 0.01 && Math.abs(L.y) < 0.01 && Math.abs(L.z) < 0.01) {
            select.value = 'x';
        } else if (Math.abs(L.x) < 0.01 && Math.abs(L.y - 1) < 0.01 && Math.abs(L.z) < 0.01) {
            select.value = 'y';
        } else if (Math.abs(L.x) < 0.01 && Math.abs(L.y) < 0.01 && Math.abs(L.z - 1) < 0.01) {
            select.value = 'z';
        }
    }

    // Update display
    const vectorEl = document.getElementById('los-vector');
    const thetaEl = document.getElementById('los-theta');
    const phiEl = document.getElementById('los-phi');
    const inclinationEl = document.getElementById('los-inclination');
    const paEl = document.getElementById('los-pa');

    if (vectorEl) vectorEl.textContent = `(${L.x.toFixed(2)}, ${L.y.toFixed(2)}, ${L.z.toFixed(2)})`;
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

// LOS direction change handler (for dropdown)
function onLOSChange(direction) {
    switch (direction) {
        case 'x':
            STATE.losVector.set(1, 0, 0);
            break;
        case 'y':
            STATE.losVector.set(0, 1, 0);
            break;
        case 'z':
            STATE.losVector.set(0, 0, 1);
            break;
    }
    updateLOSArrows();
    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

// ============================================================
// TOGGLE FUNCTIONS for visibility control
// ============================================================

// Toggle sky plane (yellow ring + N/E arrows + LOS arrow)
function toggleSkyPlane() {
    if (STATE.fixedLOSGroup) {
        STATE.fixedLOSGroup.visible = !STATE.fixedLOSGroup.visible;
        console.log('Sky plane & LOS:', STATE.fixedLOSGroup.visible ? 'visible' : 'hidden');
    }
}

// Toggle orbital pole indicator (cyan sphere + arrow + orbital plane ring)
function toggleOrbitalPole() {
    if (STATE.losSphere) {
        STATE.losSphere.visible = !STATE.losSphere.visible;
    }
    if (STATE.losArrows) {
        STATE.losArrows.visible = !STATE.losArrows.visible;
    }
    console.log('Orbital pole:', STATE.losSphere?.visible ? 'visible' : 'hidden');
}

// Toggle analytic orbit line
function toggleOrbitLine() {
    if (STATE.orbitLine) {
        STATE.orbitLine.visible = !STATE.orbitLine.visible;
    }
    if (STATE.periMarker) {
        STATE.periMarker.visible = !STATE.periMarker.visible;
    }
    console.log('Orbit line:', STATE.orbitLine?.visible ? 'visible' : 'hidden');
}

// Toggle grid helper
function toggleGrid() {
    const grid = STATE.scene.children.find(c => c.type === 'GridHelper');
    if (grid) {
        grid.visible = !grid.visible;
        console.log('Grid:', grid.visible ? 'visible' : 'hidden');
    }
}
