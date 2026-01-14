// ============================================================
// IMBH-Cloud Visualization - Camera View Controls
// ============================================================

// Store initial camera position for reset
const INITIAL_CAMERA = {
    position: new THREE.Vector3(0, 0, 25),
    target: new THREE.Vector3(0, 0, 0)
};

function setEdgeOnView() {
    // Edge-on view: looking along Z-axis (X-Y plane view)
    // Camera on Z-axis looking at origin
    animateCameraTo(
        new THREE.Vector3(0, 0, 25),
        new THREE.Vector3(0, 0, 0)
    );
}

function setFaceOnView() {
    // Face-on view: looking along Y-axis (X-Z plane view)
    // Camera on Y-axis looking at origin
    animateCameraTo(
        new THREE.Vector3(0, 25, 0),
        new THREE.Vector3(0, 0, 0)
    );
}

function setTopView() {
    // Top view: looking along X-axis (Y-Z plane view)
    // Camera on X-axis looking at origin
    animateCameraTo(
        new THREE.Vector3(25, 0, 0),
        new THREE.Vector3(0, 0, 0)
    );
}

function resetView() {
    // Reset to initial isometric-like view
    animateCameraTo(
        new THREE.Vector3(15, 10, 20),
        new THREE.Vector3(-5, 0, 0)
    );
}

function animateCameraTo(targetPos, lookAt, duration = 500) {
    const startPos = STATE.camera.position.clone();
    const startTarget = STATE.controls.target.clone();

    const startTime = Date.now();

    function updateCamera() {
        const elapsed = Date.now() - startTime;
        const t = Math.min(elapsed / duration, 1);

        // Smooth easing
        const ease = t < 0.5
            ? 2 * t * t
            : 1 - Math.pow(-2 * t + 2, 2) / 2;

        // Interpolate position
        STATE.camera.position.lerpVectors(startPos, targetPos, ease);

        // Interpolate target
        STATE.controls.target.lerpVectors(startTarget, lookAt, ease);

        STATE.controls.update();

        if (t < 1) {
            requestAnimationFrame(updateCamera);
        }
    }

    updateCamera();
}

// Set camera to follow cloud COM (optional feature)
function setCameraFollowCOM() {
    if (!STATE.comMarker) return;

    const com = STATE.comMarker.position;
    const offset = new THREE.Vector3(0, 0, 15);

    animateCameraTo(
        com.clone().add(offset),
        com.clone()
    );
}

// Toggle coordinate labels visibility
function toggleCoordinateLabels() {
    if (!STATE.coordinateLabels) return;

    STATE.coordLabelsVisible = !STATE.coordLabelsVisible;

    // Toggle visibility of all coordinate labels
    Object.values(STATE.coordinateLabels).forEach(obj => {
        if (obj && obj.visible !== undefined) {
            obj.visible = STATE.coordLabelsVisible;
        }
    });

    console.log('Coordinate labels:', STATE.coordLabelsVisible ? 'visible' : 'hidden');
}

// Initialize coordinate labels visibility state
STATE.coordLabelsVisible = true;

// ============================================================
// Orbital Pole Presets (for inclination relative to fixed LOS to Sun)
// ============================================================

// Set orbital pole to face-on view (i=0°)
// Orbital pole aligned with fixed LOS to Sun - we look straight down the pole
function setOrbitalPoleToFaceOn() {
    // Copy the fixed LOS direction to the orbital pole
    STATE.losVector.copy(STATE.fixedLOStoSun);
    updateLOSArrows();
    updatePVDiagram();
    console.log('Orbital pole set to face-on (i=0°): pole aligned with LOS to Sun');
}

// Set orbital pole to edge-on view (i=90°)
// Orbital pole perpendicular to fixed LOS to Sun - we look in the orbital plane
function setOrbitalPoleToEdgeOn() {
    // Find a direction perpendicular to the fixed LOS
    const n_sun = STATE.fixedLOStoSun;

    // Use cross product with a reference vector to get perpendicular direction
    // Choose reference that's not parallel to n_sun
    let ref = new THREE.Vector3(0, 0, 1);
    if (Math.abs(n_sun.z) > 0.9) {
        ref = new THREE.Vector3(1, 0, 0);  // Use X if LOS is near Z
    }

    // Perpendicular direction = n_sun × ref, normalized
    const perp = new THREE.Vector3().crossVectors(n_sun, ref).normalize();

    STATE.losVector.copy(perp);
    updateLOSArrows();
    updatePVDiagram();
    console.log('Orbital pole set to edge-on (i=90°): pole perpendicular to LOS to Sun');
}
