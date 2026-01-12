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
