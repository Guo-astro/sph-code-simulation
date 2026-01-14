// ============================================================
// IMBH-Cloud Visualization - Orbital Plane Rotation Functions
// Handles rotation of orbit and particles when orbital pole changes
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

        // Initialize the user-controlled pole to match the original
        STATE.losVector.copy(STATE.originalOrbitalPole);

        console.log('Original orbital pole set from L_vec:', STATE.originalOrbitalPole);
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
    updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    console.log('Orbital pole reset to original:', STATE.originalOrbitalPole);
}
