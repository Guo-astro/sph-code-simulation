// ============================================================
// IMBH-Cloud Visualization - Particle System Rendering
// Updates particle positions, colors, and UI displays
// ============================================================

function updateVisualization(frameIndex) {
    if (frameIndex < 0 || frameIndex >= STATE.snapshots.length) return;

    STATE.currentFrame = frameIndex;
    const data = STATE.snapshots[frameIndex];

    // Sync selected particles from tracked IDs BEFORE rendering
    if (typeof syncSelectedParticlesBeforeRender === 'function') {
        syncSelectedParticlesBeforeRender();
    }

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

    // COM velocity for color mapping
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

    // Update Lagrangian tracking plots if computed
    if (typeof onFrameChangeUpdateTracking === 'function') {
        onFrameChangeUpdateTracking();
    }
}
