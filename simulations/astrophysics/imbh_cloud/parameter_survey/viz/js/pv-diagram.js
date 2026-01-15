// ============================================================
// IMBH-Cloud Visualization - Position-Velocity Diagram
// Uses arbitrary LOS direction (draggable)
// Supports LSR velocity frame transformation
// ============================================================

// LSR frame state
STATE.useLSRFrame = false;
STATE.lsrVelocityOffset = -150;  // V_IMBH in LSR frame (km/s)

// Toggle LSR frame display
function toggleLSRFrame() {
    STATE.useLSRFrame = document.getElementById('lsr-toggle').checked;
    const sliderContainer = document.getElementById('lsr-slider-container');
    const offsetDisplay = document.getElementById('lsr-offset-display');

    sliderContainer.style.display = STATE.useLSRFrame ? 'block' : 'none';
    if (offsetDisplay) {
        offsetDisplay.style.display = STATE.useLSRFrame ? 'inline' : 'none';
    }

    // Update PV diagram
    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

// Update LSR velocity offset
function updateLSROffset(value) {
    STATE.lsrVelocityOffset = parseFloat(value);
    document.getElementById('lsr-offset-value').textContent = value;

    // Update PV diagram
    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

function updatePVDiagram(data) {
    const canvas = document.getElementById('pv-canvas');
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;

    // Clear canvas
    ctx.fillStyle = '#1a1a2a';
    ctx.fillRect(0, 0, width, height);

    if (!data || !data.particles || data.particles.length === 0) return;

    // Get real particles only
    const realParticles = data.particles.filter(p => p.is_ghost !== 1);

    // LOS direction (unit vector)
    const los = STATE.losVector;

    // Create perpendicular axes for position projection
    // perp1 = any vector perpendicular to LOS
    let perp1 = new THREE.Vector3();
    if (Math.abs(los.x) < 0.9) {
        perp1.crossVectors(los, new THREE.Vector3(1, 0, 0)).normalize();
    } else {
        perp1.crossVectors(los, new THREE.Vector3(0, 1, 0)).normalize();
    }

    // Compute position along perp1 and velocity along LOS for each particle
    let posValues = [];

    // Get LSR offset if enabled
    const lsrOffset = STATE.useLSRFrame ? STATE.lsrVelocityOffset : 0;

    realParticles.forEach(p => {
        const pos3d = new THREE.Vector3(p.x, p.y, p.z);
        const vel3d = new THREE.Vector3(p.vx, p.vy, p.vz);

        // Position projected onto perp1 axis
        const pos = pos3d.dot(perp1);

        // Velocity projected onto LOS (v_LOS)
        // In LSR frame: V_LSR = V_sim + V_IMBH_LSR (where V_IMBH_LSR is the IMBH bulk velocity)
        const vel = vel3d.dot(los) + lsrOffset;

        posValues.push({ pos, vel, particle: p });
    });

    // Find ranges
    let minPos = Infinity, maxPos = -Infinity;
    let minVel = Infinity, maxVel = -Infinity;

    posValues.forEach(({ pos, vel }) => {
        minPos = Math.min(minPos, pos);
        maxPos = Math.max(maxPos, pos);
        minVel = Math.min(minVel, vel);
        maxVel = Math.max(maxVel, vel);
    });

    // Add padding
    const posPad = (maxPos - minPos) * 0.1 || 1;
    const velPad = (maxVel - minVel) * 0.1 || 1;
    minPos -= posPad; maxPos += posPad;
    minVel -= velPad; maxVel += velPad;

    // Apply zoom (zoom into center)
    const zoom = STATE.pvZoom || 1.0;
    if (zoom !== 1.0) {
        const posMid = (minPos + maxPos) / 2;
        const velMid = (minVel + maxVel) / 2;
        const posRange = (maxPos - minPos) / zoom;
        const velRange = (maxVel - minVel) / zoom;
        minPos = posMid - posRange / 2;
        maxPos = posMid + posRange / 2;
        minVel = velMid - velRange / 2;
        maxVel = velMid + velRange / 2;
    }

    // Store ranges for selection coordinate conversion
    STATE.pvRanges = { minPos, maxPos, minVel, maxVel };

    // Margins for axes
    const margin = { left: 42, right: 10, top: 10, bottom: 28 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    // Draw axes
    ctx.strokeStyle = '#4a4a6a';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, height - margin.bottom);
    ctx.lineTo(width - margin.right, height - margin.bottom);
    ctx.stroke();

    // Format LOS direction for label
    const losLabel = `LOS: (${los.x.toFixed(2)}, ${los.y.toFixed(2)}, ${los.z.toFixed(2)})`;

    // Axis labels
    ctx.fillStyle = '#aaa';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';

    // X-axis: position perpendicular to LOS
    ctx.fillText('Offset (pc)', margin.left + plotWidth/2, height - 3);

    // Y-axis: V_LOS or V_LSR
    ctx.save();
    ctx.translate(12, margin.top + plotHeight/2);
    ctx.rotate(-Math.PI/2);
    if (STATE.useLSRFrame) {
        ctx.fillStyle = '#88ddff';  // Light blue for LSR
        ctx.fillText('V_LSR (km/s)', 0, 0);
    } else {
        ctx.fillText('V_LOS (km/s)', 0, 0);
    }
    ctx.restore();

    // Tick marks
    ctx.fillStyle = '#888';
    ctx.font = '9px monospace';

    // X ticks
    ctx.textAlign = 'center';
    for (let i = 0; i <= 4; i++) {
        const val = minPos + (maxPos - minPos) * i / 4;
        const x = margin.left + plotWidth * i / 4;
        ctx.fillText(val.toFixed(1), x, height - margin.bottom + 12);

        // Grid line
        ctx.strokeStyle = '#2a2a3a';
        ctx.beginPath();
        ctx.moveTo(x, margin.top);
        ctx.lineTo(x, height - margin.bottom);
        ctx.stroke();
    }

    // Y ticks
    ctx.textAlign = 'right';
    for (let i = 0; i <= 4; i++) {
        const val = minVel + (maxVel - minVel) * i / 4;
        const y = height - margin.bottom - plotHeight * i / 4;
        ctx.fillText(val.toFixed(1), margin.left - 4, y + 3);

        // Grid line
        ctx.strokeStyle = '#2a2a3a';
        ctx.beginPath();
        ctx.moveTo(margin.left, y);
        ctx.lineTo(width - margin.right, y);
        ctx.stroke();
    }

    // Draw zero velocity line if in range
    if (minVel < 0 && maxVel > 0) {
        const y0 = height - margin.bottom - (0 - minVel) / (maxVel - minVel) * plotHeight;
        ctx.strokeStyle = '#666';
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(margin.left, y0);
        ctx.lineTo(width - margin.right, y0);
        ctx.stroke();
        ctx.setLineDash([]);
    }

    // Create set of selected particle IDs for fast lookup
    const selectedSet = new Set(STATE.selectedParticles.map(i => data.particles[i]?.id));

    // Plot particles with SAME Jet colormap as main visualization
    posValues.forEach(({ pos, vel, particle }) => {
        const x = margin.left + (pos - minPos) / (maxPos - minPos) * plotWidth;
        const y = height - margin.bottom - (vel - minVel) / (maxVel - minVel) * plotHeight;

        // Skip particles outside zoom range
        if (x < margin.left || x > width - margin.right ||
            y < margin.top || y > height - margin.bottom) {
            return;
        }

        // Check if selected
        const isSelected = selectedSet.has(particle.id);

        if (isSelected) {
            // Highlight selected particles
            ctx.fillStyle = '#ff8866';
            ctx.globalAlpha = 1.0;
            ctx.beginPath();
            ctx.arc(x, y, 2.5, 0, Math.PI * 2);
            ctx.fill();
        } else {
            // Use same Jet colormap as main viz
            const t = getColorForValue(particle);
            const color = valueToColor(t);
            ctx.fillStyle = `rgb(${Math.round(color.r*255)}, ${Math.round(color.g*255)}, ${Math.round(color.b*255)})`;
            ctx.globalAlpha = 0.5;
            ctx.beginPath();
            ctx.arc(x, y, 1.5, 0, Math.PI * 2);
            ctx.fill();
        }
    });

    ctx.globalAlpha = 1.0;

    // Draw selection rectangle if selecting
    if (STATE.pvSelection && STATE.pvSelecting) {
        const sel = STATE.pvSelection;
        ctx.strokeStyle = '#ff8866';
        ctx.lineWidth = 2;
        ctx.setLineDash([4, 4]);
        ctx.strokeRect(
            Math.min(sel.x1, sel.x2),
            Math.min(sel.y1, sel.y2),
            Math.abs(sel.x2 - sel.x1),
            Math.abs(sel.y2 - sel.y1)
        );
        ctx.setLineDash([]);
    }

    // Draw LSR frame indicator if active
    if (STATE.useLSRFrame) {
        ctx.fillStyle = 'rgba(100, 170, 255, 0.2)';
        ctx.fillRect(margin.left, margin.top, plotWidth, 18);

        ctx.fillStyle = '#88ddff';
        ctx.font = 'bold 10px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(`LSR Frame (V_IMBH = ${STATE.lsrVelocityOffset} km/s)`, margin.left + 4, margin.top + 12);
    }

    // Update axis labels in HTML
    document.getElementById('pv-x-label').textContent = 'Offset (pc)';
    document.getElementById('pv-v-label').textContent = STATE.useLSRFrame ? 'V_LSR (km/s)' : 'V_LOS (km/s)';
}

function onLOSChange(direction) {
    // Set LOS to preset direction
    switch (direction) {
        case 'x':
            STATE.losVector.set(1, 0, 0);
            break;
        case 'y':
            STATE.losVector.set(0, 1, 0);
            break;
        case 'z':
        default:
            STATE.losVector.set(0, 0, 1);
            break;
    }

    // Update LOS arrows in 3D view
    updateLOSArrows();

    // Update PV diagram
    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

// PV Diagram zoom functions
function pvZoom(factor) {
    STATE.pvZoom = Math.max(0.5, Math.min(10, STATE.pvZoom * factor));
    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

function pvZoomReset() {
    STATE.pvZoom = 1.0;
    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
    }
}

// ============================================================
// PV Selection functionality
// ============================================================

function initPVSelection() {
    const canvas = document.getElementById('pv-canvas');

    canvas.addEventListener('mousedown', (e) => {
        if (e.button !== 0) return;  // Left click only
        const rect = canvas.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;

        STATE.pvSelecting = true;
        STATE.pvSelection = { x1: x, y1: y, x2: x, y2: y };
    });

    canvas.addEventListener('mousemove', (e) => {
        if (!STATE.pvSelecting) return;
        const rect = canvas.getBoundingClientRect();
        STATE.pvSelection.x2 = e.clientX - rect.left;
        STATE.pvSelection.y2 = e.clientY - rect.top;

        // Redraw with selection box
        if (STATE.snapshots.length > 0) {
            updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
        }
    });

    canvas.addEventListener('mouseup', (e) => {
        if (!STATE.pvSelecting) return;
        STATE.pvSelecting = false;

        // Finalize selection
        if (STATE.pvSelection && STATE.pvRanges) {
            selectParticlesInPVBox();
        }
    });

    canvas.addEventListener('mouseleave', (e) => {
        if (STATE.pvSelecting) {
            STATE.pvSelecting = false;
            if (STATE.pvSelection && STATE.pvRanges) {
                selectParticlesInPVBox();
            }
        }
    });
}

function selectParticlesInPVBox() {
    if (!STATE.pvSelection || !STATE.pvRanges || STATE.snapshots.length === 0) return;

    const canvas = document.getElementById('pv-canvas');
    const margin = { left: 42, right: 10, top: 10, bottom: 28 };
    const plotWidth = canvas.width - margin.left - margin.right;
    const plotHeight = canvas.height - margin.top - margin.bottom;

    const { minPos, maxPos, minVel, maxVel } = STATE.pvRanges;

    // Convert canvas coords to PV coords
    const sel = STATE.pvSelection;
    const x1 = Math.min(sel.x1, sel.x2);
    const x2 = Math.max(sel.x1, sel.x2);
    const y1 = Math.min(sel.y1, sel.y2);
    const y2 = Math.max(sel.y1, sel.y2);

    // Min selection size
    if (Math.abs(x2 - x1) < 5 && Math.abs(y2 - y1) < 5) {
        clearPVSelection();
        return;
    }

    // Convert to PV space
    const posMin = minPos + (x1 - margin.left) / plotWidth * (maxPos - minPos);
    const posMax = minPos + (x2 - margin.left) / plotWidth * (maxPos - minPos);
    const velMax = maxVel - (y1 - margin.top) / plotHeight * (maxVel - minVel);
    const velMin = maxVel - (y2 - margin.top) / plotHeight * (maxVel - minVel);

    // Find particles in selection
    const data = STATE.snapshots[STATE.currentFrame];
    const los = STATE.losVector;

    // Perpendicular axis
    let perp1 = new THREE.Vector3();
    if (Math.abs(los.x) < 0.9) {
        perp1.crossVectors(los, new THREE.Vector3(1, 0, 0)).normalize();
    } else {
        perp1.crossVectors(los, new THREE.Vector3(0, 1, 0)).normalize();
    }

    STATE.selectedParticles = [];

    // Get LSR offset if enabled
    const lsrOffset = STATE.useLSRFrame ? STATE.lsrVelocityOffset : 0;

    data.particles.forEach((p, idx) => {
        if (p.is_ghost === 1) return;

        const pos3d = new THREE.Vector3(p.x, p.y, p.z);
        const vel3d = new THREE.Vector3(p.vx, p.vy, p.vz);

        const pos = pos3d.dot(perp1);
        const vel = vel3d.dot(los) + lsrOffset;  // Apply same offset as in diagram

        if (pos >= posMin && pos <= posMax && vel >= velMin && vel <= velMax) {
            STATE.selectedParticles.push(idx);
        }
    });

    // Update analysis panel
    updateSelectionAnalysis();

    // Sync with Lagrangian tracking system
    if (typeof syncTrackingFromPVSelection === 'function') {
        syncTrackingFromPVSelection();
    }

    // Redraw to show selection
    updatePVDiagram(data);

    // Update 3D visualization to highlight selected particles
    updateVisualization(STATE.currentFrame);
}

function clearPVSelection() {
    STATE.pvSelection = null;
    STATE.selectedParticles = [];
    document.getElementById('selection-panel').style.display = 'none';

    // Clear profile selection if exists
    if (typeof clearProfileSelection === 'function') {
        clearProfileSelection();
    }

    if (STATE.snapshots.length > 0) {
        updatePVDiagram(STATE.snapshots[STATE.currentFrame]);
        updateVisualization(STATE.currentFrame);
    }
}

function updateSelectionAnalysis() {
    const panel = document.getElementById('selection-panel');

    if (STATE.selectedParticles.length === 0) {
        panel.style.display = 'none';
        return;
    }

    panel.style.display = 'block';

    const data = STATE.snapshots[STATE.currentFrame];
    const particles = STATE.selectedParticles.map(i => data.particles[i]);

    // Compute statistics
    let totalMass = 0;
    let sumX = 0, sumY = 0, sumZ = 0;
    let sumX2 = 0, sumY2 = 0, sumZ2 = 0;
    let sumDens = 0, maxDens = 0;
    let sumTemp = 0, maxTemp = 0;
    let sumVel = 0, maxVel = 0;
    let sumMach = 0, maxMach = 0;

    // Get COM velocity for Mach calculation
    let comVx = 0, comVy = 0, comVz = 0, comMass = 0;
    particles.forEach(p => {
        comVx += p.vx * p.mass;
        comVy += p.vy * p.mass;
        comVz += p.vz * p.mass;
        comMass += p.mass;
    });
    comVx /= comMass; comVy /= comMass; comVz /= comMass;

    particles.forEach(p => {
        totalMass += p.mass;
        sumX += p.x; sumY += p.y; sumZ += p.z;
        sumX2 += p.x * p.x; sumY2 += p.y * p.y; sumZ2 += p.z * p.z;

        const n_H2 = p.dens * CONFIG.densityToNH2;
        sumDens += n_H2;
        maxDens = Math.max(maxDens, n_H2);

        sumTemp += p.temp;
        maxTemp = Math.max(maxTemp, p.temp);

        const dvx = p.vx - comVx, dvy = p.vy - comVy, dvz = p.vz - comVz;
        const vRel = Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        sumVel += vRel;
        maxVel = Math.max(maxVel, vRel);

        if (p.sound > 0) {
            const mach = vRel / p.sound;
            sumMach += mach;
            maxMach = Math.max(maxMach, mach);
        }
    });

    const n = particles.length;
    const meanX = sumX / n, meanY = sumY / n, meanZ = sumZ / n;
    const stdX = Math.sqrt(sumX2 / n - meanX * meanX);
    const stdY = Math.sqrt(sumY2 / n - meanY * meanY);
    const stdZ = Math.sqrt(sumZ2 / n - meanZ * meanZ);

    // Update panel
    document.getElementById('sel-count').textContent = n;
    document.getElementById('sel-mass').textContent = totalMass.toFixed(2) + ' M☉';
    document.getElementById('sel-x').textContent = `${meanX.toFixed(2)} ± ${stdX.toFixed(2)} pc`;
    document.getElementById('sel-y').textContent = `${meanY.toFixed(2)} ± ${stdY.toFixed(2)} pc`;
    document.getElementById('sel-z').textContent = `${meanZ.toFixed(2)} ± ${stdZ.toFixed(2)} pc`;
    document.getElementById('sel-dens').textContent = `${(sumDens/n).toExponential(1)} / ${maxDens.toExponential(1)}`;
    document.getElementById('sel-temp').textContent = `${(sumTemp/n).toFixed(0)} / ${maxTemp.toFixed(0)} K`;
    document.getElementById('sel-vel').textContent = `${(sumVel/n).toFixed(2)} / ${maxVel.toFixed(2)} km/s`;
    document.getElementById('sel-mach').textContent = `${(sumMach/n).toFixed(2)} / ${maxMach.toFixed(2)}`;

    // Draw phase diagram
    drawPhaseDiagram(particles);
}

function drawPhaseDiagram(particles) {
    const canvas = document.getElementById('phase-canvas');
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;

    ctx.fillStyle = '#1a1a2a';
    ctx.fillRect(0, 0, width, height);

    if (particles.length === 0) return;

    // Compute ranges
    let minLogRho = Infinity, maxLogRho = -Infinity;
    let minLogT = Infinity, maxLogT = -Infinity;

    particles.forEach(p => {
        const logRho = Math.log10(Math.max(p.dens * CONFIG.densityToNH2, 1));
        const logT = Math.log10(Math.max(p.temp, 1));
        minLogRho = Math.min(minLogRho, logRho);
        maxLogRho = Math.max(maxLogRho, logRho);
        minLogT = Math.min(minLogT, logT);
        maxLogT = Math.max(maxLogT, logT);
    });

    // Add padding
    const pad = 0.2;
    minLogRho -= pad; maxLogRho += pad;
    minLogT -= pad; maxLogT += pad;

    const margin = { left: 30, right: 5, top: 5, bottom: 20 };
    const plotW = width - margin.left - margin.right;
    const plotH = height - margin.top - margin.bottom;

    // Axes
    ctx.strokeStyle = '#4a4a6a';
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, height - margin.bottom);
    ctx.lineTo(width - margin.right, height - margin.bottom);
    ctx.stroke();

    // Labels
    ctx.fillStyle = '#888';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('log ρ (cm⁻³)', margin.left + plotW/2, height - 2);

    ctx.save();
    ctx.translate(10, margin.top + plotH/2);
    ctx.rotate(-Math.PI/2);
    ctx.fillText('log T (K)', 0, 0);
    ctx.restore();

    // Plot points
    ctx.globalAlpha = 0.7;
    particles.forEach(p => {
        const logRho = Math.log10(Math.max(p.dens * CONFIG.densityToNH2, 1));
        const logT = Math.log10(Math.max(p.temp, 1));

        const x = margin.left + (logRho - minLogRho) / (maxLogRho - minLogRho) * plotW;
        const y = height - margin.bottom - (logT - minLogT) / (maxLogT - minLogT) * plotH;

        ctx.fillStyle = '#ff8866';
        ctx.beginPath();
        ctx.arc(x, y, 2, 0, Math.PI * 2);
        ctx.fill();
    });
    ctx.globalAlpha = 1.0;
}

// Initialize selection on page load
document.addEventListener('DOMContentLoaded', initPVSelection);
