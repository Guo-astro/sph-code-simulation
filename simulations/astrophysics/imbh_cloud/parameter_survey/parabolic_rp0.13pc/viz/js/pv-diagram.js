// ============================================================
// IMBH-Cloud Visualization - Position-Velocity Diagram
// Uses arbitrary LOS direction (draggable)
// ============================================================

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

    realParticles.forEach(p => {
        const pos3d = new THREE.Vector3(p.x, p.y, p.z);
        const vel3d = new THREE.Vector3(p.vx, p.vy, p.vz);

        // Position projected onto perp1 axis
        const pos = pos3d.dot(perp1);

        // Velocity projected onto LOS (v_LOS)
        const vel = vel3d.dot(los);

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

    // Y-axis: V_LOS
    ctx.save();
    ctx.translate(12, margin.top + plotHeight/2);
    ctx.rotate(-Math.PI/2);
    ctx.fillText('V_LOS (km/s)', 0, 0);
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

    // Plot particles with SAME color as main visualization
    posValues.forEach(({ pos, vel, particle }) => {
        const x = margin.left + (pos - minPos) / (maxPos - minPos) * plotWidth;
        const y = height - margin.bottom - (vel - minVel) / (maxVel - minVel) * plotHeight;

        // Use same color function as main viz
        const t = getColorForValue(particle);
        const hue = 0.7 * (1 - t) * 360;

        ctx.fillStyle = `hsl(${hue}, 100%, 50%)`;
        ctx.globalAlpha = 0.7;
        ctx.beginPath();
        ctx.arc(x, y, 1.5, 0, Math.PI * 2);
        ctx.fill();
    });

    ctx.globalAlpha = 1.0;

    // Update axis labels in HTML
    document.getElementById('pv-x-label').textContent = 'Offset (pc)';
    document.getElementById('pv-v-label').textContent = 'V_LOS (km/s)';
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
