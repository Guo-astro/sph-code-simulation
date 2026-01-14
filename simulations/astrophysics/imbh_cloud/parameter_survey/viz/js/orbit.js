// ============================================================
// IMBH-Cloud Visualization - Analytic Orbit Visualization
// Parabolic orbit computation and rendering
// ============================================================

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
