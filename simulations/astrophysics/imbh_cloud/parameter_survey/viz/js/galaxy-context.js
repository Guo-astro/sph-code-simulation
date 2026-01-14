// ============================================================
// IMBH-Cloud Visualization - Galaxy Context
// Shows HVCC location relative to Sgr A*, Sun, and galactic plane
// Using SIMULATION SCALE (1 unit = 1 pc)
// ============================================================

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

    // Show HVCC offset from galactic plane - simple dashed line
    const hvccOffsetLineGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(0, sgrAY, 0)
    ]);
    const hvccOffsetLine = new THREE.Line(hvccOffsetLineGeometry, new THREE.LineDashedMaterial({
        color: 0x00ffff,
        dashSize: 2,
        gapSize: 1,
        transparent: true,
        opacity: 0.7
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

    const sunDistLabel = createTextSprite('~8 kpc (schematic)', '#aaaa00', 18);
    sunDistLabel.position.set(sgrAX, sgrAY + 5, sunZ);
    STATE.galaxyGroup.add(sunDistLabel);

    // LOS lines removed - they added visual clutter
    // The yellow sky plane in the main view shows the actual LOS direction

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

    // l coordinate axis (along galactic plane at b=0°) - simple dashed line
    const lAxisGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(sgrAX - 70, sgrAY, 0),
        new THREE.Vector3(sgrAX + 70, sgrAY, 0)
    ]);
    const lAxis = new THREE.Line(lAxisGeometry, new THREE.LineDashedMaterial({
        color: gridAxisColor,
        transparent: true,
        opacity: 0.6,
        dashSize: 2,
        gapSize: 1
    }));
    lAxis.computeLineDistances();
    STATE.galaxyGroup.add(lAxis);

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

    // b coordinate axis (perpendicular to galactic plane at l=0°, through Sgr A*) - simple dashed line
    const bAxisGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(sgrAX, sgrAY - 40, 0),
        new THREE.Vector3(sgrAX, sgrAY + 40, 0)
    ]);
    const bAxis = new THREE.Line(bAxisGeometry, new THREE.LineDashedMaterial({
        color: gridAxisColor,
        transparent: true,
        opacity: 0.6,
        dashSize: 2,
        gapSize: 1
    }));
    bAxis.computeLineDistances();
    STATE.galaxyGroup.add(bAxis);

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

    // Mark HVCC position on coordinate grid - simple dashed lines only
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
        opacity: 0.6
    }));
    hvccGridLine.computeLineDistances();
    STATE.galaxyGroup.add(hvccGridLine);

    // === GALACTIC ROTATION VISUALIZATION ===
    // Show differential rotation around Sgr A*
    // Circular velocity V_c(R) ~ 220 km/s (roughly flat rotation curve)
    STATE.galacticRotationGroup = new THREE.Group();

    // Color for rotation visualization
    const rotationColor = 0x66aaff;  // Light blue for motion vectors

    // Draw rotation arrows at different radii in galactic plane
    const rotationRadii = [20, 40, 60];  // pc from Sgr A*
    const arrowLength = 8;  // Visual length

    rotationRadii.forEach(radius => {
        // Create 8 arrows at each radius (45° apart) for clearer pattern
        // Includes 90° which is the direction toward the Sun
        for (let i = 0; i < 8; i++) {
            const angle = (i * Math.PI / 4);  // 0°, 45°, 90°, 135°, 180°, 225°, 270°, 315°

            // Position on circle (in galactic plane at Y = sgrAY)
            const px = sgrAX + radius * Math.cos(angle);
            const pz = radius * Math.sin(angle);

            // Rotation direction: CLOCKWISE when viewed from +Y (North Galactic Pole)
            // This is the actual direction of Galactic rotation
            // At angle=0 (+X from GC): velocity toward -Z (away from Sun)
            // At angle=π/2 (+Z from GC, toward Sun): velocity toward +X
            const tangentX = Math.sin(angle);
            const tangentZ = -Math.cos(angle);

            // Highlight the arrow at 90° (toward Sun direction)
            const isTowardSun = (Math.abs(angle - Math.PI/2) < 0.1);
            const arrowColor = isTowardSun ? 0xffff00 : rotationColor;  // Yellow if toward Sun

            const arrow = new THREE.ArrowHelper(
                new THREE.Vector3(tangentX, 0, tangentZ),
                new THREE.Vector3(px, sgrAY, pz),
                arrowLength,
                arrowColor,
                arrowLength * 0.3,
                arrowLength * 0.15
            );
            STATE.galacticRotationGroup.add(arrow);
        }
    });

    // Add rotation direction label (clockwise when viewed from North)
    const rotLabel = createTextSprite('Galactic Rotation (CW from N)', '#66aaff', 20);
    rotLabel.position.set(sgrAX - 55, sgrAY + 5, 30);
    STATE.galacticRotationGroup.add(rotLabel);

    const rotVelLabel = createTextSprite('V_c ≈ 220 km/s', '#66aaff', 18);
    rotVelLabel.position.set(sgrAX - 55, sgrAY + 2, 30);
    STATE.galacticRotationGroup.add(rotVelLabel);

    // Sun's orbital velocity arrow (at Sun position)
    // Sun moves at ~220 km/s in the direction of Galactic rotation
    // At the Sun's position, the rotation is toward positive l
    const sunVelArrow = new THREE.ArrowHelper(
        new THREE.Vector3(1, 0, 0),  // Sun moves toward +l (Galactic East)
        new THREE.Vector3(sgrAX - 5, sgrAY, sunZ),
        15,  // Longer arrow for emphasis
        0xffff00,  // Yellow to match Sun
        5,
        2.5
    );
    STATE.galacticRotationGroup.add(sunVelArrow);

    const sunVelLabel = createTextSprite('V_☉ ≈ 220 km/s', '#ffff88', 20);
    sunVelLabel.position.set(sgrAX + 15, sgrAY + 4, sunZ);
    STATE.galacticRotationGroup.add(sunVelLabel);

    // LSR explanation box
    const lsrLabel1 = createTextSprite('LSR (Local Standard of Rest):', '#88ddff', 22);
    lsrLabel1.position.set(sgrAX, sgrAY - 50, sunZ - 10);
    STATE.galacticRotationGroup.add(lsrLabel1);

    const lsrLabel2 = createTextSprite('Reference frame co-rotating with', '#88ddff', 18);
    lsrLabel2.position.set(sgrAX, sgrAY - 54, sunZ - 10);
    STATE.galacticRotationGroup.add(lsrLabel2);

    const lsrLabel3 = createTextSprite('the Galaxy at the Sun\'s position', '#88ddff', 18);
    lsrLabel3.position.set(sgrAX, sgrAY - 58, sunZ - 10);
    STATE.galacticRotationGroup.add(lsrLabel3);

    // Key insight about HVCC/IMBH velocity
    const hvccVelLabel1 = createTextSprite('HVCC observed: V_LSR ≈ −100 to −200 km/s', '#00ffff', 20);
    hvccVelLabel1.position.set(0, -50, 0);
    STATE.galacticRotationGroup.add(hvccVelLabel1);

    const hvccVelLabel2 = createTextSprite('(Simulation frame = IMBH rest frame)', '#00aaaa', 18);
    hvccVelLabel2.position.set(0, -54, 0);
    STATE.galacticRotationGroup.add(hvccVelLabel2);

    STATE.galaxyGroup.add(STATE.galacticRotationGroup);

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

    if (STATE.galaxyGroup.visible) {
        // Zoom out to see the Sun (at Z~70) and Sgr A* (at X~52, Y~29)
        animateCameraTo(
            new THREE.Vector3(30, 50, 120),  // Camera position to see Sun
            new THREE.Vector3(25, 15, 30)    // Look at point between HVCC and Sun
        );
        console.log('Galaxy context: visible (zoomed out to see Sun)');
    } else {
        // Zoom back to simulation scale
        animateCameraTo(
            new THREE.Vector3(15, 10, 20),
            new THREE.Vector3(-5, 0, 0)
        );
        console.log('Galaxy context: hidden (zoomed back to simulation)');
    }
}

// Toggle galactic rotation arrows visibility (within galaxy context)
function toggleGalacticRotation() {
    if (!STATE.galacticRotationGroup) return;
    STATE.galacticRotationGroup.visible = !STATE.galacticRotationGroup.visible;
    console.log('Galactic rotation:', STATE.galacticRotationGroup.visible ? 'visible' : 'hidden');
}
