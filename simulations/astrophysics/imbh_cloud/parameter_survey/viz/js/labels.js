// ============================================================
// IMBH-Cloud Visualization - Text Labels and Coordinate System
// Clean, minimal labels for clarity
// ============================================================

function createTextSprite(text, color = '#ffffff', fontSize = 48) {
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    canvas.width = 256;
    canvas.height = 128;

    context.fillStyle = 'transparent';
    context.fillRect(0, 0, canvas.width, canvas.height);

    context.font = `Bold ${fontSize}px Arial`;
    context.fillStyle = color;
    context.textAlign = 'center';
    context.textBaseline = 'middle';
    context.fillText(text, canvas.width / 2, canvas.height / 2);

    const texture = new THREE.CanvasTexture(canvas);
    const material = new THREE.SpriteMaterial({ map: texture, transparent: true });
    const sprite = new THREE.Sprite(material);
    sprite.scale.set(2, 1, 1);
    return sprite;
}

// Create labeled coordinate system
// Coordinate convention:
//   X = +l direction (Galactic East)
//   Y = +b direction (Galactic North)
//   Z = away from Galactic Center (NOT the same as LOS to Sun!)
function createLabeledCoordinates() {
    const axisLength = 8;
    const labelOffset = 1.2;

    // Create groups for toggling
    STATE.coordAxesGroup = new THREE.Group();
    STATE.coordLabelsGroup = new THREE.Group();

    // === AXIS ARROWS ===
    const arrowSize = 0.3;

    // X-axis arrow (red)
    const xDir = new THREE.Vector3(1, 0, 0);
    const xArrow = new THREE.ArrowHelper(xDir, new THREE.Vector3(0, 0, 0), axisLength, 0xff4444, arrowSize * 2, arrowSize);
    STATE.coordAxesGroup.add(xArrow);

    // Y-axis arrow (green)
    const yDir = new THREE.Vector3(0, 1, 0);
    const yArrow = new THREE.ArrowHelper(yDir, new THREE.Vector3(0, 0, 0), axisLength, 0x44ff44, arrowSize * 2, arrowSize);
    STATE.coordAxesGroup.add(yArrow);

    // Z-axis arrow (blue)
    const zDir = new THREE.Vector3(0, 0, 1);
    const zArrow = new THREE.ArrowHelper(zDir, new THREE.Vector3(0, 0, 0), axisLength, 0x4444ff, arrowSize * 2, arrowSize);
    STATE.coordAxesGroup.add(zArrow);

    // === AXIS LABELS (minimal) ===
    const xLabel = createTextSprite('X', '#ff6666', 36);
    xLabel.position.set(axisLength * labelOffset, 0, 0);
    STATE.coordLabelsGroup.add(xLabel);

    const yLabel = createTextSprite('Y', '#66ff66', 36);
    yLabel.position.set(0, axisLength * labelOffset, 0);
    STATE.coordLabelsGroup.add(yLabel);

    const zLabel = createTextSprite('Z', '#6666ff', 36);
    zLabel.position.set(0, 0, axisLength * labelOffset);
    STATE.coordLabelsGroup.add(zLabel);

    // Origin label
    const originLabel = createTextSprite('IMBH', '#ff4444', 28);
    originLabel.position.set(0.8, 0.8, 0.8);
    STATE.coordLabelsGroup.add(originLabel);

    // Add groups to scene
    STATE.scene.add(STATE.coordAxesGroup);
    STATE.scene.add(STATE.coordLabelsGroup);

    // Store references for toggling
    STATE.coordinateLabels = {
        axesGroup: STATE.coordAxesGroup,
        labelsGroup: STATE.coordLabelsGroup
    };

    console.log('Coordinate system created (X=+l, Y=+b, Z=away from GC)');
}

// Toggle coordinate axes visibility
function toggleCoordAxes() {
    if (STATE.coordAxesGroup) {
        STATE.coordAxesGroup.visible = !STATE.coordAxesGroup.visible;
        console.log('Coordinate axes:', STATE.coordAxesGroup.visible ? 'visible' : 'hidden');
    }
}

// Toggle coordinate labels visibility
function toggleCoordLabels() {
    if (STATE.coordLabelsGroup) {
        STATE.coordLabelsGroup.visible = !STATE.coordLabelsGroup.visible;
        console.log('Coordinate labels:', STATE.coordLabelsGroup.visible ? 'visible' : 'hidden');
    }
}

// Legacy function for compatibility
function toggleCoordinateLabels() {
    toggleCoordAxes();
    toggleCoordLabels();
}
