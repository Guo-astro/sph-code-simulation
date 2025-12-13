'use client'

import * as THREE from 'three'

/**
 * Galactic Coordinate System Markers
 * 
 * Based on the IMBH-Cloud encounter scenario from Oka et al. (2017) Nature Astronomy:
 * 
 * Physical Parameters (Oka et al. 2017):
 * =========================================
 * Black Hole:
 *   - M_BH = 10^5 M☉
 *   - Distance from Galactic Center: ~60 pc (projected)
 *   - Point source size upper limit: < 0.022 pc
 * 
 * Molecular Cloud (CO-0.40-0.22):
 *   - Cloud size: ~5 pc (extended), dense clump ~0.3 pc
 *   - Galactic coordinates: (l, b) = (-0.398°, -0.224°)
 *   - HCN clump half-size: S = 0.15 pc
 *   - Velocity dispersion: σ_v = 22 km/s
 *   - V_LSR range: -105 to -5 km/s (main body: -80 to -40 km/s)
 *   - Temperature: 60 K
 *   - Dense clump mass: ~40 M☉
 *   - N-body model cloud mass: 1000 M☉
 * 
 * N-body Model (from Methods section):
 *   - Initial cloud σ_r = 0.2 pc
 *   - Initial position: (X, Y) = (9.8 pc, -0.65 pc)
 *   - Initial velocity: (vX, vY) = (-8.19 km/s, 0.4 km/s)
 *   - Best-fit time: t = 7.2 × 10^5 yr (just after pericentre)
 * 
 * Coordinate System:
 * - Origin: Black Hole position
 * - X: Points toward Galactic longitude l = 0° (roughly toward Sgr A*)
 * - Y: Points toward l = 90° (direction of Galactic rotation)
 * - Z: Points toward Galactic North Pole (b = 90°)
 * 
 * Observer Geometry:
 * - Inclination i = 70°: Angle between orbital plane normal and Line of Sight
 * - Position Angle PA = 41.6°: Rotation of the orbit projection on the sky
 * - Distance to Earth: 8.0 kpc (distance to Galactic Center)
 */

export interface GalacticMarkersConfig {
  // Physical distances (in pc for scene units)
  distanceToGC_pc?: number        // Distance from system to Galactic Center (~60 pc)
  distanceToGC?: number           // Alias for distanceToGC_pc (for backward compat)
  distanceToEarth_kpc?: number    // Distance from Earth to Galactic Center (8 kpc)

  // Source position
  galacticLongitude?: number   // Galactic longitude of the source (degrees)
  galacticLatitude?: number    // Galactic latitude of the source (degrees)

  // Viewing angles (from Oka et al. 2017)
  inclination?: number         // Orbital plane inclination (degrees)
  positionAngle?: number       // Position angle on sky (degrees)

  // LSR velocity
  lsrVelocity?: number         // V_LSR in km/s (negative = approaching)

  // HVCC (High-Velocity Compact Cloud) parameters from Oka et al. (2017)
  hvccSize_pc?: number         // HCN clump half-size S = 0.15 pc (full size ~0.3 pc)
  hvccSize?: number            // Alias for hvccSize_pc
  hvccPosition?: [number, number, number]  // Position offset of HVCC from BH
  cloudExtent_pc?: number      // Extended cloud size ~5 pc
  denseClumpOffset_pc?: number // Displacement between HCN peak and continuum ~0.2 pc

  // N-body model parameters (for scale reference)
  nbodyCloudSigma_pc?: number  // Initial cloud σ_r = 0.2 pc
  bhMass_Msun?: number         // BH mass 10^5 M☉
  cloudMass_Msun?: number      // Cloud mass ~1000 M☉

  // Simulation cloud parameters (from preset config)
  cloudRadius_pc?: number      // Cloud radius from simulation preset (e.g., 1.13 pc)

  // Display scale factor (scene units are pc)
  displayScale?: number        // Scale factor for visualization (1.0 = true scale)

  // Galactic visualization options
  showGalaxyDisk?: boolean     // Show low-opacity Milky Way disk
  showSolarSystem?: boolean    // Show Solar System demo at Sun position (enlarged for visibility)
  galaxyRotationSpeed?: number // Rotation speed (rad/s) for animation
}

// Default configuration based on Oka et al. (2017) Nature Astronomy
export const DEFAULT_GALACTIC_CONFIG: GalacticMarkersConfig = {
  // Physical distances
  distanceToGC_pc: 60,           // ~60 pc from Galactic Center (projected)
  distanceToGC: 60,              // Alias
  distanceToEarth_kpc: 8.0,      // 8 kpc to Galactic Center

  // Source position
  galacticLongitude: -0.398,     // l = -0.398° (or 359.602°)
  galacticLatitude: -0.224,      // b = -0.224°

  // Viewing geometry
  inclination: 70,               // 70° inclination
  positionAngle: 41.6,           // Position angle 41.6°

  // Kinematics
  lsrVelocity: -60,              // V_LSR ~ -60 km/s (mean of -80 to -40 km/s range)

  // Cloud sizes (TRUE physical values from Oka et al. 2017)
  hvccSize_pc: 0.15,             // HCN clump half-size S = 0.15 pc
  hvccSize: 0.15,                // Alias
  hvccPosition: [0.2, -0.1, 0],  // Position offset from BH
  cloudExtent_pc: 5.0,           // Extended cloud ~5 pc
  denseClumpOffset_pc: 0.2,      // Displacement ~0.2 pc

  // N-body model parameters
  nbodyCloudSigma_pc: 0.2,       // σ_r = 0.2 pc
  bhMass_Msun: 100000,           // 10^5 M☉
  cloudMass_Msun: 1000,          // 1000 M☉

  // Simulation cloud parameters (default from CAT_OKA preset)
  cloudRadius_pc: 1.13,          // Cloud radius from preset (1.13 pc)

  // Display scale (1.0 = true scale in pc)
  displayScale: 1.0,

  // Galactic visualization
  showGalaxyDisk: true,
  showSolarSystem: false,        // Solar System demo off by default
  galaxyRotationSpeed: 0.0,      // Static by default
}

// LARGER default text sprite options for better visibility
const LARGE_LABEL_OPTIONS = {
  fontSize: 72,           // Increased from 48
  scale: [12, 6, 1] as [number, number, number],  // Increased from [8, 4, 1]
}

/**
 * Create a text sprite for 3D labels
 */
export function createTextSprite(
  text: string,
  options: {
    color?: string
    fontSize?: number
    fontWeight?: string
    backgroundColor?: string
    padding?: number
    scale?: [number, number, number]
  } = {}
): THREE.Sprite {
  const {
    color = '#ffffff',
    fontSize = 72,        // Increased from 48 for better visibility
    fontWeight = 'Bold',
    backgroundColor,
    padding = 10,
    scale = [12, 6, 1]    // Increased from [8, 4, 1] for better visibility
  } = options

  const canvas = document.createElement('canvas')
  const size = 512
  canvas.width = size
  canvas.height = size / 2
  const ctx = canvas.getContext('2d')!

  // Background
  if (backgroundColor) {
    ctx.fillStyle = backgroundColor
    ctx.fillRect(0, 0, size, size / 2)
  }

  // Text
  ctx.fillStyle = color
  ctx.font = `${fontWeight} ${fontSize}px Arial, sans-serif`
  ctx.textAlign = 'center'
  ctx.textBaseline = 'middle'
  ctx.fillText(text, size / 2, size / 4)

  const texture = new THREE.CanvasTexture(canvas)
  texture.needsUpdate = true

  const material = new THREE.SpriteMaterial({
    map: texture,
    transparent: true,
    depthTest: false,
    depthWrite: false,
  })

  const sprite = new THREE.Sprite(material)
  sprite.scale.set(...scale)
  sprite.renderOrder = 999  // Render on top

  return sprite
}

/**
 * Create an arrow with cone head
 */
export function createArrow(
  from: THREE.Vector3,
  to: THREE.Vector3,
  color: number,
  options: {
    headLength?: number
    headWidth?: number
    shaftRadius?: number
    opacity?: number
  } = {}
): THREE.Group {
  const { headLength = 1.0, headWidth = 0.6, shaftRadius = 0.15, opacity = 1.0 } = options  // Thicker for visibility
  
  const group = new THREE.Group()
  const direction = to.clone().sub(from).normalize()
  const length = from.distanceTo(to)
  
  if (length <= headLength) return group

  // Shaft
  const shaftLength = length - headLength
  const shaftGeo = new THREE.CylinderGeometry(shaftRadius, shaftRadius, shaftLength, 12)
  const shaftMat = new THREE.MeshBasicMaterial({ color, transparent: opacity < 1, opacity })
  const shaft = new THREE.Mesh(shaftGeo, shaftMat)
  
  // Position shaft at midpoint between from and the base of the head
  const shaftMidpoint = from.clone().add(direction.clone().multiplyScalar(shaftLength / 2))
  shaft.position.copy(shaftMidpoint)
  shaft.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
  group.add(shaft)

  // Head (cone)
  const headGeo = new THREE.ConeGeometry(headWidth, headLength, 16)
  const headMat = new THREE.MeshBasicMaterial({ color, transparent: opacity < 1, opacity })
  const head = new THREE.Mesh(headGeo, headMat)
  
  // Position head at tip
  const headPosition = to.clone().sub(direction.clone().multiplyScalar(headLength / 2))
  head.position.copy(headPosition)
  head.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
  group.add(head)

  return group
}

/**
 * Create a dashed circle in 3D
 */
export function createDashedCircle(
  radius: number,
  color: number,
  options: {
    segments?: number
    dashSize?: number
    gapSize?: number
    opacity?: number
    normal?: THREE.Vector3
    center?: THREE.Vector3
  } = {}
): THREE.Line {
  const {
    segments = 64,
    dashSize = 0.5,
    gapSize = 0.25,
    opacity = 0.8,
    normal = new THREE.Vector3(0, 0, 1),
    center = new THREE.Vector3(0, 0, 0)
  } = options

  const points: THREE.Vector3[] = []
  
  // Create rotation to align circle with desired normal
  const defaultNormal = new THREE.Vector3(0, 0, 1)
  const quaternion = new THREE.Quaternion().setFromUnitVectors(defaultNormal, normal.clone().normalize())
  
  for (let i = 0; i <= segments; i++) {
    const theta = (i / segments) * Math.PI * 2
    const point = new THREE.Vector3(
      radius * Math.cos(theta),
      radius * Math.sin(theta),
      0
    )
    point.applyQuaternion(quaternion)
    point.add(center)
    points.push(point)
  }

  const geometry = new THREE.BufferGeometry().setFromPoints(points)
  const material = new THREE.LineDashedMaterial({
    color,
    dashSize,
    gapSize,
    transparent: true,
    opacity,
  })

  const line = new THREE.Line(geometry, material)
  line.computeLineDistances()

  return line
}

/**
 * Create the Sun marker (yellow sphere with corona effect)
 */
export function createSunMarker(position: THREE.Vector3, scale: number = 1): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)

  // Corona glow (outer)
  const coronaGeo = new THREE.SphereGeometry(1.2 * scale, 32, 32)
  const coronaMat = new THREE.MeshBasicMaterial({
    color: 0xffdd44,
    transparent: true,
    opacity: 0.2,
  })
  group.add(new THREE.Mesh(coronaGeo, coronaMat))

  // Main sun body
  const sunGeo = new THREE.SphereGeometry(0.8 * scale, 32, 32)
  const sunMat = new THREE.MeshBasicMaterial({ color: 0xffff00 })
  group.add(new THREE.Mesh(sunGeo, sunMat))

  // Core (white hot center)
  const coreGeo = new THREE.SphereGeometry(0.3 * scale, 32, 32)
  const coreMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
  group.add(new THREE.Mesh(coreGeo, coreMat))

  // Sun label
  const label = createTextSprite('☉ Sun', { color: '#ffff00', fontSize: 36 })
  label.position.set(0, 2 * scale, 0)
  label.scale.set(6, 3, 1)
  group.add(label)

  return group
}

/**
 * Create the Earth marker (blue sphere with green patches)
 */
export function createEarthMarker(position: THREE.Vector3, scale: number = 1): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)

  // Atmosphere glow
  const atmoGeo = new THREE.SphereGeometry(1.0 * scale, 32, 32)
  const atmoMat = new THREE.MeshBasicMaterial({
    color: 0x88ccff,
    transparent: true,
    opacity: 0.15,
  })
  group.add(new THREE.Mesh(atmoGeo, atmoMat))

  // Ocean (blue sphere)
  const oceanGeo = new THREE.SphereGeometry(0.7 * scale, 32, 32)
  const oceanMat = new THREE.MeshBasicMaterial({ color: 0x2266aa })
  group.add(new THREE.Mesh(oceanGeo, oceanMat))

  // Land patches (simplified)
  const landGeo = new THREE.SphereGeometry(0.72 * scale, 32, 32, 0, Math.PI * 0.8, Math.PI * 0.3, Math.PI * 0.4)
  const landMat = new THREE.MeshBasicMaterial({ color: 0x228822 })
  const land = new THREE.Mesh(landGeo, landMat)
  land.rotation.y = Math.PI / 4
  group.add(land)

  // Earth label
  const label = createTextSprite('⊕ Earth/ALMA', { color: '#88ccff', fontSize: 32 })
  label.position.set(0, 2 * scale, 0)
  label.scale.set(6, 3, 1)
  group.add(label)

  return group
}

/**
 * Create the Galactic Center marker (Sgr A*)
 */
export function createGalacticCenterMarker(position: THREE.Vector3, scale: number = 1): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)

  // Accretion disk representation (flat ring)
  const diskGeo = new THREE.RingGeometry(0.8 * scale, 2.0 * scale, 32)
  const diskMat = new THREE.MeshBasicMaterial({
    color: 0xff6600,
    transparent: true,
    opacity: 0.4,
    side: THREE.DoubleSide,
  })
  const disk = new THREE.Mesh(diskGeo, diskMat)
  disk.rotation.x = Math.PI / 2  // Lay flat in x-y plane
  group.add(disk)

  // Inner hot region
  const innerGeo = new THREE.RingGeometry(0.3 * scale, 0.8 * scale, 32)
  const innerMat = new THREE.MeshBasicMaterial({
    color: 0xffaa00,
    transparent: true,
    opacity: 0.6,
    side: THREE.DoubleSide,
  })
  const inner = new THREE.Mesh(innerGeo, innerMat)
  inner.rotation.x = Math.PI / 2
  group.add(inner)

  // Central black hole
  const bhGeo = new THREE.SphereGeometry(0.3 * scale, 32, 32)
  const bhMat = new THREE.MeshBasicMaterial({ color: 0x000000 })
  group.add(new THREE.Mesh(bhGeo, bhMat))

  // Event horizon glow
  const glowGeo = new THREE.SphereGeometry(0.35 * scale, 32, 32)
  const glowMat = new THREE.MeshBasicMaterial({
    color: 0xff4400,
    transparent: true,
    opacity: 0.5,
  })
  group.add(new THREE.Mesh(glowGeo, glowMat))

  // Sgr A* label
  const label = createTextSprite('Sgr A* (GC)', { color: '#ff8800', fontSize: 36 })
  label.position.set(0, 3 * scale, 0)
  label.scale.set(6, 3, 1)
  group.add(label)

  return group
}

/**
 * Create the HVCC (High-Velocity Compact Cloud) marker
 * CO-0.40-0.22 from Oka et al. (2017)
 * 
 * Physical parameters:
 * - HCN clump half-size: S = 0.15 pc → diameter ~0.3 pc
 * - Dense clump mass: ~40 M☉
 * - Number density: n(H₂) ~ 10^6.5 cm⁻³
 * - Temperature: T_k = 60 K
 * 
 * @param position Position of the marker
 * @param size Half-size in pc (default 0.15 pc from Oka et al.)
 * @param scale Display scale factor
 */
export function createHVCCMarker(
  position: THREE.Vector3,
  size: number = 0.15,  // HCN clump half-size S = 0.15 pc
  scale: number = 1
): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)
  
  const displaySize = size * scale

  // Molecular cloud representation (fuzzy sphere showing the HCN clump)
  const cloudGeo = new THREE.SphereGeometry(displaySize, 24, 24)
  const cloudMat = new THREE.MeshBasicMaterial({
    color: 0x4488ff,
    transparent: true,
    opacity: 0.3,
  })
  group.add(new THREE.Mesh(cloudGeo, cloudMat))

  // Denser core (innermost region)
  const coreGeo = new THREE.SphereGeometry(displaySize * 0.5, 24, 24)
  const coreMat = new THREE.MeshBasicMaterial({
    color: 0x6699ff,
    transparent: true,
    opacity: 0.5,
  })
  group.add(new THREE.Mesh(coreGeo, coreMat))

  // Cloud boundary (wireframe showing extent)
  const wireGeo = new THREE.SphereGeometry(displaySize, 16, 12)
  const wireMat = new THREE.MeshBasicMaterial({
    color: 0x88aaff,
    wireframe: true,
    transparent: true,
    opacity: 0.4,
  })
  group.add(new THREE.Mesh(wireGeo, wireMat))

  // HVCC label
  const label = createTextSprite('Dense Clump (HCN)', { color: '#88aaff', fontSize: 24 })
  label.position.set(0, displaySize + 0.8 * scale, 0)
  label.scale.set(5, 2.5, 1)
  group.add(label)

  // Sublabel with name and size
  const subLabel = createTextSprite(`CO-0.40-0.22 (S=${size} pc)`, { color: '#aaccff', fontSize: 18 })
  subLabel.position.set(0, displaySize + 0.3 * scale, 0)
  subLabel.scale.set(5, 2.5, 1)
  group.add(subLabel)

  return group
}

/**
 * Create the Line of Sight arrow and visualization
 * 
 * From Oka et al. (2017):
 * - Inclination i = 70° (angle between orbital plane normal and LoS)
 * - Position Angle PA = 41.6° (rotation of orbit projection on sky)
 * - Distance to observer: 8 kpc (shown symbolically, not to scale)
 * 
 * COORDINATE SYSTEM (Oka et al. 2017 Methods):
 * =============================================
 * The Oka paper uses a specific coordinate setup:
 * 1. Y-axis was originally set parallel to the line of sight
 * 2. Then XY plane was rotated about Z-axis by 45°
 * 3. Position angle PA = 41.6° (on the sky plane)
 * 4. Inclination i = 70° (between orbital plane normal and LoS)
 * 
 * GALACTIC COORDINATES at HVCC location (l=-0.4°, b=-0.2°):
 * - The HVCC is very close to the Galactic Center (~60 pc)
 * - Line of sight from HVCC to Earth is approximately toward l=180° from GC
 * - Since HVCC is at l=-0.4°, the LoS is roughly opposite, toward the Sun
 * 
 * For visualization:
 * - We show an arrow pointing toward Earth at a symbolic distance
 * - The orbital plane is tilted by inclination angle from LoS
 * - Position angle rotates the projection on the sky
 * 
 * @param origin Position of the BH (origin of coordinate system)
 * @param config Configuration with viewing angles
 * @param length Display length of LoS arrow (symbolic, not physical 8 kpc)
 */
export function createLineOfSight(
  origin: THREE.Vector3,
  config: GalacticMarkersConfig,
  length: number = 20
): THREE.Group {
  const group = new THREE.Group()

  // Convert viewing angles to direction vector
  // Inclination: angle from z-axis (LoS direction)
  // Position angle: rotation around z-axis
  const inc = (config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination!) * Math.PI / 180
  const pa = (config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle!) * Math.PI / 180

  // Line of sight direction (from system toward observer)
  // i=70° means LoS is 70° from vertical (orbital plane normal)
  // PA=41.6° rotates the projection on the sky
  const losDir = new THREE.Vector3(
    Math.sin(inc) * Math.cos(pa),
    Math.sin(inc) * Math.sin(pa),
    Math.cos(inc)
  ).normalize()

  // Arrow from origin toward observer (scaled for visibility)
  const arrowScale = length / 20  // Relative scale
  const losEnd = origin.clone().add(losDir.clone().multiplyScalar(length))
  const losArrow = createArrow(origin, losEnd, 0xff44ff, {
    headLength: 1.5 * arrowScale,
    headWidth: 0.6 * arrowScale,
    shaftRadius: 0.08 * arrowScale,
    opacity: 0.8,
  })
  group.add(losArrow)

  // LoS label
  const label = createTextSprite('Line of Sight → Earth', { color: '#ff88ff', fontSize: 22 })
  label.position.copy(losDir.clone().multiplyScalar(length * 0.6).add(new THREE.Vector3(2 * arrowScale, 0, 1 * arrowScale)))
  label.scale.set(5, 2.5, 1)
  group.add(label)

  // Direction indicator with angles from Oka et al.
  const angleLabel = createTextSprite(`i=${config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination}° PA=${config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle}°`, { 
    color: '#ffaaff', 
    fontSize: 18 
  })
  angleLabel.position.copy(losDir.clone().multiplyScalar(length * 0.35).add(new THREE.Vector3(2 * arrowScale, -0.5 * arrowScale, 0.5 * arrowScale)))
  angleLabel.scale.set(4, 2, 1)
  group.add(angleLabel)

  return group
}

/**
 * Create the Local Standard of Rest velocity indicator
 * 
 * V_LSR is the radial velocity of the HVCC relative to the Local Standard of Rest.
 * The LSR is a reference frame centered on the Sun but moving with the average
 * velocity of stars in the solar neighborhood.
 * 
 * V_LSR direction:
 * - NEGATIVE V_LSR = cloud is APPROACHING Earth (blueshifted)
 * - POSITIVE V_LSR = cloud is RECEDING from Earth (redshifted)
 * 
 * The arrow should point FROM the HVCC TOWARD Earth for negative V_LSR
 * (indicating the cloud's motion toward the observer).
 */
export function createLSRVelocity(
  position: THREE.Vector3,
  config: GalacticMarkersConfig,
  scale: number = 1
): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)
  group.name = 'vlsrIndicator'

  const lsrVelocity = config.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity!
  
  // V_LSR arrow direction:
  // - Negative V_LSR means cloud is moving TOWARD observer (approaching)
  // - The arrow shows the velocity vector of the cloud relative to LSR
  // - For V_LSR = -120 km/s, the cloud is approaching, so arrow points
  //   from cloud position toward Earth (but we draw it at Earth position
  //   pointing in the direction the cloud is moving relative to observer)
  //
  // Since we're drawing at Earth position, the arrow should point
  // in the direction we would SEE the cloud moving:
  // - Negative V_LSR → arrow points TOWARD us (into the screen / -Z generally)
  // - At Earth, this represents "the cloud is coming this way"
  //
  // Simplified: arrow points along -X (toward GC where cloud is)
  // with length proportional to velocity
  
  const velMag = Math.abs(lsrVelocity) / 25  // Scale factor for display
  
  // Arrow direction: for negative V_LSR (approaching), 
  // show arrow pointing toward Earth (the observer sees it coming)
  // We represent this as an arrow at Earth pointing "inward" from cloud direction
  const arrowDir = new THREE.Vector3(-1, 0, 0)  // From cloud (toward GC) toward Earth
  if (lsrVelocity > 0) {
    arrowDir.negate()  // Receding: away from Earth
  }
  
  const arrowEnd = arrowDir.multiplyScalar(velMag * scale)
  const arrow = createArrow(new THREE.Vector3(0, 0, 0), arrowEnd, 0xff6688, {
    headLength: 0.6 * scale,
    headWidth: 0.35 * scale,
    shaftRadius: 0.06 * scale,
    opacity: 0.9,
  })
  group.add(arrow)

  // V_LSR label with value
  const label = createTextSprite(`V_LSR = ${lsrVelocity} km/s`, { 
    color: '#ff8899', 
    fontSize: 20 
  })
  label.position.set(0, -0.8 * scale, 0)
  label.scale.set(5, 2.5, 1)
  group.add(label)

  // Approach/Recede indicator with explanation
  const motionText = lsrVelocity < 0 
    ? '← approaching (blueshift)' 
    : '→ receding (redshift)'
  const motionLabel = createTextSprite(motionText, {
    color: lsrVelocity < 0 ? '#88ccff' : '#ffaa88',
    fontSize: 16,
  })
  motionLabel.position.set(0, -1.5 * scale, 0)
  motionLabel.scale.set(5, 2, 1)
  group.add(motionLabel)

  return group
}

/**
 * Create local Galactic coordinate axes
 */
export function createGalacticAxes(
  origin: THREE.Vector3,
  size: number = 15
): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(origin)

  // X-axis: Toward Galactic Center (red)
  const xArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(size, 0, 0),
    0xff4444,
    { headLength: 1, headWidth: 0.4 }
  )
  group.add(xArrow)
  
  const xLabel = createTextSprite('→ GC (l=0°)', { color: '#ff6666', fontSize: 24 })
  xLabel.position.set(size + 2, 0, 0)
  xLabel.scale.set(5, 2.5, 1)
  group.add(xLabel)

  // Y-axis: Galactic rotation direction (green)
  const yArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, size, 0),
    0x44ff44,
    { headLength: 1, headWidth: 0.4 }
  )
  group.add(yArrow)
  
  const yLabel = createTextSprite('→ l=90°', { color: '#66ff66', fontSize: 24 })
  yLabel.position.set(0, size + 2, 0)
  yLabel.scale.set(5, 2.5, 1)
  group.add(yLabel)

  // Z-axis: Galactic North Pole (blue)
  const zArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, 0, size),
    0x4488ff,
    { headLength: 1, headWidth: 0.4 }
  )
  group.add(zArrow)
  
  const zLabel = createTextSprite('→ NGP (b=90°)', { color: '#6699ff', fontSize: 24 })
  zLabel.position.set(0, 0, size + 2)
  zLabel.scale.set(5, 2.5, 1)
  group.add(zLabel)

  // Origin label
  const originLabel = createTextSprite('Local Frame', { color: '#ffffff', fontSize: 28 })
  originLabel.position.set(0, -2, 0)
  originLabel.scale.set(5, 2.5, 1)
  group.add(originLabel)

  return group
}

/**
 * Create the Galactic plane reference
 */
export function createGalacticPlane(
  center: THREE.Vector3,
  radius: number = 30
): THREE.Group {
  const group = new THREE.Group()

  // Semi-transparent galactic plane
  const planeGeo = new THREE.CircleGeometry(radius, 64)
  const planeMat = new THREE.MeshBasicMaterial({
    color: 0x334488,
    transparent: true,
    opacity: 0.1,
    side: THREE.DoubleSide,
  })
  const plane = new THREE.Mesh(planeGeo, planeMat)
  plane.rotation.x = -Math.PI / 2  // Lay flat
  plane.position.copy(center)
  group.add(plane)

  // Plane boundary circle
  const boundary = createDashedCircle(radius, 0x4466aa, {
    segments: 64,
    dashSize: 1,
    gapSize: 0.5,
    opacity: 0.4,
    normal: new THREE.Vector3(0, 1, 0),
    center: center,
  })
  group.add(boundary)

  // "Galactic Plane" label
  const label = createTextSprite('Galactic Plane (b=0°)', { color: '#6688bb', fontSize: 24 })
  label.position.copy(center.clone().add(new THREE.Vector3(radius - 5, 0.5, 0)))
  label.scale.set(6, 3, 1)
  group.add(label)

  return group
}

/**
 * Create a rotating Milky Way disk visualization (schematic)
 * 
 * The Milky Way is a barred spiral galaxy:
 * - Disk diameter: ~100,000 ly (~30 kpc)
 * - Central bulge: ~10,000 ly diameter
 * - Sun is ~8 kpc from center
 * - HVCC/BH system is ~60 pc from center (in inner region)
 * 
 * For visualization, we show this at a very reduced scale with transparency
 * to give context without obscuring the simulation.
 * 
 * @param center Center of the galaxy (Galactic Center)
 * @param scale Display scale (1 = 1 pc per scene unit)
 * @param config Galaxy visualization options
 */
export function createGalaxyDisk(
  center: THREE.Vector3,
  scale: number = 1.0,
  config: {
    showSpiralArms?: boolean
    opacity?: number
    rotationSpeed?: number  // rad/s for animation
  } = {}
): THREE.Group {
  const {
    showSpiralArms = true,
    opacity = 0.4,  // Cloud-like opacity for visibility
  } = config
  
  const group = new THREE.Group()
  group.position.copy(center)
  group.name = 'galaxyDisk'
  
  // Galactic parameters (scaled down for visualization)
  // TRUE SCALE REFERENCE:
  //   - HVCC/BH system is ~60 pc from Galactic Center
  //   - Sun is 8 kpc (8000 pc) from GC, i.e., ~133× farther than HVCC
  //   - True disk radius: ~15 kpc
  //
  // For visualization, we use a SYMBOLIC representation:
  //   - Scene is centered on BH/HVCC system (which is ~60 pc from GC)
  //   - We show GC nearby at a display distance (NOT true 60 pc)
  //   - Sun is shown at a symbolic display distance along its galactic direction
  //   - The galaxy disk is purely illustrative, not to scale
  //
  // Compression: 1 display unit ≈ 300 pc (so 8 kpc Sun ≈ 27 units from GC marker)
  const displayScale = 1.0 / 300  // 1 display unit = 300 pc
  const diskRadius = 15000 * displayScale * scale  // ~50 display units
  const bulgeRadius = 2500 * displayScale * scale  // ~8 display units
  // Sun distance: 8 kpc = 8000 pc → ~27 display units from GC
  // NOTE: This is a SYMBOLIC representation for context, not true scale
  const sunDistance = 8000 * displayScale * scale  // 8 kpc → ~27 display units
  
  // Create the main disk (circular gradient)
  const diskSegments = 128
  const diskGeo = new THREE.CircleGeometry(diskRadius, diskSegments)
  
  // Create vertex colors for disk gradient
  const diskColors = new Float32Array((diskSegments + 2) * 3)
  // Center is brighter
  diskColors[0] = 1.0
  diskColors[1] = 0.95
  diskColors[2] = 0.8
  // Outer vertices fade to dark
  for (let i = 1; i <= diskSegments + 1; i++) {
    const t = 1 - (i / diskSegments) * 0.3  // Gradual fade
    diskColors[i * 3] = 0.2 * t
    diskColors[i * 3 + 1] = 0.25 * t
    diskColors[i * 3 + 2] = 0.4 * t
  }
  diskGeo.setAttribute('color', new THREE.BufferAttribute(diskColors, 3))
  
  const diskMat = new THREE.MeshBasicMaterial({
    vertexColors: true,
    transparent: true,
    opacity: opacity,
    side: THREE.DoubleSide,
  })
  const disk = new THREE.Mesh(diskGeo, diskMat)
  disk.rotation.x = -Math.PI / 2  // Lay flat in galactic plane
  group.add(disk)
  
  // Central bulge (brighter ellipsoid-like structure)
  const bulgeGeo = new THREE.CircleGeometry(bulgeRadius, 64)
  const bulgeMat = new THREE.MeshBasicMaterial({
    color: 0xffeecc,
    transparent: true,
    opacity: opacity * 2,
    side: THREE.DoubleSide,
  })
  const bulge = new THREE.Mesh(bulgeGeo, bulgeMat)
  bulge.rotation.x = -Math.PI / 2
  bulge.position.y = 0.01  // Slightly above disk
  group.add(bulge)
  
  // Spiral arms (simplified log-spiral pattern)
  if (showSpiralArms) {
    const armColor = 0x4488cc
    const numArms = 4  // Milky Way has ~4 major arms
    const armTurns = 1.5
    const armPoints = 100
    
    for (let arm = 0; arm < numArms; arm++) {
      const armOffset = (arm / numArms) * Math.PI * 2
      const points: THREE.Vector3[] = []
      
      for (let i = 0; i < armPoints; i++) {
        const t = i / armPoints
        const r = (bulgeRadius + t * (diskRadius - bulgeRadius)) * 0.9
        const theta = armOffset + t * armTurns * Math.PI * 2
        
        // Log-spiral with some variation
        const x = r * Math.cos(theta)
        const z = r * Math.sin(theta)
        points.push(new THREE.Vector3(x, 0.02, z))
      }
      
      const armGeo = new THREE.BufferGeometry().setFromPoints(points)
      const armMat = new THREE.LineBasicMaterial({
        color: armColor,
        transparent: true,
        opacity: opacity * 1.5,
      })
      const armLine = new THREE.Line(armGeo, armMat)
      group.add(armLine)
    }
  }
  
  // Galactic Center marker (Sgr A*) at origin of galaxy
  const gcGeo = new THREE.SphereGeometry(0.5 * scale, 32, 32)
  const gcMat = new THREE.MeshBasicMaterial({
    color: 0xffaa00,
    transparent: true,
    opacity: 0.3,
  })
  const gcMarker = new THREE.Mesh(gcGeo, gcMat)
  group.add(gcMarker)
  
  // Sun position indicator (8 kpc from center, toward negative X in our coords)
  // In galactic coordinates looking down from NGP:
  // - GC is at origin  
  // - Sun is toward l=180° from GC's perspective
  const sunPosInGalaxy = new THREE.Vector3(-sunDistance, 0.05, 0)
  
  // NOTE: We only show Sun position as a small marker here.
  // The detailed Sun/Earth markers are created in createAllGalacticMarkers()
  // to avoid duplication.
  
  // Small Sun position marker (just a dot to show location)
  const sunDotGeo = new THREE.SphereGeometry(0.4 * scale, 16, 16)
  const sunDotMat = new THREE.MeshBasicMaterial({
    color: 0xffff44,
    transparent: true,
    opacity: 0.4,
  })
  const sunDot = new THREE.Mesh(sunDotGeo, sunDotMat)
  sunDot.position.copy(sunPosInGalaxy)
  group.add(sunDot)
  
  // Small label showing Sun's galactic position
  const sunPosLabel = createTextSprite('☉ (8 kpc)', { color: '#ffff88', fontSize: 18 })
  sunPosLabel.position.copy(sunPosInGalaxy.clone().add(new THREE.Vector3(0, 1.2 * scale, 0)))
  sunPosLabel.scale.set(4, 2, 1)
  group.add(sunPosLabel)
  
  // Dashed circle showing Sun's orbit around GC
  const sunOrbit = createDashedCircle(sunDistance, 0xffff44, {
    segments: 64,
    dashSize: 1.0,
    gapSize: 0.5,
    opacity: opacity * 2,
    normal: new THREE.Vector3(0, 1, 0),
    center: new THREE.Vector3(0, 0, 0),
  })
  group.add(sunOrbit)
  
  // HVCC/BH location marker (very close to GC at ~60 pc)
  // With compression: 60 pc / 300 = 0.2 display units from GC center
  const hvccDistance = 60 * displayScale * scale
  const hvccDir = new THREE.Vector3(0.5, 0, -0.866)  // Approximate position toward l=-0.4°
  const hvccPosInGalaxy = hvccDir.multiplyScalar(hvccDistance)
  
  const hvccIndicator = new THREE.SphereGeometry(0.3 * scale, 16, 16)
  const hvccMat = new THREE.MeshBasicMaterial({
    color: 0x00ff88,
    transparent: true,
    opacity: 0.5,
  })
  const hvccMarker = new THREE.Mesh(hvccIndicator, hvccMat)
  hvccMarker.position.copy(hvccPosInGalaxy)
  group.add(hvccMarker)
  
  // HVCC label
  const hvccLabel = createTextSprite('HVCC (~60pc from GC)', { color: '#00ff88', fontSize: 16 })
  hvccLabel.position.copy(hvccPosInGalaxy.clone().add(new THREE.Vector3(1.5 * scale, 0.5 * scale, 0)))
  hvccLabel.scale.set(5, 2, 1)
  group.add(hvccLabel)
  
  // Scale legend showing the compression
  const legendPos = new THREE.Vector3(diskRadius * 0.7, 0.1, diskRadius * 0.5)
  const scaleNote = createTextSprite('Scale: Compressed ~300× for visualization', { color: '#888888', fontSize: 12 })
  scaleNote.position.copy(legendPos)
  scaleNote.scale.set(5, 2, 1)
  group.add(scaleNote)
  
  // Distance scale bar
  const distNote1 = createTextSprite('• Sun ↔ GC: 8 kpc (~26.7 display units)', { color: '#ffff88', fontSize: 11 })
  distNote1.position.copy(legendPos.clone().add(new THREE.Vector3(0, -1 * scale, 0)))
  distNote1.scale.set(6, 2, 1)
  group.add(distNote1)
  
  const distNote2 = createTextSprite('• HVCC ↔ GC: 60 pc (~0.2 display units)', { color: '#00ff88', fontSize: 11 })
  distNote2.position.copy(legendPos.clone().add(new THREE.Vector3(0, -2 * scale, 0)))
  distNote2.scale.set(6, 2, 1)
  group.add(distNote2)
  
  // Galaxy label
  const galaxyLabel = createTextSprite('Milky Way (schematic)', { color: '#6688aa', fontSize: 20 })
  galaxyLabel.position.set(0, 3 * scale, diskRadius - 5)
  galaxyLabel.scale.set(6, 3, 1)
  group.add(galaxyLabel)
  
  // Store reference for animation updates
  group.userData = {
    disk,
    bulge,
    rotationSpeed: config.rotationSpeed || 0,
  }
  
  return group
}

/**
 * Create a simple Sun + Earth orbit representation with LOCAL COORDINATE SYSTEM
 * 
 * This shows the Solar System with its local galactic coordinate axes:
 * - Sun at center (yellow)
 * - Earth orbit circle with Earth marker
 * - LOCAL coordinate axes:
 *   - X: toward Galactic Center (l=0°)
 *   - Y: direction of galactic rotation (l=90°)  
 *   - Z: toward North Galactic Pole
 * - V_LSR arrow attached to Earth position in local coordinates
 * 
 * At galactic scales the solar system is invisible (~0.0003 pc),
 * so we enlarge it ~10^6 times for visibility.
 * 
 * Earth's orbit (1 year) is much faster than galactic rotation (220 Myr),
 * so we don't animate Earth's orbit - it would be a blur.
 * 
 * @param position Position in the scene (Sun's galactic position)
 * @param scale Display scale factor
 * @param initialEarthAngle Initial angle of Earth on orbit (radians)
 * @param lsrVelocity V_LSR in km/s (negative = source approaching)
 * @param losDirection Line of sight direction to the source (normalized)
 */
export function createSolarSystemDemo(
  position: THREE.Vector3,
  scale: number = 1.0,
  initialEarthAngle: number = 0,
  lsrVelocity: number = -60,
  losDirection?: THREE.Vector3
): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)
  group.name = 'solarSystemDemo'
  
  // Scale: 1 AU = 1.5 scene units for visibility
  const auScale = 1.5 * scale
  const earthOrbitRadius = 1.0 * auScale
  
  // ========== LOCAL GALACTIC COORDINATE AXES ==========
  // These show the standard galactic coordinate frame at Sun's position:
  // X: toward Galactic Center (l=0°, b=0°)
  // Y: direction of galactic rotation (l=90°, b=0°)
  // Z: toward North Galactic Pole (b=90°)
  const axisLength = 2.5 * scale
  const axisOpacity = 0.7
  
  // X-axis: toward GC (red)
  const xAxisGeo = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(axisLength, 0, 0)
  ])
  const xAxisMat = new THREE.LineBasicMaterial({ color: 0xff4444, transparent: true, opacity: axisOpacity })
  group.add(new THREE.Line(xAxisGeo, xAxisMat))
  
  // X-axis arrow head
  const xArrow = new THREE.ConeGeometry(0.12 * scale, 0.3 * scale, 8)
  const xArrowMesh = new THREE.Mesh(xArrow, new THREE.MeshBasicMaterial({ color: 0xff4444 }))
  xArrowMesh.position.set(axisLength, 0, 0)
  xArrowMesh.rotation.z = -Math.PI / 2
  group.add(xArrowMesh)
  
  const xLabel = createTextSprite('→ GC (l=0°)', { color: '#ff6666', fontSize: 14 })
  xLabel.position.set(axisLength + 0.8 * scale, 0.3 * scale, 0)
  xLabel.scale.set(3, 1.5, 1)
  group.add(xLabel)
  
  // Y-axis: galactic rotation direction (green)
  const yAxisGeo = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, 0, axisLength)  // Z in THREE = Y galactic (rotation dir)
  ])
  const yAxisMat = new THREE.LineBasicMaterial({ color: 0x44ff44, transparent: true, opacity: axisOpacity })
  group.add(new THREE.Line(yAxisGeo, yAxisMat))
  
  // Y-axis arrow head
  const yArrow = new THREE.ConeGeometry(0.12 * scale, 0.3 * scale, 8)
  const yArrowMesh = new THREE.Mesh(yArrow, new THREE.MeshBasicMaterial({ color: 0x44ff44 }))
  yArrowMesh.position.set(0, 0, axisLength)
  yArrowMesh.rotation.x = Math.PI / 2
  group.add(yArrowMesh)
  
  const yLabel = createTextSprite('V_rot (l=90°)', { color: '#66ff66', fontSize: 14 })
  yLabel.position.set(0, 0.3 * scale, axisLength + 0.8 * scale)
  yLabel.scale.set(3.2, 1.5, 1)
  group.add(yLabel)
  
  // Z-axis: toward NGP (blue)  
  const zAxisGeo = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, axisLength, 0)  // Y in THREE = Z galactic (NGP)
  ])
  const zAxisMat = new THREE.LineBasicMaterial({ color: 0x4488ff, transparent: true, opacity: axisOpacity })
  group.add(new THREE.Line(zAxisGeo, zAxisMat))
  
  // Z-axis arrow head
  const zArrow = new THREE.ConeGeometry(0.12 * scale, 0.3 * scale, 8)
  const zArrowMesh = new THREE.Mesh(zArrow, new THREE.MeshBasicMaterial({ color: 0x4488ff }))
  zArrowMesh.position.set(0, axisLength, 0)
  group.add(zArrowMesh)
  
  const zLabel = createTextSprite('NGP (b=90°)', { color: '#66aaff', fontSize: 14 })
  zLabel.position.set(0.5 * scale, axisLength + 0.5 * scale, 0)
  zLabel.scale.set(3, 1.5, 1)
  group.add(zLabel)
  
  // ========== Sun ==========
  const sunGeo = new THREE.SphereGeometry(0.4 * scale, 32, 32)
  const sunMat = new THREE.MeshBasicMaterial({ color: 0xffff00 })
  const sun = new THREE.Mesh(sunGeo, sunMat)
  sun.name = 'sun'
  group.add(sun)
  
  // Sun glow
  const glowGeo = new THREE.SphereGeometry(0.6 * scale, 32, 32)
  const glowMat = new THREE.MeshBasicMaterial({
    color: 0xffee66,
    transparent: true,
    opacity: 0.25,
  })
  group.add(new THREE.Mesh(glowGeo, glowMat))
  
  // Sun label
  const sunLabel = createTextSprite('☉ Sun (LSR)', { color: '#ffff88', fontSize: 18 })
  sunLabel.position.set(0, -0.8 * scale, 0)
  sunLabel.scale.set(3.5, 1.8, 1)
  group.add(sunLabel)
  
  // ========== Earth Orbit ==========
  const orbitPoints: THREE.Vector3[] = []
  for (let i = 0; i <= 64; i++) {
    const theta = (i / 64) * Math.PI * 2
    orbitPoints.push(new THREE.Vector3(
      earthOrbitRadius * Math.cos(theta),
      0,
      earthOrbitRadius * Math.sin(theta)
    ))
  }
  const orbitGeo = new THREE.BufferGeometry().setFromPoints(orbitPoints)
  const orbitMat = new THREE.LineBasicMaterial({
    color: 0x4488cc,
    transparent: true,
    opacity: 0.5,
  })
  const orbit = new THREE.Line(orbitGeo, orbitMat)
  group.add(orbit)
  
  // ========== Earth ==========
  const earthAngle = initialEarthAngle
  const earthGeo = new THREE.SphereGeometry(0.18 * scale, 16, 16)
  const earthMat = new THREE.MeshBasicMaterial({ color: 0x2288dd })
  const earth = new THREE.Mesh(earthGeo, earthMat)
  earth.name = 'earth'
  earth.position.set(
    earthOrbitRadius * Math.cos(earthAngle),
    0,
    earthOrbitRadius * Math.sin(earthAngle)
  )
  group.add(earth)
  
  // Earth label
  const earthLabel = createTextSprite('🌍 Earth', { color: '#88ccff', fontSize: 16 })
  earthLabel.position.copy(earth.position.clone().add(new THREE.Vector3(0, 0.5 * scale, 0)))
  earthLabel.scale.set(3.5, 1.8, 1)
  earthLabel.name = 'earthLabel'
  group.add(earthLabel)
  
  // ========== V_LSR Arrow (attached to Earth) ==========
  // V_LSR = radial velocity of source relative to Local Standard of Rest
  // Negative V_LSR = source approaching Earth (blueshift)
  // Arrow points FROM Earth TOWARD source direction (LoS direction)
  // Arrow length scales with |V_LSR|
  const lsrArrowLength = Math.min(Math.abs(lsrVelocity) / 50, 3) * scale  // Scale by velocity magnitude
  
  // Line of sight direction (default: toward -X if not specified)
  // For HVCC: it's along the viewing geometry angles
  const losDir = losDirection ? losDirection.clone().normalize() : new THREE.Vector3(-0.7, 0.3, -0.5).normalize()
  
  // Arrow starts at Earth position
  const lsrStart = earth.position.clone()
  // If V_LSR < 0 (approaching), arrow points toward source (along +LoS from source to Earth, so -LoS from Earth)
  // Actually V_LSR is defined as positive = receding, negative = approaching
  // For display: arrow shows the velocity vector of the source relative to LSR
  // Negative V_LSR means source is moving TOWARD us, so we show arrow pointing toward Earth (from source direction)
  const lsrDir = lsrVelocity < 0 ? losDir.clone().negate() : losDir.clone()
  const lsrEnd = lsrStart.clone().add(lsrDir.multiplyScalar(lsrArrowLength))
  
  // Create V_LSR arrow
  const lsrLineGeo = new THREE.BufferGeometry().setFromPoints([lsrStart, lsrEnd])
  const lsrColor = lsrVelocity < 0 ? 0xff66ff : 0xff9966  // Magenta for approaching, orange for receding
  const lsrLineMat = new THREE.LineBasicMaterial({ color: lsrColor, linewidth: 3 })
  group.add(new THREE.Line(lsrLineGeo, lsrLineMat))
  
  // Arrow head for V_LSR
  const lsrArrowHead = new THREE.ConeGeometry(0.15 * scale, 0.35 * scale, 8)
  const lsrArrowMesh = new THREE.Mesh(lsrArrowHead, new THREE.MeshBasicMaterial({ color: lsrColor }))
  lsrArrowMesh.position.copy(lsrEnd)
  // Orient cone along arrow direction
  const arrowQuat = new THREE.Quaternion()
  arrowQuat.setFromUnitVectors(new THREE.Vector3(0, 1, 0), lsrDir.clone().normalize())
  lsrArrowMesh.setRotationFromQuaternion(arrowQuat)
  group.add(lsrArrowMesh)
  
  // V_LSR label
  const lsrLabel = createTextSprite(`V_LSR = ${lsrVelocity} km/s`, { color: '#ff88ff', fontSize: 16 })
  const labelOffset = lsrDir.clone().normalize().multiplyScalar(lsrArrowLength * 0.5)
  lsrLabel.position.copy(lsrStart.clone().add(labelOffset).add(new THREE.Vector3(0, 0.6 * scale, 0)))
  lsrLabel.scale.set(4, 2, 1)
  group.add(lsrLabel)
  
  // Approaching/receding indicator
  const motionLabel = createTextSprite(
    lsrVelocity < 0 ? '(approaching)' : '(receding)', 
    { color: lsrVelocity < 0 ? '#88ffff' : '#ffaa88', fontSize: 12 }
  )
  motionLabel.position.copy(lsrLabel.position.clone().add(new THREE.Vector3(0, -0.5 * scale, 0)))
  motionLabel.scale.set(3, 1.5, 1)
  group.add(motionLabel)
  
  // Store Earth position for external reference
  group.userData = {
    earthPosition: earth.position.clone(),
    earthOrbitRadius,
    earthAngle,
    lsrVelocity,
  }
  
  return group
}

/**
 * Update galaxy disk rotation for animation
 * @param galaxyGroup The group returned by createGalaxyDisk
 * @param deltaTime Time elapsed since last frame (seconds)
 */
export function updateGalaxyRotation(
  galaxyGroup: THREE.Group,
  deltaTime: number
): void {
  const { rotationSpeed } = galaxyGroup.userData || {}
  if (rotationSpeed) {
    // The galaxy disk is laid flat with rotation.x = -PI/2
    // This means the disk lies in the XZ plane with Y as the vertical axis (NGP direction)
    // So we rotate around Y-axis (the disk's normal after the X rotation)
    // Counter-clockwise when viewed from above (+Y direction, looking down at the disk)
    galaxyGroup.rotation.y += rotationSpeed * deltaTime
  }
}

/**
 * Update the V_circular animation marker in the observer geometry indicator
 * This animates a small marker moving along the Sun's orbital arc to visualize
 * the 220 km/s circular velocity of the LSR around the Galactic Center.
 * 
 * @param observerGroup The group returned by createObserverDirectionIndicator
 * @param deltaTime Time elapsed since last frame (seconds)
 */
export function updateVCircularAnimation(
  observerGroup: THREE.Group,
  deltaTime: number
): void {
  // Find the vCircularIndicator group
  const vCircularGroup = observerGroup.getObjectByName('vCircularIndicator') as THREE.Group | undefined
  if (!vCircularGroup) return
  
  const { orbitRadius, animationSpeed } = vCircularGroup.userData || {}
  if (!orbitRadius) return
  
  // Update animation angle - continuous circular motion (counter-clockwise when viewed from above)
  let angle = (vCircularGroup.userData.animationAngle || 0) + animationSpeed * deltaTime
  
  // Keep angle in [0, 2π] range
  angle = angle % (Math.PI * 2)
  
  vCircularGroup.userData.animationAngle = angle
  
  // Find and update the orbit animation marker - continuous circular motion
  const orbitMarker = vCircularGroup.getObjectByName('orbitAnimationMarker') as THREE.Mesh | undefined
  if (orbitMarker) {
    // Circular motion: x = r*cos(θ), z = r*sin(θ)
    // Start from +X (toward GC) and move counter-clockwise (direction of galactic rotation)
    orbitMarker.position.set(
      orbitRadius * Math.cos(angle),
      0,
      orbitRadius * Math.sin(angle)
    )
  }
  
  // Update trail markers (follow behind the main marker)
  for (let i = 1; i <= 3; i++) {
    const trailMarker = vCircularGroup.getObjectByName(`orbitTrail${i}`) as THREE.Mesh | undefined
    if (trailMarker) {
      const trailAngle = angle - i * 0.15  // Trail behind by small angle
      trailMarker.position.set(
        orbitRadius * Math.cos(trailAngle),
        0,
        orbitRadius * Math.sin(trailAngle)
      )
    }
  }
}

/**
 * Update Solar System / LSR frame rotation around Galactic Center
 * The LSR rotates with the galaxy at the Sun's position (V_circ = 220 km/s)
 * 
 * @param solarSystemGroup The group returned by createSolarSystemDemo
 * @param gcPosition Position of the Galactic Center
 * @param rotationSpeed Rotation speed in rad/s (same as galaxy disk)
 * @param deltaTime Time elapsed since last frame (seconds)
 */
export function updateLSRRotation(
  solarSystemGroup: THREE.Group,
  gcPosition: THREE.Vector3,
  rotationSpeed: number,
  deltaTime: number
): void {
  if (!rotationSpeed) return

  // Rotate the Solar System around the Galactic Center
  // The LSR is at the Sun's position which orbits GC at ~8 kpc
  const angle = rotationSpeed * deltaTime

  // Get current position relative to GC
  const pos = solarSystemGroup.position.clone().sub(gcPosition)

  // The galaxy disk lies in the XZ plane (after rotation.x = -PI/2)
  // So we rotate around Y-axis (the vertical axis, NGP direction)
  // This rotates in the XZ plane
  const cosA = Math.cos(angle)
  const sinA = Math.sin(angle)
  const newX = pos.x * cosA - pos.z * sinA
  const newZ = pos.x * sinA + pos.z * cosA

  // Update position (Y stays constant as we rotate in XZ plane)
  solarSystemGroup.position.set(gcPosition.x + newX, solarSystemGroup.position.y, gcPosition.z + newZ)

  // Also rotate the local coordinate axes to stay aligned with galactic frame
  solarSystemGroup.rotation.y += angle
}

/**
 * Create galactic markers for the IMBH scene
 *
 * Physical scales from Oka et al. (2017):
 * - Scene units are in parsecs (pc)
 * - BH is at origin (or bhPosition)
 * - HVCC cloud is represented by the simulation particles
 * - Distance to GC: ~60 pc
 * - Observer (Earth) is at 8 kpc from GC (shown symbolically)
 * 
 * Observer Geometry (Oka et al. 2017):
 * - Inclination i = 70° (angle between orbital plane normal and line of sight)
 * - Position Angle PA = 41.6° (rotation of orbit projection on sky plane)
 * - V_LSR = -120 km/s (cloud approaching observer)
 */
export function createAllGalacticMarkers(
  bhPosition: THREE.Vector3,
  config: GalacticMarkersConfig = DEFAULT_GALACTIC_CONFIG
): THREE.Group {
  const group = new THREE.Group()

  // Physical scale factors
  const scale = config.displayScale || 1.0
  const distanceToGC_pc = config.distanceToGC_pc ?? DEFAULT_GALACTIC_CONFIG.distanceToGC_pc!

  // ========== Galaxy disk visualization ==========
  // Shows the Milky Way structure with Sun position, GC, and rotation
  if (config.showGalaxyDisk) {
    // Position Galactic Center (it's ~60 pc from BH in +X direction)
    const gcDisplayDist = Math.min(distanceToGC_pc, 25) * scale
    const gcPosition = bhPosition.clone().add(new THREE.Vector3(gcDisplayDist, 0, 0))

    // Galaxy rotation speed:
    // Default: 0.02 rad/s gives ~5 min per full rotation (sped up ~10^10 times)
    const defaultRotationSpeed = 0.02

    // Create the galaxy disk centered on GC
    const galaxyDisk = createGalaxyDisk(gcPosition, scale, {
      showSpiralArms: true,
      opacity: 0.4,
      rotationSpeed: config.galaxyRotationSpeed ?? defaultRotationSpeed,
    })
    galaxyDisk.name = 'galaxyDisk'  // For animation updates
    group.add(galaxyDisk)

    // Add a simple label showing distance to GC
    const gcLabel = createTextSprite(`GC (~${distanceToGC_pc} pc)`, {
      color: '#ffaa44',
      fontSize: 16
    })
    gcLabel.position.copy(gcPosition.clone().add(new THREE.Vector3(0, 3 * scale, 0)))
    gcLabel.scale.set(4, 2, 1)
    group.add(gcLabel)
  }

  // ========== Observer Direction Indicator ==========
  // Shows the line of sight toward Earth with proper inclination and position angle
  // This replaces the solar system demo with a more accurate geometric representation
  if (config.showSolarSystem) {
    const observerGroup = createObserverDirectionIndicator(
      bhPosition,
      config,
      15 * scale  // Symbolic display distance for the indicator
    )
    observerGroup.name = 'observerIndicator'
    group.add(observerGroup)
  }

  return group
}

/**
 * Create an Observer Direction Indicator showing the line of sight geometry
 * 
 * This properly visualizes:
 * 1. Line of sight direction toward Earth (at symbolic distance)
 * 2. Inclination angle (i = 70°) - angle between LoS and orbital plane normal
 * 3. Position angle (PA = 41.6°) - rotation on the sky plane
 * 4. V_LSR direction (negative = approaching)
 * 
 * Coordinate System (following Oka et al. 2017):
 * - Orbital plane is defined by X-Y axes in the simulation
 * - Z is perpendicular to orbital plane (angular momentum direction)
 * - Line of sight makes angle i with Z-axis
 * - Position angle rotates the projected orbit on the sky
 * 
 * @param origin Position of BH (coordinate system origin)
 * @param config Galactic markers configuration
 * @param displayDistance Symbolic distance for the observer marker
 */
export function createObserverDirectionIndicator(
  origin: THREE.Vector3,
  config: GalacticMarkersConfig,
  displayDistance: number = 15
): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(origin)
  group.name = 'observerDirectionIndicator'
  
  const scale = config.displayScale || 1.0
  const inc = ((config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination!) * Math.PI) / 180
  const pa = ((config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle!) * Math.PI) / 180
  const lsrVelocity = config.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity!
  
  // ========== LINE OF SIGHT DIRECTION ==========
  // From Oka et al. (2017) Methods:
  // - Inclination i: angle between orbital plane normal (Z) and line of sight
  // - Position angle PA: rotation of projected orbit on sky plane
  // 
  // The LoS direction in orbital plane coordinates:
  // - i = 0° means face-on (LoS along Z, perpendicular to orbital plane)
  // - i = 90° means edge-on (LoS in XY plane)
  // - i = 70° means we're looking 70° from pole
  //
  // The PA rotates which part of the orbit is "up" on the sky
  // PA is measured from North through East on the celestial sphere
  // In our local frame, we interpret this as rotation around the LoS direction
  
  // Line of sight direction (from BH toward observer/Earth)
  // Using spherical coordinates where:
  // - inc is angle from +Z axis
  // - pa is azimuthal angle from +X axis
  const losDir = new THREE.Vector3(
    Math.sin(inc) * Math.cos(pa),
    Math.sin(inc) * Math.sin(pa),
    Math.cos(inc)
  ).normalize()
  
  // Arrow pointing toward Earth (observer)
  const losEnd = losDir.clone().multiplyScalar(displayDistance)
  const losArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    losEnd,
    0x44aaff,  // Blue for observer direction
    {
      headLength: 1.5 * scale,
      headWidth: 0.6 * scale,
      shaftRadius: 0.1 * scale,
      opacity: 0.9,
    }
  )
  group.add(losArrow)
  
  // Observer marker at the end (Earth symbol)
  const observerMarker = new THREE.Group()
  const earthGeo = new THREE.SphereGeometry(0.8 * scale, 32, 32)
  const earthMat = new THREE.MeshBasicMaterial({ color: 0x2288dd })
  const earthMesh = new THREE.Mesh(earthGeo, earthMat)
  observerMarker.add(earthMesh)
  
  // Earth atmosphere glow
  const atmoGeo = new THREE.SphereGeometry(1.0 * scale, 32, 32)
  const atmoMat = new THREE.MeshBasicMaterial({
    color: 0x88ccff,
    transparent: true,
    opacity: 0.25,
  })
  observerMarker.add(new THREE.Mesh(atmoGeo, atmoMat))
  observerMarker.position.copy(losEnd)
  group.add(observerMarker)
  
  // Observer label
  const observerLabel = createTextSprite('⊕ Observer (Earth)', { color: '#88ccff', fontSize: 20 })
  observerLabel.position.copy(losEnd.clone().add(new THREE.Vector3(0, 2 * scale, 0)))
  observerLabel.scale.set(5, 2.5, 1)
  group.add(observerLabel)
  
  // Distance annotation (not to scale)
  const distLabel = createTextSprite('~8 kpc (not to scale)', { color: '#aaddff', fontSize: 14 })
  distLabel.position.copy(losEnd.clone().add(new THREE.Vector3(0, 1.2 * scale, 0)))
  distLabel.scale.set(4, 2, 1)
  group.add(distLabel)
  
  // ========== INCLINATION ANGLE VISUALIZATION ==========
  // Show the angle between Z-axis (orbital normal) and LoS
  const incArcRadius = displayDistance * 0.3
  const incArcPoints: THREE.Vector3[] = []
  const numArcPoints = 20
  for (let i = 0; i <= numArcPoints; i++) {
    const t = i / numArcPoints
    const angle = t * inc
    // Arc from Z-axis toward LoS direction
    const arcPoint = new THREE.Vector3(
      Math.sin(angle) * Math.cos(pa),
      Math.sin(angle) * Math.sin(pa),
      Math.cos(angle)
    ).multiplyScalar(incArcRadius)
    incArcPoints.push(arcPoint)
  }
  const incArcGeo = new THREE.BufferGeometry().setFromPoints(incArcPoints)
  const incArcMat = new THREE.LineBasicMaterial({
    color: 0xffaa44,
    transparent: true,
    opacity: 0.8,
  })
  group.add(new THREE.Line(incArcGeo, incArcMat))
  
  // Inclination label
  const incLabelPos = new THREE.Vector3(
    Math.sin(inc / 2) * Math.cos(pa),
    Math.sin(inc / 2) * Math.sin(pa),
    Math.cos(inc / 2)
  ).multiplyScalar(incArcRadius + 2 * scale)
  const incLabel = createTextSprite(`i = ${config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination}°`, { 
    color: '#ffcc66', 
    fontSize: 18 
  })
  incLabel.position.copy(incLabelPos)
  incLabel.scale.set(3.5, 1.8, 1)
  group.add(incLabel)
  
  // ========== ORBITAL PLANE NORMAL (Z-axis reference) ==========
  const zAxisLength = displayDistance * 0.5
  const zAxisEnd = new THREE.Vector3(0, 0, zAxisLength)
  const zAxisArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    zAxisEnd,
    0x44ff88,  // Green for orbital plane normal
    {
      headLength: 0.8 * scale,
      headWidth: 0.35 * scale,
      shaftRadius: 0.06 * scale,
      opacity: 0.7,
    }
  )
  group.add(zAxisArrow)
  
  const zLabel = createTextSprite('L (orbital normal)', { color: '#66ffaa', fontSize: 16 })
  zLabel.position.copy(zAxisEnd.clone().add(new THREE.Vector3(2 * scale, 0.5 * scale, 0)))
  zLabel.scale.set(4, 2, 1)
  group.add(zLabel)
  
  // ========== ORBITAL PLANE SKETCH ==========
  // Small reference circle showing the orbital plane orientation
  const orbitalPlaneRadius = displayDistance * 0.25
  const orbitalPlane = createDashedCircle(orbitalPlaneRadius, 0x88ff88, {
    segments: 32,
    dashSize: 0.5 * scale,
    gapSize: 0.3 * scale,
    opacity: 0.5,
    normal: new THREE.Vector3(0, 0, 1),  // Normal is Z
    center: new THREE.Vector3(0, 0, 0),
  })
  group.add(orbitalPlane)
  
  const planeLabel = createTextSprite('Orbital Plane', { color: '#aaffaa', fontSize: 14 })
  planeLabel.position.set(orbitalPlaneRadius + 2 * scale, 0, 0)
  planeLabel.scale.set(3.5, 1.8, 1)
  group.add(planeLabel)
  
  // ========== V_LSR INDICATOR ==========
  // Show the radial velocity component (approaching/receding)
  const vlsrLength = Math.min(Math.abs(lsrVelocity) / 20, displayDistance * 0.3)
  const vlsrDir = lsrVelocity < 0 ? losDir.clone().negate() : losDir.clone()  // Toward BH if approaching
  const vlsrStart = losEnd.clone().add(losDir.clone().multiplyScalar(-2 * scale))
  const vlsrEndPoint = vlsrStart.clone().add(vlsrDir.multiplyScalar(vlsrLength))
  
  const vlsrColor = lsrVelocity < 0 ? 0xff66ff : 0xffaa66  // Magenta approaching, orange receding
  const vlsrArrow = createArrow(vlsrStart, vlsrEndPoint, vlsrColor, {
    headLength: 0.6 * scale,
    headWidth: 0.25 * scale,
    shaftRadius: 0.05 * scale,
    opacity: 0.8,
  })
  group.add(vlsrArrow)
  
  const vlsrLabel = createTextSprite(`V_LSR = ${lsrVelocity} km/s`, { 
    color: lsrVelocity < 0 ? '#ff88ff' : '#ffcc88', 
    fontSize: 16 
  })
  vlsrLabel.position.copy(vlsrEndPoint.clone().add(new THREE.Vector3(1.5 * scale, 0.5 * scale, 0)))
  vlsrLabel.scale.set(4, 2, 1)
  group.add(vlsrLabel)
  
  const motionLabel = createTextSprite(
    lsrVelocity < 0 ? '(approaching)' : '(receding)',
    { color: lsrVelocity < 0 ? '#88ffff' : '#ffddaa', fontSize: 12 }
  )
  motionLabel.position.copy(vlsrLabel.position.clone().add(new THREE.Vector3(0, -0.8 * scale, 0)))
  motionLabel.scale.set(3, 1.5, 1)
  group.add(motionLabel)
  
  // ========== V_CIRCULAR INDICATOR (Sun's orbital velocity around GC) ==========
  // The LSR moves with V_circular = 220 km/s around the Galactic Center
  // This is tangential to the Sun's orbit, pointing in the direction of galactic rotation (l=90°)
  // 
  // At the observer position, we show:
  // 1. A full circular orbit representing Sun's path around GC
  // 2. An arrow showing V_circular direction (tangent to orbit)
  // 3. Animated marker moving along the full circular orbit
  const vCircular = 220  // km/s - Sun's circular velocity
  const vCircularGroup = new THREE.Group()
  vCircularGroup.name = 'vCircularIndicator'
  
  // Position the V_circular indicator near the observer (Earth)
  vCircularGroup.position.copy(losEnd.clone().add(new THREE.Vector3(0, -4 * scale, 0)))
  
  // Full circular orbit (representing Sun's orbit around GC)
  const orbitRadius = 3 * scale
  const orbitCirclePoints: THREE.Vector3[] = []
  const numSegments = 64
  for (let i = 0; i <= numSegments; i++) {
    const angle = (i / numSegments) * Math.PI * 2
    orbitCirclePoints.push(new THREE.Vector3(
      orbitRadius * Math.cos(angle),
      0,
      orbitRadius * Math.sin(angle)
    ))
  }
  const orbitCircleGeo = new THREE.BufferGeometry().setFromPoints(orbitCirclePoints)
  const orbitCircleMat = new THREE.LineBasicMaterial({
    color: 0x88aaff,
    transparent: true,
    opacity: 0.4,
  })
  const orbitCircle = new THREE.Line(orbitCircleGeo, orbitCircleMat)
  vCircularGroup.add(orbitCircle)
  
  // GC marker at center of orbit
  const gcMarkerGeo = new THREE.SphereGeometry(0.25 * scale, 16, 16)
  const gcMarkerMat = new THREE.MeshBasicMaterial({ color: 0xffaa44 })
  const gcMarker = new THREE.Mesh(gcMarkerGeo, gcMarkerMat)
  gcMarker.position.set(0, 0, 0)
  vCircularGroup.add(gcMarker)
  
  // GC label at center
  const gcCenterLabel = createTextSprite('Galactic Center', { color: '#ffcc66', fontSize: 12 })
  gcCenterLabel.position.set(0, 0.8 * scale, 0)
  gcCenterLabel.scale.set(3.5, 1.8, 1)
  vCircularGroup.add(gcCenterLabel)
  
  // Sun marker on the orbit (fixed position showing current location)
  const sunPosOnOrbit = new THREE.Vector3(orbitRadius, 0, 0)  // Right side (toward GC direction)
  const sunMarkerGeo = new THREE.SphereGeometry(0.35 * scale, 16, 16)
  const sunMarkerMat = new THREE.MeshBasicMaterial({ color: 0xffff44 })
  const sunMarker = new THREE.Mesh(sunMarkerGeo, sunMarkerMat)
  sunMarker.position.copy(sunPosOnOrbit)
  sunMarker.name = 'sunOrbitMarker'
  vCircularGroup.add(sunMarker)
  
  // Sun label
  const sunOrbitLabel = createTextSprite('☉ Sun (8 kpc from GC)', { color: '#ffff88', fontSize: 11 })
  sunOrbitLabel.position.copy(sunPosOnOrbit.clone().add(new THREE.Vector3(0, 0.7 * scale, 0)))
  sunOrbitLabel.scale.set(4, 2, 1)
  vCircularGroup.add(sunOrbitLabel)
  
  // V_circular arrow (tangent to orbit at Sun position, direction of galactic rotation)
  // At Sun's position (+X), tangent direction is +Z (counter-clockwise motion)
  const vCircularDir = new THREE.Vector3(0, 0, 1)  // Tangent direction (l=90°, galactic rotation)
  const vCircularArrowLength = 2.5 * scale
  const vCircularArrowEnd = sunPosOnOrbit.clone().add(vCircularDir.clone().multiplyScalar(vCircularArrowLength))
  
  const vCircularArrow = createArrow(
    sunPosOnOrbit,
    vCircularArrowEnd,
    0x44ff88,  // Green for V_circular
    {
      headLength: 0.5 * scale,
      headWidth: 0.25 * scale,
      shaftRadius: 0.04 * scale,
      opacity: 0.9,
    }
  )
  vCircularGroup.add(vCircularArrow)
  
  // V_circular velocity label
  const vCircularLabel = createTextSprite('V_circular = 220 km/s', { 
    color: '#66ffaa', 
    fontSize: 14 
  })
  vCircularLabel.position.copy(vCircularArrowEnd.clone().add(new THREE.Vector3(0, 0.6 * scale, 0.5 * scale)))
  vCircularLabel.scale.set(4, 2, 1)
  vCircularGroup.add(vCircularLabel)
  
  // Direction label
  const directionLabel = createTextSprite('(Galactic rotation, l = 90°)', { 
    color: '#aaffcc', 
    fontSize: 10 
  })
  directionLabel.position.copy(vCircularLabel.position.clone().add(new THREE.Vector3(0, -0.6 * scale, 0)))
  directionLabel.scale.set(4, 2, 1)
  vCircularGroup.add(directionLabel)
  
  // Animated orbit marker (moves continuously around the full circle)
  const orbitMarkerGeo = new THREE.SphereGeometry(0.2 * scale, 12, 12)
  const orbitMarkerMat = new THREE.MeshBasicMaterial({ 
    color: 0x88ffff,
    transparent: true,
    opacity: 0.9,
  })
  const orbitMarker = new THREE.Mesh(orbitMarkerGeo, orbitMarkerMat)
  orbitMarker.name = 'orbitAnimationMarker'
  orbitMarker.position.set(orbitRadius, 0, 0)  // Start at Sun position
  vCircularGroup.add(orbitMarker)
  
  // Trail effect - add a few trailing markers for visual effect
  for (let i = 1; i <= 3; i++) {
    const trailGeo = new THREE.SphereGeometry((0.2 - i * 0.04) * scale, 8, 8)
    const trailMat = new THREE.MeshBasicMaterial({ 
      color: 0x88ffff,
      transparent: true,
      opacity: 0.6 - i * 0.15,
    })
    const trailMarker = new THREE.Mesh(trailGeo, trailMat)
    trailMarker.name = `orbitTrail${i}`
    // Position slightly behind the main marker
    const trailAngle = -i * 0.15  // Behind by small angle
    trailMarker.position.set(
      orbitRadius * Math.cos(trailAngle),
      0,
      orbitRadius * Math.sin(trailAngle)
    )
    vCircularGroup.add(trailMarker)
  }
  
  // Store animation parameters in userData
  vCircularGroup.userData = {
    orbitRadius,
    animationAngle: 0,  // Current animation angle
    animationSpeed: 0.4,  // rad/s for smooth visualization
  }
  
  // Orbital period annotation at bottom
  const periodLabel = createTextSprite('Orbital Period T ≈ 220 Myr', { 
    color: '#88aacc', 
    fontSize: 11 
  })
  periodLabel.position.set(0, -orbitRadius - 1.0 * scale, 0)
  periodLabel.scale.set(4, 2, 1)
  vCircularGroup.add(periodLabel)
  
  group.add(vCircularGroup)
  
  // ========== POSITION ANGLE REFERENCE ==========
  // Show PA on a small "sky plane" perpendicular to LoS
  const paIndicatorGroup = new THREE.Group()
  
  // Sky plane circle (perpendicular to LoS)
  const skyPlaneRadius = displayDistance * 0.15
  const skyPlane = createDashedCircle(skyPlaneRadius, 0xffff88, {
    segments: 24,
    dashSize: 0.3 * scale,
    gapSize: 0.2 * scale,
    opacity: 0.4,
    normal: losDir,
    center: new THREE.Vector3(0, 0, 0),
  })
  paIndicatorGroup.add(skyPlane)
  
  // PA arrow on sky plane
  // PA is measured from North (celestial) through East
  // Here we show it as rotation from a reference direction in the sky plane
  const northInSky = new THREE.Vector3(0, 0, 1)  // Approximate "up" direction
  const eastInSky = new THREE.Vector3().crossVectors(losDir, northInSky).normalize()
  const northProj = new THREE.Vector3().crossVectors(eastInSky, losDir).normalize()
  
  // PA direction in sky plane
  const paDirection = northProj.clone()
    .multiplyScalar(Math.cos(pa))
    .add(eastInSky.clone().multiplyScalar(Math.sin(pa)))
    .normalize()
  
  const paArrowEnd = paDirection.multiplyScalar(skyPlaneRadius)
  const paArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    paArrowEnd,
    0xffff44,
    {
      headLength: 0.4 * scale,
      headWidth: 0.2 * scale,
      shaftRadius: 0.03 * scale,
      opacity: 0.7,
    }
  )
  paIndicatorGroup.add(paArrow)
  
  // Position the PA indicator near the observer
  paIndicatorGroup.position.copy(losEnd.clone().add(losDir.clone().multiplyScalar(-3 * scale)))
  group.add(paIndicatorGroup)
  
  const paLabel = createTextSprite(`PA = ${config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle}°`, { 
    color: '#ffff88', 
    fontSize: 14 
  })
  paLabel.position.copy(paIndicatorGroup.position.clone().add(new THREE.Vector3(skyPlaneRadius + 2 * scale, 0, 0)))
  paLabel.scale.set(3.5, 1.8, 1)
  group.add(paLabel)
  
  // ========== COORDINATE REFERENCE ==========
  // Brief annotation about the geometry
  const geometryNote = createTextSprite('Observer Geometry (Oka et al. 2017)', { 
    color: '#888888', 
    fontSize: 12 
  })
  geometryNote.position.set(0, -3 * scale, displayDistance * 0.5)
  geometryNote.scale.set(5, 2, 1)
  group.add(geometryNote)
  
  return group
}

/**
 * Create a scale bar showing 1 pc
 */
function createScaleBar(
  position: THREE.Vector3,
  length_pc: number,
  scale: number
): THREE.Group {
  const group = new THREE.Group()
  group.position.copy(position)
  
  const barLength = length_pc * scale
  
  // Horizontal bar
  const barGeo = new THREE.CylinderGeometry(0.03 * scale, 0.03 * scale, barLength, 8)
  const barMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
  const bar = new THREE.Mesh(barGeo, barMat)
  bar.rotation.z = Math.PI / 2
  group.add(bar)
  
  // End caps
  const capSize = 0.15 * scale
  const capGeo = new THREE.CylinderGeometry(0.02 * scale, 0.02 * scale, capSize, 8)
  
  const leftCap = new THREE.Mesh(capGeo, barMat)
  leftCap.position.x = -barLength / 2
  group.add(leftCap)
  
  const rightCap = new THREE.Mesh(capGeo, barMat)
  rightCap.position.x = barLength / 2
  group.add(rightCap)
  
  // Label
  const label = createTextSprite(`${length_pc} pc`, { color: '#ffffff', fontSize: 20 })
  label.position.set(0, 0.4 * scale, 0)
  label.scale.set(3, 1.5, 1)
  group.add(label)
  
  return group
}

export default { 
  createAllGalacticMarkers,
  createObserverDirectionIndicator,
  createSunMarker,
  createEarthMarker,
  createGalacticCenterMarker,
  createHVCCMarker,
  createLineOfSight,
  createLSRVelocity,
  createGalacticAxes,
  createGalacticPlane,
  createGalaxyDisk,
  createSolarSystemDemo,
  updateGalaxyRotation,
  updateVCircularAnimation,
  updateLSRRotation,
  createScaleBar,
  DEFAULT_GALACTIC_CONFIG,
}
