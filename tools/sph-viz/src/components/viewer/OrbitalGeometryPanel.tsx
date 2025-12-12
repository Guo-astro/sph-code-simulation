'use client'

import { useRef, useEffect, useMemo } from 'react'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'
import {
  createTextSprite as createTextSpriteImported,
  createArrow as createArrowImported,
} from './GalacticMarkers'

/**
 * Physical Constants and Galactic Parameters
 * Based on Oka et al. (2017) and standard Galactic parameters
 */
const GALACTIC_CONSTANTS = {
  // Distances
  DISTANCE_SUN_TO_GC_kpc: 8.0,        // Distance from Sun to Galactic Center (kpc)
  DISTANCE_HVCC_TO_GC_pc: 60,         // CO-0.40-0.22 distance from GC (~60 pc projected)
  
  // Velocities
  V_CIRCULAR_SUN_kms: 220,            // Solar circular velocity (km/s)
  V_SUN_PECULIAR_U_kms: 11.1,         // U component (radially inward)
  V_SUN_PECULIAR_V_kms: 12.24,        // V component (direction of rotation)
  V_SUN_PECULIAR_W_kms: 7.25,         // W component (toward NGP)
  
  // HVCC parameters from Oka et al. (2017)
  V_LSR_HVCC_kms: -120,               // V_LSR of HVCC (-120 km/s bulk motion)
  V_DISPERSION_HVCC_kms: 22,          // Velocity dispersion σ_v
  
  // BH orbital velocity around GC (assuming circular at 60 pc)
  // v = sqrt(GM/r), for Sgr A* mass ~4×10^6 M☉
  V_ORBITAL_60PC_kms: 100,            // Approximate orbital velocity at 60 pc
  
  // Physical sizes
  CLOUD_RADIUS_PC: 1.13,              // Cloud radius from N-body model
  HVCC_SIZE_PC: 0.15,                 // HCN clump half-size
}

/**
 * Orbital Geometry Configuration
 * Based on Oka et al. (2017) IMBH-Cloud Scattering Model
 */
export interface OrbitalGeometryConfig {
  // Black hole parameters
  bhMass: number              // M_BH in units of 10^5 M_sun
  
  // Cloud parameters
  cloudRadius: number         // R_cloud in pc
  cloudMass: number           // M_cloud in units of 1000 M_sun
  
  // Orbital parameters (in orbital plane)
  impactParameter: number     // b in pc
  pericentre: number          // r_peri in pc
  eccentricity: number        // e (>1 for hyperbolic)
  
  // Initial conditions
  initialPosition: [number, number, number]  // (X₀, Y₀, Z₀) in pc
  initialVelocity: [number, number, number]  // (vx, vy, vz) in km/s
  
  // Projection angles (for observer)
  inclination: number         // i in degrees (0° = face-on, 90° = edge-on)
  positionAngle: number       // PA in degrees (rotation on sky)
  
  // Observer reference
  lsrVelocity: number         // V_LSR in km/s (bulk motion toward observer)
  distanceToGC: number        // Distance from Earth to Galactic Center in kpc
}

// Default values from Oka et al. (2017) / CAT_OKA
const DEFAULT_CONFIG: OrbitalGeometryConfig = {
  bhMass: 1.0,                // 10^5 M_sun
  cloudRadius: 1.13,          // pc
  cloudMass: 1.0,             // 1000 M_sun
  impactParameter: 5.17,      // pc (updated from config)
  pericentre: 1.69,           // pc
  eccentricity: 1.241,        // hyperbolic
  initialPosition: [20.0, -5.17, 0],  // pc
  initialVelocity: [-9.35, 4.48, 0],  // km/s
  inclination: 70,            // degrees
  positionAngle: 41.6,        // degrees
  lsrVelocity: -120,          // km/s
  distanceToGC: 8.0,          // kpc
}

export interface OrbitalGeometryPanelProps {
  config?: Partial<OrbitalGeometryConfig>
  className?: string
  width?: number
  height?: number
  showLabels?: boolean
  viewMode?: 'orbital-plane' | 'observer' | '3d' | 'galactic'
  // Synchronization with simulation
  currentTime?: number        // Current simulation time (code units)
  currentCloudPosition?: [number, number, number]  // Current cloud centroid (pc)
  currentCloudVelocity?: [number, number, number]  // Current cloud velocity (km/s)
}

/**
 * Interactive panel showing orbital geometry for IMBH-cloud encounter
 * 
 * Features:
 * - Sun-centered galactic coordinate system
 * - Impact parameter visualization
 * - Pericentre marker
 * - Orbital plane with hyperbolic trajectory
 * - Line of sight direction (Earth to HVCC)
 * - Inclination and position angle indicators
 * - Earth/Observer reference frame
 * - LSR velocity vector
 * - Velocity arrows with proper scales
 * - Synchronization with simulation time
 */
export function OrbitalGeometryPanel({
  config: configOverride,
  className = '',
  width = 600,
  height = 500,
  showLabels = true,
  viewMode = '3d',
  currentTime = 0,
  currentCloudPosition,
  currentCloudVelocity,
}: OrbitalGeometryPanelProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null)
  const sceneRef = useRef<THREE.Scene | null>(null)
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null)
  const cloudMarkerRef = useRef<THREE.Group | null>(null)
  const velocityArrowRef = useRef<THREE.Group | null>(null)
  const timeRef = useRef<{ label: THREE.Sprite | null, value: THREE.Sprite | null }>({ label: null, value: null })
  
  // Merge with defaults
  const config = useMemo(() => ({
    ...DEFAULT_CONFIG,
    ...configOverride,
  }), [configOverride])

  useEffect(() => {
    if (!containerRef.current) return

    const container = containerRef.current
    
    // Scene setup
    const scene = new THREE.Scene()
    scene.background = new THREE.Color(0x0a0a1a)
    sceneRef.current = scene
    
    // Camera
    const camera = new THREE.PerspectiveCamera(50, width / height, 0.1, 2000)
    cameraRef.current = camera
    
    // Position camera based on view mode
    // In galactic view, we scale down: 1 unit = 100 pc for overview, or 1 unit = 1 pc for local
    const isGalacticView = viewMode === 'galactic'
    const sceneScale = isGalacticView ? 0.01 : 1.0  // 1 kpc = 10 units in galactic view
    
    switch (viewMode) {
      case 'orbital-plane':
        camera.position.set(0, 0, 50)
        break
      case 'observer':
        const inc = config.inclination * Math.PI / 180
        const pa = config.positionAngle * Math.PI / 180
        camera.position.set(
          50 * Math.sin(inc) * Math.cos(pa),
          50 * Math.sin(inc) * Math.sin(pa),
          50 * Math.cos(inc)
        )
        break
      case 'galactic':
        // View of entire Milky Way scale
        camera.position.set(100, 60, 100)
        break
      default:
        camera.position.set(35, 25, 40)
    }
    camera.lookAt(0, 0, 0)
    
    // Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true })
    renderer.setSize(width, height)
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    container.appendChild(renderer.domElement)
    rendererRef.current = renderer
    
    // Controls
    const controls = new OrbitControls(camera, renderer.domElement)
    controls.enableDamping = true
    controls.dampingFactor = 0.05
    controls.target.set(0, 0, 0)
    
    // ========== HELPER FUNCTIONS ==========
    // Wrappers around imported functions to match local usage patterns
    
    const createTextSprite = (
      text: string, 
      color: string = '#ffffff',
      fontSize: number = 32,
      bgColor?: string
    ): THREE.Sprite => {
      return createTextSpriteImported(text, { color, fontSize, backgroundColor: bgColor })
    }
    
    const createArrow = (
      from: THREE.Vector3,
      to: THREE.Vector3,
      color: number,
      headLength: number = 0.5,
      headWidth: number = 0.3
    ): THREE.Group => {
      return createArrowImported(from, to, color, { headLength, headWidth })
    }
    
    // Local helper for simple circles (not available in GalacticMarkers)
    const createCircle = (
      radius: number,
      color: number,
      dashed: boolean = false,
      segments: number = 64
    ): THREE.Line => {
      const points: THREE.Vector3[] = []
      for (let i = 0; i <= segments; i++) {
        const theta = (i / segments) * Math.PI * 2
        points.push(new THREE.Vector3(radius * Math.cos(theta), radius * Math.sin(theta), 0))
      }
      const geo = new THREE.BufferGeometry().setFromPoints(points)
      const mat = dashed 
        ? new THREE.LineDashedMaterial({ color, dashSize: 0.3, gapSize: 0.15 })
        : new THREE.LineBasicMaterial({ color })
      const line = new THREE.Line(geo, mat)
      if (dashed) line.computeLineDistances()
      return line
    }
    
    // ========== COORDINATE SYSTEM ==========
    // 
    // We use a BH-centered coordinate system for the orbital dynamics:
    //   - Origin: IMBH (black hole)
    //   - X: Toward pericentre (in orbital plane)
    //   - Y: Perpendicular to X in orbital plane (impact parameter direction)
    //   - Z: Orbital angular momentum direction (perpendicular to orbital plane)
    //
    // OBSERVER GEOMETRY (Oka et al. 2017):
    //   - Inclination (i): Angle between line-of-sight and orbital angular momentum
    //     * i = 0°  → face-on (observer looks along +Z, sees full orbital plane)
    //     * i = 90° → edge-on (observer is in XY plane)
    //   - Position Angle (PA): Rotation of the ascending node on the sky
    //     * Measured from celestial/galactic North through East
    //     * PA = 0° → ascending node points North
    //
    // For visualization, the camera represents the OBSERVER (Sun/Earth):
    //   - Camera looks AT the BH-cloud system (origin)
    //   - Camera position is determined by inclination and PA
    
    const incRad = config.inclination * Math.PI / 180  // i in radians
    const paRad = config.positionAngle * Math.PI / 180 // PA in radians
    
    // LINE OF SIGHT: Direction FROM observer TO system
    // Using standard astronomical spherical coordinates:
    //   - i = 0° means looking along +Z (face-on to XY orbital plane)
    //   - i = 90° means looking in XY plane
    //   - PA rotates the viewing direction around Z
    //
    // Observer position (camera) is OPPOSITE to line-of-sight direction:
    //   - Camera at position looking toward origin
    const observerDist = 50  // Display distance for camera/observer marker
    const observerPosition = new THREE.Vector3(
      observerDist * Math.sin(incRad) * Math.sin(paRad),   // X component
      -observerDist * Math.sin(incRad) * Math.cos(paRad),  // Y component (PA from North = -Y)
      observerDist * Math.cos(incRad)                       // Z component
    )
    
    // Set camera to observer position if in observer view mode
    if (viewMode === 'observer') {
      camera.position.copy(observerPosition)
      camera.lookAt(0, 0, 0)
    }
    
    // ========== SUN/OBSERVER MARKER ==========
    // Show where the observer (Sun/Earth) is relative to the system
    
    const sunGroup = new THREE.Group()
    
    // Corona glow
    const coronaGeo = new THREE.SphereGeometry(1.0, 32, 32)
    const coronaMat = new THREE.MeshBasicMaterial({
      color: 0xffdd44,
      transparent: true,
      opacity: 0.25,
    })
    sunGroup.add(new THREE.Mesh(coronaGeo, coronaMat))
    
    // Main sun body
    const sunGeo = new THREE.SphereGeometry(0.7, 32, 32)
    const sunMat = new THREE.MeshBasicMaterial({ color: 0xffff00 })
    sunGroup.add(new THREE.Mesh(sunGeo, sunMat))
    
    // Core
    const sunCoreGeo = new THREE.SphereGeometry(0.25, 32, 32)
    const sunCoreMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
    sunGroup.add(new THREE.Mesh(sunCoreGeo, sunCoreMat))
    
    // Position Sun marker at scaled observer position (actual distance is 8 kpc)
    const sunDisplayDist = 25
    const sunPosition = observerPosition.clone().normalize().multiplyScalar(sunDisplayDist)
    sunGroup.position.copy(sunPosition)
    scene.add(sunGroup)
    
    if (showLabels) {
      const sunLabel = createTextSprite('☉ Sun', '#ffff00', 28)
      sunLabel.position.copy(sunGroup.position.clone().add(new THREE.Vector3(0, 2, 0)))
      scene.add(sunLabel)
    }
    
    // ========== EARTH MARKER (very close to Sun) ==========
    
    const earthGroup = new THREE.Group()
    
    // Atmosphere
    const atmoGeo = new THREE.SphereGeometry(0.5, 32, 32)
    const atmoMat = new THREE.MeshBasicMaterial({
      color: 0x88ccff,
      transparent: true,
      opacity: 0.2,
    })
    earthGroup.add(new THREE.Mesh(atmoGeo, atmoMat))
    
    // Ocean
    const oceanGeo = new THREE.SphereGeometry(0.35, 32, 32)
    const oceanMat = new THREE.MeshBasicMaterial({ color: 0x2266bb })
    earthGroup.add(new THREE.Mesh(oceanGeo, oceanMat))
    
    // Land patches
    const landGeo = new THREE.SphereGeometry(0.37, 16, 16, 0, Math.PI * 0.7, Math.PI * 0.25, Math.PI * 0.35)
    const landMat = new THREE.MeshBasicMaterial({ color: 0x228833 })
    const land = new THREE.Mesh(landGeo, landMat)
    land.rotation.y = 0.5
    earthGroup.add(land)
    
    // Position Earth slightly offset from Sun (in tangent direction)
    const tangentOffset = new THREE.Vector3(1.5, 0.5, 0)
    earthGroup.position.copy(sunGroup.position.clone().add(tangentOffset))
    scene.add(earthGroup)
    
    if (showLabels) {
      const earthLabel = createTextSprite('⊕ Earth (Observer)', '#88ccff', 24)
      earthLabel.position.copy(earthGroup.position.clone().add(new THREE.Vector3(0, 1.2, 0)))
      scene.add(earthLabel)
      
      // Distance label  
      const distLabel = createTextSprite(`d ≈ ${config.distanceToGC} kpc from GC`, '#aaaaff', 18)
      distLabel.position.copy(earthGroup.position.clone().add(new THREE.Vector3(0, -1.0, 0)))
      scene.add(distLabel)
    }
    
    // ========== LINE OF SIGHT ARROW ==========
    // Arrow from observer toward the BH-cloud system
    
    const losArrow = createArrow(
      sunGroup.position.clone(),
      sunGroup.position.clone().normalize().multiplyScalar(5),  // Point toward origin
      0xffff00,
      1.0,
      0.5
    )
    scene.add(losArrow)
    
    if (showLabels) {
      const losLabel = createTextSprite('Line of Sight', '#ffff88', 18)
      losLabel.position.copy(sunGroup.position.clone().normalize().multiplyScalar(15))
      scene.add(losLabel)
    }
    
    // ========== COORDINATE AXES (BH-centered orbital frame) ==========
    
    const axesGroup = new THREE.Group()
    
    // X-axis (toward pericentre in orbital plane)
    const xAxis = createArrow(
      new THREE.Vector3(-15, 0, 0),
      new THREE.Vector3(15, 0, 0),
      0xff4444
    )
    axesGroup.add(xAxis)
    
    // Y-axis (perpendicular in orbital plane, impact parameter direction)
    const yAxis = createArrow(
      new THREE.Vector3(0, -10, 0),
      new THREE.Vector3(0, 15, 0),
      0x44ff44
    )
    axesGroup.add(yAxis)
    
    // Z-axis (orbital angular momentum direction, perpendicular to orbital plane)
    const zAxis = createArrow(
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, 0, 12),
      0x4488ff
    )
    axesGroup.add(zAxis)
    
    if (showLabels) {
      const xLabel = createTextSprite('X → peri', '#ff6666')
      xLabel.position.set(16, 0, 0)
      axesGroup.add(xLabel)
      
      const yLabel = createTextSprite('Y (b)', '#66ff66')
      yLabel.position.set(0, 16, 0)
      axesGroup.add(yLabel)
      
      const zLabel = createTextSprite('Z (L)', '#6688ff')
      zLabel.position.set(0, 0, 13)
      axesGroup.add(zLabel)
    }
    
    scene.add(axesGroup)
    
    // ========== ORBITAL PLANE ==========
    
    const planeGeo = new THREE.PlaneGeometry(50, 50)
    const planeMat = new THREE.MeshBasicMaterial({
      color: 0x3366ff,
      transparent: true,
      opacity: 0.08,
      side: THREE.DoubleSide,
    })
    const orbitalPlane = new THREE.Mesh(planeGeo, planeMat)
    scene.add(orbitalPlane)
    
    // Plane border
    const planeBorder = createCircle(20, 0x3366ff, true)
    scene.add(planeBorder)
    
    if (showLabels) {
      const planeLabel = createTextSprite('Orbital Plane', '#6699ff', 22)
      planeLabel.position.set(12, -12, 0.1)
      scene.add(planeLabel)
    }
    
    // ========== BLACK HOLE (IMBH at origin) ==========
    
    const bhGroup = new THREE.Group()
    
    // Glow
    const glowGeo = new THREE.SphereGeometry(0.8, 32, 32)
    const glowMat = new THREE.MeshBasicMaterial({
      color: 0xff4444,
      transparent: true,
      opacity: 0.3,
    })
    bhGroup.add(new THREE.Mesh(glowGeo, glowMat))
    
    // Core
    const coreGeo = new THREE.SphereGeometry(0.3, 32, 32)
    const coreMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
    bhGroup.add(new THREE.Mesh(coreGeo, coreMat))
    
    // Event horizon
    const horizonGeo = new THREE.SphereGeometry(0.15, 32, 32)
    const horizonMat = new THREE.MeshBasicMaterial({ color: 0x000000 })
    bhGroup.add(new THREE.Mesh(horizonGeo, horizonMat))
    
    scene.add(bhGroup)
    
    if (showLabels) {
      const bhLabel = createTextSprite('IMBH', '#ff8888', 28)
      bhLabel.position.set(0, -1.5, 0)
      scene.add(bhLabel)
      
      const bhMassLabel = createTextSprite(`${(config.bhMass * 100).toFixed(0)}×10³ M☉`, '#ffaaaa', 20)
      bhMassLabel.position.set(0, -2.5, 0)
      scene.add(bhMassLabel)
    }
    
    // ========== DIRECTION TO GALACTIC CENTER ==========
    
    // From IMBH, GC is about 60 pc away
    // Direction: approximately -X (toward lower Galactic longitude)
    const gcDir = new THREE.Vector3(-1, 0, 0).normalize()
    const gcArrow = createArrow(
      new THREE.Vector3(0, 0, 0),
      gcDir.clone().multiplyScalar(8),
      0xff8800,
      0.6,
      0.35
    )
    scene.add(gcArrow)
    
    if (showLabels) {
      const gcLabel = createTextSprite('→ Sgr A* (60 pc)', '#ff9944', 20)
      gcLabel.position.copy(gcDir.clone().multiplyScalar(10))
      scene.add(gcLabel)
    }
    
    // ========== GALACTIC ROTATION DIRECTION ==========
    
    // At the Sun's position, rotation is toward +Y (l = 90°)
    // Near GC, rotation direction varies - approximate as +Y
    const rotArrow = createArrow(
      new THREE.Vector3(-8, -8, 0),
      new THREE.Vector3(-8, -8 + 5, 0),
      0x44ffaa,
      0.4,
      0.25
    )
    scene.add(rotArrow)
    
    if (showLabels) {
      const rotLabel = createTextSprite('Galactic rotation', '#66ffcc', 18)
      rotLabel.position.set(-8, -5.5, 0)
      scene.add(rotLabel)
      
      const rotVal = createTextSprite(`~${GALACTIC_CONSTANTS.V_CIRCULAR_SUN_kms} km/s at Sun`, '#88ffcc', 14)
      rotVal.position.set(-8, -6.8, 0)
      scene.add(rotVal)
    }
    
    // ========== HYPERBOLIC ORBIT ==========
    //
    // CORRECT ORBITAL MECHANICS - Computed from initial state vectors:
    // The orbit is computed using the actual initial position and velocity,
    // NOT the preset parameters which may be inconsistent.
    //
    // Formulas:
    //   ε = v²/2 - GM/r  (specific orbital energy)
    //   h = |r × v|      (specific angular momentum)
    //   p = h² / (GM)    (semi-latus rectum)
    //   e = √(1 + 2εh²/(GM)²)  (eccentricity)
    //   a = GM / (2ε)    (semi-major axis)
    //   rp = a(e - 1)    (pericentre)

    const [x0, y0, z0] = config.initialPosition
    const [vx0, vy0, vz0] = config.initialVelocity
    const GM = config.bhMass  // G=1 in code units

    // Compute from state vectors
    const r0 = Math.sqrt(x0 * x0 + y0 * y0)
    const v0 = Math.sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0)

    // Specific orbital energy
    const epsilon = v0 * v0 / 2 - GM / r0

    // Specific angular momentum (z-component for 2D orbit)
    const h = Math.abs(x0 * vy0 - y0 * vx0)

    // Semi-latus rectum
    const p = h * h / GM

    // Eccentricity
    const e = Math.sqrt(1 + 2 * epsilon * h * h / (GM * GM))

    // Semi-major axis and pericentre (for display)
    const a = GM / (2 * epsilon)
    const rp = a * (e - 1)

    // ========== Orbit Orientation ==========
    const theta0 = Math.atan2(y0, x0)
    const cosNu0 = Math.max(-0.9999, Math.min(0.9999, (p / r0 - 1) / e))

    // Determine sign of ν from radial velocity (r·v < 0 means approaching)
    const rdotv = x0 * vx0 + y0 * vy0
    const nu0 = rdotv < 0 ? -Math.acos(cosNu0) : Math.acos(cosNu0)

    // Argument of pericentre: ω = θ0 - ν0
    const omega = theta0 - nu0

    // Asymptotic true anomaly for hyperbola: ν_asymp = acos(-1/e)
    const nuAsymp = Math.acos(-1 / e)
    const nuMax = nuAsymp - 0.05
    const nuMin = -nuMax

    const orbitPoints: THREE.Vector3[] = []
    const numPoints = 100
    for (let i = 0; i < numPoints; i++) {
      const nu = nuMin + (nuMax - nuMin) * (i / (numPoints - 1))
      const denom = 1 + e * Math.cos(nu)
      if (denom <= 0) continue
      const r = p / denom

      // Convert to Cartesian: θ = ω + ν
      const theta = omega + nu
      const x = r * Math.cos(theta)
      const y = r * Math.sin(theta)

      if (Math.abs(x) < 30 && Math.abs(y) < 30) {
        orbitPoints.push(new THREE.Vector3(x, y, 0))
      }
    }
    
    const orbitGeo = new THREE.BufferGeometry().setFromPoints(orbitPoints)
    const orbitMat = new THREE.LineBasicMaterial({ color: 0x00ff88, linewidth: 2 })
    const orbitLine = new THREE.Line(orbitGeo, orbitMat)
    scene.add(orbitLine)
    
    // ========== PERICENTRE MARKER ==========
    // Pericentre is at distance rp from BH, at angle omega

    const periX = rp * Math.cos(omega)
    const periY = rp * Math.sin(omega)

    const periMarker = new THREE.Group()
    const periGeo = new THREE.SphereGeometry(0.2, 16, 16)
    const periMat = new THREE.MeshBasicMaterial({ color: 0xffff00 })
    periMarker.add(new THREE.Mesh(periGeo, periMat))

    // Line from origin to pericentre
    const periLineGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(periX, periY, 0)
    ])
    const periLineMat = new THREE.LineDashedMaterial({
      color: 0xffff00,
      dashSize: 0.2,
      gapSize: 0.1
    })
    const periLine = new THREE.Line(periLineGeo, periLineMat)
    periLine.computeLineDistances()
    scene.add(periLine)
    periMarker.position.set(periX, periY, 0)
    scene.add(periMarker)

    if (showLabels) {
      const periLabel = createTextSprite('Pericentre', '#ffff66', 20)
      periLabel.position.set(periX + 0.5, periY + 1.2, 0)
      scene.add(periLabel)

      const periValue = createTextSprite(`r_p = ${rp.toFixed(2)} pc`, '#ffff99', 16)
      periValue.position.set(periX + 0.5, periY + 0.3, 0)
      scene.add(periValue)
    }

    // ========== IMPACT PARAMETER ==========
    // Impact parameter b is perpendicular distance from BH to incoming asymptote
    // The asymptote direction is determined by the orbit geometry
    // Incoming asymptote angle: omega + (π - ν_asymp) = omega + π - acos(-1/e)

    const b = config.impactParameter
    const asymptoteAngle = omega + Math.PI - nuAsymp  // Direction of incoming asymptote
    const asymptoteDir = new THREE.Vector3(Math.cos(asymptoteAngle), Math.sin(asymptoteAngle), 0)

    // Perpendicular to asymptote (for impact parameter line)
    const perpDir = new THREE.Vector3(-Math.sin(asymptoteAngle), Math.cos(asymptoteAngle), 0)

    // Draw asymptote line (extended far in both directions)
    // Point on asymptote closest to origin is at distance b along perpendicular
    const asymptoteBase = perpDir.clone().multiplyScalar(b)
    const asymptoteStart = asymptoteBase.clone().add(asymptoteDir.clone().multiplyScalar(-25))
    const asymptoteEnd = asymptoteBase.clone().add(asymptoteDir.clone().multiplyScalar(15))

    const asymptoteGeo = new THREE.BufferGeometry().setFromPoints([asymptoteStart, asymptoteEnd])
    const asymptoteMat = new THREE.LineDashedMaterial({
      color: 0xff8800,
      dashSize: 0.5,
      gapSize: 0.25,
    })
    const asymptoteLine = new THREE.Line(asymptoteGeo, asymptoteMat)
    asymptoteLine.computeLineDistances()
    scene.add(asymptoteLine)

    // Impact parameter line (from origin perpendicular to asymptote)
    const impactLineGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      asymptoteBase.clone()
    ])
    const impactLineMat = new THREE.LineBasicMaterial({ color: 0xff8800 })
    const impactLine = new THREE.Line(impactLineGeo, impactLineMat)
    scene.add(impactLine)

    if (showLabels) {
      const bLabel = createTextSprite('b (impact)', '#ffaa44', 18)
      const labelPos = perpDir.clone().multiplyScalar(b * 0.5)
      bLabel.position.set(labelPos.x - 1.5, labelPos.y, 0)
      scene.add(bLabel)

      const bValue = createTextSprite(`= ${b.toFixed(2)} pc`, '#ffcc66', 16)
      bValue.position.set(labelPos.x - 1.5, labelPos.y - 1, 0)
      scene.add(bValue)
    }
    
    // ========== CLOUD (INITIAL POSITION - will be updated by second useEffect) ==========

    // Start with initial position from config - the second useEffect will update it in real-time
    // x0, y0, z0, vx0, vy0, vz0 already extracted above for orbit calculation

    const cloudGroup = new THREE.Group()
    const cloudGeo = new THREE.SphereGeometry(config.cloudRadius, 32, 32)
    const cloudMat = new THREE.MeshBasicMaterial({
      color: 0x4488ff,
      transparent: true,
      opacity: 0.4,
    })
    cloudGroup.add(new THREE.Mesh(cloudGeo, cloudMat))
    
    const cloudOutline = createCircle(config.cloudRadius, 0x66aaff)
    cloudGroup.add(cloudOutline)
    cloudGroup.position.set(x0, y0, z0)
    scene.add(cloudGroup)
    cloudMarkerRef.current = cloudGroup
    
    // Cloud velocity vector
    const velMag = Math.sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0)
    const velScale = 0.5  // Scale for visualization
    const velArrow = createArrow(
      new THREE.Vector3(x0, y0, z0),
      new THREE.Vector3(x0 + vx0 * velScale, y0 + vy0 * velScale, z0 + vz0 * velScale),
      0x44aaff,
      0.4,
      0.25
    )
    scene.add(velArrow)
    velocityArrowRef.current = velArrow
    
    if (showLabels) {
      const cloudLabel = createTextSprite(currentCloudPosition ? 'Cloud (t)' : 'Cloud (t=0)', '#66aaff', 20)
      cloudLabel.position.set(x0, y0 + config.cloudRadius + 1.5, z0)
      scene.add(cloudLabel)
      
      const posLabel = createTextSprite(
        `(${x0.toFixed(1)}, ${y0.toFixed(1)}, ${z0.toFixed(1)}) pc`,
        '#88ccff', 14
      )
      posLabel.position.set(x0, y0 - config.cloudRadius - 1.5, z0)
      scene.add(posLabel)
      
      const velLabel = createTextSprite(`v = ${velMag.toFixed(1)} km/s`, '#44aaff', 16)
      velLabel.position.set(
        x0 + vx0 * velScale + 1.5,
        y0 + vy0 * velScale + 0.5,
        z0
      )
      scene.add(velLabel)
    }
    
    // ========== V_LSR VECTOR ==========
    
    // Line of sight direction (from origin toward observer)
    const losDir = observerPosition.clone().normalize()
    
    // V_LSR is the cloud's velocity relative to Local Standard of Rest
    // Negative = approaching Earth (along line of sight toward observer)
    const lsrArrowScale = Math.abs(config.lsrVelocity) / 30  // Scale for visualization
    const lsrStart = new THREE.Vector3(-10, -10, 0)
    // Negative V_LSR means approaching, so arrow points toward observer
    const lsrEnd = lsrStart.clone().add(losDir.clone().multiplyScalar(lsrArrowScale * (config.lsrVelocity < 0 ? 1 : -1)))
    const lsrArrow = createArrow(lsrStart, lsrEnd, 0xff6688, 0.3, 0.2)
    scene.add(lsrArrow)
    
    if (showLabels) {
      const lsrLabel = createTextSprite('V_LSR', '#ff8899', 18)
      lsrLabel.position.set(-10, -11.5, 0)
      scene.add(lsrLabel)
      
      const lsrValue = createTextSprite(`= ${config.lsrVelocity} km/s`, '#ffaabb', 14)
      lsrValue.position.set(-10, -12.8, 0)
      scene.add(lsrValue)
    }
    
    // ========== INCLINATION ARC ==========
    // Arc showing the inclination angle from Z-axis (face-on) to observer direction
    // i = 0° means face-on (observer along +Z), i = 90° means edge-on
    
    const incArcPoints: THREE.Vector3[] = []
    const arcRadius = 8
    for (let i = 0; i <= 20; i++) {
      const angle = (i / 20) * incRad
      // Arc from +Z axis toward the observer direction
      incArcPoints.push(new THREE.Vector3(
        arcRadius * Math.sin(angle) * Math.sin(paRad),
        -arcRadius * Math.sin(angle) * Math.cos(paRad),
        arcRadius * Math.cos(angle)
      ))
    }
    const incArcGeo = new THREE.BufferGeometry().setFromPoints(incArcPoints)
    const incArcMat = new THREE.LineBasicMaterial({ color: 0xffaa00 })
    const incArc = new THREE.Line(incArcGeo, incArcMat)
    scene.add(incArc)
    
    if (showLabels) {
      const incLabel = createTextSprite(`i = ${config.inclination}°`, '#ffcc44', 20)
      const midAngle = incRad / 2
      incLabel.position.set(
        (arcRadius + 2) * Math.sin(midAngle) * Math.sin(paRad),
        -(arcRadius + 2) * Math.sin(midAngle) * Math.cos(paRad),
        (arcRadius + 2) * Math.cos(midAngle)
      )
      scene.add(incLabel)
    }
    
    // ========== POSITION ANGLE ARC ==========
    // Arc in the plane perpendicular to observer direction, showing PA from North
    // PA measured from -Y (galactic North proxy) through +X (East)
    
    const paArcPoints: THREE.Vector3[] = []
    const paRadius = 6
    for (let i = 0; i <= 20; i++) {
      const angle = (i / 20) * paRad
      // Arc at a small height above XY plane to show PA rotation
      paArcPoints.push(new THREE.Vector3(
        paRadius * Math.sin(angle),
        -paRadius * Math.cos(angle),
        2  // Slightly above XY plane for visibility
      ))
    }
    const paArcGeo = new THREE.BufferGeometry().setFromPoints(paArcPoints)
    const paArcMat = new THREE.LineDashedMaterial({ color: 0x44ffaa, dashSize: 0.2, gapSize: 0.1 })
    const paArc = new THREE.Line(paArcGeo, paArcMat)
    paArc.computeLineDistances()
    scene.add(paArc)
    
    if (showLabels) {
      const paLabel = createTextSprite(`PA = ${config.positionAngle}°`, '#66ffcc', 18)
      paLabel.position.set(paRadius + 2, 0, 3)
      scene.add(paLabel)
    }
    
    // ========== SCALE REFERENCE ==========
    
    // Add a scale bar to verify distances
    const scaleBarStart = new THREE.Vector3(-20, -15, 0)
    const scaleBarEnd = new THREE.Vector3(-10, -15, 0)
    const scaleBarGeo = new THREE.BufferGeometry().setFromPoints([scaleBarStart, scaleBarEnd])
    const scaleBarMat = new THREE.LineBasicMaterial({ color: 0xffffff })
    const scaleBar = new THREE.Line(scaleBarGeo, scaleBarMat)
    scene.add(scaleBar)
    
    // Scale bar end caps
    const capGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(-20, -15.3, 0),
      new THREE.Vector3(-20, -14.7, 0)
    ])
    scene.add(new THREE.Line(capGeo, scaleBarMat))
    const capGeo2 = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(-10, -15.3, 0),
      new THREE.Vector3(-10, -14.7, 0)
    ])
    scene.add(new THREE.Line(capGeo2, scaleBarMat))
    
    if (showLabels) {
      const scaleLabel = createTextSprite('10 pc', '#ffffff', 20)
      scaleLabel.position.set(-15, -16.5, 0)
      scene.add(scaleLabel)
    }
    
    // ========== LEGEND ==========
    
    const legendGroup = new THREE.Group()
    legendGroup.position.set(-24, 12, 0)
    
    const legendBg = new THREE.Mesh(
      new THREE.PlaneGeometry(13, 16),
      new THREE.MeshBasicMaterial({ color: 0x0a0a1a, transparent: true, opacity: 0.9 })
    )
    legendGroup.add(legendBg)
    
    const legendTitle = createTextSprite('LEGEND', '#ffffff', 26)
    legendTitle.position.set(0, 6, 0.1)
    legendTitle.scale.set(4, 2, 1)
    legendGroup.add(legendTitle)
    
    const items = [
      { color: '#ffff00', text: '☉ Sun (Observer)' },
      { color: '#88ccff', text: '⊕ Earth/ALMA' },
      { color: '#00ff88', text: 'Hyperbolic orbit' },
      { color: '#ff8800', text: 'Impact param b' },
      { color: '#ffff00', text: 'Pericentre r_p' },
      { color: '#ff44ff', text: 'Line of Sight' },
      { color: '#4488ff', text: 'Cloud position' },
      { color: '#ff8800', text: '→ Galactic Center' },
    ]
    
    items.forEach((item, i) => {
      const dot = new THREE.Mesh(
        new THREE.SphereGeometry(0.2, 8, 8),
        new THREE.MeshBasicMaterial({ color: item.color })
      )
      dot.position.set(-5, 4.5 - i * 1.2, 0.1)
      legendGroup.add(dot)
      
      const label = createTextSprite(item.text, item.color, 16)
      label.position.set(1.5, 4.5 - i * 1.2, 0.1)
      label.scale.set(4, 2, 1)
      legendGroup.add(label)
    })
    
    scene.add(legendGroup)
    
    // ========== SIMULATION TIME DISPLAY ==========
    
    const timeGroup = new THREE.Group()
    timeGroup.position.set(20, 14, 0)
    
    const timeBg = new THREE.Mesh(
      new THREE.PlaneGeometry(10, 3),
      new THREE.MeshBasicMaterial({ color: 0x0a0a1a, transparent: true, opacity: 0.9 })
    )
    timeGroup.add(timeBg)
    
    const timeLabel = createTextSprite('Simulation Time', '#88ffff', 18)
    timeLabel.position.set(0, 0.8, 0.1)
    timeLabel.scale.set(4, 2, 1)
    timeGroup.add(timeLabel)
    timeRef.current.label = timeLabel
    
    const timeMyr = (currentTime || 0) * 0.978
    const timeValue = createTextSprite(`t = ${timeMyr.toFixed(3)} Myr`, '#aaffff', 22)
    timeValue.position.set(0, -0.3, 0.1)
    timeValue.scale.set(4, 2, 1)
    timeGroup.add(timeValue)
    timeRef.current.value = timeValue
    
    scene.add(timeGroup)
    
    // ========== RENDER LOOP ==========
    
    let animationId: number
    const animate = () => {
      animationId = requestAnimationFrame(animate)
      controls.update()
      renderer.render(scene, camera)
    }
    animate()
    
    // Cleanup
    return () => {
      cancelAnimationFrame(animationId)
      renderer.dispose()
      container.removeChild(renderer.domElement)
    }
  }, [config, width, height, showLabels, viewMode])  // Removed currentTime, currentCloudPosition, currentCloudVelocity to prevent full scene rebuild on every frame

  // Update cloud position when simulation time changes
  useEffect(() => {
    if (!sceneRef.current || !cloudMarkerRef.current) {
      return
    }
    
    // Update cloud marker position if we have current position from simulation
    if (currentCloudPosition) {
      cloudMarkerRef.current.position.set(
        currentCloudPosition[0],
        currentCloudPosition[1],
        currentCloudPosition[2]
      )
    }
    
    // Update velocity arrow
    if (velocityArrowRef.current && currentCloudVelocity && currentCloudPosition) {
      // Remove old velocity arrow
      sceneRef.current.remove(velocityArrowRef.current)
      
      // Calculate velocity magnitude and direction
      const velMag = Math.sqrt(
        currentCloudVelocity[0] ** 2 +
        currentCloudVelocity[1] ** 2 +
        currentCloudVelocity[2] ** 2
      )
      
      // Create new arrow (scale velocity for visibility)
      const velScale = 0.5  // Scale factor for visualization
      const velArrow = createArrowHelper(
        new THREE.Vector3(currentCloudPosition[0], currentCloudPosition[1], currentCloudPosition[2]),
        new THREE.Vector3(
          currentCloudPosition[0] + currentCloudVelocity[0] * velScale,
          currentCloudPosition[1] + currentCloudVelocity[1] * velScale,
          currentCloudPosition[2] + currentCloudVelocity[2] * velScale
        ),
        0x44aaff,
        0.4,
        0.25
      )
      velocityArrowRef.current = velArrow
      sceneRef.current.add(velArrow)
    }
    
    // Update time display
    if (timeRef.current.value && currentTime !== undefined) {
      // Convert code time to Myr (1 code time ≈ 0.978 Myr)
      const timeMyr = currentTime * 0.978
      // Update canvas texture for time value
      const sprite = timeRef.current.value
      if (sprite.material.map) {
        const canvas = document.createElement('canvas')
        canvas.width = 256
        canvas.height = 64
        const ctx = canvas.getContext('2d')!
        ctx.fillStyle = '#88ffff'
        ctx.font = 'Bold 28px Arial'
        ctx.textAlign = 'left'
        ctx.fillText(`t = ${timeMyr.toFixed(2)} Myr`, 10, 40)
        sprite.material.map.dispose()
        sprite.material.map = new THREE.CanvasTexture(canvas)
        sprite.material.map.needsUpdate = true
      }
    }
  }, [currentTime, currentCloudPosition, currentCloudVelocity])

  // Helper function for arrow creation (declared outside useEffect for reuse)
  function createArrowHelper(
    from: THREE.Vector3,
    to: THREE.Vector3,
    color: number,
    headLength: number = 0.5,
    headWidth: number = 0.3
  ): THREE.Group {
    const group = new THREE.Group()
    const direction = to.clone().sub(from).normalize()
    const length = from.distanceTo(to)
    
    if (length <= headLength) return group
    
    // Shaft
    const shaftGeo = new THREE.CylinderGeometry(0.05, 0.05, length - headLength, 8)
    const shaftMat = new THREE.MeshBasicMaterial({ color })
    const shaft = new THREE.Mesh(shaftGeo, shaftMat)
    shaft.position.copy(from.clone().lerp(to, (length - headLength) / 2 / length))
    shaft.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
    group.add(shaft)
    
    // Head
    const headGeo = new THREE.ConeGeometry(headWidth, headLength, 8)
    const headMat = new THREE.MeshBasicMaterial({ color })
    const head = new THREE.Mesh(headGeo, headMat)
    head.position.copy(to.clone().sub(direction.clone().multiplyScalar(headLength / 2)))
    head.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
    group.add(head)
    
    return group
  }

  // Calculate velocities for display
  const velocityInfo = useMemo(() => {
    // Local Standard of Rest velocity (Sun's motion)
    const v_LSR = Math.sqrt(
      GALACTIC_CONSTANTS.V_SUN_PECULIAR_U_kms ** 2 +
      GALACTIC_CONSTANTS.V_SUN_PECULIAR_V_kms ** 2 +
      GALACTIC_CONSTANTS.V_SUN_PECULIAR_W_kms ** 2
    )
    
    // Galaxy rotation at Sun's position
    const v_galactic_rotation = GALACTIC_CONSTANTS.V_CIRCULAR_SUN_kms
    
    // HVCC velocity (from V_LSR)
    const v_hvcc = Math.abs(config.lsrVelocity)
    
    // Current cloud velocity magnitude
    const v_current = currentCloudVelocity
      ? Math.sqrt(
          currentCloudVelocity[0] ** 2 +
          currentCloudVelocity[1] ** 2 +
          currentCloudVelocity[2] ** 2
        )
      : Math.sqrt(
          config.initialVelocity[0] ** 2 +
          config.initialVelocity[1] ** 2 +
          config.initialVelocity[2] ** 2
        )
    
    return {
      v_LSR: v_LSR.toFixed(1),
      v_galactic_rotation: v_galactic_rotation.toFixed(0),
      v_hvcc: v_hvcc.toFixed(0),
      v_current: v_current.toFixed(2),
    }
  }, [config, currentCloudVelocity])

  return (
    <div className={`relative ${className}`}>
      <div 
        ref={containerRef} 
        style={{ width, height }}
        className="rounded-lg overflow-hidden"
      />
      
      {/* Info panel */}
      <div className="absolute top-2 left-2 bg-black/80 text-white text-xs p-2 rounded">
        <div className="font-bold mb-1 text-cyan-300">● IMBH-Centered Orbital View</div>
        <div className="text-gray-300">Oka et al. (2017) / CAT_OKA</div>
        <div className="mt-1 border-t border-gray-600 pt-1">
          <div className="text-yellow-300">☉ Sun/⊕ Earth = Observer</div>
          <div className="text-gray-400">d(GC) ≈ {config.distanceToGC} kpc</div>
        </div>
        <div className="mt-1 border-t border-gray-600 pt-1 text-[10px]">
          <div className="text-orange-300">Viewing Geometry:</div>
          <div className="text-gray-400">i=0° face-on, i=90° edge-on</div>
          <div className="text-gray-400">PA from N through E</div>
        </div>
      </div>
      
      {/* Velocity Summary */}
      <div className="absolute top-2 right-2 bg-black/80 text-white text-xs p-2 rounded">
        <div className="font-bold mb-1 text-green-300">Velocities (km/s)</div>
        <table className="w-full text-[10px]">
          <tbody>
            <tr>
              <td className="text-cyan-400">V_circular(Sun):</td>
              <td className="text-right">{velocityInfo.v_galactic_rotation}</td>
            </tr>
            <tr>
              <td className="text-yellow-400">V_peculiar(Sun):</td>
              <td className="text-right">{velocityInfo.v_LSR}</td>
            </tr>
            <tr>
              <td className="text-pink-400">V_LSR(HVCC):</td>
              <td className="text-right">{velocityInfo.v_hvcc}</td>
            </tr>
            <tr>
              <td className="text-blue-400">|v_cloud|(t):</td>
              <td className="text-right font-bold">{velocityInfo.v_current}</td>
            </tr>
          </tbody>
        </table>
      </div>
      
      {/* Parameter summary */}
      <div className="absolute bottom-2 right-2 bg-black/80 text-white text-xs p-2 rounded max-w-xs">
        <div className="font-bold mb-1 text-orange-300">Orbital Parameters</div>
        <table className="w-full">
          <tbody>
            <tr>
              <td className="text-yellow-400">r_peri:</td>
              <td className="text-right">{config.pericentre.toFixed(2)} pc</td>
            </tr>
            <tr>
              <td className="text-orange-400">b:</td>
              <td className="text-right">{config.impactParameter.toFixed(2)} pc</td>
            </tr>
            <tr>
              <td className="text-green-400">e:</td>
              <td className="text-right">{config.eccentricity.toFixed(3)}</td>
            </tr>
            <tr>
              <td className="text-purple-400">i:</td>
              <td className="text-right">{config.inclination}°</td>
            </tr>
            <tr>
              <td className="text-teal-400">PA:</td>
              <td className="text-right">{config.positionAngle}°</td>
            </tr>
            <tr>
              <td className="text-pink-400">V_LSR:</td>
              <td className="text-right">{config.lsrVelocity} km/s</td>
            </tr>
          </tbody>
        </table>
      </div>
      
      {/* Simulation Time */}
      <div className="absolute bottom-2 left-2 bg-black/80 text-white text-xs p-2 rounded">
        <div className="text-cyan-300 font-bold">
          t = {(currentTime * 0.978).toFixed(3)} Myr
        </div>
        <div className="text-gray-400 text-[10px]">
          ({currentTime.toFixed(3)} code)
        </div>
      </div>
    </div>
  )
}

export default OrbitalGeometryPanel
