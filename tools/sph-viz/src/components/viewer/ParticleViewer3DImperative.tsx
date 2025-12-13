'use client'

import { useRef, useEffect, useCallback, useState } from 'react'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'
import type { ParsedFrame, ColorMap } from '~/types/sph'
import {
  createAllGalacticMarkers,
  updateGalaxyRotation,
  updateVCircularAnimation,
  updateLSRRotation,
  type GalacticMarkersConfig,
  DEFAULT_GALACTIC_CONFIG,
} from './GalacticMarkers'
// Shared utilities (SSOT)
import { getCircleTexture, createAxesWithLabels } from '~/utils/three-helpers'
import { sampleColorMap, COLOR_MAP_DATA } from '~/utils/color-interpolation'

// ============== Physical Constants & Unit Conversions ==============
// Units: IMBH_ENCOUNTER ([L]=pc, [M]=1000 M_sun, [V]=km/s, [T]=0.978 Myr)
const UNITS = {
  TIME: 0.978,          // 1 code time = 0.978 Myr
  VELOCITY: 1.0,        // km/s
  LENGTH: 1.0,          // pc
  MASS: 1000.0,         // M_sun (1 code mass = 1000 M_sun)
  // G = 1 in code units: [G] = [L]³/([M][T]²) where T = L/V
  // Verified: GM = 100 (code) gives v_circ ≈ 10 km/s at r = 1 pc
  G: 1.0,
}

/**
 * Compute hyperbolic orbit trajectory for the cloud's center of mass.
 *
 * IMPORTANT: Orbital parameters are computed directly from the initial state
 * (position and velocity) using orbital mechanics formulas. This ensures the
 * analytical orbit matches the actual simulation trajectory.
 *
 * GEOMETRY SETUP (IMBH_ENCOUNTER coordinates):
 * ============================================
 * - BH is at origin (0, 0, 0)
 * - Cloud starts at position (x0, y0, z0), typically (20, -5.17, 0) pc
 * - Cloud has initial velocity (vx, vy, vz), typically (-10.18, 5.05, 0) km/s
 * - Orbit lies in the x-y plane (z ≈ 0)
 *
 * ORBITAL MECHANICS FORMULAS:
 * ===========================
 * Specific orbital energy:    ε = v²/2 - GM/r
 * Specific angular momentum:  h = |r × v|
 * Semi-latus rectum:          p = h² / (GM)
 * Eccentricity:               e = √(1 + 2εh²/(GM)²)
 * Semi-major axis:            a = GM / (2ε)
 * Pericentre:                 rp = a(e - 1)
 */
function computeHyperbolicOrbit(
  config: {
    bhPosition: [number, number, number]
    bhMass: number  // in code units (1 = 1000 M_sun), so GM = bhMass since G=1
    cloudInitialPosition: [number, number, number]
    cloudInitialVelocity: [number, number, number]
  },
  numPoints: number = 200
): [number, number, number][] {
  const { bhPosition, bhMass, cloudInitialPosition, cloudInitialVelocity } = config

  const [x0, y0, z0] = cloudInitialPosition
  const [vx0, vy0, vz0] = cloudInitialVelocity
  const [bh_x, bh_y, bh_z] = bhPosition

  // Position and velocity relative to BH
  const rx = x0 - bh_x
  const ry = y0 - bh_y
  const r0 = Math.sqrt(rx * rx + ry * ry)
  const v0 = Math.sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0)

  // GM in code units (G=1, M=bhMass)
  const GM = bhMass

  // ========== Compute Orbital Parameters from State Vectors ==========
  // Specific orbital energy: ε = v²/2 - GM/r
  const epsilon = v0 * v0 / 2 - GM / r0

  // Check for hyperbolic orbit (ε > 0)
  if (epsilon <= 0) {
    console.warn('[computeHyperbolicOrbit] Non-hyperbolic orbit (ε <= 0), using straight line')
    return computeStraightLineTrajectory({ cloudInitialPosition, cloudInitialVelocity }, -0.5, 5, numPoints)
  }

  // Specific angular momentum (z-component for 2D orbit): h = x*vy - y*vx
  const h = Math.abs(rx * vy0 - ry * vx0)

  // Semi-latus rectum: p = h² / (GM)
  const p = h * h / GM

  // Eccentricity: e = √(1 + 2εh²/(GM)²)
  const e = Math.sqrt(1 + 2 * epsilon * h * h / (GM * GM))

  // Semi-major axis: a = GM / (2ε)
  const a = GM / (2 * epsilon)

  // Pericentre: rp = a(e - 1)
  const r_peri = a * (e - 1)

  // ========== Orbit Orientation ==========
  // Find the angle of the initial position
  const theta0 = Math.atan2(ry, rx)

  // From the orbit equation r = p / (1 + e·cos(ν)), solve for ν at r = r0
  const cosNu0 = Math.max(-0.9999, Math.min(0.9999, (p / r0 - 1) / e))

  // Determine sign of ν from radial velocity (r·v)
  // If r·v < 0, cloud is approaching (ν < 0)
  const rdotv = rx * vx0 + ry * vy0
  const nu0 = rdotv < 0 ? -Math.acos(cosNu0) : Math.acos(cosNu0)

  // Argument of pericentre: ω = θ0 - ν0
  const omega = theta0 - nu0

  // ========== Generate Orbit Points ==========
  // For hyperbola: ν ranges from -(π - acos(1/e)) to +(π - acos(1/e))
  const nu_asymp = Math.acos(-1 / e)
  const nu_max = nu_asymp - 0.05
  const nu_min = -nu_max

  const orbit: [number, number, number][] = []

  for (let i = 0; i < numPoints; i++) {
    const nu = nu_min + (nu_max - nu_min) * (i / (numPoints - 1))
    const denom = 1 + e * Math.cos(nu)
    if (denom <= 0) continue
    const r = p / denom
    const theta = omega + nu
    const x = bh_x + r * Math.cos(theta)
    const y = bh_y + r * Math.sin(theta)
    const z = bh_z

    if (Math.abs(x) < 100 && Math.abs(y) < 100 && r > 0) {
      orbit.push([x, y, z])
    }
  }

  return orbit
}

/**
 * Compute straight-line approximation trajectory.
 * Valid when v∞ >> v_escape (the cloud isn't strongly bound).
 * 
 * x(t) = x₀ + vₓ·t
 * y(t) = y₀ (constant for straight line)
 * z(t) = z₀
 */
function computeStraightLineTrajectory(
  config: {
    cloudInitialPosition: [number, number, number]
    cloudInitialVelocity: [number, number, number]
  },
  t_start: number,
  t_end: number,
  numPoints: number = 50
): [number, number, number][] {
  const { cloudInitialPosition, cloudInitialVelocity } = config
  const [x0, y0, z0] = cloudInitialPosition
  const [vx, vy, vz] = cloudInitialVelocity
  
  const trajectory: [number, number, number][] = []
  
  for (let i = 0; i < numPoints; i++) {
    const t = t_start + (t_end - t_start) * (i / (numPoints - 1))
    trajectory.push([
      x0 + vx * t,
      y0 + vy * t,
      z0 + vz * t
    ])
  }
  
  return trajectory
}

// Note: createCircleTexture, getCircleTexture, COLOR_MAP_DATA, and sampleColorMap
// are now imported from shared utilities (~/utils/three-helpers and ~/utils/color-interpolation)

/** IMBH physics configuration for visualization */
export interface IMBHPhysicsConfig {
  enabled: boolean
  /** Black hole position [x, y, z] in code units */
  bhPosition: [number, number, number]
  /** Black hole mass in code units (100 = 1e5 Msun) */
  bhMass: number
  /** Cloud initial position [x, y, z] in code units */
  cloudInitialPosition: [number, number, number]
  /** Cloud initial velocity [vx, vy, vz] in code units (km/s) */
  cloudInitialVelocity: [number, number, number]
  /** Cloud mass in code units (1 = 1000 Msun) */
  cloudMass: number
  /** Cloud radius in pc */
  cloudRadius: number
  /** Tidal radius in pc (pre-computed or computed from masses) */
  tidalRadius: number
  /** Time unit in Myr */
  timeUnit: number
  /** Impact parameter in pc */
  impactParameter?: number
  /** Pericentre distance in pc (from preset config) */
  pericentre?: number
  /** Orbital eccentricity (>1 for hyperbolic, from preset config) */
  eccentricity?: number
}

/** Galactic coordinate configuration for markers */
export interface GalacticConfig {
  /** Distance from Earth/Sun to Galactic Center (kpc) */
  distanceToGC: number
  /** Galactic longitude of the source (degrees) */
  galacticLongitude: number
  /** Galactic latitude of the source (degrees) */
  galacticLatitude: number
  /** Orbital plane inclination (degrees) */
  inclination: number
  /** Position angle on sky (degrees) */
  positionAngle: number
  /** V_LSR in km/s (negative = approaching) */
  lsrVelocity: number
  /** Show low-opacity Milky Way disk visualization */
  showGalaxyDisk?: boolean
  /** Show Solar System demo at Sun position (enlarged for visibility) */
  showSolarSystem?: boolean
  /** Rotation speed for galaxy animation (rad/s, 0 = static) */
  galaxyRotationSpeed?: number
}

export interface ParticleViewer3DImperativeProps {
  /** Map of all loaded frames - passed once, not on every frame change */
  framesRef: React.RefObject<Map<number, ParsedFrame>>
  /** Initial frame index */
  initialFrameIndex?: number
  /** Callback to get current frame index imperatively */
  frameIndexRef: React.RefObject<number>
  colorField?: string
  colorMapName?: string
  pointSize?: number
  opacity?: number
  logScale?: boolean
  showAxes?: boolean
  showBoundingBox?: boolean
  boundingBox?: { min: [number, number, number]; max: [number, number, number] }
  className?: string
  /** Called when FPS changes */
  onFpsUpdate?: (fps: number) => void
  /** Global color range for consistent coloring across frames [min, max] */
  globalColorRange?: [number, number]
  /** IMBH physics configuration for black hole and orbital features */
  imbhPhysics?: IMBHPhysicsConfig
  /** Show the black hole */
  showBlackHole?: boolean
  /** Show orbital trajectory */
  showTrajectory?: boolean
  /** Show tidal and Hill radii circles */
  showRadii?: boolean
  /** Show galactic coordinate markers (Sun, Earth, GC, LoS, HVCC) */
  showGalacticMarkers?: boolean
  /** Show text labels on markers and arrows */
  showLabels?: boolean
  /** Galactic coordinate configuration */
  galacticConfig?: GalacticConfig
  /** Camera view mode: 'free' = orbit controls, 'earth' = view from Earth toward BH */
  cameraMode?: 'free' | 'earth'
  /** Callback when camera mode changes (for UI sync) */
  onCameraModeChange?: (mode: 'free' | 'earth') => void
}

/**
 * High-performance 3D particle viewer using imperative Three.js
 * 
 * Key differences from React Three Fiber version:
 * 1. Three.js owns the render loop - no React re-renders during animation
 * 2. Frame data accessed via refs, not props
 * 3. Buffer updates happen directly, bypassing React's reconciliation
 * 4. Can achieve 120+ FPS on capable hardware
 */
export function ParticleViewer3DImperative({
  framesRef,
  initialFrameIndex = 0,
  frameIndexRef,
  colorField = 'density',
  colorMapName = 'viridis',
  pointSize = 2.0,
  opacity = 0.9,
  logScale = false,
  showAxes = true,
  showBoundingBox = true,
  boundingBox,
  className = '',
  onFpsUpdate,
  globalColorRange,
  imbhPhysics,
  showBlackHole = true,
  showTrajectory = true,
  showRadii = true,
  showGalacticMarkers = true,
  showLabels = true,
  galacticConfig,
  cameraMode = 'free',
  onCameraModeChange,
}: ParticleViewer3DImperativeProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null)
  const sceneRef = useRef<THREE.Scene | null>(null)
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null)
  const controlsRef = useRef<OrbitControls | null>(null)
  const particlesRef = useRef<THREE.Points | null>(null)
  const geometryRef = useRef<THREE.BufferGeometry | null>(null)
  const materialRef = useRef<THREE.PointsMaterial | null>(null)
  
  // IMBH visualization objects
  const bhGroupRef = useRef<THREE.Group | null>(null)
  const trajectoryLineRef = useRef<THREE.Line | null>(null)
  const trajectoryPointsRef = useRef<[number, number, number][]>([])
  const comMarkerRef = useRef<THREE.Mesh | null>(null)
  const tidalCircleRef = useRef<THREE.Line | null>(null)
  const hillCircleRef = useRef<THREE.Line | null>(null)
  
  // Galactic markers group
  const galacticMarkersGroupRef = useRef<THREE.Group | null>(null)
  
  // Track the last rendered frame to avoid redundant updates
  const lastFrameIndexRef = useRef<number>(-1)
  const lastColorFieldRef = useRef<string>(colorField)
  const lastColorMapRef = useRef<string>(colorMapName)
  const lastLogScaleRef = useRef<boolean>(logScale)
  
  // Store global color range in ref for render loop access
  const globalColorRangeRef = useRef<[number, number] | undefined>(globalColorRange)
  useEffect(() => {
    globalColorRangeRef.current = globalColorRange
  }, [globalColorRange])
  
  // Pre-computed stats for color normalization
  const statsRef = useRef<{
    density: [number, number]
    velocity: [number, number]
    pressure: [number, number]
    energy: [number, number]
  }>({ density: [0, 1], velocity: [0, 1], pressure: [0, 1], energy: [0, 1] })
  
  // FPS tracking
  const fpsRef = useRef({ frames: 0, lastTime: performance.now(), fps: 0 })

  // Compute global stats from all loaded frames
  const computeStats = useCallback(() => {
    const frames = framesRef.current
    if (!frames || frames.size === 0) return

    const allDensity: number[] = []
    const allVelocity: number[] = []
    const allPressure: number[] = []
    const allEnergy: number[] = []

    // Sample every 10th frame and every 100th particle for efficiency
    frames.forEach((frame, idx) => {
      if (idx % 10 !== 0) return
      for (let i = 0; i < frame.particleCount; i += 100) {
        allDensity.push(frame.density[i])
        const vx = frame.velocities[i * 3]
        const vy = frame.velocities[i * 3 + 1]
        const vz = frame.velocities[i * 3 + 2]
        allVelocity.push(Math.sqrt(vx * vx + vy * vy + vz * vz))
        allPressure.push(frame.pressure[i])
        allEnergy.push(frame.energy[i])
      }
    })

    const percentile = (arr: number[], p: number) => {
      const sorted = arr.slice().sort((a, b) => a - b)
      return sorted[Math.floor(sorted.length * p / 100)] || 0
    }

    statsRef.current = {
      density: [percentile(allDensity, 1), percentile(allDensity, 99)],
      velocity: [percentile(allVelocity, 1), percentile(allVelocity, 99)],
      pressure: [percentile(allPressure, 1), percentile(allPressure, 99)],
      energy: [percentile(allEnergy, 1), percentile(allEnergy, 99)],
    }
  }, [framesRef])

  // Update particle buffers - called from render loop, NOT from React
  const updateParticles = useCallback(() => {
    const frameIndex = frameIndexRef.current ?? 0
    const frames = framesRef.current
    if (!frames || !geometryRef.current) return

    const frame = frames.get(frameIndex)
    if (!frame) return

    // Check if we need to update
    const needsUpdate = 
      frameIndex !== lastFrameIndexRef.current ||
      colorField !== lastColorFieldRef.current ||
      colorMapName !== lastColorMapRef.current ||
      logScale !== lastLogScaleRef.current

    if (!needsUpdate) return

    lastFrameIndexRef.current = frameIndex
    lastColorFieldRef.current = colorField
    lastColorMapRef.current = colorMapName
    lastLogScaleRef.current = logScale

    const geometry = geometryRef.current
    const posAttr = geometry.attributes.position as THREE.BufferAttribute
    const colorAttr = geometry.attributes.color as THREE.BufferAttribute

    // Check if we need to resize buffers
    if (posAttr.count !== frame.particleCount) {
      const positions = new Float32Array(frame.particleCount * 3)
      const colors = new Float32Array(frame.particleCount * 3)
      geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage))
      geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3).setUsage(THREE.DynamicDrawUsage))
    }

    // Update positions
    const positions = (geometry.attributes.position as THREE.BufferAttribute).array as Float32Array
    positions.set(frame.positions)
    ;(geometry.attributes.position as THREE.BufferAttribute).needsUpdate = true

    // Get field data for coloring
    let fieldData: Float32Array
    let vMin: number, vMax: number

    // Use global color range if provided, otherwise fall back to internal stats
    const useGlobalRange = globalColorRangeRef.current && globalColorRangeRef.current[0] !== globalColorRangeRef.current[1]
    
    switch (colorField) {
      case 'velocity': {
        fieldData = new Float32Array(frame.particleCount)
        for (let i = 0; i < frame.particleCount; i++) {
          const vx = frame.velocities[i * 3]
          const vy = frame.velocities[i * 3 + 1]
          const vz = frame.velocities[i * 3 + 2]
          fieldData[i] = Math.sqrt(vx * vx + vy * vy + vz * vz)
        }
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current!
        } else {
          [vMin, vMax] = statsRef.current.velocity
        }
        break
      }
      case 'pressure':
        fieldData = frame.pressure
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current!
        } else {
          [vMin, vMax] = statsRef.current.pressure
        }
        break
      case 'energy':
        fieldData = frame.energy
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current!
        } else {
          [vMin, vMax] = statsRef.current.energy
        }
        break
      default:
        fieldData = frame.density
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current!
        } else {
          [vMin, vMax] = statsRef.current.density
        }
    }

    // Update colors
    const colors = (geometry.attributes.color as THREE.BufferAttribute).array as Float32Array
    const logMin = logScale && vMin > 0 ? Math.log10(vMin) : 0
    const logRange = logScale && vMin > 0 ? Math.log10(vMax) - logMin : 1
    const range = vMax - vMin || 1

    for (let i = 0; i < frame.particleCount; i++) {
      let val = fieldData[i]
      if (!isFinite(val)) val = vMin

      let t: number
      if (logScale && vMin > 0) {
        t = (Math.log10(Math.max(val, vMin)) - logMin) / logRange
      } else {
        t = (val - vMin) / range
      }
      t = Math.max(0, Math.min(1, t))

      const [r, g, b] = sampleColorMap(colorMapName, t)
      colors[i * 3] = r
      colors[i * 3 + 1] = g
      colors[i * 3 + 2] = b
    }

    ;(geometry.attributes.color as THREE.BufferAttribute).needsUpdate = true
    geometry.computeBoundingSphere()
  }, [framesRef, frameIndexRef, colorField, colorMapName, logScale])

  // Update IMBH trajectory and COM marker
  const updateTrajectory = useCallback(() => {
    if (!imbhPhysics?.enabled) return
    const frames = framesRef?.current
    const frameIndex = frameIndexRef?.current ?? 0
    if (!frames || frames.size === 0) return

    const frame = frames.get(frameIndex)
    if (!frame) return

    // Reset trajectory if we're at frame 0
    if (frameIndex === 0 && trajectoryPointsRef.current.length > 1) {
      trajectoryPointsRef.current = []
      if (trajectoryLineRef.current) {
        trajectoryLineRef.current.geometry.setDrawRange(0, 0)
      }
    }

    // Compute center of mass
    let totalMass = 0
    let comX = 0, comY = 0, comZ = 0
    
    for (let i = 0; i < frame.particleCount; i++) {
      const mass = frame.mass?.[i] ?? 1.0
      totalMass += mass
      comX += mass * frame.positions[i * 3]
      comY += mass * frame.positions[i * 3 + 1]
      comZ += mass * frame.positions[i * 3 + 2]
    }
    
    if (totalMass > 0) {
      comX /= totalMass
      comY /= totalMass
      comZ /= totalMass
    }

    // Update COM marker position
    if (comMarkerRef.current) {
      comMarkerRef.current.position.set(comX, comY, comZ)
      comMarkerRef.current.visible = showTrajectory
    }

    // Add to trajectory if this is a new frame
    const lastIdx = trajectoryPointsRef.current.length > 0 
      ? trajectoryPointsRef.current[trajectoryPointsRef.current.length - 1]
      : null
    
    // Only add if position changed significantly (avoid duplicate points)
    const threshold = 0.01
    if (!lastIdx || 
        Math.abs(comX - lastIdx[0]) > threshold || 
        Math.abs(comY - lastIdx[1]) > threshold || 
        Math.abs(comZ - lastIdx[2]) > threshold) {
      
      trajectoryPointsRef.current.push([comX, comY, comZ])
      
      // Limit trajectory length
      const maxPoints = 500
      if (trajectoryPointsRef.current.length > maxPoints) {
        trajectoryPointsRef.current = trajectoryPointsRef.current.slice(-maxPoints)
      }

      // Update trajectory line geometry
      if (trajectoryLineRef.current && trajectoryPointsRef.current.length >= 2) {
        const positions = new Float32Array(trajectoryPointsRef.current.length * 3)
        for (let i = 0; i < trajectoryPointsRef.current.length; i++) {
          positions[i * 3] = trajectoryPointsRef.current[i][0]
          positions[i * 3 + 1] = trajectoryPointsRef.current[i][1]
          positions[i * 3 + 2] = trajectoryPointsRef.current[i][2]
        }
        trajectoryLineRef.current.geometry.setAttribute(
          'position',
          new THREE.BufferAttribute(positions, 3)
        )
        trajectoryLineRef.current.geometry.setDrawRange(0, trajectoryPointsRef.current.length)
        ;(trajectoryLineRef.current.geometry.attributes.position as THREE.BufferAttribute).needsUpdate = true
      }
    }

    // Update Hill radius circle position to follow COM
    if (hillCircleRef.current) {
      hillCircleRef.current.position.set(comX, comY, 0)
    }
  }, [framesRef, frameIndexRef, imbhPhysics, showTrajectory])

  // Initialize Three.js scene
  useEffect(() => {
    if (!containerRef.current) return

    const container = containerRef.current
    const width = container.clientWidth
    const height = container.clientHeight

    // Scene
    const scene = new THREE.Scene()
    scene.background = new THREE.Color(0x0a0a0f)
    sceneRef.current = scene

    // Camera
    const camera = new THREE.PerspectiveCamera(60, width / height, 0.1, 10000)
    camera.position.set(50, 50, 100)
    cameraRef.current = camera

    // Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true })
    renderer.setSize(width, height)
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    container.appendChild(renderer.domElement)
    rendererRef.current = renderer

    // Controls
    const controls = new OrbitControls(camera, renderer.domElement)
    controls.enableDamping = true
    controls.dampingFactor = 0.05
    controls.screenSpacePanning = true
    controls.minDistance = 1
    controls.maxDistance = 5000
    controlsRef.current = controls

    // Grid helper with physical units (positions in pc)
    const gridHelper = new THREE.GridHelper(200, 20, 0x444444, 0x222222)
    scene.add(gridHelper)

    // Axes helper with physical unit labels (using shared utility with LARGER labels)
    if (showAxes) {
      const axesGroup = createAxesWithLabels({
        size: 30,
        showTickMarks: true,
        tickInterval: 10,
        labelUnit: 'pc'
      })
      scene.add(axesGroup)
    }

    // Initialize particle system with empty buffers
    const geometry = new THREE.BufferGeometry()
    const initialCount = 64000 // Pre-allocate for expected particle count
    const positions = new Float32Array(initialCount * 3)
    const colors = new Float32Array(initialCount * 3)
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage))
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3).setUsage(THREE.DynamicDrawUsage))
    geometryRef.current = geometry

    const material = new THREE.PointsMaterial({
      size: pointSize,
      vertexColors: true,
      sizeAttenuation: true,
      transparent: true,
      opacity: opacity,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      map: getCircleTexture(),
      alphaTest: 0.01,
    })
    materialRef.current = material

    const particles = new THREE.Points(geometry, material)
    scene.add(particles)
    particlesRef.current = particles

    // =========================================================================
    // IMBH PHYSICS VISUALIZATION
    // =========================================================================
    if (imbhPhysics?.enabled) {
      const bhPos = imbhPhysics.bhPosition
      
      // --- Black Hole Group ---
      if (showBlackHole) {
        const bhGroup = new THREE.Group()
        bhGroup.position.set(bhPos[0], bhPos[1], bhPos[2])
        
        // Glow effect - multiple concentric spheres
        const glowColors = [0xff4444, 0xff6666, 0xff8888, 0xffaaaa]
        const glowSizes = [0.8, 0.6, 0.4, 0.25]
        glowSizes.forEach((size, i) => {
          const glowGeo = new THREE.SphereGeometry(size, 32, 32)
          const glowMat = new THREE.MeshBasicMaterial({
            color: glowColors[i],
            transparent: true,
            opacity: 0.15 - i * 0.03,
          })
          const glowMesh = new THREE.Mesh(glowGeo, glowMat)
          bhGroup.add(glowMesh)
        })
        
        // Core - bright white/yellow center
        const coreGeo = new THREE.SphereGeometry(0.15, 32, 32)
        const coreMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
        const coreMesh = new THREE.Mesh(coreGeo, coreMat)
        bhGroup.add(coreMesh)
        
        // Event horizon representation - dark sphere
        const horizonGeo = new THREE.SphereGeometry(0.08, 32, 32)
        const horizonMat = new THREE.MeshBasicMaterial({ color: 0x000000 })
        const horizonMesh = new THREE.Mesh(horizonGeo, horizonMat)
        bhGroup.add(horizonMesh)
        
        scene.add(bhGroup)
        bhGroupRef.current = bhGroup
      }
      
      // --- Tidal Radius Circle ---
      if (showRadii && imbhPhysics.tidalRadius > 0) {
        const tidalRadius = imbhPhysics.tidalRadius
        const circlePoints: THREE.Vector3[] = []
        for (let i = 0; i <= 64; i++) {
          const theta = (i / 64) * Math.PI * 2
          circlePoints.push(new THREE.Vector3(
            bhPos[0] + tidalRadius * Math.cos(theta),
            bhPos[1] + tidalRadius * Math.sin(theta),
            bhPos[2]
          ))
        }
        const tidalGeo = new THREE.BufferGeometry().setFromPoints(circlePoints)
        const tidalMat = new THREE.LineDashedMaterial({ 
          color: 0x00ffff, 
          dashSize: 0.5, 
          gapSize: 0.25,
          transparent: true,
          opacity: 0.6 
        })
        const tidalLine = new THREE.Line(tidalGeo, tidalMat)
        tidalLine.computeLineDistances()
        scene.add(tidalLine)
        tidalCircleRef.current = tidalLine
        
        // Hill radius (roughly 0.4 * tidal radius for typical parameters)
        const hillRadius = tidalRadius * 0.4
        const hillPoints: THREE.Vector3[] = []
        for (let i = 0; i <= 64; i++) {
          const theta = (i / 64) * Math.PI * 2
          hillPoints.push(new THREE.Vector3(
            bhPos[0] + hillRadius * Math.cos(theta),
            bhPos[1] + hillRadius * Math.sin(theta),
            bhPos[2]
          ))
        }
        const hillGeo = new THREE.BufferGeometry().setFromPoints(hillPoints)
        const hillMat = new THREE.LineDashedMaterial({ 
          color: 0xff8800, 
          dashSize: 0.3, 
          gapSize: 0.15,
          transparent: true,
          opacity: 0.5 
        })
        const hillLine = new THREE.Line(hillGeo, hillMat)
        hillLine.computeLineDistances()
        scene.add(hillLine)
        hillCircleRef.current = hillLine
      }
      
      // --- Analytical Trajectory (Hyperbolic Orbit from Orbital Mechanics) ---
      // Orbit equation: r(θ) = a(e²-1) / (1 + e·cos(θ))
      // where: a = GM/v∞², e = √(1 + (b·v∞²/GM)²)
      // When pericentre and eccentricity are provided from preset, use them directly
      if (showTrajectory) {
        // Compute hyperbolic orbit using analytical solution
        // Orbital parameters are computed from initial position and velocity
        const orbitPoints = computeHyperbolicOrbit({
          bhPosition: imbhPhysics.bhPosition,
          bhMass: imbhPhysics.bhMass,
          cloudInitialPosition: imbhPhysics.cloudInitialPosition,
          cloudInitialVelocity: imbhPhysics.cloudInitialVelocity,
        }, 100)
        
        // Create orbit line
        const orbitVectors = orbitPoints.map(p => new THREE.Vector3(p[0], p[1], p[2]))
        const orbitGeo = new THREE.BufferGeometry().setFromPoints(orbitVectors)
        const orbitMat = new THREE.LineDashedMaterial({
          color: 0x00ff88,  // Green for analytical solution
          dashSize: 0.8,
          gapSize: 0.4,
          transparent: true,
          opacity: 0.6,
        })
        const orbitLine = new THREE.Line(orbitGeo, orbitMat)
        orbitLine.computeLineDistances()
        scene.add(orbitLine)
        
        // Also draw straight-line approximation for comparison (gray, fainter)
        const straightLine = computeStraightLineTrajectory(
          {
            cloudInitialPosition: imbhPhysics.cloudInitialPosition,
            cloudInitialVelocity: imbhPhysics.cloudInitialVelocity,
          },
          0, 4.0, 50  // t from 0 to 4 code time units
        )
        const straightVectors = straightLine.map(p => new THREE.Vector3(p[0], p[1], p[2]))
        const straightGeo = new THREE.BufferGeometry().setFromPoints(straightVectors)
        const straightMat = new THREE.LineDashedMaterial({
          color: 0x888888,
          dashSize: 0.5,
          gapSize: 0.3,
          transparent: true,
          opacity: 0.3,
        })
        const straightLineObj = new THREE.Line(straightGeo, straightMat)
        straightLineObj.computeLineDistances()
        scene.add(straightLineObj)
      }
      
      // --- Trajectory Line (will be updated each frame) ---
      if (showTrajectory) {
        const trajectoryGeo = new THREE.BufferGeometry()
        // Pre-allocate for up to 200 trajectory points
        const trajectoryPositions = new Float32Array(200 * 3)
        trajectoryGeo.setAttribute('position', new THREE.BufferAttribute(trajectoryPositions, 3).setUsage(THREE.DynamicDrawUsage))
        trajectoryGeo.setDrawRange(0, 0)
        
        const trajectoryMat = new THREE.LineBasicMaterial({ 
          color: 0xffaa00, 
          transparent: true, 
          opacity: 0.8,
          linewidth: 2 
        })
        const trajectoryLine = new THREE.Line(trajectoryGeo, trajectoryMat)
        scene.add(trajectoryLine)
        trajectoryLineRef.current = trajectoryLine
        
        // Center of mass marker
        const comGeo = new THREE.SphereGeometry(0.2, 16, 16)
        const comMat = new THREE.MeshBasicMaterial({ color: 0xffaa00 })
        const comMarker = new THREE.Mesh(comGeo, comMat)
        comMarker.visible = false
        scene.add(comMarker)
        comMarkerRef.current = comMarker
      }
    }
    // =========================================================================

    // =========================================================================
    // GALACTIC COORDINATE MARKERS (Sun, Earth, LoS, GC, HVCC)
    // =========================================================================
    if (showGalacticMarkers && imbhPhysics?.enabled) {
      const bhPos = new THREE.Vector3(
        imbhPhysics.bhPosition[0],
        imbhPhysics.bhPosition[1],
        imbhPhysics.bhPosition[2]
      )
      
      // Build galactic config from props or use defaults
      const galConfig: GalacticMarkersConfig = {
        distanceToGC: galacticConfig?.distanceToGC ?? DEFAULT_GALACTIC_CONFIG.distanceToGC,
        galacticLongitude: galacticConfig?.galacticLongitude ?? DEFAULT_GALACTIC_CONFIG.galacticLongitude,
        galacticLatitude: galacticConfig?.galacticLatitude ?? DEFAULT_GALACTIC_CONFIG.galacticLatitude,
        inclination: galacticConfig?.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination,
        positionAngle: galacticConfig?.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle,
        lsrVelocity: galacticConfig?.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity,
        hvccPosition: DEFAULT_GALACTIC_CONFIG.hvccPosition,
        hvccSize: DEFAULT_GALACTIC_CONFIG.hvccSize,
        // Pass cloud radius from imbhPhysics if available
        cloudRadius_pc: imbhPhysics?.cloudRadius ?? DEFAULT_GALACTIC_CONFIG.cloudRadius_pc,
        // Galaxy disk visualization options
        showGalaxyDisk: galacticConfig?.showGalaxyDisk ?? false,  // Off by default for performance
        showSolarSystem: galacticConfig?.showSolarSystem ?? false,  // Off by default
        galaxyRotationSpeed: galacticConfig?.galaxyRotationSpeed ?? 0.1,  // Slow rotation for demo
      }

      const galacticMarkers = createAllGalacticMarkers(bhPos, galConfig)
      galacticMarkers.visible = showGalacticMarkers
      scene.add(galacticMarkers)
      galacticMarkersGroupRef.current = galacticMarkers
    }
    // =========================================================================

    // Compute initial stats
    computeStats()

    // Fit camera to bounding box
    if (boundingBox) {
      const center = new THREE.Vector3(
        (boundingBox.min[0] + boundingBox.max[0]) / 2,
        (boundingBox.min[1] + boundingBox.max[1]) / 2,
        (boundingBox.min[2] + boundingBox.max[2]) / 2
      )
      const size = Math.max(
        boundingBox.max[0] - boundingBox.min[0],
        boundingBox.max[1] - boundingBox.min[1],
        boundingBox.max[2] - boundingBox.min[2]
      )
      camera.position.set(center.x + size * 2, center.y + size, center.z + size * 2)
      camera.lookAt(center)
      controls.target.copy(center)
    }

    // Animation loop - THIS IS KEY: runs independently of React
    let animationId: number
    let lastTime = performance.now()
    
    const animate = () => {
      animationId = requestAnimationFrame(animate)
      
      // Calculate delta time for animations
      const now = performance.now()
      const deltaTime = (now - lastTime) / 1000  // Convert to seconds
      lastTime = now

      // FPS calculation
      fpsRef.current.frames++
      if (now - fpsRef.current.lastTime >= 1000) {
        fpsRef.current.fps = fpsRef.current.frames
        fpsRef.current.frames = 0
        fpsRef.current.lastTime = now
        onFpsUpdate?.(fpsRef.current.fps)
      }

      // Update particles from current frame ref
      updateParticles()
      
      // Update IMBH trajectory
      updateTrajectory()
      
      // Update galaxy disk rotation if present
      if (galacticMarkersGroupRef.current) {
        const galaxyDisk = galacticMarkersGroupRef.current.getObjectByName('galaxyDisk')
        if (galaxyDisk) {
          updateGalaxyRotation(galaxyDisk as THREE.Group, deltaTime)
          
          // Also rotate the Solar System / LSR frame around GC if animation is enabled
          const solarSystemDemo = galacticMarkersGroupRef.current.getObjectByName('solarSystemDemo')
          if (solarSystemDemo && (galaxyDisk as THREE.Group).userData?.rotationSpeed) {
            // Get GC position (center of galaxy disk)
            const gcPosition = (galaxyDisk as THREE.Group).position.clone()
            updateLSRRotation(
              solarSystemDemo as THREE.Group,
              gcPosition,
              (galaxyDisk as THREE.Group).userData.rotationSpeed,
              deltaTime
            )
          }
        }
        
        // Update V_circular animation in the observer direction indicator
        const observerIndicator = galacticMarkersGroupRef.current.getObjectByName('observerIndicator')
        if (observerIndicator) {
          updateVCircularAnimation(observerIndicator as THREE.Group, deltaTime)
        }
      }

      // Update controls
      controls.update()

      // Render
      renderer.render(scene, camera)
    }
    animate()

    // Handle resize
    const handleResize = () => {
      const width = container.clientWidth
      const height = container.clientHeight
      camera.aspect = width / height
      camera.updateProjectionMatrix()
      renderer.setSize(width, height)
    }
    window.addEventListener('resize', handleResize)

    // Cleanup
    return () => {
      cancelAnimationFrame(animationId)
      window.removeEventListener('resize', handleResize)
      renderer.dispose()
      geometry.dispose()
      material.dispose()
      container.removeChild(renderer.domElement)
    }
  }, []) // Empty deps - only run once on mount

  // Update material when props change (these are infrequent)
  useEffect(() => {
    if (materialRef.current) {
      materialRef.current.size = pointSize
      materialRef.current.opacity = opacity
      materialRef.current.needsUpdate = true
    }
  }, [pointSize, opacity])
  
  // Update IMBH visualization visibility
  useEffect(() => {
    if (bhGroupRef.current) {
      bhGroupRef.current.visible = showBlackHole
    }
    if (trajectoryLineRef.current) {
      trajectoryLineRef.current.visible = showTrajectory
    }
    if (comMarkerRef.current) {
      comMarkerRef.current.visible = showTrajectory
    }
    if (tidalCircleRef.current) {
      tidalCircleRef.current.visible = showRadii
    }
    if (hillCircleRef.current) {
      hillCircleRef.current.visible = showRadii
    }
    if (galacticMarkersGroupRef.current) {
      galacticMarkersGroupRef.current.visible = showGalacticMarkers
    }

    // Toggle visibility of all text labels (Sprite objects) in the scene
    if (sceneRef.current) {
      sceneRef.current.traverse((object) => {
        if (object instanceof THREE.Sprite) {
          object.visible = showLabels
        }
      })
    }
  }, [showBlackHole, showTrajectory, showRadii, showGalacticMarkers, showLabels])

  // Recreate galactic markers when galacticConfig changes (showGalaxyDisk, showSolarSystem)
  useEffect(() => {
    if (!sceneRef.current || !imbhPhysics?.enabled || !showGalacticMarkers) return
    
    // Remove existing galactic markers
    if (galacticMarkersGroupRef.current) {
      sceneRef.current.remove(galacticMarkersGroupRef.current)
      galacticMarkersGroupRef.current = null
    }
    
    // Recreate with new config
    const bhPos = new THREE.Vector3(
      imbhPhysics.bhPosition[0],
      imbhPhysics.bhPosition[1],
      imbhPhysics.bhPosition[2]
    )
    
    const galConfig: GalacticMarkersConfig = {
      distanceToGC: galacticConfig?.distanceToGC ?? DEFAULT_GALACTIC_CONFIG.distanceToGC,
      galacticLongitude: galacticConfig?.galacticLongitude ?? DEFAULT_GALACTIC_CONFIG.galacticLongitude,
      galacticLatitude: galacticConfig?.galacticLatitude ?? DEFAULT_GALACTIC_CONFIG.galacticLatitude,
      inclination: galacticConfig?.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination,
      positionAngle: galacticConfig?.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle,
      lsrVelocity: galacticConfig?.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity,
      hvccPosition: DEFAULT_GALACTIC_CONFIG.hvccPosition,
      hvccSize: DEFAULT_GALACTIC_CONFIG.hvccSize,
      // Pass cloud radius from imbhPhysics if available
      cloudRadius_pc: imbhPhysics?.cloudRadius ?? DEFAULT_GALACTIC_CONFIG.cloudRadius_pc,
      showGalaxyDisk: galacticConfig?.showGalaxyDisk ?? false,
      showSolarSystem: galacticConfig?.showSolarSystem ?? false,
      galaxyRotationSpeed: galacticConfig?.galaxyRotationSpeed ?? 0.1,
    }

    const galacticMarkers = createAllGalacticMarkers(bhPos, galConfig)
    galacticMarkers.visible = showGalacticMarkers
    sceneRef.current.add(galacticMarkers)
    galacticMarkersGroupRef.current = galacticMarkers
  }, [galacticConfig?.showGalaxyDisk, galacticConfig?.showSolarSystem, galacticConfig?.galaxyRotationSpeed, imbhPhysics, showGalacticMarkers])

  // Recompute stats when frames change
  useEffect(() => {
    computeStats()
  }, [framesRef.current?.size])

  // Handle camera mode changes (free orbit vs Earth view)
  useEffect(() => {
    if (!cameraRef.current || !controlsRef.current || !imbhPhysics?.enabled) return
    
    const camera = cameraRef.current
    const controls = controlsRef.current
    
    if (cameraMode === 'earth') {
      // Position camera at Earth's location, looking toward BH
      // Use same geometry as createAllGalacticMarkers
      const inclination = galacticConfig?.inclination ?? 70
      const positionAngle = galacticConfig?.positionAngle ?? 41.6
      
      const inc = inclination * Math.PI / 180
      const pa = positionAngle * Math.PI / 180
      
      // Line of sight direction (from BH toward Earth)
      const losDir = new THREE.Vector3(
        Math.sin(inc) * Math.cos(pa),
        Math.sin(inc) * Math.sin(pa),
        Math.cos(inc)
      ).normalize()
      
      // Position camera at Earth's symbolic location (along LoS direction)
      const earthDisplayDist = 50  // Symbolic distance for camera
      const bhPos = new THREE.Vector3(
        imbhPhysics.bhPosition[0],
        imbhPhysics.bhPosition[1],
        imbhPhysics.bhPosition[2]
      )
      
      const earthPos = bhPos.clone().add(losDir.clone().multiplyScalar(earthDisplayDist))
      
      // Position camera at Earth, looking at BH
      camera.position.copy(earthPos)
      controls.target.copy(bhPos)
      camera.lookAt(bhPos)
      
      // Reduce controls damping for smoother viewing
      controls.enableDamping = true
      controls.dampingFactor = 0.02
      
      console.log(`[Camera] Earth view: i=${inclination}°, PA=${positionAngle}°, pos=(${earthPos.x.toFixed(1)}, ${earthPos.y.toFixed(1)}, ${earthPos.z.toFixed(1)})`)
    } else {
      // Free mode - reset to default orbit controls
      controls.enableDamping = true
      controls.dampingFactor = 0.05
      
      // Don't force camera position - let user control it
      console.log('[Camera] Free orbit mode')
    }
    
    controls.update()
  }, [cameraMode, galacticConfig?.inclination, galacticConfig?.positionAngle, imbhPhysics])

  return (
    <div ref={containerRef} className={`w-full h-full ${className}`} />
  )
}

export default ParticleViewer3DImperative