'use client'

import { useRef, useEffect, useCallback, useState } from 'react'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'
import type { ParsedFrame, ColorMap } from '~/types/sph'

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
 * GEOMETRY SETUP:
 * ===============
 * - BH is at origin (0, 0, 0)
 * - Cloud approaches from -x direction with velocity (vx, 0, 0)
 * - Impact parameter b = y-coordinate of initial position (perpendicular distance to asymptote)
 * - Orbit lies in the x-y plane
 * 
 * HYPERBOLIC ORBIT EQUATIONS:
 * ===========================
 * Semi-major axis:     a = GM / v∞²
 * Eccentricity:        e = √(1 + (b v∞² / GM)²)
 * Semi-latus rectum:   p = a(e² - 1) = b² v∞² / GM
 * Perihelion:          r_p = a(e - 1)
 * 
 * Polar equation:      r = p / (1 + e cos(θ - θ_peri))
 * 
 * where θ_peri is the angle to perihelion from BH.
 * 
 * For approach from -x with y > 0:
 *   - Asymptotic direction: along -x
 *   - Perihelion is at angle θ_peri = π/2 - δ/2 where δ = deflection angle
 *   - Deflection: tan(δ/2) = GM / (b v∞²) = 1/e
 */
function computeHyperbolicOrbit(
  config: {
    bhPosition: [number, number, number]
    bhMass: number  // in code units (1 = 1000 M_sun)
    cloudInitialPosition: [number, number, number]
    cloudInitialVelocity: [number, number, number]
  },
  numPoints: number = 200
): [number, number, number][] {
  const { bhPosition, bhMass, cloudInitialPosition, cloudInitialVelocity } = config
  
  const [x0, y0, z0] = cloudInitialPosition
  const [vx0, vy0, vz0] = cloudInitialVelocity
  const [bh_x, bh_y, bh_z] = bhPosition
  
  // ========== Orbital Parameters ==========
  // Velocity at infinity (approach velocity)
  const v_inf = Math.sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0)
  
  // Impact parameter (perpendicular distance from BH to asymptotic velocity vector)
  // For motion in +x direction starting at (x0, y0, z0), b = |y0 - bh_y|
  const b = Math.abs(y0 - bh_y)
  
  // Gravitational parameter
  const GM = UNITS.G * bhMass
  
  // Semi-major axis (positive for hyperbola calculation)
  const a = GM / (v_inf * v_inf)
  
  // Eccentricity (e > 1 for hyperbolic orbit)
  const e = Math.sqrt(1 + Math.pow(b * v_inf * v_inf / GM, 2))
  
  // Semi-latus rectum
  const p = a * (e * e - 1)  // = b² v∞² / GM
  
  // Perihelion distance (closest approach)
  const r_peri = a * (e - 1)
  
  // ========== Orbit Orientation ==========
  // Half deflection angle
  const half_deflection = Math.atan(1 / e)  // = atan(GM / (b v∞²))
  
  // Full deflection angle
  const deflection = 2 * half_deflection
  
  // Asymptotic true anomaly (angle where r → ∞)
  // cos(θ_asymp) = -1/e, so θ_asymp = π - acos(1/e)
  const theta_asymp = Math.PI - Math.acos(1 / e)
  
  // Perihelion angle in our coordinate system:
  // Cloud approaches from -x (θ = π), so perihelion is at θ_peri
  // For y0 > 0 (cloud above x-axis): perihelion is in upper half
  const sign_y = y0 >= bh_y ? 1 : -1
  const theta_peri = Math.PI / 2 * sign_y  // Perihelion at +y (or -y if cloud below)
  
  // ========== Generate Orbit Points ==========
  // True anomaly ranges from incoming asymptote to outgoing asymptote
  // ν ranges from -(π - acos(1/e)) to +(π - acos(1/e))
  const nu_max = Math.PI - Math.acos(1 / e) - 0.02  // Slightly less than asymptotic
  const nu_min = -nu_max
  
  const orbit: [number, number, number][] = []
  
  for (let i = 0; i < numPoints; i++) {
    // True anomaly ν (measured from perihelion)
    const nu = nu_min + (nu_max - nu_min) * (i / (numPoints - 1))
    
    // Radial distance
    const r = p / (1 + e * Math.cos(nu))
    
    // Angle in x-y plane (θ = θ_peri + ν)
    // But we need to rotate so that incoming asymptote aligns with -x direction
    // 
    // At ν = -nu_max (incoming), the position should be at large -x, y ≈ b
    // At ν = 0 (perihelion), position is at (0, r_peri) for y > 0 case
    // At ν = +nu_max (outgoing), position is at large x with deflected y
    //
    // For incoming from -x with impact parameter b (y-offset):
    // We rotate the standard orbit by angle (π - theta_peri + nu_max)
    const rotation = Math.PI - theta_asymp * sign_y
    const theta = nu + rotation
    
    // Cartesian coordinates
    const x = bh_x + r * Math.cos(theta)
    const y = bh_y + r * Math.sin(theta) * sign_y
    const z = bh_z + z0 * (1 - Math.abs(nu) / nu_max * 0.1)  // Small z variation
    
    // Only include points within reasonable bounds and on the correct side of approach
    if (Math.abs(x) < 100 && Math.abs(y) < 100) {
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

/** Create a circular particle texture for smooth round particles */
function createCircleTexture(): THREE.Texture {
  const size = 64
  const canvas = document.createElement('canvas')
  canvas.width = size
  canvas.height = size
  const ctx = canvas.getContext('2d')!
  
  // Create radial gradient for smooth circular particle
  const gradient = ctx.createRadialGradient(
    size / 2, size / 2, 0,
    size / 2, size / 2, size / 2
  )
  gradient.addColorStop(0, 'rgba(255, 255, 255, 1)')
  gradient.addColorStop(0.3, 'rgba(255, 255, 255, 0.8)')
  gradient.addColorStop(0.7, 'rgba(255, 255, 255, 0.3)')
  gradient.addColorStop(1, 'rgba(255, 255, 255, 0)')
  
  ctx.fillStyle = gradient
  ctx.fillRect(0, 0, size, size)
  
  const texture = new THREE.CanvasTexture(canvas)
  texture.needsUpdate = true
  return texture
}

// Singleton texture for all particles
let circleTexture: THREE.Texture | null = null
function getCircleTexture(): THREE.Texture {
  if (!circleTexture) {
    circleTexture = createCircleTexture()
  }
  return circleTexture
}

// ============== Color Maps ==============
const COLOR_MAP_DATA: Record<string, number[][]> = {
  viridis: [[0.267,0.004,0.329],[0.282,0.140,0.458],[0.253,0.265,0.530],[0.206,0.372,0.553],[0.163,0.471,0.558],[0.127,0.566,0.551],[0.134,0.658,0.518],[0.266,0.749,0.441],[0.477,0.821,0.318],[0.741,0.873,0.150],[0.993,0.906,0.144]],
  plasma: [[0.050,0.030,0.528],[0.295,0.012,0.615],[0.492,0.012,0.658],[0.665,0.138,0.618],[0.798,0.280,0.470],[0.899,0.396,0.301],[0.966,0.530,0.128],[0.988,0.680,0.063],[0.961,0.850,0.298],[0.940,0.975,0.131]],
  inferno: [[0.001,0.000,0.014],[0.122,0.047,0.282],[0.304,0.063,0.420],[0.499,0.086,0.397],[0.680,0.144,0.295],[0.833,0.253,0.160],[0.937,0.405,0.049],[0.981,0.588,0.068],[0.987,0.772,0.264],[0.988,0.998,0.645]],
  turbo: [[0.190,0.072,0.232],[0.235,0.318,0.860],[0.137,0.572,0.996],[0.140,0.780,0.820],[0.376,0.920,0.512],[0.670,0.979,0.280],[0.924,0.904,0.145],[0.996,0.724,0.132],[0.994,0.472,0.122],[0.881,0.200,0.102],[0.528,0.055,0.052]]
}

function sampleColorMap(name: string, t: number): [number, number, number] {
  const map = COLOR_MAP_DATA[name] || COLOR_MAP_DATA.viridis
  t = Math.max(0, Math.min(1, t))
  const idx = t * (map.length - 1)
  const i = Math.floor(idx)
  const f = idx - i
  if (i >= map.length - 1) return map[map.length - 1] as [number, number, number]
  return [
    map[i][0] + f * (map[i+1][0] - map[i][0]),
    map[i][1] + f * (map[i+1][1] - map[i][1]),
    map[i][2] + f * (map[i+1][2] - map[i][2])
  ]
}

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

    // Create text sprite for axis labels
    const createTextSprite = (text: string, color: string = '#ffffff'): THREE.Sprite => {
      const canvas = document.createElement('canvas')
      const size = 256
      canvas.width = size
      canvas.height = size / 2
      const ctx = canvas.getContext('2d')!
      ctx.fillStyle = color
      ctx.font = 'Bold 48px Arial'
      ctx.textAlign = 'center'
      ctx.textBaseline = 'middle'
      ctx.fillText(text, size / 2, size / 4)
      
      const texture = new THREE.CanvasTexture(canvas)
      const material = new THREE.SpriteMaterial({ map: texture, transparent: true })
      const sprite = new THREE.Sprite(material)
      sprite.scale.set(8, 4, 1)
      return sprite
    }

    // Axes helper with physical unit labels
    if (showAxes) {
      const axesHelper = new THREE.AxesHelper(30)
      scene.add(axesHelper)
      
      // Add axis labels with physical units
      const xLabel = createTextSprite('X (pc)', '#ff4444')
      xLabel.position.set(32, 0, 0)
      scene.add(xLabel)
      
      const yLabel = createTextSprite('Y (pc)', '#44ff44')
      yLabel.position.set(0, 32, 0)
      scene.add(yLabel)
      
      const zLabel = createTextSprite('Z (pc)', '#4444ff')
      zLabel.position.set(0, 0, 32)
      scene.add(zLabel)
      
      // Add tick marks at 10 pc intervals
      const tickMaterial = new THREE.LineBasicMaterial({ color: 0x666666 })
      for (let i = -20; i <= 20; i += 10) {
        if (i === 0) continue
        // X-axis ticks
        const xTickGeom = new THREE.BufferGeometry().setFromPoints([
          new THREE.Vector3(i, -0.5, 0),
          new THREE.Vector3(i, 0.5, 0)
        ])
        scene.add(new THREE.Line(xTickGeom, tickMaterial))
        
        // Y-axis ticks  
        const yTickGeom = new THREE.BufferGeometry().setFromPoints([
          new THREE.Vector3(-0.5, i, 0),
          new THREE.Vector3(0.5, i, 0)
        ])
        scene.add(new THREE.Line(yTickGeom, tickMaterial))
        
        // Add tick labels
        if (i % 10 === 0) {
          const tickLabel = createTextSprite(`${i}`, '#888888')
          tickLabel.position.set(i, -2, 0)
          tickLabel.scale.set(4, 2, 1)
          scene.add(tickLabel)
        }
      }
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
      if (showTrajectory) {
        // Compute hyperbolic orbit using analytical solution
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
    const animate = () => {
      animationId = requestAnimationFrame(animate)

      // FPS calculation
      fpsRef.current.frames++
      const now = performance.now()
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
  }, [showBlackHole, showTrajectory, showRadii])

  // Recompute stats when frames change
  useEffect(() => {
    computeStats()
  }, [framesRef.current?.size])

  return (
    <div ref={containerRef} className={`w-full h-full ${className}`} />
  )
}

export default ParticleViewer3DImperative
