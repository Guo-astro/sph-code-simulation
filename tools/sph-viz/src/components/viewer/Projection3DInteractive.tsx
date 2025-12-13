'use client'

import { useRef, useEffect, useCallback, useState, forwardRef, useImperativeHandle } from 'react'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'
import type { ParsedFrame, ColorMap } from '~/types/sph'
import { COLOR_MAPS } from '~/types/sph'
import { getCircleTexture } from '~/utils/three-helpers'
import { sampleColorMap } from '~/utils/color-interpolation'

export interface Projection3DInteractiveProps {
  /** Map of all loaded frames - passed as ref to avoid React re-renders */
  framesRef: React.RefObject<Map<number, ParsedFrame>>
  /** Current frame index ref - read imperatively during animation */
  frameIndexRef: React.RefObject<number>
  /** Initial projection view: xy (top-down), xz (front), yz (side) */
  projection: 'xy' | 'xz' | 'yz'
  colorField?: string
  colorMap?: ColorMap
  width?: number
  height?: number
  className?: string
  /** Show shock tracer sampling regions */
  showShockSampling?: boolean
  /** Shock sampling parameters */
  shockSamplingParams?: {
    columnRadius: number
    sliceThickness: number
  }
  /** Global color range for consistent coloring across views */
  globalColorRange?: [number, number]
  /** Log scale for color mapping */
  logScale?: boolean
  /** Particle size multiplier */
  particleSize?: number
}

export interface Projection3DInteractiveHandle {
  resetCamera: () => void
}

/** Default camera positions for each projection */
const CAMERA_POSITIONS: Record<string, { position: THREE.Vector3; up: THREE.Vector3 }> = {
  xy: { position: new THREE.Vector3(0, 0, 50), up: new THREE.Vector3(0, 1, 0) },  // Looking down Z axis at XY plane
  xz: { position: new THREE.Vector3(0, -50, 0), up: new THREE.Vector3(0, 0, 1) }, // Looking from -Y at XZ plane
  yz: { position: new THREE.Vector3(50, 0, 0), up: new THREE.Vector3(0, 0, 1) },  // Looking from +X at YZ plane
}

/** Default color map */
const defaultColorMap: ColorMap = {
  name: 'Cosmic Dawn',
  colors: COLOR_MAPS.cosmicDawn?.colors || [
    '#1a0533', '#2d1b69', '#3d4fc7', '#00b4d8',
    '#48cae4', '#90e0ef', '#ffd166', '#ffeb99', '#ffffff',
  ],
  min: 0,
  max: 1,
  logScale: false,
}

/**
 * Interactive 3D Projection Viewer
 * 
 * Uses Three.js with OrbitControls for interactive viewing of particle data.
 * Starts with a camera position that mimics a 2D projection but allows
 * full 3D interaction (zoom, pan, rotate).
 */
export const Projection3DInteractive = forwardRef<Projection3DInteractiveHandle, Projection3DInteractiveProps>(({
  framesRef,
  frameIndexRef,
  projection,
  colorField = 'density',
  colorMap = defaultColorMap,
  width = 400,
  height = 400,
  className = '',
  showShockSampling = false,
  shockSamplingParams = { columnRadius: 0.15, sliceThickness: 0.15 },
  globalColorRange,
  logScale = false,
  particleSize = 1.5,
}, ref) => {
  const containerRef = useRef<HTMLDivElement>(null)
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null)
  const sceneRef = useRef<THREE.Scene | null>(null)
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null)
  const controlsRef = useRef<OrbitControls | null>(null)
  const particlesRef = useRef<THREE.Points | null>(null)
  const geometryRef = useRef<THREE.BufferGeometry | null>(null)
  const materialRef = useRef<THREE.PointsMaterial | null>(null)
  
  // Shock sampling visualization objects
  const shockGroupRef = useRef<THREE.Group | null>(null)
  
  // Track the last rendered frame
  const lastFrameIndexRef = useRef<number>(-1)
  const lastColorFieldRef = useRef<string>(colorField)
  const lastLogScaleRef = useRef<boolean>(logScale)
  
  // Animation frame ref
  const animationFrameRef = useRef<number | null>(null)
  
  // Store global color range in ref
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

  // Default camera position for this projection
  const defaultCameraConfig = CAMERA_POSITIONS[projection] || CAMERA_POSITIONS.xy
  
  // State for UI
  const [isRotated, setIsRotated] = useState(false)

  // Reset camera to default position for this projection
  const resetCamera = useCallback(() => {
    const camera = cameraRef.current
    const controls = controlsRef.current
    if (!camera || !controls) return
    
    const config = CAMERA_POSITIONS[projection] || CAMERA_POSITIONS.xy
    camera.position.copy(config.position)
    camera.up.copy(config.up)
    camera.lookAt(0, 0, 0)
    controls.target.set(0, 0, 0)
    controls.update()
    setIsRotated(false)
  }, [projection])

  // Expose resetCamera via ref
  useImperativeHandle(ref, () => ({
    resetCamera
  }), [resetCamera])

  // Compute stats from frames
  const computeStats = useCallback(() => {
    const frames = framesRef.current
    if (!frames || frames.size === 0) return

    const allDensity: number[] = []
    const allVelocity: number[] = []
    const allPressure: number[] = []
    const allEnergy: number[] = []

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

  // Compute Center of Mass
  const computeCoM = useCallback((frame: ParsedFrame): THREE.Vector3 => {
    let totalMass = 0
    let cx = 0, cy = 0, cz = 0
    for (let i = 0; i < frame.particleCount; i++) {
      const m = frame.mass?.[i] ?? 1
      cx += m * frame.positions[i * 3]
      cy += m * frame.positions[i * 3 + 1]
      cz += m * frame.positions[i * 3 + 2]
      totalMass += m
    }
    if (totalMass > 0) {
      cx /= totalMass
      cy /= totalMass
      cz /= totalMass
    }
    return new THREE.Vector3(cx, cy, cz)
  }, [])

  // Update shock sampling visualization
  const updateShockSampling = useCallback((frame: ParsedFrame) => {
    if (!shockGroupRef.current || !showShockSampling) return
    
    // Clear existing objects
    while (shockGroupRef.current.children.length > 0) {
      const child = shockGroupRef.current.children[0]
      shockGroupRef.current.remove(child)
      if (child instanceof THREE.Mesh || child instanceof THREE.Line) {
        child.geometry?.dispose()
        if (child.material instanceof THREE.Material) {
          child.material.dispose()
        }
      }
    }
    
    const com = computeCoM(frame)
    const { columnRadius, sliceThickness } = shockSamplingParams
    
    // CoM marker (white sphere)
    const comGeom = new THREE.SphereGeometry(0.08, 16, 16)
    const comMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
    const comMesh = new THREE.Mesh(comGeom, comMat)
    comMesh.position.copy(com)
    shockGroupRef.current.add(comMesh)
    
    // Z-profile cylinder (red, dashed)
    const cylinderGeom = new THREE.CylinderGeometry(columnRadius, columnRadius, 20, 32, 1, true)
    const cylinderMat = new THREE.MeshBasicMaterial({
      color: 0xff6b6b,
      transparent: true,
      opacity: 0.3,
      side: THREE.DoubleSide,
      wireframe: true,
    })
    const cylinder = new THREE.Mesh(cylinderGeom, cylinderMat)
    cylinder.position.set(com.x, com.y, 0)
    cylinder.rotation.x = Math.PI / 2  // Align along Z axis
    shockGroupRef.current.add(cylinder)
    
    // Z-profile circle outline at CoM z
    const circleGeom = new THREE.RingGeometry(columnRadius - 0.01, columnRadius + 0.01, 64)
    const circleMat = new THREE.MeshBasicMaterial({ color: 0xff6b6b, side: THREE.DoubleSide })
    const circle = new THREE.Mesh(circleGeom, circleMat)
    circle.position.copy(com)
    shockGroupRef.current.add(circle)
    
    // X-profile slice (teal box)
    const sliceGeom = new THREE.BoxGeometry(40, sliceThickness * 2, sliceThickness * 2)
    const sliceMat = new THREE.MeshBasicMaterial({
      color: 0x4ecdc4,
      transparent: true,
      opacity: 0.2,
      side: THREE.DoubleSide,
    })
    const slice = new THREE.Mesh(sliceGeom, sliceMat)
    slice.position.copy(com)
    shockGroupRef.current.add(slice)
    
    // Slice outline
    const sliceEdges = new THREE.EdgesGeometry(sliceGeom)
    const sliceLineMat = new THREE.LineBasicMaterial({ color: 0x4ecdc4 })
    const sliceLine = new THREE.LineSegments(sliceEdges, sliceLineMat)
    sliceLine.position.copy(com)
    shockGroupRef.current.add(sliceLine)
  }, [computeCoM, showShockSampling, shockSamplingParams])

  // Update particle buffers
  const updateParticles = useCallback(() => {
    const frameIndex = frameIndexRef.current ?? 0
    const frames = framesRef.current
    if (!frames || !geometryRef.current) return

    const frame = frames.get(frameIndex)
    if (!frame) return

    const needsUpdate = 
      frameIndex !== lastFrameIndexRef.current ||
      colorField !== lastColorFieldRef.current ||
      logScale !== lastLogScaleRef.current

    if (!needsUpdate) return

    lastFrameIndexRef.current = frameIndex
    lastColorFieldRef.current = colorField
    lastLogScaleRef.current = logScale

    const geometry = geometryRef.current
    const posAttr = geometry.attributes.position as THREE.BufferAttribute
    const colorAttr = geometry.attributes.color as THREE.BufferAttribute

    // Resize buffers if needed
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
      default: // density
        fieldData = frame.density
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current!
        } else {
          [vMin, vMax] = statsRef.current.density
        }
    }

    // Apply colors
    const colors = (geometry.attributes.color as THREE.BufferAttribute).array as Float32Array
    // Map color map name from colorMap.name to COLOR_MAP_DATA key
    const colorMapKey = colorMap.name || 'Cosmic Dawn'

    for (let i = 0; i < frame.particleCount; i++) {
      let value = fieldData[i]
      
      // Apply log scale if enabled
      if (logScale && value > 0 && vMin > 0) {
        const logMin = Math.log10(vMin)
        const logMax = Math.log10(vMax)
        value = (Math.log10(value) - logMin) / (logMax - logMin)
      } else {
        value = (value - vMin) / (vMax - vMin)
      }
      
      value = Math.max(0, Math.min(1, value))
      const rgb = sampleColorMap(colorMapKey, value)
      colors[i * 3] = rgb[0]
      colors[i * 3 + 1] = rgb[1]
      colors[i * 3 + 2] = rgb[2]
    }
    ;(geometry.attributes.color as THREE.BufferAttribute).needsUpdate = true

    // Update shock sampling visualization
    if (showShockSampling) {
      updateShockSampling(frame)
    }
  }, [framesRef, frameIndexRef, colorField, colorMap, logScale, showShockSampling, updateShockSampling])

  // Initialize Three.js scene
  useEffect(() => {
    if (!containerRef.current) return

    // Create renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true })
    renderer.setSize(width, height)
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    renderer.setClearColor(0x0a0a14, 1)
    containerRef.current.appendChild(renderer.domElement)
    rendererRef.current = renderer

    // Create scene
    const scene = new THREE.Scene()
    sceneRef.current = scene

    // Create camera
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000)
    const config = CAMERA_POSITIONS[projection] || CAMERA_POSITIONS.xy
    camera.position.copy(config.position)
    camera.up.copy(config.up)
    camera.lookAt(0, 0, 0)
    cameraRef.current = camera

    // Create OrbitControls
    const controls = new OrbitControls(camera, renderer.domElement)
    controls.enableDamping = true
    controls.dampingFactor = 0.1
    controls.rotateSpeed = 0.5
    controls.zoomSpeed = 1.2
    controls.panSpeed = 0.8
    controls.minDistance = 5
    controls.maxDistance = 200
    controls.target.set(0, 0, 0)
    controlsRef.current = controls
    
    // Track when user rotates the view
    controls.addEventListener('change', () => {
      const currentPos = camera.position.clone().normalize()
      const defaultPos = config.position.clone().normalize()
      const dotProduct = currentPos.dot(defaultPos)
      setIsRotated(dotProduct < 0.99)
    })

    // Create particle geometry and material
    const geometry = new THREE.BufferGeometry()
    const initialPositions = new Float32Array(100000 * 3)
    const initialColors = new Float32Array(100000 * 3)
    geometry.setAttribute('position', new THREE.BufferAttribute(initialPositions, 3).setUsage(THREE.DynamicDrawUsage))
    geometry.setAttribute('color', new THREE.BufferAttribute(initialColors, 3).setUsage(THREE.DynamicDrawUsage))
    geometryRef.current = geometry

    const texture = getCircleTexture()
    const material = new THREE.PointsMaterial({
      size: particleSize * 0.08,
      vertexColors: true,
      transparent: true,
      opacity: 0.9,
      sizeAttenuation: true,
      map: texture,
      alphaTest: 0.01,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    })
    materialRef.current = material

    const particles = new THREE.Points(geometry, material)
    scene.add(particles)
    particlesRef.current = particles

    // Add simple grid for reference
    const gridHelper = new THREE.GridHelper(40, 20, 0x444444, 0x222222)
    if (projection === 'xy') {
      gridHelper.rotation.x = Math.PI / 2  // XY plane
    } else if (projection === 'xz') {
      // Default orientation is XZ
    } else if (projection === 'yz') {
      gridHelper.rotation.z = Math.PI / 2  // YZ plane
    }
    scene.add(gridHelper)

    // Add axes helper
    const axesHelper = new THREE.AxesHelper(10)
    scene.add(axesHelper)
    
    // Create shock sampling group
    const shockGroup = new THREE.Group()
    scene.add(shockGroup)
    shockGroupRef.current = shockGroup

    // Compute initial stats
    computeStats()

    // Animation loop
    const animate = () => {
      animationFrameRef.current = requestAnimationFrame(animate)
      
      updateParticles()
      controls.update()
      renderer.render(scene, camera)
    }
    animate()

    // Cleanup
    return () => {
      if (animationFrameRef.current) {
        cancelAnimationFrame(animationFrameRef.current)
      }
      controls.dispose()
      renderer.dispose()
      geometry.dispose()
      material.dispose()
      if (containerRef.current && renderer.domElement) {
        containerRef.current.removeChild(renderer.domElement)
      }
    }
  }, [projection, computeStats, updateParticles, width, height, particleSize])

  // Update renderer size when dimensions change
  useEffect(() => {
    const renderer = rendererRef.current
    const camera = cameraRef.current
    if (renderer && camera) {
      renderer.setSize(width, height)
      camera.aspect = width / height
      camera.updateProjectionMatrix()
    }
  }, [width, height])

  // Update material when particle size changes
  useEffect(() => {
    const material = materialRef.current
    if (material) {
      material.size = particleSize * 0.08
    }
  }, [particleSize])

  return (
    <div className={`relative ${className}`}>
      <div
        ref={containerRef}
        style={{ width, height }}
        className="rounded border border-gray-700 overflow-hidden"
      />
      {/* Projection label */}
      <div className="absolute top-1 left-1 text-cyan-300 text-sm font-bold bg-black/70 px-2 py-0.5 rounded">
        {projection.toUpperCase()}
      </div>
      {/* Controls overlay */}
      <div className="absolute top-1 right-1 flex gap-1">
        {/* Reset camera button */}
        <button
          onClick={resetCamera}
          className={`
            px-1.5 py-0.5 text-[10px] rounded transition-all
            ${isRotated 
              ? 'bg-cyan-600 hover:bg-cyan-500 text-white' 
              : 'bg-gray-700/80 text-gray-400 hover:bg-gray-600'
            }
          `}
          title="Reset camera to default view"
        >
          ⟲ Reset
        </button>
      </div>
      {/* Interaction hint */}
      <div className="absolute bottom-1 left-1 text-gray-500 text-[9px] bg-black/50 px-1 rounded">
        Drag: rotate | Scroll: zoom | Shift+drag: pan
      </div>
      {/* Legend for shock sampling */}
      {showShockSampling && (
        <div className="absolute bottom-1 right-1 text-[9px] bg-black/70 px-1 rounded">
          <span className="text-red-400">● Z-col</span>
          <span className="mx-1 text-gray-500">|</span>
          <span className="text-teal-400">■ X-slice</span>
        </div>
      )}
    </div>
  )
})

Projection3DInteractive.displayName = 'Projection3DInteractive'

export default Projection3DInteractive
