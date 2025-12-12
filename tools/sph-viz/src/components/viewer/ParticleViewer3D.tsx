'use client'

import { useRef, useMemo, useEffect, useState, useCallback } from 'react'
import { Canvas, useThree } from '@react-three/fiber'
import { OrbitControls, PerspectiveCamera, Stats } from '@react-three/drei'
import * as THREE from 'three'
import type { ParsedFrame, ColorMap } from '~/types/sph'
// Shared utilities (SSOT)
import { getCircleTexture } from '~/utils/three-helpers'
import { interpolateColorHex, hexToRgbCached } from '~/utils/color-interpolation'

interface ParticleCloudProps {
  frame: ParsedFrame | null
  colorField: string
  colorMap: ColorMap
  pointSize: number
  opacity: number
}

/** GPU-accelerated particle cloud - optimized to avoid recreating geometry */
function ParticleCloud({ frame, colorField, colorMap, pointSize, opacity }: ParticleCloudProps) {
  const pointsRef = useRef<THREE.Points>(null)
  const geometryRef = useRef<THREE.BufferGeometry | null>(null)
  const materialRef = useRef<THREE.PointsMaterial | null>(null)
  
  // Pre-allocate arrays for reuse (avoid GC pressure)
  const colorsArrayRef = useRef<Float32Array | null>(null)

  // Update positions and colors when frame changes - MUTATE existing buffers
  useEffect(() => {
    if (!frame || !pointsRef.current) return

    const count = frame.particleCount
    let geometry = geometryRef.current

    // Create geometry only once, or if particle count changes
    if (!geometry || (geometry.attributes.position && geometry.attributes.position.count !== count)) {
      geometry = new THREE.BufferGeometry()
      const positions = new THREE.BufferAttribute(frame.positions, 3)
      positions.setUsage(THREE.DynamicDrawUsage)
      geometry.setAttribute('position', positions)
      
      colorsArrayRef.current = new Float32Array(count * 3)
      const colors = new THREE.BufferAttribute(colorsArrayRef.current, 3)
      colors.setUsage(THREE.DynamicDrawUsage)
      geometry.setAttribute('color', colors)
      
      geometryRef.current = geometry
      pointsRef.current.geometry = geometry
    } else {
      // Just update the position buffer data
      const posAttr = geometry.attributes.position as THREE.BufferAttribute
      posAttr.array = frame.positions
      posAttr.needsUpdate = true
    }

    // Get the color field data
    let fieldData: Float32Array | undefined
    switch (colorField) {
      case 'density':
        fieldData = frame.density
        break
      case 'pressure':
        fieldData = frame.pressure
        break
      case 'energy':
        fieldData = frame.energy
        break
      case 'velocity':
        // Compute velocity magnitude - reuse array if possible
        if (!colorsArrayRef.current || colorsArrayRef.current.length !== count * 3) {
          colorsArrayRef.current = new Float32Array(count * 3)
        }
        const velMag = new Float32Array(count)
        for (let i = 0; i < count; i++) {
          const vx = frame.velocities[i * 3]
          const vy = frame.velocities[i * 3 + 1]
          const vz = frame.velocities[i * 3 + 2]
          velMag[i] = Math.sqrt(vx * vx + vy * vy + vz * vz)
        }
        fieldData = velMag
        break
      case 'machNumber':
        fieldData = frame.machNumber
        break
      default:
        fieldData = frame.density
    }

    // Find min/max for normalization
    let min = colorMap.min ?? Infinity
    let max = colorMap.max ?? -Infinity
    
    if (min === Infinity || max === -Infinity) {
      if (fieldData) {
        for (let i = 0; i < count; i++) {
          const val = fieldData[i]
          if (isFinite(val)) {
            if (val < min) min = val
            if (val > max) max = val
          }
        }
      }
    }

    if (min === max) max = min + 1

    // Update colors directly in the buffer
    const colorAttr = geometry.attributes.color as THREE.BufferAttribute
    const colors = colorAttr.array as Float32Array
    
    const logMin = colorMap.logScale && min > 0 ? Math.log10(min) : 0
    const logRange = colorMap.logScale && min > 0 ? Math.log10(max) - logMin : 1
    const range = max - min

    for (let i = 0; i < count; i++) {
      let val = fieldData ? fieldData[i] : 0
      if (!isFinite(val)) val = min

      let t: number
      if (colorMap.logScale && min > 0) {
        t = (Math.log10(val) - logMin) / logRange
      } else {
        t = (val - min) / range
      }
      t = Math.max(0, Math.min(1, t))

      const color = interpolateColorHex(colorMap.colors, t)
      colors[i * 3] = color.r
      colors[i * 3 + 1] = color.g
      colors[i * 3 + 2] = color.b
    }

    colorAttr.needsUpdate = true
    geometry.computeBoundingSphere()
  }, [frame, colorField, colorMap])

  // Update material properties
  useEffect(() => {
    if (materialRef.current) {
      materialRef.current.size = pointSize
      materialRef.current.opacity = opacity
      materialRef.current.needsUpdate = true
    }
  }, [pointSize, opacity])

  if (!frame) return null

  return (
    <points ref={pointsRef}>
      <bufferGeometry ref={geometryRef as any} />
      <pointsMaterial
        ref={materialRef as any}
        size={pointSize}
        vertexColors
        transparent
        opacity={opacity}
        sizeAttenuation
        depthWrite={false}
        blending={THREE.AdditiveBlending}
        map={getCircleTexture()}
        alphaTest={0.01}
      />
    </points>
  )
}

// Note: Color interpolation functions (interpolateColorHex, hexToRgbCached, hexToRgb)
// are now imported from ~/utils/color-interpolation (SSOT)

/** Axes helper */
function AxesHelper({ size = 1 }: { size?: number }) {
  return <axesHelper args={[size]} />
}

/** Grid helper */
function GridHelper({ size = 10, divisions = 10 }: { size?: number; divisions?: number }) {
  return <gridHelper args={[size, divisions, '#444444', '#222222']} rotation={[Math.PI / 2, 0, 0]} />
}

/** Bounding box helper */
function BoundingBox({
  min,
  max,
}: {
  min: [number, number, number]
  max: [number, number, number]
}) {
  const geometry = useMemo(() => {
    const box = new THREE.BoxGeometry(max[0] - min[0], max[1] - min[1], max[2] - min[2])
    box.translate((max[0] + min[0]) / 2, (max[1] + min[1]) / 2, (max[2] + min[2]) / 2)
    return new THREE.EdgesGeometry(box)
  }, [min, max])

  return (
    <lineSegments geometry={geometry}>
      <lineBasicMaterial color="#666666" />
    </lineSegments>
  )
}

/** Camera controller with auto-fit */
function CameraController({
  boundingBox,
  resetKey,
}: {
  boundingBox?: { min: [number, number, number]; max: [number, number, number] }
  resetKey?: number
}) {
  const { camera } = useThree()

  useEffect(() => {
    if (boundingBox) {
      const center = new THREE.Vector3(
        (boundingBox.min[0] + boundingBox.max[0]) / 2,
        (boundingBox.min[1] + boundingBox.max[1]) / 2,
        (boundingBox.min[2] + boundingBox.max[2]) / 2
      )
      const size = new THREE.Vector3(
        boundingBox.max[0] - boundingBox.min[0],
        boundingBox.max[1] - boundingBox.min[1],
        boundingBox.max[2] - boundingBox.min[2]
      )
      const maxDim = Math.max(size.x, size.y, size.z)
      const distance = maxDim * 2

      camera.position.set(center.x + distance, center.y + distance * 0.5, center.z + distance)
      camera.lookAt(center)
    }
  }, [boundingBox, camera, resetKey])

  return null
}

export interface ParticleViewer3DProps {
  frame: ParsedFrame | null
  colorField?: string
  colorMap?: ColorMap
  pointSize?: number
  opacity?: number
  showAxes?: boolean
  showGrid?: boolean
  showBoundingBox?: boolean
  boundingBox?: { min: [number, number, number]; max: [number, number, number] }
  showStats?: boolean
  className?: string
}

/** Default color map */
const defaultColorMap: ColorMap = {
  name: 'Viridis',
  colors: [
    '#440154',
    '#482878',
    '#3e4989',
    '#31688e',
    '#26828e',
    '#1f9e89',
    '#35b779',
    '#6ece58',
    '#b5de2b',
    '#fde725',
  ],
  min: 0,
  max: 1,
  logScale: false,
}

/** Main 3D particle viewer component */
export function ParticleViewer3D({
  frame,
  colorField = 'density',
  colorMap = defaultColorMap,
  pointSize = 0.02,
  opacity = 0.8,
  showAxes = true,
  showGrid = false,
  showBoundingBox = true,
  boundingBox,
  showStats = false,
  className = '',
}: ParticleViewer3DProps) {
  const [resetKey, setResetKey] = useState(0)

  const handleResetCamera = () => {
    setResetKey((k) => k + 1)
  }

  return (
    <div className={`relative w-full h-full bg-gray-900 ${className}`}>
      <Canvas>
        <PerspectiveCamera makeDefault fov={60} near={0.001} far={1000} />
        <OrbitControls enableDamping dampingFactor={0.1} rotateSpeed={0.5} />
        <CameraController boundingBox={boundingBox} resetKey={resetKey} />

        <ambientLight intensity={0.5} />
        <directionalLight position={[10, 10, 5]} intensity={1} />

        {showAxes && <AxesHelper size={boundingBox ? Math.max(...boundingBox.max) * 0.5 : 1} />}
        {showGrid && <GridHelper size={10} divisions={10} />}
        {showBoundingBox && boundingBox && <BoundingBox min={boundingBox.min} max={boundingBox.max} />}

        <ParticleCloud
          frame={frame}
          colorField={colorField}
          colorMap={colorMap}
          pointSize={pointSize}
          opacity={opacity}
        />

        {showStats && <Stats />}
      </Canvas>

      {/* Overlay controls */}
      <div className="absolute top-2 right-2 flex gap-2">
        <button
          onClick={handleResetCamera}
          className="px-3 py-1 bg-gray-700 text-white text-sm rounded hover:bg-gray-600"
        >
          Reset Camera
        </button>
      </div>

      {/* Info overlay */}
      {frame && (
        <div className="absolute bottom-2 left-2 text-white text-xs bg-black/50 px-2 py-1 rounded">
          Frame {frame.frameIndex} | t = {frame.time.toFixed(4)} | {frame.particleCount.toLocaleString()} particles
        </div>
      )}
    </div>
  )
}

export default ParticleViewer3D
