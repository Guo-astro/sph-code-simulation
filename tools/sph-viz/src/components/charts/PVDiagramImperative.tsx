'use client'

import { useRef, useEffect, useCallback, useState } from 'react'
import type { ParsedFrame } from '~/types/sph'

// ============================================================================
// OBSERVER GEOMETRY FOR P-V DIAGRAM
// ============================================================================

interface ObserverGeometry {
  /** Orbital plane inclination (degrees) - 0 = face-on, 90 = edge-on */
  inclination: number
  /** Position angle on sky (degrees) */
  positionAngle: number
  /** Systemic velocity V_LSR (km/s) - negative = approaching */
  vLSR: number
  /** Slit position angle (degrees) - direction along which to slice */
  slitPA: number
  /** Slit offset from center (pc) */
  slitOffset: number
  /** Slit width (pc) */
  slitWidth: number
}

const DEFAULT_OBSERVER: ObserverGeometry = {
  inclination: 70,
  positionAngle: 41.6,
  vLSR: -120,
  slitPA: 0,
  slitOffset: 0,
  slitWidth: 2.0,
}

// ============================================================================
// P-V DIAGRAM COMPUTATION
// ============================================================================

interface PVData {
  positionBins: number[]
  velocityBins: number[]
  histogram: Float32Array // 2D histogram flattened
  numPosBins: number
  numVelBins: number
  maxCount: number
}

/**
 * Compute rotation matrix for observer transformation
 */
function computeRotationMatrix(inclination: number, positionAngle: number): number[][] {
  const i = inclination * Math.PI / 180
  const pa = positionAngle * Math.PI / 180

  // Rotation about X by inclination, then about Z by position angle
  const ci = Math.cos(i), si = Math.sin(i)
  const cp = Math.cos(pa), sp = Math.sin(pa)

  // Combined rotation: R_PA * R_inc
  return [
    [cp, -sp * ci, sp * si],
    [sp, cp * ci, -cp * si],
    [0, si, ci],
  ]
}

/**
 * Transform position and velocity to observer frame
 */
function transformToObserver(
  pos: [number, number, number],
  vel: [number, number, number],
  R: number[][]
): { posObs: [number, number, number]; velObs: [number, number, number] } {
  const posObs: [number, number, number] = [
    R[0][0] * pos[0] + R[0][1] * pos[1] + R[0][2] * pos[2],
    R[1][0] * pos[0] + R[1][1] * pos[1] + R[1][2] * pos[2],
    R[2][0] * pos[0] + R[2][1] * pos[1] + R[2][2] * pos[2],
  ]
  const velObs: [number, number, number] = [
    R[0][0] * vel[0] + R[0][1] * vel[1] + R[0][2] * vel[2],
    R[1][0] * vel[0] + R[1][1] * vel[1] + R[1][2] * vel[2],
    R[2][0] * vel[0] + R[2][1] * vel[1] + R[2][2] * vel[2],
  ]
  return { posObs, velObs }
}

/**
 * Compute P-V diagram from particle data
 */
function computePVDiagram(
  frame: ParsedFrame,
  observer: ObserverGeometry,
  posRange: [number, number] = [-5, 5],
  velRange: [number, number] = [-150, 0],
  numPosBins: number = 80,
  numVelBins: number = 80
): PVData {
  const R = computeRotationMatrix(observer.inclination, observer.positionAngle)
  const slitAngle = observer.slitPA * Math.PI / 180
  const cosSlit = Math.cos(slitAngle)
  const sinSlit = Math.sin(slitAngle)

  const histogram = new Float32Array(numPosBins * numVelBins)
  let maxCount = 0

  const posMin = posRange[0]
  const posMax = posRange[1]
  const velMin = velRange[0]
  const velMax = velRange[1]
  const posStep = (posMax - posMin) / numPosBins
  const velStep = (velMax - velMin) / numVelBins

  for (let i = 0; i < frame.particleCount; i++) {
    const x = frame.positions[i * 3]
    const y = frame.positions[i * 3 + 1]
    const z = frame.positions[i * 3 + 2]
    const vx = frame.velocities[i * 3]
    const vy = frame.velocities[i * 3 + 1]
    const vz = frame.velocities[i * 3 + 2]
    const density = frame.density[i]

    if (!isFinite(density) || density <= 0) continue

    // Transform to observer frame
    const { posObs, velObs } = transformToObserver([x, y, z], [vx, vy, vz], R)

    // Position along slit (projection onto slit direction)
    const posSlit = posObs[0] * cosSlit + posObs[1] * sinSlit

    // Check if particle is within slit width
    const perpDist = Math.abs(-posObs[0] * sinSlit + posObs[1] * cosSlit - observer.slitOffset)
    if (perpDist > observer.slitWidth / 2) continue

    // Line-of-sight velocity (z-component in observer frame + V_LSR)
    const vLOS = velObs[2] + observer.vLSR

    // Bin indices
    const posBin = Math.floor((posSlit - posMin) / posStep)
    const velBin = Math.floor((vLOS - velMin) / velStep)

    if (posBin >= 0 && posBin < numPosBins && velBin >= 0 && velBin < numVelBins) {
      // Weight by density (or density^2 for optically thin emission proxy)
      const weight = Math.sqrt(density) // Compromise between density and density^2
      histogram[velBin * numPosBins + posBin] += weight
      maxCount = Math.max(maxCount, histogram[velBin * numPosBins + posBin])
    }
  }

  // Create bin centers
  const positionBins = Array.from({ length: numPosBins }, (_, i) => posMin + (i + 0.5) * posStep)
  const velocityBins = Array.from({ length: numVelBins }, (_, i) => velMin + (i + 0.5) * velStep)

  return { positionBins, velocityBins, histogram, numPosBins, numVelBins, maxCount }
}

// ============================================================================
// COMPONENT
// ============================================================================

export interface PVDiagramImperativeProps {
  /** Ref to frames map - access imperatively */
  framesRef: React.RefObject<Map<number, ParsedFrame>>
  /** Ref to current frame index */
  frameIndexRef: React.RefObject<number>
  /** Initial observer geometry */
  initialObserver?: Partial<ObserverGeometry>
  /** Position range in pc */
  posRange?: [number, number]
  /** Velocity range in km/s */
  velRange?: [number, number]
  /** Number of position bins */
  numPosBins?: number
  /** Number of velocity bins */
  numVelBins?: number
  className?: string
  /** Callback when observer geometry changes */
  onObserverChange?: (observer: ObserverGeometry) => void
}

/**
 * Interactive Position-Velocity Diagram
 *
 * Uses canvas for fast rendering with draggable observer controls.
 * Computes p-v diagram from SPH particle data in real-time.
 */
export function PVDiagramImperative({
  framesRef,
  frameIndexRef,
  initialObserver = {},
  posRange = [-5, 5],
  velRange = [-150, 0],
  numPosBins = 80,
  numVelBins = 80,
  className = '',
  onObserverChange,
}: PVDiagramImperativeProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const observerRef = useRef<ObserverGeometry>({ ...DEFAULT_OBSERVER, ...initialObserver })
  const pvDataRef = useRef<PVData | null>(null)
  const lastFrameIndexRef = useRef<number>(-1)
  const lastUpdateTimeRef = useRef<number>(0)
  const animationIdRef = useRef<number>(0)
  const isDraggingRef = useRef<boolean>(false)
  const dragStartRef = useRef<{ x: number; y: number; inc: number; pa: number } | null>(null)

  // Store prop refs internally to avoid stale closures
  const internalFramesRef = useRef(framesRef)
  const internalFrameIndexRef = useRef(frameIndexRef)

  // Keep internal refs synced with props
  useEffect(() => {
    internalFramesRef.current = framesRef
    internalFrameIndexRef.current = frameIndexRef
  }, [framesRef, frameIndexRef])

  // Expose observer state for UI
  const [observer, setObserver] = useState<ObserverGeometry>({ ...DEFAULT_OBSERVER, ...initialObserver })

  // Update observer ref when state changes
  useEffect(() => {
    observerRef.current = observer
    onObserverChange?.(observer)
  }, [observer, onObserverChange])

  // Color map for p-v diagram (plasma-like)
  const getColor = useCallback((t: number): [number, number, number] => {
    // Plasma colormap
    const r = Math.min(1, 0.5 + 1.5 * t)
    const g = Math.max(0, Math.min(1, 1.5 * t - 0.5))
    const b = Math.max(0, 1 - 2 * t)
    return [r, g, b]
  }, [])

  // Render P-V diagram to canvas
  const renderPVDiagram = useCallback(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    const width = canvas.width
    const height = canvas.height
    const dpr = window.devicePixelRatio || 1

    // Clear background
    ctx.fillStyle = '#0a0a0f'
    ctx.fillRect(0, 0, width, height)

    const pvData = pvDataRef.current
    if (!pvData || pvData.maxCount === 0) {
      ctx.fillStyle = '#6b7280'
      ctx.font = `${12 * dpr}px system-ui`
      ctx.textAlign = 'center'

      // Check why there's no data
      const framesRefCurrent = internalFramesRef.current?.current
      const frameIdx = internalFrameIndexRef.current?.current ?? 0

      if (!framesRefCurrent || framesRefCurrent.size === 0) {
        ctx.fillText('Waiting for frames...', width / 2, height / 2)
      } else if (!framesRefCurrent.get(frameIdx)) {
        ctx.fillText(`Frame ${frameIdx} not loaded`, width / 2, height / 2)
      } else {
        ctx.fillText('Computing...', width / 2, height / 2)
      }
      return
    }

    // Margins for axes
    const marginLeft = 50 * dpr
    const marginRight = 20 * dpr
    const marginTop = 30 * dpr
    const marginBottom = 40 * dpr
    const plotWidth = width - marginLeft - marginRight
    const plotHeight = height - marginTop - marginBottom

    // Draw p-v image
    const cellWidth = plotWidth / pvData.numPosBins
    const cellHeight = plotHeight / pvData.numVelBins

    for (let j = 0; j < pvData.numVelBins; j++) {
      for (let i = 0; i < pvData.numPosBins; i++) {
        const val = pvData.histogram[j * pvData.numPosBins + i]
        if (val <= 0) continue

        // Log scale for better dynamic range
        const t = Math.log10(1 + val) / Math.log10(1 + pvData.maxCount)
        const [r, g, b] = getColor(t)

        ctx.fillStyle = `rgb(${Math.floor(r * 255)}, ${Math.floor(g * 255)}, ${Math.floor(b * 255)})`
        ctx.fillRect(
          marginLeft + i * cellWidth,
          marginTop + (pvData.numVelBins - 1 - j) * cellHeight, // Flip Y
          cellWidth + 0.5,
          cellHeight + 0.5
        )
      }
    }

    // Draw axes
    ctx.strokeStyle = '#4b5563'
    ctx.lineWidth = dpr
    ctx.beginPath()
    ctx.moveTo(marginLeft, marginTop)
    ctx.lineTo(marginLeft, marginTop + plotHeight)
    ctx.lineTo(marginLeft + plotWidth, marginTop + plotHeight)
    ctx.stroke()

    // Axis labels
    ctx.fillStyle = '#9ca3af'
    ctx.font = `${10 * dpr}px system-ui`

    // X-axis (position)
    ctx.textAlign = 'center'
    ctx.fillText('Position (pc)', marginLeft + plotWidth / 2, height - 5 * dpr)
    const posMin = posRange[0], posMax = posRange[1]
    for (let p = posMin; p <= posMax; p += (posMax - posMin) / 4) {
      const x = marginLeft + ((p - posMin) / (posMax - posMin)) * plotWidth
      ctx.fillText(p.toFixed(1), x, marginTop + plotHeight + 15 * dpr)
      ctx.beginPath()
      ctx.moveTo(x, marginTop + plotHeight)
      ctx.lineTo(x, marginTop + plotHeight + 5 * dpr)
      ctx.stroke()
    }

    // Y-axis (velocity)
    ctx.save()
    ctx.translate(12 * dpr, marginTop + plotHeight / 2)
    ctx.rotate(-Math.PI / 2)
    ctx.textAlign = 'center'
    ctx.fillText('V_LSR (km/s)', 0, 0)
    ctx.restore()

    ctx.textAlign = 'right'
    const velMin = velRange[0], velMax = velRange[1]
    for (let v = velMin; v <= velMax; v += (velMax - velMin) / 4) {
      const y = marginTop + ((velMax - v) / (velMax - velMin)) * plotHeight
      ctx.fillText(v.toFixed(0), marginLeft - 5 * dpr, y + 3 * dpr)
      ctx.beginPath()
      ctx.moveTo(marginLeft - 5 * dpr, y)
      ctx.lineTo(marginLeft, y)
      ctx.stroke()
    }

    // Title with observer parameters
    ctx.fillStyle = '#d1d5db'
    ctx.font = `bold ${11 * dpr}px system-ui`
    ctx.textAlign = 'left'
    ctx.fillText('P-V Diagram', marginLeft, 15 * dpr)

    ctx.fillStyle = '#9ca3af'
    ctx.font = `${9 * dpr}px system-ui`
    const obs = observerRef.current
    ctx.fillText(
      `i=${obs.inclination.toFixed(0)}° PA=${obs.positionAngle.toFixed(1)}° V_LSR=${obs.vLSR.toFixed(0)} km/s`,
      marginLeft + 80 * dpr,
      15 * dpr
    )

    // Drag hint
    ctx.fillStyle = '#6b7280'
    ctx.font = `${8 * dpr}px system-ui`
    ctx.textAlign = 'right'
    ctx.fillText('Drag to change viewing angle', width - 10 * dpr, 15 * dpr)

  }, [posRange, velRange, getColor])

  // Compute P-V diagram from current frame
  const updatePVDiagram = useCallback(() => {
    // Access internal refs to get latest values (avoid stale closures)
    const framesRefCurrent = internalFramesRef.current?.current
    const frameIndexRefCurrent = internalFrameIndexRef.current?.current ?? 0

    // Check if frames map exists and has data
    if (!framesRefCurrent || framesRefCurrent.size === 0) {
      return
    }

    const frame = framesRefCurrent.get(frameIndexRefCurrent)
    if (!frame || !frame.positions || frame.particleCount === 0) {
      return
    }

    // Throttle computation (but always compute on observer change)
    const now = performance.now()
    const observerChanged = lastUpdateTimeRef.current === 0
    if (!observerChanged && now - lastUpdateTimeRef.current < 100 && frameIndex === lastFrameIndexRef.current) {
      return
    }
    lastUpdateTimeRef.current = now
    lastFrameIndexRef.current = frameIndex

    // Compute P-V diagram
    pvDataRef.current = computePVDiagram(
      frame,
      observerRef.current,
      posRange,
      velRange,
      numPosBins,
      numVelBins
    )

    // Render
    renderPVDiagram()
  }, [posRange, velRange, numPosBins, numVelBins, renderPVDiagram])

  // Animation loop
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    // Set canvas size with DPR
    const dpr = window.devicePixelRatio || 1
    const updateCanvasSize = () => {
      const rect = canvas.getBoundingClientRect()
      canvas.width = rect.width * dpr
      canvas.height = rect.height * dpr
    }
    updateCanvasSize()

    // Animation loop
    const animate = () => {
      animationIdRef.current = requestAnimationFrame(animate)
      updatePVDiagram()
    }
    animate()

    // Handle resize
    const handleResize = () => {
      updateCanvasSize()
      renderPVDiagram()
    }
    window.addEventListener('resize', handleResize)

    return () => {
      cancelAnimationFrame(animationIdRef.current)
      window.removeEventListener('resize', handleResize)
    }
  }, [updatePVDiagram, renderPVDiagram])

  // Mouse drag to change observer geometry
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const handleMouseDown = (e: MouseEvent) => {
      isDraggingRef.current = true
      dragStartRef.current = {
        x: e.clientX,
        y: e.clientY,
        inc: observerRef.current.inclination,
        pa: observerRef.current.positionAngle,
      }
      canvas.style.cursor = 'grabbing'
    }

    const handleMouseMove = (e: MouseEvent) => {
      if (!isDraggingRef.current || !dragStartRef.current) return

      const dx = e.clientX - dragStartRef.current.x
      const dy = e.clientY - dragStartRef.current.y

      // Horizontal drag -> position angle
      // Vertical drag -> inclination
      const newPA = (dragStartRef.current.pa + dx * 0.5 + 360) % 360
      const newInc = Math.max(0, Math.min(90, dragStartRef.current.inc - dy * 0.3))

      setObserver(prev => ({
        ...prev,
        inclination: newInc,
        positionAngle: newPA,
      }))

      // Force immediate recompute
      lastUpdateTimeRef.current = 0
    }

    const handleMouseUp = () => {
      isDraggingRef.current = false
      dragStartRef.current = null
      canvas.style.cursor = 'grab'
    }

    canvas.addEventListener('mousedown', handleMouseDown)
    window.addEventListener('mousemove', handleMouseMove)
    window.addEventListener('mouseup', handleMouseUp)
    canvas.style.cursor = 'grab'

    return () => {
      canvas.removeEventListener('mousedown', handleMouseDown)
      window.removeEventListener('mousemove', handleMouseMove)
      window.removeEventListener('mouseup', handleMouseUp)
    }
  }, [])

  return (
    <div className={`flex flex-col ${className}`}>
      <canvas
        ref={canvasRef}
        className="w-full flex-1 rounded"
        style={{ minHeight: '250px' }}
      />
      {/* Observer controls */}
      <div className="mt-2 grid grid-cols-3 gap-2 text-xs">
        <label className="flex flex-col">
          <span className="text-gray-400">Inclination</span>
          <input
            type="range"
            min="0"
            max="90"
            step="1"
            value={observer.inclination}
            onChange={(e) => setObserver(prev => ({ ...prev, inclination: parseFloat(e.target.value) }))}
            className="w-full"
          />
          <span className="text-gray-500">{observer.inclination.toFixed(0)}°</span>
        </label>
        <label className="flex flex-col">
          <span className="text-gray-400">Position Angle</span>
          <input
            type="range"
            min="0"
            max="360"
            step="1"
            value={observer.positionAngle}
            onChange={(e) => setObserver(prev => ({ ...prev, positionAngle: parseFloat(e.target.value) }))}
            className="w-full"
          />
          <span className="text-gray-500">{observer.positionAngle.toFixed(1)}°</span>
        </label>
        <label className="flex flex-col">
          <span className="text-gray-400">V_LSR</span>
          <input
            type="range"
            min="-200"
            max="0"
            step="5"
            value={observer.vLSR}
            onChange={(e) => setObserver(prev => ({ ...prev, vLSR: parseFloat(e.target.value) }))}
            className="w-full"
          />
          <span className="text-gray-500">{observer.vLSR.toFixed(0)} km/s</span>
        </label>
      </div>
    </div>
  )
}

export default PVDiagramImperative
