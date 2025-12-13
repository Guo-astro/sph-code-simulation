'use client'

import { useRef, useEffect, useCallback } from 'react'
import type { ParsedFrame } from '~/types/sph'

// ============== Constants ==============
const GAMMA = 5 / 3  // Monatomic ideal gas

const COLORS = {
  background: '#0f0f14',
  panel: '#1a1a24',
  grid: '#2a2a3a',
  text: '#e0e0e0',
  accent: '#00b4d8',
  density: '#ffd166',     // Gold
  pressure: '#ef476f',    // Pink-red
  temperature: '#06d6a0', // Mint
  pAdiabatic: '#ffd166',  // Gold (dashed)
  tAdiabatic: '#06d6a0',  // Mint (dashed)
}

// ============== Types ==============
interface ProfileData {
  positions: number[]
  density: number[]
  pressure: number[]
  temperature: number[]
  pAdiabatic: number[]
  tAdiabatic: number[]
}

interface ComputedData {
  com: { x: number; y: number; z: number }
  temperature: Float32Array
  tRatio: Float32Array
  initialValues: { T: number; rho: number; P: number }
  verticalProfile: ProfileData
  horizontalProfile: ProfileData
}

export interface ShockDiagnosticsPanelImperativeProps {
  /** Map of all loaded frames - passed as ref to avoid React re-renders */
  framesRef: React.RefObject<Map<number, ParsedFrame>>
  /** Current frame index ref - read imperatively during animation */
  frameIndexRef: React.RefObject<number>
  /** Initial frame for baseline values (pass frame 0) */
  initialFrame: ParsedFrame | null
  /** Canvas width in pixels */
  width?: number
  /** Canvas height in pixels */
  height?: number
  /** Additional CSS classes */
  className?: string
}

// ============== Physics Functions ==============

/**
 * Compute Center of Mass of the cloud
 */
function computeCoM(frame: ParsedFrame): { x: number; y: number; z: number } {
  let totalMass = 0
  let cx = 0, cy = 0, cz = 0
  
  for (let i = 0; i < frame.particleCount; i++) {
    const m = frame.mass?.[i] ?? 1
    const x = frame.positions[i * 3]
    const y = frame.positions[i * 3 + 1]
    const z = frame.positions[i * 3 + 2]
    
    totalMass += m
    cx += m * x
    cy += m * y
    cz += m * z
  }
  
  if (totalMass > 0) {
    cx /= totalMass
    cy /= totalMass
    cz /= totalMass
  }
  
  return { x: cx, y: cy, z: cz }
}

/**
 * Compute temperature from pressure and density using ideal gas law
 * T = P / (ρ * k_B / (μ * m_H)) ∝ P/ρ
 * For monatomic gas: T = P/ρ * μ*m_H/k_B
 * We normalize to initial values, so T/T₀ = (P/P₀) / (ρ/ρ₀)
 */
function computeTemperature(pressure: Float32Array, density: Float32Array): Float32Array {
  const n = pressure.length
  const temperature = new Float32Array(n)
  
  for (let i = 0; i < n; i++) {
    const rho = density[i]
    const p = pressure[i]
    // T ∝ P/ρ for ideal gas - use arbitrary normalization
    temperature[i] = rho > 0 ? p / rho : 0
  }
  
  return temperature
}

/**
 * Extract 1D vertical profile (Z direction) - shows vertical compression shock
 */
function extractVerticalProfile(
  frame: ParsedFrame,
  com: { x: number; y: number; z: number },
  temperature: Float32Array,
  initial: { T: number; rho: number; P: number },
  columnRadius: number = 0.15,
  nBins: number = 40
): ProfileData {
  const positions: number[] = []
  const density: number[] = []
  const pressure: number[] = []
  const temp: number[] = []
  const pAdiabatic: number[] = []
  const tAdiabatic: number[] = []
  
  // Select particles in a vertical column around CoM
  const particlesInColumn: { z: number; dens: number; pres: number; temp: number; mass: number }[] = []
  
  for (let i = 0; i < frame.particleCount; i++) {
    const dx = frame.positions[i * 3] - com.x
    const dy = frame.positions[i * 3 + 1] - com.y
    const rXY = Math.sqrt(dx * dx + dy * dy)
    
    if (rXY < columnRadius) {
      particlesInColumn.push({
        z: frame.positions[i * 3 + 2] - com.z,
        dens: frame.density?.[i] ?? 1,
        pres: frame.pressure?.[i] ?? 1,
        temp: temperature[i],
        mass: frame.mass?.[i] ?? 1,
      })
    }
  }
  
  if (particlesInColumn.length < 10) {
    return { positions: [], density: [], pressure: [], temperature: [], pAdiabatic: [], tAdiabatic: [] }
  }
  
  // Bin the data
  const zMin = Math.min(...particlesInColumn.map(p => p.z))
  const zMax = Math.max(...particlesInColumn.map(p => p.z))
  const binWidth = (zMax - zMin) / nBins
  
  for (let bin = 0; bin < nBins; bin++) {
    const binStart = zMin + bin * binWidth
    const binEnd = binStart + binWidth
    const binCenter = (binStart + binEnd) / 2
    
    const inBin = particlesInColumn.filter(p => p.z >= binStart && p.z < binEnd)
    if (inBin.length === 0) continue
    
    const totalMass = inBin.reduce((sum, p) => sum + p.mass, 0)
    const avgDens = inBin.reduce((sum, p) => sum + p.mass * p.dens, 0) / totalMass
    const avgPres = inBin.reduce((sum, p) => sum + p.mass * p.pres, 0) / totalMass
    const avgTemp = inBin.reduce((sum, p) => sum + p.mass * p.temp, 0) / totalMass
    
    const rhoRatio = avgDens / initial.rho
    
    positions.push(binCenter)
    density.push(rhoRatio)
    pressure.push(avgPres / initial.P)
    temp.push(avgTemp / initial.T)
    pAdiabatic.push(Math.pow(rhoRatio, GAMMA))  // P/P₀ = (ρ/ρ₀)^γ
    tAdiabatic.push(Math.pow(rhoRatio, GAMMA - 1))  // T/T₀ = (ρ/ρ₀)^(γ-1)
  }
  
  return { positions, density, pressure, temperature: temp, pAdiabatic, tAdiabatic }
}

/**
 * Extract 1D horizontal profile (X direction) - shows tidal stretching
 */
function extractHorizontalProfile(
  frame: ParsedFrame,
  com: { x: number; y: number; z: number },
  temperature: Float32Array,
  initial: { T: number; rho: number; P: number },
  sliceThickness: number = 0.15,
  nBins: number = 60
): ProfileData {
  const positions: number[] = []
  const density: number[] = []
  const pressure: number[] = []
  const temp: number[] = []
  const pAdiabatic: number[] = []
  const tAdiabatic: number[] = []
  
  // Select particles in a horizontal slice
  const particlesInSlice: { x: number; dens: number; pres: number; temp: number; mass: number }[] = []
  
  for (let i = 0; i < frame.particleCount; i++) {
    const dz = Math.abs(frame.positions[i * 3 + 2] - com.z)
    const dy = Math.abs(frame.positions[i * 3 + 1] - com.y)
    
    if (dz < sliceThickness && dy < sliceThickness) {
      particlesInSlice.push({
        x: frame.positions[i * 3] - com.x,
        dens: frame.density?.[i] ?? 1,
        pres: frame.pressure?.[i] ?? 1,
        temp: temperature[i],
        mass: frame.mass?.[i] ?? 1,
      })
    }
  }
  
  if (particlesInSlice.length < 10) {
    return { positions: [], density: [], pressure: [], temperature: [], pAdiabatic: [], tAdiabatic: [] }
  }
  
  // Bin the data
  const xMin = Math.min(...particlesInSlice.map(p => p.x))
  const xMax = Math.max(...particlesInSlice.map(p => p.x))
  const binWidth = (xMax - xMin) / nBins
  
  for (let bin = 0; bin < nBins; bin++) {
    const binStart = xMin + bin * binWidth
    const binEnd = binStart + binWidth
    const binCenter = (binStart + binEnd) / 2
    
    const inBin = particlesInSlice.filter(p => p.x >= binStart && p.x < binEnd)
    if (inBin.length === 0) continue
    
    const totalMass = inBin.reduce((sum, p) => sum + p.mass, 0)
    const avgDens = inBin.reduce((sum, p) => sum + p.mass * p.dens, 0) / totalMass
    const avgPres = inBin.reduce((sum, p) => sum + p.mass * p.pres, 0) / totalMass
    const avgTemp = inBin.reduce((sum, p) => sum + p.mass * p.temp, 0) / totalMass
    
    const rhoRatio = avgDens / initial.rho
    
    positions.push(binCenter)
    density.push(rhoRatio)
    pressure.push(avgPres / initial.P)
    temp.push(avgTemp / initial.T)
    pAdiabatic.push(Math.pow(rhoRatio, GAMMA))
    tAdiabatic.push(Math.pow(rhoRatio, GAMMA - 1))
  }
  
  return { positions, density, pressure, temperature: temp, pAdiabatic, tAdiabatic }
}

/**
 * Draw a 1D profile chart
 */
function drawProfile(
  ctx: CanvasRenderingContext2D,
  profile: ProfileData,
  x: number,
  y: number,
  width: number,
  height: number,
  title: string,
  xLabel: string,
  yRange: [number, number] = [0, 5]
) {
  const margin = { left: 45, right: 10, top: 25, bottom: 30 }
  const plotWidth = width - margin.left - margin.right
  const plotHeight = height - margin.top - margin.bottom
  
  // Background
  ctx.fillStyle = COLORS.panel
  ctx.fillRect(x, y, width, height)
  
  // Title
  ctx.fillStyle = COLORS.text
  ctx.font = 'bold 11px sans-serif'
  ctx.textAlign = 'center'
  ctx.fillText(title, x + width / 2, y + 15)
  
  if (profile.positions.length === 0) {
    ctx.fillStyle = '#666'
    ctx.font = '10px sans-serif'
    ctx.fillText('No data', x + width / 2, y + height / 2)
    return
  }
  
  // Grid
  ctx.strokeStyle = COLORS.grid
  ctx.lineWidth = 0.5
  ctx.setLineDash([2, 2])
  
  for (let i = 0; i <= 4; i++) {
    const gy = y + margin.top + (i / 4) * plotHeight
    ctx.beginPath()
    ctx.moveTo(x + margin.left, gy)
    ctx.lineTo(x + margin.left + plotWidth, gy)
    ctx.stroke()
  }
  ctx.setLineDash([])
  
  // Scale functions
  const xMin = Math.min(...profile.positions)
  const xMax = Math.max(...profile.positions)
  const scaleX = (v: number) => x + margin.left + ((v - xMin) / (xMax - xMin)) * plotWidth
  const scaleY = (v: number) => y + margin.top + plotHeight - ((v - yRange[0]) / (yRange[1] - yRange[0])) * plotHeight
  
  // Reference line at 1.0
  ctx.strokeStyle = '#666'
  ctx.lineWidth = 1
  ctx.setLineDash([4, 4])
  ctx.beginPath()
  ctx.moveTo(x + margin.left, scaleY(1))
  ctx.lineTo(x + margin.left + plotWidth, scaleY(1))
  ctx.stroke()
  ctx.setLineDash([])
  
  // Draw lines
  const drawLine = (data: number[], color: string, dashed: boolean = false) => {
    if (data.length === 0) return
    ctx.strokeStyle = color
    ctx.lineWidth = dashed ? 1.5 : 2
    if (dashed) ctx.setLineDash([4, 2])
    ctx.beginPath()
    for (let i = 0; i < data.length; i++) {
      const px = scaleX(profile.positions[i])
      const py = scaleY(Math.max(yRange[0], Math.min(yRange[1], data[i])))
      if (i === 0) ctx.moveTo(px, py)
      else ctx.lineTo(px, py)
    }
    ctx.stroke()
    ctx.setLineDash([])
  }
  
  // Adiabatic references (dashed)
  drawLine(profile.pAdiabatic, COLORS.pAdiabatic, true)
  drawLine(profile.tAdiabatic, COLORS.tAdiabatic, true)
  
  // Actual values (solid)
  drawLine(profile.density, COLORS.density)
  drawLine(profile.pressure, COLORS.pressure)
  drawLine(profile.temperature, COLORS.temperature)
  
  // Axes
  ctx.strokeStyle = COLORS.accent
  ctx.lineWidth = 1
  ctx.beginPath()
  ctx.moveTo(x + margin.left, y + margin.top)
  ctx.lineTo(x + margin.left, y + height - margin.bottom)
  ctx.lineTo(x + width - margin.right, y + height - margin.bottom)
  ctx.stroke()
  
  // Labels
  ctx.fillStyle = COLORS.text
  ctx.font = '9px sans-serif'
  ctx.textAlign = 'center'
  ctx.fillText(xLabel, x + width / 2, y + height - 5)
  
  ctx.textAlign = 'right'
  ctx.fillText(yRange[1].toFixed(0), x + margin.left - 5, y + margin.top + 5)
  ctx.fillText(yRange[0].toFixed(0), x + margin.left - 5, y + height - margin.bottom)
}

/**
 * Compute initial values for normalization (cached once)
 */
function computeInitialValues(initialFrame: ParsedFrame | null, frame: ParsedFrame | null): { T: number; rho: number; P: number } {
  const refFrame = initialFrame ?? frame
  if (!refFrame) return { T: 50, rho: 1e-3, P: 1e-6 }
  
  // Compute median values
  const densities = Array.from(refFrame.density ?? []).filter(v => v > 0)
  const pressures = Array.from(refFrame.pressure ?? []).filter(v => v > 0)
  
  densities.sort((a, b) => a - b)
  pressures.sort((a, b) => a - b)
  
  const medianRho = densities[Math.floor(densities.length / 2)] || 1e-3
  const medianP = pressures[Math.floor(pressures.length / 2)] || 1e-6
  
  // Compute temperature
  const temp = computeTemperature(
    new Float32Array([medianP]),
    new Float32Array([medianRho])
  )
  
  return { T: temp[0], rho: medianRho, P: medianP }
}

/**
 * Compute all derived data for a frame
 */
function computeFrameData(
  frame: ParsedFrame,
  initialValues: { T: number; rho: number; P: number }
): ComputedData {
  const com = computeCoM(frame)
  const temperature = computeTemperature(
    frame.pressure ?? new Float32Array(frame.particleCount),
    frame.density ?? new Float32Array(frame.particleCount)
  )
  
  const tRatio = new Float32Array(frame.particleCount)
  for (let i = 0; i < frame.particleCount; i++) {
    tRatio[i] = temperature[i] / initialValues.T
  }
  
  const verticalProfile = extractVerticalProfile(frame, com, temperature, initialValues)
  const horizontalProfile = extractHorizontalProfile(frame, com, temperature, initialValues)
  
  return {
    com,
    temperature,
    tRatio,
    initialValues,
    verticalProfile,
    horizontalProfile,
  }
}

/**
 * High-performance Shock Diagnostics Panel (Imperative Mode)
 * 
 * Uses requestAnimationFrame to read frame data from refs, avoiding
 * React re-renders during playback. This enables smooth animation.
 */
export function ShockDiagnosticsPanelImperative({
  framesRef,
  frameIndexRef,
  initialFrame,
  width = 500,
  height = 300,
  className = '',
}: ShockDiagnosticsPanelImperativeProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const animationFrameRef = useRef<number | null>(null)
  const lastFrameIndexRef = useRef<number>(-1)
  const initialValuesRef = useRef<{ T: number; rho: number; P: number } | null>(null)

  // Main draw function
  const draw = useCallback((frame: ParsedFrame) => {
    const canvas = canvasRef.current
    if (!canvas) return

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    // Compute initial values once
    if (!initialValuesRef.current) {
      initialValuesRef.current = computeInitialValues(initialFrame, frame)
    }
    const initialValues = initialValuesRef.current

    // Compute all derived data
    const computedData = computeFrameData(frame, initialValues)

    const dpr = window.devicePixelRatio || 1
    canvas.width = width * dpr
    canvas.height = height * dpr
    ctx.scale(dpr, dpr)

    // Background
    ctx.fillStyle = COLORS.background
    ctx.fillRect(0, 0, width, height)

    // Layout: two profile charts side by side
    const chartWidth = (width - 30) / 2
    const chartHeight = height - 60

    // CoM info header
    ctx.fillStyle = COLORS.text
    ctx.font = 'bold 10px monospace'
    ctx.textAlign = 'left'
    ctx.fillText(
      `CoM: (${computedData.com.x.toFixed(2)}, ${computedData.com.y.toFixed(2)}, ${computedData.com.z.toFixed(2)}) pc`,
      10, 15
    )

    // IMBH distance
    const distImbh = Math.sqrt(
      computedData.com.x ** 2 + computedData.com.y ** 2 + computedData.com.z ** 2
    )
    ctx.fillText(`IMBH→CoM: ${distImbh.toFixed(2)} pc`, 10, 28)

    // Extrema info
    let maxTRatio = 0, maxRhoRatio = 0, maxPRatio = 0
    for (let i = 0; i < frame.particleCount; i++) {
      if (computedData.tRatio[i] > maxTRatio) maxTRatio = computedData.tRatio[i]
      const rhoRatio = (frame.density?.[i] ?? 0) / initialValues.rho
      if (rhoRatio > maxRhoRatio) maxRhoRatio = rhoRatio
      const pRatio = (frame.pressure?.[i] ?? 0) / initialValues.P
      if (pRatio > maxPRatio) maxPRatio = pRatio
    }

    ctx.textAlign = 'right'
    ctx.font = '9px monospace'
    ctx.fillStyle = COLORS.temperature
    ctx.fillText(`T_max/T₀=${maxTRatio.toFixed(1)}`, width - 10, 15)
    ctx.fillStyle = COLORS.density
    ctx.fillText(`ρ_max/ρ₀=${maxRhoRatio.toFixed(1)}`, width - 10, 26)
    ctx.fillStyle = COLORS.pressure
    ctx.fillText(`P_max/P₀=${maxPRatio.toFixed(1)}`, width - 10, 37)

    // Draw profiles
    drawProfile(
      ctx,
      computedData.verticalProfile,
      10, 45, chartWidth, chartHeight,
      'VERTICAL SHOCK (Z)',
      'Z - Z_CoM (pc)',
      [0, 5]
    )

    drawProfile(
      ctx,
      computedData.horizontalProfile,
      20 + chartWidth, 45, chartWidth, chartHeight,
      'HORIZONTAL STRETCH (X)',
      'X - X_CoM (pc)',
      [0, 5]
    )

    // Legend
    const legendY = height - 12
    ctx.font = '8px sans-serif'
    ctx.textAlign = 'left'

    const legendItems = [
      { color: COLORS.density, label: 'ρ/ρ₀' },
      { color: COLORS.pressure, label: 'P/P₀' },
      { color: COLORS.temperature, label: 'T/T₀' },
      { color: COLORS.pAdiabatic, label: 'P_ad', dashed: true },
      { color: COLORS.tAdiabatic, label: 'T_ad', dashed: true },
    ]

    let legendX = 10
    for (const item of legendItems) {
      ctx.fillStyle = item.color
      ctx.fillRect(legendX, legendY - 5, item.dashed ? 15 : 12, 2)
      if (item.dashed) {
        ctx.fillStyle = COLORS.background
        ctx.fillRect(legendX + 4, legendY - 5, 3, 2)
        ctx.fillRect(legendX + 10, legendY - 5, 2, 2)
      }
      ctx.fillStyle = COLORS.text
      ctx.fillText(item.label, legendX + (item.dashed ? 18 : 15), legendY)
      legendX += 55
    }
  }, [initialFrame, width, height])

  // Imperative animation loop
  useEffect(() => {
    const tick = () => {
      const frames = framesRef.current
      const frameIndex = frameIndexRef.current ?? 0

      // Only redraw if frame index changed
      if (frames && frameIndex !== lastFrameIndexRef.current) {
        const frame = frames.get(frameIndex)
        if (frame) {
          draw(frame)
          lastFrameIndexRef.current = frameIndex
        }
      }

      animationFrameRef.current = requestAnimationFrame(tick)
    }

    animationFrameRef.current = requestAnimationFrame(tick)

    return () => {
      if (animationFrameRef.current !== null) {
        cancelAnimationFrame(animationFrameRef.current)
      }
    }
  }, [framesRef, frameIndexRef, draw])

  // Initial draw when initialFrame changes
  useEffect(() => {
    initialValuesRef.current = null // Reset initial values
    const frames = framesRef.current
    const frameIndex = frameIndexRef.current ?? 0
    if (frames) {
      const frame = frames.get(frameIndex)
      if (frame) {
        draw(frame)
      }
    }
  }, [initialFrame, framesRef, frameIndexRef, draw])

  // Redraw when dimensions change
  useEffect(() => {
    const frames = framesRef.current
    const frameIndex = frameIndexRef.current ?? 0
    if (frames) {
      const frame = frames.get(frameIndex)
      if (frame) {
        draw(frame)
      }
    }
  }, [width, height, framesRef, frameIndexRef, draw])

  return (
    <div className={`relative ${className}`}>
      <canvas
        ref={canvasRef}
        style={{ width, height }}
        className="rounded border border-gray-700"
      />
    </div>
  )
}

export default ShockDiagnosticsPanelImperative
