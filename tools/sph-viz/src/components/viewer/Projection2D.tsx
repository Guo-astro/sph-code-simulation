'use client'

import { useRef, useEffect, useMemo, useCallback } from 'react'
import type { ParsedFrame, ColorMap } from '~/types/sph'

export interface Projection2DProps {
  frame: ParsedFrame | null
  projection: 'xy' | 'xz' | 'yz'
  colorField?: string
  colorMap?: ColorMap
  width?: number
  height?: number
  showAxes?: boolean
  showColorbar?: boolean
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

/** Interpolate between colors in a color map */
function interpolateColorHex(colors: string[], t: number): string {
  if (colors.length === 0) return '#ffffff'
  if (colors.length === 1) return colors[0]

  t = Math.max(0, Math.min(1, t))
  const index = t * (colors.length - 1)
  const lower = Math.floor(index)
  const upper = Math.min(lower + 1, colors.length - 1)
  const localT = index - lower

  const c1 = hexToRgb(colors[lower])
  const c2 = hexToRgb(colors[upper])

  const r = Math.round(c1.r + (c2.r - c1.r) * localT)
  const g = Math.round(c1.g + (c2.g - c1.g) * localT)
  const b = Math.round(c1.b + (c2.b - c1.b) * localT)

  return `rgb(${r},${g},${b})`
}

function hexToRgb(hex: string): { r: number; g: number; b: number } {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex)
  if (!result) return { r: 255, g: 255, b: 255 }
  return {
    r: parseInt(result[1], 16),
    g: parseInt(result[2], 16),
    b: parseInt(result[3], 16),
  }
}

/** 2D Projection Canvas Component */
export function Projection2D({
  frame,
  projection,
  colorField = 'density',
  colorMap = defaultColorMap,
  width = 400,
  height = 400,
  showAxes = true,
  showColorbar = true,
  className = '',
}: Projection2DProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null)

  // Get axis indices based on projection
  const [axisX, axisY, axisLabel] = useMemo(() => {
    switch (projection) {
      case 'xy':
        return [0, 1, ['X', 'Y']]
      case 'xz':
        return [0, 2, ['X', 'Z']]
      case 'yz':
        return [1, 2, ['Y', 'Z']]
      default:
        return [0, 1, ['X', 'Y']]
    }
  }, [projection])

  // Extract field data and compute bounds
  const { fieldData, bounds, range } = useMemo(() => {
    if (!frame) {
      return {
        fieldData: null,
        bounds: { minX: -1, maxX: 1, minY: -1, maxY: 1 },
        range: { min: 0, max: 1 },
      }
    }

    // Get color field data
    let data: Float32Array | undefined
    switch (colorField) {
      case 'density':
        data = frame.density
        break
      case 'pressure':
        data = frame.pressure
        break
      case 'energy':
        data = frame.energy
        break
      case 'velocity':
        data = new Float32Array(frame.particleCount)
        for (let i = 0; i < frame.particleCount; i++) {
          const vx = frame.velocities[i * 3]
          const vy = frame.velocities[i * 3 + 1]
          const vz = frame.velocities[i * 3 + 2]
          data[i] = Math.sqrt(vx * vx + vy * vy + vz * vz)
        }
        break
      case 'machNumber':
        data = frame.machNumber
        break
      default:
        data = frame.density
    }

    // Compute bounds
    let minX = Infinity,
      maxX = -Infinity
    let minY = Infinity,
      maxY = -Infinity
    let minVal = Infinity,
      maxVal = -Infinity

    for (let i = 0; i < frame.particleCount; i++) {
      const x = frame.positions[i * 3 + axisX]
      const y = frame.positions[i * 3 + axisY]
      const val = data ? data[i] : 0

      if (isFinite(x)) {
        if (x < minX) minX = x
        if (x > maxX) maxX = x
      }
      if (isFinite(y)) {
        if (y < minY) minY = y
        if (y > maxY) maxY = y
      }
      if (isFinite(val)) {
        if (val < minVal) minVal = val
        if (val > maxVal) maxVal = val
      }
    }

    // Add padding
    const padX = (maxX - minX) * 0.05
    const padY = (maxY - minY) * 0.05
    minX -= padX
    maxX += padX
    minY -= padY
    maxY += padY

    return {
      fieldData: data,
      bounds: { minX, maxX, minY, maxY },
      range: {
        min: colorMap.min !== undefined ? colorMap.min : minVal,
        max: colorMap.max !== undefined ? colorMap.max : maxVal,
      },
    }
  }, [frame, colorField, colorMap, axisX, axisY])

  // Draw function
  const draw = useCallback(() => {
    const canvas = canvasRef.current
    if (!canvas || !frame) return

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    const dpr = window.devicePixelRatio || 1
    canvas.width = width * dpr
    canvas.height = height * dpr
    ctx.scale(dpr, dpr)

    // Clear
    ctx.fillStyle = '#1a1a2e'
    ctx.fillRect(0, 0, width, height)

    const margin = showAxes ? 40 : 10
    const plotWidth = width - margin * 2
    const plotHeight = height - margin * 2

    // Scale functions
    const scaleX = (val: number) =>
      margin + ((val - bounds.minX) / (bounds.maxX - bounds.minX)) * plotWidth
    const scaleY = (val: number) =>
      height - margin - ((val - bounds.minY) / (bounds.maxY - bounds.minY)) * plotHeight

    // Draw grid
    if (showAxes) {
      ctx.strokeStyle = '#333'
      ctx.lineWidth = 0.5
      ctx.beginPath()

      // Vertical grid lines
      const numGridLines = 5
      for (let i = 0; i <= numGridLines; i++) {
        const x = margin + (i / numGridLines) * plotWidth
        ctx.moveTo(x, margin)
        ctx.lineTo(x, height - margin)
      }

      // Horizontal grid lines
      for (let i = 0; i <= numGridLines; i++) {
        const y = margin + (i / numGridLines) * plotHeight
        ctx.moveTo(margin, y)
        ctx.lineTo(width - margin, y)
      }
      ctx.stroke()
    }

    // Draw particles
    const { logScale } = colorMap
    const { min, max } = range

    for (let i = 0; i < frame.particleCount; i++) {
      const x = frame.positions[i * 3 + axisX]
      const y = frame.positions[i * 3 + axisY]
      const val = fieldData ? fieldData[i] : 0

      if (!isFinite(x) || !isFinite(y)) continue

      // Normalize value
      let t: number
      if (logScale && min > 0 && val > 0) {
        t = (Math.log10(val) - Math.log10(min)) / (Math.log10(max) - Math.log10(min))
      } else {
        t = (val - min) / (max - min)
      }
      t = Math.max(0, Math.min(1, isFinite(t) ? t : 0))

      const color = interpolateColorHex(colorMap.colors, t)
      const px = scaleX(x)
      const py = scaleY(y)

      ctx.fillStyle = color
      ctx.beginPath()
      ctx.arc(px, py, 1.5, 0, Math.PI * 2)
      ctx.fill()
    }

    // Draw axes
    if (showAxes) {
      ctx.strokeStyle = '#888'
      ctx.lineWidth = 1
      ctx.beginPath()
      ctx.moveTo(margin, height - margin)
      ctx.lineTo(width - margin, height - margin)
      ctx.moveTo(margin, height - margin)
      ctx.lineTo(margin, margin)
      ctx.stroke()

      // Labels
      ctx.fillStyle = '#ccc'
      ctx.font = '12px sans-serif'
      ctx.textAlign = 'center'

      // X axis label
      ctx.fillText(axisLabel[0], width / 2, height - 5)

      // Y axis label
      ctx.save()
      ctx.translate(12, height / 2)
      ctx.rotate(-Math.PI / 2)
      ctx.fillText(axisLabel[1], 0, 0)
      ctx.restore()

      // Tick labels
      ctx.font = '10px sans-serif'
      ctx.fillText(bounds.minX.toFixed(2), margin, height - 25)
      ctx.fillText(bounds.maxX.toFixed(2), width - margin, height - 25)
      ctx.textAlign = 'right'
      ctx.fillText(bounds.minY.toFixed(2), margin - 5, height - margin)
      ctx.fillText(bounds.maxY.toFixed(2), margin - 5, margin + 5)
    }

    // Draw colorbar
    if (showColorbar) {
      const barWidth = 15
      const barHeight = plotHeight * 0.6
      const barX = width - margin + 10
      const barY = margin + (plotHeight - barHeight) / 2

      // Gradient
      for (let i = 0; i < barHeight; i++) {
        const t = 1 - i / barHeight
        ctx.fillStyle = interpolateColorHex(colorMap.colors, t)
        ctx.fillRect(barX, barY + i, barWidth, 1)
      }

      // Border
      ctx.strokeStyle = '#888'
      ctx.strokeRect(barX, barY, barWidth, barHeight)

      // Labels
      ctx.fillStyle = '#ccc'
      ctx.font = '9px sans-serif'
      ctx.textAlign = 'left'

      const formatVal = (v: number) => {
        if (Math.abs(v) < 0.001 || Math.abs(v) > 10000) {
          return v.toExponential(1)
        }
        return v.toFixed(3)
      }

      ctx.fillText(formatVal(range.max), barX + barWidth + 3, barY + 8)
      ctx.fillText(formatVal(range.min), barX + barWidth + 3, barY + barHeight)
    }
  }, [
    frame,
    fieldData,
    bounds,
    range,
    colorMap,
    width,
    height,
    showAxes,
    showColorbar,
    axisX,
    axisY,
    axisLabel,
  ])

  // Redraw when dependencies change
  useEffect(() => {
    draw()
  }, [draw])

  return (
    <div className={`relative ${className}`}>
      <canvas
        ref={canvasRef}
        style={{ width, height }}
        className="rounded"
      />
      <div className="absolute top-1 left-1 text-white text-xs bg-black/50 px-1 rounded">
        {projection.toUpperCase()} Projection
      </div>
    </div>
  )
}

export default Projection2D
