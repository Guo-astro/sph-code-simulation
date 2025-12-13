'use client'

import { useRef, useEffect, useMemo, useCallback } from 'react'
import type { ParsedFrame } from '~/types/sph'

// Physical constants in CGS
const k_B = 1.38e-16 // erg/K
const m_H = 1.67e-24 // g
const G = 6.67e-8 // cm^3/(g s^2)
const mu = 2.8 // mean molecular weight for H2 + He

// Unit conversions
const M_SUN = 2e33 // g
const PC_TO_CM = 3.086e18 // cm
const KMS_TO_CMS = 1e5 // cm/s
const YR_TO_S = 3.15e7 // s
const MYR_TO_S = 1e6 * YR_TO_S

// Default IMBH mass (code units: 100 = 1e5 Msun)
const M_BH_CODE = 100
const M_BH_MSUN = M_BH_CODE * 1000

interface TimescaleStats {
  totalParticles: number
  equilibriumCount: number
  equilibriumFraction: number
  medianCoolRatio: number
  medianChemRatio: number
  coolRatioHistogram: { center: number; count: number }[]
  chemRatioHistogram: { center: number; count: number }[]
}

export interface TimescaleDiagnosticsProps {
  /** Ref to current frame data - updated imperatively */
  frameRef: React.RefObject<ParsedFrame | null>
  /** BH position in code units (pc) */
  bhPosition?: [number, number, number]
  className?: string
  /** Callback when stats are computed */
  onStatsUpdate?: (stats: TimescaleStats) => void
}

/**
 * Compute timescales for a subset of particles (for performance)
 */
function computeTimescaleStats(
  frame: ParsedFrame,
  bhPosition: [number, number, number],
  sampleRate: number = 10
): TimescaleStats {
  const n = frame.particleCount
  const coolRatios: number[] = []
  const chemRatios: number[] = []
  let eqCount = 0

  for (let i = 0; i < n; i += sampleRate) {
    const x = frame.positions[i * 3]
    const y = frame.positions[i * 3 + 1]
    const z = frame.positions[i * 3 + 2]
    const density = frame.density[i]
    const energy = frame.energy[i]

    if (!isFinite(density) || density <= 0 || !isFinite(energy) || energy <= 0) {
      continue
    }

    // Convert to physical units
    // Code units: length = 1 pc, mass = 1000 Msun, velocity = 1 km/s
    const rho_cgs = density * (1000 * M_SUN) / Math.pow(PC_TO_CM, 3)

    // Temperature from internal energy
    const gamma = 5/3
    const u_cgs = energy * Math.pow(KMS_TO_CMS, 2)
    const T = Math.max(10, (gamma - 1) * u_cgs * mu * m_H / k_B)

    // Distance from BH
    const dx = x - bhPosition[0]
    const dy = y - bhPosition[1]
    const dz = z - bhPosition[2]
    const distance = Math.max(0.1, Math.sqrt(dx*dx + dy*dy + dz*dz))

    // Number density
    const n_H = rho_cgs / (mu * m_H)

    // Cooling timescale (simplified)
    let Lambda_V: number
    if (T > 1e4) {
      Lambda_V = 1e-23 * n_H * n_H * Math.sqrt(T / 1e4)
    } else if (T > 100) {
      Lambda_V = 1e-26 * n_H * n_H * Math.sqrt(T / 100)
    } else {
      Lambda_V = 1e-27 * n_H * n_H * (T / 10)
    }
    const U_V = rho_cgs * u_cgs
    const tau_cool = Lambda_V > 0 ? (U_V / Lambda_V) / MYR_TO_S : 1e6

    // Dynamical timescale (orbital around BH)
    const r_cm = distance * PC_TO_CM
    const M_BH_g = M_BH_MSUN * M_SUN
    const tau_dyn = (2 * Math.PI * Math.sqrt(Math.pow(r_cm, 3) / (G * M_BH_g))) / MYR_TO_S

    // Chemical timescale for H2 formation
    const k_H2 = 3e-17 // cm^3/s
    const tau_chem = (1 / (k_H2 * Math.max(n_H, 1))) / MYR_TO_S

    // Ratios
    const coolRatio = tau_cool / Math.max(tau_dyn, 1e-6)
    const chemRatio = tau_chem / Math.max(tau_dyn, 1e-6)

    if (isFinite(coolRatio) && coolRatio > 0) coolRatios.push(coolRatio)
    if (isFinite(chemRatio) && chemRatio > 0) chemRatios.push(chemRatio)

    if (coolRatio < 0.1 && chemRatio < 0.1) eqCount++
  }

  // Create histograms (log scale)
  const createLogHistogram = (data: number[], bins: number = 12) => {
    if (data.length === 0) return []
    const logData = data.map(d => Math.log10(Math.max(d, 1e-10)))
    const min = Math.min(...logData)
    const max = Math.max(...logData)
    const binWidth = (max - min) / bins

    return Array.from({ length: bins }, (_, i) => {
      const binMin = min + i * binWidth
      const binMax = min + (i + 1) * binWidth
      const center = Math.pow(10, (binMin + binMax) / 2)
      const count = logData.filter(d => d >= binMin && d < binMax).length
      return { center, count }
    })
  }

  // Compute medians
  const sortedCool = [...coolRatios].sort((a, b) => a - b)
  const sortedChem = [...chemRatios].sort((a, b) => a - b)
  const medianCool = sortedCool[Math.floor(sortedCool.length / 2)] || 0
  const medianChem = sortedChem[Math.floor(sortedChem.length / 2)] || 0

  const sampledCount = Math.floor(n / sampleRate)

  return {
    totalParticles: sampledCount,
    equilibriumCount: eqCount,
    equilibriumFraction: sampledCount > 0 ? eqCount / sampledCount : 0,
    medianCoolRatio: medianCool,
    medianChemRatio: medianChem,
    coolRatioHistogram: createLogHistogram(coolRatios),
    chemRatioHistogram: createLogHistogram(chemRatios),
  }
}

/**
 * Timescale Diagnostics Panel - Uses Canvas for fast rendering
 *
 * Displays timescale comparisons to assess whether equilibrium
 * chemistry (KI2000) is applicable for post-processing.
 */
export function TimescaleDiagnostics({
  frameRef,
  bhPosition = [0, 0, 0],
  className = '',
  onStatsUpdate,
}: TimescaleDiagnosticsProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const statsRef = useRef<TimescaleStats | null>(null)
  const lastFrameTimeRef = useRef<number>(0)

  // Compute stats and render - called imperatively
  const updateDiagnostics = useCallback(() => {
    const frame = frameRef.current
    const canvas = canvasRef.current
    if (!frame || !canvas) return

    // Throttle updates to avoid excessive computation
    const now = performance.now()
    if (now - lastFrameTimeRef.current < 200) return
    lastFrameTimeRef.current = now

    // Compute stats
    const stats = computeTimescaleStats(frame, bhPosition, 20)
    statsRef.current = stats
    onStatsUpdate?.(stats)

    // Render to canvas
    const ctx = canvas.getContext('2d')
    if (!ctx) return

    const width = canvas.width
    const height = canvas.height
    const dpr = window.devicePixelRatio || 1

    // Clear
    ctx.fillStyle = '#1f2937'
    ctx.fillRect(0, 0, width, height)

    // Title
    ctx.fillStyle = '#9ca3af'
    ctx.font = `${11 * dpr}px system-ui`
    ctx.fillText('Timescale Diagnostics', 10 * dpr, 15 * dpr)

    // Equilibrium status
    const eqPct = (stats.equilibriumFraction * 100).toFixed(1)
    const isGood = stats.equilibriumFraction > 0.9
    ctx.fillStyle = isGood ? '#22c55e' : stats.equilibriumFraction > 0.5 ? '#eab308' : '#ef4444'
    ctx.font = `bold ${14 * dpr}px system-ui`
    ctx.fillText(`${eqPct}% in equilibrium`, 10 * dpr, 35 * dpr)

    ctx.fillStyle = '#6b7280'
    ctx.font = `${10 * dpr}px system-ui`
    ctx.fillText(
      isGood ? 'KI2000 applicable' : 'Consider time-dependent chemistry',
      10 * dpr,
      50 * dpr
    )

    // Draw histograms
    const histY = 70 * dpr
    const histHeight = 50 * dpr
    const histWidth = (width - 30 * dpr) / 2

    // Cooling ratio histogram
    ctx.fillStyle = '#4b5563'
    ctx.font = `${9 * dpr}px system-ui`
    ctx.fillText('τ_cool / τ_dyn', 10 * dpr, histY - 5 * dpr)

    const coolHist = stats.coolRatioHistogram
    if (coolHist.length > 0) {
      const maxCount = Math.max(...coolHist.map(h => h.count), 1)
      const barWidth = histWidth / coolHist.length
      coolHist.forEach((bin, i) => {
        const barHeight = (bin.count / maxCount) * histHeight
        const x = 10 * dpr + i * barWidth
        const y = histY + histHeight - barHeight
        ctx.fillStyle = bin.center < 0.1 ? '#22c55e' : '#ef4444'
        ctx.fillRect(x, y, barWidth - 1, barHeight)
      })
      // Reference line at 0.1
      ctx.strokeStyle = '#22c55e'
      ctx.setLineDash([3, 3])
      ctx.beginPath()
      const refX = 10 * dpr + (coolHist.findIndex(h => h.center >= 0.1) / coolHist.length) * histWidth
      ctx.moveTo(refX, histY)
      ctx.lineTo(refX, histY + histHeight)
      ctx.stroke()
      ctx.setLineDash([])
    }

    // Chemistry ratio histogram
    const chemHistX = 20 * dpr + histWidth
    ctx.fillStyle = '#4b5563'
    ctx.fillText('τ_chem / τ_dyn', chemHistX, histY - 5 * dpr)

    const chemHist = stats.chemRatioHistogram
    if (chemHist.length > 0) {
      const maxCount = Math.max(...chemHist.map(h => h.count), 1)
      const barWidth = histWidth / chemHist.length
      chemHist.forEach((bin, i) => {
        const barHeight = (bin.count / maxCount) * histHeight
        const x = chemHistX + i * barWidth
        const y = histY + histHeight - barHeight
        ctx.fillStyle = bin.center < 0.1 ? '#22c55e' : '#f59e0b'
        ctx.fillRect(x, y, barWidth - 1, barHeight)
      })
    }

    // Median values
    const medianY = histY + histHeight + 20 * dpr
    ctx.fillStyle = '#9ca3af'
    ctx.font = `${9 * dpr}px system-ui`
    ctx.fillText(`Median cool ratio: ${stats.medianCoolRatio.toExponential(1)}`, 10 * dpr, medianY)
    ctx.fillText(`Median chem ratio: ${stats.medianChemRatio.toExponential(1)}`, chemHistX, medianY)

  }, [frameRef, bhPosition, onStatsUpdate])

  // Setup canvas and animation
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    // Set canvas size with DPR
    const dpr = window.devicePixelRatio || 1
    const rect = canvas.getBoundingClientRect()
    canvas.width = rect.width * dpr
    canvas.height = rect.height * dpr

    // Initial render
    updateDiagnostics()

    // Update periodically (not on every frame - expensive)
    const intervalId = setInterval(updateDiagnostics, 500)

    return () => clearInterval(intervalId)
  }, [updateDiagnostics])

  return (
    <div className={`bg-gray-800 rounded overflow-hidden ${className}`}>
      <canvas
        ref={canvasRef}
        className="w-full h-full"
        style={{ minHeight: '180px' }}
      />
    </div>
  )
}

export default TimescaleDiagnostics
