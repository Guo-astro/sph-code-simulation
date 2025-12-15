'use client'

import { useState, useCallback, useEffect } from 'react'
import { createFileRoute, useLoaderData } from '@tanstack/react-router'
import { Dashboard } from '~/components/dashboard/Dashboard'
import type { SimulationMetadata, ParsedFrame, FrameStatistics } from '~/types/sph'
import { Folder, ChevronRight, RefreshCw } from 'lucide-react'

// Loader to fetch simulations on the server
export const Route = createFileRoute('/viz/')({
  component: VisualizationPage,
  loader: async () => {
    try {
      // Direct server-side import - no HTTP round-trip needed during SSR
      const { findSimulationsServer } = await import('~/lib/server/findSimulations')
      const simulations = await findSimulationsServer()
      console.log('[Loader] Found simulations:', simulations.length)
      return { simulations }
    } catch (err) {
      console.error('[Loader] Failed:', err)
      return { simulations: [] }
    }
  },
})

/** Parse binary frame data directly from ArrayBuffer - NO base64 overhead! */
function parseBinaryFrame(
  buffer: ArrayBuffer,
  frameIndex: number,
  time: number,
  particleCount: number,
  stride: number,
  fieldOffsets: Record<string, number>
): ParsedFrame {
  const floats = new Float32Array(buffer)

  // Extract arrays directly - much faster than base64 decode
  const positions = new Float32Array(particleCount * 3)
  const velocities = new Float32Array(particleCount * 3)
  const mass = new Float32Array(particleCount)
  const density = new Float32Array(particleCount)
  const pressure = new Float32Array(particleCount)
  const energy = new Float32Array(particleCount)
  const smoothingLength = new Float32Array(particleCount)

  for (let i = 0; i < particleCount; i++) {
    const base = i * stride

    positions[i * 3] = floats[base + (fieldOffsets.x ?? 0)]
    positions[i * 3 + 1] = floats[base + (fieldOffsets.y ?? 1)]
    positions[i * 3 + 2] = floats[base + (fieldOffsets.z ?? 2)]

    velocities[i * 3] = floats[base + (fieldOffsets.vx ?? 3)]
    velocities[i * 3 + 1] = floats[base + (fieldOffsets.vy ?? 4)]
    velocities[i * 3 + 2] = floats[base + (fieldOffsets.vz ?? 5)]

    mass[i] = floats[base + (fieldOffsets.mass ?? 6)]
    density[i] = floats[base + (fieldOffsets.density ?? 7)]
    pressure[i] = floats[base + (fieldOffsets.pressure ?? 8)]
    energy[i] = floats[base + (fieldOffsets.energy ?? 9)]
    smoothingLength[i] = floats[base + (fieldOffsets.smoothing_length ?? 10)]
  }

  return {
    frameIndex,
    time,
    particleCount,
    positions,
    velocities,
    mass,
    density,
    pressure,
    energy,
    smoothingLength,
  }
}

function VisualizationPage() {
  // Get initial data from loader
  const loaderData = useLoaderData({ from: '/viz/' }) as { simulations: SimulationMetadata[] }
  const [simulations, setSimulations] = useState<SimulationMetadata[]>(loaderData?.simulations || [])
  const [selectedSimulation, setSelectedSimulation] = useState<SimulationMetadata | null>(null)
  const [frames, setFrames] = useState<Map<number, ParsedFrame>>(new Map())
  const [statistics, setStatistics] = useState<FrameStatistics[]>([])
  const [isLoading, setIsLoading] = useState(false)
  const [loadingProgress, setLoadingProgress] = useState({ loaded: 0, total: 0 })
  const [error, setError] = useState<string | null>(null)
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false)
  const [isLoadingSimulations, setIsLoadingSimulations] = useState(false)

  // Update simulations when loader data changes
  useEffect(() => {
    if (loaderData?.simulations) {
      console.log('[useEffect] Setting simulations from loader:', loaderData.simulations.length)
      setSimulations(loaderData.simulations)
    }
  }, [loaderData])

  // Load available simulations (for refresh button)
  const loadSimulations = useCallback(async () => {
    console.log('[loadSimulations] Starting fetch...')
    setIsLoadingSimulations(true)
    try {
      const response = await fetch('/api/simulations')
      console.log('[loadSimulations] Response status:', response.status)
      const data = await response.json()
      console.log('[loadSimulations] Data received:', data)
      console.log('[loadSimulations] Simulations count:', data.simulations?.length || 0)
      setSimulations(data.simulations || [])
    } catch (err) {
      console.error('Failed to load simulations:', err)
      setError('Failed to load simulations')
    } finally {
      setIsLoadingSimulations(false)
    }
  }, [])

  // Note: No longer using useEffect to load simulations - using loader instead

  // Load a single frame from binary endpoint with retry support
  const loadSingleFrame = useCallback(async (
    simId: string,
    frameIndex: number,
    retries = 2
  ): Promise<{ frame: ParsedFrame; stats: FrameStatistics } | null> => {
    for (let attempt = 0; attempt <= retries; attempt++) {
      try {
        // Use format=binary query param for raw binary response
        const response = await fetch(`/api/simulations/${simId}/frames/${frameIndex}?format=binary`)
        if (!response.ok) {
          console.error(`Failed to load frame ${frameIndex}: HTTP ${response.status}`)
          if (attempt < retries) {
            await new Promise(r => setTimeout(r, 100 * (attempt + 1)))
            continue
          }
          return null
        }

        const frameIdx = parseInt(response.headers.get('X-Frame-Index') || '0')
        const time = parseFloat(response.headers.get('X-Frame-Time') || '0')
        const particleCount = parseInt(response.headers.get('X-Particle-Count') || '0')
        const stride = parseInt(response.headers.get('X-Stride') || '11')
        const fieldOffsetsHeader = response.headers.get('X-Field-Offsets')
        const fieldOffsets = fieldOffsetsHeader ? JSON.parse(fieldOffsetsHeader) : {}

        const buffer = await response.arrayBuffer()
        
        // Validate buffer size before parsing
        const expectedBytes = particleCount * stride * 4
        if (buffer.byteLength === 0) {
          console.error(`Frame ${frameIndex}: empty buffer received`)
          if (attempt < retries) {
            await new Promise(r => setTimeout(r, 100 * (attempt + 1)))
            continue
          }
          return null
        }
        if (buffer.byteLength % 4 !== 0) {
          console.error(`Frame ${frameIndex}: buffer size ${buffer.byteLength} not divisible by 4 (attempt ${attempt + 1})`)
          if (attempt < retries) {
            await new Promise(r => setTimeout(r, 100 * (attempt + 1)))
            continue
          }
          return null
        }
        if (particleCount > 0 && buffer.byteLength !== expectedBytes) {
          console.warn(`Frame ${frameIndex}: buffer size mismatch. Expected ${expectedBytes}, got ${buffer.byteLength}`)
        }
        
        const frame = parseBinaryFrame(buffer, frameIdx, time, particleCount, stride, fieldOffsets)
        const stats = computeFrameStatistics(frame)
        
        return { frame, stats }
      } catch (err) {
        console.error(`Failed to load frame ${frameIndex} (attempt ${attempt + 1}):`, err)
        if (attempt < retries) {
          await new Promise(r => setTimeout(r, 100 * (attempt + 1)))
          continue
        }
        return null
      }
    }
    return null
  }, [])

  // Select simulation and LOAD ALL FRAMES INTO MEMORY for instant playback
  const selectSimulation = useCallback(async (sim: SimulationMetadata) => {
    console.log('[selectSimulation] Called with:', sim.id, sim.name)
    setSelectedSimulation(sim)
    setFrames(new Map())
    setStatistics([])
    setError(null)
    setIsLoading(true)
    setLoadingProgress({ loaded: 0, total: sim.totalFrames })

    const simId = encodeURIComponent(sim.id.replace(/\//g, '|'))
    console.log('[selectSimulation] Encoded simId:', simId)
    const totalFrames = sim.totalFrames
    const newFrames = new Map<number, ParsedFrame>()
    const newStats: FrameStatistics[] = []

    console.log(`🚀 Loading all ${totalFrames} frames into memory...`)
    const startTime = performance.now()

    // Small delay to let the server stabilize before parallel loading
    await new Promise(resolve => setTimeout(resolve, 100))

    // Load frames in parallel batches for speed
    const batchSize = 10
    let loadedCount = 0
    for (let batchStart = 0; batchStart < totalFrames; batchStart += batchSize) {
      const batchEnd = Math.min(batchStart + batchSize, totalFrames)
      const batchPromises: Promise<{ idx: number; result: { frame: ParsedFrame; stats: FrameStatistics } | null }>[] = []

      for (let i = batchStart; i < batchEnd; i++) {
        batchPromises.push(
          loadSingleFrame(simId, i).then(result => ({ idx: i, result }))
        )
      }

      const results = await Promise.all(batchPromises)
      for (const { idx, result } of results) {
        if (result) {
          newFrames.set(idx, result.frame)
          newStats[idx] = result.stats
        }
        loadedCount++
      }

      // Update progress for UI
      setLoadingProgress({ loaded: loadedCount, total: totalFrames })
      console.log(`   Loaded ${loadedCount}/${totalFrames} frames`)
    }

    const elapsed = ((performance.now() - startTime) / 1000).toFixed(2)
    console.log(`✅ All ${totalFrames} frames loaded in ${elapsed}s - ready for instant playback!`)

    setFrames(newFrames)
    setStatistics(newStats)
    setIsLoading(false)
  }, [loadSingleFrame])

  // Load a frame (now just returns cached data since all frames are preloaded)
  const loadFrame = useCallback(async (frameIndex: number) => {
    // All frames should already be in memory from selectSimulation
    // This is just a fallback in case something wasn't loaded
    if (!selectedSimulation) return
    if (frames.has(frameIndex)) return
    
    // If bulk loading is in progress, skip - don't compete with it
    if (isLoading) return

    try {
      const simId = encodeURIComponent(selectedSimulation.id.replace(/\//g, '|'))
      
      // Use binary format for maximum performance
      const response = await fetch(`/api/simulations/${simId}/frames/${frameIndex}?format=binary`)
      
      if (!response.ok) {
        console.error(`Failed to load frame ${frameIndex}: HTTP ${response.status}`)
        return
      }

      // Get metadata from headers
      const frameIdx = parseInt(response.headers.get('X-Frame-Index') || '0')
      const time = parseFloat(response.headers.get('X-Frame-Time') || '0')
      const particleCount = parseInt(response.headers.get('X-Particle-Count') || '0')
      const stride = parseInt(response.headers.get('X-Stride') || '11')
      const fieldOffsets = JSON.parse(response.headers.get('X-Field-Offsets') || '{}')

      // Get raw binary data - no base64 decode needed!
      const buffer = await response.arrayBuffer()
      const parsed = parseBinaryFrame(buffer, frameIdx, time, particleCount, stride, fieldOffsets)

      setFrames((prev) => {
        const next = new Map(prev)
        next.set(frameIndex, parsed)
        return next
      })

      // Compute statistics for this frame
      const stats = computeFrameStatistics(parsed)
      setStatistics((prev) => {
        const next = [...prev]
        next[frameIndex] = stats
        return next
      })
    } catch (err) {
      console.error('Failed to load frame:', err)
      // Don't set error for fallback loads - the bulk load should have handled it
    }
  }, [selectedSimulation, frames, isLoading])

  return (
    <div className="flex h-screen bg-gray-900">
      {/* Sidebar */}
      <div
        className={`${
          sidebarCollapsed ? 'w-12' : 'w-64'
        } shrink-0 bg-gray-800 border-r border-gray-700 flex flex-col transition-all duration-200`}
      >
        {/* Sidebar header */}
        <div className="p-3 border-b border-gray-700 flex items-center justify-between">
          {!sidebarCollapsed && (
            <h2 className="text-sm font-semibold text-gray-300">Simulations</h2>
          )}
          <button
            onClick={() => setSidebarCollapsed(!sidebarCollapsed)}
            className="p-1 hover:bg-gray-700 rounded text-gray-400"
          >
            <ChevronRight
              size={16}
              className={`transform transition-transform ${sidebarCollapsed ? '' : 'rotate-180'}`}
            />
          </button>
        </div>

        {/* Simulation list */}
        {!sidebarCollapsed && (
          <div className="flex-1 overflow-y-auto">
            {isLoadingSimulations ? (
              <div className="p-4 text-center text-gray-500 text-sm">
                <div className="mb-2">Loading simulations...</div>
                <RefreshCw size={16} className="animate-spin mx-auto" />
              </div>
            ) : simulations.length === 0 ? (
              <div className="p-4 text-center text-gray-500 text-sm">
                <div className="mb-2">No simulations found</div>
                <div className="text-xs mb-2">
                  Run the data exporter to prepare simulation data
                </div>
                <div className="text-xs text-yellow-500">
                  Debug: Check browser console (F12) for fetch logs
                </div>
                {error && <div className="text-xs text-red-500 mt-2">Error: {error}</div>}
              </div>
            ) : (
              <div className="p-2 space-y-1">
                {simulations.map((sim) => (
                  <button
                    key={sim.id}
                    onClick={() => {
                      console.log('[Button Click] Simulation clicked:', sim.id)
                      alert(`Loading simulation: ${sim.name}`)
                      selectSimulation(sim)
                    }}
                    className={`w-full text-left p-2 rounded text-sm ${
                      selectedSimulation?.id === sim.id
                        ? 'bg-blue-600 text-white'
                        : 'text-gray-300 hover:bg-gray-700'
                    }`}
                  >
                    <div className="flex items-center gap-2">
                      <Folder size={14} />
                      <span className="truncate">{sim.name}</span>
                    </div>
                    <div className="text-xs opacity-60 ml-5">
                      {sim.totalFrames} frames • {sim.particleCount.toLocaleString()} particles
                    </div>
                  </button>
                ))}
              </div>
            )}
          </div>
        )}

        {/* Refresh button */}
        {!sidebarCollapsed && (
          <div className="p-2 border-t border-gray-700">
            <button
              onClick={loadSimulations}
              className="w-full flex items-center justify-center gap-2 p-2 text-sm text-gray-400 hover:bg-gray-700 rounded"
            >
              <RefreshCw size={14} />
              Refresh
            </button>
          </div>
        )}
      </div>

      {/* Main content */}
      <div className="flex-1 overflow-hidden relative">
        {/* Loading overlay when loading all frames */}
        {isLoading && loadingProgress.total > 0 && (
          <div className="absolute inset-0 bg-gray-900/90 flex items-center justify-center z-50">
            <div className="text-center">
              <div className="text-xl text-white mb-4">Loading Frames into Memory</div>
              <div className="w-64 bg-gray-700 rounded-full h-3 mb-2">
                <div 
                  className="bg-blue-500 h-3 rounded-full transition-all duration-100"
                  style={{ width: `${(loadingProgress.loaded / loadingProgress.total) * 100}%` }}
                />
              </div>
              <div className="text-gray-400">
                {loadingProgress.loaded} / {loadingProgress.total} frames
              </div>
              <div className="text-gray-500 text-sm mt-2">
                {((loadingProgress.loaded / loadingProgress.total) * 100).toFixed(0)}% complete
              </div>
            </div>
          </div>
        )}
        <Dashboard
          simulation={selectedSimulation}
          frames={frames}
          statistics={statistics}
          onLoadFrame={loadFrame}
          isLoading={isLoading}
          error={error}
        />
      </div>
    </div>
  )
}

/** Compute statistics for a frame */
function computeFrameStatistics(frame: ParsedFrame): FrameStatistics {
  let totalKinetic = 0
  let totalInternal = 0
  let totalMass = 0
  const momentum: [number, number, number] = [0, 0, 0]
  const com: [number, number, number] = [0, 0, 0]
  let minDensity = Infinity
  let maxDensity = -Infinity
  let minPressure = Infinity
  let maxPressure = -Infinity

  for (let i = 0; i < frame.particleCount; i++) {
    const m = frame.mass[i]
    const vx = frame.velocities[i * 3]
    const vy = frame.velocities[i * 3 + 1]
    const vz = frame.velocities[i * 3 + 2]
    const u = frame.energy[i]
    const rho = frame.density[i]
    const P = frame.pressure[i]

    // Kinetic energy
    const v2 = vx * vx + vy * vy + vz * vz
    totalKinetic += 0.5 * m * v2

    // Internal energy
    totalInternal += m * u

    // Momentum
    momentum[0] += m * vx
    momentum[1] += m * vy
    momentum[2] += m * vz

    // Center of mass
    com[0] += m * frame.positions[i * 3]
    com[1] += m * frame.positions[i * 3 + 1]
    com[2] += m * frame.positions[i * 3 + 2]

    totalMass += m

    // Ranges
    if (isFinite(rho)) {
      if (rho < minDensity) minDensity = rho
      if (rho > maxDensity) maxDensity = rho
    }
    if (isFinite(P)) {
      if (P < minPressure) minPressure = P
      if (P > maxPressure) maxPressure = P
    }
  }

  // Normalize center of mass
  if (totalMass > 0) {
    com[0] /= totalMass
    com[1] /= totalMass
    com[2] /= totalMass
  }

  return {
    frameIndex: frame.frameIndex,
    time: frame.time,
    totalMass,
    totalKineticEnergy: totalKinetic,
    totalInternalEnergy: totalInternal,
    totalEnergy: totalKinetic + totalInternal,
    momentum,
    centerOfMass: com,
    densityRange: [minDensity, maxDensity],
    pressureRange: [minPressure, maxPressure],
  }
}

export default VisualizationPage
