'use client'

import { useState, useCallback, useEffect, useMemo, useRef } from 'react'
import { Panel, PanelGroup, PanelResizeHandle } from 'react-resizable-panels'
import { ParticleViewer3DImperative } from '~/components/viewer/ParticleViewer3DImperative'
import { Projection3DInteractive, type Projection3DInteractiveHandle } from '~/components/viewer/Projection3DInteractive'
import { ShockDiagnosticsPanelImperative } from '~/components/viewer/ShockDiagnosticsPanelImperative'
import { EnergyChart, MomentumChart } from '~/components/charts/Charts'
import { TimescaleDiagnostics } from '~/components/charts/TimescaleDiagnostics'
import { PVDiagramImperative } from '~/components/charts/PVDiagramImperative'
import { PlaybackControls } from '~/components/controls/PlaybackControls'
import { VisualizationSettings } from '~/components/controls/VisualizationSettings'
import { COLOR_MAPS, type ParsedFrame, type FrameStatistics, type SimulationMetadata, type ColorMap, type IMBHPhysicsConfig } from '~/types/sph'
import { DEFAULT_IMBH_CONFIG, PRESET_SCENARIOS, presetToIMBHConfig, type SimulationPreset } from '~/utils/preset-loader'

/** Resize handle component for panels */
function ResizeHandle({ direction = 'horizontal' }: { direction?: 'horizontal' | 'vertical' }) {
  return (
    <PanelResizeHandle
      className={`
        ${direction === 'horizontal' ? 'w-1.5 cursor-col-resize' : 'h-1.5 cursor-row-resize'}
        bg-gray-700 hover:bg-cyan-500 active:bg-cyan-400 transition-colors
        flex items-center justify-center group
      `}
    >
      <div className={`
        ${direction === 'horizontal' ? 'w-0.5 h-8' : 'h-0.5 w-8'}
        bg-gray-500 group-hover:bg-cyan-300 rounded-full transition-colors
      `} />
    </PanelResizeHandle>
  )
}

/** Compute global color range statistics across all frames */
function computeGlobalColorStats(frames: Map<number, ParsedFrame>, colorField: string): [number, number] {
  if (frames.size === 0) return [0, 1]
  
  let globalMin = Infinity
  let globalMax = -Infinity
  let sampleCount = 0
  
  // Sample frames for efficiency (every 5th frame based on total)
  const frameArray = Array.from(frames.values())
  const step = Math.max(1, Math.floor(frameArray.length / 20)) // Sample ~20 frames max
  
  for (let frameIdx = 0; frameIdx < frameArray.length; frameIdx += step) {
    const frame = frameArray[frameIdx]
    if (!frame || frame.particleCount === 0) continue
    
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
        // Compute velocity magnitudes
        for (let i = 0; i < frame.particleCount; i += 50) {
          const vx = frame.velocities[i * 3]
          const vy = frame.velocities[i * 3 + 1]
          const vz = frame.velocities[i * 3 + 2]
          const vmag = Math.sqrt(vx * vx + vy * vy + vz * vz)
          if (isFinite(vmag) && vmag > 0) {
            if (vmag < globalMin) globalMin = vmag
            if (vmag > globalMax) globalMax = vmag
            sampleCount++
          }
        }
        continue // Skip fieldData processing for velocity
    }
    
    if (fieldData && fieldData.length > 0) {
      // Sample every 50th particle for speed
      for (let i = 0; i < fieldData.length; i += 50) {
        const v = fieldData[i]
        if (isFinite(v) && v > 0) {
          if (v < globalMin) globalMin = v
          if (v > globalMax) globalMax = v
          sampleCount++
        }
      }
    }
  }
  
  // Fallback to sensible defaults if no valid data found
  if (!isFinite(globalMin) || !isFinite(globalMax) || sampleCount === 0) {
    console.warn(`[computeGlobalColorStats] No valid data for ${colorField}, using defaults`)
    return [0.001, 1.0]
  }
  
  if (globalMin === globalMax) {
    globalMax = globalMin * 1.1 + 0.001
  }
  
  // Add small padding to range (5%)
  const range = globalMax - globalMin
  const paddedMin = Math.max(globalMin - range * 0.05, globalMin * 0.9)
  const paddedMax = globalMax + range * 0.05
  
  console.log(`[computeGlobalColorStats] ${colorField}: [${paddedMin.toExponential(2)}, ${paddedMax.toExponential(2)}] from ${sampleCount} samples`)
  
  return [paddedMin, paddedMax]
}

export interface DashboardProps {
  simulation: SimulationMetadata | null
  frames: Map<number, ParsedFrame>
  statistics: FrameStatistics[]
  onLoadFrame: (frameIndex: number) => Promise<void>
  isLoading: boolean
  error?: string | null
}

/** Main visualization dashboard */
export function Dashboard({
  simulation,
  frames,
  statistics,
  onLoadFrame,
  isLoading,
  error,
}: DashboardProps) {
  // Playback state
  const [currentFrameIndex, setCurrentFrameIndex] = useState(0)
  const [isPlaying, setIsPlaying] = useState(false)

  // Visualization settings
  const [colorField, setColorField] = useState('density')
  const [colorMapName, setColorMapName] = useState('cosmicDawn')
  const [pointSize, setPointSize] = useState(0.02)
  const [opacity, setOpacity] = useState(0.8)
  const [showAxes, setShowAxes] = useState(true)
  const [showBoundingBox, setShowBoundingBox] = useState(true)
  const [colorRange, setColorRange] = useState<[number, number]>([0, 0]) // 0,0 = auto
  const [useLogScale, setUseLogScale] = useState(false)
  
  // Per-projection color field settings (allows showing different fields in each 2D projection)
  const [projectionColorFields, setProjectionColorFields] = useState<{
    xy: string
    xz: string
    yz: string
  }>({
    xy: 'density',
    xz: 'pressure',
    yz: 'velocity',
  })
  const [useMultiColorField, setUseMultiColorField] = useState(false)

  // Layout state
  const [showProjections, setShowProjections] = useState(true)
  const [showCharts, setShowCharts] = useState(true)
  const [showShockDiagnostics, setShowShockDiagnostics] = useState(true)  // Default ON for equal importance
  const [showPVDiagram, setShowPVDiagram] = useState(true)  // P-V diagram for Oka comparison
  const [showTimescales, setShowTimescales] = useState(true)  // Equilibrium timescale diagnostics
  
  // Panel dimension tracking for responsive canvas sizing
  const [panelDimensions, setPanelDimensions] = useState({
    projection: { width: 300, height: 120 },
    shock: { width: 600, height: 250 },
  })
  const projectionPanelRef = useRef<HTMLDivElement>(null)
  const shockPanelRef = useRef<HTMLDivElement>(null)
  
  // Refs for interactive projection components (for reset camera)
  const projectionXYRef = useRef<Projection3DInteractiveHandle>(null)
  const projectionXZRef = useRef<Projection3DInteractiveHandle>(null)
  const projectionYZRef = useRef<Projection3DInteractiveHandle>(null)
  
  // Reset all projection cameras
  const resetAllProjectionCameras = useCallback(() => {
    projectionXYRef.current?.resetCamera()
    projectionXZRef.current?.resetCamera()
    projectionYZRef.current?.resetCamera()
  }, [])
  
  // Update panel dimensions on resize
  useEffect(() => {
    const updateDimensions = () => {
      if (projectionPanelRef.current) {
        const rect = projectionPanelRef.current.getBoundingClientRect()
        // Account for: padding (8px top/bottom), label (~24px with mb-1), gaps (8px = 2 * gap-1)
        // Total deduction: 8 + 24 + 8 = 40px
        const availableHeight = rect.height - 40
        const projectionHeight = Math.floor(availableHeight / 3)
        setPanelDimensions(prev => ({
          ...prev,
          projection: { 
            width: Math.floor(rect.width) - 16, 
            height: Math.max(60, projectionHeight)  // minimum 60px per projection
          }
        }))
      }
      if (shockPanelRef.current) {
        const rect = shockPanelRef.current.getBoundingClientRect()
        setPanelDimensions(prev => ({
          ...prev,
          shock: { width: Math.floor(rect.width) - 16, height: Math.floor(rect.height) - 60 }  // Extra space for description
        }))
      }
    }
    
    // Initial update
    updateDimensions()
    
    // Update on window resize
    window.addEventListener('resize', updateDimensions)
    
    // Use ResizeObserver for panel resizes
    const observer = new ResizeObserver(updateDimensions)
    if (projectionPanelRef.current) observer.observe(projectionPanelRef.current)
    if (shockPanelRef.current) observer.observe(shockPanelRef.current)
    
    return () => {
      window.removeEventListener('resize', updateDimensions)
      observer.disconnect()
    }
  }, [])
  
  // Performance monitoring
  const [currentFps, setCurrentFps] = useState(0)
  
  // IMBH visualization settings
  const [showBlackHole, setShowBlackHole] = useState(true)
  const [showTrajectory, setShowTrajectory] = useState(true)
  const [showRadii, setShowRadii] = useState(true)
  const [showGalacticMarkers, setShowGalacticMarkers] = useState(true)
  const [showLabels, setShowLabels] = useState(true)  // Text labels on markers
  const [showGalaxyDisk, setShowGalaxyDisk] = useState(true)  // Low-opacity Milky Way disk
  const [showSolarSystem, setShowSolarSystem] = useState(true)  // Observer direction indicator (i, PA, V_LSR)
  const [animateGalaxy, setAnimateGalaxy] = useState(false)  // Animate galaxy rotation
  const [galaxyAnimationSpeed, setGalaxyAnimationSpeed] = useState(1.0)  // Animation speed multiplier
  const [cameraMode, setCameraMode] = useState<'free' | 'earth'>('free')  // Camera view mode
  
  // Preset configuration state
  const [selectedPresetId, setSelectedPresetId] = useState<string | null>(null)
  const [loadedPresetConfig, setLoadedPresetConfig] = useState<IMBHPhysicsConfig | null>(null)
  
  // Load preset config when selection changes
  useEffect(() => {
    if (!selectedPresetId) {
      setLoadedPresetConfig(null)
      return
    }
    
    const preset = PRESET_SCENARIOS.find(p => p.id === selectedPresetId)
    if (!preset) return
    
    // Load preset from local file via fetch
    const loadPreset = async () => {
      try {
        // Try loading through API first
        const response = await fetch(`http://localhost:3001/api/preset?path=${encodeURIComponent(preset.path)}`)
        if (response.ok) {
          const presetData: SimulationPreset = await response.json()
          const config = presetToIMBHConfig(presetData)
          if (config) {
            console.log('[Dashboard] Loaded preset config:', preset.name, config)
            setLoadedPresetConfig(config)
          }
        } else {
          console.warn('[Dashboard] Failed to load preset from API, using default')
        }
      } catch (error) {
        console.warn('[Dashboard] Error loading preset:', error)
      }
    }
    
    loadPreset()
  }, [selectedPresetId])
  
  // IMBH physics configuration - priority: loadedPreset > simulation metadata > defaults
  const imbhPhysics = useMemo(() => {
    // Priority 1: Use loaded preset config if available
    if (loadedPresetConfig) {
      console.log('[Dashboard] Using IMBH physics from loaded preset')
      return loadedPresetConfig
    }
    
    // Priority 2: Use simulation metadata if available
    if (simulation?.imbhPhysics) {
      console.log('[Dashboard] Using IMBH physics from simulation config:', simulation.imbhPhysics)
      return simulation.imbhPhysics
    }
    
    // Priority 3: Default fallback
    console.log('[Dashboard] Using default IMBH physics')
    return DEFAULT_IMBH_CONFIG
  }, [loadedPresetConfig, simulation])

  // Refs for imperative viewer (avoids React re-renders)
  const framesRef = useRef<Map<number, ParsedFrame>>(frames)
  const frameIndexRef = useRef<number>(currentFrameIndex)
  const currentFrameRef = useRef<ParsedFrame | null>(null)
  
  // Imperative playback state (bypasses React)
  const isPlayingRef = useRef(false)
  const playbackSpeedRef = useRef(30) // FPS
  const animationFrameIdRef = useRef<number | null>(null)
  
  // Keep refs in sync (but don't trigger viewer re-renders)
  useEffect(() => {
    framesRef.current = frames
  }, [frames])

  // Keep current frame ref in sync
  useEffect(() => {
    currentFrameRef.current = frames.get(currentFrameIndex) ?? null
  }, [frames, currentFrameIndex])
  
  // Sync frame index ref (for non-playing state)
  useEffect(() => {
    if (!isPlayingRef.current) {
      frameIndexRef.current = currentFrameIndex
    }
  }, [currentFrameIndex])

  // Compute global color stats when frames or colorField changes
  const globalColorStats = useMemo(() => {
    return computeGlobalColorStats(frames, colorField)
  }, [frames, colorField])

  // Current frame data
  const currentFrame = frames.get(currentFrameIndex) || null

  // Color map with range
  const colorMap: ColorMap = useMemo(() => {
    const baseMap = COLOR_MAPS[colorMapName] || COLOR_MAPS.cosmicDawn
    
    // Use manual range if set, otherwise use global stats
    let min = colorRange[0]
    let max = colorRange[1]
    
    if (min === 0 && max === 0) {
      // Auto mode - use global stats computed across all frames
      [min, max] = globalColorStats
    }
    
    return {
      ...baseMap,
      min,
      max,
      logScale: useLogScale,
    }
  }, [colorMapName, colorRange, useLogScale, globalColorStats])

  // Load frame when index changes (fallback - frames should already be loaded)
  useEffect(() => {
    if (!frames.has(currentFrameIndex)) {
      onLoadFrame(currentFrameIndex)
    }
  }, [currentFrameIndex, frames, onLoadFrame])

  // Cache status
  const preloadedCount = frames.size
  const totalFrameCount = simulation?.totalFrames || 0
  const allFramesLoaded = preloadedCount === totalFrameCount

  // Check if a frame is ready (for smooth playback)
  const isFrameReady = useCallback((frameIndex: number) => {
    return frames.has(frameIndex)
  }, [frames])

  // Frame change handler
  const handleFrameChange = useCallback((frame: number) => {
    const clampedFrame = Math.max(0, Math.min(frame, (simulation?.totalFrames || 1) - 1))
    frameIndexRef.current = clampedFrame
    // Only update React state if not playing (to update UI controls)
    if (!isPlayingRef.current) {
      setCurrentFrameIndex(clampedFrame)
    }
  }, [simulation?.totalFrames])

  // Imperative playback system - runs outside React render cycle
  const startImperativePlayback = useCallback(() => {
    if (animationFrameIdRef.current !== null) return // Already playing
    
    isPlayingRef.current = true
    setIsPlaying(true)
    
    let lastTime = 0
    const totalFrames = simulation?.totalFrames || 1
    
    const tick = (timestamp: number) => {
      if (!isPlayingRef.current) return
      
      if (!lastTime) lastTime = timestamp
      const elapsed = timestamp - lastTime
      const frameInterval = 1000 / playbackSpeedRef.current
      
      if (elapsed >= frameInterval) {
        const currentIdx = frameIndexRef.current
        const nextFrame = currentIdx + 1
        
        if (nextFrame >= totalFrames) {
          // Reached end - auto-replay from beginning
          frameIndexRef.current = 0
          setCurrentFrameIndex(0)
          lastTime = timestamp
        } else {
          // Only advance if frame is loaded
          if (framesRef.current.has(nextFrame)) {
            frameIndexRef.current = nextFrame
            // Sync React state for 2D projections and charts
            setCurrentFrameIndex(nextFrame)
            lastTime = timestamp
          }
        }
      }
      
      animationFrameIdRef.current = requestAnimationFrame(tick)
    }
    
    animationFrameIdRef.current = requestAnimationFrame(tick)
  }, [simulation?.totalFrames])

  const stopImperativePlayback = useCallback(() => {
    isPlayingRef.current = false
    setIsPlaying(false)
    if (animationFrameIdRef.current !== null) {
      cancelAnimationFrame(animationFrameIdRef.current)
      animationFrameIdRef.current = null
    }
    // Sync React state with current ref value
    setCurrentFrameIndex(frameIndexRef.current)
  }, [])

  const handlePlayPauseChange = useCallback((playing: boolean) => {
    if (playing) {
      startImperativePlayback()
    } else {
      stopImperativePlayback()
    }
  }, [startImperativePlayback, stopImperativePlayback])

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      if (animationFrameIdRef.current !== null) {
        cancelAnimationFrame(animationFrameIdRef.current)
      }
    }
  }, [])

  if (!simulation) {
    return (
      <div className="flex items-center justify-center h-full bg-gray-900 text-gray-400">
        <div className="text-center">
          <div className="text-2xl mb-2">📊</div>
          <div>No simulation loaded</div>
          <div className="text-sm mt-1">Select a simulation from the sidebar</div>
        </div>
      </div>
    )
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-full bg-gray-900 text-red-400">
        <div className="text-center">
          <div className="text-2xl mb-2">❌</div>
          <div>Error loading simulation</div>
          <div className="text-sm mt-1">{error}</div>
        </div>
      </div>
    )
  }

  const currentTime = currentFrame?.time || 0

  return (
    <div className="flex flex-col h-full bg-gray-900">
      {/* Header */}
      <div className="flex items-center justify-between px-4 py-2 bg-gray-800 border-b border-gray-700">
        <div>
          <h1 className="text-lg font-semibold text-white">{simulation.name}</h1>
          <div className="text-xs text-gray-400">
            {simulation.method} • {simulation.kernel} • {simulation.particleCount.toLocaleString()} particles
          </div>
        </div>
        <div className="flex items-center gap-3">
          {/* Cache status */}
          <div className="text-xs">
            <span className={allFramesLoaded ? 'text-green-400' : 'text-yellow-400'}>
              {allFramesLoaded ? '✓ All frames in memory' : `${preloadedCount}/${totalFrameCount} frames`}
            </span>
          </div>
          <button
            onClick={() => setShowProjections(!showProjections)}
            className={`px-2 py-1 text-xs rounded ${showProjections ? 'bg-blue-600 text-white' : 'bg-gray-700 text-gray-400'}`}
          >
            2D Views
          </button>
          <button
            onClick={() => setShowShockDiagnostics(!showShockDiagnostics)}
            className={`px-2 py-1 text-xs rounded ${showShockDiagnostics ? 'bg-orange-600 text-white' : 'bg-gray-700 text-gray-400'}`}
          >
            Shock
          </button>
          <button
            onClick={() => setShowCharts(!showCharts)}
            className={`px-2 py-1 text-xs rounded ${showCharts ? 'bg-blue-600 text-white' : 'bg-gray-700 text-gray-400'}`}
          >
            Charts
          </button>
          <button
            onClick={() => setShowPVDiagram(!showPVDiagram)}
            className={`px-2 py-1 text-xs rounded ${showPVDiagram ? 'bg-purple-600 text-white' : 'bg-gray-700 text-gray-400'}`}
            title="Position-Velocity diagram for radio observation comparison"
          >
            P-V
          </button>
          <button
            onClick={() => setShowTimescales(!showTimescales)}
            className={`px-2 py-1 text-xs rounded ${showTimescales ? 'bg-green-600 text-white' : 'bg-gray-700 text-gray-400'}`}
            title="Timescale diagnostics for equilibrium chemistry validation"
          >
            Timescales
          </button>
        </div>
      </div>

      {/* Main content */}
      <div className="flex-1 flex overflow-hidden min-h-0">
        {/* Left sidebar - Settings */}
        <div className="w-56 shrink-0 overflow-y-auto border-r border-gray-700 p-2">
          <VisualizationSettings
            colorField={colorField}
            onColorFieldChange={setColorField}
            colorMapName={colorMapName}
            onColorMapChange={setColorMapName}
            pointSize={pointSize}
            onPointSizeChange={setPointSize}
            opacity={opacity}
            onOpacityChange={setOpacity}
            showAxes={showAxes}
            onShowAxesChange={setShowAxes}
            showBoundingBox={showBoundingBox}
            onShowBoundingBoxChange={setShowBoundingBox}
            colorRange={colorRange}
            onColorRangeChange={setColorRange}
            useLogScale={useLogScale}
            onLogScaleChange={setUseLogScale}
          />
          
          {/* Multi-Field 2D Projections Settings */}
          <div className="mt-2 bg-gray-800 rounded p-3">
            <div className="flex items-center justify-between mb-2">
              <h3 className="text-sm font-medium text-gray-300">2D Projections</h3>
              <label className="flex items-center gap-1 cursor-pointer">
                <input
                  type="checkbox"
                  checked={useMultiColorField}
                  onChange={(e) => setUseMultiColorField(e.target.checked)}
                  className="rounded"
                />
                <span className="text-xs text-gray-400">Multi-field</span>
              </label>
            </div>
            
            {useMultiColorField && (
              <div className="space-y-2 text-xs">
                {(['xy', 'xz', 'yz'] as const).map(proj => (
                  <div key={proj} className="flex items-center gap-2">
                    <span className="text-gray-400 w-8 uppercase font-mono">{proj}:</span>
                    <select
                      value={projectionColorFields[proj]}
                      onChange={(e) => setProjectionColorFields(prev => ({
                        ...prev,
                        [proj]: e.target.value
                      }))}
                      className="flex-1 bg-gray-700 text-gray-200 rounded px-1 py-0.5 border border-gray-600"
                    >
                      <option value="density">Density</option>
                      <option value="pressure">Pressure</option>
                      <option value="velocity">Velocity</option>
                      <option value="energy">Energy</option>
                      <option value="machNumber">Mach #</option>
                    </select>
                  </div>
                ))}
              </div>
            )}
          </div>
          
          {/* IMBH Visualization Settings */}
          <div className="mt-2 bg-gray-800 rounded p-3">
              <h3 className="text-sm font-medium text-gray-300 mb-2">IMBH Physics</h3>
              
              {/* Preset Selector */}
              <div className="mb-3">
                <label className="block text-xs text-gray-400 mb-1">Preset Config</label>
                <select
                  value={selectedPresetId ?? ''}
                  onChange={(e) => setSelectedPresetId(e.target.value || null)}
                  className="w-full bg-gray-700 text-gray-200 text-xs rounded px-2 py-1 border border-gray-600"
                >
                  <option value="">Default (from simulation)</option>
                  {PRESET_SCENARIOS.map(preset => (
                    <option key={preset.id} value={preset.id}>
                      {preset.name}
                    </option>
                  ))}
                </select>
              </div>
              
              <div className="space-y-2 text-xs">
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={showBlackHole}
                    onChange={(e) => setShowBlackHole(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Black Hole</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={showTrajectory}
                    onChange={(e) => setShowTrajectory(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Trajectory</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={showRadii}
                    onChange={(e) => setShowRadii(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Tidal/Hill Radii</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={showGalacticMarkers}
                    onChange={(e) => setShowGalacticMarkers(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Galactic Markers</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer" title="Show text labels on all markers and arrows">
                  <input
                    type="checkbox"
                    checked={showLabels}
                    onChange={(e) => setShowLabels(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Labels</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer" title="Show schematic Milky Way disk with Sun at 8 kpc from GC">
                  <input
                    type="checkbox"
                    checked={showGalaxyDisk}
                    onChange={(e) => setShowGalaxyDisk(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Galaxy Disk</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer" title="Show observer direction with inclination (70°), position angle (41.6°), and V_LSR indicator">
                  <input
                    type="checkbox"
                    checked={showSolarSystem}
                    onChange={(e) => setShowSolarSystem(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-gray-300">Show Observer Geometry</span>
                </label>
                {/* Galaxy Animation Controls */}
                {showGalaxyDisk && (
                  <div className="mt-2 pt-2 border-t border-gray-600">
                    <label className="flex items-center gap-2 cursor-pointer" title="Animate galactic rotation (Sun orbits GC in ~220 Myr)">
                      <input
                        type="checkbox"
                        checked={animateGalaxy}
                        onChange={(e) => setAnimateGalaxy(e.target.checked)}
                        className="rounded"
                      />
                      <span className="text-gray-300">Animate Galaxy</span>
                    </label>
                    {animateGalaxy && (
                      <div className="mt-2 ml-4">
                        <div className="text-gray-400 text-[10px] mb-1">Speed: {galaxyAnimationSpeed.toFixed(1)}× (1× = 220 Myr period)</div>
                        <input
                          type="range"
                          min="0.1"
                          max="10"
                          step="0.1"
                          value={galaxyAnimationSpeed}
                          onChange={(e) => setGalaxyAnimationSpeed(parseFloat(e.target.value))}
                          className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer"
                        />
                        <div className="flex justify-between text-[9px] text-gray-500 mt-1">
                          <span>Slow</span>
                          <span>~{(220 / galaxyAnimationSpeed).toFixed(0)} Myr/rev</span>
                          <span>Fast</span>
                        </div>
                      </div>
                    )}
                  </div>
                )}
                {/* Camera Mode Toggle */}
                <div className="mt-3 pt-2 border-t border-gray-700">
                  <div className="text-gray-400 text-[10px] mb-2">Camera View</div>
                  <div className="flex gap-1">
                    <button
                      onClick={() => setCameraMode('free')}
                      className={`flex-1 px-2 py-1 text-xs rounded ${
                        cameraMode === 'free' 
                          ? 'bg-blue-600 text-white' 
                          : 'bg-gray-700 text-gray-300 hover:bg-gray-600'
                      }`}
                      title="Free orbit camera controls"
                    >
                      🎥 Free
                    </button>
                    <button
                      onClick={() => setCameraMode('earth')}
                      className={`flex-1 px-2 py-1 text-xs rounded ${
                        cameraMode === 'earth' 
                          ? 'bg-green-600 text-white' 
                          : 'bg-gray-700 text-gray-300 hover:bg-gray-600'
                      }`}
                      title="View from Earth (observer's perspective)"
                    >
                      🌍 Earth View
                    </button>
                  </div>
                </div>
                <div className="mt-2 pt-2 border-t border-gray-700">
                  <div className="text-gray-400 text-[10px]">
                    <div>BH: {imbhPhysics.bhMass}×10⁵ M☉</div>
                    <div>Tidal: {imbhPhysics.tidalRadius.toFixed(2)} pc</div>
                    <div>Impact: {imbhPhysics.impactParameter.toFixed(2)} pc</div>
                  </div>
                </div>
              </div>
            </div>
          
          {/* Statistics panel */}
          {currentFrame && (
            <div className="mt-2 bg-gray-800 rounded p-3">
              <h3 className="text-sm font-medium text-gray-300 mb-2">Frame Statistics</h3>
              <div className="space-y-1 text-xs text-gray-400">
                <div className="flex justify-between">
                  <span>Particles:</span>
                  <span className="text-white">{currentFrame.particleCount.toLocaleString()}</span>
                </div>
                <div className="flex justify-between">
                  <span>Time:</span>
                  <span className="text-white font-mono">
                    {(currentTime * imbhPhysics.timeUnit).toFixed(3)} Myr
                    <span className="text-gray-500 text-[10px] ml-1">({currentTime.toFixed(3)} code)</span>
                  </span>
                </div>
                {statistics[currentFrameIndex] && (
                  <>
                    <div className="flex justify-between">
                      <span>Total Energy:</span>
                      <span className="text-white font-mono">
                        {statistics[currentFrameIndex].totalEnergy.toExponential(3)}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span>Max Density:</span>
                      <span className="text-white font-mono">
                        {statistics[currentFrameIndex].densityRange[1].toExponential(3)}
                      </span>
                    </div>
                  </>
                )}
              </div>
            </div>
          )}
        </div>

        {/* Main viewer area - Optimized 2-row layout */}
        <div className="flex-1 overflow-y-auto overflow-x-hidden">
          <PanelGroup direction="vertical" className="min-h-[900px]">
            {/* Row 1: 3D Viewer (full width, dominant) */}
            <Panel defaultSize={55} minSize={35}>
              <div className="h-full min-h-[400px] relative bg-gray-900">
                {isLoading && !currentFrame ? (
                  <div className="flex items-center justify-center h-full text-gray-400">
                    <div className="text-center">
                      <div className="animate-spin text-2xl mb-2">⏳</div>
                      <div>Loading frame {currentFrameIndex}...</div>
                    </div>
                  </div>
                ) : (
                  <ParticleViewer3DImperative
                    framesRef={framesRef}
                    frameIndexRef={frameIndexRef}
                    colorField={colorField}
                    colorMapName={colorMapName}
                    pointSize={pointSize * 100}
                    opacity={opacity}
                    logScale={useLogScale}
                    showAxes={showAxes}
                    showBoundingBox={showBoundingBox}
                    boundingBox={simulation.boundingBox}
                    className="h-full w-full"
                    onFpsUpdate={setCurrentFps}
                    globalColorRange={globalColorStats}
                    imbhPhysics={imbhPhysics}
                    showBlackHole={showBlackHole}
                    showTrajectory={showTrajectory}
                    showRadii={showRadii}
                    showGalacticMarkers={showGalacticMarkers}
                    showLabels={showLabels}
                    cameraMode={cameraMode}
                    galacticConfig={{
                      distanceToGC: 60,
                      galacticLongitude: -0.398,
                      galacticLatitude: -0.224,
                      inclination: simulation?.imbhPhysics?.inclination ?? 70,
                      positionAngle: simulation?.imbhPhysics?.positionAngle ?? 41.6,
                      lsrVelocity: simulation?.imbhPhysics?.lsrVelocity ?? -120,
                      showGalaxyDisk,
                      showSolarSystem,
                      galaxyRotationSpeed: animateGalaxy ? 0.1 * galaxyAnimationSpeed : 0,
                    }}
                  />
                )}
                {/* FPS overlay */}
                <div className="absolute bottom-2 right-2 text-green-400 text-xs font-mono bg-black/50 px-2 py-1 rounded">
                  FPS: {currentFps}
                </div>
                {/* Panel label */}
                <div className="absolute top-2 left-2 text-cyan-400 text-xs font-bold bg-black/60 px-2 py-1 rounded">
                  3D VIEW
                </div>
              </div>
            </Panel>

            {/* Row 2: Analysis panels (P-V, 2D Projections, Diagnostics) */}
            {(showProjections || showPVDiagram || showShockDiagnostics || showCharts || showTimescales) && (
              <>
                <ResizeHandle direction="vertical" />
                <Panel defaultSize={45} minSize={30}>
                  <PanelGroup direction="horizontal" className="h-full min-h-[350px]">
                    {/* P-V Diagram Panel - Square aspect ratio */}
                    {showPVDiagram && (
                      <Panel defaultSize={30} minSize={20}>
                        <div className="h-full bg-gray-900 p-2 overflow-hidden flex flex-col">
                          <div className="shrink-0 mb-1">
                            <div className="text-purple-400 text-xs font-bold">P-V DIAGRAM</div>
                            <div className="text-gray-500 text-[10px]">
                              Drag to rotate • i={simulation?.imbhPhysics?.inclination ?? 70}° PA={simulation?.imbhPhysics?.positionAngle ?? 41.6}°
                            </div>
                          </div>
                          <div className="flex-1 min-h-0 flex items-center justify-center">
                            <div className="w-full h-full max-w-[400px] aspect-square">
                              <PVDiagramImperative
                                framesRef={framesRef}
                                frameIndexRef={frameIndexRef}
                                initialObserver={{
                                  inclination: simulation?.imbhPhysics?.inclination ?? 70,
                                  positionAngle: simulation?.imbhPhysics?.positionAngle ?? 41.6,
                                  vLSR: simulation?.imbhPhysics?.lsrVelocity ?? -120,
                                }}
                                posRange={simulation?.pvDiagram?.positionRange ?? [-5, 5]}
                                velRange={simulation?.pvDiagram?.velocityRange ?? [-180, -60]}
                                className="h-full w-full"
                              />
                            </div>
                          </div>
                        </div>
                      </Panel>
                    )}

                    {/* 2D Projections Panel */}
                    {showProjections && (
                      <>
                        {showPVDiagram && <ResizeHandle direction="horizontal" />}
                        <Panel defaultSize={25} minSize={15}>
                          <div ref={projectionPanelRef} className="h-full bg-gray-900 p-2 flex flex-col overflow-hidden">
                            <div className="flex items-center justify-between mb-1 shrink-0">
                              <div className="text-cyan-400 text-xs font-bold">2D PROJECTIONS</div>
                              <button
                                onClick={resetAllProjectionCameras}
                                className="px-1 py-0.5 text-[9px] bg-gray-700 hover:bg-gray-600 text-gray-300 rounded"
                                title="Reset cameras"
                              >
                                ⟲
                              </button>
                            </div>
                            <div className="flex-1 flex flex-col gap-1 min-h-0 overflow-hidden">
                              <div className="flex-1 min-h-0">
                                <Projection3DInteractive
                                  ref={projectionXYRef}
                                  framesRef={framesRef}
                                  frameIndexRef={frameIndexRef}
                                  projection="xy"
                                  colorField={useMultiColorField ? projectionColorFields.xy : colorField}
                                  colorMap={colorMap}
                                  width={panelDimensions.projection.width}
                                  height={panelDimensions.projection.height}
                                  showShockSampling={showShockDiagnostics}
                                  shockSamplingParams={{ columnRadius: 0.15, sliceThickness: 0.15 }}
                                  globalColorRange={globalColorStats}
                                  logScale={useLogScale}
                                  particleSize={pointSize * 100}
                                />
                              </div>
                              <div className="flex-1 min-h-0">
                                <Projection3DInteractive
                                  ref={projectionXZRef}
                                  framesRef={framesRef}
                                  frameIndexRef={frameIndexRef}
                                  projection="xz"
                                  colorField={useMultiColorField ? projectionColorFields.xz : colorField}
                                  colorMap={colorMap}
                                  width={panelDimensions.projection.width}
                                  height={panelDimensions.projection.height}
                                  showShockSampling={showShockDiagnostics}
                                  shockSamplingParams={{ columnRadius: 0.15, sliceThickness: 0.15 }}
                                  globalColorRange={globalColorStats}
                                  logScale={useLogScale}
                                  particleSize={pointSize * 100}
                                />
                              </div>
                              <div className="flex-1 min-h-0">
                                <Projection3DInteractive
                                  ref={projectionYZRef}
                                  framesRef={framesRef}
                                  frameIndexRef={frameIndexRef}
                                  projection="yz"
                                  colorField={useMultiColorField ? projectionColorFields.yz : colorField}
                                  colorMap={colorMap}
                                  width={panelDimensions.projection.width}
                                  height={panelDimensions.projection.height}
                                  showShockSampling={showShockDiagnostics}
                                  shockSamplingParams={{ columnRadius: 0.15, sliceThickness: 0.15 }}
                                  globalColorRange={globalColorStats}
                                  logScale={useLogScale}
                                  particleSize={pointSize * 100}
                                />
                              </div>
                            </div>
                          </div>
                        </Panel>
                      </>
                    )}

                    {/* Shock Diagnostics Panel */}
                    {showShockDiagnostics && (
                      <>
                        {(showPVDiagram || showProjections) && <ResizeHandle direction="horizontal" />}
                        <Panel defaultSize={25} minSize={15}>
                          <div ref={shockPanelRef} className="h-full bg-gray-900 p-2 overflow-hidden flex flex-col">
                            <div className="shrink-0 mb-1">
                              <div className="text-orange-400 text-xs font-bold">SHOCK PROFILES</div>
                              <div className="text-gray-500 text-[9px]">
                                <span className="text-red-400">●</span> Z <span className="text-teal-400">━</span> X
                              </div>
                            </div>
                            <div className="flex-1 min-h-0">
                              <ShockDiagnosticsPanelImperative
                                framesRef={framesRef}
                                frameIndexRef={frameIndexRef}
                                initialFrame={frames.get(0) ?? null}
                                width={panelDimensions.shock.width}
                                height={panelDimensions.shock.height}
                                className="w-full h-full"
                              />
                            </div>
                          </div>
                        </Panel>
                      </>
                    )}

                    {/* Charts + Timescales Panel */}
                    {(showCharts || showTimescales) && (
                      <>
                        {(showPVDiagram || showProjections || showShockDiagnostics) && <ResizeHandle direction="horizontal" />}
                        <Panel defaultSize={20} minSize={15}>
                          <div className="h-full bg-gray-900 p-2 flex flex-col gap-1 overflow-hidden">
                            {/* Charts Section */}
                            {showCharts && statistics.length > 0 && (
                              <div className={`flex flex-col gap-1 ${showTimescales ? 'flex-1' : 'flex-1'} min-h-0`}>
                                <div className="text-cyan-400 text-[10px] font-bold shrink-0">ENERGY</div>
                                <div className="flex-1 min-h-0">
                                  <EnergyChart
                                    statistics={statistics}
                                    currentFrame={currentFrameIndex}
                                    className="h-full"
                                  />
                                </div>
                              </div>
                            )}

                            {/* Timescale Diagnostics Section */}
                            {showTimescales && (
                              <div className="flex-1 min-h-0 flex flex-col">
                                <div className="text-green-400 text-[10px] font-bold shrink-0">EQUILIBRIUM</div>
                                <div className="flex-1 min-h-0">
                                  <TimescaleDiagnostics
                                    frameRef={currentFrameRef}
                                    bhPosition={imbhPhysics.bhPosition}
                                    className="h-full"
                                  />
                                </div>
                              </div>
                            )}
                          </div>
                        </Panel>
                      </>
                    )}
                  </PanelGroup>
                </Panel>
              </>
            )}
          </PanelGroup>
        </div>
      </div>

      {/* Playback controls */}
      <div className="shrink-0 border-t border-gray-700">
        <PlaybackControls
          currentFrame={currentFrameIndex}
          totalFrames={simulation.totalFrames}
          time={currentTime}
          onFrameChange={handleFrameChange}
          onPlayPauseChange={handlePlayPauseChange}
          isFrameReady={(frame) => frames.has(frame)}
          imperativeMode={true}
          frameIndexRef={frameIndexRef}
          playbackSpeedRef={playbackSpeedRef}
        />
      </div>
    </div>
  )
}

export default Dashboard
