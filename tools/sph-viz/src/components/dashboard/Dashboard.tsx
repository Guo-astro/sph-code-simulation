'use client'

import { useState, useCallback, useEffect, useMemo, useRef } from 'react'
import { ParticleViewer3D } from '~/components/viewer/ParticleViewer3D'
import { ParticleViewer3DImperative } from '~/components/viewer/ParticleViewer3DImperative'
import { Projection2D } from '~/components/viewer/Projection2D'
import { EnergyChart, MomentumChart } from '~/components/charts/Charts'
import { PlaybackControls } from '~/components/controls/PlaybackControls'
import { VisualizationSettings } from '~/components/controls/VisualizationSettings'
import { COLOR_MAPS, type ParsedFrame, type FrameStatistics, type SimulationMetadata, type ColorMap } from '~/types/sph'

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
  const [colorMapName, setColorMapName] = useState('viridis')
  const [pointSize, setPointSize] = useState(0.02)
  const [opacity, setOpacity] = useState(0.8)
  const [showAxes, setShowAxes] = useState(true)
  const [showBoundingBox, setShowBoundingBox] = useState(true)
  const [colorRange, setColorRange] = useState<[number, number]>([0, 0]) // 0,0 = auto
  const [useLogScale, setUseLogScale] = useState(false)

  // Layout state
  const [showProjections, setShowProjections] = useState(true)
  const [showCharts, setShowCharts] = useState(true)
  
  // High-performance mode toggle
  const [useImperativeMode, setUseImperativeMode] = useState(true)
  const [currentFps, setCurrentFps] = useState(0)
  
  // IMBH visualization settings
  const [showBlackHole, setShowBlackHole] = useState(true)
  const [showTrajectory, setShowTrajectory] = useState(true)
  const [showRadii, setShowRadii] = useState(true)
  const [showGalacticMarkers, setShowGalacticMarkers] = useState(true)
  const [showGalaxyDisk, setShowGalaxyDisk] = useState(true)  // Low-opacity Milky Way disk
  const [showSolarSystem, setShowSolarSystem] = useState(true)  // Observer direction indicator (i, PA, V_LSR)
  const [animateGalaxy, setAnimateGalaxy] = useState(false)  // Animate galaxy rotation
  const [galaxyAnimationSpeed, setGalaxyAnimationSpeed] = useState(1.0)  // Animation speed multiplier
  const [cameraMode, setCameraMode] = useState<'free' | 'earth'>('free')  // Camera view mode
  
  // IMBH physics configuration - use from simulation metadata or defaults
  const imbhPhysics = useMemo(() => {
    // Use simulation metadata if available
    if (simulation?.imbhPhysics) {
      console.log('[Dashboard] Using IMBH physics from simulation config:', simulation.imbhPhysics)
      return simulation.imbhPhysics
    }
    
    // Default fallback for IMBH simulations (Oka et al. 2017 / CAT_OKA A_61k preset values)
    console.log('[Dashboard] Using default IMBH physics (CAT_OKA/A_61k/oka.json preset values)')
    return {
      enabled: true,  // Enable by default for visualizing galactic markers
      bhPosition: [0, 0, 0] as [number, number, number],
      bhMass: 100,  // 10^5 M_sun in code units (1000 M_sun)
      cloudInitialPosition: [20.0, -5.17, 0] as [number, number, number],  // From preset
      cloudInitialVelocity: [-10.18, 5.05, 0] as [number, number, number], // From preset (km/s)
      cloudMass: 1,  // 1000 M_sun in code units
      cloudRadius: 1.13,  // From preset (pc)
      tidalRadius: 5.24,  // From preset (pc)
      impactParameter: 5.17,  // From preset (pc)
      pericentre: 2.217,  // From preset (pc)
      eccentricity: 1.4504,  // From preset (hyperbolic)
      timeUnit: 0.978,  // Myr
    }
  }, [simulation])

  // Refs for imperative viewer (avoids React re-renders)
  const framesRef = useRef<Map<number, ParsedFrame>>(frames)
  const frameIndexRef = useRef<number>(currentFrameIndex)
  
  // Imperative playback state (bypasses React)
  const isPlayingRef = useRef(false)
  const playbackSpeedRef = useRef(30) // FPS
  const animationFrameIdRef = useRef<number | null>(null)
  
  // Keep refs in sync (but don't trigger viewer re-renders)
  useEffect(() => {
    framesRef.current = frames
  }, [frames])
  
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
    const baseMap = COLOR_MAPS[colorMapName] || COLOR_MAPS.viridis
    
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
    if (useImperativeMode) {
      if (playing) {
        startImperativePlayback()
      } else {
        stopImperativePlayback()
      }
    } else {
      setIsPlaying(playing)
    }
  }, [useImperativeMode, startImperativePlayback, stopImperativePlayback])

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
          {/* Performance mode toggle */}
          <button
            onClick={() => setUseImperativeMode(!useImperativeMode)}
            className={`px-2 py-1 text-xs rounded ${useImperativeMode ? 'bg-green-600 text-white' : 'bg-gray-700 text-gray-400'}`}
            title={useImperativeMode ? 'High-performance mode (120+ FPS)' : 'Standard React mode'}
          >
            {useImperativeMode ? '🚀 Fast' : '⚛️ React'}
          </button>
          <button
            onClick={() => setShowProjections(!showProjections)}
            className={`px-2 py-1 text-xs rounded ${showProjections ? 'bg-blue-600 text-white' : 'bg-gray-700 text-gray-400'}`}
          >
            2D Views
          </button>
          <button
            onClick={() => setShowCharts(!showCharts)}
            className={`px-2 py-1 text-xs rounded ${showCharts ? 'bg-blue-600 text-white' : 'bg-gray-700 text-gray-400'}`}
          >
            Charts
          </button>
        </div>
      </div>

      {/* Main content */}
      <div className="flex-1 flex overflow-hidden">
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
          
          {/* IMBH Visualization Settings */}
          {useImperativeMode && (
            <div className="mt-2 bg-gray-800 rounded p-3">
              <h3 className="text-sm font-medium text-gray-300 mb-2">IMBH Physics</h3>
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
          )}
          
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

        {/* Main viewer area */}
        <div className="flex-1 flex flex-col overflow-hidden">
          <div className="flex-1 flex overflow-hidden">
            {/* 3D Viewer */}
            <div className={`flex-1 ${showProjections ? '' : 'w-full'} relative`}>
              {isLoading && !currentFrame ? (
                <div className="flex items-center justify-center h-full bg-gray-900 text-gray-400">
                  <div className="text-center">
                    <div className="animate-spin text-2xl mb-2">⏳</div>
                    <div>Loading frame {currentFrameIndex}...</div>
                  </div>
                </div>
              ) : useImperativeMode ? (
                <ParticleViewer3DImperative
                  framesRef={framesRef}
                  frameIndexRef={frameIndexRef}
                  colorField={colorField}
                  colorMapName={colorMapName}
                  pointSize={pointSize * 100} // Scale for imperative viewer
                  opacity={opacity}
                  logScale={useLogScale}
                  showAxes={showAxes}
                  showBoundingBox={showBoundingBox}
                  boundingBox={simulation.boundingBox}
                  className="h-full"
                  onFpsUpdate={setCurrentFps}
                  globalColorRange={globalColorStats}
                  imbhPhysics={imbhPhysics}
                  showBlackHole={showBlackHole}
                  showTrajectory={showTrajectory}
                  showRadii={showRadii}
                  showGalacticMarkers={showGalacticMarkers}
                  cameraMode={cameraMode}
                  galacticConfig={{
                    distanceToGC: 60,  // ~60 pc from Galactic Center
                    galacticLongitude: -0.398,
                    galacticLatitude: -0.224,
                    inclination: simulation?.imbhPhysics?.inclination ?? 70,
                    positionAngle: simulation?.imbhPhysics?.positionAngle ?? 41.6,
                    lsrVelocity: simulation?.imbhPhysics?.lsrVelocity ?? -120,
                    showGalaxyDisk,
                    showSolarSystem,
                    // Physical rotation: Sun orbits GC at V_circ = 220 km/s, period ~220 Myr
                    // For visualization: 2π rad / (220 Myr) ≈ 2.86e-17 rad/s (too slow to see)
                    // We use a visual speed: 0.1 rad/s = ~63 seconds per revolution
                    // Speed multiplier adjusts this: higher = faster animation
                    galaxyRotationSpeed: animateGalaxy ? 0.1 * galaxyAnimationSpeed : 0,
                  }}
                />
              ) : (
                <ParticleViewer3D
                  frame={currentFrame}
                  colorField={colorField}
                  colorMap={colorMap}
                  pointSize={pointSize}
                  opacity={opacity}
                  showAxes={showAxes}
                  showBoundingBox={showBoundingBox}
                  boundingBox={simulation.boundingBox}
                  className="h-full"
                />
              )}
              
              {/* FPS overlay for imperative mode */}
              {useImperativeMode && (
                <div className="absolute bottom-2 right-2 text-green-400 text-xs font-mono bg-black/50 px-2 py-1 rounded">
                  FPS: {currentFps}
                </div>
              )}
            </div>

            {/* 2D Projections */}
            {showProjections && (
              <div className="w-64 shrink-0 flex flex-col gap-1 p-1 overflow-y-auto border-l border-gray-700">
                <Projection2D
                  frame={currentFrame}
                  projection="xy"
                  colorField={colorField}
                  colorMap={colorMap}
                  width={240}
                  height={150}
                />
                <Projection2D
                  frame={currentFrame}
                  projection="xz"
                  colorField={colorField}
                  colorMap={colorMap}
                  width={240}
                  height={150}
                />
                <Projection2D
                  frame={currentFrame}
                  projection="yz"
                  colorField={colorField}
                  colorMap={colorMap}
                  width={240}
                  height={150}
                />
              </div>
            )}
          </div>

          {/* Charts */}
          {showCharts && statistics.length > 0 && (
            <div className="h-40 shrink-0 border-t border-gray-700 flex gap-2 p-1 overflow-x-auto">
              <EnergyChart
                statistics={statistics}
                currentFrame={currentFrameIndex}
                className="flex-1 min-w-[250px]"
              />
              <MomentumChart
                statistics={statistics}
                currentFrame={currentFrameIndex}
                className="flex-1 min-w-[250px]"
              />
            </div>
          )}
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
          imperativeMode={useImperativeMode}
          frameIndexRef={frameIndexRef}
          playbackSpeedRef={playbackSpeedRef}
        />
      </div>
    </div>
  )
}

export default Dashboard
