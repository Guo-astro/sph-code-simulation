'use client'

import { useState, useEffect, useCallback, useRef, type RefObject } from 'react'
import { Play, Pause, SkipBack, SkipForward, ChevronLeft, ChevronRight } from 'lucide-react'

export interface PlaybackControlsProps {
  currentFrame: number
  totalFrames: number
  time: number
  onFrameChange: (frame: number) => void
  onPlayPauseChange?: (playing: boolean) => void
  isFrameReady?: (frame: number) => boolean
  defaultSpeed?: number
  className?: string
  /** If true, playback is controlled externally (no internal animation loop) */
  imperativeMode?: boolean
  /** Ref to frame index for imperative mode slider sync */
  frameIndexRef?: RefObject<number>
  /** Ref to playback speed for imperative mode */
  playbackSpeedRef?: RefObject<number>
}

/** Playback controls for animation */
export function PlaybackControls({
  currentFrame,
  totalFrames,
  time,
  onFrameChange,
  onPlayPauseChange,
  isFrameReady,
  defaultSpeed = 30,
  className = '',
  imperativeMode = false,
  frameIndexRef,
  playbackSpeedRef,
}: PlaybackControlsProps) {
  const [isPlaying, setIsPlaying] = useState(false)
  const [speed, setSpeed] = useState(defaultSpeed)
  const [displayFrame, setDisplayFrame] = useState(currentFrame)
  const intervalRef = useRef<NodeJS.Timeout | null>(null)
  
  // In imperative mode, poll the frameIndexRef for slider display
  useEffect(() => {
    if (imperativeMode && frameIndexRef) {
      const interval = setInterval(() => {
        if (frameIndexRef.current !== undefined) {
          setDisplayFrame(frameIndexRef.current)
        }
      }, 50) // 20 FPS UI update rate
      return () => clearInterval(interval)
    }
  }, [imperativeMode, frameIndexRef])
  
  // Sync displayFrame from prop in non-imperative mode
  useEffect(() => {
    if (!imperativeMode) {
      setDisplayFrame(currentFrame)
    }
  }, [currentFrame, imperativeMode])
  
  // Update playback speed ref when speed changes
  useEffect(() => {
    if (playbackSpeedRef && 'current' in playbackSpeedRef) {
      (playbackSpeedRef as { current: number }).current = speed
    }
  }, [speed, playbackSpeedRef])

  const play = useCallback(() => {
    setIsPlaying(true)
    onPlayPauseChange?.(true)
  }, [onPlayPauseChange])

  const pause = useCallback(() => {
    setIsPlaying(false)
    onPlayPauseChange?.(false)
  }, [onPlayPauseChange])

  const togglePlayPause = useCallback(() => {
    if (isPlaying) {
      pause()
    } else {
      play()
    }
  }, [isPlaying, play, pause])

  const goToStart = useCallback(() => {
    pause()
    onFrameChange(0)
  }, [pause, onFrameChange])

  const goToEnd = useCallback(() => {
    pause()
    onFrameChange(totalFrames - 1)
  }, [pause, onFrameChange, totalFrames])

  const stepForward = useCallback(() => {
    const frame = imperativeMode && frameIndexRef ? frameIndexRef.current : currentFrame
    if (frame < totalFrames - 1) {
      onFrameChange(frame + 1)
    }
  }, [currentFrame, totalFrames, onFrameChange, imperativeMode, frameIndexRef])

  const stepBackward = useCallback(() => {
    const frame = imperativeMode && frameIndexRef ? frameIndexRef.current : currentFrame
    if (frame > 0) {
      onFrameChange(frame - 1)
    }
  }, [currentFrame, onFrameChange, imperativeMode, frameIndexRef])

  // Animation loop - use ref to track current frame to avoid stale closure
  const currentFrameRef = useRef(currentFrame)
  useEffect(() => {
    currentFrameRef.current = currentFrame
  }, [currentFrame])

  // Store isFrameReady in a ref to avoid stale closures
  const isFrameReadyRef = useRef(isFrameReady)
  useEffect(() => {
    isFrameReadyRef.current = isFrameReady
  }, [isFrameReady])

  // Internal animation loop - only runs in non-imperative mode
  useEffect(() => {
    // Skip internal animation loop in imperative mode (parent handles it)
    if (imperativeMode) return
    
    if (isPlaying) {
      // Use requestAnimationFrame-based loop for smoother playback
      let lastTime = 0
      const frameInterval = 1000 / speed
      let animationId: number | null = null
      
      const tick = (timestamp: number) => {
        if (!lastTime) lastTime = timestamp
        const elapsed = timestamp - lastTime
        
        if (elapsed >= frameInterval) {
          const nextFrame = currentFrameRef.current + 1
          
          if (nextFrame >= totalFrames) {
            pause()
            onFrameChange(totalFrames - 1)
            return
          }
          
          // Only advance if next frame is ready (skip if not)
          const frameReady = isFrameReadyRef.current ? isFrameReadyRef.current(nextFrame) : true
          if (frameReady) {
            onFrameChange(nextFrame)
            lastTime = timestamp
          } else {
            // Frame not ready - wait a bit and retry
            // Don't update lastTime so we check again soon
          }
        }
        
        animationId = requestAnimationFrame(tick)
      }
      
      animationId = requestAnimationFrame(tick)
      
      return () => {
        if (animationId !== null) {
          cancelAnimationFrame(animationId)
        }
      }
    }
  }, [isPlaying, speed, totalFrames, pause, onFrameChange, imperativeMode])

  // Keyboard shortcuts
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.target instanceof HTMLInputElement) return

      switch (e.key) {
        case ' ':
          e.preventDefault()
          togglePlayPause()
          break
        case 'ArrowLeft':
          e.preventDefault()
          stepBackward()
          break
        case 'ArrowRight':
          e.preventDefault()
          stepForward()
          break
        case 'Home':
          e.preventDefault()
          goToStart()
          break
        case 'End':
          e.preventDefault()
          goToEnd()
          break
      }
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [togglePlayPause, stepBackward, stepForward, goToStart, goToEnd])

  const formatTime = (t: number) => {
    if (Math.abs(t) < 0.001) return t.toExponential(2)
    if (Math.abs(t) > 1000) return t.toExponential(2)
    return t.toFixed(4)
  }

  return (
    <div className={`bg-gray-800 p-3 rounded ${className}`}>
      <div className="flex items-center gap-4">
        {/* Transport controls */}
        <div className="flex items-center gap-1">
          <button
            onClick={goToStart}
            className="p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white"
            title="Go to start (Home)"
          >
            <SkipBack size={16} />
          </button>
          <button
            onClick={stepBackward}
            className="p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white"
            title="Step backward (←)"
            disabled={displayFrame === 0}
          >
            <ChevronLeft size={16} />
          </button>
          <button
            onClick={togglePlayPause}
            className="p-2 bg-blue-600 hover:bg-blue-500 rounded text-white"
            title={isPlaying ? 'Pause (Space)' : 'Play (Space)'}
          >
            {isPlaying ? <Pause size={18} /> : <Play size={18} />}
          </button>
          <button
            onClick={stepForward}
            className="p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white"
            title="Step forward (→)"
            disabled={displayFrame === totalFrames - 1}
          >
            <ChevronRight size={16} />
          </button>
          <button
            onClick={goToEnd}
            className="p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white"
            title="Go to end (End)"
          >
            <SkipForward size={16} />
          </button>
        </div>

        {/* Timeline slider */}
        <div className="flex-1">
          <input
            type="range"
            min={0}
            max={totalFrames - 1}
            value={displayFrame}
            onChange={(e) => onFrameChange(parseInt(e.target.value))}
            className="w-full h-2 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
          />
        </div>

        {/* Frame info */}
        <div className="text-xs text-gray-400 whitespace-nowrap min-w-[120px]">
          <span className="font-mono">{displayFrame + 1}</span>
          <span className="text-gray-600"> / </span>
          <span className="font-mono">{totalFrames}</span>
          <span className="text-gray-600 ml-2">t = </span>
          <span className="font-mono">{formatTime(time)}</span>
        </div>

        {/* Speed control */}
        <div className="flex items-center gap-2">
          <span className="text-xs text-gray-500">Speed:</span>
          <select
            value={speed}
            onChange={(e) => setSpeed(parseInt(e.target.value))}
            className="bg-gray-700 text-white text-xs rounded px-1 py-0.5"
          >
            <option value={5}>5 fps</option>
            <option value={10}>10 fps</option>
            <option value={15}>15 fps</option>
            <option value={20}>20 fps</option>
            <option value={30}>30 fps</option>
            <option value={45}>45 fps</option>
            <option value={60}>60 fps</option>
          </select>
        </div>
      </div>
    </div>
  )
}

export default PlaybackControls
