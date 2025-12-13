/**
 * React Hooks for Effect Integration
 * 
 * Provides React hooks to interact with Effect services,
 * bridging the functional Effect system with React's component model.
 */

'use client'

import { createContext, useContext, useEffect, useState, useCallback, useRef, type ReactNode } from 'react'
import { Effect, Exit } from 'effect'
import {
  type VisualizationState,
  type StateAction,
  defaultVisualizationState,
  stateReducer,
} from '../services/state'
import {
  type SPHApiClient,
  makeSPHApiClient,
} from '../services/api'
import type { SimulationMetadata, ParsedFrame, FrameStatistics } from '../schemas/sph'

// ============================================================================
// CONTEXT TYPES
// ============================================================================

interface EffectContextValue {
  // State
  state: VisualizationState
  
  // State Actions
  dispatch: (action: StateAction) => void
  
  // Simulation Actions
  loadSimulations: () => Promise<SimulationMetadata[]>
  selectSimulation: (sim: SimulationMetadata) => Promise<void>
  loadAllFrames: (sim: SimulationMetadata) => Promise<void>
  
  // Playback Actions
  setCurrentFrame: (frame: number) => void
  setPlaying: (playing: boolean) => void
  setPlaybackSpeed: (speed: number) => void
  
  // Visualization Settings
  setColorField: (field: string) => void
  setColorMapName: (name: string) => void
  setPointSize: (size: number) => void
  setOpacity: (opacity: number) => void
  setShowAxes: (show: boolean) => void
  setShowBoundingBox: (show: boolean) => void
  setColorRange: (range: [number, number]) => void
  setUseLogScale: (log: boolean) => void
  setShowProjections: (show: boolean) => void
  setShowCharts: (show: boolean) => void
  setUseImperativeMode: (use: boolean) => void
  
  // IMBH Settings
  setImbhSettings: (settings: Partial<VisualizationState['imbhSettings']>) => void
  
  // Loading State
  isLoading: boolean
  error: string | null
}

const EffectContext = createContext<EffectContextValue | null>(null)

// ============================================================================
// PROVIDER COMPONENT
// ============================================================================

interface EffectProviderProps {
  children: ReactNode
}

/**
 * Provider component that sets up Effect runtime and state management
 * Uses a simple reducer-based approach for state management
 */
export function EffectProvider({ children }: EffectProviderProps) {
  const [state, setState] = useState<VisualizationState>(defaultVisualizationState)
  const apiClientRef = useRef<SPHApiClient | null>(null)

  // Initialize API client on mount
  useEffect(() => {
    apiClientRef.current = makeSPHApiClient()
  }, [])

  // Dispatch action using the reducer
  const dispatch = useCallback((action: StateAction) => {
    setState(currentState => stateReducer(currentState, action))
  }, [])

  // Load simulations
  const loadSimulations = useCallback(async (): Promise<SimulationMetadata[]> => {
    if (!apiClientRef.current) {
      return []
    }

    dispatch({ type: 'SET_LOADING', payload: true })
    dispatch({ type: 'SET_ERROR', payload: null })

    try {
      const result = await Effect.runPromiseExit(apiClientRef.current.fetchSimulations)
      
      if (Exit.isSuccess(result)) {
        dispatch({ type: 'SET_LOADING', payload: false })
        return result.value
      } else {
        const error = result.cause
        dispatch({ type: 'SET_ERROR', payload: 'Failed to load simulations' })
        dispatch({ type: 'SET_LOADING', payload: false })
        console.error('Failed to load simulations:', error)
        return []
      }
    } catch (error) {
      dispatch({ type: 'SET_ERROR', payload: 'Failed to load simulations' })
      dispatch({ type: 'SET_LOADING', payload: false })
      console.error('Failed to load simulations:', error)
      return []
    }
  }, [dispatch])

  // Select simulation
  const selectSimulation = useCallback(async (sim: SimulationMetadata) => {
    dispatch({ type: 'SELECT_SIMULATION', payload: sim })
  }, [dispatch])

  // Load all frames for a simulation
  const loadAllFrames = useCallback(async (sim: SimulationMetadata) => {
    if (!apiClientRef.current) return

    dispatch({ type: 'SET_LOADING', payload: true })
    dispatch({ type: 'SET_LOADING_PROGRESS', payload: { loaded: 0, total: sim.totalFrames } })
    dispatch({ type: 'SET_ERROR', payload: null })

    try {
      const result = await Effect.runPromiseExit(
        apiClientRef.current.fetchAllFrames(sim, (loaded, total) => {
          dispatch({ type: 'SET_LOADING_PROGRESS', payload: { loaded, total } })
        })
      )

      if (Exit.isSuccess(result)) {
        dispatch({ type: 'SET_FRAMES', payload: result.value.frames })
        dispatch({ type: 'SET_STATISTICS', payload: result.value.statistics })
        dispatch({ type: 'SET_LOADING', payload: false })
      } else {
        dispatch({ type: 'SET_ERROR', payload: 'Failed to load frames' })
        dispatch({ type: 'SET_LOADING', payload: false })
        console.error('Failed to load frames:', result.cause)
      }
    } catch (error) {
      dispatch({ type: 'SET_ERROR', payload: 'Failed to load frames' })
      dispatch({ type: 'SET_LOADING', payload: false })
      console.error('Failed to load frames:', error)
    }
  }, [dispatch])

  // Simple state setters
  const setCurrentFrame = useCallback((frame: number) => {
    dispatch({ type: 'SET_CURRENT_FRAME', payload: frame })
  }, [dispatch])

  const setPlaying = useCallback((playing: boolean) => {
    dispatch({ type: 'SET_PLAYING', payload: playing })
  }, [dispatch])

  const setPlaybackSpeed = useCallback((speed: number) => {
    dispatch({ type: 'SET_PLAYBACK_SPEED', payload: speed })
  }, [dispatch])

  const setColorField = useCallback((field: string) => {
    dispatch({ type: 'SET_COLOR_FIELD', payload: field })
  }, [dispatch])

  const setColorMapName = useCallback((name: string) => {
    dispatch({ type: 'SET_COLOR_MAP_NAME', payload: name })
  }, [dispatch])

  const setPointSize = useCallback((size: number) => {
    dispatch({ type: 'SET_POINT_SIZE', payload: size })
  }, [dispatch])

  const setOpacity = useCallback((opacity: number) => {
    dispatch({ type: 'SET_OPACITY', payload: opacity })
  }, [dispatch])

  const setShowAxes = useCallback((show: boolean) => {
    dispatch({ type: 'SET_SHOW_AXES', payload: show })
  }, [dispatch])

  const setShowBoundingBox = useCallback((show: boolean) => {
    dispatch({ type: 'SET_SHOW_BOUNDING_BOX', payload: show })
  }, [dispatch])

  const setColorRange = useCallback((range: [number, number]) => {
    dispatch({ type: 'SET_COLOR_RANGE', payload: range })
  }, [dispatch])

  const setUseLogScale = useCallback((log: boolean) => {
    dispatch({ type: 'SET_USE_LOG_SCALE', payload: log })
  }, [dispatch])

  const setShowProjections = useCallback((show: boolean) => {
    dispatch({ type: 'SET_SHOW_PROJECTIONS', payload: show })
  }, [dispatch])

  const setShowCharts = useCallback((show: boolean) => {
    dispatch({ type: 'SET_SHOW_CHARTS', payload: show })
  }, [dispatch])

  const setUseImperativeMode = useCallback((use: boolean) => {
    dispatch({ type: 'SET_USE_IMPERATIVE_MODE', payload: use })
  }, [dispatch])

  const setImbhSettings = useCallback((settings: Partial<VisualizationState['imbhSettings']>) => {
    dispatch({ type: 'SET_IMBH_SETTINGS', payload: settings })
  }, [dispatch])

  const contextValue: EffectContextValue = {
    state,
    dispatch,
    loadSimulations,
    selectSimulation,
    loadAllFrames,
    setCurrentFrame,
    setPlaying,
    setPlaybackSpeed,
    setColorField,
    setColorMapName,
    setPointSize,
    setOpacity,
    setShowAxes,
    setShowBoundingBox,
    setColorRange,
    setUseLogScale,
    setShowProjections,
    setShowCharts,
    setUseImperativeMode,
    setImbhSettings,
    isLoading: state.isLoading,
    error: state.error,
  }

  return (
    <EffectContext.Provider value={contextValue}>
      {children}
    </EffectContext.Provider>
  )
}

// ============================================================================
// HOOKS
// ============================================================================

/**
 * Hook to access the Effect context
 */
export function useEffect_() {
  const context = useContext(EffectContext)
  if (!context) {
    throw new Error('useEffect_ must be used within an EffectProvider')
  }
  return context
}

/**
 * Hook to access visualization state
 */
export function useVisualizationState() {
  const { state } = useEffect_()
  return state
}

/**
 * Hook to access selected simulation
 */
export function useSelectedSimulation() {
  const { state } = useEffect_()
  return state.selectedSimulation
}

/**
 * Hook to access frames
 */
export function useFrames() {
  const { state } = useEffect_()
  return state.frames
}

/**
 * Hook to access current frame
 */
export function useCurrentFrame() {
  const { state } = useEffect_()
  return {
    frameIndex: state.currentFrameIndex,
    frame: state.frames.get(state.currentFrameIndex) ?? null,
  }
}

/**
 * Hook to access statistics
 */
export function useStatistics() {
  const { state } = useEffect_()
  return state.statistics
}

/**
 * Hook to access playback state
 */
export function usePlayback() {
  const { state, setCurrentFrame, setPlaying, setPlaybackSpeed } = useEffect_()
  return {
    currentFrame: state.currentFrameIndex,
    isPlaying: state.isPlaying,
    playbackSpeed: state.playbackSpeed,
    setCurrentFrame,
    setPlaying,
    setPlaybackSpeed,
  }
}

/**
 * Hook to access visualization settings
 */
export function useVisualizationSettings() {
  const {
    state,
    setColorField,
    setColorMapName,
    setPointSize,
    setOpacity,
    setShowAxes,
    setShowBoundingBox,
    setColorRange,
    setUseLogScale,
    setShowProjections,
    setShowCharts,
    setUseImperativeMode,
  } = useEffect_()

  return {
    colorField: state.colorField,
    colorMapName: state.colorMapName,
    pointSize: state.pointSize,
    opacity: state.opacity,
    showAxes: state.showAxes,
    showBoundingBox: state.showBoundingBox,
    colorRange: state.colorRange,
    useLogScale: state.useLogScale,
    showProjections: state.showProjections,
    showCharts: state.showCharts,
    useImperativeMode: state.useImperativeMode,
    setColorField,
    setColorMapName,
    setPointSize,
    setOpacity,
    setShowAxes,
    setShowBoundingBox,
    setColorRange,
    setUseLogScale,
    setShowProjections,
    setShowCharts,
    setUseImperativeMode,
  }
}

/**
 * Hook to access IMBH settings
 */
export function useImbhSettings() {
  const { state, setImbhSettings } = useEffect_()
  return {
    settings: state.imbhSettings,
    setSettings: setImbhSettings,
  }
}

/**
 * Hook to access loading state
 */
export function useLoadingState() {
  const { state } = useEffect_()
  return {
    isLoading: state.isLoading,
    loadingProgress: state.loadingProgress,
    error: state.error,
  }
}

/**
 * Hook for simulation management
 */
export function useSimulations() {
  const { loadSimulations, selectSimulation, loadAllFrames, state } = useEffect_()
  
  return {
    simulations: [] as SimulationMetadata[], // Would need separate state
    selectedSimulation: state.selectedSimulation,
    isLoading: state.isLoading,
    error: state.error,
    loadSimulations,
    selectSimulation,
    loadAllFrames,
  }
}
