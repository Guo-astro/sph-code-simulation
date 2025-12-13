/**
 * Simulation State Service using Effect
 * 
 * Provides a centralized state management layer for simulation data
 * using Effect's Context and Layer system for dependency injection.
 */

import { Context, Effect, Layer, Ref, Stream, pipe } from 'effect'
import type { SimulationMetadata, ParsedFrame, FrameStatistics, ColorMap } from '../schemas/sph'

// ============================================================================
// STATE TYPES
// ============================================================================

/**
 * Complete visualization state
 */
export interface VisualizationState {
  /** Currently selected simulation */
  selectedSimulation: SimulationMetadata | null
  /** Map of loaded frames (frameIndex -> ParsedFrame) */
  frames: Map<number, ParsedFrame>
  /** Statistics for each frame */
  statistics: FrameStatistics[]
  /** Current playback frame index */
  currentFrameIndex: number
  /** Whether playback is active */
  isPlaying: boolean
  /** Playback speed in FPS */
  playbackSpeed: number
  /** Loading state */
  isLoading: boolean
  /** Loading progress */
  loadingProgress: { loaded: number; total: number }
  /** Error message if any */
  error: string | null
  /** Color field for visualization */
  colorField: string
  /** Color map name */
  colorMapName: string
  /** Point size for particles */
  pointSize: number
  /** Opacity of particles */
  opacity: number
  /** Whether to show axes */
  showAxes: boolean
  /** Whether to show bounding box */
  showBoundingBox: boolean
  /** Color range (min, max) - [0,0] means auto */
  colorRange: [number, number]
  /** Whether to use log scale */
  useLogScale: boolean
  /** Whether to show 2D projections */
  showProjections: boolean
  /** Whether to show charts */
  showCharts: boolean
  /** Whether to use imperative mode for performance */
  useImperativeMode: boolean
  /** Current FPS (in imperative mode) */
  currentFps: number
  /** IMBH-specific settings */
  imbhSettings: {
    showBlackHole: boolean
    showTrajectory: boolean
    showRadii: boolean
    showGalacticMarkers: boolean
    showGalaxyDisk: boolean
    showSolarSystem: boolean
    animateGalaxy: boolean
    galaxyAnimationSpeed: number
    cameraMode: 'free' | 'earth'
  }
}

/**
 * Default visualization state
 */
export const defaultVisualizationState: VisualizationState = {
  selectedSimulation: null,
  frames: new Map(),
  statistics: [],
  currentFrameIndex: 0,
  isPlaying: false,
  playbackSpeed: 30,
  isLoading: false,
  loadingProgress: { loaded: 0, total: 0 },
  error: null,
  colorField: 'density',
  colorMapName: 'cosmicDawn',
  pointSize: 0.02,
  opacity: 0.8,
  showAxes: true,
  showBoundingBox: true,
  colorRange: [0, 0],
  useLogScale: false,
  showProjections: true,
  showCharts: true,
  useImperativeMode: true,
  currentFps: 0,
  imbhSettings: {
    showBlackHole: true,
    showTrajectory: true,
    showRadii: true,
    showGalacticMarkers: true,
    showGalaxyDisk: true,
    showSolarSystem: true,
    animateGalaxy: false,
    galaxyAnimationSpeed: 1.0,
    cameraMode: 'free',
  },
}

// ============================================================================
// STATE ACTIONS
// ============================================================================

export type StateAction =
  | { type: 'SELECT_SIMULATION'; payload: SimulationMetadata }
  | { type: 'CLEAR_SIMULATION' }
  | { type: 'SET_FRAME'; payload: { index: number; frame: ParsedFrame } }
  | { type: 'SET_FRAMES'; payload: Map<number, ParsedFrame> }
  | { type: 'SET_STATISTICS'; payload: FrameStatistics[] }
  | { type: 'SET_CURRENT_FRAME'; payload: number }
  | { type: 'SET_PLAYING'; payload: boolean }
  | { type: 'SET_PLAYBACK_SPEED'; payload: number }
  | { type: 'SET_LOADING'; payload: boolean }
  | { type: 'SET_LOADING_PROGRESS'; payload: { loaded: number; total: number } }
  | { type: 'SET_ERROR'; payload: string | null }
  | { type: 'SET_COLOR_FIELD'; payload: string }
  | { type: 'SET_COLOR_MAP_NAME'; payload: string }
  | { type: 'SET_POINT_SIZE'; payload: number }
  | { type: 'SET_OPACITY'; payload: number }
  | { type: 'SET_SHOW_AXES'; payload: boolean }
  | { type: 'SET_SHOW_BOUNDING_BOX'; payload: boolean }
  | { type: 'SET_COLOR_RANGE'; payload: [number, number] }
  | { type: 'SET_USE_LOG_SCALE'; payload: boolean }
  | { type: 'SET_SHOW_PROJECTIONS'; payload: boolean }
  | { type: 'SET_SHOW_CHARTS'; payload: boolean }
  | { type: 'SET_USE_IMPERATIVE_MODE'; payload: boolean }
  | { type: 'SET_CURRENT_FPS'; payload: number }
  | { type: 'SET_IMBH_SETTINGS'; payload: Partial<VisualizationState['imbhSettings']> }

/**
 * Pure state reducer
 */
export function stateReducer(state: VisualizationState, action: StateAction): VisualizationState {
  switch (action.type) {
    case 'SELECT_SIMULATION':
      return {
        ...state,
        selectedSimulation: action.payload,
        frames: new Map(),
        statistics: [],
        currentFrameIndex: 0,
        error: null,
      }
    case 'CLEAR_SIMULATION':
      return {
        ...state,
        selectedSimulation: null,
        frames: new Map(),
        statistics: [],
        currentFrameIndex: 0,
        isPlaying: false,
      }
    case 'SET_FRAME': {
      const newFrames = new Map(state.frames)
      newFrames.set(action.payload.index, action.payload.frame)
      return { ...state, frames: newFrames }
    }
    case 'SET_FRAMES':
      return { ...state, frames: action.payload }
    case 'SET_STATISTICS':
      return { ...state, statistics: action.payload }
    case 'SET_CURRENT_FRAME':
      return { ...state, currentFrameIndex: action.payload }
    case 'SET_PLAYING':
      return { ...state, isPlaying: action.payload }
    case 'SET_PLAYBACK_SPEED':
      return { ...state, playbackSpeed: action.payload }
    case 'SET_LOADING':
      return { ...state, isLoading: action.payload }
    case 'SET_LOADING_PROGRESS':
      return { ...state, loadingProgress: action.payload }
    case 'SET_ERROR':
      return { ...state, error: action.payload }
    case 'SET_COLOR_FIELD':
      return { ...state, colorField: action.payload }
    case 'SET_COLOR_MAP_NAME':
      return { ...state, colorMapName: action.payload }
    case 'SET_POINT_SIZE':
      return { ...state, pointSize: action.payload }
    case 'SET_OPACITY':
      return { ...state, opacity: action.payload }
    case 'SET_SHOW_AXES':
      return { ...state, showAxes: action.payload }
    case 'SET_SHOW_BOUNDING_BOX':
      return { ...state, showBoundingBox: action.payload }
    case 'SET_COLOR_RANGE':
      return { ...state, colorRange: action.payload }
    case 'SET_USE_LOG_SCALE':
      return { ...state, useLogScale: action.payload }
    case 'SET_SHOW_PROJECTIONS':
      return { ...state, showProjections: action.payload }
    case 'SET_SHOW_CHARTS':
      return { ...state, showCharts: action.payload }
    case 'SET_USE_IMPERATIVE_MODE':
      return { ...state, useImperativeMode: action.payload }
    case 'SET_CURRENT_FPS':
      return { ...state, currentFps: action.payload }
    case 'SET_IMBH_SETTINGS':
      return {
        ...state,
        imbhSettings: { ...state.imbhSettings, ...action.payload },
      }
    default:
      return state
  }
}

// ============================================================================
// STATE SERVICE
// ============================================================================

/**
 * Visualization State Service interface
 */
export interface VisualizationStateService {
  /** Get current state */
  readonly getState: Effect.Effect<VisualizationState>
  /** Dispatch an action to update state */
  readonly dispatch: (action: StateAction) => Effect.Effect<void>
  /** Subscribe to state changes */
  readonly subscribe: (callback: (state: VisualizationState) => void) => Effect.Effect<() => void>
  /** Get a specific slice of state */
  readonly select: <T>(selector: (state: VisualizationState) => T) => Effect.Effect<T>
}

/**
 * Service tag for dependency injection
 */
export class VisualizationStateServiceTag extends Context.Tag('VisualizationStateService')<
  VisualizationStateServiceTag,
  VisualizationStateService
>() {}

/**
 * Create the state service implementation
 */
export const makeVisualizationStateService = Effect.gen(function* () {
  const stateRef = yield* Ref.make(defaultVisualizationState)
  const subscribersRef = yield* Ref.make<Set<(state: VisualizationState) => void>>(new Set())

  const notifySubscribers = (state: VisualizationState) =>
    Effect.gen(function* () {
      const subscribers = yield* Ref.get(subscribersRef)
      for (const callback of subscribers) {
        callback(state)
      }
    })

  return {
    getState: Ref.get(stateRef),

    dispatch: (action: StateAction) =>
      Effect.gen(function* () {
        const currentState = yield* Ref.get(stateRef)
        const newState = stateReducer(currentState, action)
        yield* Ref.set(stateRef, newState)
        yield* notifySubscribers(newState)
      }),

    subscribe: (callback: (state: VisualizationState) => void) =>
      Effect.gen(function* () {
        yield* Ref.update(subscribersRef, (subscribers) => {
          const newSet = new Set(subscribers)
          newSet.add(callback)
          return newSet
        })

        // Return unsubscribe function
        return () => {
          Effect.runSync(
            Ref.update(subscribersRef, (subscribers) => {
              const newSet = new Set(subscribers)
              newSet.delete(callback)
              return newSet
            })
          )
        }
      }),

    select: <T>(selector: (state: VisualizationState) => T) =>
      Effect.gen(function* () {
        const state = yield* Ref.get(stateRef)
        return selector(state)
      }),
  }
})

/**
 * Layer providing the state service
 */
export const VisualizationStateServiceLive = Layer.effect(
  VisualizationStateServiceTag,
  makeVisualizationStateService
)

// ============================================================================
// HELPER EFFECTS
// ============================================================================

/**
 * Select simulation effect
 */
export const selectSimulation = (sim: SimulationMetadata) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SELECT_SIMULATION', payload: sim })
  })

/**
 * Set current frame effect
 */
export const setCurrentFrame = (frameIndex: number) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_CURRENT_FRAME', payload: frameIndex })
  })

/**
 * Toggle playback effect
 */
export const togglePlayback = Effect.gen(function* () {
  const service = yield* VisualizationStateServiceTag
  const state = yield* service.getState
  yield* service.dispatch({ type: 'SET_PLAYING', payload: !state.isPlaying })
})

/**
 * Set loading state effect
 */
export const setLoading = (isLoading: boolean) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_LOADING', payload: isLoading })
  })

/**
 * Set error effect
 */
export const setError = (error: string | null) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_ERROR', payload: error })
  })

/**
 * Update loading progress effect
 */
export const updateLoadingProgress = (loaded: number, total: number) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_LOADING_PROGRESS', payload: { loaded, total } })
  })

/**
 * Store frame effect
 */
export const storeFrame = (index: number, frame: ParsedFrame) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_FRAME', payload: { index, frame } })
  })

/**
 * Set all frames effect
 */
export const setAllFrames = (frames: Map<number, ParsedFrame>) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_FRAMES', payload: frames })
  })

/**
 * Set statistics effect
 */
export const setStatistics = (statistics: FrameStatistics[]) =>
  Effect.gen(function* () {
    const service = yield* VisualizationStateServiceTag
    yield* service.dispatch({ type: 'SET_STATISTICS', payload: statistics })
  })
