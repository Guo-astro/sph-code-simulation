/**
 * SPH Simulation API Client using Effect
 * 
 * Provides type-safe API access with automatic validation using Effect Schema.
 * All errors are properly typed and handled through the Effect system.
 */

import { Context, Effect, Layer, pipe, Schedule, Duration } from 'effect'
import { Schema } from 'effect'
import {
  SimulationMetadata,
  SimulationsListResponseSchema,
  decodeSimulationsListResponse,
  decodeSimulationMetadata,
  FieldOffsets,
  DEFAULT_FIELD_OFFSETS,
  type ParsedFrame,
  type FrameStatistics,
} from '../schemas/sph'

// ============================================================================
// ERROR TYPES
// ============================================================================

/**
 * API Error types for proper error handling
 */
export class ApiError extends Error {
  readonly _tag = 'ApiError'
  constructor(
    message: string,
    public readonly statusCode?: number,
    public readonly cause?: unknown
  ) {
    super(message)
    this.name = 'ApiError'
  }
}

export class NetworkError extends Error {
  readonly _tag = 'NetworkError'
  constructor(message: string, public readonly cause?: unknown) {
    super(message)
    this.name = 'NetworkError'
  }
}

export class ParseError extends Error {
  readonly _tag = 'ParseError'
  constructor(message: string, public readonly cause?: unknown) {
    super(message)
    this.name = 'ParseError'
  }
}

export class ValidationError extends Error {
  readonly _tag = 'ValidationError'
  constructor(message: string, public readonly cause?: unknown) {
    super(message)
    this.name = 'ValidationError'
  }
}

export type SPHApiError = ApiError | NetworkError | ParseError | ValidationError

// ============================================================================
// API CLIENT CONFIG
// ============================================================================

export interface ApiClientConfig {
  baseUrl: string
  retries: number
  retryDelay: Duration.Duration
}

export const defaultApiClientConfig: ApiClientConfig = {
  baseUrl: '',
  retries: 2,
  retryDelay: Duration.millis(100),
}

// ============================================================================
// API CLIENT SERVICE
// ============================================================================

/**
 * SPH API Client Service interface
 */
export interface SPHApiClient {
  /** Fetch all available simulations */
  readonly fetchSimulations: Effect.Effect<SimulationMetadata[], SPHApiError>
  
  /** Fetch metadata for a specific simulation */
  readonly fetchSimulationMetadata: (
    simId: string
  ) => Effect.Effect<SimulationMetadata, SPHApiError>
  
  /** Fetch a single frame as binary data */
  readonly fetchFrame: (
    simId: string,
    frameIndex: number
  ) => Effect.Effect<{ frame: ParsedFrame; stats: FrameStatistics }, SPHApiError>
  
  /** Fetch all frames for a simulation */
  readonly fetchAllFrames: (
    sim: SimulationMetadata,
    onProgress?: (loaded: number, total: number) => void
  ) => Effect.Effect<{ frames: Map<number, ParsedFrame>; statistics: FrameStatistics[] }, SPHApiError>
}

/**
 * Service tag for dependency injection
 */
export class SPHApiClientTag extends Context.Tag('SPHApiClient')<
  SPHApiClientTag,
  SPHApiClient
>() {}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Parse binary frame data from ArrayBuffer
 */
function parseBinaryFrame(
  buffer: ArrayBuffer,
  frameIndex: number,
  time: number,
  particleCount: number,
  stride: number,
  fieldOffsets: FieldOffsets
): ParsedFrame {
  const floats = new Float32Array(buffer)

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

/**
 * Compute statistics for a parsed frame
 */
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

    const v2 = vx * vx + vy * vy + vz * vz
    totalKinetic += 0.5 * m * v2
    totalInternal += m * u

    momentum[0] += m * vx
    momentum[1] += m * vy
    momentum[2] += m * vz

    com[0] += m * frame.positions[i * 3]
    com[1] += m * frame.positions[i * 3 + 1]
    com[2] += m * frame.positions[i * 3 + 2]

    totalMass += m

    if (isFinite(rho)) {
      if (rho < minDensity) minDensity = rho
      if (rho > maxDensity) maxDensity = rho
    }
    if (isFinite(P)) {
      if (P < minPressure) minPressure = P
      if (P > maxPressure) maxPressure = P
    }
  }

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

/**
 * Encode simulation ID for URL
 */
function encodeSimId(simId: string): string {
  return encodeURIComponent(simId.replace(/\//g, '|'))
}

// ============================================================================
// API CLIENT IMPLEMENTATION
// ============================================================================

/**
 * Create the API client implementation
 */
export const makeSPHApiClient = (config: ApiClientConfig = defaultApiClientConfig): SPHApiClient => {
  const { baseUrl, retries, retryDelay } = config

  /**
   * Fetch with retry logic
   */
  const fetchWithRetry = <T>(
    url: string,
    options?: RequestInit,
    parser?: (response: Response) => Promise<T>
  ): Effect.Effect<T, SPHApiError> =>
    pipe(
      Effect.tryPromise({
        try: async () => {
          const response = await fetch(`${baseUrl}${url}`, options)
          if (!response.ok) {
            throw new ApiError(`HTTP ${response.status}`, response.status)
          }
          if (parser) {
            return parser(response)
          }
          return response.json() as Promise<T>
        },
        catch: (error) => {
          if (error instanceof ApiError) return error
          return new NetworkError('Network request failed', error)
        },
      }),
      Effect.retry(
        Schedule.exponential(retryDelay).pipe(
          Schedule.compose(Schedule.recurs(retries))
        )
      )
    )

  /**
   * Fetch simulations list
   */
  const fetchSimulations: Effect.Effect<SimulationMetadata[], SPHApiError> =
    pipe(
      fetchWithRetry<unknown>('/api/simulations'),
      Effect.flatMap((data) => {
        const decoded = decodeSimulationsListResponse(data)
        if (decoded === null) {
          return Effect.fail(new ValidationError('Invalid simulations response'))
        }
        return Effect.succeed(decoded.simulations)
      })
    )

  /**
   * Fetch single simulation metadata
   */
  const fetchSimulationMetadata = (
    simId: string
  ): Effect.Effect<SimulationMetadata, SPHApiError> =>
    pipe(
      fetchWithRetry<unknown>(`/api/simulations/${encodeSimId(simId)}`),
      Effect.flatMap((data) => {
        const decoded = decodeSimulationMetadata(data)
        if (decoded === null) {
          return Effect.fail(new ValidationError('Invalid simulation metadata'))
        }
        return Effect.succeed(decoded)
      })
    )

  /**
   * Fetch a single frame with binary format
   */
  const fetchFrame = (
    simId: string,
    frameIndex: number
  ): Effect.Effect<{ frame: ParsedFrame; stats: FrameStatistics }, SPHApiError> =>
    pipe(
      Effect.tryPromise({
        try: async () => {
          const url = `${baseUrl}/api/simulations/${encodeSimId(simId)}/frames/${frameIndex}?format=binary`
          const response = await fetch(url)
          
          if (!response.ok) {
            throw new ApiError(`HTTP ${response.status}`, response.status)
          }

          const frameIdx = parseInt(response.headers.get('X-Frame-Index') || '0')
          const time = parseFloat(response.headers.get('X-Frame-Time') || '0')
          const particleCount = parseInt(response.headers.get('X-Particle-Count') || '0')
          const stride = parseInt(response.headers.get('X-Stride') || '11')
          const fieldOffsetsHeader = response.headers.get('X-Field-Offsets')
          const fieldOffsets: FieldOffsets = fieldOffsetsHeader 
            ? JSON.parse(fieldOffsetsHeader) 
            : DEFAULT_FIELD_OFFSETS

          const buffer = await response.arrayBuffer()

          // Validate buffer
          if (buffer.byteLength === 0) {
            throw new ParseError(`Frame ${frameIndex}: empty buffer received`)
          }
          if (buffer.byteLength % 4 !== 0) {
            throw new ParseError(`Frame ${frameIndex}: buffer size ${buffer.byteLength} not divisible by 4`)
          }

          const frame = parseBinaryFrame(buffer, frameIdx, time, particleCount, stride, fieldOffsets)
          const stats = computeFrameStatistics(frame)

          return { frame, stats }
        },
        catch: (error) => {
          if (error instanceof ApiError || error instanceof ParseError) return error
          return new NetworkError(`Failed to fetch frame ${frameIndex}`, error)
        },
      }),
      Effect.retry(
        Schedule.exponential(retryDelay).pipe(
          Schedule.compose(Schedule.recurs(retries))
        )
      )
    )

  /**
   * Fetch all frames for a simulation in parallel batches
   */
  const fetchAllFrames = (
    sim: SimulationMetadata,
    onProgress?: (loaded: number, total: number) => void
  ): Effect.Effect<{ frames: Map<number, ParsedFrame>; statistics: FrameStatistics[] }, SPHApiError> =>
    Effect.gen(function* () {
      const totalFrames = sim.totalFrames
      const simId = sim.id
      const batchSize = 10

      const newFrames = new Map<number, ParsedFrame>()
      const newStats: FrameStatistics[] = []

      let loadedCount = 0

      for (let batchStart = 0; batchStart < totalFrames; batchStart += batchSize) {
        const batchEnd = Math.min(batchStart + batchSize, totalFrames)
        const batchEffects: Effect.Effect<{ idx: number; result: { frame: ParsedFrame; stats: FrameStatistics } | null }, SPHApiError>[] = []

        for (let i = batchStart; i < batchEnd; i++) {
          const frameEffect = pipe(
            fetchFrame(simId, i),
            Effect.map((result) => ({ idx: i, result })),
            Effect.catchAll(() => Effect.succeed({ idx: i, result: null }))
          )
          batchEffects.push(frameEffect)
        }

        const results = yield* Effect.all(batchEffects, { concurrency: 'unbounded' })

        for (const { idx, result } of results) {
          if (result) {
            newFrames.set(idx, result.frame)
            newStats[idx] = result.stats
          }
          loadedCount++
        }

        if (onProgress) {
          onProgress(loadedCount, totalFrames)
        }
      }

      return { frames: newFrames, statistics: newStats }
    })

  return {
    fetchSimulations,
    fetchSimulationMetadata,
    fetchFrame,
    fetchAllFrames,
  }
}

/**
 * Layer providing the API client with default config
 */
export const SPHApiClientLive = Layer.succeed(
  SPHApiClientTag,
  makeSPHApiClient(defaultApiClientConfig)
)

/**
 * Create a custom API client layer
 */
export const makeSPHApiClientLayer = (config: ApiClientConfig) =>
  Layer.succeed(SPHApiClientTag, makeSPHApiClient(config))

// ============================================================================
// CONVENIENCE EFFECTS
// ============================================================================

/**
 * Fetch simulations using the API client from context
 */
export const fetchSimulations = Effect.gen(function* () {
  const client = yield* SPHApiClientTag
  return yield* client.fetchSimulations
})

/**
 * Fetch simulation metadata using the API client from context
 */
export const fetchSimulationMetadata = (simId: string) =>
  Effect.gen(function* () {
    const client = yield* SPHApiClientTag
    return yield* client.fetchSimulationMetadata(simId)
  })

/**
 * Fetch frame using the API client from context
 */
export const fetchFrame = (simId: string, frameIndex: number) =>
  Effect.gen(function* () {
    const client = yield* SPHApiClientTag
    return yield* client.fetchFrame(simId, frameIndex)
  })

/**
 * Fetch all frames using the API client from context
 */
export const fetchAllFrames = (
  sim: SimulationMetadata,
  onProgress?: (loaded: number, total: number) => void
) =>
  Effect.gen(function* () {
    const client = yield* SPHApiClientTag
    return yield* client.fetchAllFrames(sim, onProgress)
  })
