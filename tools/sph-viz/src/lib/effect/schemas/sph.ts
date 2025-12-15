/**
 * SPH Simulation Data Schemas using Effect Schema
 *
 * This module provides:
 * 1. TypeScript interfaces (mutable) for use in React components
 * 2. Effect Schemas for runtime validation
 * 3. Decode functions for parsing API responses
 * 
 * Note: The interfaces are defined separately from schemas to ensure
 * mutable types for React compatibility (Effect Schema creates readonly types).
 */

import { Schema } from 'effect'

// ============================================================================
// TYPE DEFINITIONS (Mutable for React compatibility)
// ============================================================================

/**
 * 3D Vector as mutable tuple [x, y, z]
 */
export type Vector3 = [number, number, number]

/**
 * Bounding box with min and max corners
 */
export interface BoundingBox {
  min: Vector3
  max: Vector3
}

/**
 * Physical unit conversion factors (CGS system)
 */
export interface PhysicalUnits {
  /** Mass unit in CGS (grams) */
  mass: number
  /** Length unit in CGS (cm) */
  length: number
  /** Time unit in CGS (seconds) */
  time: number
  /** Derived: length/time */
  velocity: number
  /** Derived: mass/length^3 */
  density: number
  /** Derived: mass*length^2/time^2 */
  energy: number
  /** Derived: mass/(length*time^2) */
  pressure: number
}

/**
 * IMBH (Intermediate Mass Black Hole) physics configuration
 * for encounter simulations (Oka et al. 2017 style)
 */
export interface IMBHPhysicsConfig {
  /** Whether IMBH physics is enabled */
  enabled: boolean
  /** BH position in code units (pc) */
  bhPosition: Vector3
  /** BH mass in code units (1000 M_sun) */
  bhMass: number
  /** Cloud initial position (pc) */
  cloudInitialPosition: Vector3
  /** Cloud initial velocity (km/s) */
  cloudInitialVelocity: Vector3
  /** Cloud mass in code units */
  cloudMass: number
  /** Cloud radius (pc) */
  cloudRadius: number
  /** Tidal radius (pc) */
  tidalRadius: number
  /** Impact parameter (pc) */
  impactParameter: number
  /** Pericentre distance (pc) */
  pericentre: number
  /** Orbital eccentricity (>1 for hyperbolic) */
  eccentricity: number
  /** Time unit conversion (Myr) */
  timeUnit: number
  /** Orbital plane inclination (degrees) - optional */
  inclination?: number
  /** Position angle on sky (degrees) - optional */
  positionAngle?: number
  /** V_LSR in km/s (negative = approaching) - optional */
  lsrVelocity?: number
}

/**
 * SPH method type
 */
export type SPHMethod = 'GSPH' | 'SSPH' | 'DISPH' | 'GDISPH' | 'SRGSPH' | string

/**
 * SPH kernel type
 */
export type SPHKernel = 'CubicSpline' | 'WendlandC2' | 'WendlandC4' | string

/**
 * Simulation dimensions
 */
export type Dimensions = 1 | 2 | 3

/**
 * Complete metadata for a simulation
 */
export interface SimulationMetadata {
  /** Unique identifier */
  id: string
  /** Human-readable name */
  name: string
  /** Description of the simulation */
  description: string
  /** SPH method used */
  method: SPHMethod
  /** Kernel function used */
  kernel: SPHKernel
  /** Number of spatial dimensions */
  dimensions: Dimensions
  /** Total number of frames */
  totalFrames: number
  /** Number of particles */
  particleCount: number
  /** Time range [start, end] in simulation time */
  timeRange: [number, number]
  /** Spatial bounding box */
  boundingBox: BoundingBox
  /** Physical unit conversions - optional */
  units?: PhysicalUnits
  /** IMBH physics parameters - optional */
  imbhPhysics?: IMBHPhysicsConfig
  /** Pre-computed P-V diagram ranges - optional */
  pvDiagram?: {
    positionRange: [number, number]
    velocityRange: [number, number]
    velocityStats?: {
      vx: [number, number]
      vy: [number, number]
      vz: [number, number]
    }
  }
  /** Path to config file - optional */
  configPath?: string
  /** Path to simulation data */
  dataPath: string
  /** ISO timestamp of creation */
  createdAt: string
}

/**
 * Field offsets within interleaved particle data
 */
export type FieldOffsets = Record<string, number>

/**
 * Raw frame data (before parsing)
 */
export interface FrameData {
  frameIndex: number
  time: number
  /** Base64 encoded Float32Array */
  data: string
  /** Number of floats per particle */
  stride: number
  /** Offset for each field in stride */
  fieldOffsets: FieldOffsets
  particleCount: number
}

/**
 * Parsed frame with typed arrays
 */
export interface ParsedFrame {
  frameIndex: number
  time: number
  particleCount: number
  /** xyz interleaved positions */
  positions: Float32Array
  /** vxvyvz interleaved velocities */
  velocities: Float32Array
  mass: Float32Array
  density: Float32Array
  pressure: Float32Array
  energy: Float32Array
  smoothingLength: Float32Array
  soundSpeed?: Float32Array
  machNumber?: Float32Array
  divV?: Float32Array
}

/**
 * Global statistics computed for a frame
 */
export interface FrameStatistics {
  frameIndex: number
  time: number
  totalMass: number
  totalKineticEnergy: number
  totalInternalEnergy: number
  totalEnergy: number
  momentum: Vector3
  centerOfMass: Vector3
  densityRange: [number, number]
  pressureRange: [number, number]
  temperatureRange?: [number, number]
  maxMach?: number
  /** Count of particles with divV < threshold */
  particlesInShock?: number
}

/**
 * Color map configuration
 */
export interface ColorMap {
  name: string
  /** Hex color strings for interpolation */
  colors: string[]
  min: number
  max: number
  logScale: boolean
}

/**
 * Base color map (without range)
 */
export interface BaseColorMap {
  name: string
  colors: string[]
  logScale: boolean
}

/**
 * Viewer state
 */
export interface ViewerState {
  currentFrame: number
  isPlaying: boolean
  /** Frames per second */
  playbackSpeed: number
  colorField: string
  colorMap: ColorMap
  pointSize: number
  showAxes: boolean
  showBoundingBox: boolean
  cameraPosition: Vector3
  cameraTarget: Vector3
}

/**
 * Response from GET /api/simulations
 */
export interface SimulationsListResponse {
  simulations: SimulationMetadata[]
}

/**
 * Response from GET /api/simulations/:simId/frames/:frameId
 */
export interface FrameResponse {
  frameIndex: number
  time: number
  /** Base64 encoded Float32Array */
  data: string
  stride: number
  fieldOffsets: FieldOffsets
  particleCount: number
}

/**
 * Response containing statistics for multiple frames
 */
export interface StatisticsResponse {
  frames: FrameStatistics[]
}

// ============================================================================
// EFFECT SCHEMAS (for runtime validation)
// ============================================================================

/**
 * 3D Vector schema
 */
export const Vector3Schema = Schema.Tuple(Schema.Number, Schema.Number, Schema.Number)

/**
 * Bounding box schema
 */
export const BoundingBoxSchema = Schema.Struct({
  min: Vector3Schema,
  max: Vector3Schema,
})

/**
 * Physical units schema
 */
export const PhysicalUnitsSchema = Schema.Struct({
  mass: Schema.Number,
  length: Schema.Number,
  time: Schema.Number,
  velocity: Schema.Number,
  density: Schema.Number,
  energy: Schema.Number,
  pressure: Schema.Number,
})

/**
 * IMBH physics config schema
 */
export const IMBHPhysicsConfigSchema = Schema.Struct({
  enabled: Schema.Boolean,
  bhPosition: Vector3Schema,
  bhMass: Schema.Number,
  cloudInitialPosition: Vector3Schema,
  cloudInitialVelocity: Vector3Schema,
  cloudMass: Schema.Number,
  cloudRadius: Schema.Number,
  tidalRadius: Schema.Number,
  impactParameter: Schema.Number,
  pericentre: Schema.Number,
  eccentricity: Schema.Number,
  timeUnit: Schema.Number,
  inclination: Schema.optional(Schema.Number),
  positionAngle: Schema.optional(Schema.Number),
  lsrVelocity: Schema.optional(Schema.Number),
})

/**
 * SPH method schema
 */
export const SPHMethodSchema = Schema.Union(
  Schema.Literal('GSPH'),
  Schema.Literal('SSPH'),
  Schema.Literal('DISPH'),
  Schema.Literal('GDISPH'),
  Schema.Literal('SRGSPH'),
  Schema.String
)

/**
 * SPH kernel schema
 */
export const SPHKernelSchema = Schema.Union(
  Schema.Literal('CubicSpline'),
  Schema.Literal('WendlandC2'),
  Schema.Literal('WendlandC4'),
  Schema.String
)

/**
 * Dimensions schema
 */
export const DimensionsSchema = Schema.Union(
  Schema.Literal(1),
  Schema.Literal(2),
  Schema.Literal(3)
)

/**
 * Field offsets schema
 */
export const FieldOffsetsSchema = Schema.Record({
  key: Schema.String,
  value: Schema.Number,
})

/**
 * Simulation metadata schema
 */
export const SimulationMetadataSchema = Schema.Struct({
  id: Schema.String,
  name: Schema.String,
  description: Schema.String,
  method: SPHMethodSchema,
  kernel: SPHKernelSchema,
  dimensions: DimensionsSchema,
  totalFrames: Schema.Number,
  particleCount: Schema.Number,
  timeRange: Schema.Tuple(Schema.Number, Schema.Number),
  boundingBox: BoundingBoxSchema,
  units: Schema.optional(PhysicalUnitsSchema),
  imbhPhysics: Schema.optional(IMBHPhysicsConfigSchema),
  configPath: Schema.optional(Schema.String),
  dataPath: Schema.String,
  createdAt: Schema.String,
})

/**
 * Frame data schema
 */
export const FrameDataSchema = Schema.Struct({
  frameIndex: Schema.Number,
  time: Schema.Number,
  data: Schema.String,
  stride: Schema.Number,
  fieldOffsets: FieldOffsetsSchema,
  particleCount: Schema.Number,
})

/**
 * Frame statistics schema
 */
export const FrameStatisticsSchema = Schema.Struct({
  frameIndex: Schema.Number,
  time: Schema.Number,
  totalMass: Schema.Number,
  totalKineticEnergy: Schema.Number,
  totalInternalEnergy: Schema.Number,
  totalEnergy: Schema.Number,
  momentum: Vector3Schema,
  centerOfMass: Vector3Schema,
  densityRange: Schema.Tuple(Schema.Number, Schema.Number),
  pressureRange: Schema.Tuple(Schema.Number, Schema.Number),
  temperatureRange: Schema.optional(Schema.Tuple(Schema.Number, Schema.Number)),
  maxMach: Schema.optional(Schema.Number),
  particlesInShock: Schema.optional(Schema.Number),
})

/**
 * Simulations list response schema
 */
export const SimulationsListResponseSchema = Schema.Struct({
  simulations: Schema.Array(SimulationMetadataSchema),
})

/**
 * Frame response schema
 */
export const FrameResponseSchema = Schema.Struct({
  frameIndex: Schema.Number,
  time: Schema.Number,
  data: Schema.String,
  stride: Schema.Number,
  fieldOffsets: FieldOffsetsSchema,
  particleCount: Schema.Number,
})

/**
 * Statistics response schema
 */
export const StatisticsResponseSchema = Schema.Struct({
  frames: Schema.Array(FrameStatisticsSchema),
})

/**
 * Color map schema
 */
export const ColorMapSchema = Schema.Struct({
  name: Schema.String,
  colors: Schema.Array(Schema.String),
  min: Schema.Number,
  max: Schema.Number,
  logScale: Schema.Boolean,
})

/**
 * Base color map schema
 */
export const BaseColorMapSchema = Schema.Struct({
  name: Schema.String,
  colors: Schema.Array(Schema.String),
  logScale: Schema.Boolean,
})

/**
 * Viewer state schema
 */
export const ViewerStateSchema = Schema.Struct({
  currentFrame: Schema.Number,
  isPlaying: Schema.Boolean,
  playbackSpeed: Schema.Number,
  colorField: Schema.String,
  colorMap: ColorMapSchema,
  pointSize: Schema.Number,
  showAxes: Schema.Boolean,
  showBoundingBox: Schema.Boolean,
  cameraPosition: Vector3Schema,
  cameraTarget: Vector3Schema,
})

// ============================================================================
// DEFAULT VALUES
// ============================================================================

/**
 * Default field offsets for standard SPH output
 */
export const DEFAULT_FIELD_OFFSETS: FieldOffsets = {
  x: 0,
  y: 1,
  z: 2,
  vx: 3,
  vy: 4,
  vz: 5,
  mass: 6,
  density: 7,
  pressure: 8,
  energy: 9,
  smoothingLength: 10,
}

// ============================================================================
// DECODE HELPERS
// ============================================================================

/**
 * Decode and validate simulation metadata
 * Returns the parsed data as our mutable interface type
 */
export const decodeSimulationMetadata = (data: unknown): SimulationMetadata | null => {
  try {
    const decoded = Schema.decodeUnknownSync(SimulationMetadataSchema)(data)
    // Convert readonly to mutable (safe because we just decoded it)
    return decoded as unknown as SimulationMetadata
  } catch {
    return null
  }
}

/**
 * Decode and validate simulations list response
 */
export const decodeSimulationsListResponse = (data: unknown): SimulationsListResponse | null => {
  try {
    const decoded = Schema.decodeUnknownSync(SimulationsListResponseSchema)(data)
    return decoded as unknown as SimulationsListResponse
  } catch {
    return null
  }
}

/**
 * Decode and validate frame response
 */
export const decodeFrameResponse = (data: unknown): FrameResponse | null => {
  try {
    const decoded = Schema.decodeUnknownSync(FrameResponseSchema)(data)
    return decoded as unknown as FrameResponse
  } catch {
    return null
  }
}

/**
 * Decode and validate frame statistics
 */
export const decodeFrameStatistics = (data: unknown): FrameStatistics | null => {
  try {
    const decoded = Schema.decodeUnknownSync(FrameStatisticsSchema)(data)
    return decoded as unknown as FrameStatistics
  } catch {
    return null
  }
}

/**
 * Decode simulations list with Effect Either
 */
export const decodeSimulationsListResponseEither = Schema.decodeUnknownEither(SimulationsListResponseSchema)

/**
 * Decode simulation metadata with Effect Either
 */
export const decodeSimulationMetadataEither = Schema.decodeUnknownEither(SimulationMetadataSchema)
