/**
 * SPH Simulation Data Types
 *
 * These types define the structure for SPH simulation data used across
 * the visualization tool. Compatible with various SPH methods (GSPH, SSPH,
 * DISPH, GDISPH, etc.) and test cases.
 */

/** Physical constants for unit conversion */
export interface PhysicalUnits {
  mass: number // Mass unit in CGS
  length: number // Length unit in CGS
  time: number // Time unit in CGS
  velocity: number // Derived: length/time
  density: number // Derived: mass/length^3
  energy: number // Derived: mass*length^2/time^2
  pressure: number // Derived: mass/(length*time^2)
}

/** IMBH physics configuration for encounter simulations */
export interface IMBHPhysicsConfig {
  enabled: boolean
  bhPosition: [number, number, number]      // BH position in code units (pc)
  bhMass: number                            // BH mass in code units (1000 M_sun)
  cloudInitialPosition: [number, number, number]  // Cloud initial position (pc)
  cloudInitialVelocity: [number, number, number]  // Cloud initial velocity (km/s)
  cloudMass: number                         // Cloud mass in code units
  cloudRadius: number                       // Cloud radius (pc)
  tidalRadius: number                       // Tidal radius (pc)
  impactParameter: number                   // Impact parameter (pc)
  timeUnit: number                          // Time unit conversion (Myr)
}

/** Metadata for a simulation */
export interface SimulationMetadata {
  id: string
  name: string
  description: string
  method: 'GSPH' | 'SSPH' | 'DISPH' | 'GDISPH' | 'SRGSPH' | string
  kernel: 'CubicSpline' | 'WendlandC2' | 'WendlandC4' | string
  dimensions: 1 | 2 | 3
  totalFrames: number
  particleCount: number
  timeRange: [number, number] // [start, end] in simulation time
  boundingBox: {
    min: [number, number, number]
    max: [number, number, number]
  }
  units?: PhysicalUnits
  imbhPhysics?: IMBHPhysicsConfig  // IMBH encounter physics parameters
  configPath?: string
  dataPath: string
  createdAt: string
}

/** Single particle data */
export interface Particle {
  id: number
  x: number
  y: number
  z: number
  vx: number
  vy: number
  vz: number
  mass: number
  density: number
  pressure: number
  energy: number // internal energy
  smoothingLength: number
  soundSpeed?: number
  machNumber?: number
  divV?: number // velocity divergence (shock indicator)
  curlV?: number // velocity curl magnitude
}

/** Frame data - a snapshot of all particles at a time */
export interface FrameData {
  frameIndex: number
  time: number
  particles: Float32Array // Interleaved: [x,y,z,vx,vy,vz,m,rho,p,u,h, ...]
  particleCount: number
  stride: number // Number of floats per particle
  fieldOffsets: Record<string, number> // Offset for each field in stride
}

/** Parsed frame with easy access to fields */
export interface ParsedFrame {
  frameIndex: number
  time: number
  particleCount: number
  positions: Float32Array // xyz interleaved
  velocities: Float32Array // vxvyvz interleaved
  mass: Float32Array
  density: Float32Array
  pressure: Float32Array
  energy: Float32Array
  smoothingLength: Float32Array
  soundSpeed?: Float32Array
  machNumber?: Float32Array
  divV?: Float32Array
}

/** Global simulation statistics for a frame */
export interface FrameStatistics {
  frameIndex: number
  time: number
  totalMass: number
  totalKineticEnergy: number
  totalInternalEnergy: number
  totalEnergy: number
  momentum: [number, number, number]
  centerOfMass: [number, number, number]
  densityRange: [number, number]
  pressureRange: [number, number]
  temperatureRange?: [number, number]
  maxMach?: number
  particlesInShock?: number // Count of particles with divV < threshold
}

/** Profile data for 1D plots */
export interface ProfileData {
  radius: Float32Array
  density: Float32Array
  pressure: Float32Array
  velocity: Float32Array
  temperature?: Float32Array
  analyticalDensity?: Float32Array
  analyticalPressure?: Float32Array
  analyticalVelocity?: Float32Array
}

/** Color mapping configuration */
export interface ColorMap {
  name: string
  colors: string[] // Hex colors for interpolation
  min: number
  max: number
  logScale: boolean
}

/** Visualization state */
export interface ViewerState {
  currentFrame: number
  isPlaying: boolean
  playbackSpeed: number // frames per second
  colorField: keyof ParsedFrame | string
  colorMap: ColorMap
  pointSize: number
  showAxes: boolean
  showBoundingBox: boolean
  cameraPosition: [number, number, number]
  cameraTarget: [number, number, number]
}

/** Available simulations response */
export interface SimulationsListResponse {
  simulations: SimulationMetadata[]
}

/** Frame data response */
export interface FrameResponse {
  frameIndex: number
  time: number
  data: string // Base64 encoded Float32Array
  stride: number
  fieldOffsets: Record<string, number>
  particleCount: number
}

/** Statistics response */
export interface StatisticsResponse {
  frames: FrameStatistics[]
}

/** Default field offsets for standard SPH output */
export const DEFAULT_FIELD_OFFSETS: Record<string, number> = {
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

/** Default color maps */
export const COLOR_MAPS: Record<string, Omit<ColorMap, 'min' | 'max'>> = {
  viridis: {
    name: 'Viridis',
    colors: ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'],
    logScale: false,
  },
  plasma: {
    name: 'Plasma',
    colors: ['#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786', '#d8576b', '#ed7953', '#fb9f3a', '#fdca26', '#f0f921'],
    logScale: false,
  },
  inferno: {
    name: 'Inferno',
    colors: ['#000004', '#1b0c41', '#4a0c6b', '#781c6d', '#a52c60', '#cf4446', '#ed6925', '#fb9b06', '#f7d13d', '#fcffa4'],
    logScale: false,
  },
  coolwarm: {
    name: 'Cool-Warm',
    colors: ['#3b4cc0', '#6688ee', '#88bbff', '#b8d4eb', '#dddddd', '#f5c4ad', '#f49a7b', '#d6604d', '#b40426'],
    logScale: false,
  },
  density: {
    name: 'Density',
    colors: ['#000033', '#000066', '#000099', '#0033cc', '#0066ff', '#00ccff', '#66ffcc', '#ccff66', '#ffcc00', '#ff6600', '#ff0000'],
    logScale: true,
  },
}
