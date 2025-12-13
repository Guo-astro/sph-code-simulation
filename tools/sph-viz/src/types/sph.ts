/**
 * SPH Simulation Data Types
 *
 * These types define the structure for SPH simulation data used across
 * the visualization tool. Compatible with various SPH methods (GSPH, SSPH,
 * DISPH, GDISPH, etc.) and test cases.
 * 
 * NOTE: Effect Schema versions are available in '~/lib/effect/schemas/sph'
 * for runtime validation. These types are kept for backward compatibility.
 */

// Re-export types from Effect Schema for convenience
export type {
  Vector3,
  BoundingBox,
  PhysicalUnits,
  IMBHPhysicsConfig,
  SPHMethod,
  SPHKernel,
  Dimensions,
  SimulationMetadata,
  FieldOffsets,
  FrameData,
  ParsedFrame,
  FrameStatistics,
  SimulationsListResponse,
  FrameResponse,
  StatisticsResponse,
  ColorMap,
  BaseColorMap,
  ViewerState,
} from '~/lib/effect/schemas/sph'

// Re-export schemas for validation
export {
  Vector3Schema,
  BoundingBoxSchema,
  PhysicalUnitsSchema,
  IMBHPhysicsConfigSchema,
  SPHMethodSchema,
  SPHKernelSchema,
  DimensionsSchema,
  SimulationMetadataSchema,
  FieldOffsetsSchema,
  FrameDataSchema,
  FrameStatisticsSchema,
  SimulationsListResponseSchema,
  FrameResponseSchema,
  StatisticsResponseSchema,
  ColorMapSchema,
  BaseColorMapSchema,
  ViewerStateSchema,
  DEFAULT_FIELD_OFFSETS,
  decodeSimulationMetadata,
  decodeSimulationsListResponse,
  decodeFrameResponse,
  decodeFrameStatistics,
} from '~/lib/effect/schemas/sph'

// ============================================================================
// LEGACY INTERFACES (for backward compatibility)
// These are kept for existing code that uses the old interface definitions
// ============================================================================

/** 
 * @deprecated Use the type from '~/lib/effect/schemas/sph' instead
 * Physical constants for unit conversion 
 */
export interface PhysicalUnitsLegacy {
  mass: number // Mass unit in CGS
  length: number // Length unit in CGS
  time: number // Time unit in CGS
  velocity: number // Derived: length/time
  density: number // Derived: mass/length^3
  energy: number // Derived: mass*length^2/time^2
  pressure: number // Derived: mass/(length*time^2)
}

// ============================================================================
// ADDITIONAL LEGACY TYPES
// ============================================================================

/** Single particle data (used for detailed inspection) */
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

// ============================================================================
// COLOR MAPS (predefined)
// ============================================================================

/**
 * Professionally designed color maps optimized for:
 * - Dark backgrounds (bg-gray-900, #0a0a0f canvas)
 * - High visibility and contrast
 * - Scientific accuracy (perceptually uniform where applicable)
 * - Colorblind accessibility
 * 
 * Each colormap avoids:
 * - Pure red on black (hard to see)
 * - Light green/yellow on light backgrounds
 * - Low-saturation colors that disappear on dark backgrounds
 */
export const COLOR_MAPS: Record<string, Omit<import('~/lib/effect/schemas/sph').ColorMap, 'min' | 'max'>> = {
  // ═══════════════════════════════════════════════════════════════════════════
  // SEQUENTIAL COLORMAPS (for continuous data like density, velocity magnitude)
  // ═══════════════════════════════════════════════════════════════════════════
  
  /** 
   * Cosmic Dawn - Deep space theme, excellent for density
   * Gradient: Deep purple → Electric blue → Cyan → Gold → White
   * High contrast on dark backgrounds, avoids muddy mid-tones
   */
  cosmicDawn: {
    name: 'Cosmic Dawn',
    colors: [
      '#1a0533',  // Deep purple-black (visible, not pure black)
      '#2d1b69',  // Rich purple
      '#3d4fc7',  // Royal blue
      '#00b4d8',  // Bright cyan (high visibility)
      '#48cae4',  // Light cyan
      '#90e0ef',  // Pale cyan
      '#ffd166',  // Warm gold
      '#ffeb99',  // Light gold
      '#ffffff',  // Pure white (maximum values)
    ],
    logScale: false,
  },

  /**
   * Nebula Fire - Warm colormap for energy/temperature
   * Gradient: Deep magenta → Hot pink → Orange → Yellow → White
   * Avoids pure red which is hard to see on black
   */
  nebulaFire: {
    name: 'Nebula Fire',
    colors: [
      '#2d0a31',  // Deep magenta-black
      '#5c1158',  // Dark magenta
      '#9b2d7f',  // Bright magenta (visible on dark)
      '#d64292',  // Hot pink
      '#f46e6e',  // Coral (not pure red)
      '#ff9f43',  // Bright orange
      '#ffc93c',  // Golden yellow
      '#fff176',  // Light yellow
      '#ffffff',  // White hot
    ],
    logScale: false,
  },

  /**
   * Ocean Depths - Cool colormap for pressure/potential
   * Gradient: Deep teal → Electric blue → Aqua → Mint → Cream
   * Excellent contrast, no muddy greens
   */
  oceanDepths: {
    name: 'Ocean Depths',
    colors: [
      '#0d1b2a',  // Deep navy (visible, not black)
      '#1b3a4b',  // Dark teal
      '#144d7e',  // Ocean blue
      '#1e88e5',  // Bright blue (high visibility)
      '#42a5f5',  // Sky blue
      '#80deea',  // Cyan
      '#a7ffeb',  // Aqua mint
      '#e0f7fa',  // Pale cyan
      '#fffde7',  // Warm cream
    ],
    logScale: false,
  },

  /**
   * Aurora - Vivid multi-hue for general purpose
   * Gradient: Purple → Blue → Teal → Green → Yellow
   * Perceptually uniform, colorblind-safe
   */
  aurora: {
    name: 'Aurora',
    colors: [
      '#3a0ca3',  // Deep violet
      '#4361ee',  // Electric blue
      '#4cc9f0',  // Bright cyan
      '#06d6a0',  // Teal-green (not pure green)
      '#52b788',  // Forest green
      '#99d98c',  // Light green
      '#d9ed92',  // Yellow-green
      '#fcf6bd',  // Pale yellow
      '#fff3b0',  // Cream
    ],
    logScale: false,
  },

  // ═══════════════════════════════════════════════════════════════════════════
  // SCIENTIFIC COLORMAPS (perceptually uniform)
  // ═══════════════════════════════════════════════════════════════════════════

  /**
   * Viridis - Scientific standard, perceptually uniform
   * Optimized version with brighter endpoints for dark backgrounds
   */
  viridis: {
    name: 'Viridis',
    colors: [
      '#440154',  // Deep purple
      '#482878',  // Purple
      '#3e4989',  // Blue-purple
      '#31688e',  // Steel blue
      '#26828e',  // Teal
      '#1f9e89',  // Teal-green
      '#35b779',  // Green
      '#6ece58',  // Yellow-green
      '#b5de2b',  // Lime
      '#fde725',  // Bright yellow
    ],
    logScale: false,
  },

  /**
   * Plasma - High-energy scientific colormap
   * Better for dark backgrounds than inferno
   */
  plasma: {
    name: 'Plasma',
    colors: [
      '#0d0887',  // Deep blue
      '#46039f',  // Purple
      '#7201a8',  // Magenta
      '#9c179e',  // Pink-purple
      '#bd3786',  // Hot pink
      '#d8576b',  // Coral
      '#ed7953',  // Orange
      '#fb9f3a',  // Gold
      '#fdca26',  // Yellow
      '#f0f921',  // Bright yellow
    ],
    logScale: false,
  },

  /**
   * Turbo - Rainbow without the problems
   * High-contrast, colorblind-friendly rainbow alternative
   */
  turbo: {
    name: 'Turbo',
    colors: [
      '#30123b',  // Deep indigo
      '#4662d7',  // Blue
      '#35aef4',  // Cyan
      '#1ae4b6',  // Teal
      '#72fe5e',  // Green
      '#c8ef34',  // Yellow-green
      '#faba39',  // Orange
      '#f66b19',  // Red-orange
      '#ca2a04',  // Red (dark enough to see)
      '#7a0403',  // Dark red
    ],
    logScale: false,
  },

  // ═══════════════════════════════════════════════════════════════════════════
  // DIVERGING COLORMAPS (for data with meaningful center point)
  // ═══════════════════════════════════════════════════════════════════════════

  /**
   * Velocity - Blue to White to Red-Orange
   * For velocity divergence (negative/positive)
   * Center is bright white for clear midpoint
   */
  velocity: {
    name: 'Velocity',
    colors: [
      '#0077b6',  // Deep blue
      '#00b4d8',  // Bright cyan
      '#48cae4',  // Light cyan
      '#90e0ef',  // Pale cyan
      '#caf0f8',  // Very pale cyan
      '#ffffff',  // White (center)
      '#ffccd5',  // Pale pink
      '#ff8fa3',  // Pink
      '#ff5f7e',  // Coral
      '#fb6107',  // Orange
      '#e63946',  // Warm red
    ],
    logScale: false,
  },

  /**
   * Divergent Sunset - Purple to Cream to Gold
   * Alternative diverging colormap, elegant
   */
  divergentSunset: {
    name: 'Divergent Sunset',
    colors: [
      '#5e2a84',  // Deep purple
      '#8b4d9e',  // Purple
      '#b57cc0',  // Light purple
      '#d8a8d8',  // Pale purple
      '#f5e6e8',  // Cream
      '#fff8dc',  // Light cream (center)
      '#ffe066',  // Yellow
      '#ffc233',  // Gold
      '#ff9f1c',  // Orange
      '#ff6b35',  // Red-orange
      '#d62828',  // Deep red
    ],
    logScale: false,
  },

  // ═══════════════════════════════════════════════════════════════════════════
  // SPECIALIZED COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════

  /**
   * Density - Optimized for log-scale density visualization
   * Avoids pure black (invisible) and pure red (hard on dark bg)
   */
  density: {
    name: 'Density',
    colors: [
      '#1a1a2e',  // Very dark blue (visible, not black)
      '#16213e',  // Navy
      '#1e3a5f',  // Dark blue
      '#3282b8',  // Medium blue
      '#00bcd4',  // Cyan
      '#00e5ff',  // Bright cyan
      '#76ff03',  // Neon green
      '#ffeb3b',  // Yellow
      '#ff9800',  // Orange
      '#ff5722',  // Deep orange (not red)
      '#ff1744',  // Red-pink (brighter than pure red)
    ],
    logScale: true,
  },

  /**
   * Black Hole - Dramatic colormap for IMBH simulations
   * Event horizon to accretion disk aesthetic
   */
  blackHole: {
    name: 'Black Hole',
    colors: [
      '#0f0f23',  // Near-black blue
      '#1c1c4d',  // Dark purple
      '#3b2e5a',  // Purple
      '#5c4a72',  // Dusty purple
      '#ff6b35',  // Orange (accretion disk inner)
      '#ffb627',  // Gold
      '#ffe97f',  // Pale gold
      '#fff9c4',  // Cream
      '#ffffff',  // White (hottest)
    ],
    logScale: true,
  },

  /**
   * Ice & Fire - Dramatic contrast for temperature
   * Cool ice tones to hot fire tones
   */
  iceFire: {
    name: 'Ice & Fire',
    colors: [
      '#00296b',  // Deep blue
      '#003f88',  // Royal blue
      '#00509d',  // Bright blue
      '#0077b6',  // Ocean blue
      '#00b4d8',  // Cyan
      '#90e0ef',  // Ice blue
      '#ffc300',  // Gold (transition)
      '#ff9500',  // Orange
      '#ff6700',  // Deep orange
      '#ff3d00',  // Red-orange
      '#d50000',  // Deep red (darker, visible)
    ],
    logScale: false,
  },

  /**
   * Spectral Bright - Rainbow with enhanced brightness for dark bg
   * Each color chosen for maximum visibility on near-black
   */
  spectralBright: {
    name: 'Spectral Bright',
    colors: [
      '#4a148c',  // Deep purple
      '#7b1fa2',  // Purple
      '#1565c0',  // Blue
      '#0288d1',  // Light blue
      '#00acc1',  // Cyan
      '#00897b',  // Teal
      '#43a047',  // Green
      '#7cb342',  // Light green
      '#c0ca33',  // Lime
      '#fdd835',  // Yellow
      '#ffb300',  // Amber
      '#fb8c00',  // Orange
      '#f4511e',  // Deep orange
    ],
    logScale: false,
  },

  /**
   * Mono Cyan - Single-hue gradient for clean scientific viz
   * High contrast, works well with additive blending
   */
  monoCyan: {
    name: 'Mono Cyan',
    colors: [
      '#001219',  // Near-black teal
      '#005f73',  // Dark teal
      '#0a9396',  // Teal
      '#40c9c9',  // Cyan
      '#94d2bd',  // Pale cyan-green
      '#e9d8a6',  // Cream
      '#ffffff',  // White
    ],
    logScale: false,
  },

  /**
   * Mono Gold - Warm single-hue for energy visualization
   * Elegant, high contrast on dark backgrounds
   */
  monoGold: {
    name: 'Mono Gold',
    colors: [
      '#1a1a0a',  // Very dark olive
      '#3d3d00',  // Dark gold
      '#6b6b00',  // Olive
      '#b8860b',  // Dark goldenrod
      '#daa520',  // Goldenrod
      '#ffd700',  // Gold
      '#ffeb3b',  // Yellow
      '#fff59d',  // Pale yellow
      '#ffffff',  // White
    ],
    logScale: false,
  },
}
