/**
 * Unified color interpolation utilities
 *
 * SSOT: Single Source of Truth for color mapping and interpolation
 * used across all visualization components
 * 
 * Color Design Principles:
 * - Optimized for dark backgrounds (bg-gray-900, #0a0a0f canvas)
 * - High visibility and contrast
 * - Avoids pure red on black (hard to see)
 * - Avoids low-saturation colors that disappear on dark backgrounds
 */

// ============================================================================
// COLOR MAP DATA
// ============================================================================

/**
 * Pre-defined color maps as RGB arrays (0-1 range)
 * Each array is a gradient from low to high values
 * 
 * All colormaps are designed for excellent visibility on dark backgrounds
 */
export const COLOR_MAP_DATA: Record<string, number[][]> = {
  // ═══════════════════════════════════════════════════════════════════════════
  // SEQUENTIAL COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════

  /**
   * Cosmic Dawn - Deep space theme
   * Purple → Blue → Cyan → Gold → White
   */
  cosmicDawn: [
    [0.102, 0.020, 0.200],  // #1a0533
    [0.176, 0.106, 0.412],  // #2d1b69
    [0.239, 0.310, 0.780],  // #3d4fc7
    [0.000, 0.706, 0.847],  // #00b4d8
    [0.282, 0.792, 0.894],  // #48cae4
    [0.565, 0.878, 0.937],  // #90e0ef
    [1.000, 0.820, 0.400],  // #ffd166
    [1.000, 0.922, 0.600],  // #ffeb99
    [1.000, 1.000, 1.000],  // #ffffff
  ],

  /**
   * Nebula Fire - Warm energy theme
   * Magenta → Pink → Orange → Yellow → White
   */
  nebulaFire: [
    [0.176, 0.039, 0.192],  // #2d0a31
    [0.361, 0.067, 0.345],  // #5c1158
    [0.608, 0.176, 0.498],  // #9b2d7f
    [0.839, 0.259, 0.573],  // #d64292
    [0.957, 0.431, 0.431],  // #f46e6e
    [1.000, 0.624, 0.263],  // #ff9f43
    [1.000, 0.788, 0.235],  // #ffc93c
    [1.000, 0.945, 0.463],  // #fff176
    [1.000, 1.000, 1.000],  // #ffffff
  ],

  /**
   * Ocean Depths - Cool pressure theme
   * Navy → Blue → Cyan → Mint → Cream
   */
  oceanDepths: [
    [0.051, 0.106, 0.165],  // #0d1b2a
    [0.106, 0.227, 0.294],  // #1b3a4b
    [0.078, 0.302, 0.494],  // #144d7e
    [0.118, 0.533, 0.898],  // #1e88e5
    [0.259, 0.647, 0.961],  // #42a5f5
    [0.502, 0.871, 0.918],  // #80deea
    [0.655, 1.000, 0.922],  // #a7ffeb
    [0.878, 0.969, 0.980],  // #e0f7fa
    [1.000, 0.992, 0.906],  // #fffde7
  ],

  /**
   * Aurora - Vivid multi-hue
   * Purple → Blue → Cyan → Green → Yellow
   */
  aurora: [
    [0.227, 0.047, 0.639],  // #3a0ca3
    [0.263, 0.380, 0.933],  // #4361ee
    [0.298, 0.788, 0.941],  // #4cc9f0
    [0.024, 0.839, 0.627],  // #06d6a0
    [0.322, 0.718, 0.533],  // #52b788
    [0.600, 0.851, 0.549],  // #99d98c
    [0.851, 0.929, 0.573],  // #d9ed92
    [0.988, 0.965, 0.741],  // #fcf6bd
    [1.000, 0.953, 0.690],  // #fff3b0
  ],

  // ═══════════════════════════════════════════════════════════════════════════
  // SCIENTIFIC COLORMAPS (perceptually uniform)
  // ═══════════════════════════════════════════════════════════════════════════

  viridis: [
    [0.267, 0.004, 0.329],  // #440154
    [0.282, 0.157, 0.471],  // #482878
    [0.243, 0.286, 0.537],  // #3e4989
    [0.192, 0.408, 0.557],  // #31688e
    [0.149, 0.510, 0.557],  // #26828e
    [0.122, 0.620, 0.537],  // #1f9e89
    [0.208, 0.718, 0.475],  // #35b779
    [0.431, 0.808, 0.345],  // #6ece58
    [0.710, 0.871, 0.169],  // #b5de2b
    [0.992, 0.906, 0.145],  // #fde725
  ],

  plasma: [
    [0.051, 0.031, 0.529],  // #0d0887
    [0.274, 0.012, 0.624],  // #46039f
    [0.447, 0.004, 0.659],  // #7201a8
    [0.612, 0.090, 0.620],  // #9c179e
    [0.741, 0.216, 0.525],  // #bd3786
    [0.847, 0.337, 0.420],  // #d8576b
    [0.929, 0.475, 0.325],  // #ed7953
    [0.984, 0.624, 0.227],  // #fb9f3a
    [0.992, 0.792, 0.149],  // #fdca26
    [0.941, 0.976, 0.129],  // #f0f921
  ],

  turbo: [
    [0.188, 0.071, 0.231],  // #30123b
    [0.275, 0.384, 0.843],  // #4662d7
    [0.208, 0.682, 0.957],  // #35aef4
    [0.102, 0.894, 0.714],  // #1ae4b6
    [0.447, 0.996, 0.369],  // #72fe5e
    [0.784, 0.937, 0.204],  // #c8ef34
    [0.980, 0.729, 0.224],  // #faba39
    [0.965, 0.420, 0.098],  // #f66b19
    [0.792, 0.165, 0.016],  // #ca2a04
    [0.478, 0.016, 0.012],  // #7a0403
  ],

  // ═══════════════════════════════════════════════════════════════════════════
  // DIVERGING COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════

  /**
   * Velocity - Blue to White to Red-Orange
   */
  velocity: [
    [0.000, 0.467, 0.714],  // #0077b6
    [0.000, 0.706, 0.847],  // #00b4d8
    [0.282, 0.792, 0.894],  // #48cae4
    [0.565, 0.878, 0.937],  // #90e0ef
    [0.792, 0.941, 0.973],  // #caf0f8
    [1.000, 1.000, 1.000],  // #ffffff
    [1.000, 0.800, 0.835],  // #ffccd5
    [1.000, 0.561, 0.639],  // #ff8fa3
    [1.000, 0.373, 0.494],  // #ff5f7e
    [0.984, 0.380, 0.027],  // #fb6107
    [0.902, 0.224, 0.275],  // #e63946
  ],

  /**
   * Divergent Sunset - Purple to Cream to Gold
   */
  divergentSunset: [
    [0.369, 0.165, 0.518],  // #5e2a84
    [0.545, 0.302, 0.620],  // #8b4d9e
    [0.710, 0.486, 0.753],  // #b57cc0
    [0.847, 0.659, 0.847],  // #d8a8d8
    [0.961, 0.902, 0.910],  // #f5e6e8
    [1.000, 0.973, 0.863],  // #fff8dc
    [1.000, 0.875, 0.400],  // #ffe066
    [1.000, 0.761, 0.200],  // #ffc233
    [1.000, 0.624, 0.110],  // #ff9f1c
    [1.000, 0.420, 0.208],  // #ff6b35
    [0.839, 0.157, 0.157],  // #d62828
  ],

  // ═══════════════════════════════════════════════════════════════════════════
  // SPECIALIZED COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════

  /**
   * Density - Optimized for log-scale density
   */
  density: [
    [0.102, 0.102, 0.180],  // #1a1a2e
    [0.086, 0.129, 0.243],  // #16213e
    [0.118, 0.227, 0.373],  // #1e3a5f
    [0.196, 0.510, 0.722],  // #3282b8
    [0.000, 0.737, 0.831],  // #00bcd4
    [0.000, 0.898, 1.000],  // #00e5ff
    [0.463, 1.000, 0.012],  // #76ff03
    [1.000, 0.922, 0.231],  // #ffeb3b
    [1.000, 0.596, 0.000],  // #ff9800
    [1.000, 0.341, 0.133],  // #ff5722
    [1.000, 0.090, 0.267],  // #ff1744
  ],

  /**
   * Black Hole - Dramatic accretion disk theme
   */
  blackHole: [
    [0.059, 0.059, 0.137],  // #0f0f23
    [0.110, 0.110, 0.302],  // #1c1c4d
    [0.231, 0.180, 0.353],  // #3b2e5a
    [0.361, 0.290, 0.447],  // #5c4a72
    [1.000, 0.420, 0.208],  // #ff6b35
    [1.000, 0.714, 0.153],  // #ffb627
    [1.000, 0.914, 0.498],  // #ffe97f
    [1.000, 0.976, 0.769],  // #fff9c4
    [1.000, 1.000, 1.000],  // #ffffff
  ],

  /**
   * Ice & Fire - Temperature contrast
   */
  iceFire: [
    [0.000, 0.161, 0.420],  // #00296b
    [0.000, 0.247, 0.533],  // #003f88
    [0.000, 0.314, 0.616],  // #00509d
    [0.000, 0.467, 0.714],  // #0077b6
    [0.000, 0.706, 0.847],  // #00b4d8
    [0.565, 0.878, 0.937],  // #90e0ef
    [1.000, 0.761, 0.000],  // #ffc300
    [1.000, 0.584, 0.000],  // #ff9500
    [1.000, 0.404, 0.000],  // #ff6700
    [1.000, 0.239, 0.000],  // #ff3d00
    [0.835, 0.000, 0.000],  // #d50000
  ],

  /**
   * Spectral Bright - Enhanced rainbow
   */
  spectralBright: [
    [0.290, 0.078, 0.549],  // #4a148c
    [0.482, 0.122, 0.635],  // #7b1fa2
    [0.082, 0.396, 0.753],  // #1565c0
    [0.008, 0.533, 0.820],  // #0288d1
    [0.000, 0.675, 0.757],  // #00acc1
    [0.000, 0.537, 0.482],  // #00897b
    [0.263, 0.627, 0.278],  // #43a047
    [0.486, 0.702, 0.259],  // #7cb342
    [0.753, 0.792, 0.200],  // #c0ca33
    [0.992, 0.847, 0.208],  // #fdd835
    [1.000, 0.702, 0.000],  // #ffb300
    [0.984, 0.549, 0.000],  // #fb8c00
    [0.957, 0.318, 0.118],  // #f4511e
  ],

  /**
   * Mono Cyan - Clean single-hue
   */
  monoCyan: [
    [0.000, 0.071, 0.098],  // #001219
    [0.000, 0.373, 0.451],  // #005f73
    [0.039, 0.576, 0.588],  // #0a9396
    [0.251, 0.788, 0.788],  // #40c9c9
    [0.580, 0.824, 0.741],  // #94d2bd
    [0.914, 0.847, 0.651],  // #e9d8a6
    [1.000, 1.000, 1.000],  // #ffffff
  ],

  /**
   * Mono Gold - Warm single-hue
   */
  monoGold: [
    [0.102, 0.102, 0.039],  // #1a1a0a
    [0.239, 0.239, 0.000],  // #3d3d00
    [0.420, 0.420, 0.000],  // #6b6b00
    [0.722, 0.525, 0.043],  // #b8860b
    [0.855, 0.647, 0.125],  // #daa520
    [1.000, 0.843, 0.000],  // #ffd700
    [1.000, 0.922, 0.231],  // #ffeb3b
    [1.000, 0.961, 0.616],  // #fff59d
    [1.000, 1.000, 1.000],  // #ffffff
  ],

  // ═══════════════════════════════════════════════════════════════════════════
  // LEGACY COLORMAPS (kept for backward compatibility)
  // ═══════════════════════════════════════════════════════════════════════════

  inferno: [
    [0.001, 0.000, 0.014],
    [0.122, 0.047, 0.282],
    [0.304, 0.063, 0.420],
    [0.499, 0.086, 0.397],
    [0.680, 0.144, 0.295],
    [0.833, 0.253, 0.160],
    [0.937, 0.405, 0.049],
    [0.981, 0.588, 0.068],
    [0.987, 0.772, 0.264],
    [0.988, 0.998, 0.645]
  ],

  magma: [
    [0.001, 0.000, 0.014],
    [0.104, 0.047, 0.258],
    [0.259, 0.050, 0.408],
    [0.427, 0.079, 0.430],
    [0.575, 0.134, 0.397],
    [0.716, 0.215, 0.345],
    [0.848, 0.343, 0.331],
    [0.937, 0.517, 0.388],
    [0.973, 0.699, 0.530],
    [0.988, 0.998, 0.645]
  ],

  coolwarm: [
    [0.230, 0.299, 0.754],
    [0.411, 0.484, 0.845],
    [0.593, 0.669, 0.927],
    [0.775, 0.817, 0.964],
    [0.900, 0.900, 0.900],
    [0.964, 0.775, 0.692],
    [0.927, 0.593, 0.476],
    [0.845, 0.411, 0.299],
    [0.754, 0.230, 0.173]
  ]
}

// ============================================================================
// RGB COLOR TYPE
// ============================================================================

export interface RGBColor {
  r: number
  g: number
  b: number
}

// ============================================================================
// HEX TO RGB CONVERSION (with caching)
// ============================================================================

const hexToRgbCache = new Map<string, RGBColor>()

/**
 * Convert hex color string to RGB (0-1 range)
 */
export function hexToRgb(hex: string): RGBColor {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex)
  if (!result) return { r: 1, g: 1, b: 1 }
  return {
    r: parseInt(result[1], 16) / 255,
    g: parseInt(result[2], 16) / 255,
    b: parseInt(result[3], 16) / 255,
  }
}

/**
 * Convert hex color string to RGB with caching
 */
export function hexToRgbCached(hex: string): RGBColor {
  let cached = hexToRgbCache.get(hex)
  if (!cached) {
    cached = hexToRgb(hex)
    hexToRgbCache.set(hex, cached)
  }
  return cached
}

/**
 * Convert RGB (0-1 range) to hex string
 */
export function rgbToHex(r: number, g: number, b: number): string {
  const toHex = (v: number) => {
    const h = Math.round(Math.max(0, Math.min(255, v * 255))).toString(16)
    return h.length === 1 ? '0' + h : h
  }
  return `#${toHex(r)}${toHex(g)}${toHex(b)}`
}

// ============================================================================
// COLOR MAP SAMPLING
// ============================================================================

/**
 * Sample a color from a named color map at position t (0-1)
 * Uses pre-defined RGB data for fast interpolation
 */
export function sampleColorMap(
  mapName: string,
  t: number
): [number, number, number] {
  const map = COLOR_MAP_DATA[mapName] || COLOR_MAP_DATA.viridis
  t = Math.max(0, Math.min(1, t))

  const idx = t * (map.length - 1)
  const i = Math.floor(idx)
  const f = idx - i

  if (i >= map.length - 1) {
    return map[map.length - 1] as [number, number, number]
  }

  return [
    map[i][0] + f * (map[i + 1][0] - map[i][0]),
    map[i][1] + f * (map[i + 1][1] - map[i][1]),
    map[i][2] + f * (map[i + 1][2] - map[i][2])
  ]
}

/**
 * Sample a color from a hex color array at position t (0-1)
 * Used when color map is provided as hex strings
 */
export function interpolateColorHex(
  colors: string[],
  t: number
): RGBColor {
  if (colors.length === 0) return { r: 1, g: 1, b: 1 }
  if (colors.length === 1) return hexToRgbCached(colors[0])

  t = Math.max(0, Math.min(1, t))
  const index = t * (colors.length - 1)
  const lower = Math.floor(index)
  const upper = Math.min(lower + 1, colors.length - 1)
  const localT = index - lower

  const c1 = hexToRgbCached(colors[lower])
  const c2 = hexToRgbCached(colors[upper])

  return {
    r: c1.r + (c2.r - c1.r) * localT,
    g: c1.g + (c2.g - c1.g) * localT,
    b: c1.b + (c2.b - c1.b) * localT,
  }
}

// ============================================================================
// VALUE NORMALIZATION
// ============================================================================

/**
 * Normalize a value to [0, 1] range
 * Supports both linear and logarithmic scaling
 */
export function normalizeValue(
  value: number,
  min: number,
  max: number,
  logScale: boolean = false
): number {
  if (!isFinite(value)) value = min

  let t: number
  if (logScale && min > 0) {
    const logMin = Math.log10(min)
    const logMax = Math.log10(max)
    const logRange = logMax - logMin || 1
    t = (Math.log10(Math.max(value, min)) - logMin) / logRange
  } else {
    const range = max - min || 1
    t = (value - min) / range
  }

  return Math.max(0, Math.min(1, t))
}

/**
 * Compute percentile values from an array
 */
export function percentile(arr: number[], p: number): number {
  if (arr.length === 0) return 0
  const sorted = arr.slice().sort((a, b) => a - b)
  const idx = Math.floor(sorted.length * p / 100)
  return sorted[Math.min(idx, sorted.length - 1)]
}

/**
 * Compute color range (min, max) from data using percentiles
 * to exclude outliers
 */
export function computeColorRange(
  data: Float32Array | number[],
  lowPercentile: number = 1,
  highPercentile: number = 99
): [number, number] {
  const validData: number[] = []

  for (let i = 0; i < data.length; i++) {
    const v = data[i]
    if (isFinite(v) && v > 0) {
      validData.push(v)
    }
  }

  if (validData.length === 0) {
    return [0.001, 1.0]
  }

  const min = percentile(validData, lowPercentile)
  const max = percentile(validData, highPercentile)

  return [min, max === min ? min * 1.1 + 0.001 : max]
}
