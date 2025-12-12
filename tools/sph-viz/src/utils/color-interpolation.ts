/**
 * Unified color interpolation utilities
 *
 * SSOT: Single Source of Truth for color mapping and interpolation
 * used across all visualization components
 */

// ============================================================================
// COLOR MAP DATA
// ============================================================================

/**
 * Pre-defined color maps as RGB arrays (0-1 range)
 * Each array is a gradient from low to high values
 */
export const COLOR_MAP_DATA: Record<string, number[][]> = {
  viridis: [
    [0.267, 0.004, 0.329],
    [0.282, 0.140, 0.458],
    [0.253, 0.265, 0.530],
    [0.206, 0.372, 0.553],
    [0.163, 0.471, 0.558],
    [0.127, 0.566, 0.551],
    [0.134, 0.658, 0.518],
    [0.266, 0.749, 0.441],
    [0.477, 0.821, 0.318],
    [0.741, 0.873, 0.150],
    [0.993, 0.906, 0.144]
  ],
  plasma: [
    [0.050, 0.030, 0.528],
    [0.295, 0.012, 0.615],
    [0.492, 0.012, 0.658],
    [0.665, 0.138, 0.618],
    [0.798, 0.280, 0.470],
    [0.899, 0.396, 0.301],
    [0.966, 0.530, 0.128],
    [0.988, 0.680, 0.063],
    [0.961, 0.850, 0.298],
    [0.940, 0.975, 0.131]
  ],
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
  turbo: [
    [0.190, 0.072, 0.232],
    [0.235, 0.318, 0.860],
    [0.137, 0.572, 0.996],
    [0.140, 0.780, 0.820],
    [0.376, 0.920, 0.512],
    [0.670, 0.979, 0.280],
    [0.924, 0.904, 0.145],
    [0.996, 0.724, 0.132],
    [0.994, 0.472, 0.122],
    [0.881, 0.200, 0.102],
    [0.528, 0.055, 0.052]
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
