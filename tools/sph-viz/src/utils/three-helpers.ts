/**
 * Shared Three.js helper utilities
 *
 * SSOT: Single Source of Truth for common Three.js operations
 * used across ParticleViewer3D, ParticleViewer3DImperative, GalacticMarkers, etc.
 */

import * as THREE from 'three'

// ============================================================================
// CIRCLE TEXTURE (for particle rendering)
// ============================================================================

/** Singleton circle texture for all particle systems */
let circleTexture: THREE.Texture | null = null

/**
 * Create a circular particle texture with radial gradient
 * Used for smooth round particle rendering
 */
export function createCircleTexture(size: number = 64): THREE.Texture {
  const canvas = document.createElement('canvas')
  canvas.width = size
  canvas.height = size
  const ctx = canvas.getContext('2d')!

  // Create radial gradient for smooth circular particle
  const gradient = ctx.createRadialGradient(
    size / 2, size / 2, 0,
    size / 2, size / 2, size / 2
  )
  gradient.addColorStop(0, 'rgba(255, 255, 255, 1)')
  gradient.addColorStop(0.3, 'rgba(255, 255, 255, 0.8)')
  gradient.addColorStop(0.7, 'rgba(255, 255, 255, 0.3)')
  gradient.addColorStop(1, 'rgba(255, 255, 255, 0)')

  ctx.fillStyle = gradient
  ctx.fillRect(0, 0, size, size)

  const texture = new THREE.CanvasTexture(canvas)
  texture.needsUpdate = true
  return texture
}

/**
 * Get the singleton circle texture (lazy initialization)
 * Use this to avoid creating multiple textures
 */
export function getCircleTexture(): THREE.Texture {
  if (!circleTexture) {
    circleTexture = createCircleTexture()
  }
  return circleTexture
}

// ============================================================================
// TEXT SPRITES (for labels)
// ============================================================================

export interface TextSpriteOptions {
  color?: string
  fontSize?: number
  fontWeight?: string
  fontFamily?: string
  backgroundColor?: string
  padding?: number
  scale?: [number, number, number]
  canvasSize?: number
  renderOrder?: number
}

/** Default options for text sprites - LARGER for better visibility */
export const DEFAULT_TEXT_SPRITE_OPTIONS: Required<TextSpriteOptions> = {
  color: '#ffffff',
  fontSize: 72,           // Increased from 48
  fontWeight: 'Bold',
  fontFamily: 'Arial, sans-serif',
  backgroundColor: '',
  padding: 10,
  scale: [12, 6, 1],      // Increased from [8, 4, 1]
  canvasSize: 512,
  renderOrder: 999,
}

/**
 * Create a text sprite for 3D labels
 * Renders text to a canvas and creates a sprite material
 */
export function createTextSprite(
  text: string,
  options: TextSpriteOptions = {}
): THREE.Sprite {
  const opts = { ...DEFAULT_TEXT_SPRITE_OPTIONS, ...options }
  const {
    color,
    fontSize,
    fontWeight,
    fontFamily,
    backgroundColor,
    scale,
    canvasSize,
    renderOrder,
  } = opts

  const canvas = document.createElement('canvas')
  canvas.width = canvasSize
  canvas.height = canvasSize / 2
  const ctx = canvas.getContext('2d')!

  // Background (if specified)
  if (backgroundColor) {
    ctx.fillStyle = backgroundColor
    ctx.fillRect(0, 0, canvasSize, canvasSize / 2)
  }

  // Text
  ctx.fillStyle = color
  ctx.font = `${fontWeight} ${fontSize}px ${fontFamily}`
  ctx.textAlign = 'center'
  ctx.textBaseline = 'middle'
  ctx.fillText(text, canvasSize / 2, canvasSize / 4)

  const texture = new THREE.CanvasTexture(canvas)
  texture.needsUpdate = true

  const material = new THREE.SpriteMaterial({
    map: texture,
    transparent: true,
    depthTest: false,
    depthWrite: false,
  })

  const sprite = new THREE.Sprite(material)
  sprite.scale.set(...scale)
  sprite.renderOrder = renderOrder

  return sprite
}

/**
 * Create axis label sprite with standardized styling
 */
export function createAxisLabel(
  text: string,
  axisColor: string,
  position: [number, number, number]
): THREE.Sprite {
  const sprite = createTextSprite(text, {
    color: axisColor,
    fontSize: 72,
    scale: [12, 6, 1],  // Larger scale for axis labels
  })
  sprite.position.set(...position)
  return sprite
}

/**
 * Create tick label sprite (smaller than axis labels)
 */
export function createTickLabel(
  text: string,
  position: [number, number, number]
): THREE.Sprite {
  const sprite = createTextSprite(text, {
    color: '#aaaaaa',
    fontSize: 48,
    scale: [6, 3, 1],  // Smaller scale for tick labels
  })
  sprite.position.set(...position)
  return sprite
}

// ============================================================================
// GEOMETRIC PRIMITIVES
// ============================================================================

export interface DashedCircleOptions {
  segments?: number
  dashSize?: number
  gapSize?: number
  opacity?: number
  normal?: THREE.Vector3
  center?: THREE.Vector3
}

/**
 * Create a dashed circle in 3D space
 */
export function createDashedCircle(
  radius: number,
  color: number,
  options: DashedCircleOptions = {}
): THREE.Line {
  const {
    segments = 64,
    dashSize = 0.5,
    gapSize = 0.25,
    opacity = 0.8,
    normal = new THREE.Vector3(0, 0, 1),
    center = new THREE.Vector3(0, 0, 0)
  } = options

  const points: THREE.Vector3[] = []

  // Create rotation to align circle with desired normal
  const defaultNormal = new THREE.Vector3(0, 0, 1)
  const quaternion = new THREE.Quaternion().setFromUnitVectors(
    defaultNormal,
    normal.clone().normalize()
  )

  for (let i = 0; i <= segments; i++) {
    const theta = (i / segments) * Math.PI * 2
    const point = new THREE.Vector3(
      radius * Math.cos(theta),
      radius * Math.sin(theta),
      0
    )
    point.applyQuaternion(quaternion)
    point.add(center)
    points.push(point)
  }

  const geometry = new THREE.BufferGeometry().setFromPoints(points)
  const material = new THREE.LineDashedMaterial({
    color,
    dashSize,
    gapSize,
    transparent: true,
    opacity,
  })

  const line = new THREE.Line(geometry, material)
  line.computeLineDistances()

  return line
}

export interface ArrowOptions {
  headLength?: number
  headWidth?: number
  shaftRadius?: number
  opacity?: number
}

/**
 * Create an arrow with cone head
 */
export function createArrow(
  from: THREE.Vector3,
  to: THREE.Vector3,
  color: number,
  options: ArrowOptions = {}
): THREE.Group {
  const {
    headLength = 1.0,
    headWidth = 0.5,
    shaftRadius = 0.08,
    opacity = 1.0
  } = options

  const group = new THREE.Group()
  const direction = to.clone().sub(from).normalize()
  const length = from.distanceTo(to)

  if (length <= headLength) return group

  // Shaft
  const shaftLength = length - headLength
  const shaftGeo = new THREE.CylinderGeometry(shaftRadius, shaftRadius, shaftLength, 12)
  const shaftMat = new THREE.MeshBasicMaterial({
    color,
    transparent: opacity < 1,
    opacity
  })
  const shaft = new THREE.Mesh(shaftGeo, shaftMat)

  // Position shaft at midpoint between from and the base of the head
  const shaftMidpoint = from.clone().add(direction.clone().multiplyScalar(shaftLength / 2))
  shaft.position.copy(shaftMidpoint)
  shaft.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
  group.add(shaft)

  // Head (cone)
  const headGeo = new THREE.ConeGeometry(headWidth, headLength, 16)
  const headMat = new THREE.MeshBasicMaterial({
    color,
    transparent: opacity < 1,
    opacity
  })
  const head = new THREE.Mesh(headGeo, headMat)

  // Position head at tip
  const headPosition = to.clone().sub(direction.clone().multiplyScalar(headLength / 2))
  head.position.copy(headPosition)
  head.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
  group.add(head)

  return group
}

// ============================================================================
// AXES HELPER WITH LABELS
// ============================================================================

export interface AxesWithLabelsOptions {
  size?: number
  showTickMarks?: boolean
  tickInterval?: number
  labelUnit?: string
}

/**
 * Create coordinate axes with labels and optional tick marks
 * Returns a group containing axes helper and label sprites
 */
export function createAxesWithLabels(
  options: AxesWithLabelsOptions = {}
): THREE.Group {
  const {
    size = 30,
    showTickMarks = true,
    tickInterval = 10,
    labelUnit = 'pc'
  } = options

  const group = new THREE.Group()
  group.name = 'axesWithLabels'

  // Axes helper (RGB lines)
  const axesHelper = new THREE.AxesHelper(size)
  group.add(axesHelper)

  // Axis labels - positioned beyond the axes
  const labelOffset = size + 4

  const xLabel = createAxisLabel(`X (${labelUnit})`, '#ff6666', [labelOffset, 0, 0])
  group.add(xLabel)

  const yLabel = createAxisLabel(`Y (${labelUnit})`, '#66ff66', [0, labelOffset, 0])
  group.add(yLabel)

  const zLabel = createAxisLabel(`Z (${labelUnit})`, '#6666ff', [0, 0, labelOffset])
  group.add(zLabel)

  // Tick marks
  if (showTickMarks) {
    const tickMaterial = new THREE.LineBasicMaterial({ color: 0x666666 })
    const tickSize = 0.8

    for (let i = -Math.floor(size / tickInterval) * tickInterval; i <= size; i += tickInterval) {
      if (i === 0) continue

      // X-axis ticks
      const xTickGeom = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(i, -tickSize, 0),
        new THREE.Vector3(i, tickSize, 0)
      ])
      group.add(new THREE.Line(xTickGeom, tickMaterial))

      // Y-axis ticks
      const yTickGeom = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-tickSize, i, 0),
        new THREE.Vector3(tickSize, i, 0)
      ])
      group.add(new THREE.Line(yTickGeom, tickMaterial))

      // Tick labels (every tickInterval)
      if (Math.abs(i) <= size) {
        const xTickLabel = createTickLabel(`${i}`, [i, -3, 0])
        group.add(xTickLabel)

        const yTickLabel = createTickLabel(`${i}`, [-3, i, 0])
        group.add(yTickLabel)
      }
    }
  }

  return group
}
