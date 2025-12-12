'use client'

import { useRef, useEffect, useMemo } from 'react'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'

/**
 * Orbital Geometry Configuration
 * Based on Oka et al. (2017) IMBH-Cloud Scattering Model
 */
export interface OrbitalGeometryConfig {
  // Black hole parameters
  bhMass: number              // M_BH in units of 10^5 M_sun
  
  // Cloud parameters
  cloudRadius: number         // R_cloud in pc
  cloudMass: number           // M_cloud in units of 1000 M_sun
  
  // Orbital parameters (in orbital plane)
  impactParameter: number     // b in pc
  pericentre: number          // r_peri in pc
  eccentricity: number        // e (>1 for hyperbolic)
  
  // Initial conditions
  initialPosition: [number, number, number]  // (X₀, Y₀, Z₀) in pc
  initialVelocity: [number, number, number]  // (vx, vy, vz) in km/s
  
  // Projection angles (for observer)
  inclination: number         // i in degrees (0° = face-on, 90° = edge-on)
  positionAngle: number       // PA in degrees (rotation on sky)
  
  // Observer reference
  lsrVelocity: number         // V_LSR in km/s (bulk motion toward observer)
  distanceToGC: number        // Distance from Earth to Galactic Center in kpc
}

// Default values from Oka et al. (2017) / CAT_OKA
const DEFAULT_CONFIG: OrbitalGeometryConfig = {
  bhMass: 1.0,                // 10^5 M_sun
  cloudRadius: 1.13,          // pc
  cloudMass: 1.0,             // 1000 M_sun
  impactParameter: 5.16,      // pc
  pericentre: 1.69,           // pc
  eccentricity: 1.24,         // hyperbolic
  initialPosition: [20.0, -5.17, 0],  // pc
  initialVelocity: [-9.35, 4.48, 0],  // km/s
  inclination: 70,            // degrees
  positionAngle: 41.6,        // degrees
  lsrVelocity: -120,          // km/s
  distanceToGC: 8.0,          // kpc
}

export interface OrbitalGeometryPanelProps {
  config?: Partial<OrbitalGeometryConfig>
  className?: string
  width?: number
  height?: number
  showLabels?: boolean
  viewMode?: 'orbital-plane' | 'observer' | '3d'
}

/**
 * Interactive panel showing orbital geometry for IMBH-cloud encounter
 * 
 * Features:
 * - Impact parameter visualization
 * - Pericentre marker
 * - Orbital plane with hyperbolic trajectory
 * - Line of sight direction
 * - Inclination and position angle indicators
 * - Earth/Observer reference frame
 * - LSR velocity vector
 */
export function OrbitalGeometryPanel({
  config: configOverride,
  className = '',
  width = 600,
  height = 500,
  showLabels = true,
  viewMode = '3d',
}: OrbitalGeometryPanelProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null)
  
  // Merge with defaults
  const config = useMemo(() => ({
    ...DEFAULT_CONFIG,
    ...configOverride,
  }), [configOverride])

  useEffect(() => {
    if (!containerRef.current) return

    const container = containerRef.current
    
    // Scene setup
    const scene = new THREE.Scene()
    scene.background = new THREE.Color(0x1a1a2e)
    
    // Camera
    const camera = new THREE.PerspectiveCamera(50, width / height, 0.1, 1000)
    
    // Position camera based on view mode
    switch (viewMode) {
      case 'orbital-plane':
        camera.position.set(0, 0, 50)  // Looking down at orbital plane
        break
      case 'observer':
        // Observer view (from Earth, along LoS)
        const inc = config.inclination * Math.PI / 180
        const pa = config.positionAngle * Math.PI / 180
        camera.position.set(
          50 * Math.sin(inc) * Math.cos(pa),
          50 * Math.sin(inc) * Math.sin(pa),
          50 * Math.cos(inc)
        )
        break
      default:
        camera.position.set(35, 25, 40)  // 3D overview
    }
    camera.lookAt(0, 0, 0)
    
    // Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true })
    renderer.setSize(width, height)
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    container.appendChild(renderer.domElement)
    rendererRef.current = renderer
    
    // Controls
    const controls = new OrbitControls(camera, renderer.domElement)
    controls.enableDamping = true
    controls.dampingFactor = 0.05
    controls.target.set(0, 0, 0)
    
    // ========== HELPER FUNCTIONS ==========
    
    const createTextSprite = (
      text: string, 
      color: string = '#ffffff',
      fontSize: number = 32,
      bgColor?: string
    ): THREE.Sprite => {
      const canvas = document.createElement('canvas')
      const size = 256
      canvas.width = size
      canvas.height = size / 2
      const ctx = canvas.getContext('2d')!
      
      if (bgColor) {
        ctx.fillStyle = bgColor
        ctx.fillRect(0, 0, size, size / 2)
      }
      
      ctx.fillStyle = color
      ctx.font = `Bold ${fontSize}px Arial`
      ctx.textAlign = 'center'
      ctx.textBaseline = 'middle'
      ctx.fillText(text, size / 2, size / 4)
      
      const texture = new THREE.CanvasTexture(canvas)
      const material = new THREE.SpriteMaterial({ map: texture, transparent: true })
      const sprite = new THREE.Sprite(material)
      sprite.scale.set(6, 3, 1)
      return sprite
    }
    
    const createArrow = (
      from: THREE.Vector3,
      to: THREE.Vector3,
      color: number,
      headLength: number = 0.5,
      headWidth: number = 0.3
    ): THREE.Group => {
      const group = new THREE.Group()
      const direction = to.clone().sub(from).normalize()
      const length = from.distanceTo(to)
      
      // Shaft
      const shaftGeo = new THREE.CylinderGeometry(0.05, 0.05, length - headLength, 8)
      const shaftMat = new THREE.MeshBasicMaterial({ color })
      const shaft = new THREE.Mesh(shaftGeo, shaftMat)
      shaft.position.copy(from.clone().lerp(to, (length - headLength) / 2 / length))
      shaft.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
      group.add(shaft)
      
      // Head
      const headGeo = new THREE.ConeGeometry(headWidth, headLength, 8)
      const headMat = new THREE.MeshBasicMaterial({ color })
      const head = new THREE.Mesh(headGeo, headMat)
      head.position.copy(to.clone().sub(direction.clone().multiplyScalar(headLength / 2)))
      head.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)
      group.add(head)
      
      return group
    }
    
    const createCircle = (
      radius: number,
      color: number,
      dashed: boolean = false,
      segments: number = 64
    ): THREE.Line => {
      const points: THREE.Vector3[] = []
      for (let i = 0; i <= segments; i++) {
        const theta = (i / segments) * Math.PI * 2
        points.push(new THREE.Vector3(radius * Math.cos(theta), radius * Math.sin(theta), 0))
      }
      const geo = new THREE.BufferGeometry().setFromPoints(points)
      const mat = dashed 
        ? new THREE.LineDashedMaterial({ color, dashSize: 0.3, gapSize: 0.15 })
        : new THREE.LineBasicMaterial({ color })
      const line = new THREE.Line(geo, mat)
      if (dashed) line.computeLineDistances()
      return line
    }
    
    // ========== COORDINATE AXES ==========
    
    // Orbital plane axes (X-Y)
    const axesGroup = new THREE.Group()
    
    // X-axis (toward pericentre direction in orbital plane)
    const xAxis = createArrow(
      new THREE.Vector3(-15, 0, 0),
      new THREE.Vector3(15, 0, 0),
      0xff4444
    )
    axesGroup.add(xAxis)
    
    // Y-axis (perpendicular in orbital plane)
    const yAxis = createArrow(
      new THREE.Vector3(0, -10, 0),
      new THREE.Vector3(0, 15, 0),
      0x44ff44
    )
    axesGroup.add(yAxis)
    
    // Z-axis (perpendicular to orbital plane)
    const zAxis = createArrow(
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, 0, 12),
      0x4488ff
    )
    axesGroup.add(zAxis)
    
    if (showLabels) {
      const xLabel = createTextSprite('X (pc)', '#ff6666')
      xLabel.position.set(16, 0, 0)
      axesGroup.add(xLabel)
      
      const yLabel = createTextSprite('Y (pc)', '#66ff66')
      yLabel.position.set(0, 16, 0)
      axesGroup.add(yLabel)
      
      const zLabel = createTextSprite('Z (⊥ plane)', '#6688ff')
      zLabel.position.set(0, 0, 13)
      axesGroup.add(zLabel)
    }
    
    scene.add(axesGroup)
    
    // ========== ORBITAL PLANE ==========
    
    // Semi-transparent orbital plane
    const planeGeo = new THREE.PlaneGeometry(50, 50)
    const planeMat = new THREE.MeshBasicMaterial({
      color: 0x3366ff,
      transparent: true,
      opacity: 0.1,
      side: THREE.DoubleSide,
    })
    const orbitalPlane = new THREE.Mesh(planeGeo, planeMat)
    scene.add(orbitalPlane)
    
    // Plane border
    const planeBorder = createCircle(20, 0x3366ff, true)
    scene.add(planeBorder)
    
    if (showLabels) {
      const planeLabel = createTextSprite('Orbital Plane', '#6699ff', 24)
      planeLabel.position.set(12, -12, 0.1)
      scene.add(planeLabel)
    }
    
    // ========== BLACK HOLE ==========
    
    const bhGroup = new THREE.Group()
    
    // Glow effect
    const glowGeo = new THREE.SphereGeometry(0.8, 32, 32)
    const glowMat = new THREE.MeshBasicMaterial({
      color: 0xff4444,
      transparent: true,
      opacity: 0.3,
    })
    bhGroup.add(new THREE.Mesh(glowGeo, glowMat))
    
    // Core
    const coreGeo = new THREE.SphereGeometry(0.3, 32, 32)
    const coreMat = new THREE.MeshBasicMaterial({ color: 0xffffff })
    bhGroup.add(new THREE.Mesh(coreGeo, coreMat))
    
    // Event horizon
    const horizonGeo = new THREE.SphereGeometry(0.15, 32, 32)
    const horizonMat = new THREE.MeshBasicMaterial({ color: 0x000000 })
    bhGroup.add(new THREE.Mesh(horizonGeo, horizonMat))
    
    scene.add(bhGroup)
    
    if (showLabels) {
      const bhLabel = createTextSprite('IMBH', '#ff8888', 28)
      bhLabel.position.set(0, -1.5, 0)
      scene.add(bhLabel)
      
      const bhMassLabel = createTextSprite(`10⁵ M☉`, '#ffaaaa', 20)
      bhMassLabel.position.set(0, -2.5, 0)
      scene.add(bhMassLabel)
    }
    
    // ========== HYPERBOLIC ORBIT ==========
    
    // Compute hyperbolic orbit trajectory
    const { impactParameter: b, pericentre: rp, eccentricity: e } = config
    const GM = 449.8  // pc·(km/s)² for 10^5 M_sun
    const a = rp / (e - 1)  // Semi-major axis (positive for calculation)
    const p = a * (e * e - 1)  // Semi-latus rectum
    
    // True anomaly limits (where r → ∞)
    const thetaMax = Math.PI - Math.acos(1 / e)
    
    const orbitPoints: THREE.Vector3[] = []
    const numPoints = 100
    for (let i = 0; i < numPoints; i++) {
      const theta = -thetaMax + 0.05 + (2 * thetaMax - 0.1) * (i / (numPoints - 1))
      const r = p / (1 + e * Math.cos(theta))
      
      // Rotate so pericentre is along +X
      const x = r * Math.cos(theta)
      const y = r * Math.sin(theta)
      
      if (Math.abs(x) < 30 && Math.abs(y) < 30) {
        orbitPoints.push(new THREE.Vector3(x, y, 0))
      }
    }
    
    const orbitGeo = new THREE.BufferGeometry().setFromPoints(orbitPoints)
    const orbitMat = new THREE.LineBasicMaterial({ color: 0x00ff88, linewidth: 2 })
    const orbitLine = new THREE.Line(orbitGeo, orbitMat)
    scene.add(orbitLine)
    
    // ========== PERICENTRE MARKER ==========
    
    const periMarker = new THREE.Group()
    
    // Pericentre point
    const periGeo = new THREE.SphereGeometry(0.2, 16, 16)
    const periMat = new THREE.MeshBasicMaterial({ color: 0xffff00 })
    periMarker.add(new THREE.Mesh(periGeo, periMat))
    
    // Line from BH to pericentre
    const periLineGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(rp, 0, 0)
    ])
    const periLineMat = new THREE.LineDashedMaterial({ 
      color: 0xffff00, 
      dashSize: 0.2, 
      gapSize: 0.1 
    })
    const periLine = new THREE.Line(periLineGeo, periLineMat)
    periLine.computeLineDistances()
    periMarker.add(periLine)
    
    periMarker.position.set(rp, 0, 0)
    scene.add(periMarker)
    
    if (showLabels) {
      const periLabel = createTextSprite('Pericentre', '#ffff66', 22)
      periLabel.position.set(rp + 0.5, 1.2, 0)
      scene.add(periLabel)
      
      const periValue = createTextSprite(`r_p = ${rp.toFixed(2)} pc`, '#ffff99', 18)
      periValue.position.set(rp + 0.5, 0.3, 0)
      scene.add(periValue)
    }
    
    // ========== IMPACT PARAMETER ==========
    
    // Incoming asymptote (dashed line)
    const asymptoteStart = new THREE.Vector3(-25, b, 0)
    const asymptoteEnd = new THREE.Vector3(5, b, 0)
    const asymptoteGeo = new THREE.BufferGeometry().setFromPoints([asymptoteStart, asymptoteEnd])
    const asymptoteMat = new THREE.LineDashedMaterial({
      color: 0xff8800,
      dashSize: 0.5,
      gapSize: 0.25,
    })
    const asymptoteLine = new THREE.Line(asymptoteGeo, asymptoteMat)
    asymptoteLine.computeLineDistances()
    scene.add(asymptoteLine)
    
    // Impact parameter line (from BH perpendicular to asymptote)
    const impactLineGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, b, 0)
    ])
    const impactLineMat = new THREE.LineBasicMaterial({ color: 0xff8800 })
    const impactLine = new THREE.Line(impactLineGeo, impactLineMat)
    scene.add(impactLine)
    
    // Impact parameter bracket
    const bracketSize = 0.3
    const bracketGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(-bracketSize, 0, 0),
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, bracketSize, 0),
    ])
    const bracketMat = new THREE.LineBasicMaterial({ color: 0xff8800 })
    scene.add(new THREE.Line(bracketGeo, bracketMat))
    
    if (showLabels) {
      const bLabel = createTextSprite('b (impact)', '#ffaa44', 20)
      bLabel.position.set(-1.5, b / 2, 0)
      scene.add(bLabel)
      
      const bValue = createTextSprite(`= ${b.toFixed(2)} pc`, '#ffcc66', 18)
      bValue.position.set(-1.5, b / 2 - 1, 0)
      scene.add(bValue)
      
      const asymLabel = createTextSprite('Incoming asymptote', '#ff9944', 18)
      asymLabel.position.set(-15, b + 1.5, 0)
      scene.add(asymLabel)
    }
    
    // ========== INITIAL CLOUD POSITION ==========
    
    const [x0, y0, z0] = config.initialPosition
    const [vx0, vy0, vz0] = config.initialVelocity
    
    // Cloud marker
    const cloudGroup = new THREE.Group()
    const cloudGeo = new THREE.SphereGeometry(config.cloudRadius, 32, 32)
    const cloudMat = new THREE.MeshBasicMaterial({
      color: 0x4488ff,
      transparent: true,
      opacity: 0.4,
    })
    cloudGroup.add(new THREE.Mesh(cloudGeo, cloudMat))
    
    // Cloud outline
    const cloudOutline = createCircle(config.cloudRadius, 0x66aaff)
    cloudGroup.add(cloudOutline)
    
    cloudGroup.position.set(x0, y0, z0)
    scene.add(cloudGroup)
    
    // Velocity vector
    const velScale = 0.8  // Scale velocity for visualization
    const velArrow = createArrow(
      new THREE.Vector3(x0, y0, z0),
      new THREE.Vector3(x0 + vx0 * velScale, y0 + vy0 * velScale, z0 + vz0 * velScale),
      0x44aaff,
      0.4,
      0.25
    )
    scene.add(velArrow)
    
    if (showLabels) {
      const cloudLabel = createTextSprite('Cloud (t=0)', '#66aaff', 22)
      cloudLabel.position.set(x0, y0 + config.cloudRadius + 1.5, z0)
      scene.add(cloudLabel)
      
      const posLabel = createTextSprite(
        `(${x0.toFixed(1)}, ${y0.toFixed(1)}) pc`,
        '#88ccff', 16
      )
      posLabel.position.set(x0, y0 - config.cloudRadius - 1.5, z0)
      scene.add(posLabel)
      
      const velLabel = createTextSprite(`v₀`, '#44aaff', 20)
      velLabel.position.set(
        x0 + vx0 * velScale + 1,
        y0 + vy0 * velScale + 0.5,
        z0
      )
      scene.add(velLabel)
    }
    
    // ========== LINE OF SIGHT & OBSERVER ==========
    
    const inc = config.inclination * Math.PI / 180
    const pa = config.positionAngle * Math.PI / 180
    
    // Line of sight direction (from observer to system)
    const losDir = new THREE.Vector3(
      Math.sin(inc) * Math.cos(pa),
      Math.sin(inc) * Math.sin(pa),
      Math.cos(inc)
    ).normalize()
    
    // LOS arrow (pointing toward observer)
    const losLength = 18
    const losStart = losDir.clone().multiplyScalar(-5)
    const losEnd = losDir.clone().multiplyScalar(losLength)
    const losArrow = createArrow(losStart, losEnd, 0xff44ff, 0.6, 0.35)
    scene.add(losArrow)
    
    // Earth marker at end of LOS
    const earthGroup = new THREE.Group()
    const earthGeo = new THREE.SphereGeometry(0.6, 32, 32)
    const earthMat = new THREE.MeshBasicMaterial({ color: 0x4488ff })
    earthGroup.add(new THREE.Mesh(earthGeo, earthMat))
    
    // Earth land patches (simplified)
    const landGeo = new THREE.SphereGeometry(0.62, 32, 32, 0, Math.PI, 0, Math.PI * 0.6)
    const landMat = new THREE.MeshBasicMaterial({ color: 0x44aa44 })
    const landMesh = new THREE.Mesh(landGeo, landMat)
    landMesh.rotation.x = Math.PI / 4
    earthGroup.add(landMesh)
    
    earthGroup.position.copy(losEnd)
    scene.add(earthGroup)
    
    if (showLabels) {
      const losLabel = createTextSprite('Line of Sight', '#ff88ff', 20)
      losLabel.position.copy(losDir.clone().multiplyScalar(losLength / 2).add(new THREE.Vector3(2, 0, 1)))
      scene.add(losLabel)
      
      const earthLabel = createTextSprite('Earth/ALMA', '#66aaff', 20)
      earthLabel.position.copy(losEnd.clone().add(new THREE.Vector3(2, 0, 1)))
      scene.add(earthLabel)
      
      const distLabel = createTextSprite(`d = ${config.distanceToGC} kpc`, '#88ccff', 16)
      distLabel.position.copy(losEnd.clone().add(new THREE.Vector3(2, 0, -0.5)))
      scene.add(distLabel)
    }
    
    // ========== INCLINATION ANGLE ARC ==========
    
    // Arc from Z-axis to LOS
    const incArcPoints: THREE.Vector3[] = []
    const arcRadius = 8
    for (let i = 0; i <= 20; i++) {
      const angle = (i / 20) * inc
      incArcPoints.push(new THREE.Vector3(
        arcRadius * Math.sin(angle) * Math.cos(pa),
        arcRadius * Math.sin(angle) * Math.sin(pa),
        arcRadius * Math.cos(angle)
      ))
    }
    const incArcGeo = new THREE.BufferGeometry().setFromPoints(incArcPoints)
    const incArcMat = new THREE.LineBasicMaterial({ color: 0xffaa00 })
    const incArc = new THREE.Line(incArcGeo, incArcMat)
    scene.add(incArc)
    
    if (showLabels) {
      const incLabel = createTextSprite(`i = ${config.inclination}°`, '#ffcc44', 22)
      const midAngle = inc / 2
      incLabel.position.set(
        (arcRadius + 2) * Math.sin(midAngle) * Math.cos(pa),
        (arcRadius + 2) * Math.sin(midAngle) * Math.sin(pa),
        (arcRadius + 2) * Math.cos(midAngle)
      )
      scene.add(incLabel)
    }
    
    // ========== POSITION ANGLE ARC ==========
    
    // Arc in the sky plane (perpendicular to Z-axis)
    const paArcPoints: THREE.Vector3[] = []
    const paRadius = 6
    for (let i = 0; i <= 20; i++) {
      const angle = (i / 20) * pa
      paArcPoints.push(new THREE.Vector3(
        paRadius * Math.cos(angle) * Math.sin(inc * 0.3),  // Projected
        paRadius * Math.sin(angle) * Math.sin(inc * 0.3),
        paRadius * Math.cos(inc * 0.3)
      ))
    }
    const paArcGeo = new THREE.BufferGeometry().setFromPoints(paArcPoints)
    const paArcMat = new THREE.LineDashedMaterial({ color: 0x44ffaa, dashSize: 0.2, gapSize: 0.1 })
    const paArc = new THREE.Line(paArcGeo, paArcMat)
    paArc.computeLineDistances()
    scene.add(paArc)
    
    if (showLabels) {
      const paLabel = createTextSprite(`PA = ${config.positionAngle}°`, '#66ffcc', 20)
      paLabel.position.set(paRadius + 2, 0, paRadius * Math.cos(inc * 0.3) + 1)
      scene.add(paLabel)
    }
    
    // ========== LSR VELOCITY ==========
    
    // LSR velocity vector (toward observer, negative = approaching)
    const lsrArrow = createArrow(
      new THREE.Vector3(-12, -8, 0),
      new THREE.Vector3(-12 + losDir.x * 4, -8 + losDir.y * 4, losDir.z * 4),
      0xff6688,
      0.3,
      0.2
    )
    scene.add(lsrArrow)
    
    if (showLabels) {
      const lsrLabel = createTextSprite('V_LSR', '#ff8899', 20)
      lsrLabel.position.set(-12, -9.5, 0)
      scene.add(lsrLabel)
      
      const lsrValue = createTextSprite(`= ${config.lsrVelocity} km/s`, '#ffaabb', 16)
      lsrValue.position.set(-12, -10.8, 0)
      scene.add(lsrValue)
    }
    
    // ========== LEGEND BOX ==========
    
    // Create a legend in 3D space
    const legendGroup = new THREE.Group()
    legendGroup.position.set(-22, 12, 0)
    
    const legendBg = new THREE.Mesh(
      new THREE.PlaneGeometry(12, 14),
      new THREE.MeshBasicMaterial({ color: 0x1a1a2e, transparent: true, opacity: 0.9 })
    )
    legendGroup.add(legendBg)
    
    const legendTitle = createTextSprite('LEGEND', '#ffffff', 28)
    legendTitle.position.set(0, 5.5, 0.1)
    legendTitle.scale.set(4, 2, 1)
    legendGroup.add(legendTitle)
    
    // Legend items
    const items = [
      { color: '#00ff88', text: 'Hyperbolic orbit' },
      { color: '#ff8800', text: 'Impact parameter b' },
      { color: '#ffff00', text: 'Pericentre r_p' },
      { color: '#ff44ff', text: 'Line of sight' },
      { color: '#ffaa00', text: 'Inclination i' },
      { color: '#44ffaa', text: 'Position angle PA' },
      { color: '#4488ff', text: 'Cloud position' },
    ]
    
    items.forEach((item, i) => {
      const dot = new THREE.Mesh(
        new THREE.SphereGeometry(0.2, 8, 8),
        new THREE.MeshBasicMaterial({ color: item.color })
      )
      dot.position.set(-4.5, 4 - i * 1.3, 0.1)
      legendGroup.add(dot)
      
      const label = createTextSprite(item.text, item.color, 18)
      label.position.set(1, 4 - i * 1.3, 0.1)
      label.scale.set(4, 2, 1)
      legendGroup.add(label)
    })
    
    scene.add(legendGroup)
    
    // ========== RENDER LOOP ==========
    
    let animationId: number
    const animate = () => {
      animationId = requestAnimationFrame(animate)
      controls.update()
      renderer.render(scene, camera)
    }
    animate()
    
    // Cleanup
    return () => {
      cancelAnimationFrame(animationId)
      renderer.dispose()
      container.removeChild(renderer.domElement)
    }
  }, [config, width, height, showLabels, viewMode])

  return (
    <div className={`relative ${className}`}>
      <div 
        ref={containerRef} 
        style={{ width, height }}
        className="rounded-lg overflow-hidden"
      />
      
      {/* Info panel */}
      <div className="absolute top-2 left-2 bg-black/70 text-white text-xs p-2 rounded">
        <div className="font-bold mb-1">Orbital Geometry</div>
        <div>Oka et al. (2017) / CAT_OKA</div>
      </div>
      
      {/* Parameter summary */}
      <div className="absolute bottom-2 right-2 bg-black/70 text-white text-xs p-2 rounded max-w-xs">
        <table className="w-full">
          <tbody>
            <tr>
              <td className="text-yellow-400">r_peri:</td>
              <td className="text-right">{config.pericentre.toFixed(2)} pc</td>
            </tr>
            <tr>
              <td className="text-orange-400">b:</td>
              <td className="text-right">{config.impactParameter.toFixed(2)} pc</td>
            </tr>
            <tr>
              <td className="text-green-400">e:</td>
              <td className="text-right">{config.eccentricity.toFixed(3)}</td>
            </tr>
            <tr>
              <td className="text-purple-400">i:</td>
              <td className="text-right">{config.inclination}°</td>
            </tr>
            <tr>
              <td className="text-teal-400">PA:</td>
              <td className="text-right">{config.positionAngle}°</td>
            </tr>
            <tr>
              <td className="text-pink-400">V_LSR:</td>
              <td className="text-right">{config.lsrVelocity} km/s</td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  )
}

export default OrbitalGeometryPanel
