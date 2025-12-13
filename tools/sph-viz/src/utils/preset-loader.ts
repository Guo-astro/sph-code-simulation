/**
 * Preset Configuration Loader
 *
 * Loads simulation configuration from preset JSON files stored in the
 * SPH repository. Extracts IMBH physics parameters for visualization.
 */

import type { IMBHPhysicsConfig, Vector3 } from '~/types/sph'

/**
 * Raw simulation preset JSON structure
 * (matches the format in simulations/astrophysics/imbh_cloud/config/presets/)
 */
export interface SimulationPreset {
  name: string
  description?: string
  physics_summary?: {
    cloud_mass_Msun?: number
    cloud_radius_pc?: number
    bh_mass_Msun?: number
    impact_parameter_pc?: number
    approach_velocity_kms?: number
    tidal_radius_pc?: number
    b_over_r_tidal?: number
    regime?: string
    thermal_physics?: string
    particles?: number
  }
  units?: {
    type?: string
    length_pc?: number
    mass_1e3Msun?: number
    velocity_kms?: number
  }
  initialCondition?: {
    transform?: {
      translate?: [number, number, number]
      velocity_boost?: [number, number, number]
    }
  }
  imbh_parameters?: {
    enabled?: boolean
    M_BH?: number  // Mass in M_sun
    BH_initial_position?: [number, number, number]
    BH_initial_velocity?: [number, number, number]
    softening_epsilon?: number
    cloud_initial_position?: [number, number, number]
    cloud_initial_velocity?: [number, number, number]
    comment?: string
  }
  externalForces?: {
    pointMassBH?: {
      enabled?: boolean
      mass?: number  // In code units (1000 M_sun)
      position?: [number, number, number]
      softening?: number
    }
  }
}

/**
 * Compute orbital parameters from initial conditions
 */
function computeOrbitalParams(
  cloudPos: Vector3,
  cloudVel: Vector3,
  bhMass: number,  // In code units (1000 M_sun)
  G: number = 1.0
): { pericentre: number; eccentricity: number } {
  // Position and velocity magnitudes
  const r = Math.sqrt(cloudPos[0]**2 + cloudPos[1]**2 + cloudPos[2]**2)
  const v = Math.sqrt(cloudVel[0]**2 + cloudVel[1]**2 + cloudVel[2]**2)
  
  // Specific orbital energy: E = v²/2 - GM/r
  const mu = G * bhMass
  const specificEnergy = (v * v) / 2 - mu / r
  
  // Semi-major axis (negative for hyperbolic)
  const a = -mu / (2 * specificEnergy)
  
  // Angular momentum vector: L = r × v
  const Lx = cloudPos[1] * cloudVel[2] - cloudPos[2] * cloudVel[1]
  const Ly = cloudPos[2] * cloudVel[0] - cloudPos[0] * cloudVel[2]
  const Lz = cloudPos[0] * cloudVel[1] - cloudPos[1] * cloudVel[0]
  const L = Math.sqrt(Lx**2 + Ly**2 + Lz**2)
  
  // Eccentricity: e = sqrt(1 + 2EL²/(G²M²m))
  const eccentricity = Math.sqrt(1 + (2 * specificEnergy * L * L) / (mu * mu))
  
  // Pericentre distance: r_p = a(1 - e) for elliptic, a(e - 1) for hyperbolic
  const pericentre = eccentricity > 1 
    ? Math.abs(a) * (eccentricity - 1)  // Hyperbolic
    : Math.abs(a) * (1 - eccentricity)  // Elliptic
  
  return { pericentre, eccentricity }
}

/**
 * Convert a simulation preset JSON to IMBHPhysicsConfig
 */
export function presetToIMBHConfig(preset: SimulationPreset): IMBHPhysicsConfig | null {
  // Check if IMBH physics is enabled
  const imbhParams = preset.imbh_parameters
  const externalForce = preset.externalForces?.pointMassBH
  
  // Either imbh_parameters or externalForces.pointMassBH should be present
  if (!imbhParams?.enabled && !externalForce?.enabled) {
    return null
  }
  
  // Get BH mass - prefer code units from externalForces
  const bhMassCodeUnits = externalForce?.mass ?? (imbhParams?.M_BH ? imbhParams.M_BH / 1000 : 100)
  
  // Get BH position
  const bhPosition: Vector3 = externalForce?.position as Vector3 
    ?? imbhParams?.BH_initial_position as Vector3 
    ?? [0, 0, 0]
  
  // Get cloud initial conditions
  const cloudPos: Vector3 = imbhParams?.cloud_initial_position as Vector3 
    ?? preset.initialCondition?.transform?.translate as Vector3 
    ?? [-20, 5.17, 0]
  
  const cloudVel: Vector3 = imbhParams?.cloud_initial_velocity as Vector3 
    ?? preset.initialCondition?.transform?.velocity_boost as Vector3 
    ?? [10, 0, 0]
  
  // Get cloud and tidal parameters from physics_summary
  const cloudMass = preset.physics_summary?.cloud_mass_Msun 
    ? preset.physics_summary.cloud_mass_Msun / 1000  // Convert to code units
    : 1
  const cloudRadius = preset.physics_summary?.cloud_radius_pc ?? 1.13
  const tidalRadius = preset.physics_summary?.tidal_radius_pc ?? 5.24
  const impactParameter = preset.physics_summary?.impact_parameter_pc ?? Math.abs(cloudPos[1])
  
  // Compute orbital parameters from initial conditions
  const { pericentre, eccentricity } = computeOrbitalParams(cloudPos, cloudVel, bhMassCodeUnits)
  
  // Time unit (standard IMBH_ENCOUNTER units)
  const timeUnit = 0.978  // Myr
  
  return {
    enabled: true,
    bhPosition,
    bhMass: bhMassCodeUnits,
    cloudInitialPosition: cloudPos,
    cloudInitialVelocity: cloudVel,
    cloudMass,
    cloudRadius,
    tidalRadius,
    impactParameter,
    pericentre,
    eccentricity,
    timeUnit,
    // Optional galactic parameters (can be overridden)
    inclination: 70,
    positionAngle: 41.6,
    lsrVelocity: -120,
  }
}

/**
 * Load and parse a preset JSON file
 * 
 * @param presetPath - Path to the preset JSON file (relative to API base)
 * @param apiBase - Base URL for the API (default: http://localhost:3001)
 * @returns Promise<IMBHPhysicsConfig | null>
 */
export async function loadPresetConfig(
  presetPath: string,
  apiBase: string = 'http://localhost:3001'
): Promise<IMBHPhysicsConfig | null> {
  try {
    // Fetch the preset file through the API
    const response = await fetch(`${apiBase}/api/preset?path=${encodeURIComponent(presetPath)}`)
    
    if (!response.ok) {
      console.warn(`[loadPresetConfig] Failed to load preset: ${response.statusText}`)
      return null
    }
    
    const preset: SimulationPreset = await response.json()
    return presetToIMBHConfig(preset)
  } catch (error) {
    console.error('[loadPresetConfig] Error loading preset:', error)
    return null
  }
}

/**
 * Default IMBH physics configuration (fallback)
 * Based on CAT_OKA/A_61k preset values (Oka et al. 2017)
 */
export const DEFAULT_IMBH_CONFIG: IMBHPhysicsConfig = {
  enabled: true,
  bhPosition: [0, 0, 0],
  bhMass: 100,  // 10^5 M_sun in code units (1000 M_sun)
  cloudInitialPosition: [20.0, -5.17, 0],
  cloudInitialVelocity: [-10.18, 5.05, 0],
  cloudMass: 1,  // 1000 M_sun in code units
  cloudRadius: 1.13,
  tidalRadius: 5.24,
  impactParameter: 5.17,
  pericentre: 2.217,
  eccentricity: 1.4504,
  timeUnit: 0.978,
  inclination: 70,
  positionAngle: 41.6,
  lsrVelocity: -120,
}

/**
 * List of available preset scenarios for the UI
 */
export const PRESET_SCENARIOS = [
  {
    id: 'Mc1e3_Mbh1e5_b1p5_v10_adiabatic_gsph',
    name: 'Strong Disruption (b=1.5pc, Adiabatic, GSPH)',
    path: 'simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/adiabatic_61k_gsph.json',
  },
  {
    id: 'Mc1e3_Mbh1e5_b1p5_v10_adiabatic_gdisph',
    name: 'Strong Disruption (b=1.5pc, Adiabatic, GDISPH)',
    path: 'simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/adiabatic_61k_gdisph.json',
  },
  {
    id: 'Mc1e3_Mbh1e5_b1p5_v10_radiative_gsph',
    name: 'Strong Disruption (b=1.5pc, Radiative, GSPH)',
    path: 'simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/radiative_61k_gsph.json',
  },
  {
    id: 'Mc1e3_Mbh1e5_b1p5_v10_radiative_gdisph',
    name: 'Strong Disruption (b=1.5pc, Radiative, GDISPH)',
    path: 'simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/radiative_61k_gdisph.json',
  },
  {
    id: 'Mc1e3_Mbh1e5_b3_v10_adiabatic',
    name: 'Moderate Disruption (b=3pc)',
    path: 'simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json',
  },
  {
    id: 'Mc1e3_Mbh1e5_b6_v10_adiabatic',
    name: 'Weak Disruption (b=6pc)',
    path: 'simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b6_v10/adiabatic_61k_gsph.json',
  },
] as const
