import { createFileRoute } from '@tanstack/react-router'
import { json } from '@tanstack/react-start'
import * as fs from 'fs'
import * as path from 'path'
import { fileURLToPath } from 'url'
import type { SimulationMetadata, SimulationsListResponse, IMBHPhysicsConfig } from '~/types/sph'

/**
 * Simulations API - Lists available simulations
 * GET /api/simulations - Returns list of all available simulations
 */

/** ESM-compatible __dirname equivalent */
const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

/** Data directory configuration - can be overridden via environment */
const getDataRoot = (): string => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname, '../../../../../')
}

/**
 * Try to find and load IMBH physics from simulation config.
 * Searches for config files in common locations based on simulation path.
 */
function loadIMBHPhysicsFromConfig(
  dirPath: string,
  simulationId: string,
  dataRoot: string
): IMBHPhysicsConfig | undefined {
  // Try to find config file based on simulation path
  // e.g., imbh_cloud/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph
  const pathParts = simulationId.split('/')
  if (pathParts.length < 3) return undefined
  
  const [testCase, scenario, method] = pathParts
  
  // Common config locations to search
  const configPaths = [
    // Direct path: sample/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json
    path.join(dataRoot, 'sample', testCase, 'config', 'presets', 'simulation', 'scenarios', scenario, `${method}.json`),
    // Alternative: alongside results
    path.join(dirPath, '..', '..', 'config', `${method}.json`),
    // Root config
    path.join(dataRoot, 'config.json'),
  ]
  
  for (const configPath of configPaths) {
    if (fs.existsSync(configPath)) {
      try {
        const config = JSON.parse(fs.readFileSync(configPath, 'utf-8'))
        
        // Extract IMBH physics from various possible structures
        const imbhParams = config.imbh_parameters || config.externalForces?.pointMassBH
        const initialCondition = config.initialCondition?.transform
        const physicsSummary = config.physics_summary
        
        if (imbhParams || initialCondition) {
          // Extract BH mass (check multiple sources)
          const bhMassCode = imbhParams?.M_BH 
            ? imbhParams.M_BH / 1000  // Convert from M_sun to code units
            : config.externalForces?.pointMassBH?.mass || 100
          
          // Extract initial conditions
          const cloudPos = imbhParams?.cloud_initial_position 
            || initialCondition?.translate 
            || [-20, 3, 0]
          const cloudVel = imbhParams?.cloud_initial_velocity
            || initialCondition?.velocity_boost
            || [10, 0, 0]
          const bhPos = imbhParams?.BH_initial_position
            || config.externalForces?.pointMassBH?.position
            || [0, 0, 0]
          
          // Calculate impact parameter from initial y-position
          const impactParameter = Math.abs(cloudPos[1] - bhPos[1])
          
          // Tidal radius from physics summary or calculate
          const tidalRadius = physicsSummary?.tidal_radius_pc || 3.63
          
          // Cloud radius
          const cloudRadius = physicsSummary?.cloud_radius_pc || 1.13
          
          return {
            enabled: true,
            bhPosition: bhPos as [number, number, number],
            bhMass: bhMassCode,
            cloudInitialPosition: cloudPos as [number, number, number],
            cloudInitialVelocity: cloudVel as [number, number, number],
            cloudMass: (physicsSummary?.cloud_mass_Msun || 1000) / 1000,  // To code units
            cloudRadius,
            tidalRadius,
            impactParameter,
            timeUnit: 0.978,  // From unit system
          }
        }
      } catch (e) {
        console.warn(`Failed to parse config at ${configPath}:`, e)
      }
    }
  }
  
  return undefined
}

/** Scan for simulation data in a directory */
async function scanSimulationDirectory(
  dirPath: string,
  simulationId: string,
  dataRoot: string = getDataRoot()
): Promise<SimulationMetadata | null> {
  try {
    // Look for metadata.json first
    const metadataPath = path.join(dirPath, 'viz_data', 'metadata.json')
    if (fs.existsSync(metadataPath)) {
      const metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf-8'))
      
      // Try to load IMBH physics from config
      const imbhPhysics = loadIMBHPhysicsFromConfig(dirPath, simulationId, dataRoot)
      
      return {
        ...metadata,
        id: simulationId,
        dataPath: path.join(dirPath, 'viz_data'),
        imbhPhysics,
      }
    }

    // Fallback: scan for binary frame files
    const vizDataPath = path.join(dirPath, 'viz_data')
    if (!fs.existsSync(vizDataPath)) {
      return null
    }

    const files = fs.readdirSync(vizDataPath)
    const frameFiles = files.filter((f) => f.startsWith('frame_') && f.endsWith('.bin'))

    if (frameFiles.length === 0) {
      return null
    }

    // Extract basic metadata from file names
    const frameNumbers = frameFiles
      .map((f) => parseInt(f.replace('frame_', '').replace('.bin', '')))
      .sort((a, b) => a - b)

    return {
      id: simulationId,
      name: simulationId.replace(/_/g, ' '),
      description: `Simulation data from ${simulationId}`,
      method: 'Unknown',
      kernel: 'Unknown',
      dimensions: 3,
      totalFrames: frameFiles.length,
      particleCount: 0, // Will be determined when loading frame
      timeRange: [0, frameFiles.length],
      boundingBox: {
        min: [-1, -1, -1],
        max: [1, 1, 1],
      },
      dataPath: vizDataPath,
      createdAt: new Date().toISOString(),
    }
  } catch (error) {
    console.error(`Error scanning ${dirPath}:`, error)
    return null
  }
}

/** Find all simulations in the sample directory */
async function findSimulations(): Promise<SimulationMetadata[]> {
  const dataRoot = getDataRoot()
  const simulations: SimulationMetadata[] = []
  
  console.log('🔍 Starting simulation scan...')
  console.log(`   Data root: ${dataRoot}`)
  const startTime = Date.now()

  // Scan sample directory
  const sampleDir = path.join(dataRoot, 'sample')
  if (fs.existsSync(sampleDir)) {
    console.log('📁 Scanning sample/ directory...')
    const entries = fs.readdirSync(sampleDir, { withFileTypes: true })
    const totalEntries = entries.filter(e => e.isDirectory()).length
    let scanned = 0
    
    for (const entry of entries) {
      if (entry.isDirectory()) {
        scanned++
        process.stdout.write(`\r   [${scanned}/${totalEntries}] Scanning ${entry.name}...`.padEnd(60))
        
        // Check for results subdirectories (handles nested structure like results/scenario/method/)
        const resultsDir = path.join(sampleDir, entry.name, 'results')
        if (fs.existsSync(resultsDir)) {
          const resultEntries = fs.readdirSync(resultsDir, { withFileTypes: true })
          for (const resultEntry of resultEntries) {
            if (resultEntry.isDirectory()) {
              // First check if viz_data exists at this level
              const sim = await scanSimulationDirectory(
                path.join(resultsDir, resultEntry.name),
                `${entry.name}/${resultEntry.name}`,
                dataRoot
              )
              if (sim) {
                simulations.push(sim)
                console.log(`\n   ✓ Found: ${entry.name}/${resultEntry.name} (${sim.totalFrames} frames)`)
              }
              
              // Also check one level deeper (for structure like results/scenario/method/)
              const nestedDir = path.join(resultsDir, resultEntry.name)
              const nestedEntries = fs.readdirSync(nestedDir, { withFileTypes: true })
              for (const nestedEntry of nestedEntries) {
                if (nestedEntry.isDirectory()) {
                  const nestedSim = await scanSimulationDirectory(
                    path.join(nestedDir, nestedEntry.name),
                    `${entry.name}/${resultEntry.name}/${nestedEntry.name}`,
                    dataRoot
                  )
                  if (nestedSim) {
                    simulations.push(nestedSim)
                    console.log(`\n   ✓ Found: ${entry.name}/${resultEntry.name}/${nestedEntry.name} (${nestedSim.totalFrames} frames)`)
                  }
                }
              }
            }
          }
        }

        // Also check directly in the sample directory
        const sim = await scanSimulationDirectory(
          path.join(sampleDir, entry.name),
          entry.name,
          dataRoot
        )
        if (sim) {
          simulations.push(sim)
          console.log(`\n   ✓ Found: ${entry.name} (${sim.totalFrames} frames)`)
        }
      }
    }
    console.log('')  // New line after progress
  }

  // Scan lane_emden directory
  const laneEmdenDir = path.join(dataRoot, 'lane_emden', 'results')
  if (fs.existsSync(laneEmdenDir)) {
    console.log('📁 Scanning lane_emden/results/ directory...')
    const entries = fs.readdirSync(laneEmdenDir, { withFileTypes: true })
    for (const entry of entries) {
      if (entry.isDirectory()) {
        const sim = await scanSimulationDirectory(
          path.join(laneEmdenDir, entry.name),
          `lane_emden/${entry.name}`,
          dataRoot
        )
        if (sim) {
          simulations.push(sim)
          console.log(`   ✓ Found: lane_emden/${entry.name} (${sim.totalFrames} frames)`)
        }
      }
    }
  }

  const elapsed = ((Date.now() - startTime) / 1000).toFixed(2)
  console.log(`\n✅ Scan complete: Found ${simulations.length} simulations in ${elapsed}s`)
  
  return simulations
}

export const Route = createFileRoute('/api/simulations')({
  server: {
    handlers: {
      GET: async () => {
        try {
          const simulations = await findSimulations()
          const response: SimulationsListResponse = { simulations }
          return json(response)
        } catch (error) {
          console.error('Error listing simulations:', error)
          return json({ error: 'Failed to list simulations', simulations: [] }, { status: 500 })
        }
      },
    },
  },
})
