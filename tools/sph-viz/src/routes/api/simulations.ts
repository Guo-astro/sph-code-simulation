import { createFileRoute } from '@tanstack/react-router'
import { json } from '@tanstack/react-start'
import type { SimulationMetadata, SimulationsListResponse, IMBHPhysicsConfig } from '~/types/sph'

/**
 * Simulations API - Lists available simulations
 * GET /api/simulations - Returns list of all available simulations
 * 
 * Note: Node.js modules are loaded dynamically to avoid client-side bundling issues
 */

/** 
 * Server-only function to find simulations.
 * Uses dynamic imports to avoid Node.js modules being bundled for client.
 */
async function findSimulationsServer(): Promise<SimulationMetadata[]> {
  // Dynamic imports - only executed on server
  const fs = await import('fs')
  const path = await import('path')
  const { fileURLToPath } = await import('url')
  
  const __dirname = path.dirname(fileURLToPath(import.meta.url))
  const dataRoot = process.env.SPH_DATA_ROOT || path.resolve(__dirname, '../../../../../')
  
  console.log('🔍 Starting simulation scan...')
  console.log(`   Data root: ${dataRoot}`)
  
  const simulations: SimulationMetadata[] = []
  const startTime = Date.now()

  // Only scan IMBH cloud category directories
  const categoryDirs = ['CAT1', 'CAT2', 'CAT3', 'CAT_OKA']
  const imbhResultsDir = path.join(dataRoot, 'simulations', 'astrophysics', 'imbh_cloud', 'results')
  
  if (fs.existsSync(imbhResultsDir)) {
    console.log('📁 Scanning simulations/astrophysics/imbh_cloud/results/ categories...')
    
    for (const catName of categoryDirs) {
      const catDir = path.join(imbhResultsDir, catName)
      if (!fs.existsSync(catDir)) {
        console.log(`   ⏭️  Skipping ${catName} (not found)`)
        continue
      }
      
      console.log(`   📂 Scanning ${catName}...`)
      
      // Recursively scan category directory for simulations
      const scanCategoryDir = async (dir: string, prefix: string) => {
        const entries = fs.readdirSync(dir, { withFileTypes: true })
        for (const entry of entries) {
          if (entry.isDirectory()) {
            const entryPath = path.join(dir, entry.name)
            const simId = `imbh_cloud/${prefix}/${entry.name}`
            
            // Check if this directory has viz_data
            const vizDataPath = path.join(entryPath, 'viz_data')
            const metadataPath = path.join(vizDataPath, 'metadata.json')
            
            if (fs.existsSync(metadataPath)) {
              try {
                const metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf-8'))
                simulations.push({
                  ...metadata,
                  id: simId,
                  dataPath: vizDataPath,
                })
                console.log(`      ✓ Found: ${simId} (${metadata.totalFrames} frames)`)
              } catch (e) {
                console.error(`      ✗ Error reading ${metadataPath}:`, e)
              }
            } else {
              // Check one level deeper
              await scanCategoryDir(entryPath, `${prefix}/${entry.name}`)
            }
          }
        }
      }
      
      await scanCategoryDir(catDir, catName)
    }
  } else {
    console.log('⚠️  IMBH cloud results directory not found:', imbhResultsDir)
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
          const simulations = await findSimulationsServer()
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

