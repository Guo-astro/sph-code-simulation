import { createFileRoute } from '@tanstack/react-router'
import { json } from '@tanstack/react-start'
import * as fs from 'fs'
import * as path from 'path'
import { fileURLToPath } from 'url'
import type { SimulationMetadata } from '~/types/sph'

/**
 * Simulation Metadata API
 * GET /api/simulations/[simId] - Returns metadata for a specific simulation
 */

/** ESM-compatible __dirname equivalent */
const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const getDataRoot = (): string => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname, '../../../../')
}

export const Route = createFileRoute('/api/simulations/$simId')({
  server: {
    handlers: {
      GET: async ({ params }) => {
        try {
          const { simId } = params
          const dataRoot = getDataRoot()
          
          // Decode simulation path - handle nested paths
          const simPath = decodeURIComponent(simId).replace(/\|/g, '/')
          
          // Try different possible locations
          // For nested structure like imbh_cloud/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph
          // The actual path is simulations/astrophysics/imbh_cloud/results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph/viz_data
          const pathParts = simPath.split('/')
          const possiblePaths = [
            path.join(dataRoot, 'simulations', 'astrophysics', simPath, 'viz_data'),
            path.join(dataRoot, 'simulations', 'astrophysics', simPath, 'results', 'viz_data'),
            path.join(dataRoot, 'lane_emden', 'results', simPath, 'viz_data'),
            path.join(dataRoot, simPath, 'viz_data'),
          ]
          
          // Handle nested structure: testName/scenario/method -> simulations/astrophysics/testName/results/scenario/method
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts
            possiblePaths.unshift(
              path.join(dataRoot, 'simulations', 'astrophysics', testName, 'results', ...rest, 'viz_data')
            )
          }

          let dataPath: string | null = null
          for (const p of possiblePaths) {
            if (fs.existsSync(p)) {
              dataPath = p
              break
            }
          }

          if (!dataPath) {
            return json({ error: `Simulation not found: ${simId}` }, { status: 404 })
          }

          // Load metadata
          const metadataPath = path.join(dataPath, 'metadata.json')
          if (fs.existsSync(metadataPath)) {
            const metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf-8'))
            return json({
              ...metadata,
              id: simId,
              dataPath,
            } as SimulationMetadata)
          }

          // Generate basic metadata from files
          const files = fs.readdirSync(dataPath)
          const frameFiles = files.filter((f) => f.startsWith('frame_') && f.endsWith('.bin'))
          
          const metadata: SimulationMetadata = {
            id: simId,
            name: simId.replace(/[_|]/g, ' '),
            description: `Simulation: ${simId}`,
            method: 'Unknown',
            kernel: 'Unknown',
            dimensions: 3,
            totalFrames: frameFiles.length,
            particleCount: 0,
            timeRange: [0, frameFiles.length],
            boundingBox: {
              min: [-1, -1, -1],
              max: [1, 1, 1],
            },
            dataPath,
            createdAt: new Date().toISOString(),
          }

          return json(metadata)
        } catch (error) {
          console.error('Error loading simulation metadata:', error)
          return json({ error: 'Failed to load simulation metadata' }, { status: 500 })
        }
      },
    },
  },
})
