/**
 * Server-side function to find simulations.
 * Uses dynamic imports to avoid Node.js modules being bundled for client.
 */

import type { SimulationMetadata } from '~/types/sph'

export async function findSimulationsServer(): Promise<SimulationMetadata[]> {
  // Dynamic imports - only executed on server
  const fs = await import('fs')
  const path = await import('path')
  const { fileURLToPath } = await import('url')

  const __dirname = path.dirname(fileURLToPath(import.meta.url))
  const dataRoot = process.env.SPH_DATA_ROOT || path.resolve(__dirname, '../../../../../')

  const simulations: SimulationMetadata[] = []
  const startTime = Date.now()

  // Only scan IMBH cloud category directories
  const categoryDirs = ['CAT1', 'CAT2', 'CAT3', 'CAT_OKA']
  const imbhResultsDir = path.join(dataRoot, 'simulations', 'astrophysics', 'imbh_cloud', 'results')

  if (fs.existsSync(imbhResultsDir)) {
    for (const catName of categoryDirs) {
      const catDir = path.join(imbhResultsDir, catName)
      if (!fs.existsSync(catDir)) {
        continue
      }

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
              } catch (e) {
                // Silently skip invalid metadata files
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
  }

  return simulations
}
