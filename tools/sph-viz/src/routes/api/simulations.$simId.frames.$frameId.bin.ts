import { createFileRoute } from '@tanstack/react-router'
import * as fs from 'fs'
import * as path from 'path'
import { fileURLToPath } from 'url'

/**
 * Binary Frame Data API - Serves raw binary frame data for maximum performance
 * GET /api/simulations/[simId]/frames/[frameId]/bin - Returns raw binary data
 * 
 * Metadata is returned in headers:
 * - X-Frame-Index: Frame index
 * - X-Frame-Time: Simulation time
 * - X-Particle-Count: Number of particles
 * - X-Stride: Fields per particle
 * - X-Field-Offsets: JSON-encoded field offsets
 */

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const getDataRoot = (): string => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname, '../../../../../')
}

export const Route = createFileRoute('/api/simulations/$simId/frames/$frameId/bin')({
  server: {
    handlers: {
      GET: async ({ params }) => {
        try {
          const { simId, frameId } = params
          const frameIndex = parseInt(frameId)

          if (isNaN(frameIndex)) {
            return new Response('Invalid frame ID', { status: 400 })
          }

          const dataRoot = getDataRoot()
          const simPath = decodeURIComponent(simId).replace(/\|/g, '/')
          
          // Try different possible locations
          const pathParts = simPath.split('/')
          const possiblePaths = [
            path.join(dataRoot, 'sample', simPath, 'viz_data'),
            path.join(dataRoot, 'sample', simPath, 'results', 'viz_data'),
            path.join(dataRoot, 'lane_emden', 'results', simPath, 'viz_data'),
            path.join(dataRoot, simPath, 'viz_data'),
          ]
          
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts
            possiblePaths.unshift(
              path.join(dataRoot, 'sample', testName, 'results', ...rest, 'viz_data')
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
            return new Response('Simulation not found', { status: 404 })
          }

          // Load frame file
          const framePath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, '0')}.bin`)
          if (!fs.existsSync(framePath)) {
            return new Response('Frame not found', { status: 404 })
          }

          const buffer = fs.readFileSync(framePath)
          
          // Load metadata
          const metadataPath = path.join(dataPath, 'metadata.json')
          let metadata: any = {}
          if (fs.existsSync(metadataPath)) {
            metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf-8'))
          }

          // Load frame info
          const frameInfoPath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, '0')}.json`)
          let frameInfo: any = { time: frameIndex }
          if (fs.existsSync(frameInfoPath)) {
            frameInfo = JSON.parse(fs.readFileSync(frameInfoPath, 'utf-8'))
          }

          const fieldOffsets = metadata.fieldOffsets || {
            x: 0, y: 1, z: 2,
            vx: 3, vy: 4, vz: 5,
            mass: 6, density: 7, pressure: 8,
            energy: 9, smoothing_length: 10,
          }

          const stride = metadata.stride || 11
          const particleCount = Math.floor(buffer.byteLength / (stride * 4))

          // Return raw binary with metadata in headers
          return new Response(buffer, {
            status: 200,
            headers: {
              'Content-Type': 'application/octet-stream',
              'Content-Length': buffer.byteLength.toString(),
              'Cache-Control': 'public, max-age=3600',
              'X-Frame-Index': frameIndex.toString(),
              'X-Frame-Time': (frameInfo.time || frameIndex).toString(),
              'X-Particle-Count': particleCount.toString(),
              'X-Stride': stride.toString(),
              'X-Field-Offsets': JSON.stringify(fieldOffsets),
            },
          })
        } catch (error) {
          console.error('Error loading binary frame:', error)
          return new Response('Server error', { status: 500 })
        }
      },
    },
  },
})
