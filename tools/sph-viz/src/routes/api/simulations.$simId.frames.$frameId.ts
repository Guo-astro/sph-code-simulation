import { createFileRoute } from '@tanstack/react-router'
import { json } from '@tanstack/react-start'
import * as fs from 'fs'
import * as path from 'path'
import { fileURLToPath } from 'url'
import type { FrameResponse } from '~/types/sph'

/**
 * Frame Data API - Serves individual frame data
 * GET /api/simulations/[simId]/frames/[frameId] - Returns frame binary data
 */

/** ESM-compatible __dirname equivalent */
const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

/** Logging utility for consistent server-side logging */
const log = {
  info: (msg: string, ...args: any[]) => console.log(`[FRAME-API] ${new Date().toISOString()} INFO: ${msg}`, ...args),
  warn: (msg: string, ...args: any[]) => console.warn(`[FRAME-API] ${new Date().toISOString()} WARN: ${msg}`, ...args),
  error: (msg: string, ...args: any[]) => console.error(`[FRAME-API] ${new Date().toISOString()} ERROR: ${msg}`, ...args),
  debug: (msg: string, ...args: any[]) => console.log(`[FRAME-API] ${new Date().toISOString()} DEBUG: ${msg}`, ...args),
}

const getDataRoot = (): string => {
  const root = process.env.SPH_DATA_ROOT || path.resolve(__dirname, '../../../../../')
  log.debug(`getDataRoot() __dirname=${__dirname}, resolved=${root}`)
  return root
}

export const Route = createFileRoute('/api/simulations/$simId/frames/$frameId')({
  validateSearch: (search: Record<string, unknown>) => ({
    format: (search.format as string) || 'json',
  }),
  server: {
    handlers: {
      GET: async ({ params, request }) => {
        try {
          const { simId, frameId } = params
          const frameIndex = parseInt(frameId)
          
          // Check for binary format request via query param or Accept header
          const url = new URL(request.url)
          const formatParam = url.searchParams.get('format')
          const acceptHeader = request.headers.get('Accept') || ''
          const wantBinary = formatParam === 'binary' || acceptHeader.includes('application/octet-stream')

          console.log(`📥 Frame request: simId=${simId}, frameId=${frameId}, format=${wantBinary ? 'binary' : 'json'}`)

          if (isNaN(frameIndex)) {
            if (wantBinary) {
              return new Response('Invalid frame ID', { status: 400 })
            }
            return json({ error: 'Invalid frame ID' }, { status: 400 })
          }

          const dataRoot = getDataRoot()
          
          // Decode simulation path - handle nested paths
          const simPath = decodeURIComponent(simId).replace(/\|/g, '/')
          console.log(`   Decoded simPath: ${simPath}`)
          console.log(`   Data root: ${dataRoot}`)
          
          // Try different possible locations
          // For nested structure like imbh_cloud/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph
          // The actual path is sample/imbh_cloud/results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph/viz_data
          const pathParts = simPath.split('/')
          const possiblePaths = [
            path.join(dataRoot, 'sample', simPath, 'viz_data'),
            path.join(dataRoot, 'sample', simPath, 'results', 'viz_data'),
            path.join(dataRoot, 'lane_emden', 'results', simPath, 'viz_data'),
            path.join(dataRoot, simPath, 'viz_data'),
          ]
          
          // Handle nested structure: testName/scenario/method -> sample/testName/results/scenario/method
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts
            possiblePaths.unshift(
              path.join(dataRoot, 'sample', testName, 'results', ...rest, 'viz_data')
            )
          }
          
          console.log(`   Possible paths:`)
          possiblePaths.forEach((p, i) => console.log(`     ${i}: ${p}`))

          let dataPath: string | null = null
          for (const p of possiblePaths) {
            const exists = fs.existsSync(p)
            console.log(`   Checking: ${p} => ${exists}`)
            if (exists) {
              dataPath = p
              break
            }
          }

          if (!dataPath) {
            if (wantBinary) {
              return new Response(`Simulation data not found: ${simId}`, { status: 404 })
            }
            return json({ error: `Simulation data not found: ${simId}` }, { status: 404 })
          }

          // Load frame file
          const framePath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, '0')}.bin`)
          if (!fs.existsSync(framePath)) {
            if (wantBinary) {
              return new Response(`Frame ${frameIndex} not found`, { status: 404 })
            }
            return json({ error: `Frame ${frameIndex} not found` }, { status: 404 })
          }

          const buffer = fs.readFileSync(framePath)
          
          // Load metadata if available
          const metadataPath = path.join(dataPath, 'metadata.json')
          let metadata: any = {}
          if (fs.existsSync(metadataPath)) {
            metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf-8'))
          }

          // Load frame-specific info if available
          const frameInfoPath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, '0')}.json`)
          let frameInfo: any = { time: frameIndex }
          if (fs.existsSync(frameInfoPath)) {
            frameInfo = JSON.parse(fs.readFileSync(frameInfoPath, 'utf-8'))
          }

          // Default field offsets if not specified
          const fieldOffsets = metadata.fieldOffsets || {
            x: 0, y: 1, z: 2,
            vx: 3, vy: 4, vz: 5,
            mass: 6, density: 7, pressure: 8,
            energy: 9, smoothingLength: 10,
          }

          const stride = metadata.stride || 11
          const particleCount = Math.floor(buffer.byteLength / (stride * 4)) // Float32 = 4 bytes

          // Return binary response if requested
          if (wantBinary) {
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
          }

          // Otherwise return JSON with base64 data
          const response: FrameResponse = {
            frameIndex,
            time: frameInfo.time || frameIndex,
            data: buffer.toString('base64'),
            stride,
            fieldOffsets,
            particleCount,
          }

          return json(response, {
            headers: {
              'Cache-Control': 'public, max-age=3600',
            },
          })
        } catch (error) {
          console.error('Error loading frame:', error)
          return json({ error: 'Failed to load frame data' }, { status: 500 })
        }
      },
    },
  },
})
