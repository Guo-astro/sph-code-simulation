import { createFileRoute } from '@tanstack/react-router'
import { json } from '@tanstack/react-start'
import type { SimulationsListResponse } from '~/types/sph'
import { findSimulationsServer } from '~/lib/server/findSimulations'

/**
 * Simulations API - Lists available simulations
 * GET /api/simulations - Returns list of all available simulations
 */

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

