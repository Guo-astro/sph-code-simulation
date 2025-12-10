import { createFileRoute, Link } from '@tanstack/react-router'

export const Route = createFileRoute('/')({
  component: Home,
})

function Home() {
  return (
    <div className="min-h-screen bg-gradient-to-br from-gray-900 via-gray-800 to-gray-900 flex items-center justify-center">
      <div className="text-center max-w-2xl px-8">
        <h1 className="text-5xl font-bold text-white mb-4">
          🌊 SPH Viz
        </h1>
        <p className="text-xl text-gray-400 mb-8">
          Interactive visualization dashboard for Smoothed Particle Hydrodynamics simulations
        </p>
        
        <div className="space-y-4">
          <Link
            to="/viz"
            className="inline-block px-8 py-4 bg-blue-600 hover:bg-blue-500 text-white text-lg font-semibold rounded-lg shadow-lg transition-all hover:scale-105"
          >
            Open Visualization Dashboard
          </Link>
          
          <div className="text-gray-500 text-sm mt-8">
            <p>Supports multiple SPH methods: GSPH, SSPH, DISPH, GDISPH, SRGSPH</p>
            <p className="mt-1">Real-time 3D rendering with React Three Fiber</p>
          </div>
        </div>

        <div className="mt-16 grid grid-cols-3 gap-8 text-center">
          <div className="p-4 bg-gray-800/50 rounded-lg">
            <div className="text-3xl mb-2">🎬</div>
            <div className="text-gray-400 text-sm">Animation Playback</div>
          </div>
          <div className="p-4 bg-gray-800/50 rounded-lg">
            <div className="text-3xl mb-2">📊</div>
            <div className="text-gray-400 text-sm">Real-time Statistics</div>
          </div>
          <div className="p-4 bg-gray-800/50 rounded-lg">
            <div className="text-3xl mb-2">🔬</div>
            <div className="text-gray-400 text-sm">Multiple Projections</div>
          </div>
        </div>
      </div>
    </div>
  )
}

