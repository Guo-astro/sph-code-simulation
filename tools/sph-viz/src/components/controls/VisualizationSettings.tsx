'use client'

import { useState } from 'react'
import { Settings, ChevronDown, ChevronRight } from 'lucide-react'
import { COLOR_MAPS, type ColorMap } from '~/types/sph'

export interface VisualizationSettingsProps {
  colorField: string
  onColorFieldChange: (field: string) => void
  colorMapName: string
  onColorMapChange: (name: string) => void
  pointSize: number
  onPointSizeChange: (size: number) => void
  opacity: number
  onOpacityChange: (opacity: number) => void
  showAxes: boolean
  onShowAxesChange: (show: boolean) => void
  showBoundingBox: boolean
  onShowBoundingBoxChange: (show: boolean) => void
  colorRange: [number, number]
  onColorRangeChange: (range: [number, number]) => void
  useLogScale: boolean
  onLogScaleChange: (log: boolean) => void
  availableFields?: string[]
  className?: string
}

const DEFAULT_FIELDS = ['density', 'pressure', 'energy', 'velocity', 'machNumber']

/** Collapsible settings panel */
export function VisualizationSettings({
  colorField,
  onColorFieldChange,
  colorMapName,
  onColorMapChange,
  pointSize,
  onPointSizeChange,
  opacity,
  onOpacityChange,
  showAxes,
  onShowAxesChange,
  showBoundingBox,
  onShowBoundingBoxChange,
  colorRange,
  onColorRangeChange,
  useLogScale,
  onLogScaleChange,
  availableFields = DEFAULT_FIELDS,
  className = '',
}: VisualizationSettingsProps) {
  const [expanded, setExpanded] = useState(true)

  return (
    <div className={`bg-gray-800 rounded ${className}`}>
      <button
        onClick={() => setExpanded(!expanded)}
        className="w-full flex items-center gap-2 p-3 text-left hover:bg-gray-700/50"
      >
        <Settings size={16} className="text-gray-400" />
        <span className="text-sm font-medium text-gray-300 flex-1">Visualization Settings</span>
        {expanded ? (
          <ChevronDown size={16} className="text-gray-400" />
        ) : (
          <ChevronRight size={16} className="text-gray-400" />
        )}
      </button>

      {expanded && (
        <div className="px-3 pb-3 space-y-3">
          {/* Color Field */}
          <div>
            <label className="text-xs text-gray-400 block mb-1">Color Field</label>
            <select
              value={colorField}
              onChange={(e) => onColorFieldChange(e.target.value)}
              className="w-full bg-gray-700 text-white text-sm rounded px-2 py-1.5"
            >
              {availableFields.map((field) => (
                <option key={field} value={field}>
                  {fieldDisplayName(field)}
                </option>
              ))}
            </select>
          </div>

          {/* Color Map */}
          <div>
            <label className="text-xs text-gray-400 block mb-1">Color Map</label>
            <select
              value={colorMapName}
              onChange={(e) => onColorMapChange(e.target.value)}
              className="w-full bg-gray-700 text-white text-sm rounded px-2 py-1.5"
            >
              {Object.keys(COLOR_MAPS).map((name) => (
                <option key={name} value={name}>
                  {COLOR_MAPS[name].name}
                </option>
              ))}
            </select>
            {/* Color map preview */}
            <div
              className="h-3 mt-1 rounded"
              style={{
                background: `linear-gradient(to right, ${COLOR_MAPS[colorMapName]?.colors.join(', ')})`,
              }}
            />
          </div>

          {/* Color Range */}
          <div>
            <label className="text-xs text-gray-400 block mb-1">Color Range</label>
            <div className="flex gap-2 items-center">
              <input
                type="number"
                value={colorRange[0]}
                onChange={(e) => onColorRangeChange([parseFloat(e.target.value), colorRange[1]])}
                className="w-20 bg-gray-700 text-white text-xs rounded px-2 py-1"
                step="any"
              />
              <span className="text-gray-500">—</span>
              <input
                type="number"
                value={colorRange[1]}
                onChange={(e) => onColorRangeChange([colorRange[0], parseFloat(e.target.value)])}
                className="w-20 bg-gray-700 text-white text-xs rounded px-2 py-1"
                step="any"
              />
              <button
                onClick={() => onColorRangeChange([0, 0])} // Reset to auto
                className="text-xs text-blue-400 hover:text-blue-300"
              >
                Auto
              </button>
            </div>
          </div>

          {/* Log Scale */}
          <div className="flex items-center gap-2">
            <input
              type="checkbox"
              id="logScale"
              checked={useLogScale}
              onChange={(e) => onLogScaleChange(e.target.checked)}
              className="rounded bg-gray-700"
            />
            <label htmlFor="logScale" className="text-xs text-gray-400">
              Logarithmic Scale
            </label>
          </div>

          {/* Point Size */}
          <div>
            <label className="text-xs text-gray-400 block mb-1">
              Point Size: {pointSize.toFixed(3)}
            </label>
            <input
              type="range"
              min={0.001}
              max={0.1}
              step={0.001}
              value={pointSize}
              onChange={(e) => onPointSizeChange(parseFloat(e.target.value))}
              className="w-full h-1.5 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
            />
          </div>

          {/* Opacity */}
          <div>
            <label className="text-xs text-gray-400 block mb-1">
              Opacity: {(opacity * 100).toFixed(0)}%
            </label>
            <input
              type="range"
              min={0.1}
              max={1}
              step={0.05}
              value={opacity}
              onChange={(e) => onOpacityChange(parseFloat(e.target.value))}
              className="w-full h-1.5 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
            />
          </div>

          {/* Display Options */}
          <div className="space-y-2">
            <label className="text-xs text-gray-400 block">Display Options</label>
            <div className="flex items-center gap-2">
              <input
                type="checkbox"
                id="showAxes"
                checked={showAxes}
                onChange={(e) => onShowAxesChange(e.target.checked)}
                className="rounded bg-gray-700"
              />
              <label htmlFor="showAxes" className="text-xs text-gray-400">
                Show Axes
              </label>
            </div>
            <div className="flex items-center gap-2">
              <input
                type="checkbox"
                id="showBoundingBox"
                checked={showBoundingBox}
                onChange={(e) => onShowBoundingBoxChange(e.target.checked)}
                className="rounded bg-gray-700"
              />
              <label htmlFor="showBoundingBox" className="text-xs text-gray-400">
                Show Bounding Box
              </label>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

function fieldDisplayName(field: string): string {
  const names: Record<string, string> = {
    density: 'Density (ρ)',
    pressure: 'Pressure (P)',
    energy: 'Internal Energy (u)',
    velocity: 'Velocity (|v|)',
    machNumber: 'Mach Number',
    smoothingLength: 'Smoothing Length (h)',
    mass: 'Mass (m)',
    divV: 'Velocity Divergence (∇·v)',
    temperature: 'Temperature (T)',
  }
  return names[field] || field
}

export default VisualizationSettings
