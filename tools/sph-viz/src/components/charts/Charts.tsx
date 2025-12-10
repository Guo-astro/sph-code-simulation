'use client'

import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceLine,
} from 'recharts'
import type { FrameStatistics } from '~/types/sph'

export interface EnergyChartProps {
  statistics: FrameStatistics[]
  currentFrame?: number
  showKinetic?: boolean
  showInternal?: boolean
  showTotal?: boolean
  className?: string
}

/** Energy evolution chart */
export function EnergyChart({
  statistics,
  currentFrame,
  showKinetic = true,
  showInternal = true,
  showTotal = true,
  className = '',
}: EnergyChartProps) {
  // Prepare data for recharts - filter out undefined entries
  const data = statistics
    .filter((s): s is FrameStatistics => s !== undefined && s !== null)
    .map((s) => ({
      time: s.time,
      kinetic: s.totalKineticEnergy,
      internal: s.totalInternalEnergy,
      total: s.totalEnergy,
    }))

  // Find current time for reference line
  const currentTime = currentFrame !== undefined && statistics[currentFrame]
    ? statistics[currentFrame].time
    : undefined

  return (
    <div className={`bg-gray-800 p-3 rounded ${className}`}>
      <h3 className="text-sm font-medium text-gray-300 mb-2">Energy Evolution</h3>
      <ResponsiveContainer width="100%" height={200}>
        <LineChart data={data}>
          <CartesianGrid strokeDasharray="3 3" stroke="#444" />
          <XAxis
            dataKey="time"
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => typeof v === 'number' && isFinite(v) ? v.toFixed(2) : '0'}
          />
          <YAxis
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => typeof v === 'number' && isFinite(v) ? v.toExponential(1) : '0'}
          />
          <Tooltip
            contentStyle={{ backgroundColor: '#1f2937', border: 'none' }}
            labelStyle={{ color: '#9ca3af' }}
            formatter={(value: number) => typeof value === 'number' && isFinite(value) ? value.toExponential(3) : '0'}
          />
          <Legend wrapperStyle={{ fontSize: 10 }} />
          
          {showKinetic && (
            <Line
              type="monotone"
              dataKey="kinetic"
              name="Kinetic"
              stroke="#3b82f6"
              strokeWidth={1.5}
              dot={false}
            />
          )}
          {showInternal && (
            <Line
              type="monotone"
              dataKey="internal"
              name="Internal"
              stroke="#ef4444"
              strokeWidth={1.5}
              dot={false}
            />
          )}
          {showTotal && (
            <Line
              type="monotone"
              dataKey="total"
              name="Total"
              stroke="#22c55e"
              strokeWidth={2}
              dot={false}
            />
          )}
          
          {currentTime !== undefined && (
            <ReferenceLine x={currentTime} stroke="#fff" strokeDasharray="3 3" />
          )}
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}

export interface MomentumChartProps {
  statistics: FrameStatistics[]
  currentFrame?: number
  className?: string
}

/** Momentum evolution chart */
export function MomentumChart({
  statistics,
  currentFrame,
  className = '',
}: MomentumChartProps) {
  // Filter out undefined entries
  const data = statistics
    .filter((s): s is FrameStatistics => s !== undefined && s !== null)
    .map((s) => ({
      time: s.time,
      px: s.momentum[0],
      py: s.momentum[1],
      pz: s.momentum[2],
      total: Math.sqrt(s.momentum[0] ** 2 + s.momentum[1] ** 2 + s.momentum[2] ** 2),
    }))

  const currentTime = currentFrame !== undefined && statistics[currentFrame]
    ? statistics[currentFrame].time
    : undefined

  return (
    <div className={`bg-gray-800 p-3 rounded ${className}`}>
      <h3 className="text-sm font-medium text-gray-300 mb-2">Momentum Evolution</h3>
      <ResponsiveContainer width="100%" height={200}>
        <LineChart data={data}>
          <CartesianGrid strokeDasharray="3 3" stroke="#444" />
          <XAxis
            dataKey="time"
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => typeof v === 'number' && isFinite(v) ? v.toFixed(2) : '0'}
          />
          <YAxis
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => typeof v === 'number' && isFinite(v) ? v.toExponential(1) : '0'}
          />
          <Tooltip
            contentStyle={{ backgroundColor: '#1f2937', border: 'none' }}
            labelStyle={{ color: '#9ca3af' }}
            formatter={(value: number) => typeof value === 'number' && isFinite(value) ? value.toExponential(3) : '0'}
          />
          <Legend wrapperStyle={{ fontSize: 10 }} />
          
          <Line type="monotone" dataKey="px" name="Px" stroke="#3b82f6" strokeWidth={1} dot={false} />
          <Line type="monotone" dataKey="py" name="Py" stroke="#ef4444" strokeWidth={1} dot={false} />
          <Line type="monotone" dataKey="pz" name="Pz" stroke="#22c55e" strokeWidth={1} dot={false} />
          <Line type="monotone" dataKey="total" name="|P|" stroke="#f59e0b" strokeWidth={2} dot={false} />
          
          {currentTime !== undefined && (
            <ReferenceLine x={currentTime} stroke="#fff" strokeDasharray="3 3" />
          )}
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}

export interface RadialProfileChartProps {
  radius: number[]
  density?: number[]
  pressure?: number[]
  velocity?: number[]
  analyticalDensity?: number[]
  analyticalPressure?: number[]
  analyticalVelocity?: number[]
  showDensity?: boolean
  showPressure?: boolean
  showVelocity?: boolean
  logScale?: boolean
  className?: string
}

/** Radial profile chart */
export function RadialProfileChart({
  radius,
  density,
  pressure,
  velocity,
  analyticalDensity,
  analyticalPressure,
  analyticalVelocity,
  showDensity = true,
  showPressure = true,
  showVelocity = true,
  logScale = false,
  className = '',
}: RadialProfileChartProps) {
  // Prepare data
  const data = radius.map((r, i) => ({
    r,
    density: density?.[i],
    pressure: pressure?.[i],
    velocity: velocity?.[i],
    analyticalDensity: analyticalDensity?.[i],
    analyticalPressure: analyticalPressure?.[i],
    analyticalVelocity: analyticalVelocity?.[i],
  }))

  return (
    <div className={`bg-gray-800 p-3 rounded ${className}`}>
      <h3 className="text-sm font-medium text-gray-300 mb-2">Radial Profile</h3>
      <ResponsiveContainer width="100%" height={250}>
        <LineChart data={data}>
          <CartesianGrid strokeDasharray="3 3" stroke="#444" />
          <XAxis
            dataKey="r"
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => typeof v === 'number' && isFinite(v) ? v.toFixed(2) : '0'}
            label={{ value: 'Radius', position: 'bottom', fill: '#888', fontSize: 10 }}
          />
          <YAxis
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => typeof v === 'number' && isFinite(v) ? (logScale ? v.toExponential(0) : v.toFixed(2)) : '0'}
            scale={logScale ? 'log' : 'auto'}
            domain={logScale ? ['auto', 'auto'] : undefined}
          />
          <Tooltip
            contentStyle={{ backgroundColor: '#1f2937', border: 'none' }}
            labelStyle={{ color: '#9ca3af' }}
            formatter={(value: number) => typeof value === 'number' && isFinite(value) ? value.toExponential(3) : '0'}
          />
          <Legend wrapperStyle={{ fontSize: 10 }} />
          
          {showDensity && density && (
            <Line
              type="monotone"
              dataKey="density"
              name="ρ (sim)"
              stroke="#3b82f6"
              strokeWidth={2}
              dot={false}
            />
          )}
          {showDensity && analyticalDensity && (
            <Line
              type="monotone"
              dataKey="analyticalDensity"
              name="ρ (analytical)"
              stroke="#3b82f6"
              strokeWidth={1}
              strokeDasharray="5 5"
              dot={false}
            />
          )}
          
          {showPressure && pressure && (
            <Line
              type="monotone"
              dataKey="pressure"
              name="P (sim)"
              stroke="#ef4444"
              strokeWidth={2}
              dot={false}
            />
          )}
          {showPressure && analyticalPressure && (
            <Line
              type="monotone"
              dataKey="analyticalPressure"
              name="P (analytical)"
              stroke="#ef4444"
              strokeWidth={1}
              strokeDasharray="5 5"
              dot={false}
            />
          )}
          
          {showVelocity && velocity && (
            <Line
              type="monotone"
              dataKey="velocity"
              name="v (sim)"
              stroke="#22c55e"
              strokeWidth={2}
              dot={false}
            />
          )}
          {showVelocity && analyticalVelocity && (
            <Line
              type="monotone"
              dataKey="analyticalVelocity"
              name="v (analytical)"
              stroke="#22c55e"
              strokeWidth={1}
              strokeDasharray="5 5"
              dot={false}
            />
          )}
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}

export interface HistogramChartProps {
  values: number[]
  bins?: number
  label?: string
  color?: string
  logY?: boolean
  className?: string
}

/** Histogram chart for distribution visualization */
export function HistogramChart({
  values,
  bins = 50,
  label = 'Value',
  color = '#3b82f6',
  logY = false,
  className = '',
}: HistogramChartProps) {
  // Compute histogram
  const filteredValues = values.filter((v) => isFinite(v))
  if (filteredValues.length === 0) {
    return (
      <div className={`bg-gray-800 p-3 rounded ${className}`}>
        <h3 className="text-sm font-medium text-gray-300 mb-2">{label} Distribution</h3>
        <div className="text-gray-500 text-center py-8">No data</div>
      </div>
    )
  }

  const min = Math.min(...filteredValues)
  const max = Math.max(...filteredValues)
  const binWidth = (max - min) / bins

  const histogram = new Array(bins).fill(0)
  for (const v of filteredValues) {
    const binIndex = Math.min(Math.floor((v - min) / binWidth), bins - 1)
    histogram[binIndex]++
  }

  const data = histogram.map((count, i) => ({
    value: min + (i + 0.5) * binWidth,
    count,
  }))

  return (
    <div className={`bg-gray-800 p-3 rounded ${className}`}>
      <h3 className="text-sm font-medium text-gray-300 mb-2">{label} Distribution</h3>
      <ResponsiveContainer width="100%" height={150}>
        <LineChart data={data}>
          <CartesianGrid strokeDasharray="3 3" stroke="#444" />
          <XAxis
            dataKey="value"
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            tickFormatter={(v) => v.toExponential(1)}
          />
          <YAxis
            stroke="#888"
            tick={{ fill: '#888', fontSize: 10 }}
            scale={logY ? 'log' : 'auto'}
            domain={logY ? [1, 'auto'] : undefined}
          />
          <Tooltip
            contentStyle={{ backgroundColor: '#1f2937', border: 'none' }}
            labelStyle={{ color: '#9ca3af' }}
          />
          <Line
            type="stepAfter"
            dataKey="count"
            stroke={color}
            strokeWidth={2}
            fill={color}
            fillOpacity={0.3}
            dot={false}
          />
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}
