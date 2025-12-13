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
    <div className={`bg-gray-800 p-2 rounded h-full flex flex-col ${className}`}>
      <h3 className="text-xs font-medium text-gray-300 mb-1 shrink-0">Energy Evolution</h3>
      <ResponsiveContainer width="100%" height="100%">
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
    <div className={`bg-gray-800 p-2 rounded h-full flex flex-col ${className}`}>
      <h3 className="text-xs font-medium text-gray-300 mb-1 shrink-0">Momentum Evolution</h3>
      <ResponsiveContainer width="100%" height="100%">
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

