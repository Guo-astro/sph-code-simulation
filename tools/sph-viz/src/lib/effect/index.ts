/**
 * Effect Library - Main exports
 * 
 * This module provides Effect-based state management and type-safe schemas
 * for the SPH visualization application.
 * 
 * Usage:
 * 
 * ```tsx
 * import { EffectProvider, useVisualizationState } from '~/lib/effect'
 * 
 * function App() {
 *   return (
 *     <EffectProvider>
 *       <MyComponent />
 *     </EffectProvider>
 *   )
 * }
 * 
 * function MyComponent() {
 *   const state = useVisualizationState()
 *   // ...
 * }
 * ```
 */

// Schemas
export * from './schemas'

// Services
export * from './services'

// React Hooks
export * from './hooks'
