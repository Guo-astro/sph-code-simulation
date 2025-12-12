import { tanstackStart } from '@tanstack/react-start/plugin/vite'
import { defineConfig } from 'vite'
import tsConfigPaths from 'vite-tsconfig-paths'
import viteReact from '@vitejs/plugin-react'

export default defineConfig({
  server: {
    port: 3000,
  },
  plugins: [
    tsConfigPaths({
      projects: ['./tsconfig.json'],
    }),
    tanstackStart({
      srcDirectory: 'src',
    }),
    viteReact(),
  ],
  // Ensure Node.js modules are not bundled for client
  ssr: {
    noExternal: ['@tanstack/react-router', '@tanstack/react-start'],
  },
  optimizeDeps: {
    exclude: ['fs', 'path', 'url'],
  },
})
