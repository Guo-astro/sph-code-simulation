/// <reference types="vite/client" />
import {
  HeadContent,
  Link,
  Scripts,
  createRootRoute,
} from '@tanstack/react-router'
import * as React from 'react'
import { DefaultCatchBoundary } from '~/components/DefaultCatchBoundary'
import { NotFound } from '~/components/NotFound'
import appCss from '~/styles/app.css?url'
import { seo } from '~/utils/seo'

// Client-only devtools component
function DevtoolsWrapper() {
  const [DevtoolsComponent, setDevtoolsComponent] = React.useState<React.ComponentType<any> | null>(null)

  React.useEffect(() => {
    // Only load devtools on client in development
    if (process.env.NODE_ENV !== 'production') {
      import('@tanstack/react-router-devtools').then((mod) => {
        setDevtoolsComponent(() => mod.TanStackRouterDevtools as React.ComponentType<any>)
      })
    }
  }, [])

  if (!DevtoolsComponent) return null
  return <DevtoolsComponent position="bottom-right" />
}

export const Route = createRootRoute({
  head: () => ({
    meta: [
      {
        charSet: 'utf-8',
      },
      {
        name: 'viewport',
        content: 'width=device-width, initial-scale=1',
      },
      ...seo({
        title: 'SPH Visualization Tool',
        description: 'Interactive visualization dashboard for SPH simulation data',
      }),
    ],
    links: [
      { rel: 'stylesheet', href: appCss },
      { rel: 'icon', href: '/favicon.ico' },
    ],
  }),
  errorComponent: DefaultCatchBoundary,
  notFoundComponent: () => <NotFound />,
  shellComponent: RootDocument,
})

function RootDocument({ children }: { children: React.ReactNode }) {
  return (
    <html className="dark">
      <head>
        <HeadContent />
      </head>
      <body className="bg-gray-900 text-white">
        {children}
        <DevtoolsWrapper />
        <Scripts />
      </body>
    </html>
  )
}
