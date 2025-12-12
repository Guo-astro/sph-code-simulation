import { jsxs, jsx, Fragment } from "react/jsx-runtime";
import { useRouter, useMatch, rootRouteId, ErrorComponent, Link, createRootRoute, HeadContent, Scripts, createFileRoute, lazyRouteComponent, redirect, notFound, createRouter } from "@tanstack/react-router";
import { TanStackRouterDevtools } from "@tanstack/react-router-devtools";
import { T as TSS_SERVER_FUNCTION, g as getServerFnById, c as createServerFn, j as json, a as getRequestHeaders } from "../server.js";
import { useState, useEffect, useMemo, useRef, useCallback } from "react";
import { Canvas, useThree } from "@react-three/fiber";
import { PerspectiveCamera, OrbitControls, Stats } from "@react-three/drei";
import * as THREE from "three";
import { OrbitControls as OrbitControls$1 } from "three/examples/jsm/controls/OrbitControls.js";
import { ResponsiveContainer, LineChart, CartesianGrid, XAxis, YAxis, Tooltip, Legend, Line, ReferenceLine } from "recharts";
import { SkipBack, ChevronLeft, Pause, Play, ChevronRight, SkipForward, Settings, ChevronDown, Folder, RefreshCw } from "lucide-react";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const createMiddleware = (options, __opts) => {
  const resolvedOptions = {
    type: "request",
    ...__opts || options
  };
  return {
    options: resolvedOptions,
    middleware: (middleware) => {
      return createMiddleware(
        {},
        Object.assign(resolvedOptions, { middleware })
      );
    },
    inputValidator: (inputValidator) => {
      return createMiddleware(
        {},
        Object.assign(resolvedOptions, { inputValidator })
      );
    },
    client: (client) => {
      return createMiddleware(
        {},
        Object.assign(resolvedOptions, { client })
      );
    },
    server: (server) => {
      return createMiddleware(
        {},
        Object.assign(resolvedOptions, { server })
      );
    }
  };
};
function DefaultCatchBoundary({ error }) {
  const router2 = useRouter();
  const isRoot = useMatch({
    strict: false,
    select: (state) => state.id === rootRouteId
  });
  console.error("DefaultCatchBoundary Error:", error);
  return /* @__PURE__ */ jsxs("div", { className: "min-w-0 flex-1 p-4 flex flex-col items-center justify-center gap-6", children: [
    /* @__PURE__ */ jsx(ErrorComponent, { error }),
    /* @__PURE__ */ jsxs("div", { className: "flex gap-2 items-center flex-wrap", children: [
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: () => {
            router2.invalidate();
          },
          className: `px-2 py-1 bg-gray-600 dark:bg-gray-700 rounded-sm text-white uppercase font-extrabold`,
          children: "Try Again"
        }
      ),
      isRoot ? /* @__PURE__ */ jsx(
        Link,
        {
          to: "/",
          className: `px-2 py-1 bg-gray-600 dark:bg-gray-700 rounded-sm text-white uppercase font-extrabold`,
          children: "Home"
        }
      ) : /* @__PURE__ */ jsx(
        Link,
        {
          to: "/",
          className: `px-2 py-1 bg-gray-600 dark:bg-gray-700 rounded-sm text-white uppercase font-extrabold`,
          onClick: (e) => {
            e.preventDefault();
            window.history.back();
          },
          children: "Go Back"
        }
      )
    ] })
  ] });
}
function NotFound({ children }) {
  return /* @__PURE__ */ jsxs("div", { className: "space-y-2 p-2", children: [
    /* @__PURE__ */ jsx("div", { className: "text-gray-600 dark:text-gray-400", children: children || /* @__PURE__ */ jsx("p", { children: "The page you are looking for does not exist." }) }),
    /* @__PURE__ */ jsxs("p", { className: "flex items-center gap-2 flex-wrap", children: [
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: () => window.history.back(),
          className: "bg-emerald-500 text-white px-2 py-1 rounded-sm uppercase font-black text-sm",
          children: "Go back"
        }
      ),
      /* @__PURE__ */ jsx(
        Link,
        {
          to: "/",
          className: "bg-cyan-600 text-white px-2 py-1 rounded-sm uppercase font-black text-sm",
          children: "Start Over"
        }
      )
    ] })
  ] });
}
const appCss = "/assets/app-uvzNok8D.css";
const seo = ({
  title,
  description,
  keywords,
  image
}) => {
  const tags = [
    { title },
    { name: "description", content: description },
    { name: "keywords", content: keywords },
    { name: "twitter:title", content: title },
    { name: "twitter:description", content: description },
    { name: "twitter:creator", content: "@tannerlinsley" },
    { name: "twitter:site", content: "@tannerlinsley" },
    { name: "og:type", content: "website" },
    { name: "og:title", content: title },
    { name: "og:description", content: description },
    ...image ? [
      { name: "twitter:image", content: image },
      { name: "twitter:card", content: "summary_large_image" },
      { name: "og:image", content: image }
    ] : []
  ];
  return tags;
};
const Route$m = createRootRoute({
  head: () => ({
    meta: [
      {
        charSet: "utf-8"
      },
      {
        name: "viewport",
        content: "width=device-width, initial-scale=1"
      },
      ...seo({
        title: "SPH Visualization Tool",
        description: "Interactive visualization dashboard for SPH simulation data"
      })
    ],
    links: [
      { rel: "stylesheet", href: appCss },
      { rel: "icon", href: "/favicon.ico" }
    ]
  }),
  errorComponent: DefaultCatchBoundary,
  notFoundComponent: () => /* @__PURE__ */ jsx(NotFound, {}),
  shellComponent: RootDocument
});
function RootDocument({ children }) {
  return /* @__PURE__ */ jsxs("html", { className: "dark", children: [
    /* @__PURE__ */ jsx("head", { children: /* @__PURE__ */ jsx(HeadContent, {}) }),
    /* @__PURE__ */ jsxs("body", { className: "bg-gray-900 text-white", children: [
      children,
      /* @__PURE__ */ jsx(TanStackRouterDevtools, { position: "bottom-right" }),
      /* @__PURE__ */ jsx(Scripts, {})
    ] })
  ] });
}
const $$splitComponentImporter$c = () => import("./users-BsoLiE3B.js");
const Route$l = createFileRoute("/users")({
  loader: async () => {
    const res = await fetch("/api/users");
    if (!res.ok) {
      throw new Error("Unexpected status code");
    }
    const data = await res.json();
    return data;
  },
  component: lazyRouteComponent($$splitComponentImporter$c, "component")
});
const Route$k = createFileRoute("/redirect")({
  beforeLoad: () => {
    throw redirect({
      to: "/posts"
    });
  }
});
const createSsrRpc = (functionId) => {
  const url = "/_serverFn/" + functionId;
  const fn = async (...args) => {
    const serverFn = await getServerFnById(functionId);
    return serverFn(...args);
  };
  return Object.assign(fn, {
    url,
    functionId,
    [TSS_SERVER_FUNCTION]: true
  });
};
const fetchPost_createServerFn_handler = createSsrRpc("0029094260fc8f554fa3ac223696de0e9591567ec6420250e896c91244c812c5");
const fetchPost = createServerFn({
  method: "POST"
}).inputValidator((d) => d).handler(fetchPost_createServerFn_handler, async ({
  data,
  context
}) => {
  console.log("Request context:", context);
  console.info(`Fetching post with id ${data}...`);
  const res = await fetch(`https://jsonplaceholder.typicode.com/posts/${data}`);
  if (!res.ok) {
    if (res.status === 404) {
      throw notFound();
    }
    throw new Error("Failed to fetch post");
  }
  const post = await res.json();
  return post;
});
const fetchPosts_createServerFn_handler = createSsrRpc("cbb8ca69048418e62742f2c511faa56326b80ace384144a35bb3e0bf5e8124be");
const fetchPosts = createServerFn().handler(fetchPosts_createServerFn_handler, async () => {
  console.info("Fetching posts...");
  const res = await fetch("https://jsonplaceholder.typicode.com/posts");
  if (!res.ok) {
    throw new Error("Failed to fetch posts");
  }
  const posts = await res.json();
  return posts.slice(0, 10);
});
const $$splitComponentImporter$b = () => import("./posts-fzlzT6gP.js");
const Route$j = createFileRoute("/posts")({
  loader: async () => fetchPosts(),
  component: lazyRouteComponent($$splitComponentImporter$b, "component")
});
const $$splitComponentImporter$a = () => import("./deferred-BUEJ58HT.js");
const personServerFn_createServerFn_handler = createSsrRpc("f76e8f8721c12c8547a3ced6a10916f5b5076c1a10dcbeaa607360ce419d0a48");
const personServerFn = createServerFn({
  method: "GET"
}).inputValidator((d) => d).handler(personServerFn_createServerFn_handler, ({
  data: name
}) => {
  return {
    name,
    randomNumber: Math.floor(Math.random() * 100)
  };
});
const slowServerFn_createServerFn_handler = createSsrRpc("fc3988c64f434639dfd4eab3f926b87ee39cc0c14f65b4d0e852c7fd73279a3b");
const slowServerFn = createServerFn({
  method: "GET"
}).inputValidator((d) => d).handler(slowServerFn_createServerFn_handler, async ({
  data: name
}) => {
  await new Promise((r) => setTimeout(r, 1e3));
  return {
    name,
    randomNumber: Math.floor(Math.random() * 100)
  };
});
const Route$i = createFileRoute("/deferred")({
  loader: async () => {
    return {
      deferredStuff: new Promise((r) => setTimeout(() => r("Hello deferred!"), 2e3)),
      deferredPerson: slowServerFn({
        data: "Tanner Linsley"
      }),
      person: await personServerFn({
        data: "John Doe"
      })
    };
  },
  component: lazyRouteComponent($$splitComponentImporter$a, "component")
});
const Route$h = createFileRoute("/customScript.js")({
  server: {
    handlers: {
      GET: () => {
        return new Response('console.log("Hello from customScript.js!")', {
          headers: {
            "Content-Type": "application/javascript"
          }
        });
      }
    }
  }
});
const $$splitComponentImporter$9 = () => import("./_pathlessLayout-BhrcpZGS.js");
const Route$g = createFileRoute("/_pathlessLayout")({
  component: lazyRouteComponent($$splitComponentImporter$9, "component")
});
const $$splitComponentImporter$8 = () => import("./index-DrJwRs12.js");
const Route$f = createFileRoute("/")({
  component: lazyRouteComponent($$splitComponentImporter$8, "component")
});
function createCircleTexture$1() {
  const size = 64;
  const canvas = document.createElement("canvas");
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext("2d");
  const gradient = ctx.createRadialGradient(
    size / 2,
    size / 2,
    0,
    size / 2,
    size / 2,
    size / 2
  );
  gradient.addColorStop(0, "rgba(255, 255, 255, 1)");
  gradient.addColorStop(0.3, "rgba(255, 255, 255, 0.8)");
  gradient.addColorStop(0.7, "rgba(255, 255, 255, 0.3)");
  gradient.addColorStop(1, "rgba(255, 255, 255, 0)");
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  const texture = new THREE.CanvasTexture(canvas);
  texture.needsUpdate = true;
  return texture;
}
let circleTexture$1 = null;
function getCircleTexture$1() {
  if (!circleTexture$1) {
    circleTexture$1 = createCircleTexture$1();
  }
  return circleTexture$1;
}
function ParticleCloud({ frame, colorField, colorMap, pointSize, opacity }) {
  const pointsRef = useRef(null);
  const geometryRef = useRef(null);
  const materialRef = useRef(null);
  const colorsArrayRef = useRef(null);
  useEffect(() => {
    if (!frame || !pointsRef.current) return;
    const count = frame.particleCount;
    let geometry = geometryRef.current;
    if (!geometry || geometry.attributes.position && geometry.attributes.position.count !== count) {
      geometry = new THREE.BufferGeometry();
      const positions = new THREE.BufferAttribute(frame.positions, 3);
      positions.setUsage(THREE.DynamicDrawUsage);
      geometry.setAttribute("position", positions);
      colorsArrayRef.current = new Float32Array(count * 3);
      const colors2 = new THREE.BufferAttribute(colorsArrayRef.current, 3);
      colors2.setUsage(THREE.DynamicDrawUsage);
      geometry.setAttribute("color", colors2);
      geometryRef.current = geometry;
      pointsRef.current.geometry = geometry;
    } else {
      const posAttr = geometry.attributes.position;
      posAttr.array = frame.positions;
      posAttr.needsUpdate = true;
    }
    let fieldData;
    switch (colorField) {
      case "density":
        fieldData = frame.density;
        break;
      case "pressure":
        fieldData = frame.pressure;
        break;
      case "energy":
        fieldData = frame.energy;
        break;
      case "velocity":
        if (!colorsArrayRef.current || colorsArrayRef.current.length !== count * 3) {
          colorsArrayRef.current = new Float32Array(count * 3);
        }
        const velMag = new Float32Array(count);
        for (let i = 0; i < count; i++) {
          const vx = frame.velocities[i * 3];
          const vy = frame.velocities[i * 3 + 1];
          const vz = frame.velocities[i * 3 + 2];
          velMag[i] = Math.sqrt(vx * vx + vy * vy + vz * vz);
        }
        fieldData = velMag;
        break;
      case "machNumber":
        fieldData = frame.machNumber;
        break;
      default:
        fieldData = frame.density;
    }
    let min = colorMap.min ?? Infinity;
    let max = colorMap.max ?? -Infinity;
    if (min === Infinity || max === -Infinity) {
      if (fieldData) {
        for (let i = 0; i < count; i++) {
          const val = fieldData[i];
          if (isFinite(val)) {
            if (val < min) min = val;
            if (val > max) max = val;
          }
        }
      }
    }
    if (min === max) max = min + 1;
    const colorAttr = geometry.attributes.color;
    const colors = colorAttr.array;
    const logMin = colorMap.logScale && min > 0 ? Math.log10(min) : 0;
    const logRange = colorMap.logScale && min > 0 ? Math.log10(max) - logMin : 1;
    const range = max - min;
    for (let i = 0; i < count; i++) {
      let val = fieldData ? fieldData[i] : 0;
      if (!isFinite(val)) val = min;
      let t;
      if (colorMap.logScale && min > 0) {
        t = (Math.log10(val) - logMin) / logRange;
      } else {
        t = (val - min) / range;
      }
      t = Math.max(0, Math.min(1, t));
      const color = interpolateColorFast(colorMap.colors, t);
      colors[i * 3] = color.r;
      colors[i * 3 + 1] = color.g;
      colors[i * 3 + 2] = color.b;
    }
    colorAttr.needsUpdate = true;
    geometry.computeBoundingSphere();
  }, [frame, colorField, colorMap]);
  useEffect(() => {
    if (materialRef.current) {
      materialRef.current.size = pointSize;
      materialRef.current.opacity = opacity;
      materialRef.current.needsUpdate = true;
    }
  }, [pointSize, opacity]);
  if (!frame) return null;
  return /* @__PURE__ */ jsxs("points", { ref: pointsRef, children: [
    /* @__PURE__ */ jsx("bufferGeometry", { ref: geometryRef }),
    /* @__PURE__ */ jsx(
      "pointsMaterial",
      {
        ref: materialRef,
        size: pointSize,
        vertexColors: true,
        transparent: true,
        opacity,
        sizeAttenuation: true,
        depthWrite: false,
        blending: THREE.AdditiveBlending,
        map: getCircleTexture$1(),
        alphaTest: 0.01
      }
    )
  ] });
}
const colorCache = /* @__PURE__ */ new Map();
function interpolateColorFast(colors, t) {
  if (colors.length === 0) return { r: 1, g: 1, b: 1 };
  if (colors.length === 1) return hexToRgbCached(colors[0]);
  const index = t * (colors.length - 1);
  const lower = Math.floor(index);
  const upper = Math.min(lower + 1, colors.length - 1);
  const localT = index - lower;
  const c1 = hexToRgbCached(colors[lower]);
  const c2 = hexToRgbCached(colors[upper]);
  return {
    r: c1.r + (c2.r - c1.r) * localT,
    g: c1.g + (c2.g - c1.g) * localT,
    b: c1.b + (c2.b - c1.b) * localT
  };
}
function hexToRgbCached(hex) {
  let cached = colorCache.get(hex);
  if (!cached) {
    cached = hexToRgb$1(hex);
    colorCache.set(hex, cached);
  }
  return cached;
}
function hexToRgb$1(hex) {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!result) return { r: 1, g: 1, b: 1 };
  return {
    r: parseInt(result[1], 16) / 255,
    g: parseInt(result[2], 16) / 255,
    b: parseInt(result[3], 16) / 255
  };
}
function AxesHelper({ size = 1 }) {
  return /* @__PURE__ */ jsx("axesHelper", { args: [size] });
}
function GridHelper({ size = 10, divisions = 10 }) {
  return /* @__PURE__ */ jsx("gridHelper", { args: [size, divisions, "#444444", "#222222"], rotation: [Math.PI / 2, 0, 0] });
}
function BoundingBox({
  min,
  max
}) {
  const geometry = useMemo(() => {
    const box = new THREE.BoxGeometry(max[0] - min[0], max[1] - min[1], max[2] - min[2]);
    box.translate((max[0] + min[0]) / 2, (max[1] + min[1]) / 2, (max[2] + min[2]) / 2);
    return new THREE.EdgesGeometry(box);
  }, [min, max]);
  return /* @__PURE__ */ jsx("lineSegments", { geometry, children: /* @__PURE__ */ jsx("lineBasicMaterial", { color: "#666666" }) });
}
function CameraController({
  boundingBox,
  resetKey
}) {
  const { camera } = useThree();
  useEffect(() => {
    if (boundingBox) {
      const center = new THREE.Vector3(
        (boundingBox.min[0] + boundingBox.max[0]) / 2,
        (boundingBox.min[1] + boundingBox.max[1]) / 2,
        (boundingBox.min[2] + boundingBox.max[2]) / 2
      );
      const size = new THREE.Vector3(
        boundingBox.max[0] - boundingBox.min[0],
        boundingBox.max[1] - boundingBox.min[1],
        boundingBox.max[2] - boundingBox.min[2]
      );
      const maxDim = Math.max(size.x, size.y, size.z);
      const distance = maxDim * 2;
      camera.position.set(center.x + distance, center.y + distance * 0.5, center.z + distance);
      camera.lookAt(center);
    }
  }, [boundingBox, camera, resetKey]);
  return null;
}
const defaultColorMap$1 = {
  name: "Viridis",
  colors: [
    "#440154",
    "#482878",
    "#3e4989",
    "#31688e",
    "#26828e",
    "#1f9e89",
    "#35b779",
    "#6ece58",
    "#b5de2b",
    "#fde725"
  ],
  min: 0,
  max: 1,
  logScale: false
};
function ParticleViewer3D({
  frame,
  colorField = "density",
  colorMap = defaultColorMap$1,
  pointSize = 0.02,
  opacity = 0.8,
  showAxes = true,
  showGrid = false,
  showBoundingBox = true,
  boundingBox,
  showStats = false,
  className = ""
}) {
  const [resetKey, setResetKey] = useState(0);
  const handleResetCamera = () => {
    setResetKey((k) => k + 1);
  };
  return /* @__PURE__ */ jsxs("div", { className: `relative w-full h-full bg-gray-900 ${className}`, children: [
    /* @__PURE__ */ jsxs(Canvas, { children: [
      /* @__PURE__ */ jsx(PerspectiveCamera, { makeDefault: true, fov: 60, near: 1e-3, far: 1e3 }),
      /* @__PURE__ */ jsx(OrbitControls, { enableDamping: true, dampingFactor: 0.1, rotateSpeed: 0.5 }),
      /* @__PURE__ */ jsx(CameraController, { boundingBox, resetKey }),
      /* @__PURE__ */ jsx("ambientLight", { intensity: 0.5 }),
      /* @__PURE__ */ jsx("directionalLight", { position: [10, 10, 5], intensity: 1 }),
      showAxes && /* @__PURE__ */ jsx(AxesHelper, { size: boundingBox ? Math.max(...boundingBox.max) * 0.5 : 1 }),
      showGrid && /* @__PURE__ */ jsx(GridHelper, { size: 10, divisions: 10 }),
      showBoundingBox && boundingBox && /* @__PURE__ */ jsx(BoundingBox, { min: boundingBox.min, max: boundingBox.max }),
      /* @__PURE__ */ jsx(
        ParticleCloud,
        {
          frame,
          colorField,
          colorMap,
          pointSize,
          opacity
        }
      ),
      showStats && /* @__PURE__ */ jsx(Stats, {})
    ] }),
    /* @__PURE__ */ jsx("div", { className: "absolute top-2 right-2 flex gap-2", children: /* @__PURE__ */ jsx(
      "button",
      {
        onClick: handleResetCamera,
        className: "px-3 py-1 bg-gray-700 text-white text-sm rounded hover:bg-gray-600",
        children: "Reset Camera"
      }
    ) }),
    frame && /* @__PURE__ */ jsxs("div", { className: "absolute bottom-2 left-2 text-white text-xs bg-black/50 px-2 py-1 rounded", children: [
      "Frame ",
      frame.frameIndex,
      " | t = ",
      frame.time.toFixed(4),
      " | ",
      frame.particleCount.toLocaleString(),
      " particles"
    ] })
  ] });
}
const UNITS = {
  // M_sun (1 code mass = 1000 M_sun)
  // G = 1 in code units: [G] = [L]³/([M][T]²) where T = L/V
  // Verified: GM = 100 (code) gives v_circ ≈ 10 km/s at r = 1 pc
  G: 1
};
function computeHyperbolicOrbit(config, numPoints = 200) {
  const { bhPosition, bhMass, cloudInitialPosition, cloudInitialVelocity } = config;
  const [x0, y0, z0] = cloudInitialPosition;
  const [vx0, vy0, vz0] = cloudInitialVelocity;
  const [bh_x, bh_y, bh_z] = bhPosition;
  const v_inf = Math.sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0);
  const b = Math.abs(y0 - bh_y);
  const GM = UNITS.G * bhMass;
  const a = GM / (v_inf * v_inf);
  const e = Math.sqrt(1 + Math.pow(b * v_inf * v_inf / GM, 2));
  const p = a * (e * e - 1);
  const theta_asymp = Math.PI - Math.acos(1 / e);
  const sign_y = y0 >= bh_y ? 1 : -1;
  const nu_max = Math.PI - Math.acos(1 / e) - 0.02;
  const nu_min = -nu_max;
  const orbit = [];
  for (let i = 0; i < numPoints; i++) {
    const nu = nu_min + (nu_max - nu_min) * (i / (numPoints - 1));
    const r = p / (1 + e * Math.cos(nu));
    const rotation = Math.PI - theta_asymp * sign_y;
    const theta = nu + rotation;
    const x = bh_x + r * Math.cos(theta);
    const y = bh_y + r * Math.sin(theta) * sign_y;
    const z = bh_z + z0 * (1 - Math.abs(nu) / nu_max * 0.1);
    if (Math.abs(x) < 100 && Math.abs(y) < 100) {
      orbit.push([x, y, z]);
    }
  }
  return orbit;
}
function computeStraightLineTrajectory(config, t_start, t_end, numPoints = 50) {
  const { cloudInitialPosition, cloudInitialVelocity } = config;
  const [x0, y0, z0] = cloudInitialPosition;
  const [vx, vy, vz] = cloudInitialVelocity;
  const trajectory = [];
  for (let i = 0; i < numPoints; i++) {
    const t = t_start + (t_end - t_start) * (i / (numPoints - 1));
    trajectory.push([
      x0 + vx * t,
      y0 + vy * t,
      z0 + vz * t
    ]);
  }
  return trajectory;
}
function createCircleTexture() {
  const size = 64;
  const canvas = document.createElement("canvas");
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext("2d");
  const gradient = ctx.createRadialGradient(
    size / 2,
    size / 2,
    0,
    size / 2,
    size / 2,
    size / 2
  );
  gradient.addColorStop(0, "rgba(255, 255, 255, 1)");
  gradient.addColorStop(0.3, "rgba(255, 255, 255, 0.8)");
  gradient.addColorStop(0.7, "rgba(255, 255, 255, 0.3)");
  gradient.addColorStop(1, "rgba(255, 255, 255, 0)");
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  const texture = new THREE.CanvasTexture(canvas);
  texture.needsUpdate = true;
  return texture;
}
let circleTexture = null;
function getCircleTexture() {
  if (!circleTexture) {
    circleTexture = createCircleTexture();
  }
  return circleTexture;
}
const COLOR_MAP_DATA = {
  viridis: [[0.267, 4e-3, 0.329], [0.282, 0.14, 0.458], [0.253, 0.265, 0.53], [0.206, 0.372, 0.553], [0.163, 0.471, 0.558], [0.127, 0.566, 0.551], [0.134, 0.658, 0.518], [0.266, 0.749, 0.441], [0.477, 0.821, 0.318], [0.741, 0.873, 0.15], [0.993, 0.906, 0.144]],
  plasma: [[0.05, 0.03, 0.528], [0.295, 0.012, 0.615], [0.492, 0.012, 0.658], [0.665, 0.138, 0.618], [0.798, 0.28, 0.47], [0.899, 0.396, 0.301], [0.966, 0.53, 0.128], [0.988, 0.68, 0.063], [0.961, 0.85, 0.298], [0.94, 0.975, 0.131]],
  inferno: [[1e-3, 0, 0.014], [0.122, 0.047, 0.282], [0.304, 0.063, 0.42], [0.499, 0.086, 0.397], [0.68, 0.144, 0.295], [0.833, 0.253, 0.16], [0.937, 0.405, 0.049], [0.981, 0.588, 0.068], [0.987, 0.772, 0.264], [0.988, 0.998, 0.645]],
  turbo: [[0.19, 0.072, 0.232], [0.235, 0.318, 0.86], [0.137, 0.572, 0.996], [0.14, 0.78, 0.82], [0.376, 0.92, 0.512], [0.67, 0.979, 0.28], [0.924, 0.904, 0.145], [0.996, 0.724, 0.132], [0.994, 0.472, 0.122], [0.881, 0.2, 0.102], [0.528, 0.055, 0.052]]
};
function sampleColorMap(name, t) {
  const map = COLOR_MAP_DATA[name] || COLOR_MAP_DATA.viridis;
  t = Math.max(0, Math.min(1, t));
  const idx = t * (map.length - 1);
  const i = Math.floor(idx);
  const f = idx - i;
  if (i >= map.length - 1) return map[map.length - 1];
  return [
    map[i][0] + f * (map[i + 1][0] - map[i][0]),
    map[i][1] + f * (map[i + 1][1] - map[i][1]),
    map[i][2] + f * (map[i + 1][2] - map[i][2])
  ];
}
function ParticleViewer3DImperative({
  framesRef,
  initialFrameIndex = 0,
  frameIndexRef,
  colorField = "density",
  colorMapName = "viridis",
  pointSize = 2,
  opacity = 0.9,
  logScale = false,
  showAxes = true,
  showBoundingBox = true,
  boundingBox,
  className = "",
  onFpsUpdate,
  globalColorRange,
  imbhPhysics,
  showBlackHole = true,
  showTrajectory = true,
  showRadii = true
}) {
  const containerRef = useRef(null);
  const rendererRef = useRef(null);
  const sceneRef = useRef(null);
  const cameraRef = useRef(null);
  const controlsRef = useRef(null);
  const particlesRef = useRef(null);
  const geometryRef = useRef(null);
  const materialRef = useRef(null);
  const bhGroupRef = useRef(null);
  const trajectoryLineRef = useRef(null);
  const trajectoryPointsRef = useRef([]);
  const comMarkerRef = useRef(null);
  const tidalCircleRef = useRef(null);
  const hillCircleRef = useRef(null);
  const lastFrameIndexRef = useRef(-1);
  const lastColorFieldRef = useRef(colorField);
  const lastColorMapRef = useRef(colorMapName);
  const lastLogScaleRef = useRef(logScale);
  const globalColorRangeRef = useRef(globalColorRange);
  useEffect(() => {
    globalColorRangeRef.current = globalColorRange;
  }, [globalColorRange]);
  const statsRef = useRef({ density: [0, 1], velocity: [0, 1], pressure: [0, 1], energy: [0, 1] });
  const fpsRef = useRef({ frames: 0, lastTime: performance.now(), fps: 0 });
  const computeStats = useCallback(() => {
    const frames = framesRef.current;
    if (!frames || frames.size === 0) return;
    const allDensity = [];
    const allVelocity = [];
    const allPressure = [];
    const allEnergy = [];
    frames.forEach((frame, idx) => {
      if (idx % 10 !== 0) return;
      for (let i = 0; i < frame.particleCount; i += 100) {
        allDensity.push(frame.density[i]);
        const vx = frame.velocities[i * 3];
        const vy = frame.velocities[i * 3 + 1];
        const vz = frame.velocities[i * 3 + 2];
        allVelocity.push(Math.sqrt(vx * vx + vy * vy + vz * vz));
        allPressure.push(frame.pressure[i]);
        allEnergy.push(frame.energy[i]);
      }
    });
    const percentile = (arr, p) => {
      const sorted = arr.slice().sort((a, b) => a - b);
      return sorted[Math.floor(sorted.length * p / 100)] || 0;
    };
    statsRef.current = {
      density: [percentile(allDensity, 1), percentile(allDensity, 99)],
      velocity: [percentile(allVelocity, 1), percentile(allVelocity, 99)],
      pressure: [percentile(allPressure, 1), percentile(allPressure, 99)],
      energy: [percentile(allEnergy, 1), percentile(allEnergy, 99)]
    };
  }, [framesRef]);
  const updateParticles = useCallback(() => {
    const frameIndex = frameIndexRef.current ?? 0;
    const frames = framesRef.current;
    if (!frames || !geometryRef.current) return;
    const frame = frames.get(frameIndex);
    if (!frame) return;
    const needsUpdate = frameIndex !== lastFrameIndexRef.current || colorField !== lastColorFieldRef.current || colorMapName !== lastColorMapRef.current || logScale !== lastLogScaleRef.current;
    if (!needsUpdate) return;
    lastFrameIndexRef.current = frameIndex;
    lastColorFieldRef.current = colorField;
    lastColorMapRef.current = colorMapName;
    lastLogScaleRef.current = logScale;
    const geometry = geometryRef.current;
    const posAttr = geometry.attributes.position;
    geometry.attributes.color;
    if (posAttr.count !== frame.particleCount) {
      const positions2 = new Float32Array(frame.particleCount * 3);
      const colors2 = new Float32Array(frame.particleCount * 3);
      geometry.setAttribute("position", new THREE.BufferAttribute(positions2, 3).setUsage(THREE.DynamicDrawUsage));
      geometry.setAttribute("color", new THREE.BufferAttribute(colors2, 3).setUsage(THREE.DynamicDrawUsage));
    }
    const positions = geometry.attributes.position.array;
    positions.set(frame.positions);
    geometry.attributes.position.needsUpdate = true;
    let fieldData;
    let vMin, vMax;
    const useGlobalRange = globalColorRangeRef.current && globalColorRangeRef.current[0] !== globalColorRangeRef.current[1];
    switch (colorField) {
      case "velocity": {
        fieldData = new Float32Array(frame.particleCount);
        for (let i = 0; i < frame.particleCount; i++) {
          const vx = frame.velocities[i * 3];
          const vy = frame.velocities[i * 3 + 1];
          const vz = frame.velocities[i * 3 + 2];
          fieldData[i] = Math.sqrt(vx * vx + vy * vy + vz * vz);
        }
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current;
        } else {
          [vMin, vMax] = statsRef.current.velocity;
        }
        break;
      }
      case "pressure":
        fieldData = frame.pressure;
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current;
        } else {
          [vMin, vMax] = statsRef.current.pressure;
        }
        break;
      case "energy":
        fieldData = frame.energy;
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current;
        } else {
          [vMin, vMax] = statsRef.current.energy;
        }
        break;
      default:
        fieldData = frame.density;
        if (useGlobalRange) {
          [vMin, vMax] = globalColorRangeRef.current;
        } else {
          [vMin, vMax] = statsRef.current.density;
        }
    }
    const colors = geometry.attributes.color.array;
    const logMin = logScale && vMin > 0 ? Math.log10(vMin) : 0;
    const logRange = logScale && vMin > 0 ? Math.log10(vMax) - logMin : 1;
    const range = vMax - vMin || 1;
    for (let i = 0; i < frame.particleCount; i++) {
      let val = fieldData[i];
      if (!isFinite(val)) val = vMin;
      let t;
      if (logScale && vMin > 0) {
        t = (Math.log10(Math.max(val, vMin)) - logMin) / logRange;
      } else {
        t = (val - vMin) / range;
      }
      t = Math.max(0, Math.min(1, t));
      const [r, g, b] = sampleColorMap(colorMapName, t);
      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;
    }
    geometry.attributes.color.needsUpdate = true;
    geometry.computeBoundingSphere();
  }, [framesRef, frameIndexRef, colorField, colorMapName, logScale]);
  const updateTrajectory = useCallback(() => {
    if (!imbhPhysics?.enabled) return;
    const frames = framesRef?.current;
    const frameIndex = frameIndexRef?.current ?? 0;
    if (!frames || frames.size === 0) return;
    const frame = frames.get(frameIndex);
    if (!frame) return;
    if (frameIndex === 0 && trajectoryPointsRef.current.length > 1) {
      trajectoryPointsRef.current = [];
      if (trajectoryLineRef.current) {
        trajectoryLineRef.current.geometry.setDrawRange(0, 0);
      }
    }
    let totalMass = 0;
    let comX = 0, comY = 0, comZ = 0;
    for (let i = 0; i < frame.particleCount; i++) {
      const mass = frame.mass?.[i] ?? 1;
      totalMass += mass;
      comX += mass * frame.positions[i * 3];
      comY += mass * frame.positions[i * 3 + 1];
      comZ += mass * frame.positions[i * 3 + 2];
    }
    if (totalMass > 0) {
      comX /= totalMass;
      comY /= totalMass;
      comZ /= totalMass;
    }
    if (comMarkerRef.current) {
      comMarkerRef.current.position.set(comX, comY, comZ);
      comMarkerRef.current.visible = showTrajectory;
    }
    const lastIdx = trajectoryPointsRef.current.length > 0 ? trajectoryPointsRef.current[trajectoryPointsRef.current.length - 1] : null;
    const threshold = 0.01;
    if (!lastIdx || Math.abs(comX - lastIdx[0]) > threshold || Math.abs(comY - lastIdx[1]) > threshold || Math.abs(comZ - lastIdx[2]) > threshold) {
      trajectoryPointsRef.current.push([comX, comY, comZ]);
      const maxPoints = 500;
      if (trajectoryPointsRef.current.length > maxPoints) {
        trajectoryPointsRef.current = trajectoryPointsRef.current.slice(-maxPoints);
      }
      if (trajectoryLineRef.current && trajectoryPointsRef.current.length >= 2) {
        const positions = new Float32Array(trajectoryPointsRef.current.length * 3);
        for (let i = 0; i < trajectoryPointsRef.current.length; i++) {
          positions[i * 3] = trajectoryPointsRef.current[i][0];
          positions[i * 3 + 1] = trajectoryPointsRef.current[i][1];
          positions[i * 3 + 2] = trajectoryPointsRef.current[i][2];
        }
        trajectoryLineRef.current.geometry.setAttribute(
          "position",
          new THREE.BufferAttribute(positions, 3)
        );
        trajectoryLineRef.current.geometry.setDrawRange(0, trajectoryPointsRef.current.length);
        trajectoryLineRef.current.geometry.attributes.position.needsUpdate = true;
      }
    }
    if (hillCircleRef.current) {
      hillCircleRef.current.position.set(comX, comY, 0);
    }
  }, [framesRef, frameIndexRef, imbhPhysics, showTrajectory]);
  useEffect(() => {
    if (!containerRef.current) return;
    const container = containerRef.current;
    const width = container.clientWidth;
    const height = container.clientHeight;
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(657935);
    sceneRef.current = scene;
    const camera = new THREE.PerspectiveCamera(60, width / height, 0.1, 1e4);
    camera.position.set(50, 50, 100);
    cameraRef.current = camera;
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;
    const controls = new OrbitControls$1(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.screenSpacePanning = true;
    controls.minDistance = 1;
    controls.maxDistance = 5e3;
    controlsRef.current = controls;
    const gridHelper = new THREE.GridHelper(200, 20, 4473924, 2236962);
    scene.add(gridHelper);
    const createTextSprite = (text, color = "#ffffff") => {
      const canvas = document.createElement("canvas");
      const size = 256;
      canvas.width = size;
      canvas.height = size / 2;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = color;
      ctx.font = "Bold 48px Arial";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, size / 2, size / 4);
      const texture = new THREE.CanvasTexture(canvas);
      const material2 = new THREE.SpriteMaterial({ map: texture, transparent: true });
      const sprite = new THREE.Sprite(material2);
      sprite.scale.set(8, 4, 1);
      return sprite;
    };
    if (showAxes) {
      const axesHelper = new THREE.AxesHelper(30);
      scene.add(axesHelper);
      const xLabel = createTextSprite("X (pc)", "#ff4444");
      xLabel.position.set(32, 0, 0);
      scene.add(xLabel);
      const yLabel = createTextSprite("Y (pc)", "#44ff44");
      yLabel.position.set(0, 32, 0);
      scene.add(yLabel);
      const zLabel = createTextSprite("Z (pc)", "#4444ff");
      zLabel.position.set(0, 0, 32);
      scene.add(zLabel);
      const tickMaterial = new THREE.LineBasicMaterial({ color: 6710886 });
      for (let i = -20; i <= 20; i += 10) {
        if (i === 0) continue;
        const xTickGeom = new THREE.BufferGeometry().setFromPoints([
          new THREE.Vector3(i, -0.5, 0),
          new THREE.Vector3(i, 0.5, 0)
        ]);
        scene.add(new THREE.Line(xTickGeom, tickMaterial));
        const yTickGeom = new THREE.BufferGeometry().setFromPoints([
          new THREE.Vector3(-0.5, i, 0),
          new THREE.Vector3(0.5, i, 0)
        ]);
        scene.add(new THREE.Line(yTickGeom, tickMaterial));
        if (i % 10 === 0) {
          const tickLabel = createTextSprite(`${i}`, "#888888");
          tickLabel.position.set(i, -2, 0);
          tickLabel.scale.set(4, 2, 1);
          scene.add(tickLabel);
        }
      }
    }
    const geometry = new THREE.BufferGeometry();
    const initialCount = 64e3;
    const positions = new Float32Array(initialCount * 3);
    const colors = new Float32Array(initialCount * 3);
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage));
    geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3).setUsage(THREE.DynamicDrawUsage));
    geometryRef.current = geometry;
    const material = new THREE.PointsMaterial({
      size: pointSize,
      vertexColors: true,
      sizeAttenuation: true,
      transparent: true,
      opacity,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      map: getCircleTexture(),
      alphaTest: 0.01
    });
    materialRef.current = material;
    const particles = new THREE.Points(geometry, material);
    scene.add(particles);
    particlesRef.current = particles;
    if (imbhPhysics?.enabled) {
      const bhPos = imbhPhysics.bhPosition;
      if (showBlackHole) {
        const bhGroup = new THREE.Group();
        bhGroup.position.set(bhPos[0], bhPos[1], bhPos[2]);
        const glowColors = [16729156, 16737894, 16746632, 16755370];
        const glowSizes = [0.8, 0.6, 0.4, 0.25];
        glowSizes.forEach((size, i) => {
          const glowGeo = new THREE.SphereGeometry(size, 32, 32);
          const glowMat = new THREE.MeshBasicMaterial({
            color: glowColors[i],
            transparent: true,
            opacity: 0.15 - i * 0.03
          });
          const glowMesh = new THREE.Mesh(glowGeo, glowMat);
          bhGroup.add(glowMesh);
        });
        const coreGeo = new THREE.SphereGeometry(0.15, 32, 32);
        const coreMat = new THREE.MeshBasicMaterial({ color: 16777215 });
        const coreMesh = new THREE.Mesh(coreGeo, coreMat);
        bhGroup.add(coreMesh);
        const horizonGeo = new THREE.SphereGeometry(0.08, 32, 32);
        const horizonMat = new THREE.MeshBasicMaterial({ color: 0 });
        const horizonMesh = new THREE.Mesh(horizonGeo, horizonMat);
        bhGroup.add(horizonMesh);
        scene.add(bhGroup);
        bhGroupRef.current = bhGroup;
      }
      if (showRadii && imbhPhysics.tidalRadius > 0) {
        const tidalRadius = imbhPhysics.tidalRadius;
        const circlePoints = [];
        for (let i = 0; i <= 64; i++) {
          const theta = i / 64 * Math.PI * 2;
          circlePoints.push(new THREE.Vector3(
            bhPos[0] + tidalRadius * Math.cos(theta),
            bhPos[1] + tidalRadius * Math.sin(theta),
            bhPos[2]
          ));
        }
        const tidalGeo = new THREE.BufferGeometry().setFromPoints(circlePoints);
        const tidalMat = new THREE.LineDashedMaterial({
          color: 65535,
          dashSize: 0.5,
          gapSize: 0.25,
          transparent: true,
          opacity: 0.6
        });
        const tidalLine = new THREE.Line(tidalGeo, tidalMat);
        tidalLine.computeLineDistances();
        scene.add(tidalLine);
        tidalCircleRef.current = tidalLine;
        const hillRadius = tidalRadius * 0.4;
        const hillPoints = [];
        for (let i = 0; i <= 64; i++) {
          const theta = i / 64 * Math.PI * 2;
          hillPoints.push(new THREE.Vector3(
            bhPos[0] + hillRadius * Math.cos(theta),
            bhPos[1] + hillRadius * Math.sin(theta),
            bhPos[2]
          ));
        }
        const hillGeo = new THREE.BufferGeometry().setFromPoints(hillPoints);
        const hillMat = new THREE.LineDashedMaterial({
          color: 16746496,
          dashSize: 0.3,
          gapSize: 0.15,
          transparent: true,
          opacity: 0.5
        });
        const hillLine = new THREE.Line(hillGeo, hillMat);
        hillLine.computeLineDistances();
        scene.add(hillLine);
        hillCircleRef.current = hillLine;
      }
      if (showTrajectory) {
        const orbitPoints = computeHyperbolicOrbit({
          bhPosition: imbhPhysics.bhPosition,
          bhMass: imbhPhysics.bhMass,
          cloudInitialPosition: imbhPhysics.cloudInitialPosition,
          cloudInitialVelocity: imbhPhysics.cloudInitialVelocity
        }, 100);
        const orbitVectors = orbitPoints.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
        const orbitGeo = new THREE.BufferGeometry().setFromPoints(orbitVectors);
        const orbitMat = new THREE.LineDashedMaterial({
          color: 65416,
          // Green for analytical solution
          dashSize: 0.8,
          gapSize: 0.4,
          transparent: true,
          opacity: 0.6
        });
        const orbitLine = new THREE.Line(orbitGeo, orbitMat);
        orbitLine.computeLineDistances();
        scene.add(orbitLine);
        const straightLine = computeStraightLineTrajectory(
          {
            cloudInitialPosition: imbhPhysics.cloudInitialPosition,
            cloudInitialVelocity: imbhPhysics.cloudInitialVelocity
          },
          0,
          4,
          50
          // t from 0 to 4 code time units
        );
        const straightVectors = straightLine.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
        const straightGeo = new THREE.BufferGeometry().setFromPoints(straightVectors);
        const straightMat = new THREE.LineDashedMaterial({
          color: 8947848,
          dashSize: 0.5,
          gapSize: 0.3,
          transparent: true,
          opacity: 0.3
        });
        const straightLineObj = new THREE.Line(straightGeo, straightMat);
        straightLineObj.computeLineDistances();
        scene.add(straightLineObj);
      }
      if (showTrajectory) {
        const trajectoryGeo = new THREE.BufferGeometry();
        const trajectoryPositions = new Float32Array(200 * 3);
        trajectoryGeo.setAttribute("position", new THREE.BufferAttribute(trajectoryPositions, 3).setUsage(THREE.DynamicDrawUsage));
        trajectoryGeo.setDrawRange(0, 0);
        const trajectoryMat = new THREE.LineBasicMaterial({
          color: 16755200,
          transparent: true,
          opacity: 0.8,
          linewidth: 2
        });
        const trajectoryLine = new THREE.Line(trajectoryGeo, trajectoryMat);
        scene.add(trajectoryLine);
        trajectoryLineRef.current = trajectoryLine;
        const comGeo = new THREE.SphereGeometry(0.2, 16, 16);
        const comMat = new THREE.MeshBasicMaterial({ color: 16755200 });
        const comMarker = new THREE.Mesh(comGeo, comMat);
        comMarker.visible = false;
        scene.add(comMarker);
        comMarkerRef.current = comMarker;
      }
    }
    computeStats();
    if (boundingBox) {
      const center = new THREE.Vector3(
        (boundingBox.min[0] + boundingBox.max[0]) / 2,
        (boundingBox.min[1] + boundingBox.max[1]) / 2,
        (boundingBox.min[2] + boundingBox.max[2]) / 2
      );
      const size = Math.max(
        boundingBox.max[0] - boundingBox.min[0],
        boundingBox.max[1] - boundingBox.min[1],
        boundingBox.max[2] - boundingBox.min[2]
      );
      camera.position.set(center.x + size * 2, center.y + size, center.z + size * 2);
      camera.lookAt(center);
      controls.target.copy(center);
    }
    let animationId;
    const animate = () => {
      animationId = requestAnimationFrame(animate);
      fpsRef.current.frames++;
      const now = performance.now();
      if (now - fpsRef.current.lastTime >= 1e3) {
        fpsRef.current.fps = fpsRef.current.frames;
        fpsRef.current.frames = 0;
        fpsRef.current.lastTime = now;
        onFpsUpdate?.(fpsRef.current.fps);
      }
      updateParticles();
      updateTrajectory();
      controls.update();
      renderer.render(scene, camera);
    };
    animate();
    const handleResize = () => {
      const width2 = container.clientWidth;
      const height2 = container.clientHeight;
      camera.aspect = width2 / height2;
      camera.updateProjectionMatrix();
      renderer.setSize(width2, height2);
    };
    window.addEventListener("resize", handleResize);
    return () => {
      cancelAnimationFrame(animationId);
      window.removeEventListener("resize", handleResize);
      renderer.dispose();
      geometry.dispose();
      material.dispose();
      container.removeChild(renderer.domElement);
    };
  }, []);
  useEffect(() => {
    if (materialRef.current) {
      materialRef.current.size = pointSize;
      materialRef.current.opacity = opacity;
      materialRef.current.needsUpdate = true;
    }
  }, [pointSize, opacity]);
  useEffect(() => {
    if (bhGroupRef.current) {
      bhGroupRef.current.visible = showBlackHole;
    }
    if (trajectoryLineRef.current) {
      trajectoryLineRef.current.visible = showTrajectory;
    }
    if (comMarkerRef.current) {
      comMarkerRef.current.visible = showTrajectory;
    }
    if (tidalCircleRef.current) {
      tidalCircleRef.current.visible = showRadii;
    }
    if (hillCircleRef.current) {
      hillCircleRef.current.visible = showRadii;
    }
  }, [showBlackHole, showTrajectory, showRadii]);
  useEffect(() => {
    computeStats();
  }, [framesRef.current?.size]);
  return /* @__PURE__ */ jsx("div", { ref: containerRef, className: `w-full h-full ${className}` });
}
const defaultColorMap = {
  name: "Viridis",
  colors: [
    "#440154",
    "#482878",
    "#3e4989",
    "#31688e",
    "#26828e",
    "#1f9e89",
    "#35b779",
    "#6ece58",
    "#b5de2b",
    "#fde725"
  ],
  min: 0,
  max: 1,
  logScale: false
};
function interpolateColorHex(colors, t) {
  if (colors.length === 0) return "#ffffff";
  if (colors.length === 1) return colors[0];
  t = Math.max(0, Math.min(1, t));
  const index = t * (colors.length - 1);
  const lower = Math.floor(index);
  const upper = Math.min(lower + 1, colors.length - 1);
  const localT = index - lower;
  const c1 = hexToRgb(colors[lower]);
  const c2 = hexToRgb(colors[upper]);
  const r = Math.round(c1.r + (c2.r - c1.r) * localT);
  const g = Math.round(c1.g + (c2.g - c1.g) * localT);
  const b = Math.round(c1.b + (c2.b - c1.b) * localT);
  return `rgb(${r},${g},${b})`;
}
function hexToRgb(hex) {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!result) return { r: 255, g: 255, b: 255 };
  return {
    r: parseInt(result[1], 16),
    g: parseInt(result[2], 16),
    b: parseInt(result[3], 16)
  };
}
function Projection2D({
  frame,
  projection,
  colorField = "density",
  colorMap = defaultColorMap,
  width = 400,
  height = 400,
  showAxes = true,
  showColorbar = true,
  className = ""
}) {
  const canvasRef = useRef(null);
  const [axisX, axisY, axisLabel] = useMemo(() => {
    switch (projection) {
      case "xy":
        return [0, 1, ["X", "Y"]];
      case "xz":
        return [0, 2, ["X", "Z"]];
      case "yz":
        return [1, 2, ["Y", "Z"]];
      default:
        return [0, 1, ["X", "Y"]];
    }
  }, [projection]);
  const { fieldData, bounds, range } = useMemo(() => {
    if (!frame) {
      return {
        fieldData: null,
        bounds: { minX: -1, maxX: 1, minY: -1, maxY: 1 },
        range: { min: 0, max: 1 }
      };
    }
    let data;
    switch (colorField) {
      case "density":
        data = frame.density;
        break;
      case "pressure":
        data = frame.pressure;
        break;
      case "energy":
        data = frame.energy;
        break;
      case "velocity":
        data = new Float32Array(frame.particleCount);
        for (let i = 0; i < frame.particleCount; i++) {
          const vx = frame.velocities[i * 3];
          const vy = frame.velocities[i * 3 + 1];
          const vz = frame.velocities[i * 3 + 2];
          data[i] = Math.sqrt(vx * vx + vy * vy + vz * vz);
        }
        break;
      case "machNumber":
        data = frame.machNumber;
        break;
      default:
        data = frame.density;
    }
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minVal = Infinity, maxVal = -Infinity;
    for (let i = 0; i < frame.particleCount; i++) {
      const x = frame.positions[i * 3 + axisX];
      const y = frame.positions[i * 3 + axisY];
      const val = data ? data[i] : 0;
      if (isFinite(x)) {
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
      }
      if (isFinite(y)) {
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
      }
      if (isFinite(val)) {
        if (val < minVal) minVal = val;
        if (val > maxVal) maxVal = val;
      }
    }
    const padX = (maxX - minX) * 0.05;
    const padY = (maxY - minY) * 0.05;
    minX -= padX;
    maxX += padX;
    minY -= padY;
    maxY += padY;
    return {
      fieldData: data,
      bounds: { minX, maxX, minY, maxY },
      range: {
        min: colorMap.min !== void 0 ? colorMap.min : minVal,
        max: colorMap.max !== void 0 ? colorMap.max : maxVal
      }
    };
  }, [frame, colorField, colorMap, axisX, axisY]);
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas || !frame) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    ctx.scale(dpr, dpr);
    ctx.fillStyle = "#1a1a2e";
    ctx.fillRect(0, 0, width, height);
    const margin = showAxes ? 40 : 10;
    const plotWidth = width - margin * 2;
    const plotHeight = height - margin * 2;
    const scaleX = (val) => margin + (val - bounds.minX) / (bounds.maxX - bounds.minX) * plotWidth;
    const scaleY = (val) => height - margin - (val - bounds.minY) / (bounds.maxY - bounds.minY) * plotHeight;
    if (showAxes) {
      ctx.strokeStyle = "#333";
      ctx.lineWidth = 0.5;
      ctx.beginPath();
      const numGridLines = 5;
      for (let i = 0; i <= numGridLines; i++) {
        const x = margin + i / numGridLines * plotWidth;
        ctx.moveTo(x, margin);
        ctx.lineTo(x, height - margin);
      }
      for (let i = 0; i <= numGridLines; i++) {
        const y = margin + i / numGridLines * plotHeight;
        ctx.moveTo(margin, y);
        ctx.lineTo(width - margin, y);
      }
      ctx.stroke();
    }
    const { logScale } = colorMap;
    const { min, max } = range;
    for (let i = 0; i < frame.particleCount; i++) {
      const x = frame.positions[i * 3 + axisX];
      const y = frame.positions[i * 3 + axisY];
      const val = fieldData ? fieldData[i] : 0;
      if (!isFinite(x) || !isFinite(y)) continue;
      let t;
      if (logScale && min > 0 && val > 0) {
        t = (Math.log10(val) - Math.log10(min)) / (Math.log10(max) - Math.log10(min));
      } else {
        t = (val - min) / (max - min);
      }
      t = Math.max(0, Math.min(1, isFinite(t) ? t : 0));
      const color = interpolateColorHex(colorMap.colors, t);
      const px = scaleX(x);
      const py = scaleY(y);
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(px, py, 1.5, 0, Math.PI * 2);
      ctx.fill();
    }
    if (showAxes) {
      ctx.strokeStyle = "#888";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(margin, height - margin);
      ctx.lineTo(width - margin, height - margin);
      ctx.moveTo(margin, height - margin);
      ctx.lineTo(margin, margin);
      ctx.stroke();
      ctx.fillStyle = "#ccc";
      ctx.font = "12px sans-serif";
      ctx.textAlign = "center";
      ctx.fillText(axisLabel[0], width / 2, height - 5);
      ctx.save();
      ctx.translate(12, height / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(axisLabel[1], 0, 0);
      ctx.restore();
      ctx.font = "10px sans-serif";
      ctx.fillText(bounds.minX.toFixed(2), margin, height - 25);
      ctx.fillText(bounds.maxX.toFixed(2), width - margin, height - 25);
      ctx.textAlign = "right";
      ctx.fillText(bounds.minY.toFixed(2), margin - 5, height - margin);
      ctx.fillText(bounds.maxY.toFixed(2), margin - 5, margin + 5);
    }
    if (showColorbar) {
      const barWidth = 15;
      const barHeight = plotHeight * 0.6;
      const barX = width - margin + 10;
      const barY = margin + (plotHeight - barHeight) / 2;
      for (let i = 0; i < barHeight; i++) {
        const t = 1 - i / barHeight;
        ctx.fillStyle = interpolateColorHex(colorMap.colors, t);
        ctx.fillRect(barX, barY + i, barWidth, 1);
      }
      ctx.strokeStyle = "#888";
      ctx.strokeRect(barX, barY, barWidth, barHeight);
      ctx.fillStyle = "#ccc";
      ctx.font = "9px sans-serif";
      ctx.textAlign = "left";
      const formatVal = (v) => {
        if (Math.abs(v) < 1e-3 || Math.abs(v) > 1e4) {
          return v.toExponential(1);
        }
        return v.toFixed(3);
      };
      ctx.fillText(formatVal(range.max), barX + barWidth + 3, barY + 8);
      ctx.fillText(formatVal(range.min), barX + barWidth + 3, barY + barHeight);
    }
  }, [
    frame,
    fieldData,
    bounds,
    range,
    colorMap,
    width,
    height,
    showAxes,
    showColorbar,
    axisX,
    axisY,
    axisLabel
  ]);
  useEffect(() => {
    draw();
  }, [draw]);
  return /* @__PURE__ */ jsxs("div", { className: `relative ${className}`, children: [
    /* @__PURE__ */ jsx(
      "canvas",
      {
        ref: canvasRef,
        style: { width, height },
        className: "rounded"
      }
    ),
    /* @__PURE__ */ jsxs("div", { className: "absolute top-1 left-1 text-white text-xs bg-black/50 px-1 rounded", children: [
      projection.toUpperCase(),
      " Projection"
    ] })
  ] });
}
const DEFAULT_CONFIG = {
  bhMass: 1,
  // 10^5 M_sun
  cloudRadius: 1.13,
  // pc
  cloudMass: 1,
  // 1000 M_sun
  impactParameter: 5.16,
  // pc
  pericentre: 1.69,
  // pc
  eccentricity: 1.24,
  // hyperbolic
  initialPosition: [20, -5.17, 0],
  // pc
  initialVelocity: [-9.35, 4.48, 0],
  // km/s
  inclination: 70,
  // degrees
  positionAngle: 41.6,
  // degrees
  lsrVelocity: -120,
  // km/s
  distanceToGC: 8
  // kpc
};
function OrbitalGeometryPanel({
  config: configOverride,
  className = "",
  width = 600,
  height = 500,
  showLabels = true,
  viewMode = "3d"
}) {
  const containerRef = useRef(null);
  const rendererRef = useRef(null);
  const config = useMemo(() => ({
    ...DEFAULT_CONFIG,
    ...configOverride
  }), [configOverride]);
  useEffect(() => {
    if (!containerRef.current) return;
    const container = containerRef.current;
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(1710638);
    const camera = new THREE.PerspectiveCamera(50, width / height, 0.1, 1e3);
    switch (viewMode) {
      case "orbital-plane":
        camera.position.set(0, 0, 50);
        break;
      case "observer":
        const inc2 = config.inclination * Math.PI / 180;
        const pa2 = config.positionAngle * Math.PI / 180;
        camera.position.set(
          50 * Math.sin(inc2) * Math.cos(pa2),
          50 * Math.sin(inc2) * Math.sin(pa2),
          50 * Math.cos(inc2)
        );
        break;
      default:
        camera.position.set(35, 25, 40);
    }
    camera.lookAt(0, 0, 0);
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;
    const controls = new OrbitControls$1(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.target.set(0, 0, 0);
    const createTextSprite = (text, color = "#ffffff", fontSize = 32, bgColor) => {
      const canvas = document.createElement("canvas");
      const size = 256;
      canvas.width = size;
      canvas.height = size / 2;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = color;
      ctx.font = `Bold ${fontSize}px Arial`;
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, size / 2, size / 4);
      const texture = new THREE.CanvasTexture(canvas);
      const material = new THREE.SpriteMaterial({ map: texture, transparent: true });
      const sprite = new THREE.Sprite(material);
      sprite.scale.set(6, 3, 1);
      return sprite;
    };
    const createArrow = (from, to, color, headLength = 0.5, headWidth = 0.3) => {
      const group = new THREE.Group();
      const direction = to.clone().sub(from).normalize();
      const length = from.distanceTo(to);
      const shaftGeo = new THREE.CylinderGeometry(0.05, 0.05, length - headLength, 8);
      const shaftMat = new THREE.MeshBasicMaterial({ color });
      const shaft = new THREE.Mesh(shaftGeo, shaftMat);
      shaft.position.copy(from.clone().lerp(to, (length - headLength) / 2 / length));
      shaft.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);
      group.add(shaft);
      const headGeo = new THREE.ConeGeometry(headWidth, headLength, 8);
      const headMat = new THREE.MeshBasicMaterial({ color });
      const head = new THREE.Mesh(headGeo, headMat);
      head.position.copy(to.clone().sub(direction.clone().multiplyScalar(headLength / 2)));
      head.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);
      group.add(head);
      return group;
    };
    const createCircle = (radius, color, dashed = false, segments = 64) => {
      const points = [];
      for (let i = 0; i <= segments; i++) {
        const theta = i / segments * Math.PI * 2;
        points.push(new THREE.Vector3(radius * Math.cos(theta), radius * Math.sin(theta), 0));
      }
      const geo = new THREE.BufferGeometry().setFromPoints(points);
      const mat = dashed ? new THREE.LineDashedMaterial({ color, dashSize: 0.3, gapSize: 0.15 }) : new THREE.LineBasicMaterial({ color });
      const line = new THREE.Line(geo, mat);
      if (dashed) line.computeLineDistances();
      return line;
    };
    const axesGroup = new THREE.Group();
    const xAxis = createArrow(
      new THREE.Vector3(-15, 0, 0),
      new THREE.Vector3(15, 0, 0),
      16729156
    );
    axesGroup.add(xAxis);
    const yAxis = createArrow(
      new THREE.Vector3(0, -10, 0),
      new THREE.Vector3(0, 15, 0),
      4521796
    );
    axesGroup.add(yAxis);
    const zAxis = createArrow(
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, 0, 12),
      4491519
    );
    axesGroup.add(zAxis);
    if (showLabels) {
      const xLabel = createTextSprite("X (pc)", "#ff6666");
      xLabel.position.set(16, 0, 0);
      axesGroup.add(xLabel);
      const yLabel = createTextSprite("Y (pc)", "#66ff66");
      yLabel.position.set(0, 16, 0);
      axesGroup.add(yLabel);
      const zLabel = createTextSprite("Z (⊥ plane)", "#6688ff");
      zLabel.position.set(0, 0, 13);
      axesGroup.add(zLabel);
    }
    scene.add(axesGroup);
    const planeGeo = new THREE.PlaneGeometry(50, 50);
    const planeMat = new THREE.MeshBasicMaterial({
      color: 3368703,
      transparent: true,
      opacity: 0.1,
      side: THREE.DoubleSide
    });
    const orbitalPlane = new THREE.Mesh(planeGeo, planeMat);
    scene.add(orbitalPlane);
    const planeBorder = createCircle(20, 3368703, true);
    scene.add(planeBorder);
    if (showLabels) {
      const planeLabel = createTextSprite("Orbital Plane", "#6699ff", 24);
      planeLabel.position.set(12, -12, 0.1);
      scene.add(planeLabel);
    }
    const bhGroup = new THREE.Group();
    const glowGeo = new THREE.SphereGeometry(0.8, 32, 32);
    const glowMat = new THREE.MeshBasicMaterial({
      color: 16729156,
      transparent: true,
      opacity: 0.3
    });
    bhGroup.add(new THREE.Mesh(glowGeo, glowMat));
    const coreGeo = new THREE.SphereGeometry(0.3, 32, 32);
    const coreMat = new THREE.MeshBasicMaterial({ color: 16777215 });
    bhGroup.add(new THREE.Mesh(coreGeo, coreMat));
    const horizonGeo = new THREE.SphereGeometry(0.15, 32, 32);
    const horizonMat = new THREE.MeshBasicMaterial({ color: 0 });
    bhGroup.add(new THREE.Mesh(horizonGeo, horizonMat));
    scene.add(bhGroup);
    if (showLabels) {
      const bhLabel = createTextSprite("IMBH", "#ff8888", 28);
      bhLabel.position.set(0, -1.5, 0);
      scene.add(bhLabel);
      const bhMassLabel = createTextSprite(`10⁵ M☉`, "#ffaaaa", 20);
      bhMassLabel.position.set(0, -2.5, 0);
      scene.add(bhMassLabel);
    }
    const { impactParameter: b, pericentre: rp, eccentricity: e } = config;
    const a = rp / (e - 1);
    const p = a * (e * e - 1);
    const thetaMax = Math.PI - Math.acos(1 / e);
    const orbitPoints = [];
    const numPoints = 100;
    for (let i = 0; i < numPoints; i++) {
      const theta = -thetaMax + 0.05 + (2 * thetaMax - 0.1) * (i / (numPoints - 1));
      const r = p / (1 + e * Math.cos(theta));
      const x = r * Math.cos(theta);
      const y = r * Math.sin(theta);
      if (Math.abs(x) < 30 && Math.abs(y) < 30) {
        orbitPoints.push(new THREE.Vector3(x, y, 0));
      }
    }
    const orbitGeo = new THREE.BufferGeometry().setFromPoints(orbitPoints);
    const orbitMat = new THREE.LineBasicMaterial({ color: 65416, linewidth: 2 });
    const orbitLine = new THREE.Line(orbitGeo, orbitMat);
    scene.add(orbitLine);
    const periMarker = new THREE.Group();
    const periGeo = new THREE.SphereGeometry(0.2, 16, 16);
    const periMat = new THREE.MeshBasicMaterial({ color: 16776960 });
    periMarker.add(new THREE.Mesh(periGeo, periMat));
    const periLineGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(rp, 0, 0)
    ]);
    const periLineMat = new THREE.LineDashedMaterial({
      color: 16776960,
      dashSize: 0.2,
      gapSize: 0.1
    });
    const periLine = new THREE.Line(periLineGeo, periLineMat);
    periLine.computeLineDistances();
    periMarker.add(periLine);
    periMarker.position.set(rp, 0, 0);
    scene.add(periMarker);
    if (showLabels) {
      const periLabel = createTextSprite("Pericentre", "#ffff66", 22);
      periLabel.position.set(rp + 0.5, 1.2, 0);
      scene.add(periLabel);
      const periValue = createTextSprite(`r_p = ${rp.toFixed(2)} pc`, "#ffff99", 18);
      periValue.position.set(rp + 0.5, 0.3, 0);
      scene.add(periValue);
    }
    const asymptoteStart = new THREE.Vector3(-25, b, 0);
    const asymptoteEnd = new THREE.Vector3(5, b, 0);
    const asymptoteGeo = new THREE.BufferGeometry().setFromPoints([asymptoteStart, asymptoteEnd]);
    const asymptoteMat = new THREE.LineDashedMaterial({
      color: 16746496,
      dashSize: 0.5,
      gapSize: 0.25
    });
    const asymptoteLine = new THREE.Line(asymptoteGeo, asymptoteMat);
    asymptoteLine.computeLineDistances();
    scene.add(asymptoteLine);
    const impactLineGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, b, 0)
    ]);
    const impactLineMat = new THREE.LineBasicMaterial({ color: 16746496 });
    const impactLine = new THREE.Line(impactLineGeo, impactLineMat);
    scene.add(impactLine);
    const bracketSize = 0.3;
    const bracketGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(-bracketSize, 0, 0),
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(0, bracketSize, 0)
    ]);
    const bracketMat = new THREE.LineBasicMaterial({ color: 16746496 });
    scene.add(new THREE.Line(bracketGeo, bracketMat));
    if (showLabels) {
      const bLabel = createTextSprite("b (impact)", "#ffaa44", 20);
      bLabel.position.set(-1.5, b / 2, 0);
      scene.add(bLabel);
      const bValue = createTextSprite(`= ${b.toFixed(2)} pc`, "#ffcc66", 18);
      bValue.position.set(-1.5, b / 2 - 1, 0);
      scene.add(bValue);
      const asymLabel = createTextSprite("Incoming asymptote", "#ff9944", 18);
      asymLabel.position.set(-15, b + 1.5, 0);
      scene.add(asymLabel);
    }
    const [x0, y0, z0] = config.initialPosition;
    const [vx0, vy0, vz0] = config.initialVelocity;
    const cloudGroup = new THREE.Group();
    const cloudGeo = new THREE.SphereGeometry(config.cloudRadius, 32, 32);
    const cloudMat = new THREE.MeshBasicMaterial({
      color: 4491519,
      transparent: true,
      opacity: 0.4
    });
    cloudGroup.add(new THREE.Mesh(cloudGeo, cloudMat));
    const cloudOutline = createCircle(config.cloudRadius, 6728447);
    cloudGroup.add(cloudOutline);
    cloudGroup.position.set(x0, y0, z0);
    scene.add(cloudGroup);
    const velScale = 0.8;
    const velArrow = createArrow(
      new THREE.Vector3(x0, y0, z0),
      new THREE.Vector3(x0 + vx0 * velScale, y0 + vy0 * velScale, z0 + vz0 * velScale),
      4500223,
      0.4,
      0.25
    );
    scene.add(velArrow);
    if (showLabels) {
      const cloudLabel = createTextSprite("Cloud (t=0)", "#66aaff", 22);
      cloudLabel.position.set(x0, y0 + config.cloudRadius + 1.5, z0);
      scene.add(cloudLabel);
      const posLabel = createTextSprite(
        `(${x0.toFixed(1)}, ${y0.toFixed(1)}) pc`,
        "#88ccff",
        16
      );
      posLabel.position.set(x0, y0 - config.cloudRadius - 1.5, z0);
      scene.add(posLabel);
      const velLabel = createTextSprite(`v₀`, "#44aaff", 20);
      velLabel.position.set(
        x0 + vx0 * velScale + 1,
        y0 + vy0 * velScale + 0.5,
        z0
      );
      scene.add(velLabel);
    }
    const inc = config.inclination * Math.PI / 180;
    const pa = config.positionAngle * Math.PI / 180;
    const losDir = new THREE.Vector3(
      Math.sin(inc) * Math.cos(pa),
      Math.sin(inc) * Math.sin(pa),
      Math.cos(inc)
    ).normalize();
    const losLength = 18;
    const losStart = losDir.clone().multiplyScalar(-5);
    const losEnd = losDir.clone().multiplyScalar(losLength);
    const losArrow = createArrow(losStart, losEnd, 16729343, 0.6, 0.35);
    scene.add(losArrow);
    const earthGroup = new THREE.Group();
    const earthGeo = new THREE.SphereGeometry(0.6, 32, 32);
    const earthMat = new THREE.MeshBasicMaterial({ color: 4491519 });
    earthGroup.add(new THREE.Mesh(earthGeo, earthMat));
    const landGeo = new THREE.SphereGeometry(0.62, 32, 32, 0, Math.PI, 0, Math.PI * 0.6);
    const landMat = new THREE.MeshBasicMaterial({ color: 4500036 });
    const landMesh = new THREE.Mesh(landGeo, landMat);
    landMesh.rotation.x = Math.PI / 4;
    earthGroup.add(landMesh);
    earthGroup.position.copy(losEnd);
    scene.add(earthGroup);
    if (showLabels) {
      const losLabel = createTextSprite("Line of Sight", "#ff88ff", 20);
      losLabel.position.copy(losDir.clone().multiplyScalar(losLength / 2).add(new THREE.Vector3(2, 0, 1)));
      scene.add(losLabel);
      const earthLabel = createTextSprite("Earth/ALMA", "#66aaff", 20);
      earthLabel.position.copy(losEnd.clone().add(new THREE.Vector3(2, 0, 1)));
      scene.add(earthLabel);
      const distLabel = createTextSprite(`d = ${config.distanceToGC} kpc`, "#88ccff", 16);
      distLabel.position.copy(losEnd.clone().add(new THREE.Vector3(2, 0, -0.5)));
      scene.add(distLabel);
    }
    const incArcPoints = [];
    const arcRadius = 8;
    for (let i = 0; i <= 20; i++) {
      const angle = i / 20 * inc;
      incArcPoints.push(new THREE.Vector3(
        arcRadius * Math.sin(angle) * Math.cos(pa),
        arcRadius * Math.sin(angle) * Math.sin(pa),
        arcRadius * Math.cos(angle)
      ));
    }
    const incArcGeo = new THREE.BufferGeometry().setFromPoints(incArcPoints);
    const incArcMat = new THREE.LineBasicMaterial({ color: 16755200 });
    const incArc = new THREE.Line(incArcGeo, incArcMat);
    scene.add(incArc);
    if (showLabels) {
      const incLabel = createTextSprite(`i = ${config.inclination}°`, "#ffcc44", 22);
      const midAngle = inc / 2;
      incLabel.position.set(
        (arcRadius + 2) * Math.sin(midAngle) * Math.cos(pa),
        (arcRadius + 2) * Math.sin(midAngle) * Math.sin(pa),
        (arcRadius + 2) * Math.cos(midAngle)
      );
      scene.add(incLabel);
    }
    const paArcPoints = [];
    const paRadius = 6;
    for (let i = 0; i <= 20; i++) {
      const angle = i / 20 * pa;
      paArcPoints.push(new THREE.Vector3(
        paRadius * Math.cos(angle) * Math.sin(inc * 0.3),
        // Projected
        paRadius * Math.sin(angle) * Math.sin(inc * 0.3),
        paRadius * Math.cos(inc * 0.3)
      ));
    }
    const paArcGeo = new THREE.BufferGeometry().setFromPoints(paArcPoints);
    const paArcMat = new THREE.LineDashedMaterial({ color: 4521898, dashSize: 0.2, gapSize: 0.1 });
    const paArc = new THREE.Line(paArcGeo, paArcMat);
    paArc.computeLineDistances();
    scene.add(paArc);
    if (showLabels) {
      const paLabel = createTextSprite(`PA = ${config.positionAngle}°`, "#66ffcc", 20);
      paLabel.position.set(paRadius + 2, 0, paRadius * Math.cos(inc * 0.3) + 1);
      scene.add(paLabel);
    }
    const lsrArrow = createArrow(
      new THREE.Vector3(-12, -8, 0),
      new THREE.Vector3(-12 + losDir.x * 4, -8 + losDir.y * 4, losDir.z * 4),
      16737928,
      0.3,
      0.2
    );
    scene.add(lsrArrow);
    if (showLabels) {
      const lsrLabel = createTextSprite("V_LSR", "#ff8899", 20);
      lsrLabel.position.set(-12, -9.5, 0);
      scene.add(lsrLabel);
      const lsrValue = createTextSprite(`= ${config.lsrVelocity} km/s`, "#ffaabb", 16);
      lsrValue.position.set(-12, -10.8, 0);
      scene.add(lsrValue);
    }
    const legendGroup = new THREE.Group();
    legendGroup.position.set(-22, 12, 0);
    const legendBg = new THREE.Mesh(
      new THREE.PlaneGeometry(12, 14),
      new THREE.MeshBasicMaterial({ color: 1710638, transparent: true, opacity: 0.9 })
    );
    legendGroup.add(legendBg);
    const legendTitle = createTextSprite("LEGEND", "#ffffff", 28);
    legendTitle.position.set(0, 5.5, 0.1);
    legendTitle.scale.set(4, 2, 1);
    legendGroup.add(legendTitle);
    const items = [
      { color: "#00ff88", text: "Hyperbolic orbit" },
      { color: "#ff8800", text: "Impact parameter b" },
      { color: "#ffff00", text: "Pericentre r_p" },
      { color: "#ff44ff", text: "Line of sight" },
      { color: "#ffaa00", text: "Inclination i" },
      { color: "#44ffaa", text: "Position angle PA" },
      { color: "#4488ff", text: "Cloud position" }
    ];
    items.forEach((item, i) => {
      const dot = new THREE.Mesh(
        new THREE.SphereGeometry(0.2, 8, 8),
        new THREE.MeshBasicMaterial({ color: item.color })
      );
      dot.position.set(-4.5, 4 - i * 1.3, 0.1);
      legendGroup.add(dot);
      const label = createTextSprite(item.text, item.color, 18);
      label.position.set(1, 4 - i * 1.3, 0.1);
      label.scale.set(4, 2, 1);
      legendGroup.add(label);
    });
    scene.add(legendGroup);
    let animationId;
    const animate = () => {
      animationId = requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    };
    animate();
    return () => {
      cancelAnimationFrame(animationId);
      renderer.dispose();
      container.removeChild(renderer.domElement);
    };
  }, [config, width, height, showLabels, viewMode]);
  return /* @__PURE__ */ jsxs("div", { className: `relative ${className}`, children: [
    /* @__PURE__ */ jsx(
      "div",
      {
        ref: containerRef,
        style: { width, height },
        className: "rounded-lg overflow-hidden"
      }
    ),
    /* @__PURE__ */ jsxs("div", { className: "absolute top-2 left-2 bg-black/70 text-white text-xs p-2 rounded", children: [
      /* @__PURE__ */ jsx("div", { className: "font-bold mb-1", children: "Orbital Geometry" }),
      /* @__PURE__ */ jsx("div", { children: "Oka et al. (2017) / CAT_OKA" })
    ] }),
    /* @__PURE__ */ jsx("div", { className: "absolute bottom-2 right-2 bg-black/70 text-white text-xs p-2 rounded max-w-xs", children: /* @__PURE__ */ jsx("table", { className: "w-full", children: /* @__PURE__ */ jsxs("tbody", { children: [
      /* @__PURE__ */ jsxs("tr", { children: [
        /* @__PURE__ */ jsx("td", { className: "text-yellow-400", children: "r_peri:" }),
        /* @__PURE__ */ jsxs("td", { className: "text-right", children: [
          config.pericentre.toFixed(2),
          " pc"
        ] })
      ] }),
      /* @__PURE__ */ jsxs("tr", { children: [
        /* @__PURE__ */ jsx("td", { className: "text-orange-400", children: "b:" }),
        /* @__PURE__ */ jsxs("td", { className: "text-right", children: [
          config.impactParameter.toFixed(2),
          " pc"
        ] })
      ] }),
      /* @__PURE__ */ jsxs("tr", { children: [
        /* @__PURE__ */ jsx("td", { className: "text-green-400", children: "e:" }),
        /* @__PURE__ */ jsx("td", { className: "text-right", children: config.eccentricity.toFixed(3) })
      ] }),
      /* @__PURE__ */ jsxs("tr", { children: [
        /* @__PURE__ */ jsx("td", { className: "text-purple-400", children: "i:" }),
        /* @__PURE__ */ jsxs("td", { className: "text-right", children: [
          config.inclination,
          "°"
        ] })
      ] }),
      /* @__PURE__ */ jsxs("tr", { children: [
        /* @__PURE__ */ jsx("td", { className: "text-teal-400", children: "PA:" }),
        /* @__PURE__ */ jsxs("td", { className: "text-right", children: [
          config.positionAngle,
          "°"
        ] })
      ] }),
      /* @__PURE__ */ jsxs("tr", { children: [
        /* @__PURE__ */ jsx("td", { className: "text-pink-400", children: "V_LSR:" }),
        /* @__PURE__ */ jsxs("td", { className: "text-right", children: [
          config.lsrVelocity,
          " km/s"
        ] })
      ] })
    ] }) }) })
  ] });
}
function EnergyChart({
  statistics,
  currentFrame,
  showKinetic = true,
  showInternal = true,
  showTotal = true,
  className = ""
}) {
  const data = statistics.filter((s) => s !== void 0 && s !== null).map((s) => ({
    time: s.time,
    kinetic: s.totalKineticEnergy,
    internal: s.totalInternalEnergy,
    total: s.totalEnergy
  }));
  const currentTime = currentFrame !== void 0 && statistics[currentFrame] ? statistics[currentFrame].time : void 0;
  return /* @__PURE__ */ jsxs("div", { className: `bg-gray-800 p-3 rounded ${className}`, children: [
    /* @__PURE__ */ jsx("h3", { className: "text-sm font-medium text-gray-300 mb-2", children: "Energy Evolution" }),
    /* @__PURE__ */ jsx(ResponsiveContainer, { width: "100%", height: 200, children: /* @__PURE__ */ jsxs(LineChart, { data, children: [
      /* @__PURE__ */ jsx(CartesianGrid, { strokeDasharray: "3 3", stroke: "#444" }),
      /* @__PURE__ */ jsx(
        XAxis,
        {
          dataKey: "time",
          stroke: "#888",
          tick: { fill: "#888", fontSize: 10 },
          tickFormatter: (v) => typeof v === "number" && isFinite(v) ? v.toFixed(2) : "0"
        }
      ),
      /* @__PURE__ */ jsx(
        YAxis,
        {
          stroke: "#888",
          tick: { fill: "#888", fontSize: 10 },
          tickFormatter: (v) => typeof v === "number" && isFinite(v) ? v.toExponential(1) : "0"
        }
      ),
      /* @__PURE__ */ jsx(
        Tooltip,
        {
          contentStyle: { backgroundColor: "#1f2937", border: "none" },
          labelStyle: { color: "#9ca3af" },
          formatter: (value) => typeof value === "number" && isFinite(value) ? value.toExponential(3) : "0"
        }
      ),
      /* @__PURE__ */ jsx(Legend, { wrapperStyle: { fontSize: 10 } }),
      showKinetic && /* @__PURE__ */ jsx(
        Line,
        {
          type: "monotone",
          dataKey: "kinetic",
          name: "Kinetic",
          stroke: "#3b82f6",
          strokeWidth: 1.5,
          dot: false
        }
      ),
      showInternal && /* @__PURE__ */ jsx(
        Line,
        {
          type: "monotone",
          dataKey: "internal",
          name: "Internal",
          stroke: "#ef4444",
          strokeWidth: 1.5,
          dot: false
        }
      ),
      showTotal && /* @__PURE__ */ jsx(
        Line,
        {
          type: "monotone",
          dataKey: "total",
          name: "Total",
          stroke: "#22c55e",
          strokeWidth: 2,
          dot: false
        }
      ),
      currentTime !== void 0 && /* @__PURE__ */ jsx(ReferenceLine, { x: currentTime, stroke: "#fff", strokeDasharray: "3 3" })
    ] }) })
  ] });
}
function MomentumChart({
  statistics,
  currentFrame,
  className = ""
}) {
  const data = statistics.filter((s) => s !== void 0 && s !== null).map((s) => ({
    time: s.time,
    px: s.momentum[0],
    py: s.momentum[1],
    pz: s.momentum[2],
    total: Math.sqrt(s.momentum[0] ** 2 + s.momentum[1] ** 2 + s.momentum[2] ** 2)
  }));
  const currentTime = currentFrame !== void 0 && statistics[currentFrame] ? statistics[currentFrame].time : void 0;
  return /* @__PURE__ */ jsxs("div", { className: `bg-gray-800 p-3 rounded ${className}`, children: [
    /* @__PURE__ */ jsx("h3", { className: "text-sm font-medium text-gray-300 mb-2", children: "Momentum Evolution" }),
    /* @__PURE__ */ jsx(ResponsiveContainer, { width: "100%", height: 200, children: /* @__PURE__ */ jsxs(LineChart, { data, children: [
      /* @__PURE__ */ jsx(CartesianGrid, { strokeDasharray: "3 3", stroke: "#444" }),
      /* @__PURE__ */ jsx(
        XAxis,
        {
          dataKey: "time",
          stroke: "#888",
          tick: { fill: "#888", fontSize: 10 },
          tickFormatter: (v) => typeof v === "number" && isFinite(v) ? v.toFixed(2) : "0"
        }
      ),
      /* @__PURE__ */ jsx(
        YAxis,
        {
          stroke: "#888",
          tick: { fill: "#888", fontSize: 10 },
          tickFormatter: (v) => typeof v === "number" && isFinite(v) ? v.toExponential(1) : "0"
        }
      ),
      /* @__PURE__ */ jsx(
        Tooltip,
        {
          contentStyle: { backgroundColor: "#1f2937", border: "none" },
          labelStyle: { color: "#9ca3af" },
          formatter: (value) => typeof value === "number" && isFinite(value) ? value.toExponential(3) : "0"
        }
      ),
      /* @__PURE__ */ jsx(Legend, { wrapperStyle: { fontSize: 10 } }),
      /* @__PURE__ */ jsx(Line, { type: "monotone", dataKey: "px", name: "Px", stroke: "#3b82f6", strokeWidth: 1, dot: false }),
      /* @__PURE__ */ jsx(Line, { type: "monotone", dataKey: "py", name: "Py", stroke: "#ef4444", strokeWidth: 1, dot: false }),
      /* @__PURE__ */ jsx(Line, { type: "monotone", dataKey: "pz", name: "Pz", stroke: "#22c55e", strokeWidth: 1, dot: false }),
      /* @__PURE__ */ jsx(Line, { type: "monotone", dataKey: "total", name: "|P|", stroke: "#f59e0b", strokeWidth: 2, dot: false }),
      currentTime !== void 0 && /* @__PURE__ */ jsx(ReferenceLine, { x: currentTime, stroke: "#fff", strokeDasharray: "3 3" })
    ] }) })
  ] });
}
function PlaybackControls({
  currentFrame,
  totalFrames,
  time,
  onFrameChange,
  onPlayPauseChange,
  isFrameReady,
  defaultSpeed = 30,
  className = "",
  imperativeMode = false,
  frameIndexRef,
  playbackSpeedRef
}) {
  const [isPlaying, setIsPlaying] = useState(false);
  const [speed, setSpeed] = useState(defaultSpeed);
  const [displayFrame, setDisplayFrame] = useState(currentFrame);
  useRef(null);
  useEffect(() => {
    if (imperativeMode && frameIndexRef) {
      const interval = setInterval(() => {
        if (frameIndexRef.current !== void 0) {
          setDisplayFrame(frameIndexRef.current);
        }
      }, 50);
      return () => clearInterval(interval);
    }
  }, [imperativeMode, frameIndexRef]);
  useEffect(() => {
    if (!imperativeMode) {
      setDisplayFrame(currentFrame);
    }
  }, [currentFrame, imperativeMode]);
  useEffect(() => {
    if (playbackSpeedRef && "current" in playbackSpeedRef) {
      playbackSpeedRef.current = speed;
    }
  }, [speed, playbackSpeedRef]);
  const play = useCallback(() => {
    setIsPlaying(true);
    onPlayPauseChange?.(true);
  }, [onPlayPauseChange]);
  const pause = useCallback(() => {
    setIsPlaying(false);
    onPlayPauseChange?.(false);
  }, [onPlayPauseChange]);
  const togglePlayPause = useCallback(() => {
    if (isPlaying) {
      pause();
    } else {
      play();
    }
  }, [isPlaying, play, pause]);
  const goToStart = useCallback(() => {
    pause();
    onFrameChange(0);
  }, [pause, onFrameChange]);
  const goToEnd = useCallback(() => {
    pause();
    onFrameChange(totalFrames - 1);
  }, [pause, onFrameChange, totalFrames]);
  const stepForward = useCallback(() => {
    const frame = imperativeMode && frameIndexRef ? frameIndexRef.current : currentFrame;
    if (frame < totalFrames - 1) {
      onFrameChange(frame + 1);
    }
  }, [currentFrame, totalFrames, onFrameChange, imperativeMode, frameIndexRef]);
  const stepBackward = useCallback(() => {
    const frame = imperativeMode && frameIndexRef ? frameIndexRef.current : currentFrame;
    if (frame > 0) {
      onFrameChange(frame - 1);
    }
  }, [currentFrame, onFrameChange, imperativeMode, frameIndexRef]);
  const currentFrameRef = useRef(currentFrame);
  useEffect(() => {
    currentFrameRef.current = currentFrame;
  }, [currentFrame]);
  const isFrameReadyRef = useRef(isFrameReady);
  useEffect(() => {
    isFrameReadyRef.current = isFrameReady;
  }, [isFrameReady]);
  useEffect(() => {
    if (imperativeMode) return;
    if (isPlaying) {
      let lastTime = 0;
      const frameInterval = 1e3 / speed;
      let animationId = null;
      const tick = (timestamp) => {
        if (!lastTime) lastTime = timestamp;
        const elapsed = timestamp - lastTime;
        if (elapsed >= frameInterval) {
          const nextFrame = currentFrameRef.current + 1;
          if (nextFrame >= totalFrames) {
            pause();
            onFrameChange(totalFrames - 1);
            return;
          }
          const frameReady = isFrameReadyRef.current ? isFrameReadyRef.current(nextFrame) : true;
          if (frameReady) {
            onFrameChange(nextFrame);
            lastTime = timestamp;
          }
        }
        animationId = requestAnimationFrame(tick);
      };
      animationId = requestAnimationFrame(tick);
      return () => {
        if (animationId !== null) {
          cancelAnimationFrame(animationId);
        }
      };
    }
  }, [isPlaying, speed, totalFrames, pause, onFrameChange, imperativeMode]);
  useEffect(() => {
    const handleKeyDown = (e) => {
      if (e.target instanceof HTMLInputElement) return;
      switch (e.key) {
        case " ":
          e.preventDefault();
          togglePlayPause();
          break;
        case "ArrowLeft":
          e.preventDefault();
          stepBackward();
          break;
        case "ArrowRight":
          e.preventDefault();
          stepForward();
          break;
        case "Home":
          e.preventDefault();
          goToStart();
          break;
        case "End":
          e.preventDefault();
          goToEnd();
          break;
      }
    };
    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, [togglePlayPause, stepBackward, stepForward, goToStart, goToEnd]);
  const formatTime = (t) => {
    if (Math.abs(t) < 1e-3) return t.toExponential(2);
    if (Math.abs(t) > 1e3) return t.toExponential(2);
    return t.toFixed(4);
  };
  return /* @__PURE__ */ jsx("div", { className: `bg-gray-800 p-3 rounded ${className}`, children: /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-4", children: [
    /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-1", children: [
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: goToStart,
          className: "p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white",
          title: "Go to start (Home)",
          children: /* @__PURE__ */ jsx(SkipBack, { size: 16 })
        }
      ),
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: stepBackward,
          className: "p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white",
          title: "Step backward (←)",
          disabled: displayFrame === 0,
          children: /* @__PURE__ */ jsx(ChevronLeft, { size: 16 })
        }
      ),
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: togglePlayPause,
          className: "p-2 bg-blue-600 hover:bg-blue-500 rounded text-white",
          title: isPlaying ? "Pause (Space)" : "Play (Space)",
          children: isPlaying ? /* @__PURE__ */ jsx(Pause, { size: 18 }) : /* @__PURE__ */ jsx(Play, { size: 18 })
        }
      ),
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: stepForward,
          className: "p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white",
          title: "Step forward (→)",
          disabled: displayFrame === totalFrames - 1,
          children: /* @__PURE__ */ jsx(ChevronRight, { size: 16 })
        }
      ),
      /* @__PURE__ */ jsx(
        "button",
        {
          onClick: goToEnd,
          className: "p-1.5 hover:bg-gray-700 rounded text-gray-400 hover:text-white",
          title: "Go to end (End)",
          children: /* @__PURE__ */ jsx(SkipForward, { size: 16 })
        }
      )
    ] }),
    /* @__PURE__ */ jsx("div", { className: "flex-1", children: /* @__PURE__ */ jsx(
      "input",
      {
        type: "range",
        min: 0,
        max: totalFrames - 1,
        value: displayFrame,
        onChange: (e) => onFrameChange(parseInt(e.target.value)),
        className: "w-full h-2 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
      }
    ) }),
    /* @__PURE__ */ jsxs("div", { className: "text-xs text-gray-400 whitespace-nowrap min-w-[120px]", children: [
      /* @__PURE__ */ jsx("span", { className: "font-mono", children: displayFrame + 1 }),
      /* @__PURE__ */ jsx("span", { className: "text-gray-600", children: " / " }),
      /* @__PURE__ */ jsx("span", { className: "font-mono", children: totalFrames }),
      /* @__PURE__ */ jsx("span", { className: "text-gray-600 ml-2", children: "t = " }),
      /* @__PURE__ */ jsx("span", { className: "font-mono", children: formatTime(time) })
    ] }),
    /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-2", children: [
      /* @__PURE__ */ jsx("span", { className: "text-xs text-gray-500", children: "Speed:" }),
      /* @__PURE__ */ jsxs(
        "select",
        {
          value: speed,
          onChange: (e) => setSpeed(parseInt(e.target.value)),
          className: "bg-gray-700 text-white text-xs rounded px-1 py-0.5",
          children: [
            /* @__PURE__ */ jsx("option", { value: 5, children: "5 fps" }),
            /* @__PURE__ */ jsx("option", { value: 10, children: "10 fps" }),
            /* @__PURE__ */ jsx("option", { value: 15, children: "15 fps" }),
            /* @__PURE__ */ jsx("option", { value: 20, children: "20 fps" }),
            /* @__PURE__ */ jsx("option", { value: 30, children: "30 fps" }),
            /* @__PURE__ */ jsx("option", { value: 45, children: "45 fps" }),
            /* @__PURE__ */ jsx("option", { value: 60, children: "60 fps" })
          ]
        }
      )
    ] })
  ] }) });
}
const COLOR_MAPS = {
  viridis: {
    name: "Viridis",
    colors: ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"],
    logScale: false
  },
  plasma: {
    name: "Plasma",
    colors: ["#0d0887", "#46039f", "#7201a8", "#9c179e", "#bd3786", "#d8576b", "#ed7953", "#fb9f3a", "#fdca26", "#f0f921"],
    logScale: false
  },
  inferno: {
    name: "Inferno",
    colors: ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#f7d13d", "#fcffa4"],
    logScale: false
  },
  coolwarm: {
    name: "Cool-Warm",
    colors: ["#3b4cc0", "#6688ee", "#88bbff", "#b8d4eb", "#dddddd", "#f5c4ad", "#f49a7b", "#d6604d", "#b40426"],
    logScale: false
  },
  density: {
    name: "Density",
    colors: ["#000033", "#000066", "#000099", "#0033cc", "#0066ff", "#00ccff", "#66ffcc", "#ccff66", "#ffcc00", "#ff6600", "#ff0000"],
    logScale: true
  }
};
const DEFAULT_FIELDS = ["density", "pressure", "energy", "velocity", "machNumber"];
function VisualizationSettings({
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
  className = ""
}) {
  const [expanded, setExpanded] = useState(true);
  return /* @__PURE__ */ jsxs("div", { className: `bg-gray-800 rounded ${className}`, children: [
    /* @__PURE__ */ jsxs(
      "button",
      {
        onClick: () => setExpanded(!expanded),
        className: "w-full flex items-center gap-2 p-3 text-left hover:bg-gray-700/50",
        children: [
          /* @__PURE__ */ jsx(Settings, { size: 16, className: "text-gray-400" }),
          /* @__PURE__ */ jsx("span", { className: "text-sm font-medium text-gray-300 flex-1", children: "Visualization Settings" }),
          expanded ? /* @__PURE__ */ jsx(ChevronDown, { size: 16, className: "text-gray-400" }) : /* @__PURE__ */ jsx(ChevronRight, { size: 16, className: "text-gray-400" })
        ]
      }
    ),
    expanded && /* @__PURE__ */ jsxs("div", { className: "px-3 pb-3 space-y-3", children: [
      /* @__PURE__ */ jsxs("div", { children: [
        /* @__PURE__ */ jsx("label", { className: "text-xs text-gray-400 block mb-1", children: "Color Field" }),
        /* @__PURE__ */ jsx(
          "select",
          {
            value: colorField,
            onChange: (e) => onColorFieldChange(e.target.value),
            className: "w-full bg-gray-700 text-white text-sm rounded px-2 py-1.5",
            children: availableFields.map((field) => /* @__PURE__ */ jsx("option", { value: field, children: fieldDisplayName(field) }, field))
          }
        )
      ] }),
      /* @__PURE__ */ jsxs("div", { children: [
        /* @__PURE__ */ jsx("label", { className: "text-xs text-gray-400 block mb-1", children: "Color Map" }),
        /* @__PURE__ */ jsx(
          "select",
          {
            value: colorMapName,
            onChange: (e) => onColorMapChange(e.target.value),
            className: "w-full bg-gray-700 text-white text-sm rounded px-2 py-1.5",
            children: Object.keys(COLOR_MAPS).map((name) => /* @__PURE__ */ jsx("option", { value: name, children: COLOR_MAPS[name].name }, name))
          }
        ),
        /* @__PURE__ */ jsx(
          "div",
          {
            className: "h-3 mt-1 rounded",
            style: {
              background: `linear-gradient(to right, ${COLOR_MAPS[colorMapName]?.colors.join(", ")})`
            }
          }
        )
      ] }),
      /* @__PURE__ */ jsxs("div", { children: [
        /* @__PURE__ */ jsx("label", { className: "text-xs text-gray-400 block mb-1", children: "Color Range" }),
        /* @__PURE__ */ jsxs("div", { className: "flex gap-2 items-center", children: [
          /* @__PURE__ */ jsx(
            "input",
            {
              type: "number",
              value: colorRange[0],
              onChange: (e) => onColorRangeChange([parseFloat(e.target.value), colorRange[1]]),
              className: "w-20 bg-gray-700 text-white text-xs rounded px-2 py-1",
              step: "any"
            }
          ),
          /* @__PURE__ */ jsx("span", { className: "text-gray-500", children: "—" }),
          /* @__PURE__ */ jsx(
            "input",
            {
              type: "number",
              value: colorRange[1],
              onChange: (e) => onColorRangeChange([colorRange[0], parseFloat(e.target.value)]),
              className: "w-20 bg-gray-700 text-white text-xs rounded px-2 py-1",
              step: "any"
            }
          ),
          /* @__PURE__ */ jsx(
            "button",
            {
              onClick: () => onColorRangeChange([0, 0]),
              className: "text-xs text-blue-400 hover:text-blue-300",
              children: "Auto"
            }
          )
        ] })
      ] }),
      /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-2", children: [
        /* @__PURE__ */ jsx(
          "input",
          {
            type: "checkbox",
            id: "logScale",
            checked: useLogScale,
            onChange: (e) => onLogScaleChange(e.target.checked),
            className: "rounded bg-gray-700"
          }
        ),
        /* @__PURE__ */ jsx("label", { htmlFor: "logScale", className: "text-xs text-gray-400", children: "Logarithmic Scale" })
      ] }),
      /* @__PURE__ */ jsxs("div", { children: [
        /* @__PURE__ */ jsxs("label", { className: "text-xs text-gray-400 block mb-1", children: [
          "Point Size: ",
          pointSize.toFixed(3)
        ] }),
        /* @__PURE__ */ jsx(
          "input",
          {
            type: "range",
            min: 1e-3,
            max: 0.1,
            step: 1e-3,
            value: pointSize,
            onChange: (e) => onPointSizeChange(parseFloat(e.target.value)),
            className: "w-full h-1.5 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
          }
        )
      ] }),
      /* @__PURE__ */ jsxs("div", { children: [
        /* @__PURE__ */ jsxs("label", { className: "text-xs text-gray-400 block mb-1", children: [
          "Opacity: ",
          (opacity * 100).toFixed(0),
          "%"
        ] }),
        /* @__PURE__ */ jsx(
          "input",
          {
            type: "range",
            min: 0.1,
            max: 1,
            step: 0.05,
            value: opacity,
            onChange: (e) => onOpacityChange(parseFloat(e.target.value)),
            className: "w-full h-1.5 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
          }
        )
      ] }),
      /* @__PURE__ */ jsxs("div", { className: "space-y-2", children: [
        /* @__PURE__ */ jsx("label", { className: "text-xs text-gray-400 block", children: "Display Options" }),
        /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-2", children: [
          /* @__PURE__ */ jsx(
            "input",
            {
              type: "checkbox",
              id: "showAxes",
              checked: showAxes,
              onChange: (e) => onShowAxesChange(e.target.checked),
              className: "rounded bg-gray-700"
            }
          ),
          /* @__PURE__ */ jsx("label", { htmlFor: "showAxes", className: "text-xs text-gray-400", children: "Show Axes" })
        ] }),
        /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-2", children: [
          /* @__PURE__ */ jsx(
            "input",
            {
              type: "checkbox",
              id: "showBoundingBox",
              checked: showBoundingBox,
              onChange: (e) => onShowBoundingBoxChange(e.target.checked),
              className: "rounded bg-gray-700"
            }
          ),
          /* @__PURE__ */ jsx("label", { htmlFor: "showBoundingBox", className: "text-xs text-gray-400", children: "Show Bounding Box" })
        ] })
      ] })
    ] })
  ] });
}
function fieldDisplayName(field) {
  const names = {
    density: "Density (ρ)",
    pressure: "Pressure (P)",
    energy: "Internal Energy (u)",
    velocity: "Velocity (|v|)",
    machNumber: "Mach Number",
    smoothingLength: "Smoothing Length (h)",
    mass: "Mass (m)",
    divV: "Velocity Divergence (∇·v)",
    temperature: "Temperature (T)"
  };
  return names[field] || field;
}
function computeGlobalColorStats(frames, colorField) {
  if (frames.size === 0) return [0, 1];
  let globalMin = Infinity;
  let globalMax = -Infinity;
  let sampleCount = 0;
  const frameArray = Array.from(frames.values());
  const step = Math.max(1, Math.floor(frameArray.length / 20));
  for (let frameIdx = 0; frameIdx < frameArray.length; frameIdx += step) {
    const frame = frameArray[frameIdx];
    if (!frame || frame.particleCount === 0) continue;
    let fieldData;
    switch (colorField) {
      case "density":
        fieldData = frame.density;
        break;
      case "pressure":
        fieldData = frame.pressure;
        break;
      case "energy":
        fieldData = frame.energy;
        break;
      case "velocity":
        for (let i = 0; i < frame.particleCount; i += 50) {
          const vx = frame.velocities[i * 3];
          const vy = frame.velocities[i * 3 + 1];
          const vz = frame.velocities[i * 3 + 2];
          const vmag = Math.sqrt(vx * vx + vy * vy + vz * vz);
          if (isFinite(vmag) && vmag > 0) {
            if (vmag < globalMin) globalMin = vmag;
            if (vmag > globalMax) globalMax = vmag;
            sampleCount++;
          }
        }
        continue;
    }
    if (fieldData && fieldData.length > 0) {
      for (let i = 0; i < fieldData.length; i += 50) {
        const v = fieldData[i];
        if (isFinite(v) && v > 0) {
          if (v < globalMin) globalMin = v;
          if (v > globalMax) globalMax = v;
          sampleCount++;
        }
      }
    }
  }
  if (!isFinite(globalMin) || !isFinite(globalMax) || sampleCount === 0) {
    console.warn(`[computeGlobalColorStats] No valid data for ${colorField}, using defaults`);
    return [1e-3, 1];
  }
  if (globalMin === globalMax) {
    globalMax = globalMin * 1.1 + 1e-3;
  }
  const range = globalMax - globalMin;
  const paddedMin = Math.max(globalMin - range * 0.05, globalMin * 0.9);
  const paddedMax = globalMax + range * 0.05;
  console.log(`[computeGlobalColorStats] ${colorField}: [${paddedMin.toExponential(2)}, ${paddedMax.toExponential(2)}] from ${sampleCount} samples`);
  return [paddedMin, paddedMax];
}
function Dashboard({
  simulation,
  frames,
  statistics,
  onLoadFrame,
  isLoading,
  error
}) {
  const [currentFrameIndex, setCurrentFrameIndex] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [colorField, setColorField] = useState("density");
  const [colorMapName, setColorMapName] = useState("viridis");
  const [pointSize, setPointSize] = useState(0.02);
  const [opacity, setOpacity] = useState(0.8);
  const [showAxes, setShowAxes] = useState(true);
  const [showBoundingBox, setShowBoundingBox] = useState(true);
  const [colorRange, setColorRange] = useState([0, 0]);
  const [useLogScale, setUseLogScale] = useState(false);
  const [showProjections, setShowProjections] = useState(true);
  const [showCharts, setShowCharts] = useState(true);
  const [showOrbitalGeometry, setShowOrbitalGeometry] = useState(false);
  const [useImperativeMode, setUseImperativeMode] = useState(true);
  const [currentFps, setCurrentFps] = useState(0);
  const [showBlackHole, setShowBlackHole] = useState(true);
  const [showTrajectory, setShowTrajectory] = useState(true);
  const [showRadii, setShowRadii] = useState(true);
  const imbhPhysics = useMemo(() => {
    if (simulation?.imbhPhysics) {
      console.log("[Dashboard] Using IMBH physics from simulation config:", simulation.imbhPhysics);
      return simulation.imbhPhysics;
    }
    console.log("[Dashboard] Using default IMBH physics (no config found)");
    return {
      enabled: false,
      // Disabled by default if no config
      bhPosition: [0, 0, 0],
      bhMass: 100,
      cloudInitialPosition: [-20, 3, 0],
      cloudInitialVelocity: [10, 0, 0],
      cloudMass: 1,
      cloudRadius: 1,
      tidalRadius: 3.63,
      impactParameter: 3,
      timeUnit: 0.978
    };
  }, [simulation]);
  const framesRef = useRef(frames);
  const frameIndexRef = useRef(currentFrameIndex);
  const isPlayingRef = useRef(false);
  const playbackSpeedRef = useRef(30);
  const animationFrameIdRef = useRef(null);
  useEffect(() => {
    framesRef.current = frames;
  }, [frames]);
  useEffect(() => {
    if (!isPlayingRef.current) {
      frameIndexRef.current = currentFrameIndex;
    }
  }, [currentFrameIndex]);
  const globalColorStats = useMemo(() => {
    return computeGlobalColorStats(frames, colorField);
  }, [frames, colorField]);
  const currentFrame = frames.get(currentFrameIndex) || null;
  const colorMap = useMemo(() => {
    const baseMap = COLOR_MAPS[colorMapName] || COLOR_MAPS.viridis;
    let min = colorRange[0];
    let max = colorRange[1];
    if (min === 0 && max === 0) {
      [min, max] = globalColorStats;
    }
    return {
      ...baseMap,
      min,
      max,
      logScale: useLogScale
    };
  }, [colorMapName, colorRange, useLogScale, globalColorStats]);
  useEffect(() => {
    if (!frames.has(currentFrameIndex)) {
      onLoadFrame(currentFrameIndex);
    }
  }, [currentFrameIndex, frames, onLoadFrame]);
  const preloadedCount = frames.size;
  const totalFrameCount = simulation?.totalFrames || 0;
  const allFramesLoaded = preloadedCount === totalFrameCount;
  useCallback((frameIndex) => {
    return frames.has(frameIndex);
  }, [frames]);
  const handleFrameChange = useCallback((frame) => {
    const clampedFrame = Math.max(0, Math.min(frame, (simulation?.totalFrames || 1) - 1));
    frameIndexRef.current = clampedFrame;
    if (!isPlayingRef.current) {
      setCurrentFrameIndex(clampedFrame);
    }
  }, [simulation?.totalFrames]);
  const startImperativePlayback = useCallback(() => {
    if (animationFrameIdRef.current !== null) return;
    isPlayingRef.current = true;
    setIsPlaying(true);
    let lastTime = 0;
    const totalFrames = simulation?.totalFrames || 1;
    const tick = (timestamp) => {
      if (!isPlayingRef.current) return;
      if (!lastTime) lastTime = timestamp;
      const elapsed = timestamp - lastTime;
      const frameInterval = 1e3 / playbackSpeedRef.current;
      if (elapsed >= frameInterval) {
        const currentIdx = frameIndexRef.current;
        const nextFrame = currentIdx + 1;
        if (nextFrame >= totalFrames) {
          stopImperativePlayback();
          frameIndexRef.current = totalFrames - 1;
          setCurrentFrameIndex(totalFrames - 1);
          return;
        }
        if (framesRef.current.has(nextFrame)) {
          frameIndexRef.current = nextFrame;
          lastTime = timestamp;
        }
      }
      animationFrameIdRef.current = requestAnimationFrame(tick);
    };
    animationFrameIdRef.current = requestAnimationFrame(tick);
  }, [simulation?.totalFrames]);
  const stopImperativePlayback = useCallback(() => {
    isPlayingRef.current = false;
    setIsPlaying(false);
    if (animationFrameIdRef.current !== null) {
      cancelAnimationFrame(animationFrameIdRef.current);
      animationFrameIdRef.current = null;
    }
    setCurrentFrameIndex(frameIndexRef.current);
  }, []);
  const handlePlayPauseChange = useCallback((playing) => {
    if (useImperativeMode) {
      if (playing) {
        startImperativePlayback();
      } else {
        stopImperativePlayback();
      }
    } else {
      setIsPlaying(playing);
    }
  }, [useImperativeMode, startImperativePlayback, stopImperativePlayback]);
  useEffect(() => {
    return () => {
      if (animationFrameIdRef.current !== null) {
        cancelAnimationFrame(animationFrameIdRef.current);
      }
    };
  }, []);
  if (!simulation) {
    return /* @__PURE__ */ jsx("div", { className: "flex items-center justify-center h-full bg-gray-900 text-gray-400", children: /* @__PURE__ */ jsxs("div", { className: "text-center", children: [
      /* @__PURE__ */ jsx("div", { className: "text-2xl mb-2", children: "📊" }),
      /* @__PURE__ */ jsx("div", { children: "No simulation loaded" }),
      /* @__PURE__ */ jsx("div", { className: "text-sm mt-1", children: "Select a simulation from the sidebar" })
    ] }) });
  }
  if (error) {
    return /* @__PURE__ */ jsx("div", { className: "flex items-center justify-center h-full bg-gray-900 text-red-400", children: /* @__PURE__ */ jsxs("div", { className: "text-center", children: [
      /* @__PURE__ */ jsx("div", { className: "text-2xl mb-2", children: "❌" }),
      /* @__PURE__ */ jsx("div", { children: "Error loading simulation" }),
      /* @__PURE__ */ jsx("div", { className: "text-sm mt-1", children: error })
    ] }) });
  }
  const currentTime = currentFrame?.time || 0;
  return /* @__PURE__ */ jsxs("div", { className: "flex flex-col h-full bg-gray-900", children: [
    /* @__PURE__ */ jsxs("div", { className: "flex items-center justify-between px-4 py-2 bg-gray-800 border-b border-gray-700", children: [
      /* @__PURE__ */ jsxs("div", { children: [
        /* @__PURE__ */ jsx("h1", { className: "text-lg font-semibold text-white", children: simulation.name }),
        /* @__PURE__ */ jsxs("div", { className: "text-xs text-gray-400", children: [
          simulation.method,
          " • ",
          simulation.kernel,
          " • ",
          simulation.particleCount.toLocaleString(),
          " particles"
        ] })
      ] }),
      /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-3", children: [
        /* @__PURE__ */ jsx("div", { className: "text-xs", children: /* @__PURE__ */ jsx("span", { className: allFramesLoaded ? "text-green-400" : "text-yellow-400", children: allFramesLoaded ? "✓ All frames in memory" : `${preloadedCount}/${totalFrameCount} frames` }) }),
        /* @__PURE__ */ jsx(
          "button",
          {
            onClick: () => setUseImperativeMode(!useImperativeMode),
            className: `px-2 py-1 text-xs rounded ${useImperativeMode ? "bg-green-600 text-white" : "bg-gray-700 text-gray-400"}`,
            title: useImperativeMode ? "High-performance mode (120+ FPS)" : "Standard React mode",
            children: useImperativeMode ? "🚀 Fast" : "⚛️ React"
          }
        ),
        /* @__PURE__ */ jsx(
          "button",
          {
            onClick: () => setShowProjections(!showProjections),
            className: `px-2 py-1 text-xs rounded ${showProjections ? "bg-blue-600 text-white" : "bg-gray-700 text-gray-400"}`,
            children: "2D Views"
          }
        ),
        /* @__PURE__ */ jsx(
          "button",
          {
            onClick: () => setShowOrbitalGeometry(!showOrbitalGeometry),
            className: `px-2 py-1 text-xs rounded ${showOrbitalGeometry ? "bg-purple-600 text-white" : "bg-gray-700 text-gray-400"}`,
            title: "Show orbital geometry diagram with impact parameter, pericentre, inclination, etc.",
            children: "🌍 Orbital"
          }
        ),
        /* @__PURE__ */ jsx(
          "button",
          {
            onClick: () => setShowCharts(!showCharts),
            className: `px-2 py-1 text-xs rounded ${showCharts ? "bg-blue-600 text-white" : "bg-gray-700 text-gray-400"}`,
            children: "Charts"
          }
        )
      ] })
    ] }),
    /* @__PURE__ */ jsxs("div", { className: "flex-1 flex overflow-hidden", children: [
      /* @__PURE__ */ jsxs("div", { className: "w-64 shrink-0 overflow-y-auto border-r border-gray-700 p-2", children: [
        /* @__PURE__ */ jsx(
          VisualizationSettings,
          {
            colorField,
            onColorFieldChange: setColorField,
            colorMapName,
            onColorMapChange: setColorMapName,
            pointSize,
            onPointSizeChange: setPointSize,
            opacity,
            onOpacityChange: setOpacity,
            showAxes,
            onShowAxesChange: setShowAxes,
            showBoundingBox,
            onShowBoundingBoxChange: setShowBoundingBox,
            colorRange,
            onColorRangeChange: setColorRange,
            useLogScale,
            onLogScaleChange: setUseLogScale
          }
        ),
        useImperativeMode && /* @__PURE__ */ jsxs("div", { className: "mt-2 bg-gray-800 rounded p-3", children: [
          /* @__PURE__ */ jsx("h3", { className: "text-sm font-medium text-gray-300 mb-2", children: "IMBH Physics" }),
          /* @__PURE__ */ jsxs("div", { className: "space-y-2 text-xs", children: [
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showBlackHole,
                  onChange: (e) => setShowBlackHole(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Black Hole" })
            ] }),
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showTrajectory,
                  onChange: (e) => setShowTrajectory(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Trajectory" })
            ] }),
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showRadii,
                  onChange: (e) => setShowRadii(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Tidal/Hill Radii" })
            ] }),
            /* @__PURE__ */ jsx("div", { className: "mt-2 pt-2 border-t border-gray-700", children: /* @__PURE__ */ jsxs("div", { className: "text-gray-400 text-[10px]", children: [
              /* @__PURE__ */ jsxs("div", { children: [
                "BH: ",
                imbhPhysics.bhMass,
                "×10⁵ M☉"
              ] }),
              /* @__PURE__ */ jsxs("div", { children: [
                "Tidal: ",
                imbhPhysics.tidalRadius.toFixed(2),
                " pc"
              ] }),
              /* @__PURE__ */ jsxs("div", { children: [
                "Impact: ",
                imbhPhysics.impactParameter.toFixed(2),
                " pc"
              ] })
            ] }) })
          ] })
        ] }),
        currentFrame && /* @__PURE__ */ jsxs("div", { className: "mt-2 bg-gray-800 rounded p-3", children: [
          /* @__PURE__ */ jsx("h3", { className: "text-sm font-medium text-gray-300 mb-2", children: "Frame Statistics" }),
          /* @__PURE__ */ jsxs("div", { className: "space-y-1 text-xs text-gray-400", children: [
            /* @__PURE__ */ jsxs("div", { className: "flex justify-between", children: [
              /* @__PURE__ */ jsx("span", { children: "Particles:" }),
              /* @__PURE__ */ jsx("span", { className: "text-white", children: currentFrame.particleCount.toLocaleString() })
            ] }),
            /* @__PURE__ */ jsxs("div", { className: "flex justify-between", children: [
              /* @__PURE__ */ jsx("span", { children: "Time:" }),
              /* @__PURE__ */ jsxs("span", { className: "text-white font-mono", children: [
                (currentTime * imbhPhysics.timeUnit).toFixed(3),
                " Myr",
                /* @__PURE__ */ jsxs("span", { className: "text-gray-500 text-[10px] ml-1", children: [
                  "(",
                  currentTime.toFixed(3),
                  " code)"
                ] })
              ] })
            ] }),
            statistics[currentFrameIndex] && /* @__PURE__ */ jsxs(Fragment, { children: [
              /* @__PURE__ */ jsxs("div", { className: "flex justify-between", children: [
                /* @__PURE__ */ jsx("span", { children: "Total Energy:" }),
                /* @__PURE__ */ jsx("span", { className: "text-white font-mono", children: statistics[currentFrameIndex].totalEnergy.toExponential(3) })
              ] }),
              /* @__PURE__ */ jsxs("div", { className: "flex justify-between", children: [
                /* @__PURE__ */ jsx("span", { children: "Max Density:" }),
                /* @__PURE__ */ jsx("span", { className: "text-white font-mono", children: statistics[currentFrameIndex].densityRange[1].toExponential(3) })
              ] })
            ] })
          ] })
        ] })
      ] }),
      /* @__PURE__ */ jsxs("div", { className: "flex-1 flex flex-col overflow-hidden", children: [
        /* @__PURE__ */ jsxs("div", { className: "flex-1 flex overflow-hidden", children: [
          /* @__PURE__ */ jsxs("div", { className: `flex-1 ${showProjections ? "" : "w-full"} relative`, children: [
            isLoading && !currentFrame ? /* @__PURE__ */ jsx("div", { className: "flex items-center justify-center h-full bg-gray-900 text-gray-400", children: /* @__PURE__ */ jsxs("div", { className: "text-center", children: [
              /* @__PURE__ */ jsx("div", { className: "animate-spin text-2xl mb-2", children: "⏳" }),
              /* @__PURE__ */ jsxs("div", { children: [
                "Loading frame ",
                currentFrameIndex,
                "..."
              ] })
            ] }) }) : useImperativeMode ? /* @__PURE__ */ jsx(
              ParticleViewer3DImperative,
              {
                framesRef,
                frameIndexRef,
                colorField,
                colorMapName,
                pointSize: pointSize * 100,
                opacity,
                logScale: useLogScale,
                showAxes,
                showBoundingBox,
                boundingBox: simulation.boundingBox,
                className: "h-full",
                onFpsUpdate: setCurrentFps,
                globalColorRange: globalColorStats,
                imbhPhysics,
                showBlackHole,
                showTrajectory,
                showRadii
              }
            ) : /* @__PURE__ */ jsx(
              ParticleViewer3D,
              {
                frame: currentFrame,
                colorField,
                colorMap,
                pointSize,
                opacity,
                showAxes,
                showBoundingBox,
                boundingBox: simulation.boundingBox,
                className: "h-full"
              }
            ),
            useImperativeMode && /* @__PURE__ */ jsxs("div", { className: "absolute bottom-2 right-2 text-green-400 text-xs font-mono bg-black/50 px-2 py-1 rounded", children: [
              "FPS: ",
              currentFps
            ] })
          ] }),
          showProjections && /* @__PURE__ */ jsxs("div", { className: "w-80 shrink-0 flex flex-col gap-2 p-2 overflow-y-auto border-l border-gray-700", children: [
            /* @__PURE__ */ jsx(
              Projection2D,
              {
                frame: currentFrame,
                projection: "xy",
                colorField,
                colorMap,
                width: 300,
                height: 200
              }
            ),
            /* @__PURE__ */ jsx(
              Projection2D,
              {
                frame: currentFrame,
                projection: "xz",
                colorField,
                colorMap,
                width: 300,
                height: 200
              }
            ),
            /* @__PURE__ */ jsx(
              Projection2D,
              {
                frame: currentFrame,
                projection: "yz",
                colorField,
                colorMap,
                width: 300,
                height: 200
              }
            )
          ] }),
          showOrbitalGeometry && /* @__PURE__ */ jsxs("div", { className: "w-[620px] shrink-0 flex flex-col p-2 overflow-y-auto border-l border-gray-700", children: [
            /* @__PURE__ */ jsxs("div", { className: "text-sm font-medium text-gray-300 mb-2 flex items-center gap-2", children: [
              /* @__PURE__ */ jsx("span", { children: "🌍" }),
              /* @__PURE__ */ jsx("span", { children: "Orbital Geometry (Oka et al. 2017 / CAT_OKA)" })
            ] }),
            /* @__PURE__ */ jsx(
              OrbitalGeometryPanel,
              {
                config: imbhPhysics.enabled ? {
                  bhMass: imbhPhysics.bhMass / 100,
                  // Convert from code units
                  cloudRadius: imbhPhysics.cloudRadius,
                  cloudMass: imbhPhysics.cloudMass,
                  impactParameter: imbhPhysics.impactParameter || 5.16,
                  pericentre: 1.69,
                  // CAT_OKA pericentre
                  eccentricity: 1.24,
                  // CAT_OKA eccentricity
                  initialPosition: imbhPhysics.cloudInitialPosition,
                  initialVelocity: imbhPhysics.cloudInitialVelocity,
                  inclination: 70,
                  // Oka et al. inclination
                  positionAngle: 41.6,
                  // Oka et al. position angle
                  lsrVelocity: -120,
                  // Oka et al. LSR velocity
                  distanceToGC: 8
                  // Distance to GC in kpc
                } : void 0,
                width: 600,
                height: 500,
                showLabels: true,
                viewMode: "3d"
              }
            ),
            /* @__PURE__ */ jsxs("div", { className: "mt-2 text-xs text-gray-400 bg-gray-800 p-2 rounded", children: [
              /* @__PURE__ */ jsx("div", { className: "font-medium text-gray-300 mb-1", children: "Key Elements:" }),
              /* @__PURE__ */ jsxs("ul", { className: "space-y-0.5", children: [
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-green-400", children: "●" }),
                  " Hyperbolic orbit trajectory"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-orange-400", children: "●" }),
                  " Impact parameter b (perpendicular to asymptote)"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-yellow-400", children: "●" }),
                  " Pericentre r_p (closest approach)"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-purple-400", children: "●" }),
                  " Line of Sight (Earth → system)"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-amber-400", children: "●" }),
                  " Inclination i = 70° (orbital plane tilt)"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-teal-400", children: "●" }),
                  " Position Angle PA = 41.6°"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-blue-400", children: "●" }),
                  " Cloud initial position & velocity"
                ] }),
                /* @__PURE__ */ jsxs("li", { children: [
                  /* @__PURE__ */ jsx("span", { className: "text-pink-400", children: "●" }),
                  " V_LSR = -120 km/s (bulk motion)"
                ] })
              ] })
            ] })
          ] })
        ] }),
        showCharts && statistics.length > 0 && /* @__PURE__ */ jsxs("div", { className: "h-52 shrink-0 border-t border-gray-700 flex gap-2 p-2 overflow-x-auto", children: [
          /* @__PURE__ */ jsx(
            EnergyChart,
            {
              statistics,
              currentFrame: currentFrameIndex,
              className: "flex-1 min-w-[300px]"
            }
          ),
          /* @__PURE__ */ jsx(
            MomentumChart,
            {
              statistics,
              currentFrame: currentFrameIndex,
              className: "flex-1 min-w-[300px]"
            }
          )
        ] })
      ] })
    ] }),
    /* @__PURE__ */ jsx("div", { className: "shrink-0 border-t border-gray-700", children: /* @__PURE__ */ jsx(
      PlaybackControls,
      {
        currentFrame: currentFrameIndex,
        totalFrames: simulation.totalFrames,
        time: currentTime,
        onFrameChange: handleFrameChange,
        onPlayPauseChange: handlePlayPauseChange,
        isFrameReady: (frame) => frames.has(frame),
        imperativeMode: useImperativeMode,
        frameIndexRef,
        playbackSpeedRef
      }
    ) })
  ] });
}
const Route$e = createFileRoute("/viz/")({
  component: VisualizationPage
});
function parseBinaryFrame(buffer, frameIndex, time, particleCount, stride, fieldOffsets) {
  const floats = new Float32Array(buffer);
  const positions = new Float32Array(particleCount * 3);
  const velocities = new Float32Array(particleCount * 3);
  const mass = new Float32Array(particleCount);
  const density = new Float32Array(particleCount);
  const pressure = new Float32Array(particleCount);
  const energy = new Float32Array(particleCount);
  const smoothingLength = new Float32Array(particleCount);
  for (let i = 0; i < particleCount; i++) {
    const base = i * stride;
    positions[i * 3] = floats[base + (fieldOffsets.x ?? 0)];
    positions[i * 3 + 1] = floats[base + (fieldOffsets.y ?? 1)];
    positions[i * 3 + 2] = floats[base + (fieldOffsets.z ?? 2)];
    velocities[i * 3] = floats[base + (fieldOffsets.vx ?? 3)];
    velocities[i * 3 + 1] = floats[base + (fieldOffsets.vy ?? 4)];
    velocities[i * 3 + 2] = floats[base + (fieldOffsets.vz ?? 5)];
    mass[i] = floats[base + (fieldOffsets.mass ?? 6)];
    density[i] = floats[base + (fieldOffsets.density ?? 7)];
    pressure[i] = floats[base + (fieldOffsets.pressure ?? 8)];
    energy[i] = floats[base + (fieldOffsets.energy ?? 9)];
    smoothingLength[i] = floats[base + (fieldOffsets.smoothing_length ?? 10)];
  }
  return {
    frameIndex,
    time,
    particleCount,
    positions,
    velocities,
    mass,
    density,
    pressure,
    energy,
    smoothingLength
  };
}
function VisualizationPage() {
  const [simulations, setSimulations] = useState([]);
  const [selectedSimulation, setSelectedSimulation] = useState(null);
  const [frames, setFrames] = useState(/* @__PURE__ */ new Map());
  const [statistics, setStatistics] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [loadingProgress, setLoadingProgress] = useState({ loaded: 0, total: 0 });
  const [error, setError] = useState(null);
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const loadSimulations = useCallback(async () => {
    try {
      const response = await fetch("/api/simulations");
      const data = await response.json();
      setSimulations(data.simulations || []);
    } catch (err) {
      console.error("Failed to load simulations:", err);
      setError("Failed to load simulations");
    }
  }, []);
  useEffect(() => {
    loadSimulations();
  }, [loadSimulations]);
  const loadSingleFrame = useCallback(async (simId, frameIndex, retries = 2) => {
    for (let attempt = 0; attempt <= retries; attempt++) {
      try {
        const response = await fetch(`/api/simulations/${simId}/frames/${frameIndex}?format=binary`);
        if (!response.ok) {
          console.error(`Failed to load frame ${frameIndex}: HTTP ${response.status}`);
          if (attempt < retries) {
            await new Promise((r) => setTimeout(r, 100 * (attempt + 1)));
            continue;
          }
          return null;
        }
        const frameIdx = parseInt(response.headers.get("X-Frame-Index") || "0");
        const time = parseFloat(response.headers.get("X-Frame-Time") || "0");
        const particleCount = parseInt(response.headers.get("X-Particle-Count") || "0");
        const stride = parseInt(response.headers.get("X-Stride") || "11");
        const fieldOffsetsHeader = response.headers.get("X-Field-Offsets");
        const fieldOffsets = fieldOffsetsHeader ? JSON.parse(fieldOffsetsHeader) : {};
        const buffer = await response.arrayBuffer();
        const expectedBytes = particleCount * stride * 4;
        if (buffer.byteLength === 0) {
          console.error(`Frame ${frameIndex}: empty buffer received`);
          if (attempt < retries) {
            await new Promise((r) => setTimeout(r, 100 * (attempt + 1)));
            continue;
          }
          return null;
        }
        if (buffer.byteLength % 4 !== 0) {
          console.error(`Frame ${frameIndex}: buffer size ${buffer.byteLength} not divisible by 4 (attempt ${attempt + 1})`);
          if (attempt < retries) {
            await new Promise((r) => setTimeout(r, 100 * (attempt + 1)));
            continue;
          }
          return null;
        }
        if (particleCount > 0 && buffer.byteLength !== expectedBytes) {
          console.warn(`Frame ${frameIndex}: buffer size mismatch. Expected ${expectedBytes}, got ${buffer.byteLength}`);
        }
        const frame = parseBinaryFrame(buffer, frameIdx, time, particleCount, stride, fieldOffsets);
        const stats = computeFrameStatistics(frame);
        return { frame, stats };
      } catch (err) {
        console.error(`Failed to load frame ${frameIndex} (attempt ${attempt + 1}):`, err);
        if (attempt < retries) {
          await new Promise((r) => setTimeout(r, 100 * (attempt + 1)));
          continue;
        }
        return null;
      }
    }
    return null;
  }, []);
  const selectSimulation = useCallback(async (sim) => {
    setSelectedSimulation(sim);
    setFrames(/* @__PURE__ */ new Map());
    setStatistics([]);
    setError(null);
    setIsLoading(true);
    setLoadingProgress({ loaded: 0, total: sim.totalFrames });
    const simId = encodeURIComponent(sim.id.replace(/\//g, "|"));
    const totalFrames = sim.totalFrames;
    const newFrames = /* @__PURE__ */ new Map();
    const newStats = [];
    console.log(`🚀 Loading all ${totalFrames} frames into memory...`);
    const startTime = performance.now();
    await new Promise((resolve) => setTimeout(resolve, 100));
    const batchSize = 10;
    let loadedCount = 0;
    for (let batchStart = 0; batchStart < totalFrames; batchStart += batchSize) {
      const batchEnd = Math.min(batchStart + batchSize, totalFrames);
      const batchPromises = [];
      for (let i = batchStart; i < batchEnd; i++) {
        batchPromises.push(
          loadSingleFrame(simId, i).then((result) => ({ idx: i, result }))
        );
      }
      const results = await Promise.all(batchPromises);
      for (const { idx, result } of results) {
        if (result) {
          newFrames.set(idx, result.frame);
          newStats[idx] = result.stats;
        }
        loadedCount++;
      }
      setLoadingProgress({ loaded: loadedCount, total: totalFrames });
      console.log(`   Loaded ${loadedCount}/${totalFrames} frames`);
    }
    const elapsed = ((performance.now() - startTime) / 1e3).toFixed(2);
    console.log(`✅ All ${totalFrames} frames loaded in ${elapsed}s - ready for instant playback!`);
    setFrames(newFrames);
    setStatistics(newStats);
    setIsLoading(false);
  }, [loadSingleFrame]);
  const loadFrame = useCallback(async (frameIndex) => {
    if (!selectedSimulation) return;
    if (frames.has(frameIndex)) return;
    if (isLoading) return;
    try {
      const simId = encodeURIComponent(selectedSimulation.id.replace(/\//g, "|"));
      const response = await fetch(`/api/simulations/${simId}/frames/${frameIndex}?format=binary`);
      if (!response.ok) {
        console.error(`Failed to load frame ${frameIndex}: HTTP ${response.status}`);
        return;
      }
      const frameIdx = parseInt(response.headers.get("X-Frame-Index") || "0");
      const time = parseFloat(response.headers.get("X-Frame-Time") || "0");
      const particleCount = parseInt(response.headers.get("X-Particle-Count") || "0");
      const stride = parseInt(response.headers.get("X-Stride") || "11");
      const fieldOffsets = JSON.parse(response.headers.get("X-Field-Offsets") || "{}");
      const buffer = await response.arrayBuffer();
      const parsed = parseBinaryFrame(buffer, frameIdx, time, particleCount, stride, fieldOffsets);
      setFrames((prev) => {
        const next = new Map(prev);
        next.set(frameIndex, parsed);
        return next;
      });
      const stats = computeFrameStatistics(parsed);
      setStatistics((prev) => {
        const next = [...prev];
        next[frameIndex] = stats;
        return next;
      });
    } catch (err) {
      console.error("Failed to load frame:", err);
    }
  }, [selectedSimulation, frames, isLoading]);
  return /* @__PURE__ */ jsxs("div", { className: "flex h-screen bg-gray-900", children: [
    /* @__PURE__ */ jsxs(
      "div",
      {
        className: `${sidebarCollapsed ? "w-12" : "w-64"} shrink-0 bg-gray-800 border-r border-gray-700 flex flex-col transition-all duration-200`,
        children: [
          /* @__PURE__ */ jsxs("div", { className: "p-3 border-b border-gray-700 flex items-center justify-between", children: [
            !sidebarCollapsed && /* @__PURE__ */ jsx("h2", { className: "text-sm font-semibold text-gray-300", children: "Simulations" }),
            /* @__PURE__ */ jsx(
              "button",
              {
                onClick: () => setSidebarCollapsed(!sidebarCollapsed),
                className: "p-1 hover:bg-gray-700 rounded text-gray-400",
                children: /* @__PURE__ */ jsx(
                  ChevronRight,
                  {
                    size: 16,
                    className: `transform transition-transform ${sidebarCollapsed ? "" : "rotate-180"}`
                  }
                )
              }
            )
          ] }),
          !sidebarCollapsed && /* @__PURE__ */ jsx("div", { className: "flex-1 overflow-y-auto", children: simulations.length === 0 ? /* @__PURE__ */ jsxs("div", { className: "p-4 text-center text-gray-500 text-sm", children: [
            /* @__PURE__ */ jsx("div", { className: "mb-2", children: "No simulations found" }),
            /* @__PURE__ */ jsx("div", { className: "text-xs", children: "Run the data exporter to prepare simulation data" })
          ] }) : /* @__PURE__ */ jsx("div", { className: "p-2 space-y-1", children: simulations.map((sim) => /* @__PURE__ */ jsxs(
            "button",
            {
              onClick: () => selectSimulation(sim),
              className: `w-full text-left p-2 rounded text-sm ${selectedSimulation?.id === sim.id ? "bg-blue-600 text-white" : "text-gray-300 hover:bg-gray-700"}`,
              children: [
                /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-2", children: [
                  /* @__PURE__ */ jsx(Folder, { size: 14 }),
                  /* @__PURE__ */ jsx("span", { className: "truncate", children: sim.name })
                ] }),
                /* @__PURE__ */ jsxs("div", { className: "text-xs opacity-60 ml-5", children: [
                  sim.totalFrames,
                  " frames • ",
                  sim.particleCount.toLocaleString(),
                  " particles"
                ] })
              ]
            },
            sim.id
          )) }) }),
          !sidebarCollapsed && /* @__PURE__ */ jsx("div", { className: "p-2 border-t border-gray-700", children: /* @__PURE__ */ jsxs(
            "button",
            {
              onClick: loadSimulations,
              className: "w-full flex items-center justify-center gap-2 p-2 text-sm text-gray-400 hover:bg-gray-700 rounded",
              children: [
                /* @__PURE__ */ jsx(RefreshCw, { size: 14 }),
                "Refresh"
              ]
            }
          ) })
        ]
      }
    ),
    /* @__PURE__ */ jsxs("div", { className: "flex-1 overflow-hidden relative", children: [
      isLoading && loadingProgress.total > 0 && /* @__PURE__ */ jsx("div", { className: "absolute inset-0 bg-gray-900/90 flex items-center justify-center z-50", children: /* @__PURE__ */ jsxs("div", { className: "text-center", children: [
        /* @__PURE__ */ jsx("div", { className: "text-xl text-white mb-4", children: "Loading Frames into Memory" }),
        /* @__PURE__ */ jsx("div", { className: "w-64 bg-gray-700 rounded-full h-3 mb-2", children: /* @__PURE__ */ jsx(
          "div",
          {
            className: "bg-blue-500 h-3 rounded-full transition-all duration-100",
            style: { width: `${loadingProgress.loaded / loadingProgress.total * 100}%` }
          }
        ) }),
        /* @__PURE__ */ jsxs("div", { className: "text-gray-400", children: [
          loadingProgress.loaded,
          " / ",
          loadingProgress.total,
          " frames"
        ] }),
        /* @__PURE__ */ jsxs("div", { className: "text-gray-500 text-sm mt-2", children: [
          (loadingProgress.loaded / loadingProgress.total * 100).toFixed(0),
          "% complete"
        ] })
      ] }) }),
      /* @__PURE__ */ jsx(
        Dashboard,
        {
          simulation: selectedSimulation,
          frames,
          statistics,
          onLoadFrame: loadFrame,
          isLoading,
          error
        }
      )
    ] })
  ] });
}
function computeFrameStatistics(frame) {
  let totalKinetic = 0;
  let totalInternal = 0;
  let totalMass = 0;
  const momentum = [0, 0, 0];
  const com = [0, 0, 0];
  let minDensity = Infinity;
  let maxDensity = -Infinity;
  let minPressure = Infinity;
  let maxPressure = -Infinity;
  for (let i = 0; i < frame.particleCount; i++) {
    const m = frame.mass[i];
    const vx = frame.velocities[i * 3];
    const vy = frame.velocities[i * 3 + 1];
    const vz = frame.velocities[i * 3 + 2];
    const u = frame.energy[i];
    const rho = frame.density[i];
    const P = frame.pressure[i];
    const v2 = vx * vx + vy * vy + vz * vz;
    totalKinetic += 0.5 * m * v2;
    totalInternal += m * u;
    momentum[0] += m * vx;
    momentum[1] += m * vy;
    momentum[2] += m * vz;
    com[0] += m * frame.positions[i * 3];
    com[1] += m * frame.positions[i * 3 + 1];
    com[2] += m * frame.positions[i * 3 + 2];
    totalMass += m;
    if (isFinite(rho)) {
      if (rho < minDensity) minDensity = rho;
      if (rho > maxDensity) maxDensity = rho;
    }
    if (isFinite(P)) {
      if (P < minPressure) minPressure = P;
      if (P > maxPressure) maxPressure = P;
    }
  }
  if (totalMass > 0) {
    com[0] /= totalMass;
    com[1] /= totalMass;
    com[2] /= totalMass;
  }
  return {
    frameIndex: frame.frameIndex,
    time: frame.time,
    totalMass,
    totalKineticEnergy: totalKinetic,
    totalInternalEnergy: totalInternal,
    totalEnergy: totalKinetic + totalInternal,
    momentum,
    centerOfMass: com,
    densityRange: [minDensity, maxDensity],
    pressureRange: [minPressure, maxPressure]
  };
}
const $$splitComponentImporter$7 = () => import("./users.index-Bef-9o5f.js");
const Route$d = createFileRoute("/users/")({
  component: lazyRouteComponent($$splitComponentImporter$7, "component")
});
const $$splitComponentImporter$6 = () => import("./posts.index-DU8oxB5n.js");
const Route$c = createFileRoute("/posts/")({
  component: lazyRouteComponent($$splitComponentImporter$6, "component")
});
const $$splitNotFoundComponentImporter$1 = () => import("./users._userId-aN2ESlS5.js");
const $$splitComponentImporter$5 = () => import("./users._userId-mOgQ8odf.js");
const $$splitErrorComponentImporter$2 = () => import("./users._userId-CG2IqJzb.js");
const Route$b = createFileRoute("/users/$userId")({
  loader: async ({
    params: {
      userId
    }
  }) => {
    try {
      const res = await fetch("/api/users/" + userId);
      if (!res.ok) {
        throw new Error("Unexpected status code");
      }
      const data = await res.json();
      return data;
    } catch {
      throw new Error("Failed to fetch user");
    }
  },
  errorComponent: lazyRouteComponent($$splitErrorComponentImporter$2, "errorComponent"),
  component: lazyRouteComponent($$splitComponentImporter$5, "component"),
  notFoundComponent: lazyRouteComponent($$splitNotFoundComponentImporter$1, "notFoundComponent")
});
const $$splitNotFoundComponentImporter = () => import("./posts._postId-DKFHFrut.js");
const $$splitComponentImporter$4 = () => import("./posts._postId-CqQCk-x_.js");
const $$splitErrorComponentImporter$1 = () => import("./posts._postId-C9z5TBp-.js");
const Route$a = createFileRoute("/posts/$postId")({
  loader: ({
    params: {
      postId
    }
  }) => fetchPost({
    data: postId
  }),
  errorComponent: lazyRouteComponent($$splitErrorComponentImporter$1, "errorComponent"),
  component: lazyRouteComponent($$splitComponentImporter$4, "component"),
  notFoundComponent: lazyRouteComponent($$splitNotFoundComponentImporter, "notFoundComponent")
});
const userLoggerMiddleware = createMiddleware().server(async ({
  next
}) => {
  console.info("In: /users");
  console.info("Request Headers:", getRequestHeaders());
  const result = await next();
  result.response.headers.set("x-users", "true");
  console.info("Out: /users");
  return result;
});
const testParentMiddleware = createMiddleware().server(async ({
  next
}) => {
  console.info("In: testParentMiddleware");
  const result = await next();
  result.response.headers.set("x-test-parent", "true");
  console.info("Out: testParentMiddleware");
  return result;
});
const testMiddleware = createMiddleware().middleware([testParentMiddleware]).server(async ({
  next
}) => {
  console.info("In: testMiddleware");
  const result = await next();
  result.response.headers.set("x-test", "true");
  console.info("Out: testMiddleware");
  return result;
});
const Route$9 = createFileRoute("/api/users")({
  server: {
    middleware: [testMiddleware, userLoggerMiddleware],
    handlers: {
      GET: async ({
        request
      }) => {
        console.info("GET /api/users @", request.url);
        console.info("Fetching users... @", request.url);
        const res = await fetch("https://jsonplaceholder.typicode.com/users");
        if (!res.ok) {
          throw new Error("Failed to fetch users");
        }
        const data = await res.json();
        const list = data.slice(0, 10);
        return json(list.map((u) => ({
          id: u.id,
          name: u.name,
          email: u.email
        })));
      }
    }
  }
});
const __filename$4 = fileURLToPath(import.meta.url);
const __dirname$4 = path.dirname(__filename$4);
const getDataRoot$3 = () => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname$4, "../../../../../");
};
function loadIMBHPhysicsFromConfig(dirPath, simulationId, dataRoot) {
  const pathParts = simulationId.split("/");
  if (pathParts.length < 3) return void 0;
  const [testCase, scenario, method] = pathParts;
  const configPaths = [
    // Direct path: sample/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json
    path.join(dataRoot, "sample", testCase, "config", "presets", "simulation", "scenarios", scenario, `${method}.json`),
    // Alternative: alongside results
    path.join(dirPath, "..", "..", "config", `${method}.json`),
    // Root config
    path.join(dataRoot, "config.json")
  ];
  for (const configPath of configPaths) {
    if (fs.existsSync(configPath)) {
      try {
        const config = JSON.parse(fs.readFileSync(configPath, "utf-8"));
        const imbhParams = config.imbh_parameters || config.externalForces?.pointMassBH;
        const initialCondition = config.initialCondition?.transform;
        const physicsSummary = config.physics_summary;
        if (imbhParams || initialCondition) {
          const bhMassCode = imbhParams?.M_BH ? imbhParams.M_BH / 1e3 : config.externalForces?.pointMassBH?.mass || 100;
          const cloudPos = imbhParams?.cloud_initial_position || initialCondition?.translate || [-20, 3, 0];
          const cloudVel = imbhParams?.cloud_initial_velocity || initialCondition?.velocity_boost || [10, 0, 0];
          const bhPos = imbhParams?.BH_initial_position || config.externalForces?.pointMassBH?.position || [0, 0, 0];
          const impactParameter = Math.abs(cloudPos[1] - bhPos[1]);
          const tidalRadius = physicsSummary?.tidal_radius_pc || 3.63;
          const cloudRadius = physicsSummary?.cloud_radius_pc || 1.13;
          return {
            enabled: true,
            bhPosition: bhPos,
            bhMass: bhMassCode,
            cloudInitialPosition: cloudPos,
            cloudInitialVelocity: cloudVel,
            cloudMass: (physicsSummary?.cloud_mass_Msun || 1e3) / 1e3,
            // To code units
            cloudRadius,
            tidalRadius,
            impactParameter,
            timeUnit: 0.978
            // From unit system
          };
        }
      } catch (e) {
        console.warn(`Failed to parse config at ${configPath}:`, e);
      }
    }
  }
  return void 0;
}
async function scanSimulationDirectory(dirPath, simulationId, dataRoot = getDataRoot$3()) {
  try {
    const metadataPath = path.join(dirPath, "viz_data", "metadata.json");
    if (fs.existsSync(metadataPath)) {
      const metadata = JSON.parse(fs.readFileSync(metadataPath, "utf-8"));
      const imbhPhysics = loadIMBHPhysicsFromConfig(dirPath, simulationId, dataRoot);
      return {
        ...metadata,
        id: simulationId,
        dataPath: path.join(dirPath, "viz_data"),
        imbhPhysics
      };
    }
    const vizDataPath = path.join(dirPath, "viz_data");
    if (!fs.existsSync(vizDataPath)) {
      return null;
    }
    const files = fs.readdirSync(vizDataPath);
    const frameFiles = files.filter((f) => f.startsWith("frame_") && f.endsWith(".bin"));
    if (frameFiles.length === 0) {
      return null;
    }
    const frameNumbers = frameFiles.map((f) => parseInt(f.replace("frame_", "").replace(".bin", ""))).sort((a, b) => a - b);
    return {
      id: simulationId,
      name: simulationId.replace(/_/g, " "),
      description: `Simulation data from ${simulationId}`,
      method: "Unknown",
      kernel: "Unknown",
      dimensions: 3,
      totalFrames: frameFiles.length,
      particleCount: 0,
      // Will be determined when loading frame
      timeRange: [0, frameFiles.length],
      boundingBox: {
        min: [-1, -1, -1],
        max: [1, 1, 1]
      },
      dataPath: vizDataPath,
      createdAt: (/* @__PURE__ */ new Date()).toISOString()
    };
  } catch (error) {
    console.error(`Error scanning ${dirPath}:`, error);
    return null;
  }
}
async function findSimulations() {
  const dataRoot = getDataRoot$3();
  const simulations = [];
  console.log("🔍 Starting simulation scan...");
  console.log(`   Data root: ${dataRoot}`);
  const startTime = Date.now();
  const sampleDir = path.join(dataRoot, "sample");
  if (fs.existsSync(sampleDir)) {
    console.log("📁 Scanning sample/ directory...");
    const entries = fs.readdirSync(sampleDir, { withFileTypes: true });
    const totalEntries = entries.filter((e) => e.isDirectory()).length;
    let scanned = 0;
    for (const entry of entries) {
      if (entry.isDirectory()) {
        scanned++;
        process.stdout.write(`\r   [${scanned}/${totalEntries}] Scanning ${entry.name}...`.padEnd(60));
        const resultsDir = path.join(sampleDir, entry.name, "results");
        if (fs.existsSync(resultsDir)) {
          const resultEntries = fs.readdirSync(resultsDir, { withFileTypes: true });
          for (const resultEntry of resultEntries) {
            if (resultEntry.isDirectory()) {
              const sim2 = await scanSimulationDirectory(
                path.join(resultsDir, resultEntry.name),
                `${entry.name}/${resultEntry.name}`,
                dataRoot
              );
              if (sim2) {
                simulations.push(sim2);
                console.log(`
   ✓ Found: ${entry.name}/${resultEntry.name} (${sim2.totalFrames} frames)`);
              }
              const nestedDir = path.join(resultsDir, resultEntry.name);
              const nestedEntries = fs.readdirSync(nestedDir, { withFileTypes: true });
              for (const nestedEntry of nestedEntries) {
                if (nestedEntry.isDirectory()) {
                  const nestedSim = await scanSimulationDirectory(
                    path.join(nestedDir, nestedEntry.name),
                    `${entry.name}/${resultEntry.name}/${nestedEntry.name}`,
                    dataRoot
                  );
                  if (nestedSim) {
                    simulations.push(nestedSim);
                    console.log(`
   ✓ Found: ${entry.name}/${resultEntry.name}/${nestedEntry.name} (${nestedSim.totalFrames} frames)`);
                  }
                }
              }
            }
          }
        }
        const sim = await scanSimulationDirectory(
          path.join(sampleDir, entry.name),
          entry.name,
          dataRoot
        );
        if (sim) {
          simulations.push(sim);
          console.log(`
   ✓ Found: ${entry.name} (${sim.totalFrames} frames)`);
        }
      }
    }
    console.log("");
  }
  const laneEmdenDir = path.join(dataRoot, "lane_emden", "results");
  if (fs.existsSync(laneEmdenDir)) {
    console.log("📁 Scanning lane_emden/results/ directory...");
    const entries = fs.readdirSync(laneEmdenDir, { withFileTypes: true });
    for (const entry of entries) {
      if (entry.isDirectory()) {
        const sim = await scanSimulationDirectory(
          path.join(laneEmdenDir, entry.name),
          `lane_emden/${entry.name}`,
          dataRoot
        );
        if (sim) {
          simulations.push(sim);
          console.log(`   ✓ Found: lane_emden/${entry.name} (${sim.totalFrames} frames)`);
        }
      }
    }
  }
  const elapsed = ((Date.now() - startTime) / 1e3).toFixed(2);
  console.log(`
✅ Scan complete: Found ${simulations.length} simulations in ${elapsed}s`);
  return simulations;
}
const Route$8 = createFileRoute("/api/simulations")({
  server: {
    handlers: {
      GET: async () => {
        try {
          const simulations = await findSimulations();
          const response = { simulations };
          return json(response);
        } catch (error) {
          console.error("Error listing simulations:", error);
          return json({ error: "Failed to list simulations", simulations: [] }, { status: 500 });
        }
      }
    }
  }
});
const $$splitComponentImporter$3 = () => import("./_nested-layout-BocDAsiI.js");
const Route$7 = createFileRoute("/_pathlessLayout/_nested-layout")({
  component: lazyRouteComponent($$splitComponentImporter$3, "component")
});
const $$splitComponentImporter$2 = () => import("./posts_._postId.deep-BLVjA6x7.js");
const $$splitErrorComponentImporter = () => import("./posts_._postId.deep-C9z5TBp-.js");
const Route$6 = createFileRoute("/posts_/$postId/deep")({
  loader: async ({
    params: {
      postId
    }
  }) => fetchPost({
    data: postId
  }),
  errorComponent: lazyRouteComponent($$splitErrorComponentImporter, "errorComponent"),
  component: lazyRouteComponent($$splitComponentImporter$2, "component")
});
const Route$5 = createFileRoute("/api/users/$userId")({
  server: {
    handlers: {
      GET: async ({ params, request }) => {
        console.info(`Fetching users by id=${params.userId}... @`, request.url);
        try {
          const res = await fetch(
            "https://jsonplaceholder.typicode.com/users/" + params.userId
          );
          if (!res.ok) {
            throw new Error("Failed to fetch user");
          }
          const user = await res.json();
          return json({
            id: user.id,
            name: user.name,
            email: user.email
          });
        } catch (e) {
          console.error(e);
          return json({ error: "User not found" }, { status: 404 });
        }
      }
    }
  }
});
const __filename$3 = fileURLToPath(import.meta.url);
const __dirname$3 = path.dirname(__filename$3);
const getDataRoot$2 = () => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname$3, "../../../../");
};
const Route$4 = createFileRoute("/api/simulations/$simId")({
  server: {
    handlers: {
      GET: async ({ params }) => {
        try {
          const { simId } = params;
          const dataRoot = getDataRoot$2();
          const simPath = decodeURIComponent(simId).replace(/\|/g, "/");
          const pathParts = simPath.split("/");
          const possiblePaths = [
            path.join(dataRoot, "sample", simPath, "viz_data"),
            path.join(dataRoot, "sample", simPath, "results", "viz_data"),
            path.join(dataRoot, "lane_emden", "results", simPath, "viz_data"),
            path.join(dataRoot, simPath, "viz_data")
          ];
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts;
            possiblePaths.unshift(
              path.join(dataRoot, "sample", testName, "results", ...rest, "viz_data")
            );
          }
          let dataPath = null;
          for (const p of possiblePaths) {
            if (fs.existsSync(p)) {
              dataPath = p;
              break;
            }
          }
          if (!dataPath) {
            return json({ error: `Simulation not found: ${simId}` }, { status: 404 });
          }
          const metadataPath = path.join(dataPath, "metadata.json");
          if (fs.existsSync(metadataPath)) {
            const metadata2 = JSON.parse(fs.readFileSync(metadataPath, "utf-8"));
            return json({
              ...metadata2,
              id: simId,
              dataPath
            });
          }
          const files = fs.readdirSync(dataPath);
          const frameFiles = files.filter((f) => f.startsWith("frame_") && f.endsWith(".bin"));
          const metadata = {
            id: simId,
            name: simId.replace(/[_|]/g, " "),
            description: `Simulation: ${simId}`,
            method: "Unknown",
            kernel: "Unknown",
            dimensions: 3,
            totalFrames: frameFiles.length,
            particleCount: 0,
            timeRange: [0, frameFiles.length],
            boundingBox: {
              min: [-1, -1, -1],
              max: [1, 1, 1]
            },
            dataPath,
            createdAt: (/* @__PURE__ */ new Date()).toISOString()
          };
          return json(metadata);
        } catch (error) {
          console.error("Error loading simulation metadata:", error);
          return json({ error: "Failed to load simulation metadata" }, { status: 500 });
        }
      }
    }
  }
});
const $$splitComponentImporter$1 = () => import("./route-b-CsHX6n6-.js");
const Route$3 = createFileRoute("/_pathlessLayout/_nested-layout/route-b")({
  component: lazyRouteComponent($$splitComponentImporter$1, "component")
});
const $$splitComponentImporter = () => import("./route-a-xd-e2Wm0.js");
const Route$2 = createFileRoute("/_pathlessLayout/_nested-layout/route-a")({
  component: lazyRouteComponent($$splitComponentImporter, "component")
});
const __filename$2 = fileURLToPath(import.meta.url);
const __dirname$2 = path.dirname(__filename$2);
const log = {
  info: (msg, ...args) => console.log(`[FRAME-API] ${(/* @__PURE__ */ new Date()).toISOString()} INFO: ${msg}`, ...args),
  warn: (msg, ...args) => console.warn(`[FRAME-API] ${(/* @__PURE__ */ new Date()).toISOString()} WARN: ${msg}`, ...args),
  error: (msg, ...args) => console.error(`[FRAME-API] ${(/* @__PURE__ */ new Date()).toISOString()} ERROR: ${msg}`, ...args),
  debug: (msg, ...args) => console.log(`[FRAME-API] ${(/* @__PURE__ */ new Date()).toISOString()} DEBUG: ${msg}`, ...args)
};
const getDataRoot$1 = () => {
  const root = process.env.SPH_DATA_ROOT || path.resolve(__dirname$2, "../../../../../");
  log.debug(`getDataRoot() __dirname=${__dirname$2}, resolved=${root}`);
  return root;
};
const Route$1 = createFileRoute("/api/simulations/$simId/frames/$frameId")({
  validateSearch: (search) => ({
    format: search.format || "json"
  }),
  server: {
    handlers: {
      GET: async ({ params, request }) => {
        try {
          const { simId, frameId } = params;
          const frameIndex = parseInt(frameId);
          const url = new URL(request.url);
          const formatParam = url.searchParams.get("format");
          const acceptHeader = request.headers.get("Accept") || "";
          const wantBinary = formatParam === "binary" || acceptHeader.includes("application/octet-stream");
          console.log(`📥 Frame request: simId=${simId}, frameId=${frameId}, format=${wantBinary ? "binary" : "json"}`);
          if (isNaN(frameIndex)) {
            if (wantBinary) {
              return new Response("Invalid frame ID", { status: 400 });
            }
            return json({ error: "Invalid frame ID" }, { status: 400 });
          }
          const dataRoot = getDataRoot$1();
          const simPath = decodeURIComponent(simId).replace(/\|/g, "/");
          console.log(`   Decoded simPath: ${simPath}`);
          console.log(`   Data root: ${dataRoot}`);
          const pathParts = simPath.split("/");
          const possiblePaths = [
            path.join(dataRoot, "sample", simPath, "viz_data"),
            path.join(dataRoot, "sample", simPath, "results", "viz_data"),
            path.join(dataRoot, "lane_emden", "results", simPath, "viz_data"),
            path.join(dataRoot, simPath, "viz_data")
          ];
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts;
            possiblePaths.unshift(
              path.join(dataRoot, "sample", testName, "results", ...rest, "viz_data")
            );
          }
          console.log(`   Possible paths:`);
          possiblePaths.forEach((p, i) => console.log(`     ${i}: ${p}`));
          let dataPath = null;
          for (const p of possiblePaths) {
            const exists = fs.existsSync(p);
            console.log(`   Checking: ${p} => ${exists}`);
            if (exists) {
              dataPath = p;
              break;
            }
          }
          if (!dataPath) {
            if (wantBinary) {
              return new Response(`Simulation data not found: ${simId}`, { status: 404 });
            }
            return json({ error: `Simulation data not found: ${simId}` }, { status: 404 });
          }
          const framePath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, "0")}.bin`);
          if (!fs.existsSync(framePath)) {
            if (wantBinary) {
              return new Response(`Frame ${frameIndex} not found`, { status: 404 });
            }
            return json({ error: `Frame ${frameIndex} not found` }, { status: 404 });
          }
          const buffer = fs.readFileSync(framePath);
          const metadataPath = path.join(dataPath, "metadata.json");
          let metadata = {};
          if (fs.existsSync(metadataPath)) {
            metadata = JSON.parse(fs.readFileSync(metadataPath, "utf-8"));
          }
          const frameInfoPath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, "0")}.json`);
          let frameInfo = { time: frameIndex };
          if (fs.existsSync(frameInfoPath)) {
            frameInfo = JSON.parse(fs.readFileSync(frameInfoPath, "utf-8"));
          }
          const fieldOffsets = metadata.fieldOffsets || {
            x: 0,
            y: 1,
            z: 2,
            vx: 3,
            vy: 4,
            vz: 5,
            mass: 6,
            density: 7,
            pressure: 8,
            energy: 9,
            smoothingLength: 10
          };
          const stride = metadata.stride || 11;
          const particleCount = Math.floor(buffer.byteLength / (stride * 4));
          if (wantBinary) {
            return new Response(buffer, {
              status: 200,
              headers: {
                "Content-Type": "application/octet-stream",
                "Content-Length": buffer.byteLength.toString(),
                "Cache-Control": "public, max-age=3600",
                "X-Frame-Index": frameIndex.toString(),
                "X-Frame-Time": (frameInfo.time || frameIndex).toString(),
                "X-Particle-Count": particleCount.toString(),
                "X-Stride": stride.toString(),
                "X-Field-Offsets": JSON.stringify(fieldOffsets)
              }
            });
          }
          const response = {
            frameIndex,
            time: frameInfo.time || frameIndex,
            data: buffer.toString("base64"),
            stride,
            fieldOffsets,
            particleCount
          };
          return json(response, {
            headers: {
              "Cache-Control": "public, max-age=3600"
            }
          });
        } catch (error) {
          console.error("Error loading frame:", error);
          return json({ error: "Failed to load frame data" }, { status: 500 });
        }
      }
    }
  }
});
const __filename$1 = fileURLToPath(import.meta.url);
const __dirname$1 = path.dirname(__filename$1);
const getDataRoot = () => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname$1, "../../../../../");
};
const Route = createFileRoute("/api/simulations/$simId/frames/$frameId/bin")({
  server: {
    handlers: {
      GET: async ({ params }) => {
        try {
          const { simId, frameId } = params;
          const frameIndex = parseInt(frameId);
          if (isNaN(frameIndex)) {
            return new Response("Invalid frame ID", { status: 400 });
          }
          const dataRoot = getDataRoot();
          const simPath = decodeURIComponent(simId).replace(/\|/g, "/");
          const pathParts = simPath.split("/");
          const possiblePaths = [
            path.join(dataRoot, "sample", simPath, "viz_data"),
            path.join(dataRoot, "sample", simPath, "results", "viz_data"),
            path.join(dataRoot, "lane_emden", "results", simPath, "viz_data"),
            path.join(dataRoot, simPath, "viz_data")
          ];
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts;
            possiblePaths.unshift(
              path.join(dataRoot, "sample", testName, "results", ...rest, "viz_data")
            );
          }
          let dataPath = null;
          for (const p of possiblePaths) {
            if (fs.existsSync(p)) {
              dataPath = p;
              break;
            }
          }
          if (!dataPath) {
            return new Response("Simulation not found", { status: 404 });
          }
          const framePath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, "0")}.bin`);
          if (!fs.existsSync(framePath)) {
            return new Response("Frame not found", { status: 404 });
          }
          const buffer = fs.readFileSync(framePath);
          const metadataPath = path.join(dataPath, "metadata.json");
          let metadata = {};
          if (fs.existsSync(metadataPath)) {
            metadata = JSON.parse(fs.readFileSync(metadataPath, "utf-8"));
          }
          const frameInfoPath = path.join(dataPath, `frame_${frameIndex.toString().padStart(5, "0")}.json`);
          let frameInfo = { time: frameIndex };
          if (fs.existsSync(frameInfoPath)) {
            frameInfo = JSON.parse(fs.readFileSync(frameInfoPath, "utf-8"));
          }
          const fieldOffsets = metadata.fieldOffsets || {
            x: 0,
            y: 1,
            z: 2,
            vx: 3,
            vy: 4,
            vz: 5,
            mass: 6,
            density: 7,
            pressure: 8,
            energy: 9,
            smoothing_length: 10
          };
          const stride = metadata.stride || 11;
          const particleCount = Math.floor(buffer.byteLength / (stride * 4));
          return new Response(buffer, {
            status: 200,
            headers: {
              "Content-Type": "application/octet-stream",
              "Content-Length": buffer.byteLength.toString(),
              "Cache-Control": "public, max-age=3600",
              "X-Frame-Index": frameIndex.toString(),
              "X-Frame-Time": (frameInfo.time || frameIndex).toString(),
              "X-Particle-Count": particleCount.toString(),
              "X-Stride": stride.toString(),
              "X-Field-Offsets": JSON.stringify(fieldOffsets)
            }
          });
        } catch (error) {
          console.error("Error loading binary frame:", error);
          return new Response("Server error", { status: 500 });
        }
      }
    }
  }
});
const UsersRoute = Route$l.update({
  id: "/users",
  path: "/users",
  getParentRoute: () => Route$m
});
const RedirectRoute = Route$k.update({
  id: "/redirect",
  path: "/redirect",
  getParentRoute: () => Route$m
});
const PostsRoute = Route$j.update({
  id: "/posts",
  path: "/posts",
  getParentRoute: () => Route$m
});
const DeferredRoute = Route$i.update({
  id: "/deferred",
  path: "/deferred",
  getParentRoute: () => Route$m
});
const CustomScriptDotjsRoute = Route$h.update({
  id: "/customScript.js",
  path: "/customScript.js",
  getParentRoute: () => Route$m
});
const PathlessLayoutRoute = Route$g.update({
  id: "/_pathlessLayout",
  getParentRoute: () => Route$m
});
const IndexRoute = Route$f.update({
  id: "/",
  path: "/",
  getParentRoute: () => Route$m
});
const VizIndexRoute = Route$e.update({
  id: "/viz/",
  path: "/viz/",
  getParentRoute: () => Route$m
});
const UsersIndexRoute = Route$d.update({
  id: "/",
  path: "/",
  getParentRoute: () => UsersRoute
});
const PostsIndexRoute = Route$c.update({
  id: "/",
  path: "/",
  getParentRoute: () => PostsRoute
});
const UsersUserIdRoute = Route$b.update({
  id: "/$userId",
  path: "/$userId",
  getParentRoute: () => UsersRoute
});
const PostsPostIdRoute = Route$a.update({
  id: "/$postId",
  path: "/$postId",
  getParentRoute: () => PostsRoute
});
const ApiUsersRoute = Route$9.update({
  id: "/api/users",
  path: "/api/users",
  getParentRoute: () => Route$m
});
const ApiSimulationsRoute = Route$8.update({
  id: "/api/simulations",
  path: "/api/simulations",
  getParentRoute: () => Route$m
});
const PathlessLayoutNestedLayoutRoute = Route$7.update({
  id: "/_nested-layout",
  getParentRoute: () => PathlessLayoutRoute
});
const PostsPostIdDeepRoute = Route$6.update({
  id: "/posts_/$postId/deep",
  path: "/posts/$postId/deep",
  getParentRoute: () => Route$m
});
const ApiUsersUserIdRoute = Route$5.update({
  id: "/$userId",
  path: "/$userId",
  getParentRoute: () => ApiUsersRoute
});
const ApiSimulationsSimIdRoute = Route$4.update({
  id: "/$simId",
  path: "/$simId",
  getParentRoute: () => ApiSimulationsRoute
});
const PathlessLayoutNestedLayoutRouteBRoute = Route$3.update({
  id: "/route-b",
  path: "/route-b",
  getParentRoute: () => PathlessLayoutNestedLayoutRoute
});
const PathlessLayoutNestedLayoutRouteARoute = Route$2.update({
  id: "/route-a",
  path: "/route-a",
  getParentRoute: () => PathlessLayoutNestedLayoutRoute
});
const ApiSimulationsSimIdFramesFrameIdRoute = Route$1.update({
  id: "/frames/$frameId",
  path: "/frames/$frameId",
  getParentRoute: () => ApiSimulationsSimIdRoute
});
const ApiSimulationsSimIdFramesFrameIdBinRoute = Route.update({
  id: "/bin",
  path: "/bin",
  getParentRoute: () => ApiSimulationsSimIdFramesFrameIdRoute
});
const PathlessLayoutNestedLayoutRouteChildren = {
  PathlessLayoutNestedLayoutRouteARoute,
  PathlessLayoutNestedLayoutRouteBRoute
};
const PathlessLayoutNestedLayoutRouteWithChildren = PathlessLayoutNestedLayoutRoute._addFileChildren(
  PathlessLayoutNestedLayoutRouteChildren
);
const PathlessLayoutRouteChildren = {
  PathlessLayoutNestedLayoutRoute: PathlessLayoutNestedLayoutRouteWithChildren
};
const PathlessLayoutRouteWithChildren = PathlessLayoutRoute._addFileChildren(
  PathlessLayoutRouteChildren
);
const PostsRouteChildren = {
  PostsPostIdRoute,
  PostsIndexRoute
};
const PostsRouteWithChildren = PostsRoute._addFileChildren(PostsRouteChildren);
const UsersRouteChildren = {
  UsersUserIdRoute,
  UsersIndexRoute
};
const UsersRouteWithChildren = UsersRoute._addFileChildren(UsersRouteChildren);
const ApiSimulationsSimIdFramesFrameIdRouteChildren = {
  ApiSimulationsSimIdFramesFrameIdBinRoute
};
const ApiSimulationsSimIdFramesFrameIdRouteWithChildren = ApiSimulationsSimIdFramesFrameIdRoute._addFileChildren(
  ApiSimulationsSimIdFramesFrameIdRouteChildren
);
const ApiSimulationsSimIdRouteChildren = {
  ApiSimulationsSimIdFramesFrameIdRoute: ApiSimulationsSimIdFramesFrameIdRouteWithChildren
};
const ApiSimulationsSimIdRouteWithChildren = ApiSimulationsSimIdRoute._addFileChildren(ApiSimulationsSimIdRouteChildren);
const ApiSimulationsRouteChildren = {
  ApiSimulationsSimIdRoute: ApiSimulationsSimIdRouteWithChildren
};
const ApiSimulationsRouteWithChildren = ApiSimulationsRoute._addFileChildren(
  ApiSimulationsRouteChildren
);
const ApiUsersRouteChildren = {
  ApiUsersUserIdRoute
};
const ApiUsersRouteWithChildren = ApiUsersRoute._addFileChildren(
  ApiUsersRouteChildren
);
const rootRouteChildren = {
  IndexRoute,
  PathlessLayoutRoute: PathlessLayoutRouteWithChildren,
  CustomScriptDotjsRoute,
  DeferredRoute,
  PostsRoute: PostsRouteWithChildren,
  RedirectRoute,
  UsersRoute: UsersRouteWithChildren,
  ApiSimulationsRoute: ApiSimulationsRouteWithChildren,
  ApiUsersRoute: ApiUsersRouteWithChildren,
  VizIndexRoute,
  PostsPostIdDeepRoute
};
const routeTree = Route$m._addFileChildren(rootRouteChildren)._addFileTypes();
function getRouter() {
  const router2 = createRouter({
    routeTree,
    defaultPreload: "intent",
    defaultErrorComponent: DefaultCatchBoundary,
    defaultNotFoundComponent: () => /* @__PURE__ */ jsx(NotFound, {}),
    scrollRestoration: true
  });
  return router2;
}
const router = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  getRouter
}, Symbol.toStringTag, { value: "Module" }));
export {
  NotFound as N,
  Route$l as R,
  Route$j as a,
  Route$i as b,
  Route$b as c,
  Route$a as d,
  Route$6 as e,
  router as r
};
