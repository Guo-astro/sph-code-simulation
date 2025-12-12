import { jsx, jsxs } from "react/jsx-runtime";
import { L as Link } from "./router-C1SSjqoG.js";
import "../server.js";
import "node:stream/web";
import "node:stream";
import "node:async_hooks";
import "react-dom/server";
import "react";
import "react-dom";
import "@react-three/fiber";
import "@react-three/drei";
import "three";
import "three/examples/jsm/controls/OrbitControls.js";
import "recharts";
import "lucide-react";
import "fs";
import "path";
import "url";
function Home() {
  return /* @__PURE__ */ jsx("div", { className: "min-h-screen bg-gradient-to-br from-gray-900 via-gray-800 to-gray-900 flex items-center justify-center", children: /* @__PURE__ */ jsxs("div", { className: "text-center max-w-2xl px-8", children: [
    /* @__PURE__ */ jsx("h1", { className: "text-5xl font-bold text-white mb-4", children: "🌊 SPH Viz" }),
    /* @__PURE__ */ jsx("p", { className: "text-xl text-gray-400 mb-8", children: "Interactive visualization dashboard for Smoothed Particle Hydrodynamics simulations" }),
    /* @__PURE__ */ jsxs("div", { className: "space-y-4", children: [
      /* @__PURE__ */ jsx(Link, { to: "/viz", className: "inline-block px-8 py-4 bg-blue-600 hover:bg-blue-500 text-white text-lg font-semibold rounded-lg shadow-lg transition-all hover:scale-105", children: "Open Visualization Dashboard" }),
      /* @__PURE__ */ jsxs("div", { className: "text-gray-500 text-sm mt-8", children: [
        /* @__PURE__ */ jsx("p", { children: "Supports multiple SPH methods: GSPH, SSPH, DISPH, GDISPH, SRGSPH" }),
        /* @__PURE__ */ jsx("p", { className: "mt-1", children: "Real-time 3D rendering with React Three Fiber" })
      ] })
    ] }),
    /* @__PURE__ */ jsxs("div", { className: "mt-16 grid grid-cols-3 gap-8 text-center", children: [
      /* @__PURE__ */ jsxs("div", { className: "p-4 bg-gray-800/50 rounded-lg", children: [
        /* @__PURE__ */ jsx("div", { className: "text-3xl mb-2", children: "🎬" }),
        /* @__PURE__ */ jsx("div", { className: "text-gray-400 text-sm", children: "Animation Playback" })
      ] }),
      /* @__PURE__ */ jsxs("div", { className: "p-4 bg-gray-800/50 rounded-lg", children: [
        /* @__PURE__ */ jsx("div", { className: "text-3xl mb-2", children: "📊" }),
        /* @__PURE__ */ jsx("div", { className: "text-gray-400 text-sm", children: "Real-time Statistics" })
      ] }),
      /* @__PURE__ */ jsxs("div", { className: "p-4 bg-gray-800/50 rounded-lg", children: [
        /* @__PURE__ */ jsx("div", { className: "text-3xl mb-2", children: "🔬" }),
        /* @__PURE__ */ jsx("div", { className: "text-gray-400 text-sm", children: "Multiple Projections" })
      ] })
    ] })
  ] }) });
}
export {
  Home as component
};
