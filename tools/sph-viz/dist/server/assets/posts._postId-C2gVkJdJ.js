import { jsx } from "react/jsx-runtime";
import { N as NotFound } from "./router-BnacyNlB.js";
import "@tanstack/react-router";
import "@tanstack/react-router-devtools";
import "../server.js";
import "@tanstack/history";
import "@tanstack/router-core/ssr/client";
import "@tanstack/router-core";
import "node:async_hooks";
import "@tanstack/router-core/ssr/server";
import "h3-v2";
import "tiny-invariant";
import "seroval";
import "@tanstack/react-router/ssr/server";
import "react";
import "@react-three/fiber";
import "@react-three/drei";
import "three";
import "three/examples/jsm/controls/OrbitControls.js";
import "recharts";
import "lucide-react";
import "fs";
import "path";
import "url";
const SplitNotFoundComponent = () => {
  return /* @__PURE__ */ jsx(NotFound, { children: "Post not found" });
};
export {
  SplitNotFoundComponent as notFoundComponent
};
