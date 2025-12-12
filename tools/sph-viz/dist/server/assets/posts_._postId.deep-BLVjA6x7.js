import { jsxs, jsx } from "react/jsx-runtime";
import { Link } from "@tanstack/react-router";
import { e as Route } from "./router-B2fdo7UV.js";
import "@tanstack/react-router-devtools";
import "../server.js";
import "node:async_hooks";
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
function PostDeepComponent() {
  const post = Route.useLoaderData();
  return /* @__PURE__ */ jsxs("div", { className: "p-2 space-y-2", children: [
    /* @__PURE__ */ jsx(Link, { to: "/posts", className: "block py-1 text-blue-800 hover:text-blue-600", children: "← All Posts" }),
    /* @__PURE__ */ jsx("h4", { className: "text-xl font-bold underline", children: post.title }),
    /* @__PURE__ */ jsx("div", { className: "text-sm", children: post.body })
  ] });
}
export {
  PostDeepComponent as component
};
