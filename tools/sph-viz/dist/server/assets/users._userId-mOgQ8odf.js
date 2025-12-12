import { jsxs, jsx } from "react/jsx-runtime";
import { c as Route } from "./router-B2fdo7UV.js";
import "@tanstack/react-router";
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
function UserComponent() {
  const user = Route.useLoaderData();
  return /* @__PURE__ */ jsxs("div", { className: "space-y-2", children: [
    /* @__PURE__ */ jsx("h4", { className: "text-xl font-bold underline", children: user.name }),
    /* @__PURE__ */ jsx("div", { className: "text-sm", children: user.email }),
    /* @__PURE__ */ jsx("div", { children: /* @__PURE__ */ jsx("a", { href: `/api/users/${user.id}`, className: "text-blue-800 hover:text-blue-600 underline", children: "View as JSON" }) })
  ] });
}
export {
  UserComponent as component
};
