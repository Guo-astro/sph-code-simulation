import { jsxs, jsx } from "react/jsx-runtime";
import { Link, Outlet } from "@tanstack/react-router";
import { R as Route } from "./router-B2fdo7UV.js";
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
function UsersComponent() {
  const users = Route.useLoaderData();
  return /* @__PURE__ */ jsxs("div", { className: "p-2 flex gap-2", children: [
    /* @__PURE__ */ jsx("ul", { className: "list-disc pl-4", children: [...users, {
      id: "i-do-not-exist",
      name: "Non-existent User",
      email: ""
    }].map((user) => {
      return /* @__PURE__ */ jsx("li", { className: "whitespace-nowrap", children: /* @__PURE__ */ jsx(Link, { to: "/users/$userId", params: {
        userId: String(user.id)
      }, className: "block py-1 text-blue-800 hover:text-blue-600", activeProps: {
        className: "text-black font-bold"
      }, children: /* @__PURE__ */ jsx("div", { children: user.name }) }) }, user.id);
    }) }),
    /* @__PURE__ */ jsx("hr", {}),
    /* @__PURE__ */ jsx(Outlet, {})
  ] });
}
export {
  UsersComponent as component
};
