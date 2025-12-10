import { jsxs, jsx } from "react/jsx-runtime";
import { Link, Outlet } from "@tanstack/react-router";
function LayoutComponent() {
  return /* @__PURE__ */ jsxs("div", { children: [
    /* @__PURE__ */ jsx("div", { children: "I'm a nested layout" }),
    /* @__PURE__ */ jsxs("div", { className: "flex gap-2 border-b", children: [
      /* @__PURE__ */ jsx(Link, { to: "/route-a", activeProps: {
        className: "font-bold"
      }, children: "Go to route A" }),
      /* @__PURE__ */ jsx(Link, { to: "/route-b", activeProps: {
        className: "font-bold"
      }, children: "Go to route B" })
    ] }),
    /* @__PURE__ */ jsx("div", { children: /* @__PURE__ */ jsx(Outlet, {}) })
  ] });
}
export {
  LayoutComponent as component
};
