import { jsxs, jsx } from "react/jsx-runtime";
function UsersIndexComponent() {
  return /* @__PURE__ */ jsxs("div", { children: [
    "Select a user or",
    " ",
    /* @__PURE__ */ jsx("a", { href: "/api/users", className: "text-blue-800 hover:text-blue-600 underline", children: "view as JSON" })
  ] });
}
export {
  UsersIndexComponent as component
};
