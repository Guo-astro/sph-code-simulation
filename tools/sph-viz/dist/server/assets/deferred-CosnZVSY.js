import { jsxs, jsx } from "react/jsx-runtime";
import { b as createServerRpc, c as createServerFn } from "../server.js";
import { createFileRoute, Await } from "@tanstack/react-router";
import { useState, Suspense } from "react";
import "node:async_hooks";
import "@tanstack/react-router/ssr/server";
const personServerFn_createServerFn_handler = createServerRpc("f76e8f8721c12c8547a3ced6a10916f5b5076c1a10dcbeaa607360ce419d0a48", (opts, signal) => {
  return personServerFn.__executeServer(opts, signal);
});
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
const slowServerFn_createServerFn_handler = createServerRpc("fc3988c64f434639dfd4eab3f926b87ee39cc0c14f65b4d0e852c7fd73279a3b", (opts, signal) => {
  return slowServerFn.__executeServer(opts, signal);
});
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
const Route = createFileRoute("/deferred")({
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
  component: Deferred
});
function Deferred() {
  const [count, setCount] = useState(0);
  const {
    deferredStuff,
    deferredPerson,
    person
  } = Route.useLoaderData();
  return /* @__PURE__ */ jsxs("div", { className: "p-2", children: [
    /* @__PURE__ */ jsxs("div", { "data-testid": "regular-person", children: [
      person.name,
      " - ",
      person.randomNumber
    ] }),
    /* @__PURE__ */ jsx(Suspense, { fallback: /* @__PURE__ */ jsx("div", { children: "Loading person..." }), children: /* @__PURE__ */ jsx(Await, { promise: deferredPerson, children: (data) => /* @__PURE__ */ jsxs("div", { "data-testid": "deferred-person", children: [
      data.name,
      " - ",
      data.randomNumber
    ] }) }) }),
    /* @__PURE__ */ jsx(Suspense, { fallback: /* @__PURE__ */ jsx("div", { children: "Loading stuff..." }), children: /* @__PURE__ */ jsx(Await, { promise: deferredStuff, children: (data) => /* @__PURE__ */ jsx("h3", { "data-testid": "deferred-stuff", children: data }) }) }),
    /* @__PURE__ */ jsxs("div", { children: [
      "Count: ",
      count
    ] }),
    /* @__PURE__ */ jsx("div", { children: /* @__PURE__ */ jsx("button", { onClick: () => setCount(count + 1), children: "Increment" }) })
  ] });
}
export {
  personServerFn_createServerFn_handler,
  slowServerFn_createServerFn_handler
};
