import { jsx, Fragment, jsxs } from "react/jsx-runtime";
import { exactPathTest, removeTrailingSlash, deepEqual, preloadWarning, functionalUpdate, BaseRootRoute, BaseRoute, isModuleNotFoundError, RouterCore, rootRouteId } from "@tanstack/router-core";
import warning from "tiny-warning";
import * as React from "react";
import React__default, { createElement, useRef, useEffect, useCallback, forwardRef, useState, useImperativeHandle, useMemo } from "react";
import invariant from "tiny-invariant";
import { d as dummyMatchContext, m as matchContext, u as useRouterState, a as useRouter, b as useForwardedRef, c as useIntersectionObserver, E as ErrorComponent } from "../server.js";
import { flushSync } from "react-dom";
import { PanelGroup, Panel, PanelResizeHandle } from "react-resizable-panels";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { Schema } from "effect";
import { ResponsiveContainer, LineChart, CartesianGrid, XAxis, YAxis, Tooltip, Legend, Line, ReferenceLine } from "recharts";
import { SkipBack, ChevronLeft, Pause, Play, ChevronRight, SkipForward, Settings, ChevronDown, RefreshCw, Folder } from "lucide-react";
import { json } from "@tanstack/router-core/ssr/client";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
function useMatch(opts) {
  const nearestMatchId = React.useContext(
    opts.from ? dummyMatchContext : matchContext
  );
  const matchSelection = useRouterState({
    select: (state) => {
      const match = state.matches.find(
        (d) => opts.from ? opts.from === d.routeId : d.id === nearestMatchId
      );
      invariant(
        !((opts.shouldThrow ?? true) && !match),
        `Could not find ${opts.from ? `an active match from "${opts.from}"` : "a nearest match!"}`
      );
      if (match === void 0) {
        return void 0;
      }
      return opts.select ? opts.select(match) : match;
    },
    structuralSharing: opts.structuralSharing
  });
  return matchSelection;
}
function useLoaderData(opts) {
  return useMatch({
    from: opts.from,
    strict: opts.strict,
    structuralSharing: opts.structuralSharing,
    select: (s) => {
      return opts.select ? opts.select(s.loaderData) : s.loaderData;
    }
  });
}
function useLoaderDeps(opts) {
  const { select, ...rest } = opts;
  return useMatch({
    ...rest,
    select: (s) => {
      return select ? select(s.loaderDeps) : s.loaderDeps;
    }
  });
}
function useParams(opts) {
  return useMatch({
    from: opts.from,
    shouldThrow: opts.shouldThrow,
    structuralSharing: opts.structuralSharing,
    strict: opts.strict,
    select: (match) => {
      const params = opts.strict === false ? match.params : match._strictParams;
      return opts.select ? opts.select(params) : params;
    }
  });
}
function useSearch(opts) {
  return useMatch({
    from: opts.from,
    strict: opts.strict,
    shouldThrow: opts.shouldThrow,
    structuralSharing: opts.structuralSharing,
    select: (match) => {
      return opts.select ? opts.select(match.search) : match.search;
    }
  });
}
function useNavigate(_defaultOpts) {
  const router2 = useRouter();
  return React.useCallback(
    (options) => {
      return router2.navigate({
        ...options,
        from: options.from ?? _defaultOpts?.from
      });
    },
    [_defaultOpts?.from, router2]
  );
}
function useLinkProps(options, forwardedRef) {
  const router2 = useRouter();
  const [isTransitioning, setIsTransitioning] = React.useState(false);
  const hasRenderFetched = React.useRef(false);
  const innerRef = useForwardedRef(forwardedRef);
  const {
    // custom props
    activeProps,
    inactiveProps,
    activeOptions,
    to,
    preload: userPreload,
    preloadDelay: userPreloadDelay,
    hashScrollIntoView,
    replace,
    startTransition,
    resetScroll,
    viewTransition,
    // element props
    children,
    target,
    disabled,
    style,
    className,
    onClick,
    onFocus,
    onMouseEnter,
    onMouseLeave,
    onTouchStart,
    ignoreBlocker,
    // prevent these from being returned
    params: _params,
    search: _search,
    hash: _hash,
    state: _state,
    mask: _mask,
    reloadDocument: _reloadDocument,
    unsafeRelative: _unsafeRelative,
    from: _from,
    _fromLocation,
    ...propsSafeToSpread
  } = options;
  const currentSearch = useRouterState({
    select: (s) => s.location.search,
    structuralSharing: true
  });
  const from = options.from;
  const _options = React.useMemo(
    () => {
      return { ...options, from };
    },
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [
      router2,
      currentSearch,
      from,
      options._fromLocation,
      options.hash,
      options.to,
      options.search,
      options.params,
      options.state,
      options.mask,
      options.unsafeRelative
    ]
  );
  const next = React.useMemo(
    () => router2.buildLocation({ ..._options }),
    [router2, _options]
  );
  const hrefOption = React.useMemo(() => {
    if (disabled) {
      return void 0;
    }
    let href = next.maskedLocation ? next.maskedLocation.url : next.url;
    let external = false;
    if (router2.origin) {
      if (href.startsWith(router2.origin)) {
        href = router2.history.createHref(href.replace(router2.origin, "")) || "/";
      } else {
        external = true;
      }
    }
    return { href, external };
  }, [disabled, next.maskedLocation, next.url, router2.origin, router2.history]);
  const externalLink = React.useMemo(() => {
    if (hrefOption?.external) {
      return hrefOption.href;
    }
    try {
      new URL(to);
      return to;
    } catch {
    }
    return void 0;
  }, [to, hrefOption]);
  const preload = options.reloadDocument || externalLink ? false : userPreload ?? router2.options.defaultPreload;
  const preloadDelay = userPreloadDelay ?? router2.options.defaultPreloadDelay ?? 0;
  const isActive = useRouterState({
    select: (s) => {
      if (externalLink) return false;
      if (activeOptions?.exact) {
        const testExact = exactPathTest(
          s.location.pathname,
          next.pathname,
          router2.basepath
        );
        if (!testExact) {
          return false;
        }
      } else {
        const currentPathSplit = removeTrailingSlash(
          s.location.pathname,
          router2.basepath
        );
        const nextPathSplit = removeTrailingSlash(
          next.pathname,
          router2.basepath
        );
        const pathIsFuzzyEqual = currentPathSplit.startsWith(nextPathSplit) && (currentPathSplit.length === nextPathSplit.length || currentPathSplit[nextPathSplit.length] === "/");
        if (!pathIsFuzzyEqual) {
          return false;
        }
      }
      if (activeOptions?.includeSearch ?? true) {
        const searchTest = deepEqual(s.location.search, next.search, {
          partial: !activeOptions?.exact,
          ignoreUndefined: !activeOptions?.explicitUndefined
        });
        if (!searchTest) {
          return false;
        }
      }
      if (activeOptions?.includeHash) {
        return s.location.hash === next.hash;
      }
      return true;
    }
  });
  const doPreload = React.useCallback(() => {
    router2.preloadRoute({ ..._options }).catch((err) => {
      console.warn(err);
      console.warn(preloadWarning);
    });
  }, [router2, _options]);
  const preloadViewportIoCallback = React.useCallback(
    (entry) => {
      if (entry?.isIntersecting) {
        doPreload();
      }
    },
    [doPreload]
  );
  useIntersectionObserver(
    innerRef,
    preloadViewportIoCallback,
    intersectionObserverOptions,
    { disabled: !!disabled || !(preload === "viewport") }
  );
  React.useEffect(() => {
    if (hasRenderFetched.current) {
      return;
    }
    if (!disabled && preload === "render") {
      doPreload();
      hasRenderFetched.current = true;
    }
  }, [disabled, doPreload, preload]);
  const handleClick = (e) => {
    const elementTarget = e.currentTarget.getAttribute("target");
    const effectiveTarget = target !== void 0 ? target : elementTarget;
    if (!disabled && !isCtrlEvent(e) && !e.defaultPrevented && (!effectiveTarget || effectiveTarget === "_self") && e.button === 0) {
      e.preventDefault();
      flushSync(() => {
        setIsTransitioning(true);
      });
      const unsub = router2.subscribe("onResolved", () => {
        unsub();
        setIsTransitioning(false);
      });
      router2.navigate({
        ..._options,
        replace,
        resetScroll,
        hashScrollIntoView,
        startTransition,
        viewTransition,
        ignoreBlocker
      });
    }
  };
  if (externalLink) {
    return {
      ...propsSafeToSpread,
      ref: innerRef,
      href: externalLink,
      ...children && { children },
      ...target && { target },
      ...disabled && { disabled },
      ...style && { style },
      ...className && { className },
      ...onClick && { onClick },
      ...onFocus && { onFocus },
      ...onMouseEnter && { onMouseEnter },
      ...onMouseLeave && { onMouseLeave },
      ...onTouchStart && { onTouchStart }
    };
  }
  const handleFocus = (_) => {
    if (disabled) return;
    if (preload) {
      doPreload();
    }
  };
  const handleTouchStart = handleFocus;
  const handleEnter = (e) => {
    if (disabled || !preload) return;
    if (!preloadDelay) {
      doPreload();
    } else {
      const eventTarget = e.target;
      if (timeoutMap.has(eventTarget)) {
        return;
      }
      const id = setTimeout(() => {
        timeoutMap.delete(eventTarget);
        doPreload();
      }, preloadDelay);
      timeoutMap.set(eventTarget, id);
    }
  };
  const handleLeave = (e) => {
    if (disabled || !preload || !preloadDelay) return;
    const eventTarget = e.target;
    const id = timeoutMap.get(eventTarget);
    if (id) {
      clearTimeout(id);
      timeoutMap.delete(eventTarget);
    }
  };
  const resolvedActiveProps = isActive ? functionalUpdate(activeProps, {}) ?? STATIC_ACTIVE_OBJECT : STATIC_EMPTY_OBJECT;
  const resolvedInactiveProps = isActive ? STATIC_EMPTY_OBJECT : functionalUpdate(inactiveProps, {}) ?? STATIC_EMPTY_OBJECT;
  const resolvedClassName = [
    className,
    resolvedActiveProps.className,
    resolvedInactiveProps.className
  ].filter(Boolean).join(" ");
  const resolvedStyle = (style || resolvedActiveProps.style || resolvedInactiveProps.style) && {
    ...style,
    ...resolvedActiveProps.style,
    ...resolvedInactiveProps.style
  };
  return {
    ...propsSafeToSpread,
    ...resolvedActiveProps,
    ...resolvedInactiveProps,
    href: hrefOption?.href,
    ref: innerRef,
    onClick: composeHandlers([onClick, handleClick]),
    onFocus: composeHandlers([onFocus, handleFocus]),
    onMouseEnter: composeHandlers([onMouseEnter, handleEnter]),
    onMouseLeave: composeHandlers([onMouseLeave, handleLeave]),
    onTouchStart: composeHandlers([onTouchStart, handleTouchStart]),
    disabled: !!disabled,
    target,
    ...resolvedStyle && { style: resolvedStyle },
    ...resolvedClassName && { className: resolvedClassName },
    ...disabled && STATIC_DISABLED_PROPS,
    ...isActive && STATIC_ACTIVE_PROPS,
    ...isTransitioning && STATIC_TRANSITIONING_PROPS
  };
}
const STATIC_EMPTY_OBJECT = {};
const STATIC_ACTIVE_OBJECT = { className: "active" };
const STATIC_DISABLED_PROPS = { role: "link", "aria-disabled": true };
const STATIC_ACTIVE_PROPS = { "data-status": "active", "aria-current": "page" };
const STATIC_TRANSITIONING_PROPS = { "data-transitioning": "transitioning" };
const timeoutMap = /* @__PURE__ */ new WeakMap();
const intersectionObserverOptions = {
  rootMargin: "100px"
};
const composeHandlers = (handlers) => (e) => {
  for (const handler of handlers) {
    if (!handler) continue;
    if (e.defaultPrevented) return;
    handler(e);
  }
};
const Link = React.forwardRef(
  (props, ref) => {
    const { _asChild, ...rest } = props;
    const {
      type: _type,
      ref: innerRef,
      ...linkProps
    } = useLinkProps(rest, ref);
    const children = typeof rest.children === "function" ? rest.children({
      isActive: linkProps["data-status"] === "active"
    }) : rest.children;
    if (_asChild === void 0) {
      delete linkProps.disabled;
    }
    return React.createElement(
      _asChild ? _asChild : "a",
      {
        ...linkProps,
        ref: innerRef
      },
      children
    );
  }
);
function isCtrlEvent(e) {
  return !!(e.metaKey || e.altKey || e.ctrlKey || e.shiftKey);
}
let Route$7 = class Route extends BaseRoute {
  /**
   * @deprecated Use the `createRoute` function instead.
   */
  constructor(options) {
    super(options);
    this.useMatch = (opts) => {
      return useMatch({
        select: opts?.select,
        from: this.id,
        structuralSharing: opts?.structuralSharing
      });
    };
    this.useRouteContext = (opts) => {
      return useMatch({
        ...opts,
        from: this.id,
        select: (d) => opts?.select ? opts.select(d.context) : d.context
      });
    };
    this.useSearch = (opts) => {
      return useSearch({
        select: opts?.select,
        structuralSharing: opts?.structuralSharing,
        from: this.id
      });
    };
    this.useParams = (opts) => {
      return useParams({
        select: opts?.select,
        structuralSharing: opts?.structuralSharing,
        from: this.id
      });
    };
    this.useLoaderDeps = (opts) => {
      return useLoaderDeps({ ...opts, from: this.id });
    };
    this.useLoaderData = (opts) => {
      return useLoaderData({ ...opts, from: this.id });
    };
    this.useNavigate = () => {
      return useNavigate({ from: this.fullPath });
    };
    this.Link = React__default.forwardRef(
      (props, ref) => {
        return /* @__PURE__ */ jsx(Link, { ref, from: this.fullPath, ...props });
      }
    );
    this.$$typeof = Symbol.for("react.memo");
  }
};
function createRoute(options) {
  return new Route$7(
    // TODO: Help us TypeChris, you're our only hope!
    options
  );
}
class RootRoute extends BaseRootRoute {
  /**
   * @deprecated `RootRoute` is now an internal implementation detail. Use `createRootRoute()` instead.
   */
  constructor(options) {
    super(options);
    this.useMatch = (opts) => {
      return useMatch({
        select: opts?.select,
        from: this.id,
        structuralSharing: opts?.structuralSharing
      });
    };
    this.useRouteContext = (opts) => {
      return useMatch({
        ...opts,
        from: this.id,
        select: (d) => opts?.select ? opts.select(d.context) : d.context
      });
    };
    this.useSearch = (opts) => {
      return useSearch({
        select: opts?.select,
        structuralSharing: opts?.structuralSharing,
        from: this.id
      });
    };
    this.useParams = (opts) => {
      return useParams({
        select: opts?.select,
        structuralSharing: opts?.structuralSharing,
        from: this.id
      });
    };
    this.useLoaderDeps = (opts) => {
      return useLoaderDeps({ ...opts, from: this.id });
    };
    this.useLoaderData = (opts) => {
      return useLoaderData({ ...opts, from: this.id });
    };
    this.useNavigate = () => {
      return useNavigate({ from: this.fullPath });
    };
    this.Link = React__default.forwardRef(
      (props, ref) => {
        return /* @__PURE__ */ jsx(Link, { ref, from: this.fullPath, ...props });
      }
    );
    this.$$typeof = Symbol.for("react.memo");
  }
}
function createRootRoute(options) {
  return new RootRoute(options);
}
function createFileRoute(path2) {
  if (typeof path2 === "object") {
    return new FileRoute(path2, {
      silent: true
    }).createRoute(path2);
  }
  return new FileRoute(path2, {
    silent: true
  }).createRoute;
}
class FileRoute {
  constructor(path2, _opts) {
    this.path = path2;
    this.createRoute = (options) => {
      warning(
        this.silent,
        "FileRoute is deprecated and will be removed in the next major version. Use the createFileRoute(path)(options) function instead."
      );
      const route = createRoute(options);
      route.isRoot = false;
      return route;
    };
    this.silent = _opts?.silent;
  }
}
class LazyRoute {
  constructor(opts) {
    this.useMatch = (opts2) => {
      return useMatch({
        select: opts2?.select,
        from: this.options.id,
        structuralSharing: opts2?.structuralSharing
      });
    };
    this.useRouteContext = (opts2) => {
      return useMatch({
        from: this.options.id,
        select: (d) => opts2?.select ? opts2.select(d.context) : d.context
      });
    };
    this.useSearch = (opts2) => {
      return useSearch({
        select: opts2?.select,
        structuralSharing: opts2?.structuralSharing,
        from: this.options.id
      });
    };
    this.useParams = (opts2) => {
      return useParams({
        select: opts2?.select,
        structuralSharing: opts2?.structuralSharing,
        from: this.options.id
      });
    };
    this.useLoaderDeps = (opts2) => {
      return useLoaderDeps({ ...opts2, from: this.options.id });
    };
    this.useLoaderData = (opts2) => {
      return useLoaderData({ ...opts2, from: this.options.id });
    };
    this.useNavigate = () => {
      const router2 = useRouter();
      return useNavigate({ from: router2.routesById[this.options.id].fullPath });
    };
    this.options = opts;
    this.$$typeof = Symbol.for("react.memo");
  }
}
function createLazyFileRoute(id) {
  if (typeof id === "object") {
    return new LazyRoute(id);
  }
  return (opts) => new LazyRoute({ id, ...opts });
}
function lazyRouteComponent(importer, exportName) {
  let loadPromise;
  let comp;
  let error;
  let reload;
  const load = () => {
    if (!loadPromise) {
      loadPromise = importer().then((res) => {
        loadPromise = void 0;
        comp = res[exportName];
      }).catch((err) => {
        error = err;
        if (isModuleNotFoundError(error)) {
          if (error instanceof Error && typeof window !== "undefined" && typeof sessionStorage !== "undefined") {
            const storageKey = `tanstack_router_reload:${error.message}`;
            if (!sessionStorage.getItem(storageKey)) {
              sessionStorage.setItem(storageKey, "1");
              reload = true;
            }
          }
        }
      });
    }
    return loadPromise;
  };
  const lazyComp = function Lazy(props) {
    if (reload) {
      window.location.reload();
      throw new Promise(() => {
      });
    }
    if (error) {
      throw error;
    }
    if (!comp) {
      throw load();
    }
    return React.createElement(comp, props);
  };
  lazyComp.preload = load;
  return lazyComp;
}
const createRouter = (options) => {
  return new Router(options);
};
class Router extends RouterCore {
  constructor(options) {
    super(options);
  }
}
if (typeof globalThis !== "undefined") {
  globalThis.createFileRoute = createFileRoute;
  globalThis.createLazyFileRoute = createLazyFileRoute;
} else if (typeof window !== "undefined") {
  window.createFileRoute = createFileRoute;
  window.createLazyFileRoute = createLazyFileRoute;
}
function Asset({
  tag,
  attrs,
  children,
  nonce
}) {
  switch (tag) {
    case "title":
      return /* @__PURE__ */ jsx("title", { ...attrs, suppressHydrationWarning: true, children });
    case "meta":
      return /* @__PURE__ */ jsx("meta", { ...attrs, suppressHydrationWarning: true });
    case "link":
      return /* @__PURE__ */ jsx("link", { ...attrs, nonce, suppressHydrationWarning: true });
    case "style":
      return /* @__PURE__ */ jsx(
        "style",
        {
          ...attrs,
          dangerouslySetInnerHTML: { __html: children },
          nonce
        }
      );
    case "script":
      return /* @__PURE__ */ jsx(Script, { attrs, children });
    default:
      return null;
  }
}
function Script({
  attrs,
  children
}) {
  const router2 = useRouter();
  React.useEffect(() => {
    if (attrs?.src) {
      const normSrc = (() => {
        try {
          const base = document.baseURI || window.location.href;
          return new URL(attrs.src, base).href;
        } catch {
          return attrs.src;
        }
      })();
      const existingScript = Array.from(
        document.querySelectorAll("script[src]")
      ).find((el) => el.src === normSrc);
      if (existingScript) {
        return;
      }
      const script = document.createElement("script");
      for (const [key, value] of Object.entries(attrs)) {
        if (key !== "suppressHydrationWarning" && value !== void 0 && value !== false) {
          script.setAttribute(
            key,
            typeof value === "boolean" ? "" : String(value)
          );
        }
      }
      document.head.appendChild(script);
      return () => {
        if (script.parentNode) {
          script.parentNode.removeChild(script);
        }
      };
    }
    if (typeof children === "string") {
      const typeAttr = typeof attrs?.type === "string" ? attrs.type : "text/javascript";
      const nonceAttr = typeof attrs?.nonce === "string" ? attrs.nonce : void 0;
      const existingScript = Array.from(
        document.querySelectorAll("script:not([src])")
      ).find((el) => {
        if (!(el instanceof HTMLScriptElement)) return false;
        const sType = el.getAttribute("type") ?? "text/javascript";
        const sNonce = el.getAttribute("nonce") ?? void 0;
        return el.textContent === children && sType === typeAttr && sNonce === nonceAttr;
      });
      if (existingScript) {
        return;
      }
      const script = document.createElement("script");
      script.textContent = children;
      if (attrs) {
        for (const [key, value] of Object.entries(attrs)) {
          if (key !== "suppressHydrationWarning" && value !== void 0 && value !== false) {
            script.setAttribute(
              key,
              typeof value === "boolean" ? "" : String(value)
            );
          }
        }
      }
      document.head.appendChild(script);
      return () => {
        if (script.parentNode) {
          script.parentNode.removeChild(script);
        }
      };
    }
    return void 0;
  }, [attrs, children]);
  if (!router2.isServer) {
    const { src, ...rest } = attrs || {};
    return /* @__PURE__ */ jsx(
      "script",
      {
        suppressHydrationWarning: true,
        dangerouslySetInnerHTML: { __html: "" },
        ...rest
      }
    );
  }
  if (attrs?.src && typeof attrs.src === "string") {
    return /* @__PURE__ */ jsx("script", { ...attrs, suppressHydrationWarning: true });
  }
  if (typeof children === "string") {
    return /* @__PURE__ */ jsx(
      "script",
      {
        ...attrs,
        dangerouslySetInnerHTML: { __html: children },
        suppressHydrationWarning: true
      }
    );
  }
  return null;
}
const useTags = () => {
  const router2 = useRouter();
  const nonce = router2.options.ssr?.nonce;
  const routeMeta = useRouterState({
    select: (state) => {
      return state.matches.map((match) => match.meta).filter(Boolean);
    }
  });
  const meta = React.useMemo(() => {
    const resultMeta = [];
    const metaByAttribute = {};
    let title;
    for (let i = routeMeta.length - 1; i >= 0; i--) {
      const metas = routeMeta[i];
      for (let j = metas.length - 1; j >= 0; j--) {
        const m = metas[j];
        if (!m) continue;
        if (m.title) {
          if (!title) {
            title = {
              tag: "title",
              children: m.title
            };
          }
        } else {
          const attribute = m.name ?? m.property;
          if (attribute) {
            if (metaByAttribute[attribute]) {
              continue;
            } else {
              metaByAttribute[attribute] = true;
            }
          }
          resultMeta.push({
            tag: "meta",
            attrs: {
              ...m,
              nonce
            }
          });
        }
      }
    }
    if (title) {
      resultMeta.push(title);
    }
    if (nonce) {
      resultMeta.push({
        tag: "meta",
        attrs: {
          property: "csp-nonce",
          content: nonce
        }
      });
    }
    resultMeta.reverse();
    return resultMeta;
  }, [routeMeta, nonce]);
  const links = useRouterState({
    select: (state) => {
      const constructed = state.matches.map((match) => match.links).filter(Boolean).flat(1).map((link) => ({
        tag: "link",
        attrs: {
          ...link,
          nonce
        }
      }));
      const manifest = router2.ssr?.manifest;
      const assets = state.matches.map((match) => manifest?.routes[match.routeId]?.assets ?? []).filter(Boolean).flat(1).filter((asset) => asset.tag === "link").map(
        (asset) => ({
          tag: "link",
          attrs: {
            ...asset.attrs,
            suppressHydrationWarning: true,
            nonce
          }
        })
      );
      return [...constructed, ...assets];
    },
    structuralSharing: true
  });
  const preloadLinks = useRouterState({
    select: (state) => {
      const preloadLinks2 = [];
      state.matches.map((match) => router2.looseRoutesById[match.routeId]).forEach(
        (route) => router2.ssr?.manifest?.routes[route.id]?.preloads?.filter(Boolean).forEach((preload) => {
          preloadLinks2.push({
            tag: "link",
            attrs: {
              rel: "modulepreload",
              href: preload,
              nonce
            }
          });
        })
      );
      return preloadLinks2;
    },
    structuralSharing: true
  });
  const styles = useRouterState({
    select: (state) => state.matches.map((match) => match.styles).flat(1).filter(Boolean).map(({ children, ...attrs }) => ({
      tag: "style",
      attrs,
      children,
      nonce
    })),
    structuralSharing: true
  });
  const headScripts = useRouterState({
    select: (state) => state.matches.map((match) => match.headScripts).flat(1).filter(Boolean).map(({ children, ...script }) => ({
      tag: "script",
      attrs: {
        ...script,
        nonce
      },
      children
    })),
    structuralSharing: true
  });
  return uniqBy(
    [
      ...meta,
      ...preloadLinks,
      ...links,
      ...styles,
      ...headScripts
    ],
    (d) => {
      return JSON.stringify(d);
    }
  );
};
function HeadContent() {
  const tags = useTags();
  const router2 = useRouter();
  const nonce = router2.options.ssr?.nonce;
  return tags.map((tag) => /* @__PURE__ */ createElement(Asset, { ...tag, key: `tsr-meta-${JSON.stringify(tag)}`, nonce }));
}
function uniqBy(arr, fn) {
  const seen = /* @__PURE__ */ new Set();
  return arr.filter((item) => {
    const key = fn(item);
    if (seen.has(key)) {
      return false;
    }
    seen.add(key);
    return true;
  });
}
const Scripts = () => {
  const router2 = useRouter();
  const nonce = router2.options.ssr?.nonce;
  const assetScripts = useRouterState({
    select: (state) => {
      const assetScripts2 = [];
      const manifest = router2.ssr?.manifest;
      if (!manifest) {
        return [];
      }
      state.matches.map((match) => router2.looseRoutesById[match.routeId]).forEach(
        (route) => manifest.routes[route.id]?.assets?.filter((d) => d.tag === "script").forEach((asset) => {
          assetScripts2.push({
            tag: "script",
            attrs: { ...asset.attrs, nonce },
            children: asset.children
          });
        })
      );
      return assetScripts2;
    },
    structuralSharing: true
  });
  const { scripts } = useRouterState({
    select: (state) => ({
      scripts: state.matches.map((match) => match.scripts).flat(1).filter(Boolean).map(({ children, ...script }) => ({
        tag: "script",
        attrs: {
          ...script,
          suppressHydrationWarning: true,
          nonce
        },
        children
      }))
    }),
    structuralSharing: true
  });
  let serverBufferedScript = void 0;
  if (router2.serverSsr) {
    serverBufferedScript = router2.serverSsr.takeBufferedScripts();
  }
  const allScripts = [...scripts, ...assetScripts];
  if (serverBufferedScript) {
    allScripts.unshift(serverBufferedScript);
  }
  return /* @__PURE__ */ jsx(Fragment, { children: allScripts.map((asset, i) => /* @__PURE__ */ createElement(Asset, { ...asset, key: `tsr-scripts-${asset.tag}-${i}` })) });
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
const appCss = "/assets/app-dS8jA3co.css";
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
function DevtoolsWrapper() {
  const [DevtoolsComponent, setDevtoolsComponent] = React.useState(null);
  React.useEffect(() => {
    if (process.env.NODE_ENV !== "production") {
      import("@tanstack/react-router-devtools").then((mod) => {
        setDevtoolsComponent(() => mod.TanStackRouterDevtools);
      });
    }
  }, []);
  if (!DevtoolsComponent) return null;
  return /* @__PURE__ */ jsx(DevtoolsComponent, { position: "bottom-right" });
}
const Route$6 = createRootRoute({
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
      /* @__PURE__ */ jsx(DevtoolsWrapper, {}),
      /* @__PURE__ */ jsx(Scripts, {})
    ] })
  ] });
}
const $$splitComponentImporter = () => import("./index-Cm4Dx40N.js");
const Route$5 = createFileRoute("/")({
  component: lazyRouteComponent($$splitComponentImporter, "component")
});
const DEFAULT_GALACTIC_CONFIG = {
  // Physical distances
  distanceToGC_pc: 60,
  // ~60 pc from Galactic Center (projected)
  distanceToGC: 60,
  // Alias
  distanceToEarth_kpc: 8,
  // 8 kpc to Galactic Center
  // Source position
  galacticLongitude: -0.398,
  // l = -0.398° (or 359.602°)
  galacticLatitude: -0.224,
  // b = -0.224°
  // Viewing geometry
  inclination: 70,
  // 70° inclination
  positionAngle: 41.6,
  // Position angle 41.6°
  // Kinematics
  lsrVelocity: -60,
  // V_LSR ~ -60 km/s (mean of -80 to -40 km/s range)
  // Cloud sizes (TRUE physical values from Oka et al. 2017)
  hvccSize_pc: 0.15,
  // HCN clump half-size S = 0.15 pc
  hvccSize: 0.15,
  // Alias
  hvccPosition: [0.2, -0.1, 0],
  // Position offset from BH
  cloudExtent_pc: 5,
  // Extended cloud ~5 pc
  denseClumpOffset_pc: 0.2,
  // Displacement ~0.2 pc
  // N-body model parameters
  nbodyCloudSigma_pc: 0.2,
  // σ_r = 0.2 pc
  bhMass_Msun: 1e5,
  // 10^5 M☉
  cloudMass_Msun: 1e3,
  // 1000 M☉
  // Simulation cloud parameters (default from CAT_OKA preset)
  cloudRadius_pc: 1.13,
  // Cloud radius from preset (1.13 pc)
  // Display scale (1.0 = true scale in pc)
  displayScale: 1,
  // Galactic visualization
  showGalaxyDisk: true,
  showSolarSystem: false,
  // Solar System demo off by default
  galaxyRotationSpeed: 0
  // Static by default
};
function createTextSprite$1(text, options = {}) {
  const {
    color = "#ffffff",
    fontSize = 72,
    // Increased from 48 for better visibility
    fontWeight = "Bold",
    backgroundColor,
    padding = 10,
    scale = [12, 6, 1]
    // Increased from [8, 4, 1] for better visibility
  } = options;
  const canvas = document.createElement("canvas");
  const size = 512;
  canvas.width = size;
  canvas.height = size / 2;
  const ctx = canvas.getContext("2d");
  if (backgroundColor) {
    ctx.fillStyle = backgroundColor;
    ctx.fillRect(0, 0, size, size / 2);
  }
  ctx.fillStyle = color;
  ctx.font = `${fontWeight} ${fontSize}px Arial, sans-serif`;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, size / 2, size / 4);
  const texture = new THREE.CanvasTexture(canvas);
  texture.needsUpdate = true;
  const material = new THREE.SpriteMaterial({
    map: texture,
    transparent: true,
    depthTest: false,
    depthWrite: false
  });
  const sprite = new THREE.Sprite(material);
  sprite.scale.set(...scale);
  sprite.renderOrder = 999;
  return sprite;
}
function createArrow(from, to, color, options = {}) {
  const { headLength = 1, headWidth = 0.6, shaftRadius = 0.15, opacity = 1 } = options;
  const group = new THREE.Group();
  const direction = to.clone().sub(from).normalize();
  const length = from.distanceTo(to);
  if (length <= headLength) return group;
  const shaftLength = length - headLength;
  const shaftGeo = new THREE.CylinderGeometry(shaftRadius, shaftRadius, shaftLength, 12);
  const shaftMat = new THREE.MeshBasicMaterial({ color, transparent: opacity < 1, opacity });
  const shaft = new THREE.Mesh(shaftGeo, shaftMat);
  const shaftMidpoint = from.clone().add(direction.clone().multiplyScalar(shaftLength / 2));
  shaft.position.copy(shaftMidpoint);
  shaft.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);
  group.add(shaft);
  const headGeo = new THREE.ConeGeometry(headWidth, headLength, 16);
  const headMat = new THREE.MeshBasicMaterial({ color, transparent: opacity < 1, opacity });
  const head = new THREE.Mesh(headGeo, headMat);
  const headPosition = to.clone().sub(direction.clone().multiplyScalar(headLength / 2));
  head.position.copy(headPosition);
  head.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);
  group.add(head);
  return group;
}
function createDashedCircle(radius, color, options = {}) {
  const {
    segments = 64,
    dashSize = 0.5,
    gapSize = 0.25,
    opacity = 0.8,
    normal = new THREE.Vector3(0, 0, 1),
    center = new THREE.Vector3(0, 0, 0)
  } = options;
  const points = [];
  const defaultNormal = new THREE.Vector3(0, 0, 1);
  const quaternion = new THREE.Quaternion().setFromUnitVectors(defaultNormal, normal.clone().normalize());
  for (let i = 0; i <= segments; i++) {
    const theta = i / segments * Math.PI * 2;
    const point = new THREE.Vector3(
      radius * Math.cos(theta),
      radius * Math.sin(theta),
      0
    );
    point.applyQuaternion(quaternion);
    point.add(center);
    points.push(point);
  }
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const material = new THREE.LineDashedMaterial({
    color,
    dashSize,
    gapSize,
    transparent: true,
    opacity
  });
  const line = new THREE.Line(geometry, material);
  line.computeLineDistances();
  return line;
}
function createGalaxyDisk(center, scale = 1, config = {}) {
  const {
    showSpiralArms = true,
    opacity = 0.4
    // Cloud-like opacity for visibility
  } = config;
  const group = new THREE.Group();
  group.position.copy(center);
  group.name = "galaxyDisk";
  const displayScale = 1 / 300;
  const diskRadius = 15e3 * displayScale * scale;
  const bulgeRadius = 2500 * displayScale * scale;
  const sunDistance = 8e3 * displayScale * scale;
  const diskSegments = 128;
  const diskGeo = new THREE.CircleGeometry(diskRadius, diskSegments);
  const diskColors = new Float32Array((diskSegments + 2) * 3);
  diskColors[0] = 1;
  diskColors[1] = 0.95;
  diskColors[2] = 0.8;
  for (let i = 1; i <= diskSegments + 1; i++) {
    const t = 1 - i / diskSegments * 0.3;
    diskColors[i * 3] = 0.2 * t;
    diskColors[i * 3 + 1] = 0.25 * t;
    diskColors[i * 3 + 2] = 0.4 * t;
  }
  diskGeo.setAttribute("color", new THREE.BufferAttribute(diskColors, 3));
  const diskMat = new THREE.MeshBasicMaterial({
    vertexColors: true,
    transparent: true,
    opacity,
    side: THREE.DoubleSide
  });
  const disk = new THREE.Mesh(diskGeo, diskMat);
  disk.rotation.x = -Math.PI / 2;
  group.add(disk);
  const bulgeGeo = new THREE.CircleGeometry(bulgeRadius, 64);
  const bulgeMat = new THREE.MeshBasicMaterial({
    color: 16772812,
    transparent: true,
    opacity: opacity * 2,
    side: THREE.DoubleSide
  });
  const bulge = new THREE.Mesh(bulgeGeo, bulgeMat);
  bulge.rotation.x = -Math.PI / 2;
  bulge.position.y = 0.01;
  group.add(bulge);
  if (showSpiralArms) {
    const armColor = 4491468;
    const numArms = 4;
    const armTurns = 1.5;
    const armPoints = 100;
    for (let arm = 0; arm < numArms; arm++) {
      const armOffset = arm / numArms * Math.PI * 2;
      const points = [];
      for (let i = 0; i < armPoints; i++) {
        const t = i / armPoints;
        const r = (bulgeRadius + t * (diskRadius - bulgeRadius)) * 0.9;
        const theta = armOffset + t * armTurns * Math.PI * 2;
        const x = r * Math.cos(theta);
        const z = r * Math.sin(theta);
        points.push(new THREE.Vector3(x, 0.02, z));
      }
      const armGeo = new THREE.BufferGeometry().setFromPoints(points);
      const armMat = new THREE.LineBasicMaterial({
        color: armColor,
        transparent: true,
        opacity: opacity * 1.5
      });
      const armLine = new THREE.Line(armGeo, armMat);
      group.add(armLine);
    }
  }
  const gcGeo = new THREE.SphereGeometry(0.5 * scale, 32, 32);
  const gcMat = new THREE.MeshBasicMaterial({
    color: 16755200,
    transparent: true,
    opacity: 0.3
  });
  const gcMarker = new THREE.Mesh(gcGeo, gcMat);
  group.add(gcMarker);
  const sunPosInGalaxy = new THREE.Vector3(-sunDistance, 0.05, 0);
  const sunDotGeo = new THREE.SphereGeometry(0.4 * scale, 16, 16);
  const sunDotMat = new THREE.MeshBasicMaterial({
    color: 16777028,
    transparent: true,
    opacity: 0.4
  });
  const sunDot = new THREE.Mesh(sunDotGeo, sunDotMat);
  sunDot.position.copy(sunPosInGalaxy);
  group.add(sunDot);
  const sunPosLabel = createTextSprite$1("☉ (8 kpc)", { color: "#ffff88", fontSize: 18 });
  sunPosLabel.position.copy(sunPosInGalaxy.clone().add(new THREE.Vector3(0, 1.2 * scale, 0)));
  sunPosLabel.scale.set(4, 2, 1);
  group.add(sunPosLabel);
  const sunOrbit = createDashedCircle(sunDistance, 16777028, {
    segments: 64,
    dashSize: 1,
    gapSize: 0.5,
    opacity: opacity * 2,
    normal: new THREE.Vector3(0, 1, 0),
    center: new THREE.Vector3(0, 0, 0)
  });
  group.add(sunOrbit);
  const hvccDistance = 60 * displayScale * scale;
  const hvccDir = new THREE.Vector3(0.5, 0, -0.866);
  const hvccPosInGalaxy = hvccDir.multiplyScalar(hvccDistance);
  const hvccIndicator = new THREE.SphereGeometry(0.3 * scale, 16, 16);
  const hvccMat = new THREE.MeshBasicMaterial({
    color: 65416,
    transparent: true,
    opacity: 0.5
  });
  const hvccMarker = new THREE.Mesh(hvccIndicator, hvccMat);
  hvccMarker.position.copy(hvccPosInGalaxy);
  group.add(hvccMarker);
  const hvccLabel = createTextSprite$1("HVCC (~60pc from GC)", { color: "#00ff88", fontSize: 16 });
  hvccLabel.position.copy(hvccPosInGalaxy.clone().add(new THREE.Vector3(1.5 * scale, 0.5 * scale, 0)));
  hvccLabel.scale.set(5, 2, 1);
  group.add(hvccLabel);
  const legendPos = new THREE.Vector3(diskRadius * 0.7, 0.1, diskRadius * 0.5);
  const scaleNote = createTextSprite$1("Scale: Compressed ~300× for visualization", { color: "#888888", fontSize: 12 });
  scaleNote.position.copy(legendPos);
  scaleNote.scale.set(5, 2, 1);
  group.add(scaleNote);
  const distNote1 = createTextSprite$1("• Sun ↔ GC: 8 kpc (~26.7 display units)", { color: "#ffff88", fontSize: 11 });
  distNote1.position.copy(legendPos.clone().add(new THREE.Vector3(0, -1 * scale, 0)));
  distNote1.scale.set(6, 2, 1);
  group.add(distNote1);
  const distNote2 = createTextSprite$1("• HVCC ↔ GC: 60 pc (~0.2 display units)", { color: "#00ff88", fontSize: 11 });
  distNote2.position.copy(legendPos.clone().add(new THREE.Vector3(0, -2 * scale, 0)));
  distNote2.scale.set(6, 2, 1);
  group.add(distNote2);
  const galaxyLabel = createTextSprite$1("Milky Way (schematic)", { color: "#6688aa", fontSize: 20 });
  galaxyLabel.position.set(0, 3 * scale, diskRadius - 5);
  galaxyLabel.scale.set(6, 3, 1);
  group.add(galaxyLabel);
  group.userData = {
    disk,
    bulge,
    rotationSpeed: config.rotationSpeed || 0
  };
  return group;
}
function updateGalaxyRotation(galaxyGroup, deltaTime) {
  const { rotationSpeed } = galaxyGroup.userData || {};
  if (rotationSpeed) {
    galaxyGroup.rotation.y += rotationSpeed * deltaTime;
  }
}
function updateVCircularAnimation(observerGroup, deltaTime) {
  const vCircularGroup = observerGroup.getObjectByName("vCircularIndicator");
  if (!vCircularGroup) return;
  const { orbitRadius, animationSpeed } = vCircularGroup.userData || {};
  if (!orbitRadius) return;
  let angle = (vCircularGroup.userData.animationAngle || 0) + animationSpeed * deltaTime;
  angle = angle % (Math.PI * 2);
  vCircularGroup.userData.animationAngle = angle;
  const orbitMarker = vCircularGroup.getObjectByName("orbitAnimationMarker");
  if (orbitMarker) {
    orbitMarker.position.set(
      orbitRadius * Math.cos(angle),
      0,
      orbitRadius * Math.sin(angle)
    );
  }
  for (let i = 1; i <= 3; i++) {
    const trailMarker = vCircularGroup.getObjectByName(`orbitTrail${i}`);
    if (trailMarker) {
      const trailAngle = angle - i * 0.15;
      trailMarker.position.set(
        orbitRadius * Math.cos(trailAngle),
        0,
        orbitRadius * Math.sin(trailAngle)
      );
    }
  }
}
function updateLSRRotation(solarSystemGroup, gcPosition, rotationSpeed, deltaTime) {
  if (!rotationSpeed) return;
  const angle = rotationSpeed * deltaTime;
  const pos = solarSystemGroup.position.clone().sub(gcPosition);
  const cosA = Math.cos(angle);
  const sinA = Math.sin(angle);
  const newX = pos.x * cosA - pos.z * sinA;
  const newZ = pos.x * sinA + pos.z * cosA;
  solarSystemGroup.position.set(gcPosition.x + newX, solarSystemGroup.position.y, gcPosition.z + newZ);
  solarSystemGroup.rotation.y += angle;
}
function createAllGalacticMarkers(bhPosition, config = DEFAULT_GALACTIC_CONFIG) {
  const group = new THREE.Group();
  const scale = config.displayScale || 1;
  const distanceToGC_pc = config.distanceToGC_pc ?? DEFAULT_GALACTIC_CONFIG.distanceToGC_pc;
  if (config.showGalaxyDisk) {
    const gcDisplayDist = Math.min(distanceToGC_pc, 25) * scale;
    const gcPosition = bhPosition.clone().add(new THREE.Vector3(gcDisplayDist, 0, 0));
    const defaultRotationSpeed = 0.02;
    const galaxyDisk = createGalaxyDisk(gcPosition, scale, {
      showSpiralArms: true,
      opacity: 0.4,
      rotationSpeed: config.galaxyRotationSpeed ?? defaultRotationSpeed
    });
    galaxyDisk.name = "galaxyDisk";
    group.add(galaxyDisk);
    const gcLabel = createTextSprite$1(`GC (~${distanceToGC_pc} pc)`, {
      color: "#ffaa44",
      fontSize: 16
    });
    gcLabel.position.copy(gcPosition.clone().add(new THREE.Vector3(0, 3 * scale, 0)));
    gcLabel.scale.set(4, 2, 1);
    group.add(gcLabel);
  }
  if (config.showSolarSystem) {
    const observerGroup = createObserverDirectionIndicator(
      bhPosition,
      config,
      15 * scale
      // Symbolic display distance for the indicator
    );
    observerGroup.name = "observerIndicator";
    group.add(observerGroup);
  }
  return group;
}
function createObserverDirectionIndicator(origin, config, displayDistance = 15) {
  const group = new THREE.Group();
  group.position.copy(origin);
  group.name = "observerDirectionIndicator";
  const scale = config.displayScale || 1;
  const inc = (config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination) * Math.PI / 180;
  const pa = (config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle) * Math.PI / 180;
  const lsrVelocity = config.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity;
  const losDir = new THREE.Vector3(
    Math.sin(inc) * Math.cos(pa),
    Math.sin(inc) * Math.sin(pa),
    Math.cos(inc)
  ).normalize();
  const losEnd = losDir.clone().multiplyScalar(displayDistance);
  const losArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    losEnd,
    4500223,
    // Blue for observer direction
    {
      headLength: 1.5 * scale,
      headWidth: 0.6 * scale,
      shaftRadius: 0.1 * scale,
      opacity: 0.9
    }
  );
  group.add(losArrow);
  const observerMarker = new THREE.Group();
  const earthGeo = new THREE.SphereGeometry(0.8 * scale, 32, 32);
  const earthMat = new THREE.MeshBasicMaterial({ color: 2263261 });
  const earthMesh = new THREE.Mesh(earthGeo, earthMat);
  observerMarker.add(earthMesh);
  const atmoGeo = new THREE.SphereGeometry(1 * scale, 32, 32);
  const atmoMat = new THREE.MeshBasicMaterial({
    color: 8965375,
    transparent: true,
    opacity: 0.25
  });
  observerMarker.add(new THREE.Mesh(atmoGeo, atmoMat));
  observerMarker.position.copy(losEnd);
  group.add(observerMarker);
  const observerLabel = createTextSprite$1("⊕ Observer (Earth)", { color: "#88ccff", fontSize: 20 });
  observerLabel.position.copy(losEnd.clone().add(new THREE.Vector3(0, 2 * scale, 0)));
  observerLabel.scale.set(5, 2.5, 1);
  group.add(observerLabel);
  const distLabel = createTextSprite$1("~8 kpc (not to scale)", { color: "#aaddff", fontSize: 14 });
  distLabel.position.copy(losEnd.clone().add(new THREE.Vector3(0, 1.2 * scale, 0)));
  distLabel.scale.set(4, 2, 1);
  group.add(distLabel);
  const incArcRadius = displayDistance * 0.3;
  const incArcPoints = [];
  const numArcPoints = 20;
  for (let i = 0; i <= numArcPoints; i++) {
    const t = i / numArcPoints;
    const angle = t * inc;
    const arcPoint = new THREE.Vector3(
      Math.sin(angle) * Math.cos(pa),
      Math.sin(angle) * Math.sin(pa),
      Math.cos(angle)
    ).multiplyScalar(incArcRadius);
    incArcPoints.push(arcPoint);
  }
  const incArcGeo = new THREE.BufferGeometry().setFromPoints(incArcPoints);
  const incArcMat = new THREE.LineBasicMaterial({
    color: 16755268,
    transparent: true,
    opacity: 0.8
  });
  group.add(new THREE.Line(incArcGeo, incArcMat));
  const incLabelPos = new THREE.Vector3(
    Math.sin(inc / 2) * Math.cos(pa),
    Math.sin(inc / 2) * Math.sin(pa),
    Math.cos(inc / 2)
  ).multiplyScalar(incArcRadius + 2 * scale);
  const incLabel = createTextSprite$1(`i = ${config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination}°`, {
    color: "#ffcc66",
    fontSize: 18
  });
  incLabel.position.copy(incLabelPos);
  incLabel.scale.set(3.5, 1.8, 1);
  group.add(incLabel);
  const zAxisLength = displayDistance * 0.5;
  const zAxisEnd = new THREE.Vector3(0, 0, zAxisLength);
  const zAxisArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    zAxisEnd,
    4521864,
    // Green for orbital plane normal
    {
      headLength: 0.8 * scale,
      headWidth: 0.35 * scale,
      shaftRadius: 0.06 * scale,
      opacity: 0.7
    }
  );
  group.add(zAxisArrow);
  const zLabel = createTextSprite$1("L (orbital normal)", { color: "#66ffaa", fontSize: 16 });
  zLabel.position.copy(zAxisEnd.clone().add(new THREE.Vector3(2 * scale, 0.5 * scale, 0)));
  zLabel.scale.set(4, 2, 1);
  group.add(zLabel);
  const orbitalPlaneRadius = displayDistance * 0.25;
  const orbitalPlane = createDashedCircle(orbitalPlaneRadius, 8978312, {
    segments: 32,
    dashSize: 0.5 * scale,
    gapSize: 0.3 * scale,
    opacity: 0.5,
    normal: new THREE.Vector3(0, 0, 1),
    // Normal is Z
    center: new THREE.Vector3(0, 0, 0)
  });
  group.add(orbitalPlane);
  const planeLabel = createTextSprite$1("Orbital Plane", { color: "#aaffaa", fontSize: 14 });
  planeLabel.position.set(orbitalPlaneRadius + 2 * scale, 0, 0);
  planeLabel.scale.set(3.5, 1.8, 1);
  group.add(planeLabel);
  const vlsrLength = Math.min(Math.abs(lsrVelocity) / 20, displayDistance * 0.3);
  const vlsrDir = lsrVelocity < 0 ? losDir.clone().negate() : losDir.clone();
  const vlsrStart = losEnd.clone().add(losDir.clone().multiplyScalar(-2 * scale));
  const vlsrEndPoint = vlsrStart.clone().add(vlsrDir.multiplyScalar(vlsrLength));
  const vlsrColor = lsrVelocity < 0 ? 16738047 : 16755302;
  const vlsrArrow = createArrow(vlsrStart, vlsrEndPoint, vlsrColor, {
    headLength: 0.6 * scale,
    headWidth: 0.25 * scale,
    shaftRadius: 0.05 * scale,
    opacity: 0.8
  });
  group.add(vlsrArrow);
  const vlsrLabel = createTextSprite$1(`V_LSR = ${lsrVelocity} km/s`, {
    color: lsrVelocity < 0 ? "#ff88ff" : "#ffcc88",
    fontSize: 16
  });
  vlsrLabel.position.copy(vlsrEndPoint.clone().add(new THREE.Vector3(1.5 * scale, 0.5 * scale, 0)));
  vlsrLabel.scale.set(4, 2, 1);
  group.add(vlsrLabel);
  const motionLabel = createTextSprite$1(
    lsrVelocity < 0 ? "(approaching)" : "(receding)",
    { color: lsrVelocity < 0 ? "#88ffff" : "#ffddaa", fontSize: 12 }
  );
  motionLabel.position.copy(vlsrLabel.position.clone().add(new THREE.Vector3(0, -0.8 * scale, 0)));
  motionLabel.scale.set(3, 1.5, 1);
  group.add(motionLabel);
  const vCircularGroup = new THREE.Group();
  vCircularGroup.name = "vCircularIndicator";
  vCircularGroup.position.copy(losEnd.clone().add(new THREE.Vector3(0, -4 * scale, 0)));
  const orbitRadius = 3 * scale;
  const orbitCirclePoints = [];
  const numSegments = 64;
  for (let i = 0; i <= numSegments; i++) {
    const angle = i / numSegments * Math.PI * 2;
    orbitCirclePoints.push(new THREE.Vector3(
      orbitRadius * Math.cos(angle),
      0,
      orbitRadius * Math.sin(angle)
    ));
  }
  const orbitCircleGeo = new THREE.BufferGeometry().setFromPoints(orbitCirclePoints);
  const orbitCircleMat = new THREE.LineBasicMaterial({
    color: 8956671,
    transparent: true,
    opacity: 0.4
  });
  const orbitCircle = new THREE.Line(orbitCircleGeo, orbitCircleMat);
  vCircularGroup.add(orbitCircle);
  const gcMarkerGeo = new THREE.SphereGeometry(0.25 * scale, 16, 16);
  const gcMarkerMat = new THREE.MeshBasicMaterial({ color: 16755268 });
  const gcMarker = new THREE.Mesh(gcMarkerGeo, gcMarkerMat);
  gcMarker.position.set(0, 0, 0);
  vCircularGroup.add(gcMarker);
  const gcCenterLabel = createTextSprite$1("Galactic Center", { color: "#ffcc66", fontSize: 12 });
  gcCenterLabel.position.set(0, 0.8 * scale, 0);
  gcCenterLabel.scale.set(3.5, 1.8, 1);
  vCircularGroup.add(gcCenterLabel);
  const sunPosOnOrbit = new THREE.Vector3(orbitRadius, 0, 0);
  const sunMarkerGeo = new THREE.SphereGeometry(0.35 * scale, 16, 16);
  const sunMarkerMat = new THREE.MeshBasicMaterial({ color: 16777028 });
  const sunMarker = new THREE.Mesh(sunMarkerGeo, sunMarkerMat);
  sunMarker.position.copy(sunPosOnOrbit);
  sunMarker.name = "sunOrbitMarker";
  vCircularGroup.add(sunMarker);
  const sunOrbitLabel = createTextSprite$1("☉ Sun (8 kpc from GC)", { color: "#ffff88", fontSize: 11 });
  sunOrbitLabel.position.copy(sunPosOnOrbit.clone().add(new THREE.Vector3(0, 0.7 * scale, 0)));
  sunOrbitLabel.scale.set(4, 2, 1);
  vCircularGroup.add(sunOrbitLabel);
  const vCircularDir = new THREE.Vector3(0, 0, 1);
  const vCircularArrowLength = 2.5 * scale;
  const vCircularArrowEnd = sunPosOnOrbit.clone().add(vCircularDir.clone().multiplyScalar(vCircularArrowLength));
  const vCircularArrow = createArrow(
    sunPosOnOrbit,
    vCircularArrowEnd,
    4521864,
    // Green for V_circular
    {
      headLength: 0.5 * scale,
      headWidth: 0.25 * scale,
      shaftRadius: 0.04 * scale,
      opacity: 0.9
    }
  );
  vCircularGroup.add(vCircularArrow);
  const vCircularLabel = createTextSprite$1("V_circular = 220 km/s", {
    color: "#66ffaa",
    fontSize: 14
  });
  vCircularLabel.position.copy(vCircularArrowEnd.clone().add(new THREE.Vector3(0, 0.6 * scale, 0.5 * scale)));
  vCircularLabel.scale.set(4, 2, 1);
  vCircularGroup.add(vCircularLabel);
  const directionLabel = createTextSprite$1("(Galactic rotation, l = 90°)", {
    color: "#aaffcc",
    fontSize: 10
  });
  directionLabel.position.copy(vCircularLabel.position.clone().add(new THREE.Vector3(0, -0.6 * scale, 0)));
  directionLabel.scale.set(4, 2, 1);
  vCircularGroup.add(directionLabel);
  const orbitMarkerGeo = new THREE.SphereGeometry(0.2 * scale, 12, 12);
  const orbitMarkerMat = new THREE.MeshBasicMaterial({
    color: 8978431,
    transparent: true,
    opacity: 0.9
  });
  const orbitMarker = new THREE.Mesh(orbitMarkerGeo, orbitMarkerMat);
  orbitMarker.name = "orbitAnimationMarker";
  orbitMarker.position.set(orbitRadius, 0, 0);
  vCircularGroup.add(orbitMarker);
  for (let i = 1; i <= 3; i++) {
    const trailGeo = new THREE.SphereGeometry((0.2 - i * 0.04) * scale, 8, 8);
    const trailMat = new THREE.MeshBasicMaterial({
      color: 8978431,
      transparent: true,
      opacity: 0.6 - i * 0.15
    });
    const trailMarker = new THREE.Mesh(trailGeo, trailMat);
    trailMarker.name = `orbitTrail${i}`;
    const trailAngle = -i * 0.15;
    trailMarker.position.set(
      orbitRadius * Math.cos(trailAngle),
      0,
      orbitRadius * Math.sin(trailAngle)
    );
    vCircularGroup.add(trailMarker);
  }
  vCircularGroup.userData = {
    orbitRadius,
    animationAngle: 0,
    // Current animation angle
    animationSpeed: 0.4
    // rad/s for smooth visualization
  };
  const periodLabel = createTextSprite$1("Orbital Period T ≈ 220 Myr", {
    color: "#88aacc",
    fontSize: 11
  });
  periodLabel.position.set(0, -orbitRadius - 1 * scale, 0);
  periodLabel.scale.set(4, 2, 1);
  vCircularGroup.add(periodLabel);
  group.add(vCircularGroup);
  const paIndicatorGroup = new THREE.Group();
  const skyPlaneRadius = displayDistance * 0.15;
  const skyPlane = createDashedCircle(skyPlaneRadius, 16777096, {
    segments: 24,
    dashSize: 0.3 * scale,
    gapSize: 0.2 * scale,
    opacity: 0.4,
    normal: losDir,
    center: new THREE.Vector3(0, 0, 0)
  });
  paIndicatorGroup.add(skyPlane);
  const northInSky = new THREE.Vector3(0, 0, 1);
  const eastInSky = new THREE.Vector3().crossVectors(losDir, northInSky).normalize();
  const northProj = new THREE.Vector3().crossVectors(eastInSky, losDir).normalize();
  const paDirection = northProj.clone().multiplyScalar(Math.cos(pa)).add(eastInSky.clone().multiplyScalar(Math.sin(pa))).normalize();
  const paArrowEnd = paDirection.multiplyScalar(skyPlaneRadius);
  const paArrow = createArrow(
    new THREE.Vector3(0, 0, 0),
    paArrowEnd,
    16777028,
    {
      headLength: 0.4 * scale,
      headWidth: 0.2 * scale,
      shaftRadius: 0.03 * scale,
      opacity: 0.7
    }
  );
  paIndicatorGroup.add(paArrow);
  paIndicatorGroup.position.copy(losEnd.clone().add(losDir.clone().multiplyScalar(-3 * scale)));
  group.add(paIndicatorGroup);
  const paLabel = createTextSprite$1(`PA = ${config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle}°`, {
    color: "#ffff88",
    fontSize: 14
  });
  paLabel.position.copy(paIndicatorGroup.position.clone().add(new THREE.Vector3(skyPlaneRadius + 2 * scale, 0, 0)));
  paLabel.scale.set(3.5, 1.8, 1);
  group.add(paLabel);
  const geometryNote = createTextSprite$1("Observer Geometry (Oka et al. 2017)", {
    color: "#888888",
    fontSize: 12
  });
  geometryNote.position.set(0, -3 * scale, displayDistance * 0.5);
  geometryNote.scale.set(5, 2, 1);
  group.add(geometryNote);
  return group;
}
let circleTexture = null;
function createCircleTexture(size = 64) {
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
function getCircleTexture() {
  if (!circleTexture) {
    circleTexture = createCircleTexture();
  }
  return circleTexture;
}
const DEFAULT_TEXT_SPRITE_OPTIONS = {
  color: "#ffffff",
  fontSize: 72,
  // Increased from 48
  fontWeight: "Bold",
  fontFamily: "Arial, sans-serif",
  backgroundColor: "",
  padding: 10,
  scale: [12, 6, 1],
  // Increased from [8, 4, 1]
  canvasSize: 512,
  renderOrder: 999
};
function createTextSprite(text, options = {}) {
  const opts = { ...DEFAULT_TEXT_SPRITE_OPTIONS, ...options };
  const {
    color,
    fontSize,
    fontWeight,
    fontFamily,
    backgroundColor,
    scale,
    canvasSize,
    renderOrder
  } = opts;
  const canvas = document.createElement("canvas");
  canvas.width = canvasSize;
  canvas.height = canvasSize / 2;
  const ctx = canvas.getContext("2d");
  if (backgroundColor) {
    ctx.fillStyle = backgroundColor;
    ctx.fillRect(0, 0, canvasSize, canvasSize / 2);
  }
  ctx.fillStyle = color;
  ctx.font = `${fontWeight} ${fontSize}px ${fontFamily}`;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, canvasSize / 2, canvasSize / 4);
  const texture = new THREE.CanvasTexture(canvas);
  texture.needsUpdate = true;
  const material = new THREE.SpriteMaterial({
    map: texture,
    transparent: true,
    depthTest: false,
    depthWrite: false
  });
  const sprite = new THREE.Sprite(material);
  sprite.scale.set(...scale);
  sprite.renderOrder = renderOrder;
  return sprite;
}
function createAxisLabel(text, axisColor, position) {
  const sprite = createTextSprite(text, {
    color: axisColor,
    fontSize: 72,
    scale: [12, 6, 1]
    // Larger scale for axis labels
  });
  sprite.position.set(...position);
  return sprite;
}
function createTickLabel(text, position) {
  const sprite = createTextSprite(text, {
    color: "#aaaaaa",
    fontSize: 48,
    scale: [6, 3, 1]
    // Smaller scale for tick labels
  });
  sprite.position.set(...position);
  return sprite;
}
function createAxesWithLabels(options = {}) {
  const {
    size = 30,
    showTickMarks = true,
    tickInterval = 10,
    labelUnit = "pc"
  } = options;
  const group = new THREE.Group();
  group.name = "axesWithLabels";
  const axesHelper = new THREE.AxesHelper(size);
  group.add(axesHelper);
  const labelOffset = size + 4;
  const xLabel = createAxisLabel(`X (${labelUnit})`, "#ff6666", [labelOffset, 0, 0]);
  group.add(xLabel);
  const yLabel = createAxisLabel(`Y (${labelUnit})`, "#66ff66", [0, labelOffset, 0]);
  group.add(yLabel);
  const zLabel = createAxisLabel(`Z (${labelUnit})`, "#6666ff", [0, 0, labelOffset]);
  group.add(zLabel);
  if (showTickMarks) {
    const tickMaterial = new THREE.LineBasicMaterial({ color: 6710886 });
    const tickSize = 0.8;
    for (let i = -Math.floor(size / tickInterval) * tickInterval; i <= size; i += tickInterval) {
      if (i === 0) continue;
      const xTickGeom = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(i, -tickSize, 0),
        new THREE.Vector3(i, tickSize, 0)
      ]);
      group.add(new THREE.Line(xTickGeom, tickMaterial));
      const yTickGeom = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-tickSize, i, 0),
        new THREE.Vector3(tickSize, i, 0)
      ]);
      group.add(new THREE.Line(yTickGeom, tickMaterial));
      if (Math.abs(i) <= size) {
        const xTickLabel = createTickLabel(`${i}`, [i, -3, 0]);
        group.add(xTickLabel);
        const yTickLabel = createTickLabel(`${i}`, [-3, i, 0]);
        group.add(yTickLabel);
      }
    }
  }
  return group;
}
const COLOR_MAP_DATA = {
  // ═══════════════════════════════════════════════════════════════════════════
  // SEQUENTIAL COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════
  /**
   * Cosmic Dawn - Deep space theme
   * Purple → Blue → Cyan → Gold → White
   */
  cosmicDawn: [
    [0.102, 0.02, 0.2],
    // #1a0533
    [0.176, 0.106, 0.412],
    // #2d1b69
    [0.239, 0.31, 0.78],
    // #3d4fc7
    [0, 0.706, 0.847],
    // #00b4d8
    [0.282, 0.792, 0.894],
    // #48cae4
    [0.565, 0.878, 0.937],
    // #90e0ef
    [1, 0.82, 0.4],
    // #ffd166
    [1, 0.922, 0.6],
    // #ffeb99
    [1, 1, 1]
    // #ffffff
  ],
  /**
   * Nebula Fire - Warm energy theme
   * Magenta → Pink → Orange → Yellow → White
   */
  nebulaFire: [
    [0.176, 0.039, 0.192],
    // #2d0a31
    [0.361, 0.067, 0.345],
    // #5c1158
    [0.608, 0.176, 0.498],
    // #9b2d7f
    [0.839, 0.259, 0.573],
    // #d64292
    [0.957, 0.431, 0.431],
    // #f46e6e
    [1, 0.624, 0.263],
    // #ff9f43
    [1, 0.788, 0.235],
    // #ffc93c
    [1, 0.945, 0.463],
    // #fff176
    [1, 1, 1]
    // #ffffff
  ],
  /**
   * Ocean Depths - Cool pressure theme
   * Navy → Blue → Cyan → Mint → Cream
   */
  oceanDepths: [
    [0.051, 0.106, 0.165],
    // #0d1b2a
    [0.106, 0.227, 0.294],
    // #1b3a4b
    [0.078, 0.302, 0.494],
    // #144d7e
    [0.118, 0.533, 0.898],
    // #1e88e5
    [0.259, 0.647, 0.961],
    // #42a5f5
    [0.502, 0.871, 0.918],
    // #80deea
    [0.655, 1, 0.922],
    // #a7ffeb
    [0.878, 0.969, 0.98],
    // #e0f7fa
    [1, 0.992, 0.906]
    // #fffde7
  ],
  /**
   * Aurora - Vivid multi-hue
   * Purple → Blue → Cyan → Green → Yellow
   */
  aurora: [
    [0.227, 0.047, 0.639],
    // #3a0ca3
    [0.263, 0.38, 0.933],
    // #4361ee
    [0.298, 0.788, 0.941],
    // #4cc9f0
    [0.024, 0.839, 0.627],
    // #06d6a0
    [0.322, 0.718, 0.533],
    // #52b788
    [0.6, 0.851, 0.549],
    // #99d98c
    [0.851, 0.929, 0.573],
    // #d9ed92
    [0.988, 0.965, 0.741],
    // #fcf6bd
    [1, 0.953, 0.69]
    // #fff3b0
  ],
  // ═══════════════════════════════════════════════════════════════════════════
  // SCIENTIFIC COLORMAPS (perceptually uniform)
  // ═══════════════════════════════════════════════════════════════════════════
  viridis: [
    [0.267, 4e-3, 0.329],
    // #440154
    [0.282, 0.157, 0.471],
    // #482878
    [0.243, 0.286, 0.537],
    // #3e4989
    [0.192, 0.408, 0.557],
    // #31688e
    [0.149, 0.51, 0.557],
    // #26828e
    [0.122, 0.62, 0.537],
    // #1f9e89
    [0.208, 0.718, 0.475],
    // #35b779
    [0.431, 0.808, 0.345],
    // #6ece58
    [0.71, 0.871, 0.169],
    // #b5de2b
    [0.992, 0.906, 0.145]
    // #fde725
  ],
  plasma: [
    [0.051, 0.031, 0.529],
    // #0d0887
    [0.274, 0.012, 0.624],
    // #46039f
    [0.447, 4e-3, 0.659],
    // #7201a8
    [0.612, 0.09, 0.62],
    // #9c179e
    [0.741, 0.216, 0.525],
    // #bd3786
    [0.847, 0.337, 0.42],
    // #d8576b
    [0.929, 0.475, 0.325],
    // #ed7953
    [0.984, 0.624, 0.227],
    // #fb9f3a
    [0.992, 0.792, 0.149],
    // #fdca26
    [0.941, 0.976, 0.129]
    // #f0f921
  ],
  turbo: [
    [0.188, 0.071, 0.231],
    // #30123b
    [0.275, 0.384, 0.843],
    // #4662d7
    [0.208, 0.682, 0.957],
    // #35aef4
    [0.102, 0.894, 0.714],
    // #1ae4b6
    [0.447, 0.996, 0.369],
    // #72fe5e
    [0.784, 0.937, 0.204],
    // #c8ef34
    [0.98, 0.729, 0.224],
    // #faba39
    [0.965, 0.42, 0.098],
    // #f66b19
    [0.792, 0.165, 0.016],
    // #ca2a04
    [0.478, 0.016, 0.012]
    // #7a0403
  ],
  // ═══════════════════════════════════════════════════════════════════════════
  // DIVERGING COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════
  /**
   * Velocity - Blue to White to Red-Orange
   */
  velocity: [
    [0, 0.467, 0.714],
    // #0077b6
    [0, 0.706, 0.847],
    // #00b4d8
    [0.282, 0.792, 0.894],
    // #48cae4
    [0.565, 0.878, 0.937],
    // #90e0ef
    [0.792, 0.941, 0.973],
    // #caf0f8
    [1, 1, 1],
    // #ffffff
    [1, 0.8, 0.835],
    // #ffccd5
    [1, 0.561, 0.639],
    // #ff8fa3
    [1, 0.373, 0.494],
    // #ff5f7e
    [0.984, 0.38, 0.027],
    // #fb6107
    [0.902, 0.224, 0.275]
    // #e63946
  ],
  /**
   * Divergent Sunset - Purple to Cream to Gold
   */
  divergentSunset: [
    [0.369, 0.165, 0.518],
    // #5e2a84
    [0.545, 0.302, 0.62],
    // #8b4d9e
    [0.71, 0.486, 0.753],
    // #b57cc0
    [0.847, 0.659, 0.847],
    // #d8a8d8
    [0.961, 0.902, 0.91],
    // #f5e6e8
    [1, 0.973, 0.863],
    // #fff8dc
    [1, 0.875, 0.4],
    // #ffe066
    [1, 0.761, 0.2],
    // #ffc233
    [1, 0.624, 0.11],
    // #ff9f1c
    [1, 0.42, 0.208],
    // #ff6b35
    [0.839, 0.157, 0.157]
    // #d62828
  ],
  // ═══════════════════════════════════════════════════════════════════════════
  // SPECIALIZED COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════
  /**
   * Density - Optimized for log-scale density
   */
  density: [
    [0.102, 0.102, 0.18],
    // #1a1a2e
    [0.086, 0.129, 0.243],
    // #16213e
    [0.118, 0.227, 0.373],
    // #1e3a5f
    [0.196, 0.51, 0.722],
    // #3282b8
    [0, 0.737, 0.831],
    // #00bcd4
    [0, 0.898, 1],
    // #00e5ff
    [0.463, 1, 0.012],
    // #76ff03
    [1, 0.922, 0.231],
    // #ffeb3b
    [1, 0.596, 0],
    // #ff9800
    [1, 0.341, 0.133],
    // #ff5722
    [1, 0.09, 0.267]
    // #ff1744
  ],
  /**
   * Black Hole - Dramatic accretion disk theme
   */
  blackHole: [
    [0.059, 0.059, 0.137],
    // #0f0f23
    [0.11, 0.11, 0.302],
    // #1c1c4d
    [0.231, 0.18, 0.353],
    // #3b2e5a
    [0.361, 0.29, 0.447],
    // #5c4a72
    [1, 0.42, 0.208],
    // #ff6b35
    [1, 0.714, 0.153],
    // #ffb627
    [1, 0.914, 0.498],
    // #ffe97f
    [1, 0.976, 0.769],
    // #fff9c4
    [1, 1, 1]
    // #ffffff
  ],
  /**
   * Ice & Fire - Temperature contrast
   */
  iceFire: [
    [0, 0.161, 0.42],
    // #00296b
    [0, 0.247, 0.533],
    // #003f88
    [0, 0.314, 0.616],
    // #00509d
    [0, 0.467, 0.714],
    // #0077b6
    [0, 0.706, 0.847],
    // #00b4d8
    [0.565, 0.878, 0.937],
    // #90e0ef
    [1, 0.761, 0],
    // #ffc300
    [1, 0.584, 0],
    // #ff9500
    [1, 0.404, 0],
    // #ff6700
    [1, 0.239, 0],
    // #ff3d00
    [0.835, 0, 0]
    // #d50000
  ],
  /**
   * Spectral Bright - Enhanced rainbow
   */
  spectralBright: [
    [0.29, 0.078, 0.549],
    // #4a148c
    [0.482, 0.122, 0.635],
    // #7b1fa2
    [0.082, 0.396, 0.753],
    // #1565c0
    [8e-3, 0.533, 0.82],
    // #0288d1
    [0, 0.675, 0.757],
    // #00acc1
    [0, 0.537, 0.482],
    // #00897b
    [0.263, 0.627, 0.278],
    // #43a047
    [0.486, 0.702, 0.259],
    // #7cb342
    [0.753, 0.792, 0.2],
    // #c0ca33
    [0.992, 0.847, 0.208],
    // #fdd835
    [1, 0.702, 0],
    // #ffb300
    [0.984, 0.549, 0],
    // #fb8c00
    [0.957, 0.318, 0.118]
    // #f4511e
  ],
  /**
   * Mono Cyan - Clean single-hue
   */
  monoCyan: [
    [0, 0.071, 0.098],
    // #001219
    [0, 0.373, 0.451],
    // #005f73
    [0.039, 0.576, 0.588],
    // #0a9396
    [0.251, 0.788, 0.788],
    // #40c9c9
    [0.58, 0.824, 0.741],
    // #94d2bd
    [0.914, 0.847, 0.651],
    // #e9d8a6
    [1, 1, 1]
    // #ffffff
  ],
  /**
   * Mono Gold - Warm single-hue
   */
  monoGold: [
    [0.102, 0.102, 0.039],
    // #1a1a0a
    [0.239, 0.239, 0],
    // #3d3d00
    [0.42, 0.42, 0],
    // #6b6b00
    [0.722, 0.525, 0.043],
    // #b8860b
    [0.855, 0.647, 0.125],
    // #daa520
    [1, 0.843, 0],
    // #ffd700
    [1, 0.922, 0.231],
    // #ffeb3b
    [1, 0.961, 0.616],
    // #fff59d
    [1, 1, 1]
    // #ffffff
  ],
  // ═══════════════════════════════════════════════════════════════════════════
  // LEGACY COLORMAPS (kept for backward compatibility)
  // ═══════════════════════════════════════════════════════════════════════════
  inferno: [
    [1e-3, 0, 0.014],
    [0.122, 0.047, 0.282],
    [0.304, 0.063, 0.42],
    [0.499, 0.086, 0.397],
    [0.68, 0.144, 0.295],
    [0.833, 0.253, 0.16],
    [0.937, 0.405, 0.049],
    [0.981, 0.588, 0.068],
    [0.987, 0.772, 0.264],
    [0.988, 0.998, 0.645]
  ],
  magma: [
    [1e-3, 0, 0.014],
    [0.104, 0.047, 0.258],
    [0.259, 0.05, 0.408],
    [0.427, 0.079, 0.43],
    [0.575, 0.134, 0.397],
    [0.716, 0.215, 0.345],
    [0.848, 0.343, 0.331],
    [0.937, 0.517, 0.388],
    [0.973, 0.699, 0.53],
    [0.988, 0.998, 0.645]
  ],
  coolwarm: [
    [0.23, 0.299, 0.754],
    [0.411, 0.484, 0.845],
    [0.593, 0.669, 0.927],
    [0.775, 0.817, 0.964],
    [0.9, 0.9, 0.9],
    [0.964, 0.775, 0.692],
    [0.927, 0.593, 0.476],
    [0.845, 0.411, 0.299],
    [0.754, 0.23, 0.173]
  ]
};
function sampleColorMap(mapName, t) {
  const map = COLOR_MAP_DATA[mapName] || COLOR_MAP_DATA.viridis;
  t = Math.max(0, Math.min(1, t));
  const idx = t * (map.length - 1);
  const i = Math.floor(idx);
  const f = idx - i;
  if (i >= map.length - 1) {
    return map[map.length - 1];
  }
  return [
    map[i][0] + f * (map[i + 1][0] - map[i][0]),
    map[i][1] + f * (map[i + 1][1] - map[i][1]),
    map[i][2] + f * (map[i + 1][2] - map[i][2])
  ];
}
function computeHyperbolicOrbit(config, numPoints = 200) {
  const { bhPosition, bhMass, cloudInitialPosition, cloudInitialVelocity } = config;
  const [x0, y0, z0] = cloudInitialPosition;
  const [vx0, vy0, vz0] = cloudInitialVelocity;
  const [bh_x, bh_y, bh_z] = bhPosition;
  const rx = x0 - bh_x;
  const ry = y0 - bh_y;
  const r0 = Math.sqrt(rx * rx + ry * ry);
  const v0 = Math.sqrt(vx0 * vx0 + vy0 * vy0 + vz0 * vz0);
  const GM = bhMass;
  const epsilon = v0 * v0 / 2 - GM / r0;
  if (epsilon <= 0) {
    console.warn("[computeHyperbolicOrbit] Non-hyperbolic orbit (ε <= 0), using straight line");
    return computeStraightLineTrajectory({ cloudInitialPosition, cloudInitialVelocity }, -0.5, 5, numPoints);
  }
  const h = Math.abs(rx * vy0 - ry * vx0);
  const p = h * h / GM;
  const e = Math.sqrt(1 + 2 * epsilon * h * h / (GM * GM));
  const theta0 = Math.atan2(ry, rx);
  const cosNu0 = Math.max(-0.9999, Math.min(0.9999, (p / r0 - 1) / e));
  const rdotv = rx * vx0 + ry * vy0;
  const nu0 = rdotv < 0 ? -Math.acos(cosNu0) : Math.acos(cosNu0);
  const omega = theta0 - nu0;
  const nu_asymp = Math.acos(-1 / e);
  const nu_max = nu_asymp - 0.05;
  const nu_min = -nu_max;
  const orbit = [];
  for (let i = 0; i < numPoints; i++) {
    const nu = nu_min + (nu_max - nu_min) * (i / (numPoints - 1));
    const denom = 1 + e * Math.cos(nu);
    if (denom <= 0) continue;
    const r = p / denom;
    const theta = omega + nu;
    const x = bh_x + r * Math.cos(theta);
    const y = bh_y + r * Math.sin(theta);
    const z = bh_z;
    if (Math.abs(x) < 100 && Math.abs(y) < 100 && r > 0) {
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
  showRadii = true,
  showGalacticMarkers = true,
  showLabels = true,
  galacticConfig,
  cameraMode = "free",
  onCameraModeChange
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
  const galacticMarkersGroupRef = useRef(null);
  const lastFrameIndexRef = useRef(-1);
  const lastColorFieldRef = useRef(colorField);
  const lastColorMapRef = useRef(colorMapName);
  const lastLogScaleRef = useRef(logScale);
  const colorFieldRef = useRef(colorField);
  const colorMapNameRef = useRef(colorMapName);
  const logScaleRef = useRef(logScale);
  const showTrajectoryRef = useRef(showTrajectory);
  const imbhPhysicsRef = useRef(imbhPhysics);
  const globalColorRangeRef = useRef(globalColorRange);
  const updateParticlesRef = useRef(() => {
  });
  const updateTrajectoryRef = useRef(() => {
  });
  const computeStatsRef = useRef(() => {
  });
  useEffect(() => {
    colorFieldRef.current = colorField;
    if (lastColorFieldRef.current !== colorField) {
      lastColorFieldRef.current = "__force_update__";
    }
  }, [colorField]);
  useEffect(() => {
    colorMapNameRef.current = colorMapName;
    lastColorMapRef.current = "__force_update__";
  }, [colorMapName]);
  useEffect(() => {
    logScaleRef.current = logScale;
    if (lastLogScaleRef.current !== logScale) {
      lastLogScaleRef.current = !logScale;
    }
  }, [logScale]);
  useEffect(() => {
    showTrajectoryRef.current = showTrajectory;
  }, [showTrajectory]);
  useEffect(() => {
    imbhPhysicsRef.current = imbhPhysics;
  }, [imbhPhysics]);
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
    const currentColorField = colorFieldRef.current;
    const currentColorMapName = colorMapNameRef.current;
    const currentLogScale = logScaleRef.current;
    const needsUpdate = frameIndex !== lastFrameIndexRef.current || currentColorField !== lastColorFieldRef.current || currentColorMapName !== lastColorMapRef.current || currentLogScale !== lastLogScaleRef.current;
    if (!needsUpdate) return;
    lastFrameIndexRef.current = frameIndex;
    lastColorFieldRef.current = currentColorField;
    lastColorMapRef.current = currentColorMapName;
    lastLogScaleRef.current = currentLogScale;
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
    switch (currentColorField) {
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
    const logMin = currentLogScale && vMin > 0 ? Math.log10(vMin) : 0;
    const logRange = currentLogScale && vMin > 0 ? Math.log10(vMax) - logMin : 1;
    const range = vMax - vMin || 1;
    for (let i = 0; i < frame.particleCount; i++) {
      let val = fieldData[i];
      if (!isFinite(val)) val = vMin;
      let t;
      if (currentLogScale && vMin > 0) {
        t = (Math.log10(Math.max(val, vMin)) - logMin) / logRange;
      } else {
        t = (val - vMin) / range;
      }
      t = Math.max(0, Math.min(1, t));
      const [r, g, b] = sampleColorMap(currentColorMapName, t);
      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;
    }
    geometry.attributes.color.needsUpdate = true;
    geometry.computeBoundingSphere();
  }, [framesRef, frameIndexRef]);
  const updateTrajectory = useCallback(() => {
    const currentImbhPhysics = imbhPhysicsRef.current;
    if (!currentImbhPhysics?.enabled) return;
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
    const currentShowTrajectory = showTrajectoryRef.current;
    if (comMarkerRef.current) {
      comMarkerRef.current.position.set(comX, comY, comZ);
      comMarkerRef.current.visible = currentShowTrajectory;
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
  }, [framesRef, frameIndexRef]);
  useEffect(() => {
    updateParticlesRef.current = updateParticles;
  }, [updateParticles]);
  useEffect(() => {
    updateTrajectoryRef.current = updateTrajectory;
  }, [updateTrajectory]);
  useEffect(() => {
    computeStatsRef.current = computeStats;
  }, [computeStats]);
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
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.screenSpacePanning = true;
    controls.minDistance = 1;
    controls.maxDistance = 5e3;
    controlsRef.current = controls;
    const gridHelper = new THREE.GridHelper(200, 20, 4473924, 2236962);
    scene.add(gridHelper);
    if (showAxes) {
      const axesGroup = createAxesWithLabels({
        size: 30,
        showTickMarks: true,
        tickInterval: 10,
        labelUnit: "pc"
      });
      scene.add(axesGroup);
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
    if (showGalacticMarkers && imbhPhysics?.enabled) {
      const bhPos = new THREE.Vector3(
        imbhPhysics.bhPosition[0],
        imbhPhysics.bhPosition[1],
        imbhPhysics.bhPosition[2]
      );
      const galConfig = {
        distanceToGC: galacticConfig?.distanceToGC ?? DEFAULT_GALACTIC_CONFIG.distanceToGC,
        galacticLongitude: galacticConfig?.galacticLongitude ?? DEFAULT_GALACTIC_CONFIG.galacticLongitude,
        galacticLatitude: galacticConfig?.galacticLatitude ?? DEFAULT_GALACTIC_CONFIG.galacticLatitude,
        inclination: galacticConfig?.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination,
        positionAngle: galacticConfig?.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle,
        lsrVelocity: galacticConfig?.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity,
        hvccPosition: DEFAULT_GALACTIC_CONFIG.hvccPosition,
        hvccSize: DEFAULT_GALACTIC_CONFIG.hvccSize,
        // Pass cloud radius from imbhPhysics if available
        cloudRadius_pc: imbhPhysics?.cloudRadius ?? DEFAULT_GALACTIC_CONFIG.cloudRadius_pc,
        // Galaxy disk visualization options
        showGalaxyDisk: galacticConfig?.showGalaxyDisk ?? false,
        // Off by default for performance
        showSolarSystem: galacticConfig?.showSolarSystem ?? false,
        // Off by default
        galaxyRotationSpeed: galacticConfig?.galaxyRotationSpeed ?? 0.1
        // Slow rotation for demo
      };
      const galacticMarkers = createAllGalacticMarkers(bhPos, galConfig);
      galacticMarkers.visible = showGalacticMarkers;
      scene.add(galacticMarkers);
      galacticMarkersGroupRef.current = galacticMarkers;
    }
    computeStatsRef.current();
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
    let lastTime = performance.now();
    const animate = () => {
      animationId = requestAnimationFrame(animate);
      const now = performance.now();
      const deltaTime = (now - lastTime) / 1e3;
      lastTime = now;
      fpsRef.current.frames++;
      if (now - fpsRef.current.lastTime >= 1e3) {
        fpsRef.current.fps = fpsRef.current.frames;
        fpsRef.current.frames = 0;
        fpsRef.current.lastTime = now;
        onFpsUpdate?.(fpsRef.current.fps);
      }
      updateParticlesRef.current();
      updateTrajectoryRef.current();
      if (galacticMarkersGroupRef.current) {
        const galaxyDisk = galacticMarkersGroupRef.current.getObjectByName("galaxyDisk");
        if (galaxyDisk) {
          updateGalaxyRotation(galaxyDisk, deltaTime);
          const solarSystemDemo = galacticMarkersGroupRef.current.getObjectByName("solarSystemDemo");
          if (solarSystemDemo && galaxyDisk.userData?.rotationSpeed) {
            const gcPosition = galaxyDisk.position.clone();
            updateLSRRotation(
              solarSystemDemo,
              gcPosition,
              galaxyDisk.userData.rotationSpeed,
              deltaTime
            );
          }
        }
        const observerIndicator = galacticMarkersGroupRef.current.getObjectByName("observerIndicator");
        if (observerIndicator) {
          updateVCircularAnimation(observerIndicator, deltaTime);
        }
      }
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
    if (galacticMarkersGroupRef.current) {
      galacticMarkersGroupRef.current.visible = showGalacticMarkers;
    }
    if (sceneRef.current) {
      sceneRef.current.traverse((object) => {
        if (object instanceof THREE.Sprite) {
          object.visible = showLabels;
        }
      });
    }
  }, [showBlackHole, showTrajectory, showRadii, showGalacticMarkers, showLabels]);
  useEffect(() => {
    if (!sceneRef.current || !imbhPhysics?.enabled || !showGalacticMarkers) return;
    if (galacticMarkersGroupRef.current) {
      sceneRef.current.remove(galacticMarkersGroupRef.current);
      galacticMarkersGroupRef.current = null;
    }
    const bhPos = new THREE.Vector3(
      imbhPhysics.bhPosition[0],
      imbhPhysics.bhPosition[1],
      imbhPhysics.bhPosition[2]
    );
    const galConfig = {
      distanceToGC: galacticConfig?.distanceToGC ?? DEFAULT_GALACTIC_CONFIG.distanceToGC,
      galacticLongitude: galacticConfig?.galacticLongitude ?? DEFAULT_GALACTIC_CONFIG.galacticLongitude,
      galacticLatitude: galacticConfig?.galacticLatitude ?? DEFAULT_GALACTIC_CONFIG.galacticLatitude,
      inclination: galacticConfig?.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination,
      positionAngle: galacticConfig?.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle,
      lsrVelocity: galacticConfig?.lsrVelocity ?? DEFAULT_GALACTIC_CONFIG.lsrVelocity,
      hvccPosition: DEFAULT_GALACTIC_CONFIG.hvccPosition,
      hvccSize: DEFAULT_GALACTIC_CONFIG.hvccSize,
      // Pass cloud radius from imbhPhysics if available
      cloudRadius_pc: imbhPhysics?.cloudRadius ?? DEFAULT_GALACTIC_CONFIG.cloudRadius_pc,
      showGalaxyDisk: galacticConfig?.showGalaxyDisk ?? false,
      showSolarSystem: galacticConfig?.showSolarSystem ?? false,
      galaxyRotationSpeed: galacticConfig?.galaxyRotationSpeed ?? 0.1
    };
    const galacticMarkers = createAllGalacticMarkers(bhPos, galConfig);
    galacticMarkers.visible = showGalacticMarkers;
    sceneRef.current.add(galacticMarkers);
    galacticMarkersGroupRef.current = galacticMarkers;
  }, [galacticConfig?.showGalaxyDisk, galacticConfig?.showSolarSystem, galacticConfig?.galaxyRotationSpeed, imbhPhysics, showGalacticMarkers]);
  useEffect(() => {
    computeStatsRef.current();
  }, [framesRef.current?.size]);
  useEffect(() => {
    if (!cameraRef.current || !controlsRef.current || !imbhPhysics?.enabled) return;
    const camera = cameraRef.current;
    const controls = controlsRef.current;
    if (cameraMode === "earth") {
      const inclination = galacticConfig?.inclination ?? 70;
      const positionAngle = galacticConfig?.positionAngle ?? 41.6;
      const inc = inclination * Math.PI / 180;
      const pa = positionAngle * Math.PI / 180;
      const losDir = new THREE.Vector3(
        Math.sin(inc) * Math.cos(pa),
        Math.sin(inc) * Math.sin(pa),
        Math.cos(inc)
      ).normalize();
      const earthDisplayDist = 50;
      const bhPos = new THREE.Vector3(
        imbhPhysics.bhPosition[0],
        imbhPhysics.bhPosition[1],
        imbhPhysics.bhPosition[2]
      );
      const earthPos = bhPos.clone().add(losDir.clone().multiplyScalar(earthDisplayDist));
      camera.position.copy(earthPos);
      controls.target.copy(bhPos);
      camera.lookAt(bhPos);
      controls.enableDamping = true;
      controls.dampingFactor = 0.02;
      console.log(`[Camera] Earth view: i=${inclination}°, PA=${positionAngle}°, pos=(${earthPos.x.toFixed(1)}, ${earthPos.y.toFixed(1)}, ${earthPos.z.toFixed(1)})`);
    } else {
      controls.enableDamping = true;
      controls.dampingFactor = 0.05;
      console.log("[Camera] Free orbit mode");
    }
    controls.update();
  }, [cameraMode, galacticConfig?.inclination, galacticConfig?.positionAngle, imbhPhysics]);
  return /* @__PURE__ */ jsx("div", { ref: containerRef, className: `w-full h-full ${className}` });
}
const Vector3Schema = Schema.Tuple(Schema.Number, Schema.Number, Schema.Number);
const BoundingBoxSchema = Schema.Struct({
  min: Vector3Schema,
  max: Vector3Schema
});
const PhysicalUnitsSchema = Schema.Struct({
  mass: Schema.Number,
  length: Schema.Number,
  time: Schema.Number,
  velocity: Schema.Number,
  density: Schema.Number,
  energy: Schema.Number,
  pressure: Schema.Number
});
const IMBHPhysicsConfigSchema = Schema.Struct({
  enabled: Schema.Boolean,
  bhPosition: Vector3Schema,
  bhMass: Schema.Number,
  cloudInitialPosition: Vector3Schema,
  cloudInitialVelocity: Vector3Schema,
  cloudMass: Schema.Number,
  cloudRadius: Schema.Number,
  tidalRadius: Schema.Number,
  impactParameter: Schema.Number,
  pericentre: Schema.Number,
  eccentricity: Schema.Number,
  timeUnit: Schema.Number,
  inclination: Schema.optional(Schema.Number),
  positionAngle: Schema.optional(Schema.Number),
  lsrVelocity: Schema.optional(Schema.Number)
});
const SPHMethodSchema = Schema.Union(
  Schema.Literal("GSPH"),
  Schema.Literal("SSPH"),
  Schema.Literal("DISPH"),
  Schema.Literal("GDISPH"),
  Schema.Literal("SRGSPH"),
  Schema.String
);
const SPHKernelSchema = Schema.Union(
  Schema.Literal("CubicSpline"),
  Schema.Literal("WendlandC2"),
  Schema.Literal("WendlandC4"),
  Schema.String
);
const DimensionsSchema = Schema.Union(
  Schema.Literal(1),
  Schema.Literal(2),
  Schema.Literal(3)
);
const FieldOffsetsSchema = Schema.Record({
  key: Schema.String,
  value: Schema.Number
});
const SimulationMetadataSchema = Schema.Struct({
  id: Schema.String,
  name: Schema.String,
  description: Schema.String,
  method: SPHMethodSchema,
  kernel: SPHKernelSchema,
  dimensions: DimensionsSchema,
  totalFrames: Schema.Number,
  particleCount: Schema.Number,
  timeRange: Schema.Tuple(Schema.Number, Schema.Number),
  boundingBox: BoundingBoxSchema,
  units: Schema.optional(PhysicalUnitsSchema),
  imbhPhysics: Schema.optional(IMBHPhysicsConfigSchema),
  configPath: Schema.optional(Schema.String),
  dataPath: Schema.String,
  createdAt: Schema.String
});
Schema.Struct({
  frameIndex: Schema.Number,
  time: Schema.Number,
  data: Schema.String,
  stride: Schema.Number,
  fieldOffsets: FieldOffsetsSchema,
  particleCount: Schema.Number
});
const FrameStatisticsSchema = Schema.Struct({
  frameIndex: Schema.Number,
  time: Schema.Number,
  totalMass: Schema.Number,
  totalKineticEnergy: Schema.Number,
  totalInternalEnergy: Schema.Number,
  totalEnergy: Schema.Number,
  momentum: Vector3Schema,
  centerOfMass: Vector3Schema,
  densityRange: Schema.Tuple(Schema.Number, Schema.Number),
  pressureRange: Schema.Tuple(Schema.Number, Schema.Number),
  temperatureRange: Schema.optional(Schema.Tuple(Schema.Number, Schema.Number)),
  maxMach: Schema.optional(Schema.Number),
  particlesInShock: Schema.optional(Schema.Number)
});
const SimulationsListResponseSchema = Schema.Struct({
  simulations: Schema.Array(SimulationMetadataSchema)
});
Schema.Struct({
  frameIndex: Schema.Number,
  time: Schema.Number,
  data: Schema.String,
  stride: Schema.Number,
  fieldOffsets: FieldOffsetsSchema,
  particleCount: Schema.Number
});
Schema.Struct({
  frames: Schema.Array(FrameStatisticsSchema)
});
const ColorMapSchema = Schema.Struct({
  name: Schema.String,
  colors: Schema.Array(Schema.String),
  min: Schema.Number,
  max: Schema.Number,
  logScale: Schema.Boolean
});
Schema.Struct({
  name: Schema.String,
  colors: Schema.Array(Schema.String),
  logScale: Schema.Boolean
});
Schema.Struct({
  currentFrame: Schema.Number,
  isPlaying: Schema.Boolean,
  playbackSpeed: Schema.Number,
  colorField: Schema.String,
  colorMap: ColorMapSchema,
  pointSize: Schema.Number,
  showAxes: Schema.Boolean,
  showBoundingBox: Schema.Boolean,
  cameraPosition: Vector3Schema,
  cameraTarget: Vector3Schema
});
Schema.decodeUnknownEither(SimulationsListResponseSchema);
Schema.decodeUnknownEither(SimulationMetadataSchema);
const COLOR_MAPS = {
  // ═══════════════════════════════════════════════════════════════════════════
  // SEQUENTIAL COLORMAPS (for continuous data like density, velocity magnitude)
  // ═══════════════════════════════════════════════════════════════════════════
  /** 
   * Cosmic Dawn - Deep space theme, excellent for density
   * Gradient: Deep purple → Electric blue → Cyan → Gold → White
   * High contrast on dark backgrounds, avoids muddy mid-tones
   */
  cosmicDawn: {
    name: "Cosmic Dawn",
    colors: [
      "#1a0533",
      // Deep purple-black (visible, not pure black)
      "#2d1b69",
      // Rich purple
      "#3d4fc7",
      // Royal blue
      "#00b4d8",
      // Bright cyan (high visibility)
      "#48cae4",
      // Light cyan
      "#90e0ef",
      // Pale cyan
      "#ffd166",
      // Warm gold
      "#ffeb99",
      // Light gold
      "#ffffff"
      // Pure white (maximum values)
    ],
    logScale: false
  },
  /**
   * Nebula Fire - Warm colormap for energy/temperature
   * Gradient: Deep magenta → Hot pink → Orange → Yellow → White
   * Avoids pure red which is hard to see on black
   */
  nebulaFire: {
    name: "Nebula Fire",
    colors: [
      "#2d0a31",
      // Deep magenta-black
      "#5c1158",
      // Dark magenta
      "#9b2d7f",
      // Bright magenta (visible on dark)
      "#d64292",
      // Hot pink
      "#f46e6e",
      // Coral (not pure red)
      "#ff9f43",
      // Bright orange
      "#ffc93c",
      // Golden yellow
      "#fff176",
      // Light yellow
      "#ffffff"
      // White hot
    ],
    logScale: false
  },
  /**
   * Ocean Depths - Cool colormap for pressure/potential
   * Gradient: Deep teal → Electric blue → Aqua → Mint → Cream
   * Excellent contrast, no muddy greens
   */
  oceanDepths: {
    name: "Ocean Depths",
    colors: [
      "#0d1b2a",
      // Deep navy (visible, not black)
      "#1b3a4b",
      // Dark teal
      "#144d7e",
      // Ocean blue
      "#1e88e5",
      // Bright blue (high visibility)
      "#42a5f5",
      // Sky blue
      "#80deea",
      // Cyan
      "#a7ffeb",
      // Aqua mint
      "#e0f7fa",
      // Pale cyan
      "#fffde7"
      // Warm cream
    ],
    logScale: false
  },
  /**
   * Aurora - Vivid multi-hue for general purpose
   * Gradient: Purple → Blue → Teal → Green → Yellow
   * Perceptually uniform, colorblind-safe
   */
  aurora: {
    name: "Aurora",
    colors: [
      "#3a0ca3",
      // Deep violet
      "#4361ee",
      // Electric blue
      "#4cc9f0",
      // Bright cyan
      "#06d6a0",
      // Teal-green (not pure green)
      "#52b788",
      // Forest green
      "#99d98c",
      // Light green
      "#d9ed92",
      // Yellow-green
      "#fcf6bd",
      // Pale yellow
      "#fff3b0"
      // Cream
    ],
    logScale: false
  },
  // ═══════════════════════════════════════════════════════════════════════════
  // SCIENTIFIC COLORMAPS (perceptually uniform)
  // ═══════════════════════════════════════════════════════════════════════════
  /**
   * Viridis - Scientific standard, perceptually uniform
   * Optimized version with brighter endpoints for dark backgrounds
   */
  viridis: {
    name: "Viridis",
    colors: [
      "#440154",
      // Deep purple
      "#482878",
      // Purple
      "#3e4989",
      // Blue-purple
      "#31688e",
      // Steel blue
      "#26828e",
      // Teal
      "#1f9e89",
      // Teal-green
      "#35b779",
      // Green
      "#6ece58",
      // Yellow-green
      "#b5de2b",
      // Lime
      "#fde725"
      // Bright yellow
    ],
    logScale: false
  },
  /**
   * Plasma - High-energy scientific colormap
   * Better for dark backgrounds than inferno
   */
  plasma: {
    name: "Plasma",
    colors: [
      "#0d0887",
      // Deep blue
      "#46039f",
      // Purple
      "#7201a8",
      // Magenta
      "#9c179e",
      // Pink-purple
      "#bd3786",
      // Hot pink
      "#d8576b",
      // Coral
      "#ed7953",
      // Orange
      "#fb9f3a",
      // Gold
      "#fdca26",
      // Yellow
      "#f0f921"
      // Bright yellow
    ],
    logScale: false
  },
  /**
   * Turbo - Rainbow without the problems
   * High-contrast, colorblind-friendly rainbow alternative
   */
  turbo: {
    name: "Turbo",
    colors: [
      "#30123b",
      // Deep indigo
      "#4662d7",
      // Blue
      "#35aef4",
      // Cyan
      "#1ae4b6",
      // Teal
      "#72fe5e",
      // Green
      "#c8ef34",
      // Yellow-green
      "#faba39",
      // Orange
      "#f66b19",
      // Red-orange
      "#ca2a04",
      // Red (dark enough to see)
      "#7a0403"
      // Dark red
    ],
    logScale: false
  },
  // ═══════════════════════════════════════════════════════════════════════════
  // DIVERGING COLORMAPS (for data with meaningful center point)
  // ═══════════════════════════════════════════════════════════════════════════
  /**
   * Velocity - Blue to White to Red-Orange
   * For velocity divergence (negative/positive)
   * Center is bright white for clear midpoint
   */
  velocity: {
    name: "Velocity",
    colors: [
      "#0077b6",
      // Deep blue
      "#00b4d8",
      // Bright cyan
      "#48cae4",
      // Light cyan
      "#90e0ef",
      // Pale cyan
      "#caf0f8",
      // Very pale cyan
      "#ffffff",
      // White (center)
      "#ffccd5",
      // Pale pink
      "#ff8fa3",
      // Pink
      "#ff5f7e",
      // Coral
      "#fb6107",
      // Orange
      "#e63946"
      // Warm red
    ],
    logScale: false
  },
  /**
   * Divergent Sunset - Purple to Cream to Gold
   * Alternative diverging colormap, elegant
   */
  divergentSunset: {
    name: "Divergent Sunset",
    colors: [
      "#5e2a84",
      // Deep purple
      "#8b4d9e",
      // Purple
      "#b57cc0",
      // Light purple
      "#d8a8d8",
      // Pale purple
      "#f5e6e8",
      // Cream
      "#fff8dc",
      // Light cream (center)
      "#ffe066",
      // Yellow
      "#ffc233",
      // Gold
      "#ff9f1c",
      // Orange
      "#ff6b35",
      // Red-orange
      "#d62828"
      // Deep red
    ],
    logScale: false
  },
  // ═══════════════════════════════════════════════════════════════════════════
  // SPECIALIZED COLORMAPS
  // ═══════════════════════════════════════════════════════════════════════════
  /**
   * Density - Optimized for log-scale density visualization
   * Avoids pure black (invisible) and pure red (hard on dark bg)
   */
  density: {
    name: "Density",
    colors: [
      "#1a1a2e",
      // Very dark blue (visible, not black)
      "#16213e",
      // Navy
      "#1e3a5f",
      // Dark blue
      "#3282b8",
      // Medium blue
      "#00bcd4",
      // Cyan
      "#00e5ff",
      // Bright cyan
      "#76ff03",
      // Neon green
      "#ffeb3b",
      // Yellow
      "#ff9800",
      // Orange
      "#ff5722",
      // Deep orange (not red)
      "#ff1744"
      // Red-pink (brighter than pure red)
    ],
    logScale: true
  },
  /**
   * Black Hole - Dramatic colormap for IMBH simulations
   * Event horizon to accretion disk aesthetic
   */
  blackHole: {
    name: "Black Hole",
    colors: [
      "#0f0f23",
      // Near-black blue
      "#1c1c4d",
      // Dark purple
      "#3b2e5a",
      // Purple
      "#5c4a72",
      // Dusty purple
      "#ff6b35",
      // Orange (accretion disk inner)
      "#ffb627",
      // Gold
      "#ffe97f",
      // Pale gold
      "#fff9c4",
      // Cream
      "#ffffff"
      // White (hottest)
    ],
    logScale: true
  },
  /**
   * Ice & Fire - Dramatic contrast for temperature
   * Cool ice tones to hot fire tones
   */
  iceFire: {
    name: "Ice & Fire",
    colors: [
      "#00296b",
      // Deep blue
      "#003f88",
      // Royal blue
      "#00509d",
      // Bright blue
      "#0077b6",
      // Ocean blue
      "#00b4d8",
      // Cyan
      "#90e0ef",
      // Ice blue
      "#ffc300",
      // Gold (transition)
      "#ff9500",
      // Orange
      "#ff6700",
      // Deep orange
      "#ff3d00",
      // Red-orange
      "#d50000"
      // Deep red (darker, visible)
    ],
    logScale: false
  },
  /**
   * Spectral Bright - Rainbow with enhanced brightness for dark bg
   * Each color chosen for maximum visibility on near-black
   */
  spectralBright: {
    name: "Spectral Bright",
    colors: [
      "#4a148c",
      // Deep purple
      "#7b1fa2",
      // Purple
      "#1565c0",
      // Blue
      "#0288d1",
      // Light blue
      "#00acc1",
      // Cyan
      "#00897b",
      // Teal
      "#43a047",
      // Green
      "#7cb342",
      // Light green
      "#c0ca33",
      // Lime
      "#fdd835",
      // Yellow
      "#ffb300",
      // Amber
      "#fb8c00",
      // Orange
      "#f4511e"
      // Deep orange
    ],
    logScale: false
  },
  /**
   * Mono Cyan - Single-hue gradient for clean scientific viz
   * High contrast, works well with additive blending
   */
  monoCyan: {
    name: "Mono Cyan",
    colors: [
      "#001219",
      // Near-black teal
      "#005f73",
      // Dark teal
      "#0a9396",
      // Teal
      "#40c9c9",
      // Cyan
      "#94d2bd",
      // Pale cyan-green
      "#e9d8a6",
      // Cream
      "#ffffff"
      // White
    ],
    logScale: false
  },
  /**
   * Mono Gold - Warm single-hue for energy visualization
   * Elegant, high contrast on dark backgrounds
   */
  monoGold: {
    name: "Mono Gold",
    colors: [
      "#1a1a0a",
      // Very dark olive
      "#3d3d00",
      // Dark gold
      "#6b6b00",
      // Olive
      "#b8860b",
      // Dark goldenrod
      "#daa520",
      // Goldenrod
      "#ffd700",
      // Gold
      "#ffeb3b",
      // Yellow
      "#fff59d",
      // Pale yellow
      "#ffffff"
      // White
    ],
    logScale: false
  }
};
const CAMERA_POSITIONS = {
  xy: { position: new THREE.Vector3(0, 0, 50), up: new THREE.Vector3(0, 1, 0) },
  // Looking down Z axis at XY plane
  xz: { position: new THREE.Vector3(0, -50, 0), up: new THREE.Vector3(0, 0, 1) },
  // Looking from -Y at XZ plane
  yz: { position: new THREE.Vector3(50, 0, 0), up: new THREE.Vector3(0, 0, 1) }
  // Looking from +X at YZ plane
};
const defaultColorMap = {
  name: "Cosmic Dawn",
  colors: COLOR_MAPS.cosmicDawn?.colors || [
    "#1a0533",
    "#2d1b69",
    "#3d4fc7",
    "#00b4d8",
    "#48cae4",
    "#90e0ef",
    "#ffd166",
    "#ffeb99",
    "#ffffff"
  ],
  min: 0,
  max: 1,
  logScale: false
};
const Projection3DInteractive = forwardRef(({
  framesRef,
  frameIndexRef,
  projection,
  colorField = "density",
  colorMap = defaultColorMap,
  width = 400,
  height = 400,
  className = "",
  showShockSampling = false,
  shockSamplingParams = { columnRadius: 0.15, sliceThickness: 0.15 },
  globalColorRange,
  logScale = false,
  particleSize = 1.5
}, ref) => {
  const containerRef = useRef(null);
  const rendererRef = useRef(null);
  const sceneRef = useRef(null);
  const cameraRef = useRef(null);
  const controlsRef = useRef(null);
  const particlesRef = useRef(null);
  const geometryRef = useRef(null);
  const materialRef = useRef(null);
  const shockGroupRef = useRef(null);
  const lastFrameIndexRef = useRef(-1);
  const lastColorFieldRef = useRef(colorField);
  const lastLogScaleRef = useRef(logScale);
  const animationFrameRef = useRef(null);
  const updateParticlesRef = useRef(() => {
  });
  const computeStatsRef = useRef(() => {
  });
  const colorFieldRef = useRef(colorField);
  const colorMapRef = useRef(colorMap);
  const logScaleRef = useRef(logScale);
  const showShockSamplingRef = useRef(showShockSampling);
  const shockSamplingParamsRef = useRef(shockSamplingParams);
  const globalColorRangeRef = useRef(globalColorRange);
  useEffect(() => {
    colorFieldRef.current = colorField;
    if (lastColorFieldRef.current !== colorField) {
      lastColorFieldRef.current = "__force_update__";
    }
  }, [colorField]);
  useEffect(() => {
    colorMapRef.current = colorMap;
    lastColorFieldRef.current = "__force_update__";
  }, [colorMap]);
  useEffect(() => {
    logScaleRef.current = logScale;
    if (lastLogScaleRef.current !== logScale) {
      lastLogScaleRef.current = !logScale;
    }
  }, [logScale]);
  useEffect(() => {
    showShockSamplingRef.current = showShockSampling;
  }, [showShockSampling]);
  useEffect(() => {
    shockSamplingParamsRef.current = shockSamplingParams;
  }, [shockSamplingParams]);
  useEffect(() => {
    globalColorRangeRef.current = globalColorRange;
  }, [globalColorRange]);
  const statsRef = useRef({ density: [0, 1], velocity: [0, 1], pressure: [0, 1], energy: [0, 1] });
  const [isRotated, setIsRotated] = useState(false);
  const resetCamera = useCallback(() => {
    const camera = cameraRef.current;
    const controls = controlsRef.current;
    if (!camera || !controls) return;
    const config = CAMERA_POSITIONS[projection] || CAMERA_POSITIONS.xy;
    camera.position.copy(config.position);
    camera.up.copy(config.up);
    camera.lookAt(0, 0, 0);
    controls.target.set(0, 0, 0);
    controls.update();
    setIsRotated(false);
  }, [projection]);
  useImperativeHandle(ref, () => ({
    resetCamera
  }), [resetCamera]);
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
  const computeCoM2 = useCallback((frame) => {
    let totalMass = 0;
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < frame.particleCount; i++) {
      const m = frame.mass?.[i] ?? 1;
      cx += m * frame.positions[i * 3];
      cy += m * frame.positions[i * 3 + 1];
      cz += m * frame.positions[i * 3 + 2];
      totalMass += m;
    }
    if (totalMass > 0) {
      cx /= totalMass;
      cy /= totalMass;
      cz /= totalMass;
    }
    return new THREE.Vector3(cx, cy, cz);
  }, []);
  const updateShockSampling = useCallback((frame) => {
    if (!shockGroupRef.current || !showShockSamplingRef.current) return;
    while (shockGroupRef.current.children.length > 0) {
      const child = shockGroupRef.current.children[0];
      shockGroupRef.current.remove(child);
      if (child instanceof THREE.Mesh || child instanceof THREE.Line) {
        child.geometry?.dispose();
        if (child.material instanceof THREE.Material) {
          child.material.dispose();
        }
      }
    }
    const com = computeCoM2(frame);
    const { columnRadius, sliceThickness } = shockSamplingParamsRef.current;
    const comGeom = new THREE.SphereGeometry(0.08, 16, 16);
    const comMat = new THREE.MeshBasicMaterial({ color: 16777215 });
    const comMesh = new THREE.Mesh(comGeom, comMat);
    comMesh.position.copy(com);
    shockGroupRef.current.add(comMesh);
    const cylinderGeom = new THREE.CylinderGeometry(columnRadius, columnRadius, 20, 32, 1, true);
    const cylinderMat = new THREE.MeshBasicMaterial({
      color: 16739179,
      transparent: true,
      opacity: 0.3,
      side: THREE.DoubleSide,
      wireframe: true
    });
    const cylinder = new THREE.Mesh(cylinderGeom, cylinderMat);
    cylinder.position.set(com.x, com.y, 0);
    cylinder.rotation.x = Math.PI / 2;
    shockGroupRef.current.add(cylinder);
    const circleGeom = new THREE.RingGeometry(columnRadius - 0.01, columnRadius + 0.01, 64);
    const circleMat = new THREE.MeshBasicMaterial({ color: 16739179, side: THREE.DoubleSide });
    const circle = new THREE.Mesh(circleGeom, circleMat);
    circle.position.copy(com);
    shockGroupRef.current.add(circle);
    const sliceGeom = new THREE.BoxGeometry(40, sliceThickness * 2, sliceThickness * 2);
    const sliceMat = new THREE.MeshBasicMaterial({
      color: 5164484,
      transparent: true,
      opacity: 0.2,
      side: THREE.DoubleSide
    });
    const slice = new THREE.Mesh(sliceGeom, sliceMat);
    slice.position.copy(com);
    shockGroupRef.current.add(slice);
    const sliceEdges = new THREE.EdgesGeometry(sliceGeom);
    const sliceLineMat = new THREE.LineBasicMaterial({ color: 5164484 });
    const sliceLine = new THREE.LineSegments(sliceEdges, sliceLineMat);
    sliceLine.position.copy(com);
    shockGroupRef.current.add(sliceLine);
  }, [computeCoM2]);
  const updateParticles = useCallback(() => {
    const frameIndex = frameIndexRef.current ?? 0;
    const frames = framesRef.current;
    if (!frames || !geometryRef.current) return;
    const frame = frames.get(frameIndex);
    if (!frame) return;
    const currentColorField = colorFieldRef.current;
    const currentLogScale = logScaleRef.current;
    const needsUpdate = frameIndex !== lastFrameIndexRef.current || currentColorField !== lastColorFieldRef.current || currentLogScale !== lastLogScaleRef.current;
    if (!needsUpdate) return;
    lastFrameIndexRef.current = frameIndex;
    lastColorFieldRef.current = currentColorField;
    lastLogScaleRef.current = currentLogScale;
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
    switch (currentColorField) {
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
    const currentColorMap = colorMapRef.current;
    const colorMapKey = currentColorMap.name || "Cosmic Dawn";
    for (let i = 0; i < frame.particleCount; i++) {
      let value = fieldData[i];
      if (currentLogScale && value > 0 && vMin > 0) {
        const logMin = Math.log10(vMin);
        const logMax = Math.log10(vMax);
        value = (Math.log10(value) - logMin) / (logMax - logMin);
      } else {
        value = (value - vMin) / (vMax - vMin);
      }
      value = Math.max(0, Math.min(1, value));
      const rgb = sampleColorMap(colorMapKey, value);
      colors[i * 3] = rgb[0];
      colors[i * 3 + 1] = rgb[1];
      colors[i * 3 + 2] = rgb[2];
    }
    geometry.attributes.color.needsUpdate = true;
    if (showShockSamplingRef.current) {
      updateShockSampling(frame);
    }
  }, [framesRef, frameIndexRef, updateShockSampling]);
  useEffect(() => {
    updateParticlesRef.current = updateParticles;
  }, [updateParticles]);
  useEffect(() => {
    computeStatsRef.current = computeStats;
  }, [computeStats]);
  useEffect(() => {
    if (!containerRef.current) return;
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setClearColor(657940, 1);
    containerRef.current.appendChild(renderer.domElement);
    rendererRef.current = renderer;
    const scene = new THREE.Scene();
    sceneRef.current = scene;
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1e3);
    const config = CAMERA_POSITIONS[projection] || CAMERA_POSITIONS.xy;
    camera.position.copy(config.position);
    camera.up.copy(config.up);
    camera.lookAt(0, 0, 0);
    cameraRef.current = camera;
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.1;
    controls.rotateSpeed = 0.5;
    controls.zoomSpeed = 1.2;
    controls.panSpeed = 0.8;
    controls.minDistance = 5;
    controls.maxDistance = 200;
    controls.target.set(0, 0, 0);
    controlsRef.current = controls;
    controls.addEventListener("change", () => {
      const currentPos = camera.position.clone().normalize();
      const defaultPos = config.position.clone().normalize();
      const dotProduct = currentPos.dot(defaultPos);
      setIsRotated(dotProduct < 0.99);
    });
    const geometry = new THREE.BufferGeometry();
    const initialPositions = new Float32Array(1e5 * 3);
    const initialColors = new Float32Array(1e5 * 3);
    geometry.setAttribute("position", new THREE.BufferAttribute(initialPositions, 3).setUsage(THREE.DynamicDrawUsage));
    geometry.setAttribute("color", new THREE.BufferAttribute(initialColors, 3).setUsage(THREE.DynamicDrawUsage));
    geometryRef.current = geometry;
    const texture = getCircleTexture();
    const material = new THREE.PointsMaterial({
      size: particleSize * 0.08,
      vertexColors: true,
      transparent: true,
      opacity: 0.9,
      sizeAttenuation: true,
      map: texture,
      alphaTest: 0.01,
      depthWrite: false,
      blending: THREE.AdditiveBlending
    });
    materialRef.current = material;
    const particles = new THREE.Points(geometry, material);
    scene.add(particles);
    particlesRef.current = particles;
    const gridHelper = new THREE.GridHelper(40, 20, 4473924, 2236962);
    if (projection === "xy") {
      gridHelper.rotation.x = Math.PI / 2;
    } else if (projection === "xz") ;
    else if (projection === "yz") {
      gridHelper.rotation.z = Math.PI / 2;
    }
    scene.add(gridHelper);
    const axesHelper = new THREE.AxesHelper(10);
    scene.add(axesHelper);
    const shockGroup = new THREE.Group();
    scene.add(shockGroup);
    shockGroupRef.current = shockGroup;
    computeStatsRef.current();
    const animate = () => {
      animationFrameRef.current = requestAnimationFrame(animate);
      updateParticlesRef.current();
      controls.update();
      renderer.render(scene, camera);
    };
    animate();
    return () => {
      if (animationFrameRef.current) {
        cancelAnimationFrame(animationFrameRef.current);
      }
      controls.dispose();
      renderer.dispose();
      geometry.dispose();
      material.dispose();
      if (containerRef.current && renderer.domElement) {
        containerRef.current.removeChild(renderer.domElement);
      }
    };
  }, [projection]);
  useEffect(() => {
    const renderer = rendererRef.current;
    const camera = cameraRef.current;
    if (renderer && camera) {
      renderer.setSize(width, height);
      camera.aspect = width / height;
      camera.updateProjectionMatrix();
    }
  }, [width, height]);
  useEffect(() => {
    const material = materialRef.current;
    if (material) {
      material.size = particleSize * 0.08;
    }
  }, [particleSize]);
  return /* @__PURE__ */ jsxs("div", { className: `relative ${className}`, children: [
    /* @__PURE__ */ jsx(
      "div",
      {
        ref: containerRef,
        style: { width, height },
        className: "rounded border border-gray-700 overflow-hidden"
      }
    ),
    /* @__PURE__ */ jsx("div", { className: "absolute top-1 left-1 text-cyan-300 text-sm font-bold bg-black/70 px-2 py-0.5 rounded", children: projection.toUpperCase() }),
    /* @__PURE__ */ jsx("div", { className: "absolute top-1 right-1 flex gap-1", children: /* @__PURE__ */ jsx(
      "button",
      {
        onClick: resetCamera,
        className: `
            px-1.5 py-0.5 text-[10px] rounded transition-all
            ${isRotated ? "bg-cyan-600 hover:bg-cyan-500 text-white" : "bg-gray-700/80 text-gray-400 hover:bg-gray-600"}
          `,
        title: "Reset camera to default view",
        children: "⟲ Reset"
      }
    ) }),
    /* @__PURE__ */ jsx("div", { className: "absolute bottom-1 left-1 text-gray-500 text-[9px] bg-black/50 px-1 rounded", children: "Drag: rotate | Scroll: zoom | Shift+drag: pan" }),
    showShockSampling && /* @__PURE__ */ jsxs("div", { className: "absolute bottom-1 right-1 text-[9px] bg-black/70 px-1 rounded", children: [
      /* @__PURE__ */ jsx("span", { className: "text-red-400", children: "● Z-col" }),
      /* @__PURE__ */ jsx("span", { className: "mx-1 text-gray-500", children: "|" }),
      /* @__PURE__ */ jsx("span", { className: "text-teal-400", children: "■ X-slice" })
    ] })
  ] });
});
Projection3DInteractive.displayName = "Projection3DInteractive";
const GAMMA = 5 / 3;
const COLORS = {
  background: "#0f0f14",
  panel: "#1a1a24",
  grid: "#2a2a3a",
  text: "#e0e0e0",
  accent: "#00b4d8",
  density: "#ffd166",
  // Gold
  pressure: "#ef476f",
  // Pink-red
  temperature: "#06d6a0",
  // Mint
  pAdiabatic: "#ffd166",
  // Gold (dashed)
  tAdiabatic: "#06d6a0"
  // Mint (dashed)
};
function computeCoM(frame) {
  let totalMass = 0;
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < frame.particleCount; i++) {
    const m = frame.mass?.[i] ?? 1;
    const x = frame.positions[i * 3];
    const y = frame.positions[i * 3 + 1];
    const z = frame.positions[i * 3 + 2];
    totalMass += m;
    cx += m * x;
    cy += m * y;
    cz += m * z;
  }
  if (totalMass > 0) {
    cx /= totalMass;
    cy /= totalMass;
    cz /= totalMass;
  }
  return { x: cx, y: cy, z: cz };
}
function computeTemperature(pressure, density) {
  const n = pressure.length;
  const temperature = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const rho = density[i];
    const p = pressure[i];
    temperature[i] = rho > 0 ? p / rho : 0;
  }
  return temperature;
}
function extractVerticalProfile(frame, com, temperature, initial, columnRadius = 0.15, nBins = 40) {
  const positions = [];
  const density = [];
  const pressure = [];
  const temp = [];
  const pAdiabatic = [];
  const tAdiabatic = [];
  const particlesInColumn = [];
  for (let i = 0; i < frame.particleCount; i++) {
    const dx = frame.positions[i * 3] - com.x;
    const dy = frame.positions[i * 3 + 1] - com.y;
    const rXY = Math.sqrt(dx * dx + dy * dy);
    if (rXY < columnRadius) {
      particlesInColumn.push({
        z: frame.positions[i * 3 + 2] - com.z,
        dens: frame.density?.[i] ?? 1,
        pres: frame.pressure?.[i] ?? 1,
        temp: temperature[i],
        mass: frame.mass?.[i] ?? 1
      });
    }
  }
  if (particlesInColumn.length < 10) {
    return { positions: [], density: [], pressure: [], temperature: [], pAdiabatic: [], tAdiabatic: [] };
  }
  const zMin = Math.min(...particlesInColumn.map((p) => p.z));
  const zMax = Math.max(...particlesInColumn.map((p) => p.z));
  const binWidth = (zMax - zMin) / nBins;
  for (let bin = 0; bin < nBins; bin++) {
    const binStart = zMin + bin * binWidth;
    const binEnd = binStart + binWidth;
    const binCenter = (binStart + binEnd) / 2;
    const inBin = particlesInColumn.filter((p) => p.z >= binStart && p.z < binEnd);
    if (inBin.length === 0) continue;
    const totalMass = inBin.reduce((sum, p) => sum + p.mass, 0);
    const avgDens = inBin.reduce((sum, p) => sum + p.mass * p.dens, 0) / totalMass;
    const avgPres = inBin.reduce((sum, p) => sum + p.mass * p.pres, 0) / totalMass;
    const avgTemp = inBin.reduce((sum, p) => sum + p.mass * p.temp, 0) / totalMass;
    const rhoRatio = avgDens / initial.rho;
    positions.push(binCenter);
    density.push(rhoRatio);
    pressure.push(avgPres / initial.P);
    temp.push(avgTemp / initial.T);
    pAdiabatic.push(Math.pow(rhoRatio, GAMMA));
    tAdiabatic.push(Math.pow(rhoRatio, GAMMA - 1));
  }
  return { positions, density, pressure, temperature: temp, pAdiabatic, tAdiabatic };
}
function extractHorizontalProfile(frame, com, temperature, initial, sliceThickness = 0.15, nBins = 60) {
  const positions = [];
  const density = [];
  const pressure = [];
  const temp = [];
  const pAdiabatic = [];
  const tAdiabatic = [];
  const particlesInSlice = [];
  for (let i = 0; i < frame.particleCount; i++) {
    const dz = Math.abs(frame.positions[i * 3 + 2] - com.z);
    const dy = Math.abs(frame.positions[i * 3 + 1] - com.y);
    if (dz < sliceThickness && dy < sliceThickness) {
      particlesInSlice.push({
        x: frame.positions[i * 3] - com.x,
        dens: frame.density?.[i] ?? 1,
        pres: frame.pressure?.[i] ?? 1,
        temp: temperature[i],
        mass: frame.mass?.[i] ?? 1
      });
    }
  }
  if (particlesInSlice.length < 10) {
    return { positions: [], density: [], pressure: [], temperature: [], pAdiabatic: [], tAdiabatic: [] };
  }
  const xMin = Math.min(...particlesInSlice.map((p) => p.x));
  const xMax = Math.max(...particlesInSlice.map((p) => p.x));
  const binWidth = (xMax - xMin) / nBins;
  for (let bin = 0; bin < nBins; bin++) {
    const binStart = xMin + bin * binWidth;
    const binEnd = binStart + binWidth;
    const binCenter = (binStart + binEnd) / 2;
    const inBin = particlesInSlice.filter((p) => p.x >= binStart && p.x < binEnd);
    if (inBin.length === 0) continue;
    const totalMass = inBin.reduce((sum, p) => sum + p.mass, 0);
    const avgDens = inBin.reduce((sum, p) => sum + p.mass * p.dens, 0) / totalMass;
    const avgPres = inBin.reduce((sum, p) => sum + p.mass * p.pres, 0) / totalMass;
    const avgTemp = inBin.reduce((sum, p) => sum + p.mass * p.temp, 0) / totalMass;
    const rhoRatio = avgDens / initial.rho;
    positions.push(binCenter);
    density.push(rhoRatio);
    pressure.push(avgPres / initial.P);
    temp.push(avgTemp / initial.T);
    pAdiabatic.push(Math.pow(rhoRatio, GAMMA));
    tAdiabatic.push(Math.pow(rhoRatio, GAMMA - 1));
  }
  return { positions, density, pressure, temperature: temp, pAdiabatic, tAdiabatic };
}
function drawProfile(ctx, profile, x, y, width, height, title, xLabel, yRange = [0, 5]) {
  const margin = { left: 45, right: 10, top: 25, bottom: 30 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  ctx.fillStyle = COLORS.panel;
  ctx.fillRect(x, y, width, height);
  ctx.fillStyle = COLORS.text;
  ctx.font = "bold 11px sans-serif";
  ctx.textAlign = "center";
  ctx.fillText(title, x + width / 2, y + 15);
  if (profile.positions.length === 0) {
    ctx.fillStyle = "#666";
    ctx.font = "10px sans-serif";
    ctx.fillText("No data", x + width / 2, y + height / 2);
    return;
  }
  ctx.strokeStyle = COLORS.grid;
  ctx.lineWidth = 0.5;
  ctx.setLineDash([2, 2]);
  for (let i = 0; i <= 4; i++) {
    const gy = y + margin.top + i / 4 * plotHeight;
    ctx.beginPath();
    ctx.moveTo(x + margin.left, gy);
    ctx.lineTo(x + margin.left + plotWidth, gy);
    ctx.stroke();
  }
  ctx.setLineDash([]);
  const xMin = Math.min(...profile.positions);
  const xMax = Math.max(...profile.positions);
  const scaleX = (v) => x + margin.left + (v - xMin) / (xMax - xMin) * plotWidth;
  const scaleY = (v) => y + margin.top + plotHeight - (v - yRange[0]) / (yRange[1] - yRange[0]) * plotHeight;
  ctx.strokeStyle = "#666";
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(x + margin.left, scaleY(1));
  ctx.lineTo(x + margin.left + plotWidth, scaleY(1));
  ctx.stroke();
  ctx.setLineDash([]);
  const drawLine = (data, color, dashed = false) => {
    if (data.length === 0) return;
    ctx.strokeStyle = color;
    ctx.lineWidth = dashed ? 1.5 : 2;
    if (dashed) ctx.setLineDash([4, 2]);
    ctx.beginPath();
    for (let i = 0; i < data.length; i++) {
      const px = scaleX(profile.positions[i]);
      const py = scaleY(Math.max(yRange[0], Math.min(yRange[1], data[i])));
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  };
  drawLine(profile.pAdiabatic, COLORS.pAdiabatic, true);
  drawLine(profile.tAdiabatic, COLORS.tAdiabatic, true);
  drawLine(profile.density, COLORS.density);
  drawLine(profile.pressure, COLORS.pressure);
  drawLine(profile.temperature, COLORS.temperature);
  ctx.strokeStyle = COLORS.accent;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(x + margin.left, y + margin.top);
  ctx.lineTo(x + margin.left, y + height - margin.bottom);
  ctx.lineTo(x + width - margin.right, y + height - margin.bottom);
  ctx.stroke();
  ctx.fillStyle = COLORS.text;
  ctx.font = "9px sans-serif";
  ctx.textAlign = "center";
  ctx.fillText(xLabel, x + width / 2, y + height - 5);
  ctx.textAlign = "right";
  ctx.fillText(yRange[1].toFixed(0), x + margin.left - 5, y + margin.top + 5);
  ctx.fillText(yRange[0].toFixed(0), x + margin.left - 5, y + height - margin.bottom);
}
function computeInitialValues(initialFrame, frame) {
  const refFrame = initialFrame ?? frame;
  if (!refFrame) return { T: 50, rho: 1e-3, P: 1e-6 };
  const densities = Array.from(refFrame.density ?? []).filter((v) => v > 0);
  const pressures = Array.from(refFrame.pressure ?? []).filter((v) => v > 0);
  densities.sort((a, b) => a - b);
  pressures.sort((a, b) => a - b);
  const medianRho = densities[Math.floor(densities.length / 2)] || 1e-3;
  const medianP = pressures[Math.floor(pressures.length / 2)] || 1e-6;
  const temp = computeTemperature(
    new Float32Array([medianP]),
    new Float32Array([medianRho])
  );
  return { T: temp[0], rho: medianRho, P: medianP };
}
function computeFrameData(frame, initialValues) {
  const com = computeCoM(frame);
  const temperature = computeTemperature(
    frame.pressure ?? new Float32Array(frame.particleCount),
    frame.density ?? new Float32Array(frame.particleCount)
  );
  const tRatio = new Float32Array(frame.particleCount);
  for (let i = 0; i < frame.particleCount; i++) {
    tRatio[i] = temperature[i] / initialValues.T;
  }
  const verticalProfile = extractVerticalProfile(frame, com, temperature, initialValues);
  const horizontalProfile = extractHorizontalProfile(frame, com, temperature, initialValues);
  return {
    com,
    temperature,
    tRatio,
    initialValues,
    verticalProfile,
    horizontalProfile
  };
}
function ShockDiagnosticsPanelImperative({
  framesRef,
  frameIndexRef,
  initialFrame,
  width = 500,
  height = 300,
  className = ""
}) {
  const canvasRef = useRef(null);
  const animationFrameRef = useRef(null);
  const lastFrameIndexRef = useRef(-1);
  const initialValuesRef = useRef(null);
  const draw = useCallback((frame) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    if (!initialValuesRef.current) {
      initialValuesRef.current = computeInitialValues(initialFrame, frame);
    }
    const initialValues = initialValuesRef.current;
    const computedData = computeFrameData(frame, initialValues);
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    ctx.scale(dpr, dpr);
    ctx.fillStyle = COLORS.background;
    ctx.fillRect(0, 0, width, height);
    const chartWidth = (width - 30) / 2;
    const chartHeight = height - 60;
    ctx.fillStyle = COLORS.text;
    ctx.font = "bold 10px monospace";
    ctx.textAlign = "left";
    ctx.fillText(
      `CoM: (${computedData.com.x.toFixed(2)}, ${computedData.com.y.toFixed(2)}, ${computedData.com.z.toFixed(2)}) pc`,
      10,
      15
    );
    const distImbh = Math.sqrt(
      computedData.com.x ** 2 + computedData.com.y ** 2 + computedData.com.z ** 2
    );
    ctx.fillText(`IMBH→CoM: ${distImbh.toFixed(2)} pc`, 10, 28);
    let maxTRatio = 0, maxRhoRatio = 0, maxPRatio = 0;
    for (let i = 0; i < frame.particleCount; i++) {
      if (computedData.tRatio[i] > maxTRatio) maxTRatio = computedData.tRatio[i];
      const rhoRatio = (frame.density?.[i] ?? 0) / initialValues.rho;
      if (rhoRatio > maxRhoRatio) maxRhoRatio = rhoRatio;
      const pRatio = (frame.pressure?.[i] ?? 0) / initialValues.P;
      if (pRatio > maxPRatio) maxPRatio = pRatio;
    }
    ctx.textAlign = "right";
    ctx.font = "9px monospace";
    ctx.fillStyle = COLORS.temperature;
    ctx.fillText(`T_max/T₀=${maxTRatio.toFixed(1)}`, width - 10, 15);
    ctx.fillStyle = COLORS.density;
    ctx.fillText(`ρ_max/ρ₀=${maxRhoRatio.toFixed(1)}`, width - 10, 26);
    ctx.fillStyle = COLORS.pressure;
    ctx.fillText(`P_max/P₀=${maxPRatio.toFixed(1)}`, width - 10, 37);
    drawProfile(
      ctx,
      computedData.verticalProfile,
      10,
      45,
      chartWidth,
      chartHeight,
      "VERTICAL SHOCK (Z)",
      "Z - Z_CoM (pc)",
      [0, 5]
    );
    drawProfile(
      ctx,
      computedData.horizontalProfile,
      20 + chartWidth,
      45,
      chartWidth,
      chartHeight,
      "HORIZONTAL STRETCH (X)",
      "X - X_CoM (pc)",
      [0, 5]
    );
    const legendY = height - 12;
    ctx.font = "8px sans-serif";
    ctx.textAlign = "left";
    const legendItems = [
      { color: COLORS.density, label: "ρ/ρ₀" },
      { color: COLORS.pressure, label: "P/P₀" },
      { color: COLORS.temperature, label: "T/T₀" },
      { color: COLORS.pAdiabatic, label: "P_ad", dashed: true },
      { color: COLORS.tAdiabatic, label: "T_ad", dashed: true }
    ];
    let legendX = 10;
    for (const item of legendItems) {
      ctx.fillStyle = item.color;
      ctx.fillRect(legendX, legendY - 5, item.dashed ? 15 : 12, 2);
      if (item.dashed) {
        ctx.fillStyle = COLORS.background;
        ctx.fillRect(legendX + 4, legendY - 5, 3, 2);
        ctx.fillRect(legendX + 10, legendY - 5, 2, 2);
      }
      ctx.fillStyle = COLORS.text;
      ctx.fillText(item.label, legendX + (item.dashed ? 18 : 15), legendY);
      legendX += 55;
    }
  }, [initialFrame, width, height]);
  useEffect(() => {
    const tick = () => {
      const frames = framesRef.current;
      const frameIndex = frameIndexRef.current ?? 0;
      if (frames && frameIndex !== lastFrameIndexRef.current) {
        const frame = frames.get(frameIndex);
        if (frame) {
          draw(frame);
          lastFrameIndexRef.current = frameIndex;
        }
      }
      animationFrameRef.current = requestAnimationFrame(tick);
    };
    animationFrameRef.current = requestAnimationFrame(tick);
    return () => {
      if (animationFrameRef.current !== null) {
        cancelAnimationFrame(animationFrameRef.current);
      }
    };
  }, [framesRef, frameIndexRef, draw]);
  useEffect(() => {
    initialValuesRef.current = null;
    const frames = framesRef.current;
    const frameIndex = frameIndexRef.current ?? 0;
    if (frames) {
      const frame = frames.get(frameIndex);
      if (frame) {
        draw(frame);
      }
    }
  }, [initialFrame, framesRef, frameIndexRef, draw]);
  useEffect(() => {
    const frames = framesRef.current;
    const frameIndex = frameIndexRef.current ?? 0;
    if (frames) {
      const frame = frames.get(frameIndex);
      if (frame) {
        draw(frame);
      }
    }
  }, [width, height, framesRef, frameIndexRef, draw]);
  return /* @__PURE__ */ jsx("div", { className: `relative ${className}`, children: /* @__PURE__ */ jsx(
    "canvas",
    {
      ref: canvasRef,
      style: { width, height },
      className: "rounded border border-gray-700"
    }
  ) });
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
  return /* @__PURE__ */ jsxs("div", { className: `bg-gray-800 p-2 rounded h-full flex flex-col ${className}`, children: [
    /* @__PURE__ */ jsx("h3", { className: "text-xs font-medium text-gray-300 mb-1 shrink-0", children: "Energy Evolution" }),
    /* @__PURE__ */ jsx(ResponsiveContainer, { width: "100%", height: "100%", children: /* @__PURE__ */ jsxs(LineChart, { data, children: [
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
  return /* @__PURE__ */ jsxs("div", { className: `bg-gray-800 p-2 rounded h-full flex flex-col ${className}`, children: [
    /* @__PURE__ */ jsx("h3", { className: "text-xs font-medium text-gray-300 mb-1 shrink-0", children: "Momentum Evolution" }),
    /* @__PURE__ */ jsx(ResponsiveContainer, { width: "100%", height: "100%", children: /* @__PURE__ */ jsxs(LineChart, { data, children: [
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
function computeOrbitalParams(cloudPos, cloudVel, bhMass, G = 1) {
  const r = Math.sqrt(cloudPos[0] ** 2 + cloudPos[1] ** 2 + cloudPos[2] ** 2);
  const v = Math.sqrt(cloudVel[0] ** 2 + cloudVel[1] ** 2 + cloudVel[2] ** 2);
  const mu = G * bhMass;
  const specificEnergy = v * v / 2 - mu / r;
  const a = -mu / (2 * specificEnergy);
  const Lx = cloudPos[1] * cloudVel[2] - cloudPos[2] * cloudVel[1];
  const Ly = cloudPos[2] * cloudVel[0] - cloudPos[0] * cloudVel[2];
  const Lz = cloudPos[0] * cloudVel[1] - cloudPos[1] * cloudVel[0];
  const L = Math.sqrt(Lx ** 2 + Ly ** 2 + Lz ** 2);
  const eccentricity = Math.sqrt(1 + 2 * specificEnergy * L * L / (mu * mu));
  const pericentre = eccentricity > 1 ? Math.abs(a) * (eccentricity - 1) : Math.abs(a) * (1 - eccentricity);
  return { pericentre, eccentricity };
}
function presetToIMBHConfig(preset) {
  const imbhParams = preset.imbh_parameters;
  const externalForce = preset.externalForces?.pointMassBH;
  if (!imbhParams?.enabled && !externalForce?.enabled) {
    return null;
  }
  const bhMassCodeUnits = externalForce?.mass ?? (imbhParams?.M_BH ? imbhParams.M_BH / 1e3 : 100);
  const bhPosition = externalForce?.position ?? imbhParams?.BH_initial_position ?? [0, 0, 0];
  const cloudPos = imbhParams?.cloud_initial_position ?? preset.initialCondition?.transform?.translate ?? [-20, 5.17, 0];
  const cloudVel = imbhParams?.cloud_initial_velocity ?? preset.initialCondition?.transform?.velocity_boost ?? [10, 0, 0];
  const cloudMass = preset.physics_summary?.cloud_mass_Msun ? preset.physics_summary.cloud_mass_Msun / 1e3 : 1;
  const cloudRadius = preset.physics_summary?.cloud_radius_pc ?? 1.13;
  const tidalRadius = preset.physics_summary?.tidal_radius_pc ?? 5.24;
  const impactParameter = preset.physics_summary?.impact_parameter_pc ?? Math.abs(cloudPos[1]);
  const { pericentre, eccentricity } = computeOrbitalParams(cloudPos, cloudVel, bhMassCodeUnits);
  const timeUnit = 0.978;
  return {
    enabled: true,
    bhPosition,
    bhMass: bhMassCodeUnits,
    cloudInitialPosition: cloudPos,
    cloudInitialVelocity: cloudVel,
    cloudMass,
    cloudRadius,
    tidalRadius,
    impactParameter,
    pericentre,
    eccentricity,
    timeUnit,
    // Optional galactic parameters (can be overridden)
    inclination: 70,
    positionAngle: 41.6,
    lsrVelocity: -120
  };
}
const DEFAULT_IMBH_CONFIG = {
  enabled: true,
  bhPosition: [0, 0, 0],
  bhMass: 100,
  // 10^5 M_sun in code units (1000 M_sun)
  cloudInitialPosition: [20, -5.17, 0],
  cloudInitialVelocity: [-10.18, 5.05, 0],
  cloudMass: 1,
  // 1000 M_sun in code units
  cloudRadius: 1.13,
  tidalRadius: 5.24,
  impactParameter: 5.17,
  pericentre: 2.217,
  eccentricity: 1.4504,
  timeUnit: 0.978,
  inclination: 70,
  positionAngle: 41.6,
  lsrVelocity: -120
};
const PRESET_SCENARIOS = [
  {
    id: "Mc1e3_Mbh1e5_b1p5_v10_adiabatic_gsph",
    name: "Strong Disruption (b=1.5pc, Adiabatic, GSPH)",
    path: "simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/adiabatic_61k_gsph.json"
  },
  {
    id: "Mc1e3_Mbh1e5_b1p5_v10_adiabatic_gdisph",
    name: "Strong Disruption (b=1.5pc, Adiabatic, GDISPH)",
    path: "simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/adiabatic_61k_gdisph.json"
  },
  {
    id: "Mc1e3_Mbh1e5_b1p5_v10_radiative_gsph",
    name: "Strong Disruption (b=1.5pc, Radiative, GSPH)",
    path: "simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/radiative_61k_gsph.json"
  },
  {
    id: "Mc1e3_Mbh1e5_b1p5_v10_radiative_gdisph",
    name: "Strong Disruption (b=1.5pc, Radiative, GDISPH)",
    path: "simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b1p5_v10/radiative_61k_gdisph.json"
  },
  {
    id: "Mc1e3_Mbh1e5_b3_v10_adiabatic",
    name: "Moderate Disruption (b=3pc)",
    path: "simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json"
  },
  {
    id: "Mc1e3_Mbh1e5_b6_v10_adiabatic",
    name: "Weak Disruption (b=6pc)",
    path: "simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b6_v10/adiabatic_61k_gsph.json"
  }
];
function ResizeHandle({ direction = "horizontal" }) {
  return /* @__PURE__ */ jsx(
    PanelResizeHandle,
    {
      className: `
        ${direction === "horizontal" ? "w-1.5 cursor-col-resize" : "h-1.5 cursor-row-resize"}
        bg-gray-700 hover:bg-cyan-500 active:bg-cyan-400 transition-colors
        flex items-center justify-center group
      `,
      children: /* @__PURE__ */ jsx("div", { className: `
        ${direction === "horizontal" ? "w-0.5 h-8" : "h-0.5 w-8"}
        bg-gray-500 group-hover:bg-cyan-300 rounded-full transition-colors
      ` })
    }
  );
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
  const [colorMapName, setColorMapName] = useState("cosmicDawn");
  const [pointSize, setPointSize] = useState(0.02);
  const [opacity, setOpacity] = useState(0.8);
  const [showAxes, setShowAxes] = useState(true);
  const [showBoundingBox, setShowBoundingBox] = useState(true);
  const [colorRange, setColorRange] = useState([0, 0]);
  const [useLogScale, setUseLogScale] = useState(false);
  const [projectionColorFields, setProjectionColorFields] = useState({
    xy: "density",
    xz: "pressure",
    yz: "velocity"
  });
  const [useMultiColorField, setUseMultiColorField] = useState(false);
  const [showProjections, setShowProjections] = useState(true);
  const [showCharts, setShowCharts] = useState(true);
  const [showShockDiagnostics, setShowShockDiagnostics] = useState(true);
  const [panelDimensions, setPanelDimensions] = useState({
    projection: { width: 300, height: 120 },
    shock: { width: 600, height: 250 }
  });
  const projectionPanelRef = useRef(null);
  const shockPanelRef = useRef(null);
  const projectionXYRef = useRef(null);
  const projectionXZRef = useRef(null);
  const projectionYZRef = useRef(null);
  const resetAllProjectionCameras = useCallback(() => {
    projectionXYRef.current?.resetCamera();
    projectionXZRef.current?.resetCamera();
    projectionYZRef.current?.resetCamera();
  }, []);
  useEffect(() => {
    const updateDimensions = () => {
      if (projectionPanelRef.current) {
        const rect = projectionPanelRef.current.getBoundingClientRect();
        const availableHeight = rect.height - 40;
        const projectionHeight = Math.floor(availableHeight / 3);
        setPanelDimensions((prev) => ({
          ...prev,
          projection: {
            width: Math.floor(rect.width) - 16,
            height: Math.max(60, projectionHeight)
            // minimum 60px per projection
          }
        }));
      }
      if (shockPanelRef.current) {
        const rect = shockPanelRef.current.getBoundingClientRect();
        setPanelDimensions((prev) => ({
          ...prev,
          shock: { width: Math.floor(rect.width) - 16, height: Math.floor(rect.height) - 60 }
          // Extra space for description
        }));
      }
    };
    updateDimensions();
    window.addEventListener("resize", updateDimensions);
    const observer = new ResizeObserver(updateDimensions);
    if (projectionPanelRef.current) observer.observe(projectionPanelRef.current);
    if (shockPanelRef.current) observer.observe(shockPanelRef.current);
    return () => {
      window.removeEventListener("resize", updateDimensions);
      observer.disconnect();
    };
  }, []);
  const [currentFps, setCurrentFps] = useState(0);
  const [showBlackHole, setShowBlackHole] = useState(true);
  const [showTrajectory, setShowTrajectory] = useState(true);
  const [showRadii, setShowRadii] = useState(true);
  const [showGalacticMarkers, setShowGalacticMarkers] = useState(true);
  const [showLabels, setShowLabels] = useState(true);
  const [showGalaxyDisk, setShowGalaxyDisk] = useState(true);
  const [showSolarSystem, setShowSolarSystem] = useState(true);
  const [animateGalaxy, setAnimateGalaxy] = useState(false);
  const [galaxyAnimationSpeed, setGalaxyAnimationSpeed] = useState(1);
  const [cameraMode, setCameraMode] = useState("free");
  const [selectedPresetId, setSelectedPresetId] = useState(null);
  const [loadedPresetConfig, setLoadedPresetConfig] = useState(null);
  useEffect(() => {
    if (!selectedPresetId) {
      setLoadedPresetConfig(null);
      return;
    }
    const preset = PRESET_SCENARIOS.find((p) => p.id === selectedPresetId);
    if (!preset) return;
    const loadPreset = async () => {
      try {
        const response = await fetch(`http://localhost:3001/api/preset?path=${encodeURIComponent(preset.path)}`);
        if (response.ok) {
          const presetData = await response.json();
          const config = presetToIMBHConfig(presetData);
          if (config) {
            console.log("[Dashboard] Loaded preset config:", preset.name, config);
            setLoadedPresetConfig(config);
          }
        } else {
          console.warn("[Dashboard] Failed to load preset from API, using default");
        }
      } catch (error2) {
        console.warn("[Dashboard] Error loading preset:", error2);
      }
    };
    loadPreset();
  }, [selectedPresetId]);
  const imbhPhysics = useMemo(() => {
    if (loadedPresetConfig) {
      console.log("[Dashboard] Using IMBH physics from loaded preset");
      return loadedPresetConfig;
    }
    if (simulation?.imbhPhysics) {
      console.log("[Dashboard] Using IMBH physics from simulation config:", simulation.imbhPhysics);
      return simulation.imbhPhysics;
    }
    console.log("[Dashboard] Using default IMBH physics");
    return DEFAULT_IMBH_CONFIG;
  }, [loadedPresetConfig, simulation]);
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
    const baseMap = COLOR_MAPS[colorMapName] || COLOR_MAPS.cosmicDawn;
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
          frameIndexRef.current = 0;
          setCurrentFrameIndex(0);
          lastTime = timestamp;
        } else {
          if (framesRef.current.has(nextFrame)) {
            frameIndexRef.current = nextFrame;
            setCurrentFrameIndex(nextFrame);
            lastTime = timestamp;
          }
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
    if (playing) {
      startImperativePlayback();
    } else {
      stopImperativePlayback();
    }
  }, [startImperativePlayback, stopImperativePlayback]);
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
            onClick: () => setShowProjections(!showProjections),
            className: `px-2 py-1 text-xs rounded ${showProjections ? "bg-blue-600 text-white" : "bg-gray-700 text-gray-400"}`,
            children: "2D Views"
          }
        ),
        /* @__PURE__ */ jsx(
          "button",
          {
            onClick: () => setShowShockDiagnostics(!showShockDiagnostics),
            className: `px-2 py-1 text-xs rounded ${showShockDiagnostics ? "bg-orange-600 text-white" : "bg-gray-700 text-gray-400"}`,
            children: "Shock"
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
      /* @__PURE__ */ jsxs("div", { className: "w-56 shrink-0 overflow-y-auto border-r border-gray-700 p-2", children: [
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
        /* @__PURE__ */ jsxs("div", { className: "mt-2 bg-gray-800 rounded p-3", children: [
          /* @__PURE__ */ jsxs("div", { className: "flex items-center justify-between mb-2", children: [
            /* @__PURE__ */ jsx("h3", { className: "text-sm font-medium text-gray-300", children: "2D Projections" }),
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-1 cursor-pointer", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: useMultiColorField,
                  onChange: (e) => setUseMultiColorField(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-xs text-gray-400", children: "Multi-field" })
            ] })
          ] }),
          useMultiColorField && /* @__PURE__ */ jsx("div", { className: "space-y-2 text-xs", children: ["xy", "xz", "yz"].map((proj) => /* @__PURE__ */ jsxs("div", { className: "flex items-center gap-2", children: [
            /* @__PURE__ */ jsxs("span", { className: "text-gray-400 w-8 uppercase font-mono", children: [
              proj,
              ":"
            ] }),
            /* @__PURE__ */ jsxs(
              "select",
              {
                value: projectionColorFields[proj],
                onChange: (e) => setProjectionColorFields((prev) => ({
                  ...prev,
                  [proj]: e.target.value
                })),
                className: "flex-1 bg-gray-700 text-gray-200 rounded px-1 py-0.5 border border-gray-600",
                children: [
                  /* @__PURE__ */ jsx("option", { value: "density", children: "Density" }),
                  /* @__PURE__ */ jsx("option", { value: "pressure", children: "Pressure" }),
                  /* @__PURE__ */ jsx("option", { value: "velocity", children: "Velocity" }),
                  /* @__PURE__ */ jsx("option", { value: "energy", children: "Energy" }),
                  /* @__PURE__ */ jsx("option", { value: "machNumber", children: "Mach #" })
                ]
              }
            )
          ] }, proj)) })
        ] }),
        /* @__PURE__ */ jsxs("div", { className: "mt-2 bg-gray-800 rounded p-3", children: [
          /* @__PURE__ */ jsx("h3", { className: "text-sm font-medium text-gray-300 mb-2", children: "IMBH Physics" }),
          /* @__PURE__ */ jsxs("div", { className: "mb-3", children: [
            /* @__PURE__ */ jsx("label", { className: "block text-xs text-gray-400 mb-1", children: "Preset Config" }),
            /* @__PURE__ */ jsxs(
              "select",
              {
                value: selectedPresetId ?? "",
                onChange: (e) => setSelectedPresetId(e.target.value || null),
                className: "w-full bg-gray-700 text-gray-200 text-xs rounded px-2 py-1 border border-gray-600",
                children: [
                  /* @__PURE__ */ jsx("option", { value: "", children: "Default (from simulation)" }),
                  PRESET_SCENARIOS.map((preset) => /* @__PURE__ */ jsx("option", { value: preset.id, children: preset.name }, preset.id))
                ]
              }
            )
          ] }),
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
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showGalacticMarkers,
                  onChange: (e) => setShowGalacticMarkers(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Galactic Markers" })
            ] }),
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", title: "Show text labels on all markers and arrows", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showLabels,
                  onChange: (e) => setShowLabels(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Labels" })
            ] }),
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", title: "Show schematic Milky Way disk with Sun at 8 kpc from GC", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showGalaxyDisk,
                  onChange: (e) => setShowGalaxyDisk(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Galaxy Disk" })
            ] }),
            /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", title: "Show observer direction with inclination (70°), position angle (41.6°), and V_LSR indicator", children: [
              /* @__PURE__ */ jsx(
                "input",
                {
                  type: "checkbox",
                  checked: showSolarSystem,
                  onChange: (e) => setShowSolarSystem(e.target.checked),
                  className: "rounded"
                }
              ),
              /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Show Observer Geometry" })
            ] }),
            showGalaxyDisk && /* @__PURE__ */ jsxs("div", { className: "mt-2 pt-2 border-t border-gray-600", children: [
              /* @__PURE__ */ jsxs("label", { className: "flex items-center gap-2 cursor-pointer", title: "Animate galactic rotation (Sun orbits GC in ~220 Myr)", children: [
                /* @__PURE__ */ jsx(
                  "input",
                  {
                    type: "checkbox",
                    checked: animateGalaxy,
                    onChange: (e) => setAnimateGalaxy(e.target.checked),
                    className: "rounded"
                  }
                ),
                /* @__PURE__ */ jsx("span", { className: "text-gray-300", children: "Animate Galaxy" })
              ] }),
              animateGalaxy && /* @__PURE__ */ jsxs("div", { className: "mt-2 ml-4", children: [
                /* @__PURE__ */ jsxs("div", { className: "text-gray-400 text-[10px] mb-1", children: [
                  "Speed: ",
                  galaxyAnimationSpeed.toFixed(1),
                  "× (1× = 220 Myr period)"
                ] }),
                /* @__PURE__ */ jsx(
                  "input",
                  {
                    type: "range",
                    min: "0.1",
                    max: "10",
                    step: "0.1",
                    value: galaxyAnimationSpeed,
                    onChange: (e) => setGalaxyAnimationSpeed(parseFloat(e.target.value)),
                    className: "w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer"
                  }
                ),
                /* @__PURE__ */ jsxs("div", { className: "flex justify-between text-[9px] text-gray-500 mt-1", children: [
                  /* @__PURE__ */ jsx("span", { children: "Slow" }),
                  /* @__PURE__ */ jsxs("span", { children: [
                    "~",
                    (220 / galaxyAnimationSpeed).toFixed(0),
                    " Myr/rev"
                  ] }),
                  /* @__PURE__ */ jsx("span", { children: "Fast" })
                ] })
              ] })
            ] }),
            /* @__PURE__ */ jsxs("div", { className: "mt-3 pt-2 border-t border-gray-700", children: [
              /* @__PURE__ */ jsx("div", { className: "text-gray-400 text-[10px] mb-2", children: "Camera View" }),
              /* @__PURE__ */ jsxs("div", { className: "flex gap-1", children: [
                /* @__PURE__ */ jsx(
                  "button",
                  {
                    onClick: () => setCameraMode("free"),
                    className: `flex-1 px-2 py-1 text-xs rounded ${cameraMode === "free" ? "bg-blue-600 text-white" : "bg-gray-700 text-gray-300 hover:bg-gray-600"}`,
                    title: "Free orbit camera controls",
                    children: "🎥 Free"
                  }
                ),
                /* @__PURE__ */ jsx(
                  "button",
                  {
                    onClick: () => setCameraMode("earth"),
                    className: `flex-1 px-2 py-1 text-xs rounded ${cameraMode === "earth" ? "bg-green-600 text-white" : "bg-gray-700 text-gray-300 hover:bg-gray-600"}`,
                    title: "View from Earth (observer's perspective)",
                    children: "🌍 Earth View"
                  }
                )
              ] })
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
      /* @__PURE__ */ jsx("div", { className: "flex-1 overflow-hidden", children: /* @__PURE__ */ jsxs(PanelGroup, { direction: "vertical", className: "h-full", children: [
        /* @__PURE__ */ jsx(Panel, { defaultSize: 50, minSize: 25, children: /* @__PURE__ */ jsxs(PanelGroup, { direction: "horizontal", className: "h-full", children: [
          /* @__PURE__ */ jsx(Panel, { defaultSize: 50, minSize: 25, children: /* @__PURE__ */ jsxs("div", { className: "h-full relative bg-gray-900", children: [
            isLoading && !currentFrame ? /* @__PURE__ */ jsx("div", { className: "flex items-center justify-center h-full text-gray-400", children: /* @__PURE__ */ jsxs("div", { className: "text-center", children: [
              /* @__PURE__ */ jsx("div", { className: "animate-spin text-2xl mb-2", children: "⏳" }),
              /* @__PURE__ */ jsxs("div", { children: [
                "Loading frame ",
                currentFrameIndex,
                "..."
              ] })
            ] }) }) : /* @__PURE__ */ jsx(
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
                className: "h-full w-full",
                onFpsUpdate: setCurrentFps,
                globalColorRange: globalColorStats,
                imbhPhysics,
                showBlackHole,
                showTrajectory,
                showRadii,
                showGalacticMarkers,
                showLabels,
                cameraMode,
                galacticConfig: {
                  distanceToGC: 60,
                  galacticLongitude: -0.398,
                  galacticLatitude: -0.224,
                  inclination: simulation?.imbhPhysics?.inclination ?? 70,
                  positionAngle: simulation?.imbhPhysics?.positionAngle ?? 41.6,
                  lsrVelocity: simulation?.imbhPhysics?.lsrVelocity ?? -120,
                  showGalaxyDisk,
                  showSolarSystem,
                  galaxyRotationSpeed: animateGalaxy ? 0.1 * galaxyAnimationSpeed : 0
                }
              }
            ),
            /* @__PURE__ */ jsxs("div", { className: "absolute bottom-2 right-2 text-green-400 text-xs font-mono bg-black/50 px-2 py-1 rounded", children: [
              "FPS: ",
              currentFps
            ] }),
            /* @__PURE__ */ jsx("div", { className: "absolute top-2 left-2 text-cyan-400 text-xs font-bold bg-black/60 px-2 py-1 rounded", children: "3D VIEW" })
          ] }) }),
          showProjections && /* @__PURE__ */ jsxs(Fragment, { children: [
            /* @__PURE__ */ jsx(ResizeHandle, { direction: "horizontal" }),
            /* @__PURE__ */ jsx(Panel, { defaultSize: 50, minSize: 25, children: /* @__PURE__ */ jsxs("div", { ref: projectionPanelRef, className: "h-full bg-gray-900 p-2 flex flex-col overflow-hidden", children: [
              /* @__PURE__ */ jsxs("div", { className: "flex items-center justify-between mb-1 shrink-0", children: [
                /* @__PURE__ */ jsx("div", { className: "text-cyan-400 text-xs font-bold", children: "3D PROJECTIONS" }),
                /* @__PURE__ */ jsx(
                  "button",
                  {
                    onClick: resetAllProjectionCameras,
                    className: "px-1.5 py-0.5 text-[10px] bg-gray-700 hover:bg-gray-600 text-gray-300 rounded transition-colors",
                    title: "Reset all cameras to default view",
                    children: "⟲ Reset All"
                  }
                )
              ] }),
              /* @__PURE__ */ jsxs("div", { className: "flex-1 flex flex-col gap-1 min-h-0 overflow-hidden", children: [
                /* @__PURE__ */ jsx("div", { className: "flex-1 min-h-0 overflow-hidden", children: /* @__PURE__ */ jsx(
                  Projection3DInteractive,
                  {
                    ref: projectionXYRef,
                    framesRef,
                    frameIndexRef,
                    projection: "xy",
                    colorField: useMultiColorField ? projectionColorFields.xy : colorField,
                    colorMap,
                    width: panelDimensions.projection.width,
                    height: panelDimensions.projection.height,
                    showShockSampling: showShockDiagnostics,
                    shockSamplingParams: { columnRadius: 0.15, sliceThickness: 0.15 },
                    globalColorRange: globalColorStats,
                    logScale: useLogScale,
                    particleSize: pointSize * 100
                  }
                ) }),
                /* @__PURE__ */ jsx("div", { className: "flex-1 min-h-0 overflow-hidden", children: /* @__PURE__ */ jsx(
                  Projection3DInteractive,
                  {
                    ref: projectionXZRef,
                    framesRef,
                    frameIndexRef,
                    projection: "xz",
                    colorField: useMultiColorField ? projectionColorFields.xz : colorField,
                    colorMap,
                    width: panelDimensions.projection.width,
                    height: panelDimensions.projection.height,
                    showShockSampling: showShockDiagnostics,
                    shockSamplingParams: { columnRadius: 0.15, sliceThickness: 0.15 },
                    globalColorRange: globalColorStats,
                    logScale: useLogScale,
                    particleSize: pointSize * 100
                  }
                ) }),
                /* @__PURE__ */ jsx("div", { className: "flex-1 min-h-0 overflow-hidden", children: /* @__PURE__ */ jsx(
                  Projection3DInteractive,
                  {
                    ref: projectionYZRef,
                    framesRef,
                    frameIndexRef,
                    projection: "yz",
                    colorField: useMultiColorField ? projectionColorFields.yz : colorField,
                    colorMap,
                    width: panelDimensions.projection.width,
                    height: panelDimensions.projection.height,
                    showShockSampling: showShockDiagnostics,
                    shockSamplingParams: { columnRadius: 0.15, sliceThickness: 0.15 },
                    globalColorRange: globalColorStats,
                    logScale: useLogScale,
                    particleSize: pointSize * 100
                  }
                ) })
              ] })
            ] }) })
          ] })
        ] }) }),
        (showShockDiagnostics || showCharts) && /* @__PURE__ */ jsxs(Fragment, { children: [
          /* @__PURE__ */ jsx(ResizeHandle, { direction: "vertical" }),
          /* @__PURE__ */ jsx(Panel, { defaultSize: 50, minSize: 25, children: /* @__PURE__ */ jsxs(PanelGroup, { direction: "horizontal", className: "h-full", children: [
            showShockDiagnostics && /* @__PURE__ */ jsx(Panel, { defaultSize: showCharts ? 50 : 100, minSize: 25, children: /* @__PURE__ */ jsxs("div", { ref: shockPanelRef, className: "h-full bg-gray-900 p-2 overflow-hidden flex flex-col", children: [
              /* @__PURE__ */ jsxs("div", { className: "shrink-0 mb-1", children: [
                /* @__PURE__ */ jsx("div", { className: "text-cyan-400 text-xs font-bold", children: "SHOCK DIAGNOSTICS" }),
                /* @__PURE__ */ jsxs("div", { className: "text-gray-500 text-[10px] leading-tight mt-0.5", children: [
                  /* @__PURE__ */ jsx("span", { className: "text-red-400", children: "●" }),
                  " Z-profile: cylinder r<0.15 pc around CoM (vertical compression)",
                  /* @__PURE__ */ jsx("span", { className: "mx-1", children: "|" }),
                  /* @__PURE__ */ jsx("span", { className: "text-teal-400", children: "━" }),
                  " X-profile: slice |y|,|z|<0.15 pc (tidal stretching)"
                ] })
              ] }),
              /* @__PURE__ */ jsx("div", { className: "flex-1 min-h-0", children: /* @__PURE__ */ jsx(
                ShockDiagnosticsPanelImperative,
                {
                  framesRef,
                  frameIndexRef,
                  initialFrame: frames.get(0) ?? null,
                  width: panelDimensions.shock.width,
                  height: panelDimensions.shock.height,
                  className: "w-full"
                }
              ) })
            ] }) }),
            showCharts && statistics.length > 0 && /* @__PURE__ */ jsxs(Fragment, { children: [
              showShockDiagnostics && /* @__PURE__ */ jsx(ResizeHandle, { direction: "horizontal" }),
              /* @__PURE__ */ jsx(Panel, { defaultSize: showShockDiagnostics ? 50 : 100, minSize: 25, children: /* @__PURE__ */ jsxs("div", { className: "h-full bg-gray-900 p-2 flex flex-col gap-2 overflow-hidden", children: [
                /* @__PURE__ */ jsx("div", { className: "text-cyan-400 text-xs font-bold", children: "ENERGY & MOMENTUM" }),
                /* @__PURE__ */ jsxs("div", { className: "flex-1 flex flex-col gap-2 min-h-0", children: [
                  /* @__PURE__ */ jsx("div", { className: "flex-1 min-h-0", children: /* @__PURE__ */ jsx(
                    EnergyChart,
                    {
                      statistics,
                      currentFrame: currentFrameIndex,
                      className: "h-full"
                    }
                  ) }),
                  /* @__PURE__ */ jsx("div", { className: "flex-1 min-h-0", children: /* @__PURE__ */ jsx(
                    MomentumChart,
                    {
                      statistics,
                      currentFrame: currentFrameIndex,
                      className: "h-full"
                    }
                  ) })
                ] })
              ] }) })
            ] })
          ] }) })
        ] })
      ] }) })
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
        imperativeMode: true,
        frameIndexRef,
        playbackSpeedRef
      }
    ) })
  ] });
}
const Route$4 = createFileRoute("/viz/")({
  component: VisualizationPage,
  loader: async () => {
    try {
      const response = await fetch("http://localhost:3000/api/simulations");
      const data = await response.json();
      console.log("[Loader] Fetched simulations:", data.simulations?.length || 0);
      return { simulations: data.simulations || [] };
    } catch (err) {
      console.error("[Loader] Failed:", err);
      return { simulations: [] };
    }
  }
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
  const loaderData = useLoaderData({ from: "/viz/" });
  const [simulations, setSimulations] = useState(loaderData?.simulations || []);
  const [selectedSimulation, setSelectedSimulation] = useState(null);
  const [frames, setFrames] = useState(/* @__PURE__ */ new Map());
  const [statistics, setStatistics] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [loadingProgress, setLoadingProgress] = useState({ loaded: 0, total: 0 });
  const [error, setError] = useState(null);
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const [isLoadingSimulations, setIsLoadingSimulations] = useState(false);
  useEffect(() => {
    if (loaderData?.simulations) {
      console.log("[useEffect] Setting simulations from loader:", loaderData.simulations.length);
      setSimulations(loaderData.simulations);
    }
  }, [loaderData]);
  const loadSimulations = useCallback(async () => {
    console.log("[loadSimulations] Starting fetch...");
    setIsLoadingSimulations(true);
    try {
      const response = await fetch("/api/simulations");
      console.log("[loadSimulations] Response status:", response.status);
      const data = await response.json();
      console.log("[loadSimulations] Data received:", data);
      console.log("[loadSimulations] Simulations count:", data.simulations?.length || 0);
      setSimulations(data.simulations || []);
    } catch (err) {
      console.error("Failed to load simulations:", err);
      setError("Failed to load simulations");
    } finally {
      setIsLoadingSimulations(false);
    }
  }, []);
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
    console.log("[selectSimulation] Called with:", sim.id, sim.name);
    setSelectedSimulation(sim);
    setFrames(/* @__PURE__ */ new Map());
    setStatistics([]);
    setError(null);
    setIsLoading(true);
    setLoadingProgress({ loaded: 0, total: sim.totalFrames });
    const simId = encodeURIComponent(sim.id.replace(/\//g, "|"));
    console.log("[selectSimulation] Encoded simId:", simId);
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
          !sidebarCollapsed && /* @__PURE__ */ jsx("div", { className: "flex-1 overflow-y-auto", children: isLoadingSimulations ? /* @__PURE__ */ jsxs("div", { className: "p-4 text-center text-gray-500 text-sm", children: [
            /* @__PURE__ */ jsx("div", { className: "mb-2", children: "Loading simulations..." }),
            /* @__PURE__ */ jsx(RefreshCw, { size: 16, className: "animate-spin mx-auto" })
          ] }) : simulations.length === 0 ? /* @__PURE__ */ jsxs("div", { className: "p-4 text-center text-gray-500 text-sm", children: [
            /* @__PURE__ */ jsx("div", { className: "mb-2", children: "No simulations found" }),
            /* @__PURE__ */ jsx("div", { className: "text-xs mb-2", children: "Run the data exporter to prepare simulation data" }),
            /* @__PURE__ */ jsx("div", { className: "text-xs text-yellow-500", children: "Debug: Check browser console (F12) for fetch logs" }),
            error && /* @__PURE__ */ jsxs("div", { className: "text-xs text-red-500 mt-2", children: [
              "Error: ",
              error
            ] })
          ] }) : /* @__PURE__ */ jsx("div", { className: "p-2 space-y-1", children: simulations.map((sim) => /* @__PURE__ */ jsxs(
            "button",
            {
              onClick: () => {
                console.log("[Button Click] Simulation clicked:", sim.id);
                alert(`Loading simulation: ${sim.name}`);
                selectSimulation(sim);
              },
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
async function findSimulationsServer() {
  const fs2 = await import("fs");
  const path2 = await import("path");
  const { fileURLToPath: fileURLToPath2 } = await import("url");
  const __dirname = path2.dirname(fileURLToPath2(import.meta.url));
  const dataRoot = process.env.SPH_DATA_ROOT || path2.resolve(__dirname, "../../../../../");
  console.log("🔍 Starting simulation scan...");
  console.log(`   Data root: ${dataRoot}`);
  const simulations = [];
  const startTime = Date.now();
  const categoryDirs = ["CAT1", "CAT2", "CAT3", "CAT_OKA"];
  const imbhResultsDir = path2.join(dataRoot, "simulations", "astrophysics", "imbh_cloud", "results");
  if (fs2.existsSync(imbhResultsDir)) {
    console.log("📁 Scanning simulations/astrophysics/imbh_cloud/results/ categories...");
    for (const catName of categoryDirs) {
      const catDir = path2.join(imbhResultsDir, catName);
      if (!fs2.existsSync(catDir)) {
        console.log(`   ⏭️  Skipping ${catName} (not found)`);
        continue;
      }
      console.log(`   📂 Scanning ${catName}...`);
      const scanCategoryDir = async (dir, prefix) => {
        const entries = fs2.readdirSync(dir, { withFileTypes: true });
        for (const entry of entries) {
          if (entry.isDirectory()) {
            const entryPath = path2.join(dir, entry.name);
            const simId = `imbh_cloud/${prefix}/${entry.name}`;
            const vizDataPath = path2.join(entryPath, "viz_data");
            const metadataPath = path2.join(vizDataPath, "metadata.json");
            if (fs2.existsSync(metadataPath)) {
              try {
                const metadata = JSON.parse(fs2.readFileSync(metadataPath, "utf-8"));
                simulations.push({
                  ...metadata,
                  id: simId,
                  dataPath: vizDataPath
                });
                console.log(`      ✓ Found: ${simId} (${metadata.totalFrames} frames)`);
              } catch (e) {
                console.error(`      ✗ Error reading ${metadataPath}:`, e);
              }
            } else {
              await scanCategoryDir(entryPath, `${prefix}/${entry.name}`);
            }
          }
        }
      };
      await scanCategoryDir(catDir, catName);
    }
  } else {
    console.log("⚠️  IMBH cloud results directory not found:", imbhResultsDir);
  }
  const elapsed = ((Date.now() - startTime) / 1e3).toFixed(2);
  console.log(`
✅ Scan complete: Found ${simulations.length} simulations in ${elapsed}s`);
  return simulations;
}
const Route$3 = createFileRoute("/api/simulations")({
  server: {
    handlers: {
      GET: async () => {
        try {
          const simulations = await findSimulationsServer();
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
const __filename$2 = fileURLToPath(import.meta.url);
const __dirname$2 = path.dirname(__filename$2);
const getDataRoot$1 = () => {
  return process.env.SPH_DATA_ROOT || path.resolve(__dirname$2, "../../../../");
};
const Route$2 = createFileRoute("/api/simulations/$simId")({
  server: {
    handlers: {
      GET: async ({ params }) => {
        try {
          const { simId } = params;
          const dataRoot = getDataRoot$1();
          const simPath = decodeURIComponent(simId).replace(/\|/g, "/");
          const pathParts = simPath.split("/");
          const possiblePaths = [
            path.join(dataRoot, "simulations", "astrophysics", simPath, "viz_data"),
            path.join(dataRoot, "simulations", "astrophysics", simPath, "results", "viz_data"),
            path.join(dataRoot, "lane_emden", "results", simPath, "viz_data"),
            path.join(dataRoot, simPath, "viz_data")
          ];
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts;
            possiblePaths.unshift(
              path.join(dataRoot, "simulations", "astrophysics", testName, "results", ...rest, "viz_data")
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
const Route$1 = createFileRoute("/api/simulations/$simId/frames/$frameId")({
  validateSearch: (search) => ({
    format: search.format || "json"
  }),
  server: {
    handlers: {
      GET: async ({ params, request }) => {
        const fs2 = await import("fs");
        const path2 = await import("path");
        const { fileURLToPath: fileURLToPath2 } = await import("url");
        const __dirname = path2.dirname(fileURLToPath2(import.meta.url));
        const dataRoot = process.env.SPH_DATA_ROOT || path2.resolve(__dirname, "../../../../../");
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
          const simPath = decodeURIComponent(simId).replace(/\|/g, "/");
          console.log(`   Decoded simPath: ${simPath}`);
          console.log(`   Data root: ${dataRoot}`);
          const pathParts = simPath.split("/");
          const possiblePaths = [
            path2.join(dataRoot, "simulations", "astrophysics", simPath, "viz_data"),
            path2.join(dataRoot, "simulations", "astrophysics", simPath, "results", "viz_data"),
            path2.join(dataRoot, "lane_emden", "results", simPath, "viz_data"),
            path2.join(dataRoot, simPath, "viz_data")
          ];
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts;
            possiblePaths.unshift(
              path2.join(dataRoot, "simulations", "astrophysics", testName, "results", ...rest, "viz_data")
            );
          }
          console.log(`   Possible paths:`);
          possiblePaths.forEach((p, i) => console.log(`     ${i}: ${p}`));
          let dataPath = null;
          for (const p of possiblePaths) {
            const exists = fs2.existsSync(p);
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
          const framePath = path2.join(dataPath, `frame_${frameIndex.toString().padStart(5, "0")}.bin`);
          if (!fs2.existsSync(framePath)) {
            if (wantBinary) {
              return new Response(`Frame ${frameIndex} not found`, { status: 404 });
            }
            return json({ error: `Frame ${frameIndex} not found` }, { status: 404 });
          }
          const buffer = fs2.readFileSync(framePath);
          const metadataPath = path2.join(dataPath, "metadata.json");
          let metadata = {};
          if (fs2.existsSync(metadataPath)) {
            metadata = JSON.parse(fs2.readFileSync(metadataPath, "utf-8"));
          }
          const frameInfoPath = path2.join(dataPath, `frame_${frameIndex.toString().padStart(5, "0")}.json`);
          let frameInfo = { time: frameIndex };
          if (fs2.existsSync(frameInfoPath)) {
            frameInfo = JSON.parse(fs2.readFileSync(frameInfoPath, "utf-8"));
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
const Route2 = createFileRoute("/api/simulations/$simId/frames/$frameId/bin")({
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
const IndexRoute = Route$5.update({
  id: "/",
  path: "/",
  getParentRoute: () => Route$6
});
const VizIndexRoute = Route$4.update({
  id: "/viz/",
  path: "/viz/",
  getParentRoute: () => Route$6
});
const ApiSimulationsRoute = Route$3.update({
  id: "/api/simulations",
  path: "/api/simulations",
  getParentRoute: () => Route$6
});
const ApiSimulationsSimIdRoute = Route$2.update({
  id: "/$simId",
  path: "/$simId",
  getParentRoute: () => ApiSimulationsRoute
});
const ApiSimulationsSimIdFramesFrameIdRoute = Route$1.update({
  id: "/frames/$frameId",
  path: "/frames/$frameId",
  getParentRoute: () => ApiSimulationsSimIdRoute
});
const ApiSimulationsSimIdFramesFrameIdBinRoute = Route2.update({
  id: "/bin",
  path: "/bin",
  getParentRoute: () => ApiSimulationsSimIdFramesFrameIdRoute
});
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
const rootRouteChildren = {
  IndexRoute,
  ApiSimulationsRoute: ApiSimulationsRouteWithChildren,
  VizIndexRoute
};
const routeTree = Route$6._addFileChildren(rootRouteChildren)._addFileTypes();
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
  Link as L,
  router as r
};
