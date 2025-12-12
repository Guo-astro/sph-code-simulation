import { jsx, Fragment, jsxs } from "react/jsx-runtime";
import { r as rootRouteId, i as invariant, t as trimPathLeft, j as joinPaths, d as dummyMatchContext, m as matchContext, u as useRouterState, a as useRouter, b as useForwardedRef, c as useIntersectionObserver, f as functionalUpdate, e as exactPathTest, g as removeTrailingSlash, h as deepEqual, w as warning, k as isModuleNotFoundError, R as RouterCore, E as ErrorComponent, l as json } from "../server.js";
import * as React from "react";
import React__default, { createElement, useState, useEffect, useMemo, useRef, useCallback } from "react";
import { flushSync } from "react-dom";
import { Canvas, useThree } from "@react-three/fiber";
import { PerspectiveCamera, OrbitControls, Stats } from "@react-three/drei";
import * as THREE from "three";
import { OrbitControls as OrbitControls$1 } from "three/examples/jsm/controls/OrbitControls.js";
import { ResponsiveContainer, LineChart, CartesianGrid, XAxis, YAxis, Tooltip, Legend, Line, ReferenceLine } from "recharts";
import { SkipBack, ChevronLeft, Pause, Play, ChevronRight, SkipForward, Settings, ChevronDown, RefreshCw, Folder } from "lucide-react";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const preloadWarning = "Error preloading route! ☝️";
class BaseRoute {
  constructor(options) {
    this.init = (opts) => {
      this.originalIndex = opts.originalIndex;
      const options2 = this.options;
      const isRoot = !options2?.path && !options2?.id;
      this.parentRoute = this.options.getParentRoute?.();
      if (isRoot) {
        this._path = rootRouteId;
      } else if (!this.parentRoute) {
        invariant(
          false,
          `Child Route instances must pass a 'getParentRoute: () => ParentRoute' option that returns a Route instance.`
        );
      }
      let path2 = isRoot ? rootRouteId : options2?.path;
      if (path2 && path2 !== "/") {
        path2 = trimPathLeft(path2);
      }
      const customId = options2?.id || path2;
      let id = isRoot ? rootRouteId : joinPaths([
        this.parentRoute.id === rootRouteId ? "" : this.parentRoute.id,
        customId
      ]);
      if (path2 === rootRouteId) {
        path2 = "/";
      }
      if (id !== rootRouteId) {
        id = joinPaths(["/", id]);
      }
      const fullPath = id === rootRouteId ? "/" : joinPaths([this.parentRoute.fullPath, path2]);
      this._path = path2;
      this._id = id;
      this._fullPath = fullPath;
      this._to = fullPath;
    };
    this.addChildren = (children) => {
      return this._addFileChildren(children);
    };
    this._addFileChildren = (children) => {
      if (Array.isArray(children)) {
        this.children = children;
      }
      if (typeof children === "object" && children !== null) {
        this.children = Object.values(children);
      }
      return this;
    };
    this._addFileTypes = () => {
      return this;
    };
    this.updateLoader = (options2) => {
      Object.assign(this.options, options2);
      return this;
    };
    this.update = (options2) => {
      Object.assign(this.options, options2);
      return this;
    };
    this.lazy = (lazyFn) => {
      this.lazyFn = lazyFn;
      return this;
    };
    this.options = options || {};
    this.isRoot = !options?.getParentRoute;
    if (options?.id && options?.path) {
      throw new Error(`Route cannot have both an 'id' and a 'path' option.`);
    }
  }
  get to() {
    return this._to;
  }
  get id() {
    return this._id;
  }
  get path() {
    return this._path;
  }
  get fullPath() {
    return this._fullPath;
  }
}
class BaseRootRoute extends BaseRoute {
  constructor(options) {
    super(options);
  }
}
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
const appCss = "/assets/app-DtwQlxQQ.css";
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
const $$splitComponentImporter = () => import("./index-BWIiq-yf.js");
const Route$5 = createFileRoute("/")({
  component: lazyRouteComponent($$splitComponentImporter, "component")
});
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
function createTextSprite$1(text, options = {}) {
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
  const sprite = createTextSprite$1(text, {
    color: axisColor,
    fontSize: 72,
    scale: [12, 6, 1]
    // Larger scale for axis labels
  });
  sprite.position.set(...position);
  return sprite;
}
function createTickLabel(text, position) {
  const sprite = createTextSprite$1(text, {
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
  viridis: [
    [0.267, 4e-3, 0.329],
    [0.282, 0.14, 0.458],
    [0.253, 0.265, 0.53],
    [0.206, 0.372, 0.553],
    [0.163, 0.471, 0.558],
    [0.127, 0.566, 0.551],
    [0.134, 0.658, 0.518],
    [0.266, 0.749, 0.441],
    [0.477, 0.821, 0.318],
    [0.741, 0.873, 0.15],
    [0.993, 0.906, 0.144]
  ],
  plasma: [
    [0.05, 0.03, 0.528],
    [0.295, 0.012, 0.615],
    [0.492, 0.012, 0.658],
    [0.665, 0.138, 0.618],
    [0.798, 0.28, 0.47],
    [0.899, 0.396, 0.301],
    [0.966, 0.53, 0.128],
    [0.988, 0.68, 0.063],
    [0.961, 0.85, 0.298],
    [0.94, 0.975, 0.131]
  ],
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
  turbo: [
    [0.19, 0.072, 0.232],
    [0.235, 0.318, 0.86],
    [0.137, 0.572, 0.996],
    [0.14, 0.78, 0.82],
    [0.376, 0.92, 0.512],
    [0.67, 0.979, 0.28],
    [0.924, 0.904, 0.145],
    [0.996, 0.724, 0.132],
    [0.994, 0.472, 0.122],
    [0.881, 0.2, 0.102],
    [0.528, 0.055, 0.052]
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
const hexToRgbCache = /* @__PURE__ */ new Map();
function hexToRgb$1(hex) {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!result) return { r: 1, g: 1, b: 1 };
  return {
    r: parseInt(result[1], 16) / 255,
    g: parseInt(result[2], 16) / 255,
    b: parseInt(result[3], 16) / 255
  };
}
function hexToRgbCached(hex) {
  let cached = hexToRgbCache.get(hex);
  if (!cached) {
    cached = hexToRgb$1(hex);
    hexToRgbCache.set(hex, cached);
  }
  return cached;
}
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
function interpolateColorHex$1(colors, t) {
  if (colors.length === 0) return { r: 1, g: 1, b: 1 };
  if (colors.length === 1) return hexToRgbCached(colors[0]);
  t = Math.max(0, Math.min(1, t));
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
      const color = interpolateColorHex$1(colorMap.colors, t);
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
        map: getCircleTexture(),
        alphaTest: 0.01
      }
    )
  ] });
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
function createTextSprite(text, options = {}) {
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
  const { headLength = 1, headWidth = 0.5, shaftRadius = 0.08, opacity = 1 } = options;
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
  const sunPosLabel = createTextSprite("☉ (8 kpc)", { color: "#ffff88", fontSize: 18 });
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
  const hvccLabel = createTextSprite("HVCC (~60pc from GC)", { color: "#00ff88", fontSize: 16 });
  hvccLabel.position.copy(hvccPosInGalaxy.clone().add(new THREE.Vector3(1.5 * scale, 0.5 * scale, 0)));
  hvccLabel.scale.set(5, 2, 1);
  group.add(hvccLabel);
  const legendPos = new THREE.Vector3(diskRadius * 0.7, 0.1, diskRadius * 0.5);
  const scaleNote = createTextSprite("Scale: Compressed ~300× for visualization", { color: "#888888", fontSize: 12 });
  scaleNote.position.copy(legendPos);
  scaleNote.scale.set(5, 2, 1);
  group.add(scaleNote);
  const distNote1 = createTextSprite("• Sun ↔ GC: 8 kpc (~26.7 display units)", { color: "#ffff88", fontSize: 11 });
  distNote1.position.copy(legendPos.clone().add(new THREE.Vector3(0, -1 * scale, 0)));
  distNote1.scale.set(6, 2, 1);
  group.add(distNote1);
  const distNote2 = createTextSprite("• HVCC ↔ GC: 60 pc (~0.2 display units)", { color: "#00ff88", fontSize: 11 });
  distNote2.position.copy(legendPos.clone().add(new THREE.Vector3(0, -2 * scale, 0)));
  distNote2.scale.set(6, 2, 1);
  group.add(distNote2);
  const galaxyLabel = createTextSprite("Milky Way (schematic)", { color: "#6688aa", fontSize: 20 });
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
    const gcLabel = createTextSprite(`GC (~${distanceToGC_pc} pc)`, {
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
  const observerLabel = createTextSprite("⊕ Observer (Earth)", { color: "#88ccff", fontSize: 20 });
  observerLabel.position.copy(losEnd.clone().add(new THREE.Vector3(0, 2 * scale, 0)));
  observerLabel.scale.set(5, 2.5, 1);
  group.add(observerLabel);
  const distLabel = createTextSprite("~8 kpc (not to scale)", { color: "#aaddff", fontSize: 14 });
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
  const incLabel = createTextSprite(`i = ${config.inclination ?? DEFAULT_GALACTIC_CONFIG.inclination}°`, {
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
  const zLabel = createTextSprite("L (orbital normal)", { color: "#66ffaa", fontSize: 16 });
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
  const planeLabel = createTextSprite("Orbital Plane", { color: "#aaffaa", fontSize: 14 });
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
  const vlsrLabel = createTextSprite(`V_LSR = ${lsrVelocity} km/s`, {
    color: lsrVelocity < 0 ? "#ff88ff" : "#ffcc88",
    fontSize: 16
  });
  vlsrLabel.position.copy(vlsrEndPoint.clone().add(new THREE.Vector3(1.5 * scale, 0.5 * scale, 0)));
  vlsrLabel.scale.set(4, 2, 1);
  group.add(vlsrLabel);
  const motionLabel = createTextSprite(
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
  const gcCenterLabel = createTextSprite("Galactic Center", { color: "#ffcc66", fontSize: 12 });
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
  const sunOrbitLabel = createTextSprite("☉ Sun (8 kpc from GC)", { color: "#ffff88", fontSize: 11 });
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
  const vCircularLabel = createTextSprite("V_circular = 220 km/s", {
    color: "#66ffaa",
    fontSize: 14
  });
  vCircularLabel.position.copy(vCircularArrowEnd.clone().add(new THREE.Vector3(0, 0.6 * scale, 0.5 * scale)));
  vCircularLabel.scale.set(4, 2, 1);
  vCircularGroup.add(vCircularLabel);
  const directionLabel = createTextSprite("(Galactic rotation, l = 90°)", {
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
  const periodLabel = createTextSprite("Orbital Period T ≈ 220 Myr", {
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
  const paLabel = createTextSprite(`PA = ${config.positionAngle ?? DEFAULT_GALACTIC_CONFIG.positionAngle}°`, {
    color: "#ffff88",
    fontSize: 14
  });
  paLabel.position.copy(paIndicatorGroup.position.clone().add(new THREE.Vector3(skyPlaneRadius + 2 * scale, 0, 0)));
  paLabel.scale.set(3.5, 1.8, 1);
  group.add(paLabel);
  const geometryNote = createTextSprite("Observer Geometry (Oka et al. 2017)", {
    color: "#888888",
    fontSize: 12
  });
  geometryNote.position.set(0, -3 * scale, displayDistance * 0.5);
  geometryNote.scale.set(5, 2, 1);
  group.add(geometryNote);
  return group;
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
      updateParticles();
      updateTrajectory();
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
  }, [showBlackHole, showTrajectory, showRadii, showGalacticMarkers]);
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
    computeStats();
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
  const [useImperativeMode, setUseImperativeMode] = useState(true);
  const [currentFps, setCurrentFps] = useState(0);
  const [showBlackHole, setShowBlackHole] = useState(true);
  const [showTrajectory, setShowTrajectory] = useState(true);
  const [showRadii, setShowRadii] = useState(true);
  const [showGalacticMarkers, setShowGalacticMarkers] = useState(true);
  const [showGalaxyDisk, setShowGalaxyDisk] = useState(true);
  const [showSolarSystem, setShowSolarSystem] = useState(true);
  const [animateGalaxy, setAnimateGalaxy] = useState(false);
  const [galaxyAnimationSpeed, setGalaxyAnimationSpeed] = useState(1);
  const [cameraMode, setCameraMode] = useState("free");
  const imbhPhysics = useMemo(() => {
    if (simulation?.imbhPhysics) {
      console.log("[Dashboard] Using IMBH physics from simulation config:", simulation.imbhPhysics);
      return simulation.imbhPhysics;
    }
    console.log("[Dashboard] Using default IMBH physics (CAT_OKA/A_61k/oka.json preset values)");
    return {
      enabled: true,
      // Enable by default for visualizing galactic markers
      bhPosition: [0, 0, 0],
      bhMass: 100,
      // 10^5 M_sun in code units (1000 M_sun)
      cloudInitialPosition: [20, -5.17, 0],
      // From preset
      cloudInitialVelocity: [-10.18, 5.05, 0],
      // From preset (km/s)
      cloudMass: 1,
      // 1000 M_sun in code units
      cloudRadius: 1.13,
      // From preset (pc)
      tidalRadius: 5.24,
      // From preset (pc)
      impactParameter: 5.17,
      // From preset (pc)
      pericentre: 2.217,
      // From preset (pc)
      eccentricity: 1.4504,
      // From preset (hyperbolic)
      timeUnit: 0.978
      // Myr
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
          frameIndexRef.current = 0;
          setCurrentFrameIndex(0);
          lastTime = timestamp;
        } else {
          if (framesRef.current.has(nextFrame)) {
            frameIndexRef.current = nextFrame;
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
                showRadii,
                showGalacticMarkers,
                cameraMode,
                galacticConfig: {
                  distanceToGC: 60,
                  // ~60 pc from Galactic Center
                  galacticLongitude: -0.398,
                  galacticLatitude: -0.224,
                  inclination: simulation?.imbhPhysics?.inclination ?? 70,
                  positionAngle: simulation?.imbhPhysics?.positionAngle ?? 41.6,
                  lsrVelocity: simulation?.imbhPhysics?.lsrVelocity ?? -120,
                  showGalaxyDisk,
                  showSolarSystem,
                  // Physical rotation: Sun orbits GC at V_circ = 220 km/s, period ~220 Myr
                  // For visualization: 2π rad / (220 Myr) ≈ 2.86e-17 rad/s (too slow to see)
                  // We use a visual speed: 0.1 rad/s = ~63 seconds per revolution
                  // Speed multiplier adjusts this: higher = faster animation
                  galaxyRotationSpeed: animateGalaxy ? 0.1 * galaxyAnimationSpeed : 0
                }
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
          showProjections && /* @__PURE__ */ jsxs("div", { className: "w-64 shrink-0 flex flex-col gap-1 p-1 overflow-y-auto border-l border-gray-700", children: [
            /* @__PURE__ */ jsx(
              Projection2D,
              {
                frame: currentFrame,
                projection: "xy",
                colorField,
                colorMap,
                width: 240,
                height: 150
              }
            ),
            /* @__PURE__ */ jsx(
              Projection2D,
              {
                frame: currentFrame,
                projection: "xz",
                colorField,
                colorMap,
                width: 240,
                height: 150
              }
            ),
            /* @__PURE__ */ jsx(
              Projection2D,
              {
                frame: currentFrame,
                projection: "yz",
                colorField,
                colorMap,
                width: 240,
                height: 150
              }
            )
          ] })
        ] }),
        showCharts && statistics.length > 0 && /* @__PURE__ */ jsxs("div", { className: "h-40 shrink-0 border-t border-gray-700 flex gap-2 p-1 overflow-x-auto", children: [
          /* @__PURE__ */ jsx(
            EnergyChart,
            {
              statistics,
              currentFrame: currentFrameIndex,
              className: "flex-1 min-w-[250px]"
            }
          ),
          /* @__PURE__ */ jsx(
            MomentumChart,
            {
              statistics,
              currentFrame: currentFrameIndex,
              className: "flex-1 min-w-[250px]"
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
  const imbhResultsDir = path2.join(dataRoot, "sample", "imbh_cloud", "results");
  if (fs2.existsSync(imbhResultsDir)) {
    console.log("📁 Scanning sample/imbh_cloud/results/ categories...");
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
            path2.join(dataRoot, "sample", simPath, "viz_data"),
            path2.join(dataRoot, "sample", simPath, "results", "viz_data"),
            path2.join(dataRoot, "lane_emden", "results", simPath, "viz_data"),
            path2.join(dataRoot, simPath, "viz_data")
          ];
          if (pathParts.length >= 3) {
            const [testName, ...rest] = pathParts;
            possiblePaths.unshift(
              path2.join(dataRoot, "sample", testName, "results", ...rest, "viz_data")
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
