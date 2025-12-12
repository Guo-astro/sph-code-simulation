import { ReadableStream as ReadableStream$1 } from "node:stream/web";
import { Readable, PassThrough } from "node:stream";
import { AsyncLocalStorage } from "node:async_hooks";
import { jsxs, jsx, Fragment } from "react/jsx-runtime";
import ReactDOMServer from "react-dom/server";
import * as React from "react";
import React__default, { useRef } from "react";
const __storeToDerived = /* @__PURE__ */ new WeakMap();
const __derivedToStore = /* @__PURE__ */ new WeakMap();
const __depsThatHaveWrittenThisTick = {
  current: []
};
let __isFlushing = false;
let __batchDepth = 0;
const __pendingUpdates = /* @__PURE__ */ new Set();
const __initialBatchValues = /* @__PURE__ */ new Map();
function __flush_internals(relatedVals) {
  for (const derived of relatedVals) {
    if (__depsThatHaveWrittenThisTick.current.includes(derived)) {
      continue;
    }
    __depsThatHaveWrittenThisTick.current.push(derived);
    derived.recompute();
    const stores = __derivedToStore.get(derived);
    if (stores) {
      for (const store of stores) {
        const relatedLinkedDerivedVals = __storeToDerived.get(store);
        if (!(relatedLinkedDerivedVals == null ? void 0 : relatedLinkedDerivedVals.length)) continue;
        __flush_internals(relatedLinkedDerivedVals);
      }
    }
  }
}
function __notifyListeners(store) {
  const value = {
    prevVal: store.prevState,
    currentVal: store.state
  };
  for (const listener of store.listeners) {
    listener(value);
  }
}
function __notifyDerivedListeners(derived) {
  const value = {
    prevVal: derived.prevState,
    currentVal: derived.state
  };
  for (const listener of derived.listeners) {
    listener(value);
  }
}
function __flush(store) {
  if (__batchDepth > 0 && !__initialBatchValues.has(store)) {
    __initialBatchValues.set(store, store.prevState);
  }
  __pendingUpdates.add(store);
  if (__batchDepth > 0) return;
  if (__isFlushing) return;
  try {
    __isFlushing = true;
    while (__pendingUpdates.size > 0) {
      const stores = Array.from(__pendingUpdates);
      __pendingUpdates.clear();
      for (const store2 of stores) {
        const prevState = __initialBatchValues.get(store2) ?? store2.prevState;
        store2.prevState = prevState;
        __notifyListeners(store2);
      }
      for (const store2 of stores) {
        const derivedVals = __storeToDerived.get(store2);
        if (!derivedVals) continue;
        __depsThatHaveWrittenThisTick.current.push(store2);
        __flush_internals(derivedVals);
      }
      for (const store2 of stores) {
        const derivedVals = __storeToDerived.get(store2);
        if (!derivedVals) continue;
        for (const derived of derivedVals) {
          __notifyDerivedListeners(derived);
        }
      }
    }
  } finally {
    __isFlushing = false;
    __depsThatHaveWrittenThisTick.current = [];
    __initialBatchValues.clear();
  }
}
function batch(fn2) {
  __batchDepth++;
  try {
    fn2();
  } finally {
    __batchDepth--;
    if (__batchDepth === 0) {
      const pendingUpdateToFlush = __pendingUpdates.values().next().value;
      if (pendingUpdateToFlush) {
        __flush(pendingUpdateToFlush);
      }
    }
  }
}
function isUpdaterFunction(updater) {
  return typeof updater === "function";
}
class Store {
  constructor(initialState, options) {
    this.listeners = /* @__PURE__ */ new Set();
    this.subscribe = (listener) => {
      var _a, _b;
      this.listeners.add(listener);
      const unsub = (_b = (_a = this.options) == null ? void 0 : _a.onSubscribe) == null ? void 0 : _b.call(_a, listener, this);
      return () => {
        this.listeners.delete(listener);
        unsub == null ? void 0 : unsub();
      };
    };
    this.prevState = initialState;
    this.state = initialState;
    this.options = options;
  }
  setState(updater) {
    var _a, _b, _c;
    this.prevState = this.state;
    if ((_a = this.options) == null ? void 0 : _a.updateFn) {
      this.state = this.options.updateFn(this.prevState)(updater);
    } else {
      if (isUpdaterFunction(updater)) {
        this.state = updater(this.prevState);
      } else {
        this.state = updater;
      }
    }
    (_c = (_b = this.options) == null ? void 0 : _b.onUpdate) == null ? void 0 : _c.call(_b);
    __flush(this);
  }
}
const stateIndexKey = "__TSR_index";
const popStateEvent = "popstate";
const beforeUnloadEvent = "beforeunload";
function createHistory(opts) {
  let location = opts.getLocation();
  const subscribers = /* @__PURE__ */ new Set();
  const notify = (action) => {
    location = opts.getLocation();
    subscribers.forEach((subscriber) => subscriber({ location, action }));
  };
  const handleIndexChange = (action) => {
    if (opts.notifyOnIndexChange ?? true) notify(action);
    else location = opts.getLocation();
  };
  const tryNavigation = async ({
    task,
    navigateOpts,
    ...actionInfo
  }) => {
    const ignoreBlocker = navigateOpts?.ignoreBlocker ?? false;
    if (ignoreBlocker) {
      task();
      return;
    }
    const blockers = opts.getBlockers?.() ?? [];
    const isPushOrReplace = actionInfo.type === "PUSH" || actionInfo.type === "REPLACE";
    if (typeof document !== "undefined" && blockers.length && isPushOrReplace) {
      for (const blocker of blockers) {
        const nextLocation = parseHref(actionInfo.path, actionInfo.state);
        const isBlocked = await blocker.blockerFn({
          currentLocation: location,
          nextLocation,
          action: actionInfo.type
        });
        if (isBlocked) {
          opts.onBlocked?.();
          return;
        }
      }
    }
    task();
  };
  return {
    get location() {
      return location;
    },
    get length() {
      return opts.getLength();
    },
    subscribers,
    subscribe: (cb) => {
      subscribers.add(cb);
      return () => {
        subscribers.delete(cb);
      };
    },
    push: (path, state, navigateOpts) => {
      const currentIndex = location.state[stateIndexKey];
      state = assignKeyAndIndex(currentIndex + 1, state);
      tryNavigation({
        task: () => {
          opts.pushState(path, state);
          notify({ type: "PUSH" });
        },
        navigateOpts,
        type: "PUSH",
        path,
        state
      });
    },
    replace: (path, state, navigateOpts) => {
      const currentIndex = location.state[stateIndexKey];
      state = assignKeyAndIndex(currentIndex, state);
      tryNavigation({
        task: () => {
          opts.replaceState(path, state);
          notify({ type: "REPLACE" });
        },
        navigateOpts,
        type: "REPLACE",
        path,
        state
      });
    },
    go: (index, navigateOpts) => {
      tryNavigation({
        task: () => {
          opts.go(index);
          handleIndexChange({ type: "GO", index });
        },
        navigateOpts,
        type: "GO"
      });
    },
    back: (navigateOpts) => {
      tryNavigation({
        task: () => {
          opts.back(navigateOpts?.ignoreBlocker ?? false);
          handleIndexChange({ type: "BACK" });
        },
        navigateOpts,
        type: "BACK"
      });
    },
    forward: (navigateOpts) => {
      tryNavigation({
        task: () => {
          opts.forward(navigateOpts?.ignoreBlocker ?? false);
          handleIndexChange({ type: "FORWARD" });
        },
        navigateOpts,
        type: "FORWARD"
      });
    },
    canGoBack: () => location.state[stateIndexKey] !== 0,
    createHref: (str) => opts.createHref(str),
    block: (blocker) => {
      if (!opts.setBlockers) return () => {
      };
      const blockers = opts.getBlockers?.() ?? [];
      opts.setBlockers([...blockers, blocker]);
      return () => {
        const blockers2 = opts.getBlockers?.() ?? [];
        opts.setBlockers?.(blockers2.filter((b2) => b2 !== blocker));
      };
    },
    flush: () => opts.flush?.(),
    destroy: () => opts.destroy?.(),
    notify
  };
}
function assignKeyAndIndex(index, state) {
  if (!state) {
    state = {};
  }
  const key = createRandomKey();
  return {
    ...state,
    key,
    // TODO: Remove in v2 - use __TSR_key instead
    __TSR_key: key,
    [stateIndexKey]: index
  };
}
function createBrowserHistory(opts) {
  const win = typeof document !== "undefined" ? window : void 0;
  const originalPushState = win.history.pushState;
  const originalReplaceState = win.history.replaceState;
  let blockers = [];
  const _getBlockers = () => blockers;
  const _setBlockers = (newBlockers) => blockers = newBlockers;
  const createHref = ((path) => path);
  const parseLocation = (() => parseHref(
    `${win.location.pathname}${win.location.search}${win.location.hash}`,
    win.history.state
  ));
  if (!win.history.state?.__TSR_key && !win.history.state?.key) {
    const addedKey = createRandomKey();
    win.history.replaceState(
      {
        [stateIndexKey]: 0,
        key: addedKey,
        // TODO: Remove in v2 - use __TSR_key instead
        __TSR_key: addedKey
      },
      ""
    );
  }
  let currentLocation = parseLocation();
  let rollbackLocation;
  let nextPopIsGo = false;
  let ignoreNextPop = false;
  let skipBlockerNextPop = false;
  let ignoreNextBeforeUnload = false;
  const getLocation = () => currentLocation;
  let next;
  let scheduled;
  const flush = () => {
    if (!next) {
      return;
    }
    history._ignoreSubscribers = true;
    (next.isPush ? win.history.pushState : win.history.replaceState)(
      next.state,
      "",
      next.href
    );
    history._ignoreSubscribers = false;
    next = void 0;
    scheduled = void 0;
    rollbackLocation = void 0;
  };
  const queueHistoryAction = (type, destHref, state) => {
    const href = createHref(destHref);
    if (!scheduled) {
      rollbackLocation = currentLocation;
    }
    currentLocation = parseHref(destHref, state);
    next = {
      href,
      state,
      isPush: next?.isPush || type === "push"
    };
    if (!scheduled) {
      scheduled = Promise.resolve().then(() => flush());
    }
  };
  const onPushPop = (type) => {
    currentLocation = parseLocation();
    history.notify({ type });
  };
  const onPushPopEvent = async () => {
    if (ignoreNextPop) {
      ignoreNextPop = false;
      return;
    }
    const nextLocation = parseLocation();
    const delta = nextLocation.state[stateIndexKey] - currentLocation.state[stateIndexKey];
    const isForward = delta === 1;
    const isBack = delta === -1;
    const isGo = !isForward && !isBack || nextPopIsGo;
    nextPopIsGo = false;
    const action = isGo ? "GO" : isBack ? "BACK" : "FORWARD";
    const notify = isGo ? {
      type: "GO",
      index: delta
    } : {
      type: isBack ? "BACK" : "FORWARD"
    };
    if (skipBlockerNextPop) {
      skipBlockerNextPop = false;
    } else {
      const blockers2 = _getBlockers();
      if (typeof document !== "undefined" && blockers2.length) {
        for (const blocker of blockers2) {
          const isBlocked = await blocker.blockerFn({
            currentLocation,
            nextLocation,
            action
          });
          if (isBlocked) {
            ignoreNextPop = true;
            win.history.go(1);
            history.notify(notify);
            return;
          }
        }
      }
    }
    currentLocation = parseLocation();
    history.notify(notify);
  };
  const onBeforeUnload = (e) => {
    if (ignoreNextBeforeUnload) {
      ignoreNextBeforeUnload = false;
      return;
    }
    let shouldBlock = false;
    const blockers2 = _getBlockers();
    if (typeof document !== "undefined" && blockers2.length) {
      for (const blocker of blockers2) {
        const shouldHaveBeforeUnload = blocker.enableBeforeUnload ?? true;
        if (shouldHaveBeforeUnload === true) {
          shouldBlock = true;
          break;
        }
        if (typeof shouldHaveBeforeUnload === "function" && shouldHaveBeforeUnload() === true) {
          shouldBlock = true;
          break;
        }
      }
    }
    if (shouldBlock) {
      e.preventDefault();
      return e.returnValue = "";
    }
    return;
  };
  const history = createHistory({
    getLocation,
    getLength: () => win.history.length,
    pushState: (href, state) => queueHistoryAction("push", href, state),
    replaceState: (href, state) => queueHistoryAction("replace", href, state),
    back: (ignoreBlocker) => {
      if (ignoreBlocker) skipBlockerNextPop = true;
      ignoreNextBeforeUnload = true;
      return win.history.back();
    },
    forward: (ignoreBlocker) => {
      if (ignoreBlocker) skipBlockerNextPop = true;
      ignoreNextBeforeUnload = true;
      win.history.forward();
    },
    go: (n) => {
      nextPopIsGo = true;
      win.history.go(n);
    },
    createHref: (href) => createHref(href),
    flush,
    destroy: () => {
      win.history.pushState = originalPushState;
      win.history.replaceState = originalReplaceState;
      win.removeEventListener(beforeUnloadEvent, onBeforeUnload, {
        capture: true
      });
      win.removeEventListener(popStateEvent, onPushPopEvent);
    },
    onBlocked: () => {
      if (rollbackLocation && currentLocation !== rollbackLocation) {
        currentLocation = rollbackLocation;
      }
    },
    getBlockers: _getBlockers,
    setBlockers: _setBlockers,
    notifyOnIndexChange: false
  });
  win.addEventListener(beforeUnloadEvent, onBeforeUnload, { capture: true });
  win.addEventListener(popStateEvent, onPushPopEvent);
  win.history.pushState = function(...args) {
    const res = originalPushState.apply(win.history, args);
    if (!history._ignoreSubscribers) onPushPop("PUSH");
    return res;
  };
  win.history.replaceState = function(...args) {
    const res = originalReplaceState.apply(win.history, args);
    if (!history._ignoreSubscribers) onPushPop("REPLACE");
    return res;
  };
  return history;
}
function createMemoryHistory(opts = {
  initialEntries: ["/"]
}) {
  const entries = opts.initialEntries;
  let index = opts.initialIndex ? Math.min(Math.max(opts.initialIndex, 0), entries.length - 1) : entries.length - 1;
  const states = entries.map(
    (_entry, index2) => assignKeyAndIndex(index2, void 0)
  );
  const getLocation = () => parseHref(entries[index], states[index]);
  return createHistory({
    getLocation,
    getLength: () => entries.length,
    pushState: (path, state) => {
      if (index < entries.length - 1) {
        entries.splice(index + 1);
        states.splice(index + 1);
      }
      states.push(state);
      entries.push(path);
      index = Math.max(entries.length - 1, 0);
    },
    replaceState: (path, state) => {
      states[index] = state;
      entries[index] = path;
    },
    back: () => {
      index = Math.max(index - 1, 0);
    },
    forward: () => {
      index = Math.min(index + 1, entries.length - 1);
    },
    go: (n) => {
      index = Math.min(Math.max(index + n, 0), entries.length - 1);
    },
    createHref: (path) => path
  });
}
function parseHref(href, state) {
  const hashIndex = href.indexOf("#");
  const searchIndex = href.indexOf("?");
  const addedKey = createRandomKey();
  return {
    href,
    pathname: href.substring(
      0,
      hashIndex > 0 ? searchIndex > 0 ? Math.min(hashIndex, searchIndex) : hashIndex : searchIndex > 0 ? searchIndex : href.length
    ),
    hash: hashIndex > -1 ? href.substring(hashIndex) : "",
    search: searchIndex > -1 ? href.slice(searchIndex, hashIndex === -1 ? void 0 : hashIndex) : "",
    state: state || { [stateIndexKey]: 0, key: addedKey, __TSR_key: addedKey }
  };
}
function createRandomKey() {
  return (Math.random() + 1).toString(36).substring(7);
}
function last(arr) {
  return arr[arr.length - 1];
}
function isFunction(d) {
  return typeof d === "function";
}
function functionalUpdate(updater, previous) {
  if (isFunction(updater)) {
    return updater(previous);
  }
  return updater;
}
const hasOwn = Object.prototype.hasOwnProperty;
function replaceEqualDeep(prev, _next) {
  if (prev === _next) {
    return prev;
  }
  const next = _next;
  const array = isPlainArray(prev) && isPlainArray(next);
  if (!array && !(isPlainObject(prev) && isPlainObject(next))) return next;
  const prevItems = array ? prev : getEnumerableOwnKeys(prev);
  if (!prevItems) return next;
  const nextItems = array ? next : getEnumerableOwnKeys(next);
  if (!nextItems) return next;
  const prevSize = prevItems.length;
  const nextSize = nextItems.length;
  const copy = array ? new Array(nextSize) : {};
  let equalItems = 0;
  for (let i = 0; i < nextSize; i++) {
    const key = array ? i : nextItems[i];
    const p2 = prev[key];
    const n = next[key];
    if (p2 === n) {
      copy[key] = p2;
      if (array ? i < prevSize : hasOwn.call(prev, key)) equalItems++;
      continue;
    }
    if (p2 === null || n === null || typeof p2 !== "object" || typeof n !== "object") {
      copy[key] = n;
      continue;
    }
    const v2 = replaceEqualDeep(p2, n);
    copy[key] = v2;
    if (v2 === p2) equalItems++;
  }
  return prevSize === nextSize && equalItems === prevSize ? prev : copy;
}
function getEnumerableOwnKeys(o2) {
  const keys = [];
  const names = Object.getOwnPropertyNames(o2);
  for (const name of names) {
    if (!Object.prototype.propertyIsEnumerable.call(o2, name)) return false;
    keys.push(name);
  }
  const symbols = Object.getOwnPropertySymbols(o2);
  for (const symbol of symbols) {
    if (!Object.prototype.propertyIsEnumerable.call(o2, symbol)) return false;
    keys.push(symbol);
  }
  return keys;
}
function isPlainObject(o2) {
  if (!hasObjectPrototype(o2)) {
    return false;
  }
  const ctor = o2.constructor;
  if (typeof ctor === "undefined") {
    return true;
  }
  const prot = ctor.prototype;
  if (!hasObjectPrototype(prot)) {
    return false;
  }
  if (!prot.hasOwnProperty("isPrototypeOf")) {
    return false;
  }
  return true;
}
function hasObjectPrototype(o2) {
  return Object.prototype.toString.call(o2) === "[object Object]";
}
function isPlainArray(value) {
  return Array.isArray(value) && value.length === Object.keys(value).length;
}
function deepEqual(a, b2, opts) {
  if (a === b2) {
    return true;
  }
  if (typeof a !== typeof b2) {
    return false;
  }
  if (Array.isArray(a) && Array.isArray(b2)) {
    if (a.length !== b2.length) return false;
    for (let i = 0, l = a.length; i < l; i++) {
      if (!deepEqual(a[i], b2[i], opts)) return false;
    }
    return true;
  }
  if (isPlainObject(a) && isPlainObject(b2)) {
    const ignoreUndefined = opts?.ignoreUndefined ?? true;
    if (opts?.partial) {
      for (const k2 in b2) {
        if (!ignoreUndefined || b2[k2] !== void 0) {
          if (!deepEqual(a[k2], b2[k2], opts)) return false;
        }
      }
      return true;
    }
    let aCount = 0;
    if (!ignoreUndefined) {
      aCount = Object.keys(a).length;
    } else {
      for (const k2 in a) {
        if (a[k2] !== void 0) aCount++;
      }
    }
    let bCount = 0;
    for (const k2 in b2) {
      if (!ignoreUndefined || b2[k2] !== void 0) {
        bCount++;
        if (bCount > aCount || !deepEqual(a[k2], b2[k2], opts)) return false;
      }
    }
    return aCount === bCount;
  }
  return false;
}
function createControlledPromise(onResolve) {
  let resolveLoadPromise;
  let rejectLoadPromise;
  const controlledPromise = new Promise((resolve, reject) => {
    resolveLoadPromise = resolve;
    rejectLoadPromise = reject;
  });
  controlledPromise.status = "pending";
  controlledPromise.resolve = (value) => {
    controlledPromise.status = "resolved";
    controlledPromise.value = value;
    resolveLoadPromise(value);
    onResolve?.(value);
  };
  controlledPromise.reject = (e) => {
    controlledPromise.status = "rejected";
    rejectLoadPromise(e);
  };
  return controlledPromise;
}
function isModuleNotFoundError(error) {
  if (typeof error?.message !== "string") return false;
  return error.message.startsWith("Failed to fetch dynamically imported module") || error.message.startsWith("error loading dynamically imported module") || error.message.startsWith("Importing a module script failed");
}
function isPromise(value) {
  return Boolean(
    value && typeof value === "object" && typeof value.then === "function"
  );
}
function findLast(array, predicate) {
  for (let i = array.length - 1; i >= 0; i--) {
    const item = array[i];
    if (predicate(item)) return item;
  }
  return void 0;
}
function decodeSegment(segment) {
  try {
    return decodeURI(segment);
  } catch {
    return segment.replaceAll(/%[0-9A-F]{2}/gi, (match) => {
      try {
        return decodeURI(match);
      } catch {
        return match;
      }
    });
  }
}
function decodePath(path, decodeIgnore) {
  if (!path) return path;
  const re2 = /%25|%5C/gi;
  let cursor = 0;
  let result = "";
  let match;
  while (null !== (match = re2.exec(path))) {
    result += decodeSegment(path.slice(cursor, match.index)) + match[0];
    cursor = re2.lastIndex;
  }
  return result + decodeSegment(cursor ? path.slice(cursor) : path);
}
var isProduction$1 = process.env.NODE_ENV === "production";
var prefix = "Invariant failed";
function invariant(condition, message) {
  if (condition) {
    return;
  }
  if (isProduction$1) {
    throw new Error(prefix);
  }
  var provided = typeof message === "function" ? message() : message;
  var value = provided ? "".concat(prefix, ": ").concat(provided) : prefix;
  throw new Error(value);
}
function createLRUCache(max) {
  const cache = /* @__PURE__ */ new Map();
  let oldest;
  let newest;
  const touch = (entry) => {
    if (!entry.next) return;
    if (!entry.prev) {
      entry.next.prev = void 0;
      oldest = entry.next;
      entry.next = void 0;
      if (newest) {
        entry.prev = newest;
        newest.next = entry;
      }
    } else {
      entry.prev.next = entry.next;
      entry.next.prev = entry.prev;
      entry.next = void 0;
      if (newest) {
        newest.next = entry;
        entry.prev = newest;
      }
    }
    newest = entry;
  };
  return {
    get(key) {
      const entry = cache.get(key);
      if (!entry) return void 0;
      touch(entry);
      return entry.value;
    },
    set(key, value) {
      if (cache.size >= max && oldest) {
        const toDelete = oldest;
        cache.delete(toDelete.key);
        if (toDelete.next) {
          oldest = toDelete.next;
          toDelete.next.prev = void 0;
        }
        if (toDelete === newest) {
          newest = void 0;
        }
      }
      const existing = cache.get(key);
      if (existing) {
        existing.value = value;
        touch(existing);
      } else {
        const entry = { key, value, prev: newest };
        if (newest) newest.next = entry;
        newest = entry;
        if (!oldest) oldest = entry;
        cache.set(key, entry);
      }
    },
    clear() {
      cache.clear();
      oldest = void 0;
      newest = void 0;
    }
  };
}
const SEGMENT_TYPE_PATHNAME = 0;
const SEGMENT_TYPE_PARAM = 1;
const SEGMENT_TYPE_WILDCARD = 2;
const SEGMENT_TYPE_OPTIONAL_PARAM = 3;
const PARAM_W_CURLY_BRACES_RE = /^([^{]*)\{\$([a-zA-Z_$][a-zA-Z0-9_$]*)\}([^}]*)$/;
const OPTIONAL_PARAM_W_CURLY_BRACES_RE = /^([^{]*)\{-\$([a-zA-Z_$][a-zA-Z0-9_$]*)\}([^}]*)$/;
const WILDCARD_W_CURLY_BRACES_RE = /^([^{]*)\{\$\}([^}]*)$/;
function parseSegment(path, start, output = new Uint16Array(6)) {
  const next = path.indexOf("/", start);
  const end = next === -1 ? path.length : next;
  const part = path.substring(start, end);
  if (!part || !part.includes("$")) {
    output[0] = SEGMENT_TYPE_PATHNAME;
    output[1] = start;
    output[2] = start;
    output[3] = end;
    output[4] = end;
    output[5] = end;
    return output;
  }
  if (part === "$") {
    const total = path.length;
    output[0] = SEGMENT_TYPE_WILDCARD;
    output[1] = start;
    output[2] = start;
    output[3] = total;
    output[4] = total;
    output[5] = total;
    return output;
  }
  if (part.charCodeAt(0) === 36) {
    output[0] = SEGMENT_TYPE_PARAM;
    output[1] = start;
    output[2] = start + 1;
    output[3] = end;
    output[4] = end;
    output[5] = end;
    return output;
  }
  const wildcardBracesMatch = part.match(WILDCARD_W_CURLY_BRACES_RE);
  if (wildcardBracesMatch) {
    const prefix2 = wildcardBracesMatch[1];
    const pLength = prefix2.length;
    output[0] = SEGMENT_TYPE_WILDCARD;
    output[1] = start + pLength;
    output[2] = start + pLength + 1;
    output[3] = start + pLength + 2;
    output[4] = start + pLength + 3;
    output[5] = path.length;
    return output;
  }
  const optionalParamBracesMatch = part.match(OPTIONAL_PARAM_W_CURLY_BRACES_RE);
  if (optionalParamBracesMatch) {
    const prefix2 = optionalParamBracesMatch[1];
    const paramName = optionalParamBracesMatch[2];
    const suffix = optionalParamBracesMatch[3];
    const pLength = prefix2.length;
    output[0] = SEGMENT_TYPE_OPTIONAL_PARAM;
    output[1] = start + pLength;
    output[2] = start + pLength + 3;
    output[3] = start + pLength + 3 + paramName.length;
    output[4] = end - suffix.length;
    output[5] = end;
    return output;
  }
  const paramBracesMatch = part.match(PARAM_W_CURLY_BRACES_RE);
  if (paramBracesMatch) {
    const prefix2 = paramBracesMatch[1];
    const paramName = paramBracesMatch[2];
    const suffix = paramBracesMatch[3];
    const pLength = prefix2.length;
    output[0] = SEGMENT_TYPE_PARAM;
    output[1] = start + pLength;
    output[2] = start + pLength + 2;
    output[3] = start + pLength + 2 + paramName.length;
    output[4] = end - suffix.length;
    output[5] = end;
    return output;
  }
  output[0] = SEGMENT_TYPE_PATHNAME;
  output[1] = start;
  output[2] = start;
  output[3] = end;
  output[4] = end;
  output[5] = end;
  return output;
}
function parseSegments(defaultCaseSensitive, data, route, start, node, depth, onRoute) {
  onRoute?.(route);
  let cursor = start;
  {
    const path = route.fullPath ?? route.from;
    const length = path.length;
    const caseSensitive = route.options?.caseSensitive ?? defaultCaseSensitive;
    while (cursor < length) {
      const segment = parseSegment(path, cursor, data);
      let nextNode;
      const start2 = cursor;
      const end = segment[5];
      cursor = end + 1;
      depth++;
      const kind = segment[0];
      switch (kind) {
        case SEGMENT_TYPE_PATHNAME: {
          const value = path.substring(segment[2], segment[3]);
          if (caseSensitive) {
            const existingNode = node.static?.get(value);
            if (existingNode) {
              nextNode = existingNode;
            } else {
              node.static ??= /* @__PURE__ */ new Map();
              const next = createStaticNode(
                route.fullPath ?? route.from
              );
              next.parent = node;
              next.depth = depth;
              nextNode = next;
              node.static.set(value, next);
            }
          } else {
            const name = value.toLowerCase();
            const existingNode = node.staticInsensitive?.get(name);
            if (existingNode) {
              nextNode = existingNode;
            } else {
              node.staticInsensitive ??= /* @__PURE__ */ new Map();
              const next = createStaticNode(
                route.fullPath ?? route.from
              );
              next.parent = node;
              next.depth = depth;
              nextNode = next;
              node.staticInsensitive.set(name, next);
            }
          }
          break;
        }
        case SEGMENT_TYPE_PARAM: {
          const prefix_raw = path.substring(start2, segment[1]);
          const suffix_raw = path.substring(segment[4], end);
          const actuallyCaseSensitive = caseSensitive && !!(prefix_raw || suffix_raw);
          const prefix2 = !prefix_raw ? void 0 : actuallyCaseSensitive ? prefix_raw : prefix_raw.toLowerCase();
          const suffix = !suffix_raw ? void 0 : actuallyCaseSensitive ? suffix_raw : suffix_raw.toLowerCase();
          const existingNode = node.dynamic?.find(
            (s) => s.caseSensitive === actuallyCaseSensitive && s.prefix === prefix2 && s.suffix === suffix
          );
          if (existingNode) {
            nextNode = existingNode;
          } else {
            const next = createDynamicNode(
              SEGMENT_TYPE_PARAM,
              route.fullPath ?? route.from,
              actuallyCaseSensitive,
              prefix2,
              suffix
            );
            nextNode = next;
            next.depth = depth;
            next.parent = node;
            node.dynamic ??= [];
            node.dynamic.push(next);
          }
          break;
        }
        case SEGMENT_TYPE_OPTIONAL_PARAM: {
          const prefix_raw = path.substring(start2, segment[1]);
          const suffix_raw = path.substring(segment[4], end);
          const actuallyCaseSensitive = caseSensitive && !!(prefix_raw || suffix_raw);
          const prefix2 = !prefix_raw ? void 0 : actuallyCaseSensitive ? prefix_raw : prefix_raw.toLowerCase();
          const suffix = !suffix_raw ? void 0 : actuallyCaseSensitive ? suffix_raw : suffix_raw.toLowerCase();
          const existingNode = node.optional?.find(
            (s) => s.caseSensitive === actuallyCaseSensitive && s.prefix === prefix2 && s.suffix === suffix
          );
          if (existingNode) {
            nextNode = existingNode;
          } else {
            const next = createDynamicNode(
              SEGMENT_TYPE_OPTIONAL_PARAM,
              route.fullPath ?? route.from,
              actuallyCaseSensitive,
              prefix2,
              suffix
            );
            nextNode = next;
            next.parent = node;
            next.depth = depth;
            node.optional ??= [];
            node.optional.push(next);
          }
          break;
        }
        case SEGMENT_TYPE_WILDCARD: {
          const prefix_raw = path.substring(start2, segment[1]);
          const suffix_raw = path.substring(segment[4], end);
          const actuallyCaseSensitive = caseSensitive && !!(prefix_raw || suffix_raw);
          const prefix2 = !prefix_raw ? void 0 : actuallyCaseSensitive ? prefix_raw : prefix_raw.toLowerCase();
          const suffix = !suffix_raw ? void 0 : actuallyCaseSensitive ? suffix_raw : suffix_raw.toLowerCase();
          const next = createDynamicNode(
            SEGMENT_TYPE_WILDCARD,
            route.fullPath ?? route.from,
            actuallyCaseSensitive,
            prefix2,
            suffix
          );
          nextNode = next;
          next.parent = node;
          next.depth = depth;
          node.wildcard ??= [];
          node.wildcard.push(next);
        }
      }
      node = nextNode;
    }
    if ((route.path || !route.children) && !route.isRoot) {
      const isIndex = path.endsWith("/");
      if (!isIndex) node.notFound = route;
      if (!node.route || !node.isIndex && isIndex) {
        node.route = route;
        node.fullPath = route.fullPath ?? route.from;
      }
      node.isIndex ||= isIndex;
    }
  }
  if (route.children)
    for (const child of route.children) {
      parseSegments(
        defaultCaseSensitive,
        data,
        child,
        cursor,
        node,
        depth,
        onRoute
      );
    }
}
function sortDynamic(a, b2) {
  if (a.prefix && b2.prefix && a.prefix !== b2.prefix) {
    if (a.prefix.startsWith(b2.prefix)) return -1;
    if (b2.prefix.startsWith(a.prefix)) return 1;
  }
  if (a.suffix && b2.suffix && a.suffix !== b2.suffix) {
    if (a.suffix.endsWith(b2.suffix)) return -1;
    if (b2.suffix.endsWith(a.suffix)) return 1;
  }
  if (a.prefix && !b2.prefix) return -1;
  if (!a.prefix && b2.prefix) return 1;
  if (a.suffix && !b2.suffix) return -1;
  if (!a.suffix && b2.suffix) return 1;
  if (a.caseSensitive && !b2.caseSensitive) return -1;
  if (!a.caseSensitive && b2.caseSensitive) return 1;
  return 0;
}
function sortTreeNodes(node) {
  if (node.static) {
    for (const child of node.static.values()) {
      sortTreeNodes(child);
    }
  }
  if (node.staticInsensitive) {
    for (const child of node.staticInsensitive.values()) {
      sortTreeNodes(child);
    }
  }
  if (node.dynamic?.length) {
    node.dynamic.sort(sortDynamic);
    for (const child of node.dynamic) {
      sortTreeNodes(child);
    }
  }
  if (node.optional?.length) {
    node.optional.sort(sortDynamic);
    for (const child of node.optional) {
      sortTreeNodes(child);
    }
  }
  if (node.wildcard?.length) {
    node.wildcard.sort(sortDynamic);
    for (const child of node.wildcard) {
      sortTreeNodes(child);
    }
  }
}
function createStaticNode(fullPath) {
  return {
    kind: SEGMENT_TYPE_PATHNAME,
    depth: 0,
    static: null,
    staticInsensitive: null,
    dynamic: null,
    optional: null,
    wildcard: null,
    route: null,
    fullPath,
    parent: null,
    isIndex: false,
    notFound: null
  };
}
function createDynamicNode(kind, fullPath, caseSensitive, prefix2, suffix) {
  return {
    kind,
    depth: 0,
    static: null,
    staticInsensitive: null,
    dynamic: null,
    optional: null,
    wildcard: null,
    route: null,
    fullPath,
    parent: null,
    isIndex: false,
    notFound: null,
    caseSensitive,
    prefix: prefix2,
    suffix
  };
}
function processRouteMasks(routeList, processedTree) {
  const segmentTree = createStaticNode("/");
  const data = new Uint16Array(6);
  for (const route of routeList) {
    parseSegments(false, data, route, 1, segmentTree, 0);
  }
  sortTreeNodes(segmentTree);
  processedTree.masksTree = segmentTree;
  processedTree.flatCache = createLRUCache(1e3);
}
function findFlatMatch(path, processedTree) {
  path ||= "/";
  const cached = processedTree.flatCache.get(path);
  if (cached) return cached;
  const result = findMatch(path, processedTree.masksTree);
  processedTree.flatCache.set(path, result);
  return result;
}
function findSingleMatch(from, caseSensitive, fuzzy, path, processedTree) {
  from ||= "/";
  path ||= "/";
  const key = caseSensitive ? `case\0${from}` : from;
  let tree = processedTree.singleCache.get(key);
  if (!tree) {
    tree = createStaticNode("/");
    const data = new Uint16Array(6);
    parseSegments(caseSensitive, data, { from }, 1, tree, 0);
    processedTree.singleCache.set(key, tree);
  }
  return findMatch(path, tree, fuzzy);
}
function findRouteMatch(path, processedTree, fuzzy = false) {
  const key = fuzzy ? path : `nofuzz\0${path}`;
  const cached = processedTree.matchCache.get(key);
  if (cached !== void 0) return cached;
  path ||= "/";
  const result = findMatch(
    path,
    processedTree.segmentTree,
    fuzzy
  );
  if (result) result.branch = buildRouteBranch(result.route);
  processedTree.matchCache.set(key, result);
  return result;
}
function trimPathRight$1(path) {
  return path === "/" ? path : path.replace(/\/{1,}$/, "");
}
function processRouteTree(routeTree, caseSensitive = false, initRoute) {
  const segmentTree = createStaticNode(routeTree.fullPath);
  const data = new Uint16Array(6);
  const routesById = {};
  const routesByPath = {};
  let index = 0;
  parseSegments(caseSensitive, data, routeTree, 1, segmentTree, 0, (route) => {
    initRoute?.(route, index);
    invariant(
      !(route.id in routesById),
      `Duplicate routes found with id: ${String(route.id)}`
    );
    routesById[route.id] = route;
    if (index !== 0 && route.path) {
      const trimmedFullPath = trimPathRight$1(route.fullPath);
      if (!routesByPath[trimmedFullPath] || route.fullPath.endsWith("/")) {
        routesByPath[trimmedFullPath] = route;
      }
    }
    index++;
  });
  sortTreeNodes(segmentTree);
  const processedTree = {
    segmentTree,
    singleCache: createLRUCache(1e3),
    matchCache: createLRUCache(1e3),
    flatCache: null,
    masksTree: null
  };
  return {
    processedTree,
    routesById,
    routesByPath
  };
}
function findMatch(path, segmentTree, fuzzy = false) {
  const parts = path.split("/");
  const leaf = getNodeMatch(path, parts, segmentTree, fuzzy);
  if (!leaf) return null;
  const params = extractParams(path, parts, leaf);
  const isFuzzyMatch = "**" in leaf;
  if (isFuzzyMatch) params["**"] = leaf["**"];
  const route = isFuzzyMatch ? leaf.node.notFound ?? leaf.node.route : leaf.node.route;
  return {
    route,
    params
  };
}
function extractParams(path, parts, leaf) {
  const list = buildBranch(leaf.node);
  let nodeParts = null;
  const params = {};
  for (let partIndex = 0, nodeIndex = 0, pathIndex = 0; nodeIndex < list.length; partIndex++, nodeIndex++, pathIndex++) {
    const node = list[nodeIndex];
    const part = parts[partIndex];
    const currentPathIndex = pathIndex;
    if (part) pathIndex += part.length;
    if (node.kind === SEGMENT_TYPE_PARAM) {
      nodeParts ??= leaf.node.fullPath.split("/");
      const nodePart = nodeParts[nodeIndex];
      const preLength = node.prefix?.length ?? 0;
      const isCurlyBraced = nodePart.charCodeAt(preLength) === 123;
      if (isCurlyBraced) {
        const sufLength = node.suffix?.length ?? 0;
        const name = nodePart.substring(
          preLength + 2,
          nodePart.length - sufLength - 1
        );
        const value = part.substring(preLength, part.length - sufLength);
        params[name] = decodeURIComponent(value);
      } else {
        const name = nodePart.substring(1);
        params[name] = decodeURIComponent(part);
      }
    } else if (node.kind === SEGMENT_TYPE_OPTIONAL_PARAM) {
      if (leaf.skipped & 1 << nodeIndex) {
        partIndex--;
        continue;
      }
      nodeParts ??= leaf.node.fullPath.split("/");
      const nodePart = nodeParts[nodeIndex];
      const preLength = node.prefix?.length ?? 0;
      const sufLength = node.suffix?.length ?? 0;
      const name = nodePart.substring(
        preLength + 3,
        nodePart.length - sufLength - 1
      );
      const value = node.suffix || node.prefix ? part.substring(preLength, part.length - sufLength) : part;
      if (value) params[name] = decodeURIComponent(value);
    } else if (node.kind === SEGMENT_TYPE_WILDCARD) {
      const n = node;
      const value = path.substring(
        currentPathIndex + (n.prefix?.length ?? 0),
        path.length - (n.suffix?.length ?? 0)
      );
      const splat = decodeURIComponent(value);
      params["*"] = splat;
      params._splat = splat;
      break;
    }
  }
  return params;
}
function buildRouteBranch(route) {
  const list = [route];
  while (route.parentRoute) {
    route = route.parentRoute;
    list.push(route);
  }
  list.reverse();
  return list;
}
function buildBranch(node) {
  const list = Array(node.depth + 1);
  do {
    list[node.depth] = node;
    node = node.parent;
  } while (node);
  return list;
}
function getNodeMatch(path, parts, segmentTree, fuzzy) {
  const trailingSlash = !last(parts);
  const pathIsIndex = trailingSlash && path !== "/";
  const partsLength = parts.length - (trailingSlash ? 1 : 0);
  const stack = [
    {
      node: segmentTree,
      index: 1,
      skipped: 0,
      depth: 1,
      statics: 1,
      dynamics: 0,
      optionals: 0
    }
  ];
  let wildcardMatch = null;
  let bestFuzzy = null;
  let bestMatch = null;
  while (stack.length) {
    const frame = stack.pop();
    let { node, index, skipped, depth, statics, dynamics, optionals } = frame;
    if (fuzzy && node.notFound && isFrameMoreSpecific(bestFuzzy, frame)) {
      bestFuzzy = frame;
    }
    const isBeyondPath = index === partsLength;
    if (isBeyondPath) {
      if (node.route && (!pathIsIndex || node.isIndex)) {
        if (isFrameMoreSpecific(bestMatch, frame)) {
          bestMatch = frame;
        }
        if (statics === partsLength && node.isIndex) return bestMatch;
      }
      if (!node.optional && !node.wildcard) continue;
    }
    const part = isBeyondPath ? void 0 : parts[index];
    let lowerPart;
    if (node.wildcard && isFrameMoreSpecific(wildcardMatch, frame)) {
      for (const segment of node.wildcard) {
        const { prefix: prefix2, suffix } = segment;
        if (prefix2) {
          if (isBeyondPath) continue;
          const casePart = segment.caseSensitive ? part : lowerPart ??= part.toLowerCase();
          if (!casePart.startsWith(prefix2)) continue;
        }
        if (suffix) {
          if (isBeyondPath) continue;
          const end = parts.slice(index).join("/").slice(-suffix.length);
          const casePart = segment.caseSensitive ? end : end.toLowerCase();
          if (casePart !== suffix) continue;
        }
        wildcardMatch = {
          node: segment,
          index,
          skipped,
          depth,
          statics,
          dynamics,
          optionals
        };
        break;
      }
    }
    if (node.optional) {
      const nextSkipped = skipped | 1 << depth;
      const nextDepth = depth + 1;
      for (let i = node.optional.length - 1; i >= 0; i--) {
        const segment = node.optional[i];
        stack.push({
          node: segment,
          index,
          skipped: nextSkipped,
          depth: nextDepth,
          statics,
          dynamics,
          optionals
        });
      }
      if (!isBeyondPath) {
        for (let i = node.optional.length - 1; i >= 0; i--) {
          const segment = node.optional[i];
          const { prefix: prefix2, suffix } = segment;
          if (prefix2 || suffix) {
            const casePart = segment.caseSensitive ? part : lowerPart ??= part.toLowerCase();
            if (prefix2 && !casePart.startsWith(prefix2)) continue;
            if (suffix && !casePart.endsWith(suffix)) continue;
          }
          stack.push({
            node: segment,
            index: index + 1,
            skipped,
            depth: nextDepth,
            statics,
            dynamics,
            optionals: optionals + 1
          });
        }
      }
    }
    if (!isBeyondPath && node.dynamic && part) {
      for (let i = node.dynamic.length - 1; i >= 0; i--) {
        const segment = node.dynamic[i];
        const { prefix: prefix2, suffix } = segment;
        if (prefix2 || suffix) {
          const casePart = segment.caseSensitive ? part : lowerPart ??= part.toLowerCase();
          if (prefix2 && !casePart.startsWith(prefix2)) continue;
          if (suffix && !casePart.endsWith(suffix)) continue;
        }
        stack.push({
          node: segment,
          index: index + 1,
          skipped,
          depth: depth + 1,
          statics,
          dynamics: dynamics + 1,
          optionals
        });
      }
    }
    if (!isBeyondPath && node.staticInsensitive) {
      const match = node.staticInsensitive.get(
        lowerPart ??= part.toLowerCase()
      );
      if (match) {
        stack.push({
          node: match,
          index: index + 1,
          skipped,
          depth: depth + 1,
          statics: statics + 1,
          dynamics,
          optionals
        });
      }
    }
    if (!isBeyondPath && node.static) {
      const match = node.static.get(part);
      if (match) {
        stack.push({
          node: match,
          index: index + 1,
          skipped,
          depth: depth + 1,
          statics: statics + 1,
          dynamics,
          optionals
        });
      }
    }
  }
  if (bestMatch && wildcardMatch) {
    return isFrameMoreSpecific(wildcardMatch, bestMatch) ? bestMatch : wildcardMatch;
  }
  if (bestMatch) return bestMatch;
  if (wildcardMatch) return wildcardMatch;
  if (fuzzy && bestFuzzy) {
    let sliceIndex = bestFuzzy.index;
    for (let i = 0; i < bestFuzzy.index; i++) {
      sliceIndex += parts[i].length;
    }
    const splat = sliceIndex === path.length ? "/" : path.slice(sliceIndex);
    return {
      node: bestFuzzy.node,
      skipped: bestFuzzy.skipped,
      "**": decodeURIComponent(splat)
    };
  }
  return null;
}
function isFrameMoreSpecific(prev, next) {
  if (!prev) return true;
  return next.statics > prev.statics || next.statics === prev.statics && (next.dynamics > prev.dynamics || next.dynamics === prev.dynamics && (next.optionals > prev.optionals || next.optionals === prev.optionals && (next.node.isIndex > prev.node.isIndex || next.node.isIndex === prev.node.isIndex && next.depth > prev.depth)));
}
function joinPaths(paths) {
  return cleanPath(
    paths.filter((val) => {
      return val !== void 0;
    }).join("/")
  );
}
function cleanPath(path) {
  return path.replace(/\/{2,}/g, "/");
}
function trimPathLeft(path) {
  return path === "/" ? path : path.replace(/^\/{1,}/, "");
}
function trimPathRight(path) {
  const len = path.length;
  return len > 1 && path[len - 1] === "/" ? path.replace(/\/{1,}$/, "") : path;
}
function trimPath(path) {
  return trimPathRight(trimPathLeft(path));
}
function removeTrailingSlash(value, basepath) {
  if (value?.endsWith("/") && value !== "/" && value !== `${basepath}/`) {
    return value.slice(0, -1);
  }
  return value;
}
function exactPathTest(pathName1, pathName2, basepath) {
  return removeTrailingSlash(pathName1, basepath) === removeTrailingSlash(pathName2, basepath);
}
function resolvePath({
  base,
  to: to2,
  trailingSlash = "never",
  cache
}) {
  const isAbsolute = to2.startsWith("/");
  const isBase = !isAbsolute && to2 === ".";
  let key;
  if (cache) {
    key = isAbsolute ? to2 : isBase ? base : base + "\0" + to2;
    const cached = cache.get(key);
    if (cached) return cached;
  }
  let baseSegments;
  if (isBase) {
    baseSegments = base.split("/");
  } else if (isAbsolute) {
    baseSegments = to2.split("/");
  } else {
    baseSegments = base.split("/");
    while (baseSegments.length > 1 && last(baseSegments) === "") {
      baseSegments.pop();
    }
    const toSegments = to2.split("/");
    for (let index = 0, length = toSegments.length; index < length; index++) {
      const value = toSegments[index];
      if (value === "") {
        if (!index) {
          baseSegments = [value];
        } else if (index === length - 1) {
          baseSegments.push(value);
        } else ;
      } else if (value === "..") {
        baseSegments.pop();
      } else if (value === ".") ;
      else {
        baseSegments.push(value);
      }
    }
  }
  if (baseSegments.length > 1) {
    if (last(baseSegments) === "") {
      if (trailingSlash === "never") {
        baseSegments.pop();
      }
    } else if (trailingSlash === "always") {
      baseSegments.push("");
    }
  }
  let segment;
  let joined = "";
  for (let i = 0; i < baseSegments.length; i++) {
    if (i > 0) joined += "/";
    const part = baseSegments[i];
    if (!part) continue;
    segment = parseSegment(part, 0, segment);
    const kind = segment[0];
    if (kind === SEGMENT_TYPE_PATHNAME) {
      joined += part;
      continue;
    }
    const end = segment[5];
    const prefix2 = part.substring(0, segment[1]);
    const suffix = part.substring(segment[4], end);
    const value = part.substring(segment[2], segment[3]);
    if (kind === SEGMENT_TYPE_PARAM) {
      joined += prefix2 || suffix ? `${prefix2}{$${value}}${suffix}` : `$${value}`;
    } else if (kind === SEGMENT_TYPE_WILDCARD) {
      joined += prefix2 || suffix ? `${prefix2}{$}${suffix}` : "$";
    } else {
      joined += `${prefix2}{-$${value}}${suffix}`;
    }
  }
  joined = cleanPath(joined);
  const result = joined || "/";
  if (key && cache) cache.set(key, result);
  return result;
}
function encodeParam(key, params, decodeCharMap) {
  const value = params[key];
  if (typeof value !== "string") return value;
  if (key === "_splat") {
    return encodeURI(value);
  } else {
    return encodePathParam(value, decodeCharMap);
  }
}
function interpolatePath({
  path,
  params,
  decodeCharMap
}) {
  let isMissingParams = false;
  const usedParams = {};
  if (!path || path === "/")
    return { interpolatedPath: "/", usedParams, isMissingParams };
  if (!path.includes("$"))
    return { interpolatedPath: path, usedParams, isMissingParams };
  const length = path.length;
  let cursor = 0;
  let segment;
  let joined = "";
  while (cursor < length) {
    const start = cursor;
    segment = parseSegment(path, start, segment);
    const end = segment[5];
    cursor = end + 1;
    if (start === end) continue;
    const kind = segment[0];
    if (kind === SEGMENT_TYPE_PATHNAME) {
      joined += "/" + path.substring(start, end);
      continue;
    }
    if (kind === SEGMENT_TYPE_WILDCARD) {
      const splat = params._splat;
      usedParams._splat = splat;
      usedParams["*"] = splat;
      const prefix2 = path.substring(start, segment[1]);
      const suffix = path.substring(segment[4], end);
      if (!splat) {
        isMissingParams = true;
        if (prefix2 || suffix) {
          joined += "/" + prefix2 + suffix;
        }
        continue;
      }
      const value = encodeParam("_splat", params, decodeCharMap);
      joined += "/" + prefix2 + value + suffix;
      continue;
    }
    if (kind === SEGMENT_TYPE_PARAM) {
      const key = path.substring(segment[2], segment[3]);
      if (!isMissingParams && !(key in params)) {
        isMissingParams = true;
      }
      usedParams[key] = params[key];
      const prefix2 = path.substring(start, segment[1]);
      const suffix = path.substring(segment[4], end);
      const value = encodeParam(key, params, decodeCharMap) ?? "undefined";
      joined += "/" + prefix2 + value + suffix;
      continue;
    }
    if (kind === SEGMENT_TYPE_OPTIONAL_PARAM) {
      const key = path.substring(segment[2], segment[3]);
      const valueRaw = params[key];
      if (valueRaw == null) continue;
      usedParams[key] = valueRaw;
      const prefix2 = path.substring(start, segment[1]);
      const suffix = path.substring(segment[4], end);
      const value = encodeParam(key, params, decodeCharMap) ?? "";
      joined += "/" + prefix2 + value + suffix;
      continue;
    }
  }
  if (path.endsWith("/")) joined += "/";
  const interpolatedPath = joined || "/";
  return { usedParams, interpolatedPath, isMissingParams };
}
function encodePathParam(value, decodeCharMap) {
  let encoded = encodeURIComponent(value);
  if (decodeCharMap) {
    for (const [encodedChar, char] of decodeCharMap) {
      encoded = encoded.replaceAll(encodedChar, char);
    }
  }
  return encoded;
}
function isNotFound(obj) {
  return !!obj?.isNotFound;
}
function getSafeSessionStorage() {
  try {
    if (typeof window !== "undefined" && typeof window.sessionStorage === "object") {
      return window.sessionStorage;
    }
  } catch {
  }
  return void 0;
}
const storageKey = "tsr-scroll-restoration-v1_3";
const throttle = (fn2, wait) => {
  let timeout;
  return (...args) => {
    if (!timeout) {
      timeout = setTimeout(() => {
        fn2(...args);
        timeout = null;
      }, wait);
    }
  };
};
function createScrollRestorationCache() {
  const safeSessionStorage = getSafeSessionStorage();
  if (!safeSessionStorage) {
    return null;
  }
  const persistedState = safeSessionStorage.getItem(storageKey);
  let state = persistedState ? JSON.parse(persistedState) : {};
  return {
    state,
    // This setter is simply to make sure that we set the sessionStorage right
    // after the state is updated. It doesn't necessarily need to be a functional
    // update.
    set: (updater) => (state = functionalUpdate(updater, state) || state, safeSessionStorage.setItem(storageKey, JSON.stringify(state)))
  };
}
const scrollRestorationCache = createScrollRestorationCache();
const defaultGetScrollRestorationKey = (location) => {
  return location.state.__TSR_key || location.href;
};
function getCssSelector(el) {
  const path = [];
  let parent;
  while (parent = el.parentNode) {
    path.push(
      `${el.tagName}:nth-child(${Array.prototype.indexOf.call(parent.children, el) + 1})`
    );
    el = parent;
  }
  return `${path.reverse().join(" > ")}`.toLowerCase();
}
let ignoreScroll = false;
function restoreScroll({
  storageKey: storageKey2,
  key,
  behavior,
  shouldScrollRestoration,
  scrollToTopSelectors,
  location
}) {
  let byKey;
  try {
    byKey = JSON.parse(sessionStorage.getItem(storageKey2) || "{}");
  } catch (error) {
    console.error(error);
    return;
  }
  const resolvedKey = key || window.history.state?.__TSR_key;
  const elementEntries = byKey[resolvedKey];
  ignoreScroll = true;
  scroll: {
    if (shouldScrollRestoration && elementEntries && Object.keys(elementEntries).length > 0) {
      for (const elementSelector in elementEntries) {
        const entry = elementEntries[elementSelector];
        if (elementSelector === "window") {
          window.scrollTo({
            top: entry.scrollY,
            left: entry.scrollX,
            behavior
          });
        } else if (elementSelector) {
          const element = document.querySelector(elementSelector);
          if (element) {
            element.scrollLeft = entry.scrollX;
            element.scrollTop = entry.scrollY;
          }
        }
      }
      break scroll;
    }
    const hash = (location ?? window.location).hash.split("#", 2)[1];
    if (hash) {
      const hashScrollIntoViewOptions = window.history.state?.__hashScrollIntoViewOptions ?? true;
      if (hashScrollIntoViewOptions) {
        const el = document.getElementById(hash);
        if (el) {
          el.scrollIntoView(hashScrollIntoViewOptions);
        }
      }
      break scroll;
    }
    const scrollOptions = { top: 0, left: 0, behavior };
    window.scrollTo(scrollOptions);
    if (scrollToTopSelectors) {
      for (const selector of scrollToTopSelectors) {
        if (selector === "window") continue;
        const element = typeof selector === "function" ? selector() : document.querySelector(selector);
        if (element) element.scrollTo(scrollOptions);
      }
    }
  }
  ignoreScroll = false;
}
function setupScrollRestoration(router, force) {
  if (!scrollRestorationCache && !router.isServer) {
    return;
  }
  const shouldScrollRestoration = router.options.scrollRestoration ?? false;
  if (shouldScrollRestoration) {
    router.isScrollRestoring = true;
  }
  if (router.isServer || router.isScrollRestorationSetup || !scrollRestorationCache) {
    return;
  }
  router.isScrollRestorationSetup = true;
  ignoreScroll = false;
  const getKey = router.options.getScrollRestorationKey || defaultGetScrollRestorationKey;
  window.history.scrollRestoration = "manual";
  const onScroll = (event) => {
    if (ignoreScroll || !router.isScrollRestoring) {
      return;
    }
    let elementSelector = "";
    if (event.target === document || event.target === window) {
      elementSelector = "window";
    } else {
      const attrId = event.target.getAttribute(
        "data-scroll-restoration-id"
      );
      if (attrId) {
        elementSelector = `[data-scroll-restoration-id="${attrId}"]`;
      } else {
        elementSelector = getCssSelector(event.target);
      }
    }
    const restoreKey = getKey(router.state.location);
    scrollRestorationCache.set((state) => {
      const keyEntry = state[restoreKey] ||= {};
      const elementEntry = keyEntry[elementSelector] ||= {};
      if (elementSelector === "window") {
        elementEntry.scrollX = window.scrollX || 0;
        elementEntry.scrollY = window.scrollY || 0;
      } else if (elementSelector) {
        const element = document.querySelector(elementSelector);
        if (element) {
          elementEntry.scrollX = element.scrollLeft || 0;
          elementEntry.scrollY = element.scrollTop || 0;
        }
      }
      return state;
    });
  };
  if (typeof document !== "undefined") {
    document.addEventListener("scroll", throttle(onScroll, 100), true);
  }
  router.subscribe("onRendered", (event) => {
    const cacheKey = getKey(event.toLocation);
    if (!router.resetNextScroll) {
      router.resetNextScroll = true;
      return;
    }
    if (typeof router.options.scrollRestoration === "function") {
      const shouldRestore = router.options.scrollRestoration({
        location: router.latestLocation
      });
      if (!shouldRestore) {
        return;
      }
    }
    restoreScroll({
      storageKey,
      key: cacheKey,
      behavior: router.options.scrollRestorationBehavior,
      shouldScrollRestoration: router.isScrollRestoring,
      scrollToTopSelectors: router.options.scrollToTopSelectors,
      location: router.history.location
    });
    if (router.isScrollRestoring) {
      scrollRestorationCache.set((state) => {
        state[cacheKey] ||= {};
        return state;
      });
    }
  });
}
function handleHashScroll(router) {
  if (typeof document !== "undefined" && document.querySelector) {
    const hashScrollIntoViewOptions = router.state.location.state.__hashScrollIntoViewOptions ?? true;
    if (hashScrollIntoViewOptions && router.state.location.hash !== "") {
      const el = document.getElementById(router.state.location.hash);
      if (el) {
        el.scrollIntoView(hashScrollIntoViewOptions);
      }
    }
  }
}
function encode(obj, stringify = String) {
  const result = new URLSearchParams();
  for (const key in obj) {
    const val = obj[key];
    if (val !== void 0) {
      result.set(key, stringify(val));
    }
  }
  return result.toString();
}
function toValue(str) {
  if (!str) return "";
  if (str === "false") return false;
  if (str === "true") return true;
  return +str * 0 === 0 && +str + "" === str ? +str : str;
}
function decode(str) {
  const searchParams = new URLSearchParams(str);
  const result = {};
  for (const [key, value] of searchParams.entries()) {
    const previousValue = result[key];
    if (previousValue == null) {
      result[key] = toValue(value);
    } else if (Array.isArray(previousValue)) {
      previousValue.push(toValue(value));
    } else {
      result[key] = [previousValue, toValue(value)];
    }
  }
  return result;
}
const defaultParseSearch = parseSearchWith(JSON.parse);
const defaultStringifySearch = stringifySearchWith(
  JSON.stringify,
  JSON.parse
);
function parseSearchWith(parser) {
  return (searchStr) => {
    if (searchStr[0] === "?") {
      searchStr = searchStr.substring(1);
    }
    const query = decode(searchStr);
    for (const key in query) {
      const value = query[key];
      if (typeof value === "string") {
        try {
          query[key] = parser(value);
        } catch (_err) {
        }
      }
    }
    return query;
  };
}
function stringifySearchWith(stringify, parser) {
  const hasParser = typeof parser === "function";
  function stringifyValue(val) {
    if (typeof val === "object" && val !== null) {
      try {
        return stringify(val);
      } catch (_err) {
      }
    } else if (hasParser && typeof val === "string") {
      try {
        parser(val);
        return stringify(val);
      } catch (_err) {
      }
    }
    return val;
  }
  return (search) => {
    const searchStr = encode(search, stringifyValue);
    return searchStr ? `?${searchStr}` : "";
  };
}
const rootRouteId = "__root__";
function redirect(opts) {
  opts.statusCode = opts.statusCode || opts.code || 307;
  if (!opts.reloadDocument && typeof opts.href === "string") {
    try {
      new URL(opts.href);
      opts.reloadDocument = true;
    } catch {
    }
  }
  const headers = new Headers(opts.headers);
  if (opts.href && headers.get("Location") === null) {
    headers.set("Location", opts.href);
  }
  const response = new Response(null, {
    status: opts.statusCode,
    headers
  });
  response.options = opts;
  if (opts.throw) {
    throw response;
  }
  return response;
}
function isRedirect(obj) {
  return obj instanceof Response && !!obj.options;
}
function isResolvedRedirect(obj) {
  return isRedirect(obj) && !!obj.options.href;
}
const triggerOnReady = (inner) => {
  if (!inner.rendered) {
    inner.rendered = true;
    return inner.onReady?.();
  }
};
const resolvePreload = (inner, matchId) => {
  return !!(inner.preload && !inner.router.state.matches.some((d) => d.id === matchId));
};
const _handleNotFound = (inner, err) => {
  const routeCursor = inner.router.routesById[err.routeId ?? ""] ?? inner.router.routeTree;
  if (!routeCursor.options.notFoundComponent && inner.router.options?.defaultNotFoundComponent) {
    routeCursor.options.notFoundComponent = inner.router.options.defaultNotFoundComponent;
  }
  invariant(
    routeCursor.options.notFoundComponent,
    "No notFoundComponent found. Please set a notFoundComponent on your route or provide a defaultNotFoundComponent to the router."
  );
  const matchForRoute = inner.matches.find((m2) => m2.routeId === routeCursor.id);
  invariant(matchForRoute, "Could not find match for route: " + routeCursor.id);
  inner.updateMatch(matchForRoute.id, (prev) => ({
    ...prev,
    status: "notFound",
    error: err,
    isFetching: false
  }));
  if (err.routerCode === "BEFORE_LOAD" && routeCursor.parentRoute) {
    err.routeId = routeCursor.parentRoute.id;
    _handleNotFound(inner, err);
  }
};
const handleRedirectAndNotFound = (inner, match, err) => {
  if (!isRedirect(err) && !isNotFound(err)) return;
  if (isRedirect(err) && err.redirectHandled && !err.options.reloadDocument) {
    throw err;
  }
  if (match) {
    match._nonReactive.beforeLoadPromise?.resolve();
    match._nonReactive.loaderPromise?.resolve();
    match._nonReactive.beforeLoadPromise = void 0;
    match._nonReactive.loaderPromise = void 0;
    const status = isRedirect(err) ? "redirected" : "notFound";
    match._nonReactive.error = err;
    inner.updateMatch(match.id, (prev) => ({
      ...prev,
      status,
      isFetching: false,
      error: err
    }));
    if (isNotFound(err) && !err.routeId) {
      err.routeId = match.routeId;
    }
    match._nonReactive.loadPromise?.resolve();
  }
  if (isRedirect(err)) {
    inner.rendered = true;
    err.options._fromLocation = inner.location;
    err.redirectHandled = true;
    err = inner.router.resolveRedirect(err);
    throw err;
  } else {
    _handleNotFound(inner, err);
    throw err;
  }
};
const shouldSkipLoader = (inner, matchId) => {
  const match = inner.router.getMatch(matchId);
  if (!inner.router.isServer && match._nonReactive.dehydrated) {
    return true;
  }
  if (inner.router.isServer && match.ssr === false) {
    return true;
  }
  return false;
};
const handleSerialError = (inner, index, err, routerCode) => {
  const { id: matchId, routeId } = inner.matches[index];
  const route = inner.router.looseRoutesById[routeId];
  if (err instanceof Promise) {
    throw err;
  }
  err.routerCode = routerCode;
  inner.firstBadMatchIndex ??= index;
  handleRedirectAndNotFound(inner, inner.router.getMatch(matchId), err);
  try {
    route.options.onError?.(err);
  } catch (errorHandlerErr) {
    err = errorHandlerErr;
    handleRedirectAndNotFound(inner, inner.router.getMatch(matchId), err);
  }
  inner.updateMatch(matchId, (prev) => {
    prev._nonReactive.beforeLoadPromise?.resolve();
    prev._nonReactive.beforeLoadPromise = void 0;
    prev._nonReactive.loadPromise?.resolve();
    return {
      ...prev,
      error: err,
      status: "error",
      isFetching: false,
      updatedAt: Date.now(),
      abortController: new AbortController()
    };
  });
};
const isBeforeLoadSsr = (inner, matchId, index, route) => {
  const existingMatch = inner.router.getMatch(matchId);
  const parentMatchId = inner.matches[index - 1]?.id;
  const parentMatch = parentMatchId ? inner.router.getMatch(parentMatchId) : void 0;
  if (inner.router.isShell()) {
    existingMatch.ssr = route.id === rootRouteId;
    return;
  }
  if (parentMatch?.ssr === false) {
    existingMatch.ssr = false;
    return;
  }
  const parentOverride = (tempSsr2) => {
    if (tempSsr2 === true && parentMatch?.ssr === "data-only") {
      return "data-only";
    }
    return tempSsr2;
  };
  const defaultSsr = inner.router.options.defaultSsr ?? true;
  if (route.options.ssr === void 0) {
    existingMatch.ssr = parentOverride(defaultSsr);
    return;
  }
  if (typeof route.options.ssr !== "function") {
    existingMatch.ssr = parentOverride(route.options.ssr);
    return;
  }
  const { search, params } = existingMatch;
  const ssrFnContext = {
    search: makeMaybe(search, existingMatch.searchError),
    params: makeMaybe(params, existingMatch.paramsError),
    location: inner.location,
    matches: inner.matches.map((match) => ({
      index: match.index,
      pathname: match.pathname,
      fullPath: match.fullPath,
      staticData: match.staticData,
      id: match.id,
      routeId: match.routeId,
      search: makeMaybe(match.search, match.searchError),
      params: makeMaybe(match.params, match.paramsError),
      ssr: match.ssr
    }))
  };
  const tempSsr = route.options.ssr(ssrFnContext);
  if (isPromise(tempSsr)) {
    return tempSsr.then((ssr) => {
      existingMatch.ssr = parentOverride(ssr ?? defaultSsr);
    });
  }
  existingMatch.ssr = parentOverride(tempSsr ?? defaultSsr);
  return;
};
const setupPendingTimeout = (inner, matchId, route, match) => {
  if (match._nonReactive.pendingTimeout !== void 0) return;
  const pendingMs = route.options.pendingMs ?? inner.router.options.defaultPendingMs;
  const shouldPending = !!(inner.onReady && !inner.router.isServer && !resolvePreload(inner, matchId) && (route.options.loader || route.options.beforeLoad || routeNeedsPreload(route)) && typeof pendingMs === "number" && pendingMs !== Infinity && (route.options.pendingComponent ?? inner.router.options?.defaultPendingComponent));
  if (shouldPending) {
    const pendingTimeout = setTimeout(() => {
      triggerOnReady(inner);
    }, pendingMs);
    match._nonReactive.pendingTimeout = pendingTimeout;
  }
};
const preBeforeLoadSetup = (inner, matchId, route) => {
  const existingMatch = inner.router.getMatch(matchId);
  if (!existingMatch._nonReactive.beforeLoadPromise && !existingMatch._nonReactive.loaderPromise)
    return;
  setupPendingTimeout(inner, matchId, route, existingMatch);
  const then = () => {
    const match = inner.router.getMatch(matchId);
    if (match.preload && (match.status === "redirected" || match.status === "notFound")) {
      handleRedirectAndNotFound(inner, match, match.error);
    }
  };
  return existingMatch._nonReactive.beforeLoadPromise ? existingMatch._nonReactive.beforeLoadPromise.then(then) : then();
};
const executeBeforeLoad = (inner, matchId, index, route) => {
  const match = inner.router.getMatch(matchId);
  const prevLoadPromise = match._nonReactive.loadPromise;
  match._nonReactive.loadPromise = createControlledPromise(() => {
    prevLoadPromise?.resolve();
  });
  const { paramsError, searchError } = match;
  if (paramsError) {
    handleSerialError(inner, index, paramsError, "PARSE_PARAMS");
  }
  if (searchError) {
    handleSerialError(inner, index, searchError, "VALIDATE_SEARCH");
  }
  setupPendingTimeout(inner, matchId, route, match);
  const abortController = new AbortController();
  const parentMatchId = inner.matches[index - 1]?.id;
  const parentMatch = parentMatchId ? inner.router.getMatch(parentMatchId) : void 0;
  const parentMatchContext = parentMatch?.context ?? inner.router.options.context ?? void 0;
  let isPending = false;
  const pending = () => {
    if (isPending) return;
    isPending = true;
    inner.updateMatch(matchId, (prev) => ({
      ...prev,
      isFetching: "beforeLoad",
      fetchCount: prev.fetchCount + 1,
      abortController
      // Note: We intentionally don't update context here.
      // Context should only be updated after beforeLoad resolves to avoid
      // components seeing incomplete context during async beforeLoad execution.
    }));
  };
  const resolve = () => {
    match._nonReactive.beforeLoadPromise?.resolve();
    match._nonReactive.beforeLoadPromise = void 0;
    inner.updateMatch(matchId, (prev) => ({
      ...prev,
      isFetching: false
    }));
  };
  if (!route.options.beforeLoad) {
    batch(() => {
      pending();
      resolve();
    });
    return;
  }
  match._nonReactive.beforeLoadPromise = createControlledPromise();
  const { search, params, cause } = match;
  const preload = resolvePreload(inner, matchId);
  const beforeLoadFnContext = {
    search,
    abortController,
    params,
    preload,
    // Include parent's __beforeLoadContext so child routes can access it during their beforeLoad
    context: {
      ...parentMatchContext,
      ...parentMatch?.__beforeLoadContext,
      ...match.__routeContext
    },
    location: inner.location,
    navigate: (opts) => inner.router.navigate({
      ...opts,
      _fromLocation: inner.location
    }),
    buildLocation: inner.router.buildLocation,
    cause: preload ? "preload" : cause,
    matches: inner.matches,
    ...inner.router.options.additionalContext
  };
  const updateContext = (beforeLoadContext2) => {
    if (beforeLoadContext2 === void 0) {
      batch(() => {
        pending();
        resolve();
      });
      return;
    }
    if (isRedirect(beforeLoadContext2) || isNotFound(beforeLoadContext2)) {
      pending();
      handleSerialError(inner, index, beforeLoadContext2, "BEFORE_LOAD");
    }
    batch(() => {
      pending();
      inner.updateMatch(matchId, (prev) => ({
        ...prev,
        __beforeLoadContext: beforeLoadContext2
      }));
      resolve();
    });
  };
  let beforeLoadContext;
  try {
    beforeLoadContext = route.options.beforeLoad(beforeLoadFnContext);
    if (isPromise(beforeLoadContext)) {
      pending();
      return beforeLoadContext.catch((err) => {
        handleSerialError(inner, index, err, "BEFORE_LOAD");
      }).then(updateContext);
    }
  } catch (err) {
    pending();
    handleSerialError(inner, index, err, "BEFORE_LOAD");
  }
  updateContext(beforeLoadContext);
  return;
};
const handleBeforeLoad = (inner, index) => {
  const { id: matchId, routeId } = inner.matches[index];
  const route = inner.router.looseRoutesById[routeId];
  const serverSsr = () => {
    if (inner.router.isServer) {
      const maybePromise = isBeforeLoadSsr(inner, matchId, index, route);
      if (isPromise(maybePromise)) return maybePromise.then(queueExecution);
    }
    return queueExecution();
  };
  const execute = () => executeBeforeLoad(inner, matchId, index, route);
  const queueExecution = () => {
    if (shouldSkipLoader(inner, matchId)) return;
    const result = preBeforeLoadSetup(inner, matchId, route);
    return isPromise(result) ? result.then(execute) : execute();
  };
  return serverSsr();
};
const executeHead = (inner, matchId, route) => {
  const match = inner.router.getMatch(matchId);
  if (!match) {
    return;
  }
  if (!route.options.head && !route.options.scripts && !route.options.headers) {
    return;
  }
  const assetContext = {
    matches: inner.matches,
    match,
    params: match.params,
    loaderData: match.loaderData
  };
  return Promise.all([
    route.options.head?.(assetContext),
    route.options.scripts?.(assetContext),
    route.options.headers?.(assetContext)
  ]).then(([headFnContent, scripts, headers]) => {
    const meta = headFnContent?.meta;
    const links = headFnContent?.links;
    const headScripts = headFnContent?.scripts;
    const styles = headFnContent?.styles;
    return {
      meta,
      links,
      headScripts,
      headers,
      scripts,
      styles
    };
  });
};
const getLoaderContext = (inner, matchId, index, route) => {
  const parentMatchPromise = inner.matchPromises[index - 1];
  const { params, loaderDeps, abortController, cause } = inner.router.getMatch(matchId);
  let context = inner.router.options.context ?? {};
  for (let i = 0; i <= index; i++) {
    const innerMatch = inner.matches[i];
    if (!innerMatch) continue;
    const m2 = inner.router.getMatch(innerMatch.id);
    if (!m2) continue;
    context = {
      ...context,
      ...m2.__routeContext ?? {},
      ...m2.__beforeLoadContext ?? {}
    };
  }
  const preload = resolvePreload(inner, matchId);
  return {
    params,
    deps: loaderDeps,
    preload: !!preload,
    parentMatchPromise,
    abortController,
    context,
    location: inner.location,
    navigate: (opts) => inner.router.navigate({
      ...opts,
      _fromLocation: inner.location
    }),
    cause: preload ? "preload" : cause,
    route,
    ...inner.router.options.additionalContext
  };
};
const runLoader = async (inner, matchId, index, route) => {
  try {
    const match = inner.router.getMatch(matchId);
    try {
      if (!inner.router.isServer || match.ssr === true) {
        loadRouteChunk(route);
      }
      const loaderResult = route.options.loader?.(
        getLoaderContext(inner, matchId, index, route)
      );
      const loaderResultIsPromise = route.options.loader && isPromise(loaderResult);
      const willLoadSomething = !!(loaderResultIsPromise || route._lazyPromise || route._componentsPromise || route.options.head || route.options.scripts || route.options.headers || match._nonReactive.minPendingPromise);
      if (willLoadSomething) {
        inner.updateMatch(matchId, (prev) => ({
          ...prev,
          isFetching: "loader"
        }));
      }
      if (route.options.loader) {
        const loaderData = loaderResultIsPromise ? await loaderResult : loaderResult;
        handleRedirectAndNotFound(
          inner,
          inner.router.getMatch(matchId),
          loaderData
        );
        if (loaderData !== void 0) {
          inner.updateMatch(matchId, (prev) => ({
            ...prev,
            loaderData
          }));
        }
      }
      if (route._lazyPromise) await route._lazyPromise;
      const headResult = executeHead(inner, matchId, route);
      const head = headResult ? await headResult : void 0;
      const pendingPromise = match._nonReactive.minPendingPromise;
      if (pendingPromise) await pendingPromise;
      if (route._componentsPromise) await route._componentsPromise;
      inner.updateMatch(matchId, (prev) => ({
        ...prev,
        error: void 0,
        status: "success",
        isFetching: false,
        updatedAt: Date.now(),
        ...head
      }));
    } catch (e) {
      let error = e;
      const pendingPromise = match._nonReactive.minPendingPromise;
      if (pendingPromise) await pendingPromise;
      if (isNotFound(e)) {
        await route.options.notFoundComponent?.preload?.();
      }
      handleRedirectAndNotFound(inner, inner.router.getMatch(matchId), e);
      try {
        route.options.onError?.(e);
      } catch (onErrorError) {
        error = onErrorError;
        handleRedirectAndNotFound(
          inner,
          inner.router.getMatch(matchId),
          onErrorError
        );
      }
      const headResult = executeHead(inner, matchId, route);
      const head = headResult ? await headResult : void 0;
      inner.updateMatch(matchId, (prev) => ({
        ...prev,
        error,
        status: "error",
        isFetching: false,
        ...head
      }));
    }
  } catch (err) {
    const match = inner.router.getMatch(matchId);
    if (match) {
      const headResult = executeHead(inner, matchId, route);
      if (headResult) {
        const head = await headResult;
        inner.updateMatch(matchId, (prev) => ({
          ...prev,
          ...head
        }));
      }
      match._nonReactive.loaderPromise = void 0;
    }
    handleRedirectAndNotFound(inner, match, err);
  }
};
const loadRouteMatch = async (inner, index) => {
  const { id: matchId, routeId } = inner.matches[index];
  let loaderShouldRunAsync = false;
  let loaderIsRunningAsync = false;
  const route = inner.router.looseRoutesById[routeId];
  const commitContext = () => {
    const context = { ...inner.router.options.context };
    for (let i = 0; i <= index; i++) {
      const innerMatch = inner.matches[i];
      if (!innerMatch) continue;
      const m2 = inner.router.getMatch(innerMatch.id);
      if (!m2) continue;
      Object.assign(context, m2.__routeContext, m2.__beforeLoadContext);
    }
    inner.updateMatch(matchId, (prev) => ({
      ...prev,
      context
    }));
  };
  if (shouldSkipLoader(inner, matchId)) {
    if (inner.router.isServer) {
      const headResult = executeHead(inner, matchId, route);
      if (headResult) {
        const head = await headResult;
        inner.updateMatch(matchId, (prev) => ({
          ...prev,
          ...head
        }));
      }
      return inner.router.getMatch(matchId);
    }
  } else {
    const prevMatch = inner.router.getMatch(matchId);
    if (prevMatch._nonReactive.loaderPromise) {
      if (prevMatch.status === "success" && !inner.sync && !prevMatch.preload) {
        return prevMatch;
      }
      await prevMatch._nonReactive.loaderPromise;
      const match2 = inner.router.getMatch(matchId);
      const error = match2._nonReactive.error || match2.error;
      if (error) {
        handleRedirectAndNotFound(inner, match2, error);
      }
    } else {
      const age = Date.now() - prevMatch.updatedAt;
      const preload = resolvePreload(inner, matchId);
      const staleAge = preload ? route.options.preloadStaleTime ?? inner.router.options.defaultPreloadStaleTime ?? 3e4 : route.options.staleTime ?? inner.router.options.defaultStaleTime ?? 0;
      const shouldReloadOption = route.options.shouldReload;
      const shouldReload = typeof shouldReloadOption === "function" ? shouldReloadOption(getLoaderContext(inner, matchId, index, route)) : shouldReloadOption;
      const nextPreload = !!preload && !inner.router.state.matches.some((d) => d.id === matchId);
      const match2 = inner.router.getMatch(matchId);
      match2._nonReactive.loaderPromise = createControlledPromise();
      if (nextPreload !== match2.preload) {
        inner.updateMatch(matchId, (prev) => ({
          ...prev,
          preload: nextPreload
        }));
      }
      const { status, invalid } = match2;
      loaderShouldRunAsync = status === "success" && (invalid || (shouldReload ?? age > staleAge));
      if (preload && route.options.preload === false) ;
      else if (loaderShouldRunAsync && !inner.sync) {
        loaderIsRunningAsync = true;
        (async () => {
          try {
            await runLoader(inner, matchId, index, route);
            commitContext();
            const match3 = inner.router.getMatch(matchId);
            match3._nonReactive.loaderPromise?.resolve();
            match3._nonReactive.loadPromise?.resolve();
            match3._nonReactive.loaderPromise = void 0;
          } catch (err) {
            if (isRedirect(err)) {
              await inner.router.navigate(err.options);
            }
          }
        })();
      } else if (status !== "success" || loaderShouldRunAsync && inner.sync) {
        await runLoader(inner, matchId, index, route);
      } else {
        const headResult = executeHead(inner, matchId, route);
        if (headResult) {
          const head = await headResult;
          inner.updateMatch(matchId, (prev) => ({
            ...prev,
            ...head
          }));
        }
      }
    }
  }
  const match = inner.router.getMatch(matchId);
  if (!loaderIsRunningAsync) {
    match._nonReactive.loaderPromise?.resolve();
    match._nonReactive.loadPromise?.resolve();
  }
  clearTimeout(match._nonReactive.pendingTimeout);
  match._nonReactive.pendingTimeout = void 0;
  if (!loaderIsRunningAsync) match._nonReactive.loaderPromise = void 0;
  match._nonReactive.dehydrated = void 0;
  if (!loaderIsRunningAsync) {
    commitContext();
  }
  const nextIsFetching = loaderIsRunningAsync ? match.isFetching : false;
  if (nextIsFetching !== match.isFetching || match.invalid !== false) {
    inner.updateMatch(matchId, (prev) => ({
      ...prev,
      isFetching: nextIsFetching,
      invalid: false
    }));
    return inner.router.getMatch(matchId);
  } else {
    return match;
  }
};
async function loadMatches(arg) {
  const inner = Object.assign(arg, {
    matchPromises: []
  });
  if (!inner.router.isServer && inner.router.state.matches.some((d) => d._forcePending)) {
    triggerOnReady(inner);
  }
  try {
    for (let i = 0; i < inner.matches.length; i++) {
      const beforeLoad = handleBeforeLoad(inner, i);
      if (isPromise(beforeLoad)) await beforeLoad;
    }
    const max = inner.firstBadMatchIndex ?? inner.matches.length;
    for (let i = 0; i < max; i++) {
      inner.matchPromises.push(loadRouteMatch(inner, i));
    }
    await Promise.all(inner.matchPromises);
    const readyPromise = triggerOnReady(inner);
    if (isPromise(readyPromise)) await readyPromise;
  } catch (err) {
    if (isNotFound(err) && !inner.preload) {
      const readyPromise = triggerOnReady(inner);
      if (isPromise(readyPromise)) await readyPromise;
      throw err;
    }
    if (isRedirect(err)) {
      throw err;
    }
  }
  return inner.matches;
}
async function loadRouteChunk(route) {
  if (!route._lazyLoaded && route._lazyPromise === void 0) {
    if (route.lazyFn) {
      route._lazyPromise = route.lazyFn().then((lazyRoute) => {
        const { id: _id, ...options } = lazyRoute.options;
        Object.assign(route.options, options);
        route._lazyLoaded = true;
        route._lazyPromise = void 0;
      });
    } else {
      route._lazyLoaded = true;
    }
  }
  if (!route._componentsLoaded && route._componentsPromise === void 0) {
    const loadComponents = () => {
      const preloads = [];
      for (const type of componentTypes) {
        const preload = route.options[type]?.preload;
        if (preload) preloads.push(preload());
      }
      if (preloads.length)
        return Promise.all(preloads).then(() => {
          route._componentsLoaded = true;
          route._componentsPromise = void 0;
        });
      route._componentsLoaded = true;
      route._componentsPromise = void 0;
      return;
    };
    route._componentsPromise = route._lazyPromise ? route._lazyPromise.then(loadComponents) : loadComponents();
  }
  return route._componentsPromise;
}
function makeMaybe(value, error) {
  if (error) {
    return { status: "error", error };
  }
  return { status: "success", value };
}
function routeNeedsPreload(route) {
  for (const componentType of componentTypes) {
    if (route.options[componentType]?.preload) {
      return true;
    }
  }
  return false;
}
const componentTypes = [
  "component",
  "errorComponent",
  "pendingComponent",
  "notFoundComponent"
];
function composeRewrites(rewrites) {
  return {
    input: ({ url }) => {
      for (const rewrite of rewrites) {
        url = executeRewriteInput(rewrite, url);
      }
      return url;
    },
    output: ({ url }) => {
      for (let i = rewrites.length - 1; i >= 0; i--) {
        url = executeRewriteOutput(rewrites[i], url);
      }
      return url;
    }
  };
}
function rewriteBasepath(opts) {
  const trimmedBasepath = trimPath(opts.basepath);
  const normalizedBasepath = `/${trimmedBasepath}`;
  const normalizedBasepathWithSlash = `${normalizedBasepath}/`;
  const checkBasepath = opts.caseSensitive ? normalizedBasepath : normalizedBasepath.toLowerCase();
  const checkBasepathWithSlash = opts.caseSensitive ? normalizedBasepathWithSlash : normalizedBasepathWithSlash.toLowerCase();
  return {
    input: ({ url }) => {
      const pathname = opts.caseSensitive ? url.pathname : url.pathname.toLowerCase();
      if (pathname === checkBasepath) {
        url.pathname = "/";
      } else if (pathname.startsWith(checkBasepathWithSlash)) {
        url.pathname = url.pathname.slice(normalizedBasepath.length);
      }
      return url;
    },
    output: ({ url }) => {
      url.pathname = joinPaths(["/", trimmedBasepath, url.pathname]);
      return url;
    }
  };
}
function executeRewriteInput(rewrite, url) {
  const res = rewrite?.input?.({ url });
  if (res) {
    if (typeof res === "string") {
      return new URL(res);
    } else if (res instanceof URL) {
      return res;
    }
  }
  return url;
}
function executeRewriteOutput(rewrite, url) {
  const res = rewrite?.output?.({ url });
  if (res) {
    if (typeof res === "string") {
      return new URL(res);
    } else if (res instanceof URL) {
      return res;
    }
  }
  return url;
}
function getLocationChangeInfo(routerState) {
  const fromLocation = routerState.resolvedLocation;
  const toLocation = routerState.location;
  const pathChanged = fromLocation?.pathname !== toLocation.pathname;
  const hrefChanged = fromLocation?.href !== toLocation.href;
  const hashChanged = fromLocation?.hash !== toLocation.hash;
  return { fromLocation, toLocation, pathChanged, hrefChanged, hashChanged };
}
class RouterCore {
  /**
   * @deprecated Use the `createRouter` function instead
   */
  constructor(options) {
    this.tempLocationKey = `${Math.round(
      Math.random() * 1e7
    )}`;
    this.resetNextScroll = true;
    this.shouldViewTransition = void 0;
    this.isViewTransitionTypesSupported = void 0;
    this.subscribers = /* @__PURE__ */ new Set();
    this.isScrollRestoring = false;
    this.isScrollRestorationSetup = false;
    this.startTransition = (fn2) => fn2();
    this.update = (newOptions) => {
      if (newOptions.notFoundRoute) {
        console.warn(
          "The notFoundRoute API is deprecated and will be removed in the next major version. See https://tanstack.com/router/v1/docs/framework/react/guide/not-found-errors#migrating-from-notfoundroute for more info."
        );
      }
      const prevOptions = this.options;
      const prevBasepath = this.basepath ?? prevOptions?.basepath ?? "/";
      const basepathWasUnset = this.basepath === void 0;
      const prevRewriteOption = prevOptions?.rewrite;
      this.options = {
        ...prevOptions,
        ...newOptions
      };
      this.isServer = this.options.isServer ?? typeof document === "undefined";
      this.pathParamsDecodeCharMap = this.options.pathParamsAllowedCharacters ? new Map(
        this.options.pathParamsAllowedCharacters.map((char) => [
          encodeURIComponent(char),
          char
        ])
      ) : void 0;
      if (!this.history || this.options.history && this.options.history !== this.history) {
        if (!this.options.history) {
          if (!this.isServer) {
            this.history = createBrowserHistory();
          }
        } else {
          this.history = this.options.history;
        }
      }
      this.origin = this.options.origin;
      if (!this.origin) {
        if (!this.isServer && window?.origin && window.origin !== "null") {
          this.origin = window.origin;
        } else {
          this.origin = "http://localhost";
        }
      }
      if (this.history) {
        this.updateLatestLocation();
      }
      if (this.options.routeTree !== this.routeTree) {
        this.routeTree = this.options.routeTree;
        this.buildRouteTree();
      }
      if (!this.__store && this.latestLocation) {
        this.__store = new Store(getInitialRouterState(this.latestLocation), {
          onUpdate: () => {
            this.__store.state = {
              ...this.state,
              cachedMatches: this.state.cachedMatches.filter(
                (d) => !["redirected"].includes(d.status)
              )
            };
          }
        });
        setupScrollRestoration(this);
      }
      let needsLocationUpdate = false;
      const nextBasepath = this.options.basepath ?? "/";
      const nextRewriteOption = this.options.rewrite;
      const basepathChanged = basepathWasUnset || prevBasepath !== nextBasepath;
      const rewriteChanged = prevRewriteOption !== nextRewriteOption;
      if (basepathChanged || rewriteChanged) {
        this.basepath = nextBasepath;
        const rewrites = [];
        if (trimPath(nextBasepath) !== "") {
          rewrites.push(
            rewriteBasepath({
              basepath: nextBasepath
            })
          );
        }
        if (nextRewriteOption) {
          rewrites.push(nextRewriteOption);
        }
        this.rewrite = rewrites.length === 0 ? void 0 : rewrites.length === 1 ? rewrites[0] : composeRewrites(rewrites);
        if (this.history) {
          this.updateLatestLocation();
        }
        needsLocationUpdate = true;
      }
      if (needsLocationUpdate && this.__store) {
        this.__store.state = {
          ...this.state,
          location: this.latestLocation
        };
      }
      if (typeof window !== "undefined" && "CSS" in window && typeof window.CSS?.supports === "function") {
        this.isViewTransitionTypesSupported = window.CSS.supports(
          "selector(:active-view-transition-type(a)"
        );
      }
    };
    this.updateLatestLocation = () => {
      this.latestLocation = this.parseLocation(
        this.history.location,
        this.latestLocation
      );
    };
    this.buildRouteTree = () => {
      const { routesById, routesByPath, processedTree } = processRouteTree(
        this.routeTree,
        this.options.caseSensitive,
        (route, i) => {
          route.init({
            originalIndex: i
          });
        }
      );
      if (this.options.routeMasks) {
        processRouteMasks(this.options.routeMasks, processedTree);
      }
      this.routesById = routesById;
      this.routesByPath = routesByPath;
      this.processedTree = processedTree;
      const notFoundRoute = this.options.notFoundRoute;
      if (notFoundRoute) {
        notFoundRoute.init({
          originalIndex: 99999999999
        });
        this.routesById[notFoundRoute.id] = notFoundRoute;
      }
    };
    this.subscribe = (eventType, fn2) => {
      const listener = {
        eventType,
        fn: fn2
      };
      this.subscribers.add(listener);
      return () => {
        this.subscribers.delete(listener);
      };
    };
    this.emit = (routerEvent) => {
      this.subscribers.forEach((listener) => {
        if (listener.eventType === routerEvent.type) {
          listener.fn(routerEvent);
        }
      });
    };
    this.parseLocation = (locationToParse, previousLocation) => {
      const parse = ({
        href,
        state
      }) => {
        const fullUrl = new URL(href, this.origin);
        const url = executeRewriteInput(this.rewrite, fullUrl);
        const parsedSearch = this.options.parseSearch(url.search);
        const searchStr = this.options.stringifySearch(parsedSearch);
        url.search = searchStr;
        const fullPath = url.href.replace(url.origin, "");
        const { pathname, hash } = url;
        return {
          href: fullPath,
          publicHref: href,
          url: url.href,
          pathname: decodePath(pathname),
          searchStr,
          search: replaceEqualDeep(previousLocation?.search, parsedSearch),
          hash: hash.split("#").reverse()[0] ?? "",
          state: replaceEqualDeep(previousLocation?.state, state)
        };
      };
      const location = parse(locationToParse);
      const { __tempLocation, __tempKey } = location.state;
      if (__tempLocation && (!__tempKey || __tempKey === this.tempLocationKey)) {
        const parsedTempLocation = parse(__tempLocation);
        parsedTempLocation.state.key = location.state.key;
        parsedTempLocation.state.__TSR_key = location.state.__TSR_key;
        delete parsedTempLocation.state.__tempLocation;
        return {
          ...parsedTempLocation,
          maskedLocation: location
        };
      }
      return location;
    };
    this.resolvePathCache = createLRUCache(1e3);
    this.resolvePathWithBase = (from, path) => {
      const resolvedPath = resolvePath({
        base: from,
        to: cleanPath(path),
        trailingSlash: this.options.trailingSlash,
        cache: this.resolvePathCache
      });
      return resolvedPath;
    };
    this.matchRoutes = (pathnameOrNext, locationSearchOrOpts, opts) => {
      if (typeof pathnameOrNext === "string") {
        return this.matchRoutesInternal(
          {
            pathname: pathnameOrNext,
            search: locationSearchOrOpts
          },
          opts
        );
      }
      return this.matchRoutesInternal(pathnameOrNext, locationSearchOrOpts);
    };
    this.getMatchedRoutes = (pathname) => {
      return getMatchedRoutes({
        pathname,
        routesById: this.routesById,
        processedTree: this.processedTree
      });
    };
    this.cancelMatch = (id) => {
      const match = this.getMatch(id);
      if (!match) return;
      match.abortController.abort();
      clearTimeout(match._nonReactive.pendingTimeout);
      match._nonReactive.pendingTimeout = void 0;
    };
    this.cancelMatches = () => {
      const currentPendingMatches = this.state.matches.filter(
        (match) => match.status === "pending"
      );
      const currentLoadingMatches = this.state.matches.filter(
        (match) => match.isFetching === "loader"
      );
      const matchesToCancelArray = /* @__PURE__ */ new Set([
        ...this.state.pendingMatches ?? [],
        ...currentPendingMatches,
        ...currentLoadingMatches
      ]);
      matchesToCancelArray.forEach((match) => {
        this.cancelMatch(match.id);
      });
    };
    this.buildLocation = (opts) => {
      const build = (dest = {}) => {
        const currentLocation = dest._fromLocation || this.pendingBuiltLocation || this.latestLocation;
        const allCurrentLocationMatches = this.matchRoutes(currentLocation, {
          _buildLocation: true
        });
        const lastMatch = last(allCurrentLocationMatches);
        if (dest.from && process.env.NODE_ENV !== "production" && dest._isNavigate) {
          const allFromMatches = this.getMatchedRoutes(dest.from).matchedRoutes;
          const matchedFrom = findLast(allCurrentLocationMatches, (d) => {
            return comparePaths(d.fullPath, dest.from);
          });
          const matchedCurrent = findLast(allFromMatches, (d) => {
            return comparePaths(d.fullPath, lastMatch.fullPath);
          });
          if (!matchedFrom && !matchedCurrent) {
            console.warn(`Could not find match for from: ${dest.from}`);
          }
        }
        const defaultedFromPath = dest.unsafeRelative === "path" ? currentLocation.pathname : dest.from ?? lastMatch.fullPath;
        const fromPath = this.resolvePathWithBase(defaultedFromPath, ".");
        const fromSearch = lastMatch.search;
        const fromParams = { ...lastMatch.params };
        const nextTo = dest.to ? this.resolvePathWithBase(fromPath, `${dest.to}`) : this.resolvePathWithBase(fromPath, ".");
        const nextParams = dest.params === false || dest.params === null ? {} : (dest.params ?? true) === true ? fromParams : Object.assign(
          fromParams,
          functionalUpdate(dest.params, fromParams)
        );
        const interpolatedNextTo = interpolatePath({
          path: nextTo,
          params: nextParams
        }).interpolatedPath;
        const destRoutes = this.matchRoutes(interpolatedNextTo, void 0, {
          _buildLocation: true
        }).map((d) => this.looseRoutesById[d.routeId]);
        if (Object.keys(nextParams).length > 0) {
          for (const route of destRoutes) {
            const fn2 = route.options.params?.stringify ?? route.options.stringifyParams;
            if (fn2) {
              Object.assign(nextParams, fn2(nextParams));
            }
          }
        }
        const nextPathname = opts.leaveParams ? (
          // Use the original template path for interpolation
          // This preserves the original parameter syntax including optional parameters
          nextTo
        ) : decodePath(
          interpolatePath({
            path: nextTo,
            params: nextParams,
            decodeCharMap: this.pathParamsDecodeCharMap
          }).interpolatedPath
        );
        let nextSearch = fromSearch;
        if (opts._includeValidateSearch && this.options.search?.strict) {
          const validatedSearch = {};
          destRoutes.forEach((route) => {
            if (route.options.validateSearch) {
              try {
                Object.assign(
                  validatedSearch,
                  validateSearch(route.options.validateSearch, {
                    ...validatedSearch,
                    ...nextSearch
                  })
                );
              } catch {
              }
            }
          });
          nextSearch = validatedSearch;
        }
        nextSearch = applySearchMiddleware({
          search: nextSearch,
          dest,
          destRoutes,
          _includeValidateSearch: opts._includeValidateSearch
        });
        nextSearch = replaceEqualDeep(fromSearch, nextSearch);
        const searchStr = this.options.stringifySearch(nextSearch);
        const hash = dest.hash === true ? currentLocation.hash : dest.hash ? functionalUpdate(dest.hash, currentLocation.hash) : void 0;
        const hashStr = hash ? `#${hash}` : "";
        let nextState = dest.state === true ? currentLocation.state : dest.state ? functionalUpdate(dest.state, currentLocation.state) : {};
        nextState = replaceEqualDeep(currentLocation.state, nextState);
        const fullPath = `${nextPathname}${searchStr}${hashStr}`;
        const url = new URL(fullPath, this.origin);
        const rewrittenUrl = executeRewriteOutput(this.rewrite, url);
        return {
          publicHref: rewrittenUrl.pathname + rewrittenUrl.search + rewrittenUrl.hash,
          href: fullPath,
          url: rewrittenUrl.href,
          pathname: nextPathname,
          search: nextSearch,
          searchStr,
          state: nextState,
          hash: hash ?? "",
          unmaskOnReload: dest.unmaskOnReload
        };
      };
      const buildWithMatches = (dest = {}, maskedDest) => {
        const next = build(dest);
        let maskedNext = maskedDest ? build(maskedDest) : void 0;
        if (!maskedNext) {
          const params = {};
          if (this.options.routeMasks) {
            const match = findFlatMatch(
              next.pathname,
              this.processedTree
            );
            if (match) {
              Object.assign(params, match.params);
              const {
                from: _from,
                params: maskParams,
                ...maskProps
              } = match.route;
              const nextParams = maskParams === false || maskParams === null ? {} : (maskParams ?? true) === true ? params : Object.assign(params, functionalUpdate(maskParams, params));
              maskedDest = {
                from: opts.from,
                ...maskProps,
                params: nextParams
              };
              maskedNext = build(maskedDest);
            }
          }
        }
        if (maskedNext) {
          next.maskedLocation = maskedNext;
        }
        return next;
      };
      if (opts.mask) {
        return buildWithMatches(opts, {
          from: opts.from,
          ...opts.mask
        });
      }
      return buildWithMatches(opts);
    };
    this.commitLocation = ({
      viewTransition,
      ignoreBlocker,
      ...next
    }) => {
      const isSameState = () => {
        const ignoredProps = [
          "key",
          // TODO: Remove in v2 - use __TSR_key instead
          "__TSR_key",
          "__TSR_index",
          "__hashScrollIntoViewOptions"
        ];
        ignoredProps.forEach((prop) => {
          next.state[prop] = this.latestLocation.state[prop];
        });
        const isEqual = deepEqual(next.state, this.latestLocation.state);
        ignoredProps.forEach((prop) => {
          delete next.state[prop];
        });
        return isEqual;
      };
      const isSameUrl = trimPathRight(this.latestLocation.href) === trimPathRight(next.href);
      const previousCommitPromise = this.commitLocationPromise;
      this.commitLocationPromise = createControlledPromise(() => {
        previousCommitPromise?.resolve();
      });
      if (isSameUrl && isSameState()) {
        this.load();
      } else {
        let { maskedLocation, hashScrollIntoView, ...nextHistory } = next;
        if (maskedLocation) {
          nextHistory = {
            ...maskedLocation,
            state: {
              ...maskedLocation.state,
              __tempKey: void 0,
              __tempLocation: {
                ...nextHistory,
                search: nextHistory.searchStr,
                state: {
                  ...nextHistory.state,
                  __tempKey: void 0,
                  __tempLocation: void 0,
                  __TSR_key: void 0,
                  key: void 0
                  // TODO: Remove in v2 - use __TSR_key instead
                }
              }
            }
          };
          if (nextHistory.unmaskOnReload ?? this.options.unmaskOnReload ?? false) {
            nextHistory.state.__tempKey = this.tempLocationKey;
          }
        }
        nextHistory.state.__hashScrollIntoViewOptions = hashScrollIntoView ?? this.options.defaultHashScrollIntoView ?? true;
        this.shouldViewTransition = viewTransition;
        this.history[next.replace ? "replace" : "push"](
          nextHistory.publicHref,
          nextHistory.state,
          { ignoreBlocker }
        );
      }
      this.resetNextScroll = next.resetScroll ?? true;
      if (!this.history.subscribers.size) {
        this.load();
      }
      return this.commitLocationPromise;
    };
    this.buildAndCommitLocation = ({
      replace,
      resetScroll,
      hashScrollIntoView,
      viewTransition,
      ignoreBlocker,
      href,
      ...rest
    } = {}) => {
      if (href) {
        const currentIndex = this.history.location.state.__TSR_index;
        const parsed = parseHref(href, {
          __TSR_index: replace ? currentIndex : currentIndex + 1
        });
        rest.to = parsed.pathname;
        rest.search = this.options.parseSearch(parsed.search);
        rest.hash = parsed.hash.slice(1);
      }
      const location = this.buildLocation({
        ...rest,
        _includeValidateSearch: true
      });
      this.pendingBuiltLocation = location;
      const commitPromise = this.commitLocation({
        ...location,
        viewTransition,
        replace,
        resetScroll,
        hashScrollIntoView,
        ignoreBlocker
      });
      Promise.resolve().then(() => {
        if (this.pendingBuiltLocation === location) {
          this.pendingBuiltLocation = void 0;
        }
      });
      return commitPromise;
    };
    this.navigate = async ({ to: to2, reloadDocument, href, ...rest }) => {
      if (!reloadDocument && href) {
        try {
          new URL(`${href}`);
          reloadDocument = true;
        } catch {
        }
      }
      if (reloadDocument) {
        if (!href) {
          const location = this.buildLocation({ to: to2, ...rest });
          href = location.url;
        }
        if (!rest.ignoreBlocker) {
          const historyWithBlockers = this.history;
          const blockers = historyWithBlockers.getBlockers?.() ?? [];
          for (const blocker of blockers) {
            if (blocker?.blockerFn) {
              const shouldBlock = await blocker.blockerFn({
                currentLocation: this.latestLocation,
                nextLocation: this.latestLocation,
                // External URLs don't have a next location in our router
                action: "PUSH"
              });
              if (shouldBlock) {
                return Promise.resolve();
              }
            }
          }
        }
        if (rest.replace) {
          window.location.replace(href);
        } else {
          window.location.href = href;
        }
        return Promise.resolve();
      }
      return this.buildAndCommitLocation({
        ...rest,
        href,
        to: to2,
        _isNavigate: true
      });
    };
    this.beforeLoad = () => {
      this.cancelMatches();
      this.updateLatestLocation();
      if (this.isServer) {
        const nextLocation = this.buildLocation({
          to: this.latestLocation.pathname,
          search: true,
          params: true,
          hash: true,
          state: true,
          _includeValidateSearch: true
        });
        const normalizeUrl = (url) => {
          try {
            return encodeURI(decodeURI(url));
          } catch {
            return url;
          }
        };
        if (trimPath(normalizeUrl(this.latestLocation.href)) !== trimPath(normalizeUrl(nextLocation.href))) {
          let href = nextLocation.url;
          if (this.origin && href.startsWith(this.origin)) {
            href = href.replace(this.origin, "") || "/";
          }
          throw redirect({ href });
        }
      }
      const pendingMatches = this.matchRoutes(this.latestLocation);
      this.__store.setState((s) => ({
        ...s,
        status: "pending",
        statusCode: 200,
        isLoading: true,
        location: this.latestLocation,
        pendingMatches,
        // If a cached moved to pendingMatches, remove it from cachedMatches
        cachedMatches: s.cachedMatches.filter(
          (d) => !pendingMatches.some((e) => e.id === d.id)
        )
      }));
    };
    this.load = async (opts) => {
      let redirect2;
      let notFound;
      let loadPromise;
      loadPromise = new Promise((resolve) => {
        this.startTransition(async () => {
          try {
            this.beforeLoad();
            const next = this.latestLocation;
            const prevLocation = this.state.resolvedLocation;
            if (!this.state.redirect) {
              this.emit({
                type: "onBeforeNavigate",
                ...getLocationChangeInfo({
                  resolvedLocation: prevLocation,
                  location: next
                })
              });
            }
            this.emit({
              type: "onBeforeLoad",
              ...getLocationChangeInfo({
                resolvedLocation: prevLocation,
                location: next
              })
            });
            await loadMatches({
              router: this,
              sync: opts?.sync,
              matches: this.state.pendingMatches,
              location: next,
              updateMatch: this.updateMatch,
              // eslint-disable-next-line @typescript-eslint/require-await
              onReady: async () => {
                this.startTransition(() => {
                  this.startViewTransition(async () => {
                    let exitingMatches = [];
                    let enteringMatches = [];
                    let stayingMatches = [];
                    batch(() => {
                      this.__store.setState((s) => {
                        const previousMatches = s.matches;
                        const newMatches = s.pendingMatches || s.matches;
                        exitingMatches = previousMatches.filter(
                          (match) => !newMatches.some((d) => d.id === match.id)
                        );
                        enteringMatches = newMatches.filter(
                          (match) => !previousMatches.some((d) => d.id === match.id)
                        );
                        stayingMatches = newMatches.filter(
                          (match) => previousMatches.some((d) => d.id === match.id)
                        );
                        return {
                          ...s,
                          isLoading: false,
                          loadedAt: Date.now(),
                          matches: newMatches,
                          pendingMatches: void 0,
                          /**
                           * When committing new matches, cache any exiting matches that are still usable.
                           * Routes that resolved with `status: 'error'` or `status: 'notFound'` are
                           * deliberately excluded from `cachedMatches` so that subsequent invalidations
                           * or reloads re-run their loaders instead of reusing the failed/not-found data.
                           */
                          cachedMatches: [
                            ...s.cachedMatches,
                            ...exitingMatches.filter(
                              (d) => d.status !== "error" && d.status !== "notFound"
                            )
                          ]
                        };
                      });
                      this.clearExpiredCache();
                    });
                    [
                      [exitingMatches, "onLeave"],
                      [enteringMatches, "onEnter"],
                      [stayingMatches, "onStay"]
                    ].forEach(([matches, hook]) => {
                      matches.forEach((match) => {
                        this.looseRoutesById[match.routeId].options[hook]?.(
                          match
                        );
                      });
                    });
                  });
                });
              }
            });
          } catch (err) {
            if (isRedirect(err)) {
              redirect2 = err;
              if (!this.isServer) {
                this.navigate({
                  ...redirect2.options,
                  replace: true,
                  ignoreBlocker: true
                });
              }
            } else if (isNotFound(err)) {
              notFound = err;
            }
            this.__store.setState((s) => ({
              ...s,
              statusCode: redirect2 ? redirect2.status : notFound ? 404 : s.matches.some((d) => d.status === "error") ? 500 : 200,
              redirect: redirect2
            }));
          }
          if (this.latestLoadPromise === loadPromise) {
            this.commitLocationPromise?.resolve();
            this.latestLoadPromise = void 0;
            this.commitLocationPromise = void 0;
          }
          resolve();
        });
      });
      this.latestLoadPromise = loadPromise;
      await loadPromise;
      while (this.latestLoadPromise && loadPromise !== this.latestLoadPromise) {
        await this.latestLoadPromise;
      }
      let newStatusCode = void 0;
      if (this.hasNotFoundMatch()) {
        newStatusCode = 404;
      } else if (this.__store.state.matches.some((d) => d.status === "error")) {
        newStatusCode = 500;
      }
      if (newStatusCode !== void 0) {
        this.__store.setState((s) => ({
          ...s,
          statusCode: newStatusCode
        }));
      }
    };
    this.startViewTransition = (fn2) => {
      const shouldViewTransition = this.shouldViewTransition ?? this.options.defaultViewTransition;
      delete this.shouldViewTransition;
      if (shouldViewTransition && typeof document !== "undefined" && "startViewTransition" in document && typeof document.startViewTransition === "function") {
        let startViewTransitionParams;
        if (typeof shouldViewTransition === "object" && this.isViewTransitionTypesSupported) {
          const next = this.latestLocation;
          const prevLocation = this.state.resolvedLocation;
          const resolvedViewTransitionTypes = typeof shouldViewTransition.types === "function" ? shouldViewTransition.types(
            getLocationChangeInfo({
              resolvedLocation: prevLocation,
              location: next
            })
          ) : shouldViewTransition.types;
          if (resolvedViewTransitionTypes === false) {
            fn2();
            return;
          }
          startViewTransitionParams = {
            update: fn2,
            types: resolvedViewTransitionTypes
          };
        } else {
          startViewTransitionParams = fn2;
        }
        document.startViewTransition(startViewTransitionParams);
      } else {
        fn2();
      }
    };
    this.updateMatch = (id, updater) => {
      this.startTransition(() => {
        const matchesKey = this.state.pendingMatches?.some((d) => d.id === id) ? "pendingMatches" : this.state.matches.some((d) => d.id === id) ? "matches" : this.state.cachedMatches.some((d) => d.id === id) ? "cachedMatches" : "";
        if (matchesKey) {
          this.__store.setState((s) => ({
            ...s,
            [matchesKey]: s[matchesKey]?.map(
              (d) => d.id === id ? updater(d) : d
            )
          }));
        }
      });
    };
    this.getMatch = (matchId) => {
      const findFn = (d) => d.id === matchId;
      return this.state.cachedMatches.find(findFn) ?? this.state.pendingMatches?.find(findFn) ?? this.state.matches.find(findFn);
    };
    this.invalidate = (opts) => {
      const invalidate = (d) => {
        if (opts?.filter?.(d) ?? true) {
          return {
            ...d,
            invalid: true,
            ...opts?.forcePending || d.status === "error" || d.status === "notFound" ? { status: "pending", error: void 0 } : void 0
          };
        }
        return d;
      };
      this.__store.setState((s) => ({
        ...s,
        matches: s.matches.map(invalidate),
        cachedMatches: s.cachedMatches.map(invalidate),
        pendingMatches: s.pendingMatches?.map(invalidate)
      }));
      this.shouldViewTransition = false;
      return this.load({ sync: opts?.sync });
    };
    this.resolveRedirect = (redirect2) => {
      if (!redirect2.options.href) {
        const location = this.buildLocation(redirect2.options);
        let href = location.url;
        if (this.origin && href.startsWith(this.origin)) {
          href = href.replace(this.origin, "") || "/";
        }
        redirect2.options.href = location.href;
        redirect2.headers.set("Location", href);
      }
      if (!redirect2.headers.get("Location")) {
        redirect2.headers.set("Location", redirect2.options.href);
      }
      return redirect2;
    };
    this.clearCache = (opts) => {
      const filter = opts?.filter;
      if (filter !== void 0) {
        this.__store.setState((s) => {
          return {
            ...s,
            cachedMatches: s.cachedMatches.filter(
              (m2) => !filter(m2)
            )
          };
        });
      } else {
        this.__store.setState((s) => {
          return {
            ...s,
            cachedMatches: []
          };
        });
      }
    };
    this.clearExpiredCache = () => {
      const filter = (d) => {
        const route = this.looseRoutesById[d.routeId];
        if (!route.options.loader) {
          return true;
        }
        const gcTime = (d.preload ? route.options.preloadGcTime ?? this.options.defaultPreloadGcTime : route.options.gcTime ?? this.options.defaultGcTime) ?? 5 * 60 * 1e3;
        const isError = d.status === "error";
        if (isError) return true;
        const gcEligible = Date.now() - d.updatedAt >= gcTime;
        return gcEligible;
      };
      this.clearCache({ filter });
    };
    this.loadRouteChunk = loadRouteChunk;
    this.preloadRoute = async (opts) => {
      const next = this.buildLocation(opts);
      let matches = this.matchRoutes(next, {
        throwOnError: true,
        preload: true,
        dest: opts
      });
      const activeMatchIds = new Set(
        [...this.state.matches, ...this.state.pendingMatches ?? []].map(
          (d) => d.id
        )
      );
      const loadedMatchIds = /* @__PURE__ */ new Set([
        ...activeMatchIds,
        ...this.state.cachedMatches.map((d) => d.id)
      ]);
      batch(() => {
        matches.forEach((match) => {
          if (!loadedMatchIds.has(match.id)) {
            this.__store.setState((s) => ({
              ...s,
              cachedMatches: [...s.cachedMatches, match]
            }));
          }
        });
      });
      try {
        matches = await loadMatches({
          router: this,
          matches,
          location: next,
          preload: true,
          updateMatch: (id, updater) => {
            if (activeMatchIds.has(id)) {
              matches = matches.map((d) => d.id === id ? updater(d) : d);
            } else {
              this.updateMatch(id, updater);
            }
          }
        });
        return matches;
      } catch (err) {
        if (isRedirect(err)) {
          if (err.options.reloadDocument) {
            return void 0;
          }
          return await this.preloadRoute({
            ...err.options,
            _fromLocation: next
          });
        }
        if (!isNotFound(err)) {
          console.error(err);
        }
        return void 0;
      }
    };
    this.matchRoute = (location, opts) => {
      const matchLocation = {
        ...location,
        to: location.to ? this.resolvePathWithBase(
          location.from || "",
          location.to
        ) : void 0,
        params: location.params || {},
        leaveParams: true
      };
      const next = this.buildLocation(matchLocation);
      if (opts?.pending && this.state.status !== "pending") {
        return false;
      }
      const pending = opts?.pending === void 0 ? !this.state.isLoading : opts.pending;
      const baseLocation = pending ? this.latestLocation : this.state.resolvedLocation || this.state.location;
      const match = findSingleMatch(
        next.pathname,
        opts?.caseSensitive ?? false,
        opts?.fuzzy ?? false,
        baseLocation.pathname,
        this.processedTree
      );
      if (!match) {
        return false;
      }
      if (location.params) {
        if (!deepEqual(match.params, location.params, { partial: true })) {
          return false;
        }
      }
      if (opts?.includeSearch ?? true) {
        return deepEqual(baseLocation.search, next.search, { partial: true }) ? match.params : false;
      }
      return match.params;
    };
    this.hasNotFoundMatch = () => {
      return this.__store.state.matches.some(
        (d) => d.status === "notFound" || d.globalNotFound
      );
    };
    this.update({
      defaultPreloadDelay: 50,
      defaultPendingMs: 1e3,
      defaultPendingMinMs: 500,
      context: void 0,
      ...options,
      caseSensitive: options.caseSensitive ?? false,
      notFoundMode: options.notFoundMode ?? "fuzzy",
      stringifySearch: options.stringifySearch ?? defaultStringifySearch,
      parseSearch: options.parseSearch ?? defaultParseSearch
    });
    if (typeof document !== "undefined") {
      self.__TSR_ROUTER__ = this;
    }
  }
  isShell() {
    return !!this.options.isShell;
  }
  isPrerendering() {
    return !!this.options.isPrerendering;
  }
  get state() {
    return this.__store.state;
  }
  get looseRoutesById() {
    return this.routesById;
  }
  matchRoutesInternal(next, opts) {
    const matchedRoutesResult = this.getMatchedRoutes(next.pathname);
    const { foundRoute, routeParams } = matchedRoutesResult;
    let { matchedRoutes } = matchedRoutesResult;
    let isGlobalNotFound = false;
    if (
      // If we found a route, and it's not an index route and we have left over path
      foundRoute ? foundRoute.path !== "/" && routeParams["**"] : (
        // Or if we didn't find a route and we have left over path
        trimPathRight(next.pathname)
      )
    ) {
      if (this.options.notFoundRoute) {
        matchedRoutes = [...matchedRoutes, this.options.notFoundRoute];
      } else {
        isGlobalNotFound = true;
      }
    }
    const globalNotFoundRouteId = (() => {
      if (!isGlobalNotFound) {
        return void 0;
      }
      if (this.options.notFoundMode !== "root") {
        for (let i = matchedRoutes.length - 1; i >= 0; i--) {
          const route = matchedRoutes[i];
          if (route.children) {
            return route.id;
          }
        }
      }
      return rootRouteId;
    })();
    const matches = [];
    const getParentContext = (parentMatch) => {
      const parentMatchId = parentMatch?.id;
      const parentContext = !parentMatchId ? this.options.context ?? void 0 : parentMatch.context ?? this.options.context ?? void 0;
      return parentContext;
    };
    matchedRoutes.forEach((route, index) => {
      const parentMatch = matches[index - 1];
      const [preMatchSearch, strictMatchSearch, searchError] = (() => {
        const parentSearch = parentMatch?.search ?? next.search;
        const parentStrictSearch = parentMatch?._strictSearch ?? void 0;
        try {
          const strictSearch = validateSearch(route.options.validateSearch, { ...parentSearch }) ?? void 0;
          return [
            {
              ...parentSearch,
              ...strictSearch
            },
            { ...parentStrictSearch, ...strictSearch },
            void 0
          ];
        } catch (err) {
          let searchParamError = err;
          if (!(err instanceof SearchParamError)) {
            searchParamError = new SearchParamError(err.message, {
              cause: err
            });
          }
          if (opts?.throwOnError) {
            throw searchParamError;
          }
          return [parentSearch, {}, searchParamError];
        }
      })();
      const loaderDeps = route.options.loaderDeps?.({
        search: preMatchSearch
      }) ?? "";
      const loaderDepsHash = loaderDeps ? JSON.stringify(loaderDeps) : "";
      const { interpolatedPath, usedParams } = interpolatePath({
        path: route.fullPath,
        params: routeParams,
        decodeCharMap: this.pathParamsDecodeCharMap
      });
      const matchId = (
        // route.id for disambiguation
        route.id + // interpolatedPath for param changes
        interpolatedPath + // explicit deps
        loaderDepsHash
      );
      const existingMatch = this.getMatch(matchId);
      const previousMatch = this.state.matches.find(
        (d) => d.routeId === route.id
      );
      const strictParams = existingMatch?._strictParams ?? usedParams;
      let paramsError = void 0;
      if (!existingMatch) {
        const strictParseParams = route.options.params?.parse ?? route.options.parseParams;
        if (strictParseParams) {
          try {
            Object.assign(
              strictParams,
              strictParseParams(strictParams)
            );
          } catch (err) {
            if (isNotFound(err) || isRedirect(err)) {
              paramsError = err;
            } else {
              paramsError = new PathParamError(err.message, {
                cause: err
              });
            }
            if (opts?.throwOnError) {
              throw paramsError;
            }
          }
        }
      }
      Object.assign(routeParams, strictParams);
      const cause = previousMatch ? "stay" : "enter";
      let match;
      if (existingMatch) {
        match = {
          ...existingMatch,
          cause,
          params: previousMatch ? replaceEqualDeep(previousMatch.params, routeParams) : routeParams,
          _strictParams: strictParams,
          search: previousMatch ? replaceEqualDeep(previousMatch.search, preMatchSearch) : replaceEqualDeep(existingMatch.search, preMatchSearch),
          _strictSearch: strictMatchSearch
        };
      } else {
        const status = route.options.loader || route.options.beforeLoad || route.lazyFn || routeNeedsPreload(route) ? "pending" : "success";
        match = {
          id: matchId,
          ssr: this.isServer ? void 0 : route.options.ssr,
          index,
          routeId: route.id,
          params: previousMatch ? replaceEqualDeep(previousMatch.params, routeParams) : routeParams,
          _strictParams: strictParams,
          pathname: interpolatedPath,
          updatedAt: Date.now(),
          search: previousMatch ? replaceEqualDeep(previousMatch.search, preMatchSearch) : preMatchSearch,
          _strictSearch: strictMatchSearch,
          searchError: void 0,
          status,
          isFetching: false,
          error: void 0,
          paramsError,
          __routeContext: void 0,
          _nonReactive: {
            loadPromise: createControlledPromise()
          },
          __beforeLoadContext: void 0,
          context: {},
          abortController: new AbortController(),
          fetchCount: 0,
          cause,
          loaderDeps: previousMatch ? replaceEqualDeep(previousMatch.loaderDeps, loaderDeps) : loaderDeps,
          invalid: false,
          preload: false,
          links: void 0,
          scripts: void 0,
          headScripts: void 0,
          meta: void 0,
          staticData: route.options.staticData || {},
          fullPath: route.fullPath
        };
      }
      if (!opts?.preload) {
        match.globalNotFound = globalNotFoundRouteId === route.id;
      }
      match.searchError = searchError;
      const parentContext = getParentContext(parentMatch);
      match.context = {
        ...parentContext,
        ...match.__routeContext,
        ...match.__beforeLoadContext
      };
      matches.push(match);
    });
    matches.forEach((match, index) => {
      const route = this.looseRoutesById[match.routeId];
      const existingMatch = this.getMatch(match.id);
      if (!existingMatch && opts?._buildLocation !== true) {
        const parentMatch = matches[index - 1];
        const parentContext = getParentContext(parentMatch);
        if (route.options.context) {
          const contextFnContext = {
            deps: match.loaderDeps,
            params: match.params,
            context: parentContext ?? {},
            location: next,
            navigate: (opts2) => this.navigate({ ...opts2, _fromLocation: next }),
            buildLocation: this.buildLocation,
            cause: match.cause,
            abortController: match.abortController,
            preload: !!match.preload,
            matches
          };
          match.__routeContext = route.options.context(contextFnContext) ?? void 0;
        }
        match.context = {
          ...parentContext,
          ...match.__routeContext,
          ...match.__beforeLoadContext
        };
      }
    });
    return matches;
  }
}
class SearchParamError extends Error {
}
class PathParamError extends Error {
}
const normalize = (str) => str.endsWith("/") && str.length > 1 ? str.slice(0, -1) : str;
function comparePaths(a, b2) {
  return normalize(a) === normalize(b2);
}
function getInitialRouterState(location) {
  return {
    loadedAt: 0,
    isLoading: false,
    isTransitioning: false,
    status: "idle",
    resolvedLocation: void 0,
    location,
    matches: [],
    pendingMatches: [],
    cachedMatches: [],
    statusCode: 200
  };
}
function validateSearch(validateSearch2, input) {
  if (validateSearch2 == null) return {};
  if ("~standard" in validateSearch2) {
    const result = validateSearch2["~standard"].validate(input);
    if (result instanceof Promise)
      throw new SearchParamError("Async validation not supported");
    if (result.issues)
      throw new SearchParamError(JSON.stringify(result.issues, void 0, 2), {
        cause: result
      });
    return result.value;
  }
  if ("parse" in validateSearch2) {
    return validateSearch2.parse(input);
  }
  if (typeof validateSearch2 === "function") {
    return validateSearch2(input);
  }
  return {};
}
function getMatchedRoutes({
  pathname,
  routesById,
  processedTree
}) {
  const routeParams = {};
  const trimmedPath = trimPathRight(pathname);
  let foundRoute = void 0;
  const match = findRouteMatch(trimmedPath, processedTree, true);
  if (match) {
    foundRoute = match.route;
    Object.assign(routeParams, match.params);
  }
  const matchedRoutes = match?.branch || [routesById[rootRouteId]];
  return { matchedRoutes, routeParams, foundRoute };
}
function applySearchMiddleware({
  search,
  dest,
  destRoutes,
  _includeValidateSearch
}) {
  const allMiddlewares = destRoutes.reduce(
    (acc, route) => {
      const middlewares = [];
      if ("search" in route.options) {
        if (route.options.search?.middlewares) {
          middlewares.push(...route.options.search.middlewares);
        }
      } else if (route.options.preSearchFilters || route.options.postSearchFilters) {
        const legacyMiddleware = ({
          search: search2,
          next
        }) => {
          let nextSearch = search2;
          if ("preSearchFilters" in route.options && route.options.preSearchFilters) {
            nextSearch = route.options.preSearchFilters.reduce(
              (prev, next2) => next2(prev),
              search2
            );
          }
          const result = next(nextSearch);
          if ("postSearchFilters" in route.options && route.options.postSearchFilters) {
            return route.options.postSearchFilters.reduce(
              (prev, next2) => next2(prev),
              result
            );
          }
          return result;
        };
        middlewares.push(legacyMiddleware);
      }
      if (_includeValidateSearch && route.options.validateSearch) {
        const validate = ({ search: search2, next }) => {
          const result = next(search2);
          try {
            const validatedSearch = {
              ...result,
              ...validateSearch(route.options.validateSearch, result) ?? void 0
            };
            return validatedSearch;
          } catch {
            return result;
          }
        };
        middlewares.push(validate);
      }
      return acc.concat(middlewares);
    },
    []
  ) ?? [];
  const final = ({ search: search2 }) => {
    if (!dest.search) {
      return {};
    }
    if (dest.search === true) {
      return search2;
    }
    return functionalUpdate(dest.search, search2);
  };
  allMiddlewares.push(final);
  const applyNext = (index, currentSearch) => {
    if (index >= allMiddlewares.length) {
      return currentSearch;
    }
    const middleware = allMiddlewares[index];
    const next = (newSearch) => {
      return applyNext(index + 1, newSearch);
    };
    return middleware({ search: currentSearch, next });
  };
  return applyNext(0, search);
}
var K = ((s) => (s[s.AggregateError = 1] = "AggregateError", s[s.ArrowFunction = 2] = "ArrowFunction", s[s.ErrorPrototypeStack = 4] = "ErrorPrototypeStack", s[s.ObjectAssign = 8] = "ObjectAssign", s[s.BigIntTypedArray = 16] = "BigIntTypedArray", s))(K || {});
var b = Symbol.asyncIterator, lr = Symbol.hasInstance, P$1 = Symbol.isConcatSpreadable, A = Symbol.iterator, ur = Symbol.match, cr = Symbol.matchAll, fr = Symbol.replace, Sr = Symbol.search, pr = Symbol.species, dr = Symbol.split, mr = Symbol.toPrimitive, x = Symbol.toStringTag, gr = Symbol.unscopables;
var Zr = { 0: "Symbol.asyncIterator", 1: "Symbol.hasInstance", 2: "Symbol.isConcatSpreadable", 3: "Symbol.iterator", 4: "Symbol.match", 5: "Symbol.matchAll", 6: "Symbol.replace", 7: "Symbol.search", 8: "Symbol.species", 9: "Symbol.split", 10: "Symbol.toPrimitive", 11: "Symbol.toStringTag", 12: "Symbol.unscopables" }, ge = { [b]: 0, [lr]: 1, [P$1]: 2, [A]: 3, [ur]: 4, [cr]: 5, [fr]: 6, [Sr]: 7, [pr]: 8, [dr]: 9, [mr]: 10, [x]: 11, [gr]: 12 }, Jr = { 0: b, 1: lr, 2: P$1, 3: A, 4: ur, 5: cr, 6: fr, 7: Sr, 8: pr, 9: dr, 10: mr, 11: x, 12: gr }, Hr = { 2: "!0", 3: "!1", 1: "void 0", 0: "null", 4: "-0", 5: "1/0", 6: "-1/0", 7: "0/0" }, $r = { 2: true, 3: false, 1: void 0, 0: null, 4: -0, 5: Number.POSITIVE_INFINITY, 6: Number.NEGATIVE_INFINITY, 7: Number.NaN };
var ye = { 0: "Error", 1: "EvalError", 2: "RangeError", 3: "ReferenceError", 4: "SyntaxError", 5: "TypeError", 6: "URIError" }, qr = { 0: Error, 1: EvalError, 2: RangeError, 3: ReferenceError, 4: SyntaxError, 5: TypeError, 6: URIError }, o = void 0;
function c(e, r, t, n, a, s, i, l, u2, g, S, d) {
  return { t: e, i: r, s: t, l: n, c: a, m: s, p: i, e: l, a: u2, f: g, b: S, o: d };
}
function _(e) {
  return c(2, o, e, o, o, o, o, o, o, o, o, o);
}
var Z = _(2), J = _(3), Ne = _(1), Ce = _(0), Xr = _(4), Qr = _(5), et = _(6), rt = _(7);
function Qt(e) {
  switch (e) {
    case '"':
      return '\\"';
    case "\\":
      return "\\\\";
    case `
`:
      return "\\n";
    case "\r":
      return "\\r";
    case "\b":
      return "\\b";
    case "	":
      return "\\t";
    case "\f":
      return "\\f";
    case "<":
      return "\\x3C";
    case "\u2028":
      return "\\u2028";
    case "\u2029":
      return "\\u2029";
    default:
      return;
  }
}
function y(e) {
  let r = "", t = 0, n;
  for (let a = 0, s = e.length; a < s; a++) n = Qt(e[a]), n && (r += e.slice(t, a) + n, t = a + 1);
  return t === 0 ? r = e : r += e.slice(t), r;
}
function en(e) {
  switch (e) {
    case "\\\\":
      return "\\";
    case '\\"':
      return '"';
    case "\\n":
      return `
`;
    case "\\r":
      return "\r";
    case "\\b":
      return "\b";
    case "\\t":
      return "	";
    case "\\f":
      return "\f";
    case "\\x3C":
      return "<";
    case "\\u2028":
      return "\u2028";
    case "\\u2029":
      return "\u2029";
    default:
      return e;
  }
}
function z(e) {
  return e.replace(/(\\\\|\\"|\\n|\\r|\\b|\\t|\\f|\\u2028|\\u2029|\\x3C)/g, en);
}
var B = "__SEROVAL_REFS__", ie = "$R", ve = `self.${ie}`;
function rn(e) {
  return e == null ? `${ve}=${ve}||[]` : `(${ve}=${ve}||{})["${y(e)}"]=[]`;
}
var yr = /* @__PURE__ */ new Map(), V = /* @__PURE__ */ new Map();
function Nr(e) {
  return yr.has(e);
}
function nn(e) {
  return V.has(e);
}
function tt(e) {
  if (Nr(e)) return yr.get(e);
  throw new be(e);
}
function nt(e) {
  if (nn(e)) return V.get(e);
  throw new Ae(e);
}
typeof globalThis != "undefined" ? Object.defineProperty(globalThis, B, { value: V, configurable: true, writable: false, enumerable: false }) : typeof window != "undefined" ? Object.defineProperty(window, B, { value: V, configurable: true, writable: false, enumerable: false }) : typeof self != "undefined" ? Object.defineProperty(self, B, { value: V, configurable: true, writable: false, enumerable: false }) : typeof global != "undefined" && Object.defineProperty(global, B, { value: V, configurable: true, writable: false, enumerable: false });
function Re(e) {
  return e instanceof EvalError ? 1 : e instanceof RangeError ? 2 : e instanceof ReferenceError ? 3 : e instanceof SyntaxError ? 4 : e instanceof TypeError ? 5 : e instanceof URIError ? 6 : 0;
}
function on(e) {
  let r = ye[Re(e)];
  return e.name !== r ? { name: e.name } : e.constructor.name !== r ? { name: e.constructor.name } : {};
}
function H(e, r) {
  let t = on(e), n = Object.getOwnPropertyNames(e);
  for (let a = 0, s = n.length, i; a < s; a++) i = n[a], i !== "name" && i !== "message" && (i === "stack" ? r & 4 && (t = t || {}, t[i] = e[i]) : (t = t || {}, t[i] = e[i]));
  return t;
}
function Ie(e) {
  return Object.isFrozen(e) ? 3 : Object.isSealed(e) ? 2 : Object.isExtensible(e) ? 0 : 1;
}
function Ee(e) {
  switch (e) {
    case Number.POSITIVE_INFINITY:
      return Qr;
    case Number.NEGATIVE_INFINITY:
      return et;
  }
  return e !== e ? rt : Object.is(e, -0) ? Xr : c(0, o, e, o, o, o, o, o, o, o, o, o);
}
function $(e) {
  return c(1, o, y(e), o, o, o, o, o, o, o, o, o);
}
function Pe(e) {
  return c(3, o, "" + e, o, o, o, o, o, o, o, o, o);
}
function at(e) {
  return c(4, e, o, o, o, o, o, o, o, o, o, o);
}
function xe(e, r) {
  let t = r.valueOf();
  return c(5, e, t !== t ? "" : r.toISOString(), o, o, o, o, o, o, o, o, o);
}
function Te(e, r) {
  return c(6, e, o, o, y(r.source), r.flags, o, o, o, o, o, o);
}
function st(e, r) {
  return c(17, e, ge[r], o, o, o, o, o, o, o, o, o);
}
function it(e, r) {
  return c(18, e, y(tt(r)), o, o, o, o, o, o, o, o, o);
}
function le(e, r, t) {
  return c(25, e, t, o, y(r), o, o, o, o, o, o, o);
}
function Oe(e, r, t) {
  return c(9, e, o, r.length, o, o, o, o, t, o, o, Ie(r));
}
function he(e, r) {
  return c(21, e, o, o, o, o, o, o, o, r, o, o);
}
function we(e, r, t) {
  return c(15, e, o, r.length, r.constructor.name, o, o, o, o, t, r.byteOffset, o);
}
function ze(e, r, t) {
  return c(16, e, o, r.length, r.constructor.name, o, o, o, o, t, r.byteOffset, o);
}
function ke(e, r, t) {
  return c(20, e, o, r.byteLength, o, o, o, o, o, t, r.byteOffset, o);
}
function _e(e, r, t) {
  return c(13, e, Re(r), o, o, y(r.message), t, o, o, o, o, o);
}
function De(e, r, t) {
  return c(14, e, Re(r), o, o, y(r.message), t, o, o, o, o, o);
}
function Fe(e, r, t) {
  return c(7, e, o, r, o, o, o, o, t, o, o, o);
}
function Be(e, r) {
  return c(28, o, o, o, o, o, o, o, [e, r], o, o, o);
}
function Ve(e, r) {
  return c(30, o, o, o, o, o, o, o, [e, r], o, o, o);
}
function Me(e, r, t) {
  return c(31, e, o, o, o, o, o, o, t, r, o, o);
}
function je(e, r) {
  return c(32, e, o, o, o, o, o, o, o, r, o, o);
}
function Ue(e, r) {
  return c(33, e, o, o, o, o, o, o, o, r, o, o);
}
function Le(e, r) {
  return c(34, e, o, o, o, o, o, o, o, r, o, o);
}
var an = { parsing: 1, serialization: 2, deserialization: 3 };
function sn(e) {
  return `Seroval Error (step: ${an[e]})`;
}
var ln = (e, r) => sn(e), ue = class extends Error {
  constructor(t, n) {
    super(ln(t));
    this.cause = n;
  }
}, h = class extends ue {
  constructor(r) {
    super("parsing", r);
  }
}, Ye = class extends ue {
  constructor(r) {
    super("deserialization", r);
  }
};
function M(e) {
  return `Seroval Error (specific: ${e})`;
}
var E = class extends Error {
  constructor(t) {
    super(M(1));
    this.value = t;
  }
}, k = class extends Error {
  constructor(r) {
    super(M(2));
  }
}, q = class extends Error {
  constructor(r) {
    super(M(3));
  }
}, D = class extends Error {
  constructor(r) {
    super(M(4));
  }
}, be = class extends Error {
  constructor(t) {
    super(M(5));
    this.value = t;
  }
}, Ae = class extends Error {
  constructor(r) {
    super(M(6));
  }
}, We = class extends Error {
  constructor(r) {
    super(M(7));
  }
};
var j = class {
  constructor(r, t) {
    this.value = r;
    this.replacement = t;
  }
};
var X = () => {
  let e = { p: 0, s: 0, f: 0 };
  return e.p = new Promise((r, t) => {
    e.s = r, e.f = t;
  }), e;
}, un = (e, r) => {
  e.s(r), e.p.s = 1, e.p.v = r;
}, cn = (e, r) => {
  e.f(r), e.p.s = 2, e.p.v = r;
}, lt = X.toString(), ut = un.toString(), ct = cn.toString(), vr = () => {
  let e = [], r = [], t = true, n = false, a = 0, s = (u2, g, S) => {
    for (S = 0; S < a; S++) r[S] && r[S][g](u2);
  }, i = (u2, g, S, d) => {
    for (g = 0, S = e.length; g < S; g++) d = e[g], !t && g === S - 1 ? u2[n ? "return" : "throw"](d) : u2.next(d);
  }, l = (u2, g) => (t && (g = a++, r[g] = u2), i(u2), () => {
    t && (r[g] = r[a], r[a--] = void 0);
  });
  return { __SEROVAL_STREAM__: true, on: (u2) => l(u2), next: (u2) => {
    t && (e.push(u2), s(u2, "next"));
  }, throw: (u2) => {
    t && (e.push(u2), s(u2, "throw"), t = false, n = false, r.length = 0);
  }, return: (u2) => {
    t && (e.push(u2), s(u2, "return"), t = false, n = true, r.length = 0);
  } };
}, ft = vr.toString(), br = (e) => (r) => () => {
  let t = 0, n = { [e]: () => n, next: () => {
    if (t > r.d) return { done: true, value: void 0 };
    let a = t++, s = r.v[a];
    if (a === r.t) throw s;
    return { done: a === r.d, value: s };
  } };
  return n;
}, St = br.toString(), Ar = (e, r) => (t) => () => {
  let n = 0, a = -1, s = false, i = [], l = [], u2 = (S = 0, d = l.length) => {
    for (; S < d; S++) l[S].s({ done: true, value: void 0 });
  };
  t.on({ next: (S) => {
    let d = l.shift();
    d && d.s({ done: false, value: S }), i.push(S);
  }, throw: (S) => {
    let d = l.shift();
    d && d.f(S), u2(), a = i.length, s = true, i.push(S);
  }, return: (S) => {
    let d = l.shift();
    d && d.s({ done: true, value: S }), u2(), a = i.length, i.push(S);
  } });
  let g = { [e]: () => g, next: () => {
    if (a === -1) {
      let W = n++;
      if (W >= i.length) {
        let Gr = r();
        return l.push(Gr), Gr.p;
      }
      return { done: false, value: i[W] };
    }
    if (n > a) return { done: true, value: void 0 };
    let S = n++, d = i[S];
    if (S !== a) return { done: false, value: d };
    if (s) throw d;
    return { done: true, value: d };
  } };
  return g;
}, pt = Ar.toString(), Rr = (e, r) => {
  let t = atob(r), n = new Uint8Array(e);
  for (let a = 0; a < e; a++) n[a] = t.charCodeAt(a);
  return n.buffer;
}, dt = Rr.toString();
var mt = {}, gt = {};
var yt = { 0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {} }, Nt = { 0: "[]", 1: lt, 2: ut, 3: ct, 4: ft, 5: dt };
function Ke(e) {
  return "__SEROVAL_STREAM__" in e;
}
function Q() {
  return vr();
}
function Ge(e) {
  let r = Q(), t = e[b]();
  async function n() {
    try {
      let a = await t.next();
      a.done ? r.return(a.value) : (r.next(a.value), await n());
    } catch (a) {
      r.throw(a);
    }
  }
  return n().catch(() => {
  }), r;
}
var fn = Ar(b, X);
function Ct(e) {
  return fn(e);
}
function Ze(e) {
  let r = [], t = -1, n = -1, a = e[A]();
  for (; ; ) try {
    let s = a.next();
    if (r.push(s.value), s.done) {
      n = r.length - 1;
      break;
    }
  } catch (s) {
    t = r.length, r.push(s);
  }
  return { v: r, t, d: n };
}
var Sn = br(A);
function vt(e) {
  return Sn(e);
}
async function Ir(e) {
  try {
    return [1, await e];
  } catch (r) {
    return [0, r];
  }
}
function ce(e, r) {
  return { plugins: r.plugins, mode: e, marked: /* @__PURE__ */ new Set(), features: 31 ^ (r.disabledFeatures || 0), refs: r.refs || /* @__PURE__ */ new Map() };
}
function fe(e, r) {
  e.marked.add(r);
}
function Er(e, r) {
  let t = e.refs.size;
  return e.refs.set(r, t), t;
}
function Je(e, r) {
  let t = e.refs.get(r);
  return t != null ? (fe(e, t), { type: 1, value: at(t) }) : { type: 0, value: Er(e, r) };
}
function U(e, r) {
  let t = Je(e, r);
  return t.type === 1 ? t : Nr(r) ? { type: 2, value: it(t.value, r) } : t;
}
function I(e, r) {
  let t = U(e, r);
  if (t.type !== 0) return t.value;
  if (r in ge) return st(t.value, r);
  throw new E(r);
}
function w$1(e, r) {
  let t = Je(e, yt[r]);
  return t.type === 1 ? t.value : c(26, t.value, r, o, o, o, o, o, o, o, o, o);
}
function He(e) {
  let r = Je(e, mt);
  return r.type === 1 ? r.value : c(27, r.value, o, o, o, o, o, o, o, I(e, A), o, o);
}
function $e(e) {
  let r = Je(e, gt);
  return r.type === 1 ? r.value : c(29, r.value, o, o, o, o, o, o, [w$1(e, 1), I(e, b)], o, o, o);
}
function qe(e, r, t, n) {
  return c(t ? 11 : 10, e, o, o, o, o, n, o, o, o, o, Ie(r));
}
function Xe(e, r, t, n, a) {
  return c(8, r, o, o, o, o, o, { k: t, v: n, s: a }, o, w$1(e, 0), o, o);
}
function At(e, r, t) {
  return c(22, r, t, o, o, o, o, o, o, w$1(e, 1), o, o);
}
function Qe(e, r, t) {
  let n = new Uint8Array(t), a = n.length, s = "";
  for (let i = 0; i < a; i++) s += String.fromCharCode(n[i]);
  return c(19, r, y(btoa(s)), a, o, o, o, o, o, w$1(e, 5), o, o);
}
function ee$1(e, r) {
  return { base: ce(e, r), child: void 0 };
}
var xr = class {
  constructor(r) {
    this._p = r;
  }
  parse(r) {
    return N(this._p, r);
  }
};
async function mn(e, r) {
  let t = [];
  for (let n = 0, a = r.length; n < a; n++) n in r && (t[n] = await N(e, r[n]));
  return t;
}
async function gn(e, r, t) {
  return Oe(r, t, await mn(e, t));
}
async function Tr(e, r) {
  let t = Object.entries(r), n = [], a = [];
  for (let s = 0, i = t.length; s < i; s++) n.push(y(t[s][0])), a.push(await N(e, t[s][1]));
  return A in r && (n.push(I(e.base, A)), a.push(Be(He(e.base), await N(e, Ze(r))))), b in r && (n.push(I(e.base, b)), a.push(Ve($e(e.base), await N(e, Ge(r))))), x in r && (n.push(I(e.base, x)), a.push($(r[x]))), P$1 in r && (n.push(I(e.base, P$1)), a.push(r[P$1] ? Z : J)), { k: n, v: a, s: n.length };
}
async function Pr(e, r, t, n) {
  return qe(r, t, n, await Tr(e, t));
}
async function yn(e, r, t) {
  return he(r, await N(e, t.valueOf()));
}
async function Nn(e, r, t) {
  return we(r, t, await N(e, t.buffer));
}
async function Cn(e, r, t) {
  return ze(r, t, await N(e, t.buffer));
}
async function vn(e, r, t) {
  return ke(r, t, await N(e, t.buffer));
}
async function Rt(e, r, t) {
  let n = H(t, e.base.features);
  return _e(r, t, n ? await Tr(e, n) : o);
}
async function bn(e, r, t) {
  let n = H(t, e.base.features);
  return De(r, t, n ? await Tr(e, n) : o);
}
async function An(e, r, t) {
  let n = [], a = [];
  for (let [s, i] of t.entries()) n.push(await N(e, s)), a.push(await N(e, i));
  return Xe(e.base, r, n, a, t.size);
}
async function Rn(e, r, t) {
  let n = [];
  for (let a of t.keys()) n.push(await N(e, a));
  return Fe(r, t.size, n);
}
async function It(e, r, t) {
  let n = e.base.plugins;
  if (n) for (let a = 0, s = n.length; a < s; a++) {
    let i = n[a];
    if (i.parse.async && i.test(t)) return e.child == null && (e.child = new xr(e)), le(r, i.tag, await i.parse.async(t, e.child, { id: r }));
  }
  return o;
}
async function In(e, r, t) {
  let [n, a] = await Ir(t);
  return c(12, r, n, o, o, o, o, o, o, await N(e, a), o, o);
}
function En(e, r, t, n) {
  let a = [], s = r.on({ next: (i) => {
    fe(this.base, e), N(this, i).then((l) => {
      a.push(je(e, l));
    }, (l) => {
      n(l), s();
    });
  }, throw: (i) => {
    fe(this.base, e), N(this, i).then((l) => {
      a.push(Ue(e, l)), t(a), s();
    }, (l) => {
      n(l), s();
    });
  }, return: (i) => {
    fe(this.base, e), N(this, i).then((l) => {
      a.push(Le(e, l)), t(a), s();
    }, (l) => {
      n(l), s();
    });
  } });
}
async function Pn(e, r, t) {
  return Me(r, w$1(e.base, 4), await new Promise(En.bind(e, r, t)));
}
async function xn(e, r, t) {
  if (Array.isArray(t)) return gn(e, r, t);
  if (Ke(t)) return Pn(e, r, t);
  let n = t.constructor;
  if (n === j) return N(e, t.replacement);
  let a = await It(e, r, t);
  if (a) return a;
  switch (n) {
    case Object:
      return Pr(e, r, t, false);
    case o:
      return Pr(e, r, t, true);
    case Date:
      return xe(r, t);
    case RegExp:
      return Te(r, t);
    case Error:
    case EvalError:
    case RangeError:
    case ReferenceError:
    case SyntaxError:
    case TypeError:
    case URIError:
      return Rt(e, r, t);
    case Number:
    case Boolean:
    case String:
    case BigInt:
      return yn(e, r, t);
    case ArrayBuffer:
      return Qe(e.base, r, t);
    case Int8Array:
    case Int16Array:
    case Int32Array:
    case Uint8Array:
    case Uint16Array:
    case Uint32Array:
    case Uint8ClampedArray:
    case Float32Array:
    case Float64Array:
      return Nn(e, r, t);
    case DataView:
      return vn(e, r, t);
    case Map:
      return An(e, r, t);
    case Set:
      return Rn(e, r, t);
  }
  if (n === Promise || t instanceof Promise) return In(e, r, t);
  let s = e.base.features;
  if (s & 16) switch (n) {
    case BigInt64Array:
    case BigUint64Array:
      return Cn(e, r, t);
  }
  if (s & 1 && typeof AggregateError != "undefined" && (n === AggregateError || t instanceof AggregateError)) return bn(e, r, t);
  if (t instanceof Error) return Rt(e, r, t);
  if (A in t || b in t) return Pr(e, r, t, !!n);
  throw new E(t);
}
async function Tn(e, r) {
  let t = U(e.base, r);
  if (t.type !== 0) return t.value;
  let n = await It(e, t.value, r);
  if (n) return n;
  throw new E(r);
}
async function N(e, r) {
  switch (typeof r) {
    case "boolean":
      return r ? Z : J;
    case "undefined":
      return Ne;
    case "string":
      return $(r);
    case "number":
      return Ee(r);
    case "bigint":
      return Pe(r);
    case "object": {
      if (r) {
        let t = U(e.base, r);
        return t.type === 0 ? await xn(e, t.value, r) : t.value;
      }
      return Ce;
    }
    case "symbol":
      return I(e.base, r);
    case "function":
      return Tn(e, r);
    default:
      throw new E(r);
  }
}
async function re$1(e, r) {
  try {
    return await N(e, r);
  } catch (t) {
    throw t instanceof h ? t : new h(t);
  }
}
var te = ((t) => (t[t.Vanilla = 1] = "Vanilla", t[t.Cross = 2] = "Cross", t))(te || {});
function _s(e) {
  return e;
}
function Et(e, r) {
  for (let t = 0, n = r.length; t < n; t++) {
    let a = r[t];
    e.has(a) || (e.add(a), a.extends && Et(e, a.extends));
  }
}
function v(e) {
  if (e) {
    let r = /* @__PURE__ */ new Set();
    return Et(r, e), [...r];
  }
}
function Pt(e) {
  switch (e) {
    case "Int8Array":
      return Int8Array;
    case "Int16Array":
      return Int16Array;
    case "Int32Array":
      return Int32Array;
    case "Uint8Array":
      return Uint8Array;
    case "Uint16Array":
      return Uint16Array;
    case "Uint32Array":
      return Uint32Array;
    case "Uint8ClampedArray":
      return Uint8ClampedArray;
    case "Float32Array":
      return Float32Array;
    case "Float64Array":
      return Float64Array;
    case "BigInt64Array":
      return BigInt64Array;
    case "BigUint64Array":
      return BigUint64Array;
    default:
      throw new We(e);
  }
}
function xt(e, r) {
  switch (r) {
    case 3:
      return Object.freeze(e);
    case 1:
      return Object.preventExtensions(e);
    case 2:
      return Object.seal(e);
    default:
      return e;
  }
}
function Tt(e, r) {
  return { mode: e, plugins: r.plugins, refs: r.refs || /* @__PURE__ */ new Map() };
}
function Ot(e) {
  return { mode: 1, base: Tt(1, e), child: void 0, state: { marked: new Set(e.markedRefs) } };
}
var Or = class {
  constructor(r) {
    this._p = r;
  }
  deserialize(r) {
    return m(this._p, r);
  }
};
function On(e, r, t) {
  return e.state.marked.has(r) && e.base.refs.set(r, t), t;
}
function hn(e, r, t) {
  return e.base.refs.has(r) || e.base.refs.set(r, t), t;
}
function C(e, r, t) {
  return e.mode === 1 ? On(e, r, t) : hn(e, r, t);
}
function wn(e, r) {
  return C(e, r.i, nt(z(r.s)));
}
function zn(e, r) {
  let t = r.l, n = C(e, r.i, new Array(t)), a;
  for (let s = 0; s < t; s++) a = r.a[s], a && (n[s] = m(e, a));
  return xt(n, r.o), n;
}
function wt(e, r, t) {
  let n = r.s;
  if (n) {
    let a = r.k, s = r.v;
    for (let i = 0, l; i < n; i++) l = a[i], typeof l == "string" ? t[z(l)] = m(e, s[i]) : t[m(e, l)] = m(e, s[i]);
  }
  return t;
}
function kn(e, r) {
  let t = C(e, r.i, r.t === 10 ? {} : /* @__PURE__ */ Object.create(null));
  return wt(e, r.p, t), xt(t, r.o), t;
}
function _n(e, r) {
  return C(e, r.i, new Date(r.s));
}
function Dn(e, r) {
  return C(e, r.i, new RegExp(z(r.c), r.m));
}
function Fn(e, r) {
  let t = C(e, r.i, /* @__PURE__ */ new Set()), n = r.a;
  for (let a = 0, s = r.l; a < s; a++) t.add(m(e, n[a]));
  return t;
}
function Bn(e, r) {
  let t = C(e, r.i, /* @__PURE__ */ new Map()), n = r.e.k, a = r.e.v;
  for (let s = 0, i = r.e.s; s < i; s++) t.set(m(e, n[s]), m(e, a[s]));
  return t;
}
function Vn(e, r) {
  return C(e, r.i, Rr(r.l, z(r.s)));
}
function Mn(e, r) {
  let t = Pt(r.c), n = m(e, r.f);
  return C(e, r.i, new t(n, r.b, r.l));
}
function jn(e, r) {
  let t = m(e, r.f);
  return C(e, r.i, new DataView(t, r.b, r.l));
}
function zt(e, r, t) {
  if (r.p) {
    let n = wt(e, r.p, {});
    Object.assign(t, n);
  }
  return t;
}
function Un(e, r) {
  let t = C(e, r.i, new AggregateError([], z(r.m)));
  return zt(e, r, t);
}
function Ln(e, r) {
  let t = qr[r.s], n = C(e, r.i, new t(z(r.m)));
  return zt(e, r, n);
}
function Yn(e, r) {
  let t = X(), n = C(e, r.i, t.p), a = m(e, r.f);
  return r.s ? t.s(a) : t.f(a), n;
}
function Wn(e, r) {
  return C(e, r.i, Object(m(e, r.f)));
}
function Kn(e, r) {
  let t = e.base.plugins;
  if (t) {
    let n = z(r.c);
    for (let a = 0, s = t.length; a < s; a++) {
      let i = t[a];
      if (i.tag === n) return e.child == null && (e.child = new Or(e)), C(e, r.i, i.deserialize(r.s, e.child, { id: r.i }));
    }
  }
  throw new q(r.c);
}
function Gn(e, r) {
  return C(e, r.i, C(e, r.s, X()).p);
}
function Zn(e, r) {
  let t = e.base.refs.get(r.i);
  if (t) {
    t.s(m(e, r.a[1]));
    return;
  }
  throw new D("Promise");
}
function Jn(e, r) {
  let t = e.base.refs.get(r.i);
  if (t) {
    t.f(m(e, r.a[1]));
    return;
  }
  throw new D("Promise");
}
function Hn(e, r) {
  m(e, r.a[0]);
  let t = m(e, r.a[1]);
  return vt(t);
}
function $n(e, r) {
  m(e, r.a[0]);
  let t = m(e, r.a[1]);
  return Ct(t);
}
function qn(e, r) {
  let t = C(e, r.i, Q()), n = r.a.length;
  if (n) for (let a = 0; a < n; a++) m(e, r.a[a]);
  return t;
}
function Xn(e, r) {
  let t = e.base.refs.get(r.i);
  if (t) {
    t.next(m(e, r.f));
    return;
  }
  throw new D("Stream");
}
function Qn(e, r) {
  let t = e.base.refs.get(r.i);
  if (t) {
    t.throw(m(e, r.f));
    return;
  }
  throw new D("Stream");
}
function eo(e, r) {
  let t = e.base.refs.get(r.i);
  if (t) {
    t.return(m(e, r.f));
    return;
  }
  throw new D("Stream");
}
function ro(e, r) {
  m(e, r.f);
}
function to(e, r) {
  m(e, r.a[1]);
}
function m(e, r) {
  switch (r.t) {
    case 2:
      return $r[r.s];
    case 0:
      return r.s;
    case 1:
      return z(r.s);
    case 3:
      return BigInt(r.s);
    case 4:
      return e.base.refs.get(r.i);
    case 18:
      return wn(e, r);
    case 9:
      return zn(e, r);
    case 10:
    case 11:
      return kn(e, r);
    case 5:
      return _n(e, r);
    case 6:
      return Dn(e, r);
    case 7:
      return Fn(e, r);
    case 8:
      return Bn(e, r);
    case 19:
      return Vn(e, r);
    case 16:
    case 15:
      return Mn(e, r);
    case 20:
      return jn(e, r);
    case 14:
      return Un(e, r);
    case 13:
      return Ln(e, r);
    case 12:
      return Yn(e, r);
    case 17:
      return Jr[r.s];
    case 21:
      return Wn(e, r);
    case 25:
      return Kn(e, r);
    case 22:
      return Gn(e, r);
    case 23:
      return Zn(e, r);
    case 24:
      return Jn(e, r);
    case 28:
      return Hn(e, r);
    case 30:
      return $n(e, r);
    case 31:
      return qn(e, r);
    case 32:
      return Xn(e, r);
    case 33:
      return Qn(e, r);
    case 34:
      return eo(e, r);
    case 27:
      return ro(e, r);
    case 29:
      return to(e, r);
    default:
      throw new k(r);
  }
}
function er(e, r) {
  try {
    return m(e, r);
  } catch (t) {
    throw new Ye(t);
  }
}
var no = () => T, oo = no.toString(), kt = /=>/.test(oo);
function rr(e, r) {
  return kt ? (e.length === 1 ? e[0] : "(" + e.join(",") + ")") + "=>" + (r.startsWith("{") ? "(" + r + ")" : r) : "function(" + e.join(",") + "){return " + r + "}";
}
function _t(e, r) {
  return kt ? (e.length === 1 ? e[0] : "(" + e.join(",") + ")") + "=>{" + r + "}" : "function(" + e.join(",") + "){" + r + "}";
}
var Bt = "hjkmoquxzABCDEFGHIJKLNPQRTUVWXYZ$_", Dt = Bt.length, Vt = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789$_", Ft = Vt.length;
function hr(e) {
  let r = e % Dt, t = Bt[r];
  for (e = (e - r) / Dt; e > 0; ) r = e % Ft, t += Vt[r], e = (e - r) / Ft;
  return t;
}
var ao = /^[$A-Z_][0-9A-Z_$]*$/i;
function wr(e) {
  let r = e[0];
  return (r === "$" || r === "_" || r >= "A" && r <= "Z" || r >= "a" && r <= "z") && ao.test(e);
}
function pe(e) {
  switch (e.t) {
    case 0:
      return e.s + "=" + e.v;
    case 2:
      return e.s + ".set(" + e.k + "," + e.v + ")";
    case 1:
      return e.s + ".add(" + e.v + ")";
    case 3:
      return e.s + ".delete(" + e.k + ")";
  }
}
function so(e) {
  let r = [], t = e[0];
  for (let n = 1, a = e.length, s, i = t; n < a; n++) s = e[n], s.t === 0 && s.v === i.v ? t = { t: 0, s: s.s, k: o, v: pe(t) } : s.t === 2 && s.s === i.s ? t = { t: 2, s: pe(t), k: s.k, v: s.v } : s.t === 1 && s.s === i.s ? t = { t: 1, s: pe(t), k: o, v: s.v } : s.t === 3 && s.s === i.s ? t = { t: 3, s: pe(t), k: s.k, v: o } : (r.push(t), t = s), i = s;
  return r.push(t), r;
}
function Wt(e) {
  if (e.length) {
    let r = "", t = so(e);
    for (let n = 0, a = t.length; n < a; n++) r += pe(t[n]) + ",";
    return r;
  }
  return o;
}
var io = "Object.create(null)", lo = "new Set", uo = "new Map", co = "Promise.resolve", fo = "Promise.reject", So = { 3: "Object.freeze", 2: "Object.seal", 1: "Object.preventExtensions", 0: o };
function Kt(e, r) {
  return { mode: e, plugins: r.plugins, features: r.features, marked: new Set(r.markedRefs), stack: [], flags: [], assignments: [] };
}
function nr(e) {
  return { mode: 2, base: Kt(2, e), state: e, child: void 0 };
}
var zr = class {
  constructor(r) {
    this._p = r;
  }
  serialize(r) {
    return f(this._p, r);
  }
};
function mo(e, r) {
  let t = e.valid.get(r);
  t == null && (t = e.valid.size, e.valid.set(r, t));
  let n = e.vars[t];
  return n == null && (n = hr(t), e.vars[t] = n), n;
}
function go(e) {
  return ie + "[" + e + "]";
}
function p$1(e, r) {
  return e.mode === 1 ? mo(e.state, r) : go(r);
}
function O(e, r) {
  e.marked.add(r);
}
function kr(e, r) {
  return e.marked.has(r);
}
function Dr(e, r, t) {
  r !== 0 && (O(e.base, t), e.base.flags.push({ type: r, value: p$1(e, t) }));
}
function yo(e) {
  let r = "";
  for (let t = 0, n = e.flags, a = n.length; t < a; t++) {
    let s = n[t];
    r += So[s.type] + "(" + s.value + "),";
  }
  return r;
}
function Gt(e) {
  let r = Wt(e.assignments), t = yo(e);
  return r ? t ? r + t : r : t;
}
function Zt(e, r, t) {
  e.assignments.push({ t: 0, s: r, k: o, v: t });
}
function No(e, r, t) {
  e.base.assignments.push({ t: 1, s: p$1(e, r), k: o, v: t });
}
function Se(e, r, t, n) {
  e.base.assignments.push({ t: 2, s: p$1(e, r), k: t, v: n });
}
function Mt(e, r, t) {
  e.base.assignments.push({ t: 3, s: p$1(e, r), k: t, v: o });
}
function de(e, r, t, n) {
  Zt(e.base, p$1(e, r) + "[" + t + "]", n);
}
function _r(e, r, t, n) {
  Zt(e.base, p$1(e, r) + "." + t, n);
}
function F(e, r) {
  return r.t === 4 && e.stack.includes(r.i);
}
function ne(e, r, t) {
  return e.mode === 1 && !kr(e.base, r) ? t : p$1(e, r) + "=" + t;
}
function Co(e) {
  return B + '.get("' + e.s + '")';
}
function jt(e, r, t, n) {
  return t ? F(e.base, t) ? (O(e.base, r), de(e, r, n, p$1(e, t.i)), "") : f(e, t) : "";
}
function vo(e, r) {
  let t = r.i;
  if (r.l) {
    e.base.stack.push(t);
    let n = r.a, a = jt(e, t, n[0], 0), s = a === "";
    for (let i = 1, l = r.l, u2; i < l; i++) u2 = jt(e, t, n[i], i), a += "," + u2, s = u2 === "";
    return e.base.stack.pop(), Dr(e, r.o, r.i), "[" + a + (s ? ",]" : "]");
  }
  return "[]";
}
function Ut(e, r, t, n) {
  if (typeof t == "string") {
    let a = Number(t), s = a >= 0 && a.toString() === t || wr(t);
    if (F(e.base, n)) {
      let i = p$1(e, n.i);
      return O(e.base, r.i), s && a !== a ? _r(e, r.i, t, i) : de(e, r.i, s ? t : '"' + t + '"', i), "";
    }
    return (s ? t : '"' + t + '"') + ":" + f(e, n);
  }
  return "[" + f(e, t) + "]:" + f(e, n);
}
function Jt(e, r, t) {
  let n = t.s;
  if (n) {
    let a = t.k, s = t.v;
    e.base.stack.push(r.i);
    let i = Ut(e, r, a[0], s[0]);
    for (let l = 1, u2 = i; l < n; l++) u2 = Ut(e, r, a[l], s[l]), i += (u2 && i && ",") + u2;
    return e.base.stack.pop(), "{" + i + "}";
  }
  return "{}";
}
function bo(e, r) {
  return Dr(e, r.o, r.i), Jt(e, r, r.p);
}
function Ao(e, r, t, n) {
  let a = Jt(e, r, t);
  return a !== "{}" ? "Object.assign(" + n + "," + a + ")" : n;
}
function Ro(e, r, t, n, a) {
  let s = e.base, i = f(e, a), l = Number(n), u2 = l >= 0 && l.toString() === n || wr(n);
  if (F(s, a)) u2 && l !== l ? _r(e, r.i, n, i) : de(e, r.i, u2 ? n : '"' + n + '"', i);
  else {
    let g = s.assignments;
    s.assignments = t, u2 && l !== l ? _r(e, r.i, n, i) : de(e, r.i, u2 ? n : '"' + n + '"', i), s.assignments = g;
  }
}
function Io(e, r, t, n, a) {
  if (typeof n == "string") Ro(e, r, t, n, a);
  else {
    let s = e.base, i = s.stack;
    s.stack = [];
    let l = f(e, a);
    s.stack = i;
    let u2 = s.assignments;
    s.assignments = t, de(e, r.i, f(e, n), l), s.assignments = u2;
  }
}
function Eo(e, r, t) {
  let n = t.s;
  if (n) {
    let a = [], s = t.k, i = t.v;
    e.base.stack.push(r.i);
    for (let l = 0; l < n; l++) Io(e, r, a, s[l], i[l]);
    return e.base.stack.pop(), Wt(a);
  }
  return o;
}
function Fr(e, r, t) {
  if (r.p) {
    let n = e.base;
    if (n.features & 8) t = Ao(e, r, r.p, t);
    else {
      O(n, r.i);
      let a = Eo(e, r, r.p);
      if (a) return "(" + ne(e, r.i, t) + "," + a + p$1(e, r.i) + ")";
    }
  }
  return t;
}
function Po(e, r) {
  return Dr(e, r.o, r.i), Fr(e, r, io);
}
function xo(e) {
  return 'new Date("' + e.s + '")';
}
function To(e) {
  return "/" + e.c + "/" + e.m;
}
function Lt(e, r, t) {
  let n = e.base;
  return F(n, t) ? (O(n, r), No(e, r, p$1(e, t.i)), "") : f(e, t);
}
function Oo(e, r) {
  let t = lo, n = r.l, a = r.i;
  if (n) {
    let s = r.a;
    e.base.stack.push(a);
    let i = Lt(e, a, s[0]);
    for (let l = 1, u2 = i; l < n; l++) u2 = Lt(e, a, s[l]), i += (u2 && i && ",") + u2;
    e.base.stack.pop(), i && (t += "([" + i + "])");
  }
  return t;
}
function Yt(e, r, t, n, a) {
  let s = e.base;
  if (F(s, t)) {
    let i = p$1(e, t.i);
    if (O(s, r), F(s, n)) {
      let u2 = p$1(e, n.i);
      return Se(e, r, i, u2), "";
    }
    if (n.t !== 4 && n.i != null && kr(s, n.i)) {
      let u2 = "(" + f(e, n) + ",[" + a + "," + a + "])";
      return Se(e, r, i, p$1(e, n.i)), Mt(e, r, a), u2;
    }
    let l = s.stack;
    return s.stack = [], Se(e, r, i, f(e, n)), s.stack = l, "";
  }
  if (F(s, n)) {
    let i = p$1(e, n.i);
    if (O(s, r), t.t !== 4 && t.i != null && kr(s, t.i)) {
      let u2 = "(" + f(e, t) + ",[" + a + "," + a + "])";
      return Se(e, r, p$1(e, t.i), i), Mt(e, r, a), u2;
    }
    let l = s.stack;
    return s.stack = [], Se(e, r, f(e, t), i), s.stack = l, "";
  }
  return "[" + f(e, t) + "," + f(e, n) + "]";
}
function ho(e, r) {
  let t = uo, n = r.e.s, a = r.i, s = r.f, i = p$1(e, s.i), l = e.base;
  if (n) {
    let u2 = r.e.k, g = r.e.v;
    l.stack.push(a);
    let S = Yt(e, a, u2[0], g[0], i);
    for (let d = 1, W = S; d < n; d++) W = Yt(e, a, u2[d], g[d], i), S += (W && S && ",") + W;
    l.stack.pop(), S && (t += "([" + S + "])");
  }
  return s.t === 26 && (O(l, s.i), t = "(" + f(e, s) + "," + t + ")"), t;
}
function wo(e, r) {
  return L(e, r.f) + "(" + r.l + ',"' + r.s + '")';
}
function zo(e, r) {
  return "new " + r.c + "(" + f(e, r.f) + "," + r.b + "," + r.l + ")";
}
function ko(e, r) {
  return "new DataView(" + f(e, r.f) + "," + r.b + "," + r.l + ")";
}
function _o(e, r) {
  let t = r.i;
  e.base.stack.push(t);
  let n = Fr(e, r, 'new AggregateError([],"' + r.m + '")');
  return e.base.stack.pop(), n;
}
function Do(e, r) {
  return Fr(e, r, "new " + ye[r.s] + '("' + r.m + '")');
}
function Fo(e, r) {
  let t, n = r.f, a = r.i, s = r.s ? co : fo, i = e.base;
  if (F(i, n)) {
    let l = p$1(e, n.i);
    t = s + (r.s ? "().then(" + rr([], l) + ")" : "().catch(" + _t([], "throw " + l) + ")");
  } else {
    i.stack.push(a);
    let l = f(e, n);
    i.stack.pop(), t = s + "(" + l + ")";
  }
  return t;
}
function Bo(e, r) {
  return "Object(" + f(e, r.f) + ")";
}
function L(e, r) {
  let t = f(e, r);
  return r.t === 4 ? t : "(" + t + ")";
}
function Vo(e, r) {
  if (e.mode === 1) throw new k(r);
  return "(" + ne(e, r.s, L(e, r.f) + "()") + ").p";
}
function Mo(e, r) {
  if (e.mode === 1) throw new k(r);
  return L(e, r.a[0]) + "(" + p$1(e, r.i) + "," + f(e, r.a[1]) + ")";
}
function jo(e, r) {
  if (e.mode === 1) throw new k(r);
  return L(e, r.a[0]) + "(" + p$1(e, r.i) + "," + f(e, r.a[1]) + ")";
}
function Uo(e, r) {
  let t = e.base.plugins;
  if (t) for (let n = 0, a = t.length; n < a; n++) {
    let s = t[n];
    if (s.tag === r.c) return e.child == null && (e.child = new zr(e)), s.serialize(r.s, e.child, { id: r.i });
  }
  throw new q(r.c);
}
function Lo(e, r) {
  let t = "", n = false;
  return r.f.t !== 4 && (O(e.base, r.f.i), t = "(" + f(e, r.f) + ",", n = true), t += ne(e, r.i, "(" + St + ")(" + p$1(e, r.f.i) + ")"), n && (t += ")"), t;
}
function Yo(e, r) {
  return L(e, r.a[0]) + "(" + f(e, r.a[1]) + ")";
}
function Wo(e, r) {
  let t = r.a[0], n = r.a[1], a = e.base, s = "";
  t.t !== 4 && (O(a, t.i), s += "(" + f(e, t)), n.t !== 4 && (O(a, n.i), s += (s ? "," : "(") + f(e, n)), s && (s += ",");
  let i = ne(e, r.i, "(" + pt + ")(" + p$1(e, n.i) + "," + p$1(e, t.i) + ")");
  return s ? s + i + ")" : i;
}
function Ko(e, r) {
  return L(e, r.a[0]) + "(" + f(e, r.a[1]) + ")";
}
function Go(e, r) {
  let t = ne(e, r.i, L(e, r.f) + "()"), n = r.a.length;
  if (n) {
    let a = f(e, r.a[0]);
    for (let s = 1; s < n; s++) a += "," + f(e, r.a[s]);
    return "(" + t + "," + a + "," + p$1(e, r.i) + ")";
  }
  return t;
}
function Zo(e, r) {
  return p$1(e, r.i) + ".next(" + f(e, r.f) + ")";
}
function Jo(e, r) {
  return p$1(e, r.i) + ".throw(" + f(e, r.f) + ")";
}
function Ho(e, r) {
  return p$1(e, r.i) + ".return(" + f(e, r.f) + ")";
}
function $o(e, r) {
  switch (r.t) {
    case 17:
      return Zr[r.s];
    case 18:
      return Co(r);
    case 9:
      return vo(e, r);
    case 10:
      return bo(e, r);
    case 11:
      return Po(e, r);
    case 5:
      return xo(r);
    case 6:
      return To(r);
    case 7:
      return Oo(e, r);
    case 8:
      return ho(e, r);
    case 19:
      return wo(e, r);
    case 16:
    case 15:
      return zo(e, r);
    case 20:
      return ko(e, r);
    case 14:
      return _o(e, r);
    case 13:
      return Do(e, r);
    case 12:
      return Fo(e, r);
    case 21:
      return Bo(e, r);
    case 22:
      return Vo(e, r);
    case 25:
      return Uo(e, r);
    case 26:
      return Nt[r.s];
    default:
      throw new k(r);
  }
}
function f(e, r) {
  switch (r.t) {
    case 2:
      return Hr[r.s];
    case 0:
      return "" + r.s;
    case 1:
      return '"' + r.s + '"';
    case 3:
      return r.s + "n";
    case 4:
      return p$1(e, r.i);
    case 23:
      return Mo(e, r);
    case 24:
      return jo(e, r);
    case 27:
      return Lo(e, r);
    case 28:
      return Yo(e, r);
    case 29:
      return Wo(e, r);
    case 30:
      return Ko(e, r);
    case 31:
      return Go(e, r);
    case 32:
      return Zo(e, r);
    case 33:
      return Jo(e, r);
    case 34:
      return Ho(e, r);
    default:
      return ne(e, r.i, $o(e, r));
  }
}
function ar(e, r) {
  let t = f(e, r), n = r.i;
  if (n == null) return t;
  let a = Gt(e.base), s = p$1(e, n), i = e.state.scopeId, l = i == null ? "" : ie, u2 = a ? "(" + t + "," + a + s + ")" : t;
  if (l === "") return r.t === 10 && !a ? "(" + u2 + ")" : u2;
  let g = i == null ? "()" : "(" + ie + '["' + y(i) + '"])';
  return "(" + rr([l], u2) + ")" + g;
}
var Vr = class {
  constructor(r) {
    this._p = r;
  }
  parse(r) {
    return R(this._p, r);
  }
}, Mr = class {
  constructor(r) {
    this._p = r;
  }
  parse(r) {
    return R(this._p, r);
  }
  parseWithError(r) {
    return Y(this._p, r);
  }
  isAlive() {
    return this._p.state.alive;
  }
  pushPendingState() {
    Wr(this._p);
  }
  popPendingState() {
    me(this._p);
  }
  onParse(r) {
    oe(this._p, r);
  }
  onError(r) {
    Lr(this._p, r);
  }
};
function qo(e) {
  return { alive: true, pending: 0, initial: true, buffer: [], onParse: e.onParse, onError: e.onError, onDone: e.onDone };
}
function jr(e) {
  return { type: 2, base: ce(2, e), child: void 0, state: qo(e) };
}
function Xo(e, r) {
  let t = [];
  for (let n = 0, a = r.length; n < a; n++) n in r && (t[n] = R(e, r[n]));
  return t;
}
function Qo(e, r, t) {
  return Oe(r, t, Xo(e, t));
}
function Ur(e, r) {
  let t = Object.entries(r), n = [], a = [];
  for (let s = 0, i = t.length; s < i; s++) n.push(y(t[s][0])), a.push(R(e, t[s][1]));
  return A in r && (n.push(I(e.base, A)), a.push(Be(He(e.base), R(e, Ze(r))))), b in r && (n.push(I(e.base, b)), a.push(Ve($e(e.base), R(e, e.type === 1 ? Q() : Ge(r))))), x in r && (n.push(I(e.base, x)), a.push($(r[x]))), P$1 in r && (n.push(I(e.base, P$1)), a.push(r[P$1] ? Z : J)), { k: n, v: a, s: n.length };
}
function Br(e, r, t, n) {
  return qe(r, t, n, Ur(e, t));
}
function ea(e, r, t) {
  return he(r, R(e, t.valueOf()));
}
function ra(e, r, t) {
  return we(r, t, R(e, t.buffer));
}
function ta(e, r, t) {
  return ze(r, t, R(e, t.buffer));
}
function na(e, r, t) {
  return ke(r, t, R(e, t.buffer));
}
function Ht(e, r, t) {
  let n = H(t, e.base.features);
  return _e(r, t, n ? Ur(e, n) : o);
}
function oa(e, r, t) {
  let n = H(t, e.base.features);
  return De(r, t, n ? Ur(e, n) : o);
}
function aa(e, r, t) {
  let n = [], a = [];
  for (let [s, i] of t.entries()) n.push(R(e, s)), a.push(R(e, i));
  return Xe(e.base, r, n, a, t.size);
}
function sa(e, r, t) {
  let n = [];
  for (let a of t.keys()) n.push(R(e, a));
  return Fe(r, t.size, n);
}
function ia(e, r, t) {
  let n = Me(r, w$1(e.base, 4), []);
  return e.type === 1 || (Wr(e), t.on({ next: (a) => {
    if (e.state.alive) {
      let s = Y(e, a);
      s && oe(e, je(r, s));
    }
  }, throw: (a) => {
    if (e.state.alive) {
      let s = Y(e, a);
      s && oe(e, Ue(r, s));
    }
    me(e);
  }, return: (a) => {
    if (e.state.alive) {
      let s = Y(e, a);
      s && oe(e, Le(r, s));
    }
    me(e);
  } })), n;
}
function la(e, r) {
  if (this.state.alive) {
    let t = Y(this, r);
    t && oe(this, c(23, e, o, o, o, o, o, o, [w$1(this.base, 2), t], o, o, o)), me(this);
  }
}
function ua(e, r) {
  if (this.state.alive) {
    let t = Y(this, r);
    t && oe(this, c(24, e, o, o, o, o, o, o, [w$1(this.base, 3), t], o, o, o));
  }
  me(this);
}
function ca(e, r, t) {
  let n = Er(e.base, {});
  return e.type === 2 && (Wr(e), t.then(la.bind(e, n), ua.bind(e, n))), At(e.base, r, n);
}
function fa(e, r, t, n) {
  for (let a = 0, s = n.length; a < s; a++) {
    let i = n[a];
    if (i.parse.sync && i.test(t)) return e.child == null && (e.child = new Vr(e)), le(r, i.tag, i.parse.sync(t, e.child, { id: r }));
  }
}
function Sa(e, r, t, n) {
  for (let a = 0, s = n.length; a < s; a++) {
    let i = n[a];
    if (i.parse.stream && i.test(t)) return e.child == null && (e.child = new Mr(e)), le(r, i.tag, i.parse.stream(t, e.child, { id: r }));
  }
}
function $t(e, r, t) {
  let n = e.base.plugins;
  if (n) return e.type === 1 ? fa(e, r, t, n) : Sa(e, r, t, n);
}
function pa(e, r, t, n) {
  switch (n) {
    case Object:
      return Br(e, r, t, false);
    case void 0:
      return Br(e, r, t, true);
    case Date:
      return xe(r, t);
    case RegExp:
      return Te(r, t);
    case Error:
    case EvalError:
    case RangeError:
    case ReferenceError:
    case SyntaxError:
    case TypeError:
    case URIError:
      return Ht(e, r, t);
    case Number:
    case Boolean:
    case String:
    case BigInt:
      return ea(e, r, t);
    case ArrayBuffer:
      return Qe(e.base, r, t);
    case Int8Array:
    case Int16Array:
    case Int32Array:
    case Uint8Array:
    case Uint16Array:
    case Uint32Array:
    case Uint8ClampedArray:
    case Float32Array:
    case Float64Array:
      return ra(e, r, t);
    case DataView:
      return na(e, r, t);
    case Map:
      return aa(e, r, t);
    case Set:
      return sa(e, r, t);
  }
  if (n === Promise || t instanceof Promise) return ca(e, r, t);
  let a = e.base.features;
  if (a & 16) switch (n) {
    case BigInt64Array:
    case BigUint64Array:
      return ta(e, r, t);
  }
  if (a & 1 && typeof AggregateError != "undefined" && (n === AggregateError || t instanceof AggregateError)) return oa(e, r, t);
  if (t instanceof Error) return Ht(e, r, t);
  if (A in t || b in t) return Br(e, r, t, !!n);
  throw new E(t);
}
function da(e, r, t) {
  if (Array.isArray(t)) return Qo(e, r, t);
  if (Ke(t)) return ia(e, r, t);
  let n = t.constructor;
  if (n === j) return R(e, t.replacement);
  let a = $t(e, r, t);
  return a || pa(e, r, t, n);
}
function ma(e, r) {
  let t = U(e.base, r);
  if (t.type !== 0) return t.value;
  let n = $t(e, t.value, r);
  if (n) return n;
  throw new E(r);
}
function R(e, r) {
  switch (typeof r) {
    case "boolean":
      return r ? Z : J;
    case "undefined":
      return Ne;
    case "string":
      return $(r);
    case "number":
      return Ee(r);
    case "bigint":
      return Pe(r);
    case "object": {
      if (r) {
        let t = U(e.base, r);
        return t.type === 0 ? da(e, t.value, r) : t.value;
      }
      return Ce;
    }
    case "symbol":
      return I(e.base, r);
    case "function":
      return ma(e, r);
    default:
      throw new E(r);
  }
}
function oe(e, r) {
  e.state.initial ? e.state.buffer.push(r) : Yr(e, r, false);
}
function Lr(e, r) {
  if (e.state.onError) e.state.onError(r);
  else throw r instanceof h ? r : new h(r);
}
function qt(e) {
  e.state.onDone && e.state.onDone();
}
function Yr(e, r, t) {
  try {
    e.state.onParse(r, t);
  } catch (n) {
    Lr(e, n);
  }
}
function Wr(e) {
  e.state.pending++;
}
function me(e) {
  --e.state.pending <= 0 && qt(e);
}
function Y(e, r) {
  try {
    return R(e, r);
  } catch (t) {
    return Lr(e, t), o;
  }
}
function Kr(e, r) {
  let t = Y(e, r);
  t && (Yr(e, t, true), e.state.initial = false, ga(e, e.state), e.state.pending <= 0 && sr(e));
}
function ga(e, r) {
  for (let t = 0, n = r.buffer.length; t < n; t++) Yr(e, r.buffer[t], false);
}
function sr(e) {
  e.state.alive && (qt(e), e.state.alive = false);
}
async function ki(e, r = {}) {
  let t = v(r.plugins), n = ee$1(2, { plugins: t, disabledFeatures: r.disabledFeatures, refs: r.refs });
  return await re$1(n, e);
}
function Xt(e, r) {
  let t = v(r.plugins), n = jr({ plugins: t, refs: r.refs, disabledFeatures: r.disabledFeatures, onParse(a, s) {
    let i = nr({ plugins: t, features: n.base.features, scopeId: r.scopeId, markedRefs: n.base.marked }), l;
    try {
      l = ar(i, a);
    } catch (u2) {
      r.onError && r.onError(u2);
      return;
    }
    r.onSerialize(l, s);
  }, onError: r.onError, onDone: r.onDone });
  return Kr(n, e), sr.bind(null, n);
}
function _i(e, r) {
  let t = v(r.plugins), n = jr({ plugins: t, refs: r.refs, disabledFeatures: r.disabledFeatures, onParse: r.onParse, onError: r.onError, onDone: r.onDone });
  return Kr(n, e), sr.bind(null, n);
}
function Xi(e, r = {}) {
  let t = v(r.plugins), n = Ot({ plugins: t, markedRefs: e.m });
  return er(n, e.t);
}
const GLOBAL_TSR = "$_TSR";
function createSerializationAdapter(opts) {
  return opts;
}
function makeSsrSerovalPlugin(serializationAdapter, options) {
  return _s({
    tag: "$TSR/t/" + serializationAdapter.key,
    test: serializationAdapter.test,
    parse: {
      stream(value, ctx) {
        return ctx.parse(serializationAdapter.toSerializable(value));
      }
    },
    serialize(node, ctx) {
      options.didRun = true;
      return GLOBAL_TSR + '.t.get("' + serializationAdapter.key + '")(' + ctx.serialize(node) + ")";
    },
    // we never deserialize on the server during SSR
    deserialize: void 0
  });
}
function makeSerovalPlugin(serializationAdapter) {
  return _s({
    tag: "$TSR/t/" + serializationAdapter.key,
    test: serializationAdapter.test,
    parse: {
      sync(value, ctx) {
        return ctx.parse(serializationAdapter.toSerializable(value));
      },
      async async(value, ctx) {
        return await ctx.parse(serializationAdapter.toSerializable(value));
      },
      stream(value, ctx) {
        return ctx.parse(serializationAdapter.toSerializable(value));
      }
    },
    // we don't generate JS code outside of SSR (for now)
    serialize: void 0,
    deserialize(node, ctx) {
      return serializationAdapter.fromSerializable(ctx.deserialize(node));
    }
  });
}
var p = {}, P = (e) => new ReadableStream({ start: (r) => {
  e.on({ next: (a) => {
    try {
      r.enqueue(a);
    } catch (t) {
    }
  }, throw: (a) => {
    r.error(a);
  }, return: () => {
    try {
      r.close();
    } catch (a) {
    }
  } });
} }), ee = _s({ tag: "seroval-plugins/web/ReadableStreamFactory", test(e) {
  return e === p;
}, parse: { sync() {
}, async async() {
  return await Promise.resolve(void 0);
}, stream() {
} }, serialize() {
  return P.toString();
}, deserialize() {
  return p;
} });
function w(e) {
  let r = Q(), a = e.getReader();
  async function t() {
    try {
      let n = await a.read();
      n.done ? r.return(n.value) : (r.next(n.value), await t());
    } catch (n) {
      r.throw(n);
    }
  }
  return t().catch(() => {
  }), r;
}
var re = _s({ tag: "seroval/plugins/web/ReadableStream", extends: [ee], test(e) {
  return typeof ReadableStream == "undefined" ? false : e instanceof ReadableStream;
}, parse: { sync(e, r) {
  return { factory: r.parse(p), stream: r.parse(Q()) };
}, async async(e, r) {
  return { factory: await r.parse(p), stream: await r.parse(w(e)) };
}, stream(e, r) {
  return { factory: r.parse(p), stream: r.parse(w(e)) };
} }, serialize(e, r) {
  return "(" + r.serialize(e.factory) + ")(" + r.serialize(e.stream) + ")";
}, deserialize(e, r) {
  let a = r.deserialize(e.stream);
  return P(a);
} }), u = re;
const ShallowErrorPlugin = /* @__PURE__ */ _s({
  tag: "$TSR/Error",
  test(value) {
    return value instanceof Error;
  },
  parse: {
    sync(value, ctx) {
      return {
        message: ctx.parse(value.message)
      };
    },
    async async(value, ctx) {
      return {
        message: await ctx.parse(value.message)
      };
    },
    stream(value, ctx) {
      return {
        message: ctx.parse(value.message)
      };
    }
  },
  serialize(node, ctx) {
    return "new Error(" + ctx.serialize(node.message) + ")";
  },
  deserialize(node, ctx) {
    return new Error(ctx.deserialize(node.message));
  }
});
const defaultSerovalPlugins = [
  ShallowErrorPlugin,
  // ReadableStreamNode is not exported by seroval
  u
];
function CatchBoundary(props) {
  const errorComponent = props.errorComponent ?? ErrorComponent;
  return /* @__PURE__ */ jsx(
    CatchBoundaryImpl,
    {
      getResetKey: props.getResetKey,
      onCatch: props.onCatch,
      children: ({ error, reset }) => {
        if (error) {
          return React.createElement(errorComponent, {
            error,
            reset
          });
        }
        return props.children;
      }
    }
  );
}
class CatchBoundaryImpl extends React.Component {
  constructor() {
    super(...arguments);
    this.state = { error: null };
  }
  static getDerivedStateFromProps(props) {
    return { resetKey: props.getResetKey() };
  }
  static getDerivedStateFromError(error) {
    return { error };
  }
  reset() {
    this.setState({ error: null });
  }
  componentDidUpdate(prevProps, prevState) {
    if (prevState.error && prevState.resetKey !== this.state.resetKey) {
      this.reset();
    }
  }
  componentDidCatch(error, errorInfo) {
    if (this.props.onCatch) {
      this.props.onCatch(error, errorInfo);
    }
  }
  render() {
    return this.props.children({
      error: this.state.resetKey !== this.props.getResetKey() ? null : this.state.error,
      reset: () => {
        this.reset();
      }
    });
  }
}
function ErrorComponent({ error }) {
  const [show, setShow] = React.useState(process.env.NODE_ENV !== "production");
  return /* @__PURE__ */ jsxs("div", { style: { padding: ".5rem", maxWidth: "100%" }, children: [
    /* @__PURE__ */ jsxs("div", { style: { display: "flex", alignItems: "center", gap: ".5rem" }, children: [
      /* @__PURE__ */ jsx("strong", { style: { fontSize: "1rem" }, children: "Something went wrong!" }),
      /* @__PURE__ */ jsx(
        "button",
        {
          style: {
            appearance: "none",
            fontSize: ".6em",
            border: "1px solid currentColor",
            padding: ".1rem .2rem",
            fontWeight: "bold",
            borderRadius: ".25rem"
          },
          onClick: () => setShow((d) => !d),
          children: show ? "Hide Error" : "Show Error"
        }
      )
    ] }),
    /* @__PURE__ */ jsx("div", { style: { height: ".25rem" } }),
    show ? /* @__PURE__ */ jsx("div", { children: /* @__PURE__ */ jsx(
      "pre",
      {
        style: {
          fontSize: ".7em",
          border: "1px solid red",
          borderRadius: ".25rem",
          padding: ".3rem",
          color: "red",
          overflow: "auto"
        },
        children: error.message ? /* @__PURE__ */ jsx("code", { children: error.message }) : null
      }
    ) }) : null
  ] });
}
function ClientOnly({ children, fallback = null }) {
  return useHydrated() ? /* @__PURE__ */ jsx(React__default.Fragment, { children }) : /* @__PURE__ */ jsx(React__default.Fragment, { children: fallback });
}
function useHydrated() {
  return React__default.useSyncExternalStore(
    subscribe,
    () => true,
    () => false
  );
}
function subscribe() {
  return () => {
  };
}
var isProduction = process.env.NODE_ENV === "production";
function warning(condition, message) {
  if (!isProduction) {
    if (condition) {
      return;
    }
    var text = "Warning: " + message;
    if (typeof console !== "undefined") {
      console.warn(text);
    }
    try {
      throw Error(text);
    } catch (x2) {
    }
  }
}
var withSelector = { exports: {} };
var withSelector_production = {};
var shim = { exports: {} };
var useSyncExternalStoreShim_production = {};
var hasRequiredUseSyncExternalStoreShim_production;
function requireUseSyncExternalStoreShim_production() {
  if (hasRequiredUseSyncExternalStoreShim_production) return useSyncExternalStoreShim_production;
  hasRequiredUseSyncExternalStoreShim_production = 1;
  var React2 = React__default;
  function is(x2, y2) {
    return x2 === y2 && (0 !== x2 || 1 / x2 === 1 / y2) || x2 !== x2 && y2 !== y2;
  }
  var objectIs = "function" === typeof Object.is ? Object.is : is, useState = React2.useState, useEffect = React2.useEffect, useLayoutEffect2 = React2.useLayoutEffect, useDebugValue = React2.useDebugValue;
  function useSyncExternalStore$2(subscribe2, getSnapshot) {
    var value = getSnapshot(), _useState = useState({ inst: { value, getSnapshot } }), inst = _useState[0].inst, forceUpdate = _useState[1];
    useLayoutEffect2(
      function() {
        inst.value = value;
        inst.getSnapshot = getSnapshot;
        checkIfSnapshotChanged(inst) && forceUpdate({ inst });
      },
      [subscribe2, value, getSnapshot]
    );
    useEffect(
      function() {
        checkIfSnapshotChanged(inst) && forceUpdate({ inst });
        return subscribe2(function() {
          checkIfSnapshotChanged(inst) && forceUpdate({ inst });
        });
      },
      [subscribe2]
    );
    useDebugValue(value);
    return value;
  }
  function checkIfSnapshotChanged(inst) {
    var latestGetSnapshot = inst.getSnapshot;
    inst = inst.value;
    try {
      var nextValue = latestGetSnapshot();
      return !objectIs(inst, nextValue);
    } catch (error) {
      return true;
    }
  }
  function useSyncExternalStore$1(subscribe2, getSnapshot) {
    return getSnapshot();
  }
  var shim2 = "undefined" === typeof window || "undefined" === typeof window.document || "undefined" === typeof window.document.createElement ? useSyncExternalStore$1 : useSyncExternalStore$2;
  useSyncExternalStoreShim_production.useSyncExternalStore = void 0 !== React2.useSyncExternalStore ? React2.useSyncExternalStore : shim2;
  return useSyncExternalStoreShim_production;
}
var useSyncExternalStoreShim_development = {};
var hasRequiredUseSyncExternalStoreShim_development;
function requireUseSyncExternalStoreShim_development() {
  if (hasRequiredUseSyncExternalStoreShim_development) return useSyncExternalStoreShim_development;
  hasRequiredUseSyncExternalStoreShim_development = 1;
  "production" !== process.env.NODE_ENV && (function() {
    function is(x2, y2) {
      return x2 === y2 && (0 !== x2 || 1 / x2 === 1 / y2) || x2 !== x2 && y2 !== y2;
    }
    function useSyncExternalStore$2(subscribe2, getSnapshot) {
      didWarnOld18Alpha || void 0 === React2.startTransition || (didWarnOld18Alpha = true, console.error(
        "You are using an outdated, pre-release alpha of React 18 that does not support useSyncExternalStore. The use-sync-external-store shim will not work correctly. Upgrade to a newer pre-release."
      ));
      var value = getSnapshot();
      if (!didWarnUncachedGetSnapshot) {
        var cachedValue = getSnapshot();
        objectIs(value, cachedValue) || (console.error(
          "The result of getSnapshot should be cached to avoid an infinite loop"
        ), didWarnUncachedGetSnapshot = true);
      }
      cachedValue = useState({
        inst: { value, getSnapshot }
      });
      var inst = cachedValue[0].inst, forceUpdate = cachedValue[1];
      useLayoutEffect2(
        function() {
          inst.value = value;
          inst.getSnapshot = getSnapshot;
          checkIfSnapshotChanged(inst) && forceUpdate({ inst });
        },
        [subscribe2, value, getSnapshot]
      );
      useEffect(
        function() {
          checkIfSnapshotChanged(inst) && forceUpdate({ inst });
          return subscribe2(function() {
            checkIfSnapshotChanged(inst) && forceUpdate({ inst });
          });
        },
        [subscribe2]
      );
      useDebugValue(value);
      return value;
    }
    function checkIfSnapshotChanged(inst) {
      var latestGetSnapshot = inst.getSnapshot;
      inst = inst.value;
      try {
        var nextValue = latestGetSnapshot();
        return !objectIs(inst, nextValue);
      } catch (error) {
        return true;
      }
    }
    function useSyncExternalStore$1(subscribe2, getSnapshot) {
      return getSnapshot();
    }
    "undefined" !== typeof __REACT_DEVTOOLS_GLOBAL_HOOK__ && "function" === typeof __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStart && __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStart(Error());
    var React2 = React__default, objectIs = "function" === typeof Object.is ? Object.is : is, useState = React2.useState, useEffect = React2.useEffect, useLayoutEffect2 = React2.useLayoutEffect, useDebugValue = React2.useDebugValue, didWarnOld18Alpha = false, didWarnUncachedGetSnapshot = false, shim2 = "undefined" === typeof window || "undefined" === typeof window.document || "undefined" === typeof window.document.createElement ? useSyncExternalStore$1 : useSyncExternalStore$2;
    useSyncExternalStoreShim_development.useSyncExternalStore = void 0 !== React2.useSyncExternalStore ? React2.useSyncExternalStore : shim2;
    "undefined" !== typeof __REACT_DEVTOOLS_GLOBAL_HOOK__ && "function" === typeof __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStop && __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStop(Error());
  })();
  return useSyncExternalStoreShim_development;
}
var hasRequiredShim;
function requireShim() {
  if (hasRequiredShim) return shim.exports;
  hasRequiredShim = 1;
  if (process.env.NODE_ENV === "production") {
    shim.exports = requireUseSyncExternalStoreShim_production();
  } else {
    shim.exports = requireUseSyncExternalStoreShim_development();
  }
  return shim.exports;
}
var hasRequiredWithSelector_production;
function requireWithSelector_production() {
  if (hasRequiredWithSelector_production) return withSelector_production;
  hasRequiredWithSelector_production = 1;
  var React2 = React__default, shim2 = requireShim();
  function is(x2, y2) {
    return x2 === y2 && (0 !== x2 || 1 / x2 === 1 / y2) || x2 !== x2 && y2 !== y2;
  }
  var objectIs = "function" === typeof Object.is ? Object.is : is, useSyncExternalStore = shim2.useSyncExternalStore, useRef2 = React2.useRef, useEffect = React2.useEffect, useMemo = React2.useMemo, useDebugValue = React2.useDebugValue;
  withSelector_production.useSyncExternalStoreWithSelector = function(subscribe2, getSnapshot, getServerSnapshot, selector, isEqual) {
    var instRef = useRef2(null);
    if (null === instRef.current) {
      var inst = { hasValue: false, value: null };
      instRef.current = inst;
    } else inst = instRef.current;
    instRef = useMemo(
      function() {
        function memoizedSelector(nextSnapshot) {
          if (!hasMemo) {
            hasMemo = true;
            memoizedSnapshot = nextSnapshot;
            nextSnapshot = selector(nextSnapshot);
            if (void 0 !== isEqual && inst.hasValue) {
              var currentSelection = inst.value;
              if (isEqual(currentSelection, nextSnapshot))
                return memoizedSelection = currentSelection;
            }
            return memoizedSelection = nextSnapshot;
          }
          currentSelection = memoizedSelection;
          if (objectIs(memoizedSnapshot, nextSnapshot)) return currentSelection;
          var nextSelection = selector(nextSnapshot);
          if (void 0 !== isEqual && isEqual(currentSelection, nextSelection))
            return memoizedSnapshot = nextSnapshot, currentSelection;
          memoizedSnapshot = nextSnapshot;
          return memoizedSelection = nextSelection;
        }
        var hasMemo = false, memoizedSnapshot, memoizedSelection, maybeGetServerSnapshot = void 0 === getServerSnapshot ? null : getServerSnapshot;
        return [
          function() {
            return memoizedSelector(getSnapshot());
          },
          null === maybeGetServerSnapshot ? void 0 : function() {
            return memoizedSelector(maybeGetServerSnapshot());
          }
        ];
      },
      [getSnapshot, getServerSnapshot, selector, isEqual]
    );
    var value = useSyncExternalStore(subscribe2, instRef[0], instRef[1]);
    useEffect(
      function() {
        inst.hasValue = true;
        inst.value = value;
      },
      [value]
    );
    useDebugValue(value);
    return value;
  };
  return withSelector_production;
}
var withSelector_development = {};
var hasRequiredWithSelector_development;
function requireWithSelector_development() {
  if (hasRequiredWithSelector_development) return withSelector_development;
  hasRequiredWithSelector_development = 1;
  "production" !== process.env.NODE_ENV && (function() {
    function is(x2, y2) {
      return x2 === y2 && (0 !== x2 || 1 / x2 === 1 / y2) || x2 !== x2 && y2 !== y2;
    }
    "undefined" !== typeof __REACT_DEVTOOLS_GLOBAL_HOOK__ && "function" === typeof __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStart && __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStart(Error());
    var React2 = React__default, shim2 = requireShim(), objectIs = "function" === typeof Object.is ? Object.is : is, useSyncExternalStore = shim2.useSyncExternalStore, useRef2 = React2.useRef, useEffect = React2.useEffect, useMemo = React2.useMemo, useDebugValue = React2.useDebugValue;
    withSelector_development.useSyncExternalStoreWithSelector = function(subscribe2, getSnapshot, getServerSnapshot, selector, isEqual) {
      var instRef = useRef2(null);
      if (null === instRef.current) {
        var inst = { hasValue: false, value: null };
        instRef.current = inst;
      } else inst = instRef.current;
      instRef = useMemo(
        function() {
          function memoizedSelector(nextSnapshot) {
            if (!hasMemo) {
              hasMemo = true;
              memoizedSnapshot = nextSnapshot;
              nextSnapshot = selector(nextSnapshot);
              if (void 0 !== isEqual && inst.hasValue) {
                var currentSelection = inst.value;
                if (isEqual(currentSelection, nextSnapshot))
                  return memoizedSelection = currentSelection;
              }
              return memoizedSelection = nextSnapshot;
            }
            currentSelection = memoizedSelection;
            if (objectIs(memoizedSnapshot, nextSnapshot))
              return currentSelection;
            var nextSelection = selector(nextSnapshot);
            if (void 0 !== isEqual && isEqual(currentSelection, nextSelection))
              return memoizedSnapshot = nextSnapshot, currentSelection;
            memoizedSnapshot = nextSnapshot;
            return memoizedSelection = nextSelection;
          }
          var hasMemo = false, memoizedSnapshot, memoizedSelection, maybeGetServerSnapshot = void 0 === getServerSnapshot ? null : getServerSnapshot;
          return [
            function() {
              return memoizedSelector(getSnapshot());
            },
            null === maybeGetServerSnapshot ? void 0 : function() {
              return memoizedSelector(maybeGetServerSnapshot());
            }
          ];
        },
        [getSnapshot, getServerSnapshot, selector, isEqual]
      );
      var value = useSyncExternalStore(subscribe2, instRef[0], instRef[1]);
      useEffect(
        function() {
          inst.hasValue = true;
          inst.value = value;
        },
        [value]
      );
      useDebugValue(value);
      return value;
    };
    "undefined" !== typeof __REACT_DEVTOOLS_GLOBAL_HOOK__ && "function" === typeof __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStop && __REACT_DEVTOOLS_GLOBAL_HOOK__.registerInternalModuleStop(Error());
  })();
  return withSelector_development;
}
var hasRequiredWithSelector;
function requireWithSelector() {
  if (hasRequiredWithSelector) return withSelector.exports;
  hasRequiredWithSelector = 1;
  if (process.env.NODE_ENV === "production") {
    withSelector.exports = requireWithSelector_production();
  } else {
    withSelector.exports = requireWithSelector_development();
  }
  return withSelector.exports;
}
var withSelectorExports = requireWithSelector();
function useStore(store, selector = (d) => d, options = {}) {
  const equal = options.equal ?? shallow;
  const slice = withSelectorExports.useSyncExternalStoreWithSelector(
    store.subscribe,
    () => store.state,
    () => store.state,
    selector,
    equal
  );
  return slice;
}
function shallow(objA, objB) {
  if (Object.is(objA, objB)) {
    return true;
  }
  if (typeof objA !== "object" || objA === null || typeof objB !== "object" || objB === null) {
    return false;
  }
  if (objA instanceof Map && objB instanceof Map) {
    if (objA.size !== objB.size) return false;
    for (const [k2, v2] of objA) {
      if (!objB.has(k2) || !Object.is(v2, objB.get(k2))) return false;
    }
    return true;
  }
  if (objA instanceof Set && objB instanceof Set) {
    if (objA.size !== objB.size) return false;
    for (const v2 of objA) {
      if (!objB.has(v2)) return false;
    }
    return true;
  }
  if (objA instanceof Date && objB instanceof Date) {
    if (objA.getTime() !== objB.getTime()) return false;
    return true;
  }
  const keysA = getOwnKeys(objA);
  if (keysA.length !== getOwnKeys(objB).length) {
    return false;
  }
  for (let i = 0; i < keysA.length; i++) {
    if (!Object.prototype.hasOwnProperty.call(objB, keysA[i]) || !Object.is(objA[keysA[i]], objB[keysA[i]])) {
      return false;
    }
  }
  return true;
}
function getOwnKeys(obj) {
  return Object.keys(obj).concat(
    Object.getOwnPropertySymbols(obj)
  );
}
const routerContext = React.createContext(null);
function getRouterContext() {
  if (typeof document === "undefined") {
    return routerContext;
  }
  if (window.__TSR_ROUTER_CONTEXT__) {
    return window.__TSR_ROUTER_CONTEXT__;
  }
  window.__TSR_ROUTER_CONTEXT__ = routerContext;
  return routerContext;
}
function useRouter(opts) {
  const value = React.useContext(getRouterContext());
  warning(
    !((opts?.warn ?? true) && !value),
    "useRouter must be used inside a <RouterProvider> component!"
  );
  return value;
}
function useRouterState(opts) {
  const contextRouter = useRouter({
    warn: opts?.router === void 0
  });
  const router = opts?.router || contextRouter;
  const previousResult = useRef(void 0);
  return useStore(router.__store, (state) => {
    if (opts?.select) {
      if (opts.structuralSharing ?? router.options.defaultStructuralSharing) {
        const newSlice = replaceEqualDeep(
          previousResult.current,
          opts.select(state)
        );
        previousResult.current = newSlice;
        return newSlice;
      }
      return opts.select(state);
    }
    return state;
  });
}
const matchContext = React.createContext(void 0);
const dummyMatchContext = React.createContext(
  void 0
);
const useLayoutEffect = typeof window !== "undefined" ? React.useLayoutEffect : React.useEffect;
function usePrevious(value) {
  const ref = React.useRef({
    value,
    prev: null
  });
  const current = ref.current.value;
  if (value !== current) {
    ref.current = {
      value,
      prev: current
    };
  }
  return ref.current.prev;
}
function useIntersectionObserver(ref, callback, intersectionObserverOptions = {}, options = {}) {
  React.useEffect(() => {
    if (!ref.current || options.disabled || typeof IntersectionObserver !== "function") {
      return;
    }
    const observer = new IntersectionObserver(([entry]) => {
      callback(entry);
    }, intersectionObserverOptions);
    observer.observe(ref.current);
    return () => {
      observer.disconnect();
    };
  }, [callback, intersectionObserverOptions, options.disabled, ref]);
}
function useForwardedRef(ref) {
  const innerRef = React.useRef(null);
  React.useImperativeHandle(ref, () => innerRef.current, []);
  return innerRef;
}
function Transitioner() {
  const router = useRouter();
  const mountLoadForRouter = React.useRef({ router, mounted: false });
  const [isTransitioning, setIsTransitioning] = React.useState(false);
  const { hasPendingMatches, isLoading } = useRouterState({
    select: (s) => ({
      isLoading: s.isLoading,
      hasPendingMatches: s.matches.some((d) => d.status === "pending")
    }),
    structuralSharing: true
  });
  const previousIsLoading = usePrevious(isLoading);
  const isAnyPending = isLoading || isTransitioning || hasPendingMatches;
  const previousIsAnyPending = usePrevious(isAnyPending);
  const isPagePending = isLoading || hasPendingMatches;
  const previousIsPagePending = usePrevious(isPagePending);
  router.startTransition = (fn2) => {
    setIsTransitioning(true);
    React.startTransition(() => {
      fn2();
      setIsTransitioning(false);
    });
  };
  React.useEffect(() => {
    const unsub = router.history.subscribe(router.load);
    const nextLocation = router.buildLocation({
      to: router.latestLocation.pathname,
      search: true,
      params: true,
      hash: true,
      state: true,
      _includeValidateSearch: true
    });
    if (trimPathRight(router.latestLocation.href) !== trimPathRight(nextLocation.href)) {
      router.commitLocation({ ...nextLocation, replace: true });
    }
    return () => {
      unsub();
    };
  }, [router, router.history]);
  useLayoutEffect(() => {
    if (
      // if we are hydrating from SSR, loading is triggered in ssr-client
      typeof window !== "undefined" && router.ssr || mountLoadForRouter.current.router === router && mountLoadForRouter.current.mounted
    ) {
      return;
    }
    mountLoadForRouter.current = { router, mounted: true };
    const tryLoad = async () => {
      try {
        await router.load();
      } catch (err) {
        console.error(err);
      }
    };
    tryLoad();
  }, [router]);
  useLayoutEffect(() => {
    if (previousIsLoading && !isLoading) {
      router.emit({
        type: "onLoad",
        // When the new URL has committed, when the new matches have been loaded into state.matches
        ...getLocationChangeInfo(router.state)
      });
    }
  }, [previousIsLoading, router, isLoading]);
  useLayoutEffect(() => {
    if (previousIsPagePending && !isPagePending) {
      router.emit({
        type: "onBeforeRouteMount",
        ...getLocationChangeInfo(router.state)
      });
    }
  }, [isPagePending, previousIsPagePending, router]);
  useLayoutEffect(() => {
    if (previousIsAnyPending && !isAnyPending) {
      const changeInfo = getLocationChangeInfo(router.state);
      router.emit({
        type: "onResolved",
        ...changeInfo
      });
      router.__store.setState((s) => ({
        ...s,
        status: "idle",
        resolvedLocation: s.location
      }));
      if (changeInfo.hrefChanged) {
        handleHashScroll(router);
      }
    }
  }, [isAnyPending, previousIsAnyPending, router]);
  return null;
}
function CatchNotFound(props) {
  const resetKey = useRouterState({
    select: (s) => `not-found-${s.location.pathname}-${s.status}`
  });
  return /* @__PURE__ */ jsx(
    CatchBoundary,
    {
      getResetKey: () => resetKey,
      onCatch: (error, errorInfo) => {
        if (isNotFound(error)) {
          props.onCatch?.(error, errorInfo);
        } else {
          throw error;
        }
      },
      errorComponent: ({ error }) => {
        if (isNotFound(error)) {
          return props.fallback?.(error);
        } else {
          throw error;
        }
      },
      children: props.children
    }
  );
}
function DefaultGlobalNotFound() {
  return /* @__PURE__ */ jsx("p", { children: "Not Found" });
}
function SafeFragment(props) {
  return /* @__PURE__ */ jsx(Fragment, { children: props.children });
}
function renderRouteNotFound(router, route, data) {
  if (!route.options.notFoundComponent) {
    if (router.options.defaultNotFoundComponent) {
      return /* @__PURE__ */ jsx(router.options.defaultNotFoundComponent, { ...data });
    }
    if (process.env.NODE_ENV === "development") {
      warning(
        route.options.notFoundComponent,
        `A notFoundError was encountered on the route with ID "${route.id}", but a notFoundComponent option was not configured, nor was a router level defaultNotFoundComponent configured. Consider configuring at least one of these to avoid TanStack Router's overly generic defaultNotFoundComponent (<p>Not Found</p>)`
      );
    }
    return /* @__PURE__ */ jsx(DefaultGlobalNotFound, {});
  }
  return /* @__PURE__ */ jsx(route.options.notFoundComponent, { ...data });
}
function ScriptOnce({ children }) {
  const router = useRouter();
  if (!router.isServer) {
    return null;
  }
  return /* @__PURE__ */ jsx(
    "script",
    {
      nonce: router.options.ssr?.nonce,
      className: "$tsr",
      dangerouslySetInnerHTML: {
        __html: children + ';typeof $_TSR !== "undefined" && $_TSR.c()'
      }
    }
  );
}
function ScrollRestoration() {
  const router = useRouter();
  if (!router.isScrollRestoring || !router.isServer) {
    return null;
  }
  if (typeof router.options.scrollRestoration === "function") {
    const shouldRestore = router.options.scrollRestoration({
      location: router.latestLocation
    });
    if (!shouldRestore) {
      return null;
    }
  }
  const getKey = router.options.getScrollRestorationKey || defaultGetScrollRestorationKey;
  const userKey = getKey(router.latestLocation);
  const resolvedKey = userKey !== defaultGetScrollRestorationKey(router.latestLocation) ? userKey : void 0;
  const restoreScrollOptions = {
    storageKey,
    shouldScrollRestoration: true
  };
  if (resolvedKey) {
    restoreScrollOptions.key = resolvedKey;
  }
  return /* @__PURE__ */ jsx(
    ScriptOnce,
    {
      children: `(${restoreScroll.toString()})(${JSON.stringify(restoreScrollOptions)})`
    }
  );
}
const Match = React.memo(function MatchImpl({
  matchId
}) {
  const router = useRouter();
  const matchState = useRouterState({
    select: (s) => {
      const match = s.matches.find((d) => d.id === matchId);
      invariant(
        match,
        `Could not find match for matchId "${matchId}". Please file an issue!`
      );
      return {
        routeId: match.routeId,
        ssr: match.ssr,
        _displayPending: match._displayPending
      };
    },
    structuralSharing: true
  });
  const route = router.routesById[matchState.routeId];
  const PendingComponent = route.options.pendingComponent ?? router.options.defaultPendingComponent;
  const pendingElement = PendingComponent ? /* @__PURE__ */ jsx(PendingComponent, {}) : null;
  const routeErrorComponent = route.options.errorComponent ?? router.options.defaultErrorComponent;
  const routeOnCatch = route.options.onCatch ?? router.options.defaultOnCatch;
  const routeNotFoundComponent = route.isRoot ? (
    // If it's the root route, use the globalNotFound option, with fallback to the notFoundRoute's component
    route.options.notFoundComponent ?? router.options.notFoundRoute?.options.component
  ) : route.options.notFoundComponent;
  const resolvedNoSsr = matchState.ssr === false || matchState.ssr === "data-only";
  const ResolvedSuspenseBoundary = (
    // If we're on the root route, allow forcefully wrapping in suspense
    (!route.isRoot || route.options.wrapInSuspense || resolvedNoSsr) && (route.options.wrapInSuspense ?? PendingComponent ?? (route.options.errorComponent?.preload || resolvedNoSsr)) ? React.Suspense : SafeFragment
  );
  const ResolvedCatchBoundary = routeErrorComponent ? CatchBoundary : SafeFragment;
  const ResolvedNotFoundBoundary = routeNotFoundComponent ? CatchNotFound : SafeFragment;
  const resetKey = useRouterState({
    select: (s) => s.loadedAt
  });
  const parentRouteId = useRouterState({
    select: (s) => {
      const index = s.matches.findIndex((d) => d.id === matchId);
      return s.matches[index - 1]?.routeId;
    }
  });
  const ShellComponent = route.isRoot ? route.options.shellComponent ?? SafeFragment : SafeFragment;
  return /* @__PURE__ */ jsxs(ShellComponent, { children: [
    /* @__PURE__ */ jsx(matchContext.Provider, { value: matchId, children: /* @__PURE__ */ jsx(ResolvedSuspenseBoundary, { fallback: pendingElement, children: /* @__PURE__ */ jsx(
      ResolvedCatchBoundary,
      {
        getResetKey: () => resetKey,
        errorComponent: routeErrorComponent || ErrorComponent,
        onCatch: (error, errorInfo) => {
          if (isNotFound(error)) throw error;
          warning(false, `Error in route match: ${matchId}`);
          routeOnCatch?.(error, errorInfo);
        },
        children: /* @__PURE__ */ jsx(
          ResolvedNotFoundBoundary,
          {
            fallback: (error) => {
              if (!routeNotFoundComponent || error.routeId && error.routeId !== matchState.routeId || !error.routeId && !route.isRoot)
                throw error;
              return React.createElement(routeNotFoundComponent, error);
            },
            children: resolvedNoSsr || matchState._displayPending ? /* @__PURE__ */ jsx(ClientOnly, { fallback: pendingElement, children: /* @__PURE__ */ jsx(MatchInner, { matchId }) }) : /* @__PURE__ */ jsx(MatchInner, { matchId })
          }
        )
      }
    ) }) }),
    parentRouteId === rootRouteId && router.options.scrollRestoration ? /* @__PURE__ */ jsxs(Fragment, { children: [
      /* @__PURE__ */ jsx(OnRendered, {}),
      /* @__PURE__ */ jsx(ScrollRestoration, {})
    ] }) : null
  ] });
});
function OnRendered() {
  const router = useRouter();
  const prevLocationRef = React.useRef(
    void 0
  );
  return /* @__PURE__ */ jsx(
    "script",
    {
      suppressHydrationWarning: true,
      ref: (el) => {
        if (el && (prevLocationRef.current === void 0 || prevLocationRef.current.href !== router.latestLocation.href)) {
          router.emit({
            type: "onRendered",
            ...getLocationChangeInfo(router.state)
          });
          prevLocationRef.current = router.latestLocation;
        }
      }
    },
    router.latestLocation.state.__TSR_key
  );
}
const MatchInner = React.memo(function MatchInnerImpl({
  matchId
}) {
  const router = useRouter();
  const { match, key, routeId } = useRouterState({
    select: (s) => {
      const match2 = s.matches.find((d) => d.id === matchId);
      const routeId2 = match2.routeId;
      const remountFn = router.routesById[routeId2].options.remountDeps ?? router.options.defaultRemountDeps;
      const remountDeps = remountFn?.({
        routeId: routeId2,
        loaderDeps: match2.loaderDeps,
        params: match2._strictParams,
        search: match2._strictSearch
      });
      const key2 = remountDeps ? JSON.stringify(remountDeps) : void 0;
      return {
        key: key2,
        routeId: routeId2,
        match: {
          id: match2.id,
          status: match2.status,
          error: match2.error,
          _forcePending: match2._forcePending,
          _displayPending: match2._displayPending
        }
      };
    },
    structuralSharing: true
  });
  const route = router.routesById[routeId];
  const out = React.useMemo(() => {
    const Comp = route.options.component ?? router.options.defaultComponent;
    if (Comp) {
      return /* @__PURE__ */ jsx(Comp, {}, key);
    }
    return /* @__PURE__ */ jsx(Outlet, {});
  }, [key, route.options.component, router.options.defaultComponent]);
  if (match._displayPending) {
    throw router.getMatch(match.id)?._nonReactive.displayPendingPromise;
  }
  if (match._forcePending) {
    throw router.getMatch(match.id)?._nonReactive.minPendingPromise;
  }
  if (match.status === "pending") {
    const pendingMinMs = route.options.pendingMinMs ?? router.options.defaultPendingMinMs;
    if (pendingMinMs) {
      const routerMatch = router.getMatch(match.id);
      if (routerMatch && !routerMatch._nonReactive.minPendingPromise) {
        if (!router.isServer) {
          const minPendingPromise = createControlledPromise();
          routerMatch._nonReactive.minPendingPromise = minPendingPromise;
          setTimeout(() => {
            minPendingPromise.resolve();
            routerMatch._nonReactive.minPendingPromise = void 0;
          }, pendingMinMs);
        }
      }
    }
    throw router.getMatch(match.id)?._nonReactive.loadPromise;
  }
  if (match.status === "notFound") {
    invariant(isNotFound(match.error), "Expected a notFound error");
    return renderRouteNotFound(router, route, match.error);
  }
  if (match.status === "redirected") {
    invariant(isRedirect(match.error), "Expected a redirect error");
    throw router.getMatch(match.id)?._nonReactive.loadPromise;
  }
  if (match.status === "error") {
    if (router.isServer) {
      const RouteErrorComponent = (route.options.errorComponent ?? router.options.defaultErrorComponent) || ErrorComponent;
      return /* @__PURE__ */ jsx(
        RouteErrorComponent,
        {
          error: match.error,
          reset: void 0,
          info: {
            componentStack: ""
          }
        }
      );
    }
    throw match.error;
  }
  return out;
});
const Outlet = React.memo(function OutletImpl() {
  const router = useRouter();
  const matchId = React.useContext(matchContext);
  const routeId = useRouterState({
    select: (s) => s.matches.find((d) => d.id === matchId)?.routeId
  });
  const route = router.routesById[routeId];
  const parentGlobalNotFound = useRouterState({
    select: (s) => {
      const matches = s.matches;
      const parentMatch = matches.find((d) => d.id === matchId);
      invariant(
        parentMatch,
        `Could not find parent match for matchId "${matchId}"`
      );
      return parentMatch.globalNotFound;
    }
  });
  const childMatchId = useRouterState({
    select: (s) => {
      const matches = s.matches;
      const index = matches.findIndex((d) => d.id === matchId);
      return matches[index + 1]?.id;
    }
  });
  const pendingElement = router.options.defaultPendingComponent ? /* @__PURE__ */ jsx(router.options.defaultPendingComponent, {}) : null;
  if (parentGlobalNotFound) {
    return renderRouteNotFound(router, route, void 0);
  }
  if (!childMatchId) {
    return null;
  }
  const nextMatch = /* @__PURE__ */ jsx(Match, { matchId: childMatchId });
  if (routeId === rootRouteId) {
    return /* @__PURE__ */ jsx(React.Suspense, { fallback: pendingElement, children: nextMatch });
  }
  return nextMatch;
});
function Matches() {
  const router = useRouter();
  const rootRoute = router.routesById[rootRouteId];
  const PendingComponent = rootRoute.options.pendingComponent ?? router.options.defaultPendingComponent;
  const pendingElement = PendingComponent ? /* @__PURE__ */ jsx(PendingComponent, {}) : null;
  const ResolvedSuspense = router.isServer || typeof document !== "undefined" && router.ssr ? SafeFragment : React.Suspense;
  const inner = /* @__PURE__ */ jsxs(ResolvedSuspense, { fallback: pendingElement, children: [
    !router.isServer && /* @__PURE__ */ jsx(Transitioner, {}),
    /* @__PURE__ */ jsx(MatchesInner, {})
  ] });
  return router.options.InnerWrap ? /* @__PURE__ */ jsx(router.options.InnerWrap, { children: inner }) : inner;
}
function MatchesInner() {
  const router = useRouter();
  const matchId = useRouterState({
    select: (s) => {
      return s.matches[0]?.id;
    }
  });
  const resetKey = useRouterState({
    select: (s) => s.loadedAt
  });
  const matchComponent = matchId ? /* @__PURE__ */ jsx(Match, { matchId }) : null;
  return /* @__PURE__ */ jsx(matchContext.Provider, { value: matchId, children: router.options.disableGlobalCatchBoundary ? matchComponent : /* @__PURE__ */ jsx(
    CatchBoundary,
    {
      getResetKey: () => resetKey,
      errorComponent: ErrorComponent,
      onCatch: (error) => {
        warning(
          false,
          `The following error wasn't caught by any route! At the very least, consider setting an 'errorComponent' in your RootRoute!`
        );
        warning(false, error.message || error.toString());
      },
      children: matchComponent
    }
  ) });
}
function RouterContextProvider({
  router,
  children,
  ...rest
}) {
  if (Object.keys(rest).length > 0) {
    router.update({
      ...router.options,
      ...rest,
      context: {
        ...router.options.context,
        ...rest.context
      }
    });
  }
  const routerContext2 = getRouterContext();
  const provider = /* @__PURE__ */ jsx(routerContext2.Provider, { value: router, children });
  if (router.options.Wrap) {
    return /* @__PURE__ */ jsx(router.options.Wrap, { children: provider });
  }
  return provider;
}
function RouterProvider({ router, ...rest }) {
  return /* @__PURE__ */ jsx(RouterContextProvider, { router, ...rest, children: /* @__PURE__ */ jsx(Matches, {}) });
}
function StartServer(props) {
  return /* @__PURE__ */ jsx(RouterProvider, { router: props.router });
}
function splitSetCookieString$1(cookiesString) {
  if (Array.isArray(cookiesString)) {
    return cookiesString.flatMap((c2) => splitSetCookieString$1(c2));
  }
  if (typeof cookiesString !== "string") {
    return [];
  }
  const cookiesStrings = [];
  let pos = 0;
  let start;
  let ch;
  let lastComma;
  let nextStart;
  let cookiesSeparatorFound;
  const skipWhitespace = () => {
    while (pos < cookiesString.length && /\s/.test(cookiesString.charAt(pos))) {
      pos += 1;
    }
    return pos < cookiesString.length;
  };
  const notSpecialChar = () => {
    ch = cookiesString.charAt(pos);
    return ch !== "=" && ch !== ";" && ch !== ",";
  };
  while (pos < cookiesString.length) {
    start = pos;
    cookiesSeparatorFound = false;
    while (skipWhitespace()) {
      ch = cookiesString.charAt(pos);
      if (ch === ",") {
        lastComma = pos;
        pos += 1;
        skipWhitespace();
        nextStart = pos;
        while (pos < cookiesString.length && notSpecialChar()) {
          pos += 1;
        }
        if (pos < cookiesString.length && cookiesString.charAt(pos) === "=") {
          cookiesSeparatorFound = true;
          pos = nextStart;
          cookiesStrings.push(cookiesString.slice(start, lastComma));
          start = pos;
        } else {
          pos = lastComma + 1;
        }
      } else {
        pos += 1;
      }
    }
    if (!cookiesSeparatorFound || pos >= cookiesString.length) {
      cookiesStrings.push(cookiesString.slice(start));
    }
  }
  return cookiesStrings;
}
function toHeadersInstance(init) {
  if (init instanceof Headers) {
    return new Headers(init);
  } else if (Array.isArray(init)) {
    return new Headers(init);
  } else if (typeof init === "object") {
    return new Headers(init);
  } else {
    return new Headers();
  }
}
function mergeHeaders(...headers) {
  return headers.reduce((acc, header) => {
    const headersInstance = toHeadersInstance(header);
    for (const [key, value] of headersInstance.entries()) {
      if (key === "set-cookie") {
        const splitCookies = splitSetCookieString$1(value);
        splitCookies.forEach((cookie) => acc.append("set-cookie", cookie));
      } else {
        acc.set(key, value);
      }
    }
    return acc;
  }, new Headers());
}
const minifiedTsrBootStrapScript = 'self.$_TSR={c(){document.querySelectorAll(".\\\\$tsr").forEach(e=>{e.remove()}),this.hydrated&&this.streamEnd&&(delete self.$_TSR,delete self.$R.tsr)},p(e){this.initialized?e():this.buffer.push(e)},buffer:[]};\n';
function transformReadableStreamWithRouter(router, routerStream) {
  return transformStreamWithRouter(router, routerStream);
}
function transformPipeableStreamWithRouter(router, routerStream) {
  return Readable.fromWeb(
    transformStreamWithRouter(router, Readable.toWeb(routerStream))
  );
}
const TSR_SCRIPT_BARRIER_ID = "$tsr-stream-barrier";
const patternBodyEnd = /(<\/body>)/;
const patternHtmlEnd = /(<\/html>)/;
const patternClosingTag = /(<\/[a-zA-Z][\w:.-]*?>)/g;
function createPassthrough(onCancel) {
  let controller;
  const encoder = new TextEncoder();
  const stream = new ReadableStream$1({
    start(c2) {
      controller = c2;
    },
    cancel() {
      res.destroyed = true;
      onCancel();
    }
  });
  const res = {
    stream,
    write: (chunk) => {
      if (res.destroyed) return;
      if (typeof chunk === "string") {
        controller.enqueue(encoder.encode(chunk));
      } else {
        controller.enqueue(chunk);
      }
    },
    end: (chunk) => {
      if (res.destroyed) return;
      if (chunk) {
        res.write(chunk);
      }
      res.destroyed = true;
      controller.close();
    },
    destroy: (error) => {
      if (res.destroyed) return;
      res.destroyed = true;
      controller.error(error);
    },
    destroyed: false
  };
  return res;
}
async function readStream(stream, opts) {
  const reader = stream.getReader();
  try {
    let chunk;
    while (!(chunk = await reader.read()).done) {
      opts.onData?.(chunk);
    }
    opts.onEnd?.();
  } catch (error) {
    opts.onError?.(error);
  } finally {
    reader.releaseLock();
  }
}
function transformStreamWithRouter(router, appStream, opts) {
  let stopListeningToInjectedHtml = void 0;
  let timeoutHandle;
  let cleanedUp = false;
  function cleanup() {
    if (cleanedUp) return;
    cleanedUp = true;
    if (stopListeningToInjectedHtml) {
      stopListeningToInjectedHtml();
      stopListeningToInjectedHtml = void 0;
    }
    clearTimeout(timeoutHandle);
    router.serverSsr?.cleanup();
  }
  const finalPassThrough = createPassthrough(cleanup);
  const textDecoder = new TextDecoder();
  let isAppRendering = true;
  let routerStreamBuffer = "";
  let pendingClosingTags = "";
  let streamBarrierLifted = false;
  let leftover = "";
  let leftoverHtml = "";
  function getBufferedRouterStream() {
    const html = routerStreamBuffer;
    routerStreamBuffer = "";
    return html;
  }
  function decodeChunk(chunk) {
    if (chunk instanceof Uint8Array) {
      return textDecoder.decode(chunk, { stream: true });
    }
    return String(chunk);
  }
  const injectedHtmlDonePromise = createControlledPromise();
  let processingCount = 0;
  handleInjectedHtml();
  stopListeningToInjectedHtml = router.subscribe(
    "onInjectedHtml",
    handleInjectedHtml
  );
  function handleInjectedHtml() {
    if (cleanedUp) return;
    router.serverSsr.injectedHtml.forEach((promise) => {
      processingCount++;
      promise.then((html) => {
        if (cleanedUp || finalPassThrough.destroyed) {
          return;
        }
        if (isAppRendering) {
          routerStreamBuffer += html;
        } else {
          finalPassThrough.write(html);
        }
      }).catch((err) => {
        injectedHtmlDonePromise.reject(err);
      }).finally(() => {
        processingCount--;
        if (!isAppRendering && processingCount === 0) {
          injectedHtmlDonePromise.resolve();
        }
      });
    });
    router.serverSsr.injectedHtml = [];
  }
  injectedHtmlDonePromise.then(() => {
    if (cleanedUp || finalPassThrough.destroyed) {
      return;
    }
    clearTimeout(timeoutHandle);
    const finalHtml = leftover + leftoverHtml + getBufferedRouterStream() + pendingClosingTags;
    leftover = "";
    leftoverHtml = "";
    pendingClosingTags = "";
    finalPassThrough.end(finalHtml);
  }).catch((err) => {
    if (cleanedUp || finalPassThrough.destroyed) {
      return;
    }
    console.error("Error reading routerStream:", err);
    finalPassThrough.destroy(err);
  }).finally(cleanup);
  readStream(appStream, {
    onData: (chunk) => {
      if (cleanedUp || finalPassThrough.destroyed) {
        return;
      }
      const text = decodeChunk(chunk.value);
      const chunkString = leftover + text;
      const bodyEndMatch = chunkString.match(patternBodyEnd);
      const htmlEndMatch = chunkString.match(patternHtmlEnd);
      if (!streamBarrierLifted) {
        const streamBarrierIdIncluded = chunkString.includes(
          TSR_SCRIPT_BARRIER_ID
        );
        if (streamBarrierIdIncluded) {
          streamBarrierLifted = true;
          router.serverSsr.liftScriptBarrier();
        }
      }
      if (bodyEndMatch && htmlEndMatch && bodyEndMatch.index < htmlEndMatch.index) {
        const bodyEndIndex = bodyEndMatch.index;
        pendingClosingTags = chunkString.slice(bodyEndIndex);
        finalPassThrough.write(
          chunkString.slice(0, bodyEndIndex) + getBufferedRouterStream() + leftoverHtml
        );
        leftover = "";
        leftoverHtml = "";
        return;
      }
      let result;
      let lastIndex = 0;
      patternClosingTag.lastIndex = 0;
      while ((result = patternClosingTag.exec(chunkString)) !== null) {
        lastIndex = result.index + result[0].length;
      }
      if (lastIndex > 0) {
        const processed = chunkString.slice(0, lastIndex) + getBufferedRouterStream() + leftoverHtml;
        finalPassThrough.write(processed);
        leftover = chunkString.slice(lastIndex);
        leftoverHtml = "";
      } else {
        leftover = chunkString;
        leftoverHtml += getBufferedRouterStream();
      }
    },
    onEnd: () => {
      if (cleanedUp || finalPassThrough.destroyed) {
        return;
      }
      isAppRendering = false;
      router.serverSsr.setRenderFinished();
      if (processingCount === 0) {
        injectedHtmlDonePromise.resolve();
      } else {
        const timeoutMs = 6e4;
        timeoutHandle = setTimeout(() => {
          injectedHtmlDonePromise.reject(
            new Error("Injected HTML timeout after app render finished")
          );
        }, timeoutMs);
      }
    },
    onError: (error) => {
      if (cleanedUp) {
        return;
      }
      console.error("Error reading appStream:", error);
      isAppRendering = false;
      router.serverSsr.setRenderFinished();
      clearTimeout(timeoutHandle);
      leftover = "";
      leftoverHtml = "";
      routerStreamBuffer = "";
      pendingClosingTags = "";
      finalPassThrough.destroy(error);
      injectedHtmlDonePromise.reject(error);
    }
  });
  return finalPassThrough.stream;
}
const SCOPE_ID = "tsr";
function dehydrateMatch(match) {
  const dehydratedMatch = {
    i: match.id,
    u: match.updatedAt,
    s: match.status
  };
  const properties = [
    ["__beforeLoadContext", "b"],
    ["loaderData", "l"],
    ["error", "e"],
    ["ssr", "ssr"]
  ];
  for (const [key, shorthand] of properties) {
    if (match[key] !== void 0) {
      dehydratedMatch[shorthand] = match[key];
    }
  }
  return dehydratedMatch;
}
const INITIAL_SCRIPTS = [
  rn(SCOPE_ID),
  minifiedTsrBootStrapScript
];
class ScriptBuffer {
  constructor(router) {
    this._queue = [...INITIAL_SCRIPTS];
    this._scriptBarrierLifted = false;
    this._cleanedUp = false;
    this.router = router;
  }
  enqueue(script) {
    if (this._cleanedUp) return;
    if (this._scriptBarrierLifted && this._queue.length === 0) {
      queueMicrotask(() => {
        this.injectBufferedScripts();
      });
    }
    this._queue.push(script);
  }
  liftBarrier() {
    if (this._scriptBarrierLifted || this._cleanedUp) return;
    this._scriptBarrierLifted = true;
    if (this._queue.length > 0) {
      queueMicrotask(() => {
        this.injectBufferedScripts();
      });
    }
  }
  takeAll() {
    const bufferedScripts = this._queue;
    this._queue = [];
    if (bufferedScripts.length === 0) {
      return void 0;
    }
    bufferedScripts.push(`${GLOBAL_TSR}.c()`);
    const joinedScripts = bufferedScripts.join(";");
    return joinedScripts;
  }
  injectBufferedScripts() {
    if (this._cleanedUp) return;
    const scriptsToInject = this.takeAll();
    if (scriptsToInject && this.router?.serverSsr) {
      this.router.serverSsr.injectScript(() => scriptsToInject);
    }
  }
  cleanup() {
    this._cleanedUp = true;
    this._queue = [];
    this.router = void 0;
  }
}
function attachRouterServerSsrUtils({
  router,
  manifest: manifest2
}) {
  router.ssr = {
    manifest: manifest2
  };
  let _dehydrated = false;
  const listeners = [];
  const scriptBuffer = new ScriptBuffer(router);
  router.serverSsr = {
    injectedHtml: [],
    injectHtml: (getHtml) => {
      const promise = Promise.resolve().then(getHtml);
      router.serverSsr.injectedHtml.push(promise);
      router.emit({
        type: "onInjectedHtml",
        promise
      });
      return promise.then(() => {
      });
    },
    injectScript: (getScript) => {
      return router.serverSsr.injectHtml(async () => {
        const script = await getScript();
        if (!script) {
          return "";
        }
        return `<script${router.options.ssr?.nonce ? ` nonce='${router.options.ssr.nonce}'` : ""} class='$tsr'>${script}<\/script>`;
      });
    },
    dehydrate: async () => {
      invariant(!_dehydrated, "router is already dehydrated!");
      let matchesToDehydrate = router.state.matches;
      if (router.isShell()) {
        matchesToDehydrate = matchesToDehydrate.slice(0, 1);
      }
      const matches = matchesToDehydrate.map(dehydrateMatch);
      let manifestToDehydrate = void 0;
      if (manifest2) {
        const filteredRoutes = Object.fromEntries(
          router.state.matches.map((k2) => [
            k2.routeId,
            manifest2.routes[k2.routeId]
          ])
        );
        manifestToDehydrate = {
          routes: filteredRoutes
        };
      }
      const dehydratedRouter = {
        manifest: manifestToDehydrate,
        matches
      };
      const lastMatchId = matchesToDehydrate[matchesToDehydrate.length - 1]?.id;
      if (lastMatchId) {
        dehydratedRouter.lastMatchId = lastMatchId;
      }
      const dehydratedData = await router.options.dehydrate?.();
      if (dehydratedData) {
        dehydratedRouter.dehydratedData = dehydratedData;
      }
      _dehydrated = true;
      const p2 = createControlledPromise();
      const trackPlugins = { didRun: false };
      const plugins = router.options.serializationAdapters?.map((t) => makeSsrSerovalPlugin(t, trackPlugins)) ?? [];
      Xt(dehydratedRouter, {
        refs: /* @__PURE__ */ new Map(),
        plugins: [...plugins, ...defaultSerovalPlugins],
        onSerialize: (data, initial) => {
          let serialized = initial ? GLOBAL_TSR + ".router=" + data : data;
          if (trackPlugins.didRun) {
            serialized = GLOBAL_TSR + ".p(()=>" + serialized + ")";
          }
          scriptBuffer.enqueue(serialized);
        },
        scopeId: SCOPE_ID,
        onDone: () => {
          scriptBuffer.enqueue(GLOBAL_TSR + ".streamEnd=true");
          p2.resolve("");
        },
        onError: (err) => p2.reject(err)
      });
      router.serverSsr.injectHtml(() => p2);
    },
    isDehydrated() {
      return _dehydrated;
    },
    onRenderFinished: (listener) => listeners.push(listener),
    setRenderFinished: () => {
      listeners.forEach((l) => l());
      listeners.length = 0;
      scriptBuffer.liftBarrier();
    },
    takeBufferedScripts() {
      const scripts = scriptBuffer.takeAll();
      const serverBufferedScript = {
        tag: "script",
        attrs: {
          nonce: router.options.ssr?.nonce,
          className: "$tsr",
          id: TSR_SCRIPT_BARRIER_ID
        },
        children: scripts
      };
      return serverBufferedScript;
    },
    liftScriptBarrier() {
      scriptBuffer.liftBarrier();
    },
    cleanup() {
      if (!router.serverSsr) return;
      listeners.length = 0;
      scriptBuffer.cleanup();
      router.serverSsr.injectedHtml = [];
      router.serverSsr = void 0;
    }
  };
}
function getOrigin(request) {
  const originHeader = request.headers.get("Origin");
  if (originHeader) {
    try {
      new URL(originHeader);
      return originHeader;
    } catch {
    }
  }
  try {
    return new URL(request.url).origin;
  } catch {
  }
  return "http://localhost";
}
function defineHandlerCallback(handler) {
  return handler;
}
var fullPattern = " daum[ /]| deusu/|(?:^|[^g])news(?!sapphire)|(?<! (?:channel/|google/))google(?!(app|/google| pixel))|(?<! cu)bots?(?:\\b|_)|(?<!(?:lib))http|(?<![hg]m)score|(?<!cam)scan|24x7|@[a-z][\\w-]+\\.|\\(\\)|\\.com\\b|\\bperl\\b|\\btime/|\\||^[\\w \\.\\-\\(?:\\):%]+(?:/v?\\d+(?:\\.\\d+)?(?:\\.\\d{1,10})*?)?(?:,|$)|^[^ ]{50,}$|^\\d+\\b|^\\W|^\\w*search\\b|^\\w+/[\\w\\(\\)]*$|^active|^ad muncher|^amaya|^avsdevicesdk/|^azure|^biglotron|^bot|^bw/|^clamav[ /]|^client/|^cobweb/|^custom|^ddg[_-]android|^discourse|^dispatch/\\d|^downcast/|^duckduckgo|^email|^facebook|^getright/|^gozilla/|^hobbit|^hotzonu|^hwcdn/|^igetter/|^jeode/|^jetty/|^jigsaw|^microsoft bits|^movabletype|^mozilla/\\d\\.\\d\\s[\\w\\.-]+$|^mozilla/\\d\\.\\d\\s\\(compatible;?(?:\\s\\w+\\/\\d+\\.\\d+)?\\)$|^navermailapp|^netsurf|^offline|^openai/|^owler|^php|^postman|^python|^rank|^read|^reed|^rest|^rss|^snapchat|^space bison|^svn|^swcd |^taringa|^thumbor/|^track|^w3c|^webbandit/|^webcopier|^wget|^whatsapp|^wordpress|^xenu link sleuth|^yahoo|^yandex|^zdm/\\d|^zoom marketplace/|agent\\b|analyzer|archive|ask jeeves/teoma|audit|bit\\.ly/|bluecoat drtr|browsex|burpcollaborator|capture|catch|check\\b|checker|chrome-lighthouse|chromeframe|classifier|cloudflare|convertify|crawl|cypress/|dareboost|datanyze|dejaclick|detect|dmbrowser|download|evc-batch/|exaleadcloudview|feed|fetcher|firephp|functionize|grab|headless|httrack|hubspot marketing grader|hydra|ibisbrowser|infrawatch|insight|inspect|iplabel|java(?!;)|library|linkcheck|mail\\.ru/|manager|measure|neustar wpm|node\\b|nutch|offbyone|onetrust|optimize|pageburst|pagespeed|parser|phantomjs|pingdom|powermarks|preview|proxy|ptst[ /]\\d|retriever|rexx;|rigor|rss\\b|scrape|server|sogou|sparkler/|speedcurve|spider|splash|statuscake|supercleaner|synapse|synthetic|tools|torrent|transcoder|url|validator|virtuoso|wappalyzer|webglance|webkit2png|whatcms/|xtate/";
var naivePattern = /bot|crawl|http|lighthouse|scan|search|spider/i;
var pattern;
function getPattern() {
  if (pattern instanceof RegExp) {
    return pattern;
  }
  try {
    pattern = new RegExp(fullPattern, "i");
  } catch (error) {
    pattern = naivePattern;
  }
  return pattern;
}
function isbot(userAgent) {
  return Boolean(userAgent) && getPattern().test(userAgent);
}
const renderRouterToStream = async ({
  request,
  router,
  responseHeaders,
  children
}) => {
  if (typeof ReactDOMServer.renderToReadableStream === "function") {
    const stream = await ReactDOMServer.renderToReadableStream(children, {
      signal: request.signal,
      nonce: router.options.ssr?.nonce,
      progressiveChunkSize: Number.POSITIVE_INFINITY
    });
    if (isbot(request.headers.get("User-Agent"))) {
      await stream.allReady;
    }
    const responseStream = transformReadableStreamWithRouter(
      router,
      stream
    );
    return new Response(responseStream, {
      status: router.state.statusCode,
      headers: responseHeaders
    });
  }
  if (typeof ReactDOMServer.renderToPipeableStream === "function") {
    const reactAppPassthrough = new PassThrough();
    try {
      const pipeable = ReactDOMServer.renderToPipeableStream(children, {
        nonce: router.options.ssr?.nonce,
        progressiveChunkSize: Number.POSITIVE_INFINITY,
        ...isbot(request.headers.get("User-Agent")) ? {
          onAllReady() {
            pipeable.pipe(reactAppPassthrough);
          }
        } : {
          onShellReady() {
            pipeable.pipe(reactAppPassthrough);
          }
        },
        onError: (error, info) => {
          console.error("Error in renderToPipeableStream:", error, info);
          if (!reactAppPassthrough.destroyed) {
            reactAppPassthrough.destroy(
              error instanceof Error ? error : new Error(String(error))
            );
          }
        }
      });
    } catch (e) {
      console.error("Error in renderToPipeableStream:", e);
      reactAppPassthrough.destroy(e instanceof Error ? e : new Error(String(e)));
    }
    const responseStream = transformPipeableStreamWithRouter(
      router,
      reactAppPassthrough
    );
    return new Response(responseStream, {
      status: router.state.statusCode,
      headers: responseHeaders
    });
  }
  throw new Error(
    "No renderToReadableStream or renderToPipeableStream found in react-dom/server. Ensure you are using a version of react-dom that supports streaming."
  );
};
const defaultStreamHandler = defineHandlerCallback(
  ({ request, router, responseHeaders }) => renderRouterToStream({
    request,
    router,
    responseHeaders,
    children: /* @__PURE__ */ jsx(StartServer, { router })
  })
);
function json(payload, init) {
  return new Response(JSON.stringify(payload), {
    ...init,
    headers: mergeHeaders(
      { "content-type": "application/json" },
      init?.headers
    )
  });
}
const TSS_FORMDATA_CONTEXT = "__TSS_CONTEXT";
const TSS_SERVER_FUNCTION = Symbol.for("TSS_SERVER_FUNCTION");
const X_TSS_SERIALIZED = "x-tss-serialized";
const X_TSS_RAW_RESPONSE = "x-tss-raw";
const startStorage = new AsyncLocalStorage();
async function runWithStartContext(context, fn2) {
  return startStorage.run(context, fn2);
}
function getStartContext(opts) {
  const context = startStorage.getStore();
  if (!context && opts?.throwIfNotFound !== false) {
    throw new Error(
      `No Start context found in AsyncLocalStorage. Make sure you are using the function within the server runtime.`
    );
  }
  return context;
}
const getStartOptions = () => getStartContext().startOptions;
function flattenMiddlewares(middlewares) {
  const seen = /* @__PURE__ */ new Set();
  const flattened = [];
  const recurse = (middleware) => {
    middleware.forEach((m2) => {
      if (m2.options.middleware) {
        recurse(m2.options.middleware);
      }
      if (!seen.has(m2)) {
        seen.add(m2);
        flattened.push(m2);
      }
    });
  };
  recurse(middlewares);
  return flattened;
}
function getDefaultSerovalPlugins() {
  const start = getStartOptions();
  const adapters = start?.serializationAdapters;
  return [
    ...adapters?.map(makeSerovalPlugin) ?? [],
    ...defaultSerovalPlugins
  ];
}
const NullProtoObj = /* @__PURE__ */ (() => {
  const e = function() {
  };
  return e.prototype = /* @__PURE__ */ Object.create(null), Object.freeze(e.prototype), e;
})();
function splitSetCookieString(cookiesString) {
  if (Array.isArray(cookiesString)) return cookiesString.flatMap((c2) => splitSetCookieString(c2));
  if (typeof cookiesString !== "string") return [];
  const cookiesStrings = [];
  let pos = 0;
  let start;
  let ch;
  let lastComma;
  let nextStart;
  let cookiesSeparatorFound;
  const skipWhitespace = () => {
    while (pos < cookiesString.length && /\s/.test(cookiesString.charAt(pos))) pos += 1;
    return pos < cookiesString.length;
  };
  const notSpecialChar = () => {
    ch = cookiesString.charAt(pos);
    return ch !== "=" && ch !== ";" && ch !== ",";
  };
  while (pos < cookiesString.length) {
    start = pos;
    cookiesSeparatorFound = false;
    while (skipWhitespace()) {
      ch = cookiesString.charAt(pos);
      if (ch === ",") {
        lastComma = pos;
        pos += 1;
        skipWhitespace();
        nextStart = pos;
        while (pos < cookiesString.length && notSpecialChar()) pos += 1;
        if (pos < cookiesString.length && cookiesString.charAt(pos) === "=") {
          cookiesSeparatorFound = true;
          pos = nextStart;
          cookiesStrings.push(cookiesString.slice(start, lastComma));
          start = pos;
        } else pos = lastComma + 1;
      } else pos += 1;
    }
    if (!cookiesSeparatorFound || pos >= cookiesString.length) cookiesStrings.push(cookiesString.slice(start));
  }
  return cookiesStrings;
}
function lazyInherit(target, source, sourceKey) {
  for (const key of Object.getOwnPropertyNames(source)) {
    if (key === "constructor") continue;
    const targetDesc = Object.getOwnPropertyDescriptor(target, key);
    const desc = Object.getOwnPropertyDescriptor(source, key);
    let modified = false;
    if (desc.get) {
      modified = true;
      desc.get = targetDesc?.get || function() {
        return this[sourceKey][key];
      };
    }
    if (desc.set) {
      modified = true;
      desc.set = targetDesc?.set || function(value) {
        this[sourceKey][key] = value;
      };
    }
    if (typeof desc.value === "function") {
      modified = true;
      desc.value = function(...args) {
        return this[sourceKey][key](...args);
      };
    }
    if (modified) Object.defineProperty(target, key, desc);
  }
}
const FastURL = /* @__PURE__ */ (() => {
  const NativeURL = globalThis.URL;
  const FastURL$1 = class URL {
    #url;
    #href;
    #protocol;
    #host;
    #pathname;
    #search;
    #searchParams;
    #pos;
    constructor(url) {
      if (typeof url === "string") this.#href = url;
      else {
        this.#protocol = url.protocol;
        this.#host = url.host;
        this.#pathname = url.pathname;
        this.#search = url.search;
      }
    }
    get _url() {
      if (this.#url) return this.#url;
      this.#url = new NativeURL(this.href);
      this.#href = void 0;
      this.#protocol = void 0;
      this.#host = void 0;
      this.#pathname = void 0;
      this.#search = void 0;
      this.#searchParams = void 0;
      this.#pos = void 0;
      return this.#url;
    }
    get href() {
      if (this.#url) return this.#url.href;
      if (!this.#href) this.#href = `${this.#protocol || "http:"}//${this.#host || "localhost"}${this.#pathname || "/"}${this.#search || ""}`;
      return this.#href;
    }
    #getPos() {
      if (!this.#pos) {
        const url = this.href;
        const protoIndex = url.indexOf("://");
        const pathnameIndex = protoIndex === -1 ? -1 : url.indexOf("/", protoIndex + 4);
        const qIndex = pathnameIndex === -1 ? -1 : url.indexOf("?", pathnameIndex);
        this.#pos = [
          protoIndex,
          pathnameIndex,
          qIndex
        ];
      }
      return this.#pos;
    }
    get pathname() {
      if (this.#url) return this.#url.pathname;
      if (this.#pathname === void 0) {
        const [, pathnameIndex, queryIndex] = this.#getPos();
        if (pathnameIndex === -1) return this._url.pathname;
        this.#pathname = this.href.slice(pathnameIndex, queryIndex === -1 ? void 0 : queryIndex);
      }
      return this.#pathname;
    }
    get search() {
      if (this.#url) return this.#url.search;
      if (this.#search === void 0) {
        const [, pathnameIndex, queryIndex] = this.#getPos();
        if (pathnameIndex === -1) return this._url.search;
        const url = this.href;
        this.#search = queryIndex === -1 || queryIndex === url.length - 1 ? "" : url.slice(queryIndex);
      }
      return this.#search;
    }
    get searchParams() {
      if (this.#url) return this.#url.searchParams;
      if (!this.#searchParams) this.#searchParams = new URLSearchParams(this.search);
      return this.#searchParams;
    }
    get protocol() {
      if (this.#url) return this.#url.protocol;
      if (this.#protocol === void 0) {
        const [protocolIndex] = this.#getPos();
        if (protocolIndex === -1) return this._url.protocol;
        const url = this.href;
        this.#protocol = url.slice(0, protocolIndex + 1);
      }
      return this.#protocol;
    }
    toString() {
      return this.href;
    }
    toJSON() {
      return this.href;
    }
  };
  lazyInherit(FastURL$1.prototype, NativeURL.prototype, "_url");
  Object.setPrototypeOf(FastURL$1.prototype, NativeURL.prototype);
  Object.setPrototypeOf(FastURL$1, NativeURL);
  return FastURL$1;
})();
const NodeResponse = /* @__PURE__ */ (() => {
  const NativeResponse = globalThis.Response;
  const STATUS_CODES = globalThis.process?.getBuiltinModule?.("node:http")?.STATUS_CODES || {};
  class NodeResponse$1 {
    #body;
    #init;
    #headers;
    #response;
    constructor(body, init) {
      this.#body = body;
      this.#init = init;
    }
    get status() {
      return this.#response?.status || this.#init?.status || 200;
    }
    get statusText() {
      return this.#response?.statusText || this.#init?.statusText || STATUS_CODES[this.status] || "";
    }
    get headers() {
      if (this.#response) return this.#response.headers;
      if (this.#headers) return this.#headers;
      const initHeaders = this.#init?.headers;
      return this.#headers = initHeaders instanceof Headers ? initHeaders : new Headers(initHeaders);
    }
    get ok() {
      if (this.#response) return this.#response.ok;
      const status = this.status;
      return status >= 200 && status < 300;
    }
    get _response() {
      if (this.#response) return this.#response;
      this.#response = new NativeResponse(this.#body, this.#headers ? {
        ...this.#init,
        headers: this.#headers
      } : this.#init);
      this.#init = void 0;
      this.#headers = void 0;
      this.#body = void 0;
      return this.#response;
    }
    nodeResponse() {
      const status = this.status;
      const statusText = this.statusText;
      let body;
      let contentType;
      let contentLength;
      if (this.#response) body = this.#response.body;
      else if (this.#body) if (this.#body instanceof ReadableStream) body = this.#body;
      else if (typeof this.#body === "string") {
        body = this.#body;
        contentType = "text/plain; charset=UTF-8";
        contentLength = Buffer.byteLength(this.#body);
      } else if (this.#body instanceof ArrayBuffer) {
        body = Buffer.from(this.#body);
        contentLength = this.#body.byteLength;
      } else if (this.#body instanceof Uint8Array) {
        body = this.#body;
        contentLength = this.#body.byteLength;
      } else if (this.#body instanceof DataView) {
        body = Buffer.from(this.#body.buffer);
        contentLength = this.#body.byteLength;
      } else if (this.#body instanceof Blob) {
        body = this.#body.stream();
        contentType = this.#body.type;
        contentLength = this.#body.size;
      } else if (typeof this.#body.pipe === "function") body = this.#body;
      else body = this._response.body;
      const rawNodeHeaders = [];
      const initHeaders = this.#init?.headers;
      const headerEntries = this.#response?.headers || this.#headers || (initHeaders ? Array.isArray(initHeaders) ? initHeaders : initHeaders?.entries ? initHeaders.entries() : Object.entries(initHeaders).map(([k2, v2]) => [k2.toLowerCase(), v2]) : void 0);
      let hasContentTypeHeader;
      let hasContentLength;
      if (headerEntries) for (const [key, value] of headerEntries) {
        if (key === "set-cookie") {
          for (const setCookie of splitSetCookieString(value)) rawNodeHeaders.push(["set-cookie", setCookie]);
          continue;
        }
        rawNodeHeaders.push([key, value]);
        if (key === "content-type") hasContentTypeHeader = true;
        else if (key === "content-length") hasContentLength = true;
      }
      if (contentType && !hasContentTypeHeader) rawNodeHeaders.push(["content-type", contentType]);
      if (contentLength && !hasContentLength) rawNodeHeaders.push(["content-length", String(contentLength)]);
      this.#init = void 0;
      this.#headers = void 0;
      this.#response = void 0;
      this.#body = void 0;
      return {
        status,
        statusText,
        headers: rawNodeHeaders,
        body
      };
    }
  }
  lazyInherit(NodeResponse$1.prototype, NativeResponse.prototype, "_response");
  Object.setPrototypeOf(NodeResponse$1, NativeResponse);
  Object.setPrototypeOf(NodeResponse$1.prototype, NativeResponse.prototype);
  return NodeResponse$1;
})();
var H3Event = class {
  /**
  * Access to the H3 application instance.
  */
  app;
  /**
  * Incoming HTTP request info.
  *
  * [MDN Reference](https://developer.mozilla.org/en-US/docs/Web/API/Request)
  */
  req;
  /**
  * Access to the parsed request URL.
  *
  * [MDN Reference](https://developer.mozilla.org/en-US/docs/Web/API/URL)
  */
  url;
  /**
  * Event context.
  */
  context;
  /**
  * @internal
  */
  static __is_event__ = true;
  /**
  * @internal
  */
  _res;
  constructor(req, context, app) {
    this.context = context || req.context || new NullProtoObj();
    this.req = req;
    this.app = app;
    const _url = req._url;
    this.url = _url && _url instanceof URL ? _url : new FastURL(req.url);
  }
  /**
  * Prepared HTTP response.
  */
  get res() {
    if (!this._res) this._res = new H3EventResponse();
    return this._res;
  }
  /**
  * Access to runtime specific additional context.
  *
  */
  get runtime() {
    return this.req.runtime;
  }
  /**
  * Tell the runtime about an ongoing operation that shouldn't close until the promise resolves.
  */
  waitUntil(promise) {
    this.req.waitUntil?.(promise);
  }
  toString() {
    return `[${this.req.method}] ${this.req.url}`;
  }
  toJSON() {
    return this.toString();
  }
  /**
  * Access to the raw Node.js req/res objects.
  *
  * @deprecated Use `event.runtime.{node|deno|bun|...}.` instead.
  */
  get node() {
    return this.req.runtime?.node;
  }
  /**
  * Access to the incoming request headers.
  *
  * @deprecated Use `event.req.headers` instead.
  *
  */
  get headers() {
    return this.req.headers;
  }
  /**
  * Access to the incoming request url (pathname+search).
  *
  * @deprecated Use `event.url.pathname + event.url.search` instead.
  *
  * Example: `/api/hello?name=world`
  * */
  get path() {
    return this.url.pathname + this.url.search;
  }
  /**
  * Access to the incoming request method.
  *
  * @deprecated Use `event.req.method` instead.
  */
  get method() {
    return this.req.method;
  }
};
var H3EventResponse = class {
  status;
  statusText;
  _headers;
  get headers() {
    if (!this._headers) this._headers = new Headers();
    return this._headers;
  }
};
const DISALLOWED_STATUS_CHARS = /[^\u0009\u0020-\u007E]/g;
function sanitizeStatusMessage(statusMessage = "") {
  return statusMessage.replace(DISALLOWED_STATUS_CHARS, "");
}
function sanitizeStatusCode(statusCode, defaultStatusCode = 200) {
  if (!statusCode) return defaultStatusCode;
  if (typeof statusCode === "string") statusCode = +statusCode;
  if (statusCode < 100 || statusCode > 599) return defaultStatusCode;
  return statusCode;
}
var HTTPError = class HTTPError2 extends Error {
  get name() {
    return "HTTPError";
  }
  /**
  * HTTP status code in range [200...599]
  */
  status;
  /**
  * HTTP status text
  *
  * **NOTE:** This should be short (max 512 to 1024 characters).
  * Allowed characters are tabs, spaces, visible ASCII characters, and extended characters (byte value 128–255).
  *
  * **TIP:** Use `message` for longer error descriptions in JSON body.
  */
  statusText;
  /**
  * Additional HTTP headers to be sent in error response.
  */
  headers;
  /**
  * Original error object that caused this error.
  */
  cause;
  /**
  * Additional data attached in the error JSON body under `data` key.
  */
  data;
  /**
  * Additional top level JSON body properties to attach in the error JSON body.
  */
  body;
  /**
  * Flag to indicate that the error was not handled by the application.
  *
  * Unhandled error stack trace, data and message are hidden in non debug mode for security reasons.
  */
  unhandled;
  /**
  * Check if the input is an instance of HTTPError using its constructor name.
  *
  * It is safer than using `instanceof` because it works across different contexts (e.g., if the error was thrown in a different module).
  */
  static isError(input) {
    return input instanceof Error && input?.name === "HTTPError";
  }
  /**
  * Create a new HTTPError with the given status code and optional status text and details.
  *
  * @example
  *
  * HTTPError.status(404)
  * HTTPError.status(418, "I'm a teapot")
  * HTTPError.status(403, "Forbidden", { message: "Not authenticated" })
  */
  static status(status, statusText, details) {
    return new HTTPError2({
      ...details,
      statusText,
      status
    });
  }
  constructor(arg1, arg2) {
    let messageInput;
    let details;
    if (typeof arg1 === "string") {
      messageInput = arg1;
      details = arg2;
    } else details = arg1;
    const status = sanitizeStatusCode(details?.status || details?.cause?.status || details?.status || details?.statusCode, 500);
    const statusText = sanitizeStatusMessage(details?.statusText || details?.cause?.statusText || details?.statusText || details?.statusMessage);
    const message = messageInput || details?.message || details?.cause?.message || details?.statusText || details?.statusMessage || [
      "HTTPError",
      status,
      statusText
    ].filter(Boolean).join(" ");
    super(message, { cause: details });
    this.cause = details;
    Error.captureStackTrace?.(this, this.constructor);
    this.status = status;
    this.statusText = statusText || void 0;
    const rawHeaders = details?.headers || details?.cause?.headers;
    this.headers = rawHeaders ? new Headers(rawHeaders) : void 0;
    this.unhandled = details?.unhandled ?? details?.cause?.unhandled ?? void 0;
    this.data = details?.data;
    this.body = details?.body;
  }
  /**
  * @deprecated Use `status`
  */
  get statusCode() {
    return this.status;
  }
  /**
  * @deprecated Use `statusText`
  */
  get statusMessage() {
    return this.statusText;
  }
  toJSON() {
    const unhandled = this.unhandled;
    return {
      status: this.status,
      statusText: this.statusText,
      unhandled,
      message: unhandled ? "HTTPError" : this.message,
      data: unhandled ? void 0 : this.data,
      ...unhandled ? void 0 : this.body
    };
  }
};
function isJSONSerializable(value, _type) {
  if (value === null || value === void 0) return true;
  if (_type !== "object") return _type === "boolean" || _type === "number" || _type === "string";
  if (typeof value.toJSON === "function") return true;
  if (Array.isArray(value)) return true;
  if (typeof value.pipe === "function" || typeof value.pipeTo === "function") return false;
  if (value instanceof NullProtoObj) return true;
  const proto = Object.getPrototypeOf(value);
  return proto === Object.prototype || proto === null;
}
const kNotFound = /* @__PURE__ */ Symbol.for("h3.notFound");
const kHandled = /* @__PURE__ */ Symbol.for("h3.handled");
function toResponse(val, event, config = {}) {
  if (typeof val?.then === "function") return (val.catch?.((error) => error) || Promise.resolve(val)).then((resolvedVal) => toResponse(resolvedVal, event, config));
  const response = prepareResponse(val, event, config);
  if (typeof response?.then === "function") return toResponse(response, event, config);
  const { onResponse: onResponse$1 } = config;
  return onResponse$1 ? Promise.resolve(onResponse$1(response, event)).then(() => response) : response;
}
function prepareResponse(val, event, config, nested) {
  if (val === kHandled) return new NodeResponse(null);
  if (val === kNotFound) val = new HTTPError({
    status: 404,
    message: `Cannot find any route matching [${event.req.method}] ${event.url}`
  });
  if (val && val instanceof Error) {
    const isHTTPError = HTTPError.isError(val);
    const error = isHTTPError ? val : new HTTPError(val);
    if (!isHTTPError) {
      error.unhandled = true;
      if (val?.stack) error.stack = val.stack;
    }
    if (error.unhandled && !config.silent) console.error(error);
    const { onError: onError$1 } = config;
    return onError$1 && !nested ? Promise.resolve(onError$1(error, event)).catch((error$1) => error$1).then((newVal) => prepareResponse(newVal ?? val, event, config, true)) : errorResponse(error, config.debug);
  }
  const eventHeaders = event.res._headers;
  if (!(val instanceof Response)) {
    const res = prepareResponseBody(val, event, config);
    const status = event.res.status;
    return new NodeResponse(nullBody(event.req.method, status) ? null : res.body, {
      status,
      statusText: event.res.statusText,
      headers: res.headers && eventHeaders ? mergeHeaders$1(res.headers, eventHeaders) : res.headers || eventHeaders
    });
  }
  if (!eventHeaders) return val;
  return new NodeResponse(nullBody(event.req.method, val.status) ? null : val.body, {
    status: val.status,
    statusText: val.statusText,
    headers: mergeHeaders$1(eventHeaders, val.headers)
  });
}
function mergeHeaders$1(base, merge) {
  const mergedHeaders = new Headers(base);
  for (const [name, value] of merge) if (name === "set-cookie") mergedHeaders.append(name, value);
  else mergedHeaders.set(name, value);
  return mergedHeaders;
}
const emptyHeaders = /* @__PURE__ */ new Headers({ "content-length": "0" });
const jsonHeaders = /* @__PURE__ */ new Headers({ "content-type": "application/json;charset=UTF-8" });
function prepareResponseBody(val, event, config) {
  if (val === null || val === void 0) return {
    body: "",
    headers: emptyHeaders
  };
  const valType = typeof val;
  if (valType === "string") return { body: val };
  if (val instanceof Uint8Array) {
    event.res.headers.set("content-length", val.byteLength.toString());
    return { body: val };
  }
  if (isJSONSerializable(val, valType)) return {
    body: JSON.stringify(val, void 0, config.debug ? 2 : void 0),
    headers: jsonHeaders
  };
  if (valType === "bigint") return {
    body: val.toString(),
    headers: jsonHeaders
  };
  if (val instanceof Blob) {
    const headers = {
      "content-type": val.type,
      "content-length": val.size.toString()
    };
    let filename = val.name;
    if (filename) {
      filename = encodeURIComponent(filename);
      headers["content-disposition"] = `filename="${filename}"; filename*=UTF-8''${filename}`;
    }
    return {
      body: val.stream(),
      headers
    };
  }
  if (valType === "symbol") return { body: val.toString() };
  if (valType === "function") return { body: `${val.name}()` };
  return { body: val };
}
function nullBody(method, status) {
  return method === "HEAD" || status === 100 || status === 101 || status === 102 || status === 204 || status === 205 || status === 304;
}
function errorResponse(error, debug) {
  return new NodeResponse(JSON.stringify({
    ...error.toJSON(),
    stack: debug && error.stack ? error.stack.split("\n").map((l) => l.trim()) : void 0
  }, void 0, debug ? 2 : void 0), {
    status: error.status,
    statusText: error.statusText,
    headers: error.headers ? mergeHeaders$1(jsonHeaders, error.headers) : jsonHeaders
  });
}
const eventStorage = new AsyncLocalStorage();
function requestHandler(handler) {
  return (request, requestOpts) => {
    const h3Event = new H3Event(request);
    const response = eventStorage.run(
      { h3Event },
      () => handler(request, requestOpts)
    );
    return toResponse(response, h3Event);
  };
}
function getH3Event() {
  const event = eventStorage.getStore();
  if (!event) {
    throw new Error(
      `No StartEvent found in AsyncLocalStorage. Make sure you are using the function within the server runtime.`
    );
  }
  return event.h3Event;
}
function getResponse() {
  const event = getH3Event();
  return event._res;
}
async function getStartManifest() {
  const { tsrStartManifest } = await import("./assets/_tanstack-start-manifest_v-uxLmYTyG.js");
  const startManifest = tsrStartManifest();
  const rootRoute = startManifest.routes[rootRouteId] = startManifest.routes[rootRouteId] || {};
  rootRoute.assets = rootRoute.assets || [];
  let script = `import('${startManifest.clientEntry}')`;
  rootRoute.assets.push({
    tag: "script",
    attrs: {
      type: "module",
      async: true
    },
    children: script
  });
  const manifest2 = {
    routes: Object.fromEntries(
      Object.entries(startManifest.routes).map(([k2, v2]) => {
        const result = {};
        let hasData = false;
        if (v2.preloads && v2.preloads.length > 0) {
          result["preloads"] = v2.preloads;
          hasData = true;
        }
        if (v2.assets && v2.assets.length > 0) {
          result["assets"] = v2.assets;
          hasData = true;
        }
        if (!hasData) {
          return [];
        }
        return [k2, result];
      })
    )
  };
  return manifest2;
}
const manifest = {};
async function getServerFnById(id) {
  const serverFnInfo = manifest[id];
  if (!serverFnInfo) {
    throw new Error("Server function info not found for " + id);
  }
  const fnModule = await serverFnInfo.importer();
  if (!fnModule) {
    console.info("serverFnInfo", serverFnInfo);
    throw new Error("Server function module not resolved for " + id);
  }
  const action = fnModule[serverFnInfo.functionName];
  if (!action) {
    console.info("serverFnInfo", serverFnInfo);
    console.info("fnModule", fnModule);
    throw new Error(
      `Server function module export not resolved for serverFn ID: ${id}`
    );
  }
  return action;
}
let regex = void 0;
const handleServerAction = async ({
  request,
  context
}) => {
  const controller = new AbortController();
  const signal = controller.signal;
  const abort = () => controller.abort();
  request.signal.addEventListener("abort", abort);
  if (regex === void 0) {
    regex = new RegExp(`${"/_serverFn/"}([^/?#]+)`);
  }
  const method = request.method;
  const url = new URL(request.url, "http://localhost:3000");
  const match = url.pathname.match(regex);
  const serverFnId = match ? match[1] : null;
  const search = Object.fromEntries(url.searchParams.entries());
  const isCreateServerFn = "createServerFn" in search;
  if (typeof serverFnId !== "string") {
    throw new Error("Invalid server action param for serverFnId: " + serverFnId);
  }
  const action = await getServerFnById(serverFnId);
  const formDataContentTypes = [
    "multipart/form-data",
    "application/x-www-form-urlencoded"
  ];
  const contentType = request.headers.get("Content-Type");
  const serovalPlugins = getDefaultSerovalPlugins();
  function parsePayload(payload) {
    const parsedPayload = Xi(payload, { plugins: serovalPlugins });
    return parsedPayload;
  }
  const response = await (async () => {
    try {
      let result = await (async () => {
        if (formDataContentTypes.some(
          (type) => contentType && contentType.includes(type)
        )) {
          invariant(
            method.toLowerCase() !== "get",
            "GET requests with FormData payloads are not supported"
          );
          const formData = await request.formData();
          const serializedContext = formData.get(TSS_FORMDATA_CONTEXT);
          formData.delete(TSS_FORMDATA_CONTEXT);
          const params = {
            context,
            data: formData
          };
          if (typeof serializedContext === "string") {
            try {
              const parsedContext = JSON.parse(serializedContext);
              const deserializedContext = Xi(parsedContext, {
                plugins: serovalPlugins
              });
              if (typeof deserializedContext === "object" && deserializedContext) {
                params.context = { ...context, ...deserializedContext };
              }
            } catch {
            }
          }
          return await action(params, signal);
        }
        if (method.toLowerCase() === "get") {
          invariant(
            isCreateServerFn,
            "expected GET request to originate from createServerFn"
          );
          let payload = search.payload;
          payload = payload ? parsePayload(JSON.parse(payload)) : {};
          payload.context = { ...context, ...payload.context };
          return await action(payload, signal);
        }
        if (method.toLowerCase() !== "post") {
          throw new Error("expected POST method");
        }
        let jsonPayload;
        if (contentType?.includes("application/json")) {
          jsonPayload = await request.json();
        }
        if (isCreateServerFn) {
          const payload = jsonPayload ? parsePayload(jsonPayload) : {};
          payload.context = { ...payload.context, ...context };
          return await action(payload, signal);
        }
        return await action(...jsonPayload);
      })();
      if (result.result instanceof Response) {
        result.result.headers.set(X_TSS_RAW_RESPONSE, "true");
        return result.result;
      }
      if (!isCreateServerFn) {
        result = result.result;
        if (result instanceof Response) {
          return result;
        }
      }
      if (isNotFound(result)) {
        return isNotFoundResponse(result);
      }
      const response2 = getResponse();
      let nonStreamingBody = void 0;
      if (result !== void 0) {
        let done = false;
        const callbacks = {
          onParse: (value) => {
            nonStreamingBody = value;
          },
          onDone: () => {
            done = true;
          },
          onError: (error) => {
            throw error;
          }
        };
        _i(result, {
          refs: /* @__PURE__ */ new Map(),
          plugins: serovalPlugins,
          onParse(value) {
            callbacks.onParse(value);
          },
          onDone() {
            callbacks.onDone();
          },
          onError: (error) => {
            callbacks.onError(error);
          }
        });
        if (done) {
          return new Response(
            nonStreamingBody ? JSON.stringify(nonStreamingBody) : void 0,
            {
              status: response2?.status,
              statusText: response2?.statusText,
              headers: {
                "Content-Type": "application/json",
                [X_TSS_SERIALIZED]: "true"
              }
            }
          );
        }
        const encoder = new TextEncoder();
        const stream = new ReadableStream({
          start(controller2) {
            callbacks.onParse = (value) => controller2.enqueue(encoder.encode(JSON.stringify(value) + "\n"));
            callbacks.onDone = () => {
              try {
                controller2.close();
              } catch (error) {
                controller2.error(error);
              }
            };
            callbacks.onError = (error) => controller2.error(error);
            if (nonStreamingBody !== void 0) {
              callbacks.onParse(nonStreamingBody);
            }
          }
        });
        return new Response(stream, {
          status: response2?.status,
          statusText: response2?.statusText,
          headers: {
            "Content-Type": "application/x-ndjson",
            [X_TSS_SERIALIZED]: "true"
          }
        });
      }
      return new Response(void 0, {
        status: response2?.status,
        statusText: response2?.statusText
      });
    } catch (error) {
      if (error instanceof Response) {
        return error;
      }
      if (isNotFound(error)) {
        return isNotFoundResponse(error);
      }
      console.info();
      console.info("Server Fn Error!");
      console.info();
      console.error(error);
      console.info();
      const serializedError = JSON.stringify(
        await Promise.resolve(
          ki(error, {
            refs: /* @__PURE__ */ new Map(),
            plugins: serovalPlugins
          })
        )
      );
      const response2 = getResponse();
      return new Response(serializedError, {
        status: response2?.status ?? 500,
        statusText: response2?.statusText,
        headers: {
          "Content-Type": "application/json",
          [X_TSS_SERIALIZED]: "true"
        }
      });
    }
  })();
  request.signal.removeEventListener("abort", abort);
  return response;
};
function isNotFoundResponse(error) {
  const { headers, ...rest } = error;
  return new Response(JSON.stringify(rest), {
    status: 404,
    headers: {
      "Content-Type": "application/json",
      ...headers || {}
    }
  });
}
const HEADERS = {
  TSS_SHELL: "X-TSS_SHELL"
};
const createServerRpc = (functionId, splitImportFn) => {
  return Object.assign(splitImportFn, {
    functionId,
    [TSS_SERVER_FUNCTION]: true
  });
};
const ServerFunctionSerializationAdapter = createSerializationAdapter({
  key: "$TSS/serverfn",
  test: (v2) => {
    if (typeof v2 !== "function") return false;
    if (!(TSS_SERVER_FUNCTION in v2)) return false;
    return !!v2[TSS_SERVER_FUNCTION];
  },
  toSerializable: ({ functionId }) => ({ functionId }),
  fromSerializable: ({ functionId }) => {
    const fn2 = async (opts, signal) => {
      const serverFn = await getServerFnById(functionId);
      const result = await serverFn(opts ?? {}, signal);
      return result.result;
    };
    return createServerRpc(functionId, fn2);
  }
});
function getStartResponseHeaders(opts) {
  const headers = mergeHeaders(
    {
      "Content-Type": "text/html; charset=utf-8"
    },
    ...opts.router.state.matches.map((match) => {
      return match.headers;
    })
  );
  return headers;
}
function createStartHandler(cb) {
  const ROUTER_BASEPATH = "/";
  let startRoutesManifest = null;
  let startEntry = null;
  let routerEntry = null;
  const getEntries = async () => {
    if (routerEntry === null) {
      routerEntry = await import("./assets/router-C1SSjqoG.js").then((n) => n.r);
    }
    if (startEntry === null) {
      startEntry = await import("./assets/start-HYkvq4Ni.js");
    }
    return {
      startEntry,
      routerEntry
    };
  };
  const startRequestResolver = async (request, requestOpts) => {
    let router = null;
    let cbWillCleanup = false;
    try {
      const origin = getOrigin(request);
      const url = new URL(request.url);
      const href = url.href.replace(url.origin, "");
      const startOptions = await (await getEntries()).startEntry.startInstance?.getOptions() || {};
      const serializationAdapters = [
        ...startOptions.serializationAdapters || [],
        ServerFunctionSerializationAdapter
      ];
      const requestStartOptions = {
        ...startOptions,
        serializationAdapters
      };
      const getRouter = async () => {
        if (router) return router;
        router = await (await getEntries()).routerEntry.getRouter();
        const isPrerendering = process.env.TSS_PRERENDERING === "true";
        let isShell = process.env.TSS_SHELL === "true";
        if (isPrerendering && !isShell) {
          isShell = request.headers.get(HEADERS.TSS_SHELL) === "true";
        }
        const history = createMemoryHistory({
          initialEntries: [href]
        });
        router.update({
          history,
          isShell,
          isPrerendering,
          origin: router.options.origin ?? origin,
          ...{
            defaultSsr: requestStartOptions.defaultSsr,
            serializationAdapters: [
              ...requestStartOptions.serializationAdapters,
              ...router.options.serializationAdapters || []
            ]
          },
          basepath: ROUTER_BASEPATH
        });
        return router;
      };
      const requestHandlerMiddleware = handlerToMiddleware(
        async ({ context }) => {
          const response2 = await runWithStartContext(
            {
              getRouter,
              startOptions: requestStartOptions,
              contextAfterGlobalMiddlewares: context,
              request
            },
            async () => {
              try {
                if (href.startsWith("/_serverFn/")) {
                  return await handleServerAction({
                    request,
                    context: requestOpts?.context
                  });
                }
                const executeRouter = async ({
                  serverContext
                }) => {
                  const requestAcceptHeader = request.headers.get("Accept") || "*/*";
                  const splitRequestAcceptHeader = requestAcceptHeader.split(",");
                  const supportedMimeTypes = ["*/*", "text/html"];
                  const isRouterAcceptSupported = supportedMimeTypes.some(
                    (mimeType) => splitRequestAcceptHeader.some(
                      (acceptedMimeType) => acceptedMimeType.trim().startsWith(mimeType)
                    )
                  );
                  if (!isRouterAcceptSupported) {
                    return json(
                      {
                        error: "Only HTML requests are supported here"
                      },
                      {
                        status: 500
                      }
                    );
                  }
                  if (startRoutesManifest === null) {
                    startRoutesManifest = await getStartManifest();
                  }
                  const router2 = await getRouter();
                  attachRouterServerSsrUtils({
                    router: router2,
                    manifest: startRoutesManifest
                  });
                  router2.update({ additionalContext: { serverContext } });
                  await router2.load();
                  if (router2.state.redirect) {
                    return router2.state.redirect;
                  }
                  await router2.serverSsr.dehydrate();
                  const responseHeaders = getStartResponseHeaders({ router: router2 });
                  cbWillCleanup = true;
                  const response4 = await cb({
                    request,
                    router: router2,
                    responseHeaders
                  });
                  return response4;
                };
                const response3 = await handleServerRoutes({
                  getRouter,
                  request,
                  executeRouter,
                  context
                });
                return response3;
              } catch (err) {
                if (err instanceof Response) {
                  return err;
                }
                throw err;
              }
            }
          );
          return response2;
        }
      );
      const flattenedMiddlewares = startOptions.requestMiddleware ? flattenMiddlewares(startOptions.requestMiddleware) : [];
      const middlewares = flattenedMiddlewares.map((d) => d.options.server);
      const ctx = await executeMiddleware(
        [...middlewares, requestHandlerMiddleware],
        {
          request,
          context: requestOpts?.context || {}
        }
      );
      const response = ctx.response;
      if (isRedirect(response)) {
        if (isResolvedRedirect(response)) {
          if (request.headers.get("x-tsr-redirect") === "manual") {
            return json(
              {
                ...response.options,
                isSerializedRedirect: true
              },
              {
                headers: response.headers
              }
            );
          }
          return response;
        }
        if (response.options.to && typeof response.options.to === "string" && !response.options.to.startsWith("/")) {
          throw new Error(
            `Server side redirects must use absolute paths via the 'href' or 'to' options. The redirect() method's "to" property accepts an internal path only. Use the "href" property to provide an external URL. Received: ${JSON.stringify(response.options)}`
          );
        }
        if (["params", "search", "hash"].some(
          (d) => typeof response.options[d] === "function"
        )) {
          throw new Error(
            `Server side redirects must use static search, params, and hash values and do not support functional values. Received functional values for: ${Object.keys(
              response.options
            ).filter((d) => typeof response.options[d] === "function").map((d) => `"${d}"`).join(", ")}`
          );
        }
        const router2 = await getRouter();
        const redirect2 = router2.resolveRedirect(response);
        if (request.headers.get("x-tsr-redirect") === "manual") {
          return json(
            {
              ...response.options,
              isSerializedRedirect: true
            },
            {
              headers: response.headers
            }
          );
        }
        return redirect2;
      }
      return response;
    } finally {
      if (router && !cbWillCleanup) {
        router.serverSsr?.cleanup();
      }
      router = null;
    }
  };
  return requestHandler(startRequestResolver);
}
async function handleServerRoutes({
  getRouter,
  request,
  executeRouter,
  context
}) {
  const router = await getRouter();
  let url = new URL(request.url);
  url = executeRewriteInput(router.rewrite, url);
  const pathname = url.pathname;
  const { matchedRoutes, foundRoute, routeParams } = router.getMatchedRoutes(pathname);
  const middlewares = flattenMiddlewares(
    matchedRoutes.flatMap((r) => r.options.server?.middleware).filter(Boolean)
  ).map((d) => d.options.server);
  const server2 = foundRoute?.options.server;
  if (server2) {
    if (server2.handlers) {
      const handlers = typeof server2.handlers === "function" ? server2.handlers({
        createHandlers: (d) => d
      }) : server2.handlers;
      const requestMethod = request.method.toUpperCase();
      const handler = handlers[requestMethod] ?? handlers["ANY"];
      if (handler) {
        const mayDefer = !!foundRoute.options.component;
        if (typeof handler === "function") {
          middlewares.push(handlerToMiddleware(handler, mayDefer));
        } else {
          const { middleware } = handler;
          if (middleware && middleware.length) {
            middlewares.push(
              ...flattenMiddlewares(middleware).map((d) => d.options.server)
            );
          }
          if (handler.handler) {
            middlewares.push(handlerToMiddleware(handler.handler, mayDefer));
          }
        }
      }
    }
  }
  middlewares.push(
    handlerToMiddleware((ctx2) => executeRouter({ serverContext: ctx2.context }))
  );
  const ctx = await executeMiddleware(middlewares, {
    request,
    context,
    params: routeParams,
    pathname
  });
  const response = ctx.response;
  return response;
}
function throwRouteHandlerError() {
  if (process.env.NODE_ENV === "development") {
    throw new Error(
      `It looks like you forgot to return a response from your server route handler. If you want to defer to the app router, make sure to have a component set in this route.`
    );
  }
  throw new Error("Internal Server Error");
}
function throwIfMayNotDefer() {
  if (process.env.NODE_ENV === "development") {
    throw new Error(
      `You cannot defer to the app router if there is no component defined on this route.`
    );
  }
  throw new Error("Internal Server Error");
}
function handlerToMiddleware(handler, mayDefer = false) {
  if (mayDefer) {
    return handler;
  }
  return async ({ next: _next, ...rest }) => {
    const response = await handler({ ...rest, next: throwIfMayNotDefer });
    if (!response) {
      throwRouteHandlerError();
    }
    return response;
  };
}
function executeMiddleware(middlewares, ctx) {
  let index = -1;
  const next = async (ctx2) => {
    index++;
    const middleware = middlewares[index];
    if (!middleware) return ctx2;
    let result;
    try {
      result = await middleware({
        ...ctx2,
        // Allow the middleware to call the next middleware in the chain
        next: async (nextCtx) => {
          const nextResult = await next({
            ...ctx2,
            ...nextCtx,
            context: {
              ...ctx2.context,
              ...nextCtx?.context || {}
            }
          });
          return Object.assign(ctx2, handleCtxResult(nextResult));
        }
        // Allow the middleware result to extend the return context
      });
    } catch (err) {
      if (isSpecialResponse(err)) {
        result = {
          response: err
        };
      } else {
        throw err;
      }
    }
    return Object.assign(ctx2, handleCtxResult(result));
  };
  return handleCtxResult(next(ctx));
}
function handleCtxResult(result) {
  if (isSpecialResponse(result)) {
    return {
      response: result
    };
  }
  return result;
}
function isSpecialResponse(err) {
  return isResponse(err) || isRedirect(err);
}
function isResponse(response) {
  return response instanceof Response;
}
const fetch$1 = createStartHandler(defaultStreamHandler);
function createServerEntry(entry) {
  return {
    async fetch(...args) {
      return await entry.fetch(...args);
    }
  };
}
const server = createServerEntry({ fetch: fetch$1 });
export {
  ErrorComponent as E,
  RouterCore as R,
  useRouter as a,
  useForwardedRef as b,
  useIntersectionObserver as c,
  createServerEntry,
  dummyMatchContext as d,
  server as default,
  exactPathTest as e,
  functionalUpdate as f,
  removeTrailingSlash as g,
  deepEqual as h,
  invariant as i,
  joinPaths as j,
  isModuleNotFoundError as k,
  json as l,
  matchContext as m,
  rootRouteId as r,
  trimPathLeft as t,
  useRouterState as u,
  warning as w
};
