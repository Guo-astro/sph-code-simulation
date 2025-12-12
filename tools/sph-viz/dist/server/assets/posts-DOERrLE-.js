import { b as createServerRpc, c as createServerFn } from "../server.js";
import { notFound } from "@tanstack/react-router";
import "node:async_hooks";
import "react/jsx-runtime";
import "@tanstack/react-router/ssr/server";
const fetchPost_createServerFn_handler = createServerRpc("0029094260fc8f554fa3ac223696de0e9591567ec6420250e896c91244c812c5", (opts, signal) => {
  return fetchPost.__executeServer(opts, signal);
});
const fetchPost = createServerFn({
  method: "POST"
}).inputValidator((d) => d).handler(fetchPost_createServerFn_handler, async ({
  data,
  context
}) => {
  console.log("Request context:", context);
  console.info(`Fetching post with id ${data}...`);
  const res = await fetch(`https://jsonplaceholder.typicode.com/posts/${data}`);
  if (!res.ok) {
    if (res.status === 404) {
      throw notFound();
    }
    throw new Error("Failed to fetch post");
  }
  const post = await res.json();
  return post;
});
const fetchPosts_createServerFn_handler = createServerRpc("cbb8ca69048418e62742f2c511faa56326b80ace384144a35bb3e0bf5e8124be", (opts, signal) => {
  return fetchPosts.__executeServer(opts, signal);
});
const fetchPosts = createServerFn().handler(fetchPosts_createServerFn_handler, async () => {
  console.info("Fetching posts...");
  const res = await fetch("https://jsonplaceholder.typicode.com/posts");
  if (!res.ok) {
    throw new Error("Failed to fetch posts");
  }
  const posts = await res.json();
  return posts.slice(0, 10);
});
export {
  fetchPost_createServerFn_handler,
  fetchPosts_createServerFn_handler
};
