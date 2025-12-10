import { jsx } from "react/jsx-runtime";
import { ErrorComponent } from "@tanstack/react-router";
function PostErrorComponent({ error }) {
  return /* @__PURE__ */ jsx(ErrorComponent, { error });
}
export {
  PostErrorComponent as P
};
