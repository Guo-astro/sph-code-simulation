import { jsx } from "react/jsx-runtime";
import { ErrorComponent } from "@tanstack/react-router";
function UserErrorComponent({ error }) {
  return /* @__PURE__ */ jsx(ErrorComponent, { error });
}
const SplitErrorComponent = UserErrorComponent;
export {
  SplitErrorComponent as errorComponent
};
