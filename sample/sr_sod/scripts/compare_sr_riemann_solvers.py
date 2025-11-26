#!/usr/bin/env python3
"""Legacy helper left in place to explain that solver toggles vanished."""

from __future__ import annotations

import sys
from textwrap import dedent


def main() -> None:
    notice = dedent(
        """
        SR-GSPH now always uses the exact relativistic Riemann solver, matching
        the Python reference implementation. The old ITERATIVE vs EXACT
        comparison workflow has been retired, so this helper intentionally does
        not launch any simulations.

        To verify the new pipeline, run one of the SR Sod presets and compare
        the outputs (CSV or animations) against sample/sr_sod/results/sharp or
        the Python driver in src/srgsph/reference.
        """
    ).strip()
    print(notice)
    sys.exit(0)


if __name__ == "__main__":
    main()
