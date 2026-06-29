#!/usr/bin/env python3
"""Generate the deterministic distribution and manufactured target."""

from __future__ import annotations

from pathlib import Path

from manufactured_chicane import write_all


def main() -> int:
    paths = write_all(Path(__file__).resolve().parent)
    print(f"Wrote {paths.distribution}")
    print(f"Wrote {paths.opal_distribution}")
    print(f"Wrote {paths.solution_csv}")
    print(f"Wrote {paths.solution_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
