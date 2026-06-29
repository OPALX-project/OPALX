#!/usr/bin/env python3
"""Compare generated OPALX and OPAL element marker positions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re

import numpy as np


MARKER_RE = re.compile(
    r'"(?P<marker>[^"]+)"\s+'
    r"(?P<c0>[-+0-9.eE]+)\s+"
    r"(?P<c1>[-+0-9.eE]+)\s+"
    r"(?P<c2>[-+0-9.eE]+)"
)


def read_markers(path: Path, keep_midpoints: bool = False) -> dict[str, np.ndarray]:
    markers: dict[str, np.ndarray] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = MARKER_RE.search(line)
        if not match:
            continue
        marker = match.group("marker")
        if marker.startswith("MID:") and not keep_midpoints:
            continue
        markers[marker] = np.array(
            [
                float(match.group("c0")),
                float(match.group("c1")),
                float(match.group("c2")),
            ]
        )
    return markers


def marker_sort_key(marker: str) -> tuple[int, int, str]:
    element_order = {
        "DPRE": 0,
        "B1": 1,
        "D12": 2,
        "B2": 3,
        "D23": 4,
        "B3": 5,
        "D34": 6,
        "B4": 7,
        "MEXIT": 8,
    }
    marker_order = {
        "FIELD BEGIN": 0,
        "ENTRY EDGE": 1,
        "BEGIN": 2,
        "END": 3,
        "EXIT EDGE": 4,
        "FIELD END": 5,
    }
    if ": " not in marker:
        return (999, 999, marker)
    kind, element = marker.split(": ", 1)
    return (element_order.get(element, 999), marker_order.get(kind, 999), marker)


def compare_markers(
    opalx: dict[str, np.ndarray],
    opal: dict[str, np.ndarray],
    element: str | None = None,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for marker in sorted(set(opalx) & set(opal), key=marker_sort_key):
        if element is not None and not marker.endswith(f": {element}"):
            continue
        x_value = opalx[marker]
        o_value = opal[marker]
        diff = x_value - o_value
        rows.append(
            {
                "marker": marker,
                "opalx_c0_m": float(x_value[0]),
                "opalx_c1_m": float(x_value[1]),
                "opalx_c2_m": float(x_value[2]),
                "opal_c0_m": float(o_value[0]),
                "opal_c1_m": float(o_value[1]),
                "opal_c2_m": float(o_value[2]),
                "diff_c0_mm": float(1.0e3 * diff[0]),
                "diff_c1_mm": float(1.0e3 * diff[1]),
                "diff_c2_mm": float(1.0e3 * diff[2]),
                "diff_norm_mm": float(1.0e3 * np.linalg.norm(diff)),
            }
        )
    return rows


def compare_opalx_field_to_opal_body(
    opalx: dict[str, np.ndarray],
    opal: dict[str, np.ndarray],
    element: str | None = None,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    bend_elements = ["B1", "B2", "B3", "B4"] if element is None else [element]
    pairs = [
        ("FIELD BEGIN", "BEGIN"),
        ("FIELD END", "END"),
    ]
    for bend in bend_elements:
        for opalx_kind, opal_kind in pairs:
            opalx_marker = f"{opalx_kind}: {bend}"
            opal_marker = f"{opal_kind}: {bend}"
            if opalx_marker not in opalx or opal_marker not in opal:
                continue
            x_value = opalx[opalx_marker]
            o_value = opal[opal_marker]
            diff = x_value - o_value
            rows.append(
                {
                    "marker": f"OPALX {opalx_marker} vs OPAL {opal_marker}",
                    "opalx_c0_m": float(x_value[0]),
                    "opalx_c1_m": float(x_value[1]),
                    "opalx_c2_m": float(x_value[2]),
                    "opal_c0_m": float(o_value[0]),
                    "opal_c1_m": float(o_value[1]),
                    "opal_c2_m": float(o_value[2]),
                    "diff_c0_mm": float(1.0e3 * diff[0]),
                    "diff_c1_mm": float(1.0e3 * diff[1]),
                    "diff_c2_mm": float(1.0e3 * diff[2]),
                    "diff_norm_mm": float(1.0e3 * np.linalg.norm(diff)),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, float | str]]) -> str:
    headers = [
        "marker",
        "opalx_c0_m",
        "opalx_c1_m",
        "opal_c0_m",
        "opal_c1_m",
        "diff_c0_mm",
        "diff_c1_mm",
        "diff_norm_mm",
    ]
    rendered: list[list[str]] = []
    for row in rows:
        values: list[str] = []
        for header in headers:
            value = row[header]
            if isinstance(value, float):
                values.append(f"{value:.9e}")
            else:
                values.append(str(value))
        rendered.append(values)
    widths = [
        max(len(header), *(len(row[index]) for row in rendered))
        for index, header in enumerate(headers)
    ]
    lines = [
        "| " + " | ".join(header.ljust(widths[index]) for index, header in enumerate(headers)) + " |",
        "| " + " | ".join("-" * width for width in widths) + " |",
    ]
    lines.extend(
        "| " + " | ".join(value.ljust(widths[index]) for index, value in enumerate(row)) + " |"
        for row in rendered
    )
    return "\n".join(lines)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--opalx-markers",
        type=Path,
        default=default_base
        / "opalx_elemedge"
        / "data"
        / "test-chicane-distribution-2_ElementPositions.txt",
    )
    parser.add_argument(
        "--opal-markers",
        type=Path,
        default=default_base
        / "opal_2022"
        / "data"
        / "test-chicane-distribution-1_opal_ElementPositions.txt",
    )
    parser.add_argument("--element", help="Only compare markers for one element, e.g. B1.")
    parser.add_argument(
        "--field-vs-opal-body",
        action="store_true",
        help="Compare OPALX FIELD BEGIN/END markers against OPAL BEGIN/END markers for bends.",
    )
    parser.add_argument("--out-dir", type=Path, default=default_base / "marker_compare")
    args = parser.parse_args()

    opalx = read_markers(args.opalx_markers)
    opal = read_markers(args.opal_markers)
    rows = (
        compare_opalx_field_to_opal_body(opalx, opal, element=args.element)
        if args.field_vs_opal_body
        else compare_markers(opalx, opal, element=args.element)
    )
    if not rows:
        raise RuntimeError("no matching markers found")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_field_vs_opal_body" if args.field_vs_opal_body else ""
    if args.element:
        suffix += f"_{args.element}"
    csv_path = args.out_dir / f"element_marker_comparison{suffix}.csv"
    md_path = args.out_dir / f"element_marker_comparison{suffix}.md"
    write_csv(csv_path, rows)
    markdown = "\n".join(
        [
            "# Element marker comparison",
            "",
            f"OPALX markers: `{args.opalx_markers}`",
            f"OPAL markers: `{args.opal_markers}`",
            "",
            (
                "Differences are `OPALX - OPAL`.  Coordinates are the three columns written by the generated "
                "element-position files."
            ),
            "",
            markdown_table(rows),
            "",
        ]
    )
    md_path.write_text(markdown, encoding="utf-8")
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
