#!/usr/bin/env python3
"""Diagnose OPAL/OPALX bend marker and active-interval placement."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re

import numpy as np


MARKER_RE = re.compile(
    r'"(?P<marker>[^"]+)"\s+'
    r"(?P<z>[-+0-9.eE]+)\s+"
    r"(?P<x>[-+0-9.eE]+)\s+"
    r"(?P<y>[-+0-9.eE]+)"
)
SECTION_RE = re.compile(
    r"Key:\s+\((?P<begin>[-+0-9.eE]+)\s+-\s+(?P<end>[-+0-9.eE]+)\)"
)
ELEMENT_RE = re.compile(r"\*\s+(?P<element>[A-Za-z][A-Za-z0-9_]*)\s*$")


def read_markers(path: Path) -> dict[str, np.ndarray]:
    markers: dict[str, np.ndarray] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = MARKER_RE.search(line)
        if match:
            markers[match.group("marker")] = np.array(
                [float(match.group("z")), float(match.group("x")), float(match.group("y"))]
            )
    return markers


def read_section_intervals(path: Path) -> dict[str, tuple[float, float]]:
    intervals: dict[str, list[float]] = {}
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    index = 0
    while index < len(lines):
        section = SECTION_RE.search(lines[index])
        if section is None:
            index += 1
            continue
        begin = float(section.group("begin"))
        end = float(section.group("end"))
        index += 1
        while index < len(lines):
            if SECTION_RE.search(lines[index]):
                break
            element = ELEMENT_RE.search(lines[index])
            if element:
                values = intervals.setdefault(element.group("element"), [begin, end])
                values[0] = min(values[0], begin)
                values[1] = max(values[1], end)
            index += 1
    return {name: (values[0], values[1]) for name, values in intervals.items()}


def marker_interval(
    markers: dict[str, np.ndarray], element: str, begin_kind: str, end_kind: str
) -> tuple[float, float] | None:
    begin_marker = f"{begin_kind}: {element}"
    end_marker = f"{end_kind}: {element}"
    if begin_marker not in markers or end_marker not in markers:
        return None
    return float(markers[begin_marker][0]), float(markers[end_marker][0])


def tangent_from_markers(markers: dict[str, np.ndarray], element: str) -> np.ndarray:
    begin = markers[f"BEGIN: {element}"]
    end = markers[f"END: {element}"]
    tangent = end - begin
    norm = np.linalg.norm(tangent[:2])
    if norm == 0.0:
        return np.array([1.0, 0.0, 0.0])
    # Element marker coordinates are [z, x, y].  The chicane is in z-x.
    return tangent / np.linalg.norm(tangent)


def projected_offset(
    markers: dict[str, np.ndarray], start_marker: str, end_marker: str, tangent_element: str
) -> tuple[float, float]:
    delta = markers[end_marker] - markers[start_marker]
    tangent = tangent_from_markers(markers, tangent_element)
    projected = float(np.dot(delta, tangent))
    perpendicular = float(np.linalg.norm(delta - projected * tangent))
    return projected, perpendicular


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, title: str, rows: list[dict[str, object]]) -> None:
    headers = list(rows[0].keys())
    rendered = []
    for row in rows:
        rendered.append(
            [
                f"{value:.9e}" if isinstance(value, float) else str(value)
                for value in (row[header] for header in headers)
            ]
        )
    widths = [
        max(len(header), *(len(row[index]) for row in rendered))
        for index, header in enumerate(headers)
    ]
    lines = [f"# {title}", ""]
    lines.append("| " + " | ".join(header.ljust(widths[i]) for i, header in enumerate(headers)) + " |")
    lines.append("| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |")
    for row in rendered:
        lines.append("| " + " | ".join(value.ljust(widths[i]) for i, value in enumerate(row)) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--old-log",
        type=Path,
        default=script_dir / "diagnostics" / "old_opal_info" / "old_opal_info.log",
    )
    parser.add_argument(
        "--old-markers",
        type=Path,
        default=script_dir
        / "diagnostics"
        / "old_opal_info"
        / "data"
        / "test-chicane-distribution-1_opal_ElementPositions.txt",
    )
    parser.add_argument(
        "--opalx-markers",
        type=Path,
        default=script_dir
        / "diagnostics"
        / "opalx_elemedge_info"
        / "data"
        / "test-chicane-distribution-2_ElementPositions.txt",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=script_dir / "diagnostics" / "bend_interval_compare",
    )
    args = parser.parse_args()

    old_markers = read_markers(args.old_markers)
    opalx_markers = read_markers(args.opalx_markers)
    old_sections = read_section_intervals(args.old_log)

    bends = ["B1", "B2", "B3", "B4"]
    interval_rows: list[dict[str, object]] = []
    for bend in bends:
        old_marker = marker_interval(old_markers, bend, "BEGIN", "END")
        opalx_body = marker_interval(opalx_markers, bend, "BEGIN", "END")
        opalx_field = marker_interval(opalx_markers, bend, "FIELD BEGIN", "FIELD END")
        if old_marker is None or opalx_body is None or opalx_field is None:
            raise RuntimeError(f"missing marker interval for {bend}")
        old_active = old_sections[bend]
        interval_rows.append(
            {
                "bend": bend,
                "old_active_begin_m": old_active[0],
                "old_active_end_m": old_active[1],
                "old_active_len_m": old_active[1] - old_active[0],
                "old_marker_begin_m": old_marker[0],
                "old_marker_end_m": old_marker[1],
                "old_marker_len_m": old_marker[1] - old_marker[0],
                "opalx_field_begin_m": opalx_field[0],
                "opalx_field_end_m": opalx_field[1],
                "opalx_field_len_m": opalx_field[1] - opalx_field[0],
                "opalx_body_begin_m": opalx_body[0],
                "opalx_body_end_m": opalx_body[1],
                "old_active_minus_opalx_field_begin_mm": 1.0e3 * (old_active[0] - opalx_field[0]),
                "old_active_minus_opalx_field_end_mm": 1.0e3 * (old_active[1] - opalx_field[1]),
            }
        )

    next_elements = {"B1": "D12", "B2": "D23", "B3": "D34", "B4": "MEXIT"}
    interface_rows: list[dict[str, object]] = []
    for bend, next_element in next_elements.items():
        old_projected, old_perp = projected_offset(
            old_markers, f"EXIT EDGE: {bend}", f"BEGIN: {next_element}", bend
        )
        opalx_projected, opalx_perp = projected_offset(
            opalx_markers, f"EXIT EDGE: {bend}", f"BEGIN: {next_element}", bend
        )
        interface_rows.append(
            {
                "interface": f"{bend}->{next_element}",
                "old_projected_gap_mm": 1.0e3 * old_projected,
                "old_perp_gap_mm": 1.0e3 * old_perp,
                "opalx_projected_gap_mm": 1.0e3 * opalx_projected,
                "opalx_perp_gap_mm": 1.0e3 * opalx_perp,
                "projected_gap_diff_mm": 1.0e3 * (opalx_projected - old_projected),
            }
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "bend_intervals.csv", interval_rows)
    write_csv(args.out_dir / "bend_interfaces.csv", interface_rows)
    write_markdown(args.out_dir / "bend_intervals.md", "Bend Intervals", interval_rows)
    write_markdown(args.out_dir / "bend_interfaces.md", "Bend Interfaces", interface_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
