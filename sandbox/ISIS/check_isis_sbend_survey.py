#!/usr/bin/env python3
"""Check the exported ISIS SBEND placement and make a survey plot."""

from __future__ import annotations

import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-isis-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/opalx-isis-cache")

import matplotlib
import pandas as pd

from generate_isis_sbend_survey import ROOT, source_elements

matplotlib.use("Agg")
import matplotlib.pyplot as plt


POSITIONS = ROOT / "data/isis_sbend_survey_ElementPositions.txt"
SUMMARY = ROOT / "isis_sbend_survey_summary.csv"
REPORT = ROOT / "isis_sbend_survey_report.md"
PLOT_PNG = ROOT / "data/isis_sbend_survey.png"
PLOT_PDF = ROOT / "data/isis_sbend_survey.pdf"


def read_positions() -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    pattern = re.compile(
        r'^"([^"]+)"\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)'
    )
    for line in POSITIONS.read_text().splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        label, z, x, y = match.groups()
        rows.append({"label": label, "x_m": float(x), "y_m": float(y), "z_m": float(z)})
    if not rows:
        raise RuntimeError(f"no positions found in {POSITIONS}")
    return pd.DataFrame(rows)


def unique_position(frame: pd.DataFrame, label: str) -> pd.Series:
    selected = frame.loc[frame["label"] == label]
    if len(selected) != 1:
        raise RuntimeError(f"expected one '{label}' row, found {len(selected)}")
    return selected.iloc[0]


def main() -> None:
    positions = read_positions()
    elements = source_elements()
    previous = unique_position(positions, "BEGIN: MSTART")
    records: list[dict[str, float | str]] = []
    path_parts: list[pd.DataFrame] = []
    bend_centers: list[tuple[str, float, float]] = []

    for element in elements:
        name = element.name.upper()
        is_bend = element.angle != 0.0
        entry_label = f"ENTRY EDGE: {name}" if is_bend else f"BEGIN: {name}"
        exit_label = f"EXIT EDGE: {name}" if is_bend else f"END: {name}"
        entry = unique_position(positions, entry_label)
        exit_ = unique_position(positions, exit_label)
        gap = math.hypot(float(entry.x_m - previous.x_m), float(entry.z_m - previous.z_m))
        records.append(
            {
                "name": element.name,
                "source_type": element.kind,
                "survey_type": "SBEND" if is_bend else "DRIFT",
                "length_m": element.length,
                "source_elemedge_m": element.elemedge,
                "entry_x_m": float(entry.x_m),
                "entry_z_m": float(entry.z_m),
                "exit_x_m": float(exit_.x_m),
                "exit_z_m": float(exit_.z_m),
                "gap_from_previous_m": gap,
            }
        )

        labels = positions["label"]
        if is_bend:
            selected = positions.loc[
                labels.isin(
                    [f"BEGIN: {name}", f"MID: {name}", f"END: {name}"]
                )
            ]
            middle = selected.iloc[len(selected) // 2]
            bend_centers.append((name, float(middle.z_m), float(middle.x_m)))
        else:
            selected = positions.loc[labels.isin([entry_label, exit_label])]
        path_parts.append(selected)
        previous = exit_

    summary = pd.DataFrame.from_records(records)
    summary.to_csv(SUMMARY, index=False)

    start = unique_position(positions, "BEGIN: MSTART")
    close = unique_position(positions, "BEGIN: MCLOSE")
    closure_x = float(close.x_m - start.x_m)
    closure_y = float(close.y_m - start.y_m)
    closure_z = float(close.z_m - start.z_m)
    closure_norm = math.sqrt(closure_x**2 + closure_y**2 + closure_z**2)
    max_gap = float(summary["gap_from_previous_m"].max())
    final_gap = math.hypot(float(previous.x_m - close.x_m), float(previous.z_m - close.z_m))

    # The text export rounds coordinates at roughly the 10 nm level for this
    # 50 m footprint, so allow two units in its last printed decimal place.
    if closure_norm > 1.0e-9 or max_gap > 2.0e-8 or final_gap > 2.0e-8:
        raise RuntimeError(
            "placement check failed: "
            f"closure={closure_norm:.3e} m, max_gap={max_gap:.3e} m, "
            f"last_to_close={final_gap:.3e} m"
        )

    path = pd.concat(path_parts, ignore_index=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(7.4, 7.0), constrained_layout=True)
    ax.plot(path["z_m"], path["x_m"], color="#1769aa", linewidth=1.8)
    ax.scatter([start.z_m], [start.x_m], s=55, marker="o", color="#d1495b", zorder=4)
    for index, (_, z, x) in enumerate(bend_centers, start=1):
        ax.text(z, x, str(index), fontsize=8, ha="center", va="center", color="#222222")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Global Z [m]")
    ax.set_ylabel("Global X [m]")
    ax.set_title("ISIS native-SBEND placement survey")
    ax.text(
        0.02,
        0.98,
        f"192 elements, 10 bends\nclosure = {closure_norm:.1e} m",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9},
    )
    fig.savefig(PLOT_PNG, dpi=220)
    fig.savefig(PLOT_PDF)
    plt.close(fig)

    report = f"""# ISIS SBEND survey report

- Source elements: {len(elements)}
- Native SBENDs: {int((summary['survey_type'] == 'SBEND').sum())}
- Reference-path length: {summary['length_m'].sum():.8f} m
- Closure delta X: {closure_x:.3e} m
- Closure delta Y: {closure_y:.3e} m
- Closure delta Z: {closure_z:.3e} m
- Closure norm: {closure_norm:.3e} m
- Maximum interface gap: {max_gap:.3e} m
- Final bend exit to close-monitor gap: {final_gap:.3e} m
- Survey bounds X: [{path['x_m'].min():.6f}, {path['x_m'].max():.6f}] m
- Survey bounds Z: [{path['z_m'].min():.6f}, {path['z_m'].max():.6f}] m

Result: PASS (closure <= 1 nm and every component interface gap <= 20 nm,
the resolution limit of the text placement export at this footprint).
"""
    REPORT.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
