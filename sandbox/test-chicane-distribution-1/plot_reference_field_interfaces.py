#!/usr/bin/env python3
"""Plot reference-particle dipole field around chicane bend interfaces."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compare_opal_opalx import read_sdds_stat


MARKER_RE = re.compile(
    r'"(?P<marker>[^"]+)"\s+'
    r"(?P<z>[-+0-9.eE]+)\s+"
    r"(?P<x>[-+0-9.eE]+)\s+"
    r"(?P<y>[-+0-9.eE]+)"
)


def read_markers(path: Path) -> dict[str, float]:
    markers: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = MARKER_RE.search(line)
        if match:
            markers[match.group("marker")] = float(match.group("z"))
    return markers


def sorted_xy(data: dict[str, np.ndarray], field: str) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(data["s"])
    return np.asarray(data["s"])[order], np.asarray(data[field])[order]


def draw_marker_lines(ax: plt.Axes, markers: dict[str, float], code: str, bend: str) -> None:
    styles = {
        "FIELD BEGIN": (":", "tab:blue"),
        "BEGIN": ("--", "tab:blue"),
        "END": ("--", "tab:red"),
        "FIELD END": (":", "tab:red"),
        "EXIT EDGE": ("-", "0.35"),
    }
    for kind, (linestyle, color) in styles.items():
        marker = f"{kind}: {bend}"
        if marker not in markers:
            continue
        ax.axvline(
            markers[marker],
            color=color,
            linestyle=linestyle,
            linewidth=0.9,
            alpha=0.55 if code == "OPALX" else 0.35,
        )


def plot_window(
    out_path: Path,
    opalx: dict[str, np.ndarray],
    opal: dict[str, np.ndarray],
    opalx_markers: dict[str, float],
    opal_markers: dict[str, float],
    center: float,
    width: float,
    title: str,
    bend: str,
    field: str,
) -> None:
    sx, bx = sorted_xy(opalx, field)
    so, bo = sorted_xy(opal, field)
    lo = center - 0.5 * width
    hi = center + 0.5 * width
    mask_x = (sx >= lo) & (sx <= hi)
    mask_o = (so >= lo) & (so <= hi)

    fig, (ax_field, ax_diff) = plt.subplots(
        2,
        1,
        figsize=(10.5, 6.2),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
        constrained_layout=True,
    )
    component = field.split("_", maxsplit=1)[0]
    ax_field.plot(
        sx[mask_x],
        bx[mask_x],
        color="#009E73",
        lw=2.0,
        label=rf"OPALX ${component}_{{\rm ref}}$",
    )
    ax_field.plot(
        so[mask_o],
        bo[mask_o],
        color="#009E73",
        lw=1.6,
        ls="--",
        label=rf"OPAL ${component}_{{\rm ref}}$",
    )
    draw_marker_lines(ax_field, opalx_markers, "OPALX", bend)
    draw_marker_lines(ax_field, opal_markers, "OPAL", bend)
    ax_field.set_ylabel(rf"${component}_{{\rm ref}}$ [T]")
    ax_field.set_title(title)
    ax_field.grid(True, which="major", alpha=0.35)
    ax_field.grid(True, which="minor", alpha=0.18)
    ax_field.minorticks_on()
    ax_field.legend(loc="best", frameon=False)

    common = sx[mask_x]
    if common.size:
        bo_interp = np.interp(common, so, bo)
        ax_diff.plot(common, bx[mask_x] - bo_interp, color="#D55E00", lw=1.5)
    ax_diff.axhline(0.0, color="0.3", lw=0.8)
    ax_diff.set_xlabel("s [m]")
    ax_diff.set_ylabel(rf"$\Delta {component}$ [T]")
    ax_diff.grid(True, which="major", alpha=0.35)
    ax_diff.grid(True, which="minor", alpha=0.18)
    ax_diff.minorticks_on()
    ax_diff.set_xlim(lo, hi)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--opalx-stat",
        type=Path,
        default=default_base / "opalx_elemedge" / "test-chicane-distribution-2.stat",
    )
    parser.add_argument(
        "--opal-stat",
        type=Path,
        default=default_base / "opal_2022" / "test-chicane-distribution-1_opal.stat",
    )
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
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=default_base / "reference_field_interface_plots",
    )
    parser.add_argument("--width", type=float, default=0.18)
    parser.add_argument("--field", choices=("Bx_ref", "By_ref", "Bz_ref"), default="Bz_ref")
    args = parser.parse_args()

    opalx = read_sdds_stat(args.opalx_stat)
    opal = read_sdds_stat(args.opal_stat)
    opalx_markers = read_markers(args.opalx_markers)
    opal_markers = read_markers(args.opal_markers)

    windows = [
        ("B1", "exit", opalx_markers["EXIT EDGE: B1"]),
        ("B2", "entry", opalx_markers["BEGIN: B2"]),
        ("B2", "exit", opalx_markers["EXIT EDGE: B2"]),
        ("B3", "entry", opalx_markers["BEGIN: B3"]),
        ("B3", "exit", opalx_markers["EXIT EDGE: B3"]),
        ("B4", "entry", opalx_markers["BEGIN: B4"]),
        ("B4", "exit", opalx_markers["EXIT EDGE: B4"]),
    ]
    for bend, edge, center in windows:
        plot_window(
            args.out_dir / f"{args.field.lower()}_{bend.lower()}_{edge}.png",
            opalx,
            opal,
            opalx_markers,
            opal_markers,
            center,
            args.width,
            f"{args.field} around {bend} {edge}",
            bend,
            args.field,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
