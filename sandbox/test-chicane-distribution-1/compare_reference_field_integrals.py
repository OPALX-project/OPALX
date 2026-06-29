#!/usr/bin/env python3
"""Compare OPALX and old OPAL reference-field integrals on common s grids."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compare_opal_opalx import read_sdds_stat


def sorted_series(data: dict[str, np.ndarray], field: str) -> tuple[np.ndarray, np.ndarray]:
    s = np.asarray(data["s"], dtype=float)
    values = np.asarray(data[field], dtype=float)
    finite = np.isfinite(s) & np.isfinite(values)
    order = np.argsort(s[finite])
    return s[finite][order], values[finite][order]


def field_intervals(
        s: np.ndarray,
        values: np.ndarray,
        threshold: float) -> list[tuple[float, float]]:
    active = np.abs(values) > threshold
    if not active.any():
        return []
    starts = np.where(active & np.r_[True, ~active[:-1]])[0]
    ends = np.where(active & np.r_[~active[1:], True])[0]
    return [(float(s[start]), float(s[end])) for start, end in zip(starts, ends)]


def cumulative_trapezoid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x)
    if x.size > 1:
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def markdown_table(frame: pd.DataFrame) -> str:
    headers = list(frame.columns)
    rendered = []
    for _, row in frame.iterrows():
        rendered.append(
            [
                f"{value:.9e}" if isinstance(value, float) else str(value)
                for value in row.tolist()
            ]
        )
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


def plot_bend(
        path: Path,
        bend: str,
        s_grid: np.ndarray,
        bx: np.ndarray,
        bo: np.ndarray,
        ix: np.ndarray,
        io: np.ndarray,
        field: str) -> None:
    component = field.split("_", maxsplit=1)[0]
    fig, (ax_field, ax_int, ax_diff) = plt.subplots(
        3,
        1,
        figsize=(10.5, 7.4),
        sharex=True,
        gridspec_kw={"height_ratios": [1.35, 1.0, 0.8]},
        constrained_layout=True,
    )
    ax_field.plot(s_grid, bx, color="#0072B2", lw=1.8, label=f"OPALX {field}")
    ax_field.plot(s_grid, bo, color="#D55E00", lw=1.5, ls="--", label=f"OPAL {field}")
    ax_field.set_ylabel(f"{component} [T]")
    ax_field.set_title(f"{field} and cumulative integral on common grid: {bend}")
    ax_field.grid(True, which="major", alpha=0.35)
    ax_field.grid(True, which="minor", alpha=0.18)
    ax_field.minorticks_on()
    ax_field.ticklabel_format(axis="y", useOffset=False)
    ax_field.legend(frameon=False, loc="best")

    ax_int.plot(s_grid, ix, color="#0072B2", lw=1.8, label=r"OPALX $\int B_y ds$")
    ax_int.plot(s_grid, io, color="#D55E00", lw=1.5, ls="--", label=r"OPAL $\int B_y ds$")
    ax_int.axhline(0.0, color="0.4", lw=0.8)
    ax_int.set_ylabel(r"$\int B\,ds$ [T m]")
    ax_int.grid(True, which="major", alpha=0.35)
    ax_int.grid(True, which="minor", alpha=0.18)
    ax_int.minorticks_on()
    ax_int.ticklabel_format(axis="y", useOffset=False)
    ax_int.legend(frameon=False, loc="best")

    ax_diff.plot(s_grid, ix - io, color="0.15", lw=1.4)
    ax_diff.axhline(0.0, color="0.4", lw=0.8)
    ax_diff.set_xlabel("s [m]")
    ax_diff.set_ylabel(r"$\Delta\int B\,ds$ [T m]")
    ax_diff.grid(True, which="major", alpha=0.35)
    ax_diff.grid(True, which="minor", alpha=0.18)
    ax_diff.minorticks_on()
    ax_diff.ticklabel_format(axis="y", useOffset=False)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
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
    parser.add_argument("--field", default="By_ref")
    parser.add_argument("--threshold", type=float, default=1.0e-6)
    parser.add_argument("--samples-per-bend", type=int, default=4001)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=default_base / "reference_field_integrals",
    )
    args = parser.parse_args()

    opalx = read_sdds_stat(args.opalx_stat)
    opal = read_sdds_stat(args.opal_stat)
    if args.field not in opalx or args.field not in opal:
        print(f"missing field column {args.field}", file=sys.stderr)
        return 2

    sx, bx = sorted_series(opalx, args.field)
    so, bo = sorted_series(opal, args.field)
    x_intervals = field_intervals(sx, bx, args.threshold)
    o_intervals = field_intervals(so, bo, args.threshold)
    count = min(len(x_intervals), len(o_intervals))
    if count == 0:
        print(f"no intervals found above {args.threshold:g} T", file=sys.stderr)
        return 1

    rows = []
    for index in range(count):
        bend = f"B{index + 1}"
        lo = min(x_intervals[index][0], o_intervals[index][0])
        hi = max(x_intervals[index][1], o_intervals[index][1])
        s_grid = np.linspace(lo, hi, args.samples_per_bend)
        bx_grid = np.interp(s_grid, sx, bx)
        bo_grid = np.interp(s_grid, so, bo)
        ix = cumulative_trapezoid(s_grid, bx_grid)
        io = cumulative_trapezoid(s_grid, bo_grid)
        diff = ix - io
        rows.append(
            {
                "bend": bend,
                "s_start_m": lo,
                "s_end_m": hi,
                "opalx_integral_Tm": ix[-1],
                "opal_integral_Tm": io[-1],
                "integral_diff_Tm": diff[-1],
                "integral_rel_diff": diff[-1] / io[-1] if io[-1] != 0.0 else np.nan,
                "max_abs_cumulative_diff_Tm": float(np.max(np.abs(diff))),
                "max_abs_field_diff_T": float(np.max(np.abs(bx_grid - bo_grid))),
            }
        )
        plot_bend(
            args.out_dir / f"{args.field.lower()}_{bend.lower()}_integral.png",
            bend,
            s_grid,
            bx_grid,
            bo_grid,
            ix,
            io,
            args.field,
        )

    frame = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out_dir / "summary.csv", index=False)

    total_x = float(frame["opalx_integral_Tm"].sum())
    total_o = float(frame["opal_integral_Tm"].sum())
    total_diff = total_x - total_o
    summary = [
        "# Reference Field Integral Comparison",
        "",
        f"Field column: `{args.field}`",
        f"Threshold for interval detection: `{args.threshold:g} T`",
        f"Samples per bend: `{args.samples_per_bend}`",
        "",
        markdown_table(frame),
        "",
        "## Total Over Four Bend Intervals",
        "",
        f"- OPALX: `{total_x:.12e} T m`",
        f"- OPAL: `{total_o:.12e} T m`",
        f"- Difference: `{total_diff:.12e} T m`",
        f"- Relative difference: `{total_diff / total_o:.12e}`" if total_o != 0.0 else "- Relative difference: `nan`",
        "",
        "This diagnostic compares stat-file reference fields only. It does not yet",
        "force both codes to evaluate fields on an externally prescribed 3D trajectory.",
    ]
    (args.out_dir / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
