#!/usr/bin/env python3
"""Compare longitudinal-slice shear at truncated chicane stop points."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compare_longitudinal_slices import (
    linear_fit_r2,
    make_rows,
    read_ordered_monitor,
    variance_decomposition,
)
from compare_results import read_distribution


STOPS = ["after_b1", "after_b2", "after_b3", "after_b4"]


def analyze_stop(base_dir: Path, stop: str, distribution: Path) -> tuple[dict[str, float | str], pd.DataFrame]:
    stop_dir = base_dir / f"section_stop_{stop}"
    initial = read_distribution(distribution)
    n_particles = len(initial)
    opalx_step, opalx_spos, opalx = read_ordered_monitor(
        stop_dir / "opalx_elemedge" / "test-chicane-distribution-2.h5",
        n_particles,
    )
    opal_step, opal_spos, opal = read_ordered_monitor(
        stop_dir / "opal_2022" / "test-chicane-distribution-1_opal.h5",
        n_particles,
    )
    frame = pd.DataFrame(make_rows(initial, opalx, opal))
    opalx_decomp = variance_decomposition(frame, "opalx")
    opal_decomp = variance_decomposition(frame, "opal")
    slope, intercept, r2 = linear_fit_r2(
        frame["z0_m"].to_numpy(dtype=float),
        frame["diff_mean_px"].to_numpy(dtype=float),
    )
    summary = {
        "stop": stop,
        "opalx_step": opalx_step,
        "opal_step": opal_step,
        "opalx_spos_m": opalx_spos,
        "opal_spos_m": opal_spos,
        "dpx_slope_per_m": slope,
        "dpx_intercept": intercept,
        "dpx_r2": r2,
        "opalx_sigma_x_m": opalx_decomp["sigma_total"],
        "opal_sigma_x_m": opal_decomp["sigma_total"],
        "sigma_x_diff_m": opalx_decomp["sigma_total"] - opal_decomp["sigma_total"],
        "opalx_within_sigma_x_m": opalx_decomp["sigma_within"],
        "opal_within_sigma_x_m": opal_decomp["sigma_within"],
        "opalx_centroid_sigma_x_m": opalx_decomp["sigma_centroid"],
        "opal_centroid_sigma_x_m": opal_decomp["sigma_centroid"],
        "centroid_sigma_x_diff_m": opalx_decomp["sigma_centroid"] - opal_decomp["sigma_centroid"],
    }
    return summary, frame


def write_plots(out_dir: Path, summary: pd.DataFrame) -> None:
    x = np.arange(len(summary))
    labels = [item.replace("after_", "") for item in summary["stop"]]

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.plot(x, summary["opalx_centroid_sigma_x_m"], marker="o", lw=1.8, label="OPALX")
    ax.plot(x, summary["opal_centroid_sigma_x_m"], marker="o", lw=1.8, label="OPAL 2022.1")
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"slice-centroid contribution to $\sigma_x$ [m]")
    ax.set_title(r"Horizontal slice-centroid shear by stop point")
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(True, which="minor", alpha=0.18)
    ax.minorticks_on()
    ax.legend(frameon=False)
    fig.savefig(out_dir / "section_centroid_sigma_x.png", dpi=220)
    fig.savefig(out_dir / "section_centroid_sigma_x.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.plot(x, summary["dpx_slope_per_m"], marker="o", lw=1.8, color="#D55E00")
    ax.set_xticks(x, labels)
    ax.axhline(0.0, color="0.35", lw=0.9)
    ax.set_ylabel(r"slope of mean $\Delta p_x$ vs $z_0$ [1/m]")
    ax.set_title(r"Growth of longitudinal-slice horizontal kick shear")
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(True, which="minor", alpha=0.18)
    ax.minorticks_on()
    fig.savefig(out_dir / "section_dpx_slope.png", dpi=220)
    fig.savefig(out_dir / "section_dpx_slope.pdf")
    plt.close(fig)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=default_base)
    parser.add_argument(
        "--distribution",
        type=Path,
        default=default_base / "opalx_elemedge" / "test-chicane-distribution-1_distribution.txt",
    )
    parser.add_argument("--out-dir", type=Path, default=default_base / "section_stop_slice_compare")
    args = parser.parse_args()

    summaries = []
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for stop in STOPS:
        summary, frame = analyze_stop(args.base_dir, stop, args.distribution)
        summaries.append(summary)
        frame.to_csv(args.out_dir / f"{stop}_slice_summary.csv", index=False)

    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(args.out_dir / "summary.csv", index=False)
    write_plots(args.out_dir, summary_frame)

    lines = [
        "# Section Stop Longitudinal Slice Comparison",
        "",
        "The original lattice is unchanged. Each run only changes `ZSTOP` and compares",
        "the normal H5 bunch output at the truncated stop.",
        "",
        "| stop | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x OPALX [m] | sigma_x OPAL [m] | centroid sigma OPALX [m] | centroid sigma OPAL [m] |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in summary_frame.iterrows():
        lines.append(
            "| "
            f"{row['stop']} | "
            f"{row['opalx_spos_m']:.9e} | "
            f"{row['opal_spos_m']:.9e} | "
            f"{row['dpx_slope_per_m']:.9e} | "
            f"{row['dpx_r2']:.6f} | "
            f"{row['opalx_sigma_x_m']:.9e} | "
            f"{row['opal_sigma_x_m']:.9e} | "
            f"{row['opalx_centroid_sigma_x_m']:.9e} | "
            f"{row['opal_centroid_sigma_x_m']:.9e} |"
        )
    lines.extend(
        [
            "",
            "## Plots",
            "",
            "- `section_centroid_sigma_x.png`",
            "- `section_dpx_slope.png`",
        ]
    )
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
