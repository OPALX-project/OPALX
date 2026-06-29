#!/usr/bin/env python3
"""Compare longitudinal-slice shear at terminal temporal section monitors."""

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


SECTIONS = ["after_b1", "after_b2", "after_b3", "after_b4"]


def analyze_section(base_dir: Path, section: str, distribution: Path) -> tuple[dict[str, float | str], pd.DataFrame]:
    section_dir = base_dir / f"section_monitor_{section}"
    initial = read_distribution(distribution)
    n_particles = len(initial)
    opalx_step, opalx_spos, opalx = read_ordered_monitor(
        section_dir / "opalx_elemedge" / f"test-chicane-distribution-2_{section}.h5",
        n_particles,
    )
    opal_step, opal_spos, opal = read_ordered_monitor(
        section_dir / "opal_2022" / f"test-chicane-distribution-1_opal_{section}.h5",
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
        "section": section,
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
    labels = [item.replace("after_", "") for item in summary["section"]]

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.plot(x, summary["opalx_centroid_sigma_x_m"], marker="o", lw=1.8, label="OPALX")
    ax.plot(x, summary["opal_centroid_sigma_x_m"], marker="o", lw=1.8, label="OPAL 2022.1")
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"slice-centroid contribution to monitor $\sigma_x$ [m]")
    ax.set_title(r"Horizontal slice-centroid shear at section temporal monitors")
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(True, which="minor", alpha=0.18)
    ax.minorticks_on()
    ax.legend(frameon=False)
    fig.savefig(out_dir / "monitor_centroid_sigma_x.png", dpi=220)
    fig.savefig(out_dir / "monitor_centroid_sigma_x.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.plot(x, summary["dpx_slope_per_m"], marker="o", lw=1.8, color="#D55E00")
    ax.set_xticks(x, labels)
    ax.axhline(0.0, color="0.35", lw=0.9)
    ax.set_ylabel(r"slope of monitor mean $\Delta p_x$ vs $z_0$ [1/m]")
    ax.set_title(r"Longitudinal-slice horizontal kick shear at temporal monitors")
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(True, which="minor", alpha=0.18)
    ax.minorticks_on()
    fig.savefig(out_dir / "monitor_dpx_slope.png", dpi=220)
    fig.savefig(out_dir / "monitor_dpx_slope.pdf")
    plt.close(fig)


def write_slice_plots(out_dir: Path, frames: dict[str, pd.DataFrame]) -> None:
    specs = [
        (
            "slice_mean_x_by_section",
            "opalx_mean_x",
            "opal_mean_x",
            r"slice mean $x$ [m]",
            r"Longitudinal-slice horizontal centroids at section monitors",
        ),
        (
            "slice_mean_px_by_section",
            "opalx_mean_px",
            "opal_mean_px",
            r"slice mean $p_x$",
            r"Longitudinal-slice horizontal momentum centroids at section monitors",
        ),
        (
            "slice_rms_x_by_section",
            "opalx_rms_x",
            "opal_rms_x",
            r"within-slice $\sigma_x$ [m]",
            r"Longitudinal-slice horizontal widths at section monitors",
        ),
    ]
    for stem, opalx_column, opal_column, ylabel, title in specs:
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True, constrained_layout=True)
        for ax, section in zip(axes.ravel(), SECTIONS):
            frame = frames[section]
            z_mm = 1.0e3 * frame["z0_m"].to_numpy(dtype=float)
            ax.plot(z_mm, frame[opalx_column], marker="o", lw=1.5, label="OPALX")
            ax.plot(z_mm, frame[opal_column], marker="s", lw=1.5, label="OPAL 2022.1")
            ax.set_title(section.replace("after_", "after "))
            ax.set_ylabel(ylabel)
            ax.grid(True, which="major", alpha=0.35)
            ax.grid(True, which="minor", alpha=0.18)
            ax.minorticks_on()
        for ax in axes[-1, :]:
            ax.set_xlabel(r"initial $z_0$ [mm]")
        axes[0, 0].legend(frameon=False)
        fig.suptitle(title)
        fig.savefig(out_dir / f"{stem}.png", dpi=220)
        fig.savefig(out_dir / f"{stem}.pdf")
        plt.close(fig)

    diff_specs = [
        (
            "slice_diff_mean_x_by_section",
            "diff_mean_x",
            r"mean $\Delta x$ [m]",
            r"OPALX - OPAL longitudinal-slice horizontal centroid difference",
        ),
        (
            "slice_diff_mean_px_by_section",
            "diff_mean_px",
            r"mean $\Delta p_x$",
            r"OPALX - OPAL longitudinal-slice horizontal momentum-centroid difference",
        ),
        (
            "slice_diff_rms_x_by_section",
            "diff_rms_x",
            r"centered RMS of $\Delta x$ [m]",
            r"Within-slice centered OPALX - OPAL horizontal difference",
        ),
    ]
    for stem, column, ylabel, title in diff_specs:
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True, constrained_layout=True)
        for ax, section in zip(axes.ravel(), SECTIONS):
            frame = frames[section]
            z_mm = 1.0e3 * frame["z0_m"].to_numpy(dtype=float)
            ax.plot(z_mm, frame[column], marker="o", lw=1.5, color="#D55E00")
            ax.axhline(0.0, color="0.35", lw=0.9)
            ax.set_title(section.replace("after_", "after "))
            ax.set_ylabel(ylabel)
            ax.grid(True, which="major", alpha=0.35)
            ax.grid(True, which="minor", alpha=0.18)
            ax.minorticks_on()
        for ax in axes[-1, :]:
            ax.set_xlabel(r"initial $z_0$ [mm]")
        fig.suptitle(title)
        fig.savefig(out_dir / f"{stem}.png", dpi=220)
        fig.savefig(out_dir / f"{stem}.pdf")
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
    parser.add_argument("--out-dir", type=Path, default=default_base / "section_monitor_slice_compare")
    args = parser.parse_args()

    summaries = []
    frames = {}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for section in SECTIONS:
        summary, frame = analyze_section(args.base_dir, section, args.distribution)
        summaries.append(summary)
        frames[section] = frame
        frame.to_csv(args.out_dir / f"{section}_slice_summary.csv", index=False)

    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(args.out_dir / "summary.csv", index=False)
    write_plots(args.out_dir, summary_frame)
    write_slice_plots(args.out_dir, frames)

    lines = [
        "# Section Temporal Monitor Longitudinal Slice Comparison",
        "",
        "Each run uses the original upstream lattice truncated at one section and adds",
        "a single terminal temporal monitor at that section.  The reported coordinates",
        "therefore use the same monitor-local convention as the final-exit diagnostic.",
        "",
        "| section | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x OPALX [m] | sigma_x OPAL [m] | centroid sigma OPALX [m] | centroid sigma OPAL [m] |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in summary_frame.iterrows():
        lines.append(
            "| "
            f"{row['section']} | "
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
            "- `monitor_centroid_sigma_x.png`",
            "- `monitor_dpx_slope.png`",
            "- `slice_mean_x_by_section.png`",
            "- `slice_diff_mean_x_by_section.png`",
            "- `slice_mean_px_by_section.png`",
            "- `slice_diff_mean_px_by_section.png`",
            "- `slice_rms_x_by_section.png`",
            "- `slice_diff_rms_x_by_section.png`",
        ]
    )
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
