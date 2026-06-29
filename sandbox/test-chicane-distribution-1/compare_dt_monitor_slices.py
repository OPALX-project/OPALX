#!/usr/bin/env python3
"""Summarize exit-monitor longitudinal-slice diagnostics across timestep scans."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-mplconfig")

import h5py
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


def parse_dt(path: Path) -> float:
    match = re.fullmatch(r"dt_([0-9]+)p([0-9]+)em([0-9]+)", path.name)
    if not match:
        raise ValueError(f"cannot parse timestep directory name: {path.name}")
    mantissa = float(f"{match.group(1)}.{match.group(2)}")
    return mantissa * 10.0 ** (-int(match.group(3)))


def step_index(name: str) -> int:
    match = re.fullmatch(r"Step#(\d+)", name)
    return int(match.group(1)) if match else -1


def h5_attr_scalar(group: h5py.Group, name: str) -> float:
    if name not in group.attrs:
        return float("nan")
    value = np.asarray(group.attrs[name], dtype=float)
    if value.size == 0:
        return float("nan")
    return float(value.reshape(-1)[0])


def h5_steps(path: Path) -> str:
    parts = []
    with h5py.File(path, "r") as h5:
        for name in sorted(h5.keys(), key=step_index):
            group = h5[name]
            n_particles = len(group["id"]) if "id" in group else -1
            spos = h5_attr_scalar(group, "SPOS")
            parts.append(f"{name}:n={n_particles}:s={spos:.12e}")
    return "; ".join(parts)


def analyze_dt(dt_dir: Path) -> dict[str, float | str]:
    distribution = dt_dir / "opalx_elemedge" / "test-chicane-distribution-1_distribution.txt"
    opalx_h5 = dt_dir / "opalx_elemedge" / "test-chicane-distribution-2_exit.h5"
    opal_h5 = dt_dir / "opal_2022" / "test-chicane-distribution-1_opal_exit.h5"
    for path in (distribution, opalx_h5, opal_h5):
        if not path.exists():
            raise FileNotFoundError(path)

    initial = read_distribution(distribution)
    n_particles = len(initial)
    opalx_step, opalx_spos, opalx = read_ordered_monitor(opalx_h5, n_particles)
    opal_step, opal_spos, opal = read_ordered_monitor(opal_h5, n_particles)
    frame = pd.DataFrame(make_rows(initial, opalx, opal))

    slope, intercept, r2 = linear_fit_r2(
        frame["z0_m"].to_numpy(dtype=float),
        frame["diff_mean_px"].to_numpy(dtype=float),
    )
    dx_slope, dx_intercept, dx_r2 = linear_fit_r2(
        frame["z0_m"].to_numpy(dtype=float),
        frame["diff_mean_x"].to_numpy(dtype=float),
    )
    opalx_decomp = variance_decomposition(frame, "opalx")
    opal_decomp = variance_decomposition(frame, "opal")

    return {
        "dt_s": parse_dt(dt_dir),
        "dt_dir": dt_dir.name,
        "n_particles": n_particles,
        "opalx_step": opalx_step,
        "opal_step": opal_step,
        "opalx_spos_m": opalx_spos,
        "opal_spos_m": opal_spos,
        "spos_diff_m": opalx_spos - opal_spos,
        "dpx_slope_per_m": slope,
        "dpx_intercept": intercept,
        "dpx_r2": r2,
        "dx_slope": dx_slope,
        "dx_intercept": dx_intercept,
        "dx_r2": dx_r2,
        "opalx_sigma_x_m": opalx_decomp["sigma_total"],
        "opal_sigma_x_m": opal_decomp["sigma_total"],
        "sigma_x_rel_diff_percent": 100.0
        * (opalx_decomp["sigma_total"] - opal_decomp["sigma_total"])
        / opal_decomp["sigma_total"],
        "opalx_centroid_sigma_x_m": opalx_decomp["sigma_centroid"],
        "opal_centroid_sigma_x_m": opal_decomp["sigma_centroid"],
        "centroid_sigma_x_rel_diff_percent": 100.0
        * (opalx_decomp["sigma_centroid"] - opal_decomp["sigma_centroid"])
        / opal_decomp["sigma_centroid"],
        "opalx_h5_steps": h5_steps(opalx_h5),
        "opal_h5_steps": h5_steps(opal_h5),
    }


def write_plots(out_dir: Path, summary: pd.DataFrame) -> None:
    ordered = summary.sort_values("dt_s")

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.6), sharex=True, constrained_layout=True)
    axes[0].plot(ordered["dt_s"], ordered["dpx_slope_per_m"], marker="o", lw=1.8)
    axes[0].axhline(0.0, color="0.35", lw=0.9)
    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"slope of mean $\Delta p_x$ vs $z_0$ [1/m]")
    axes[0].grid(True, which="both", alpha=0.3)

    axes[1].plot(ordered["dt_s"], ordered["sigma_x_rel_diff_percent"], marker="o", lw=1.8)
    axes[1].axhline(0.0, color="0.35", lw=0.9)
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"DT [s]")
    axes[1].set_ylabel(r"$\sigma_x$ relative difference [%]")
    axes[1].grid(True, which="both", alpha=0.3)
    fig.savefig(out_dir / "dt_monitor_slice_summary.png", dpi=220)
    fig.savefig(out_dir / "dt_monitor_slice_summary.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.8), constrained_layout=True)
    ax.plot(ordered["dt_s"], 1.0e3 * ordered["spos_diff_m"], marker="o", lw=1.8)
    ax.axhline(0.0, color="0.35", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlabel(r"DT [s]")
    ax.set_ylabel(r"OPALX - OPAL monitor SPOS [mm]")
    ax.grid(True, which="both", alpha=0.3)
    fig.savefig(out_dir / "dt_monitor_spos_difference.png", dpi=220)
    fig.savefig(out_dir / "dt_monitor_spos_difference.pdf")
    plt.close(fig)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_scan = script_dir / "comparison" / "dt_scan" / "finite"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", type=Path, default=default_scan)
    parser.add_argument("--out-dir", type=Path, default=default_scan / "dt_monitor_slice_compare")
    args = parser.parse_args()

    rows = []
    skipped = []
    for dt_dir in sorted(args.scan_dir.glob("dt_*"), key=parse_dt):
        try:
            rows.append(analyze_dt(dt_dir))
        except FileNotFoundError as exc:
            skipped.append(f"{dt_dir.name}: missing {exc.filename}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows).sort_values("dt_s")
    summary.to_csv(args.out_dir / "summary.csv", index=False)
    write_plots(args.out_dir, summary)

    lines = [
        "# Exit Monitor Longitudinal-Slice DT Summary",
        "",
        "This uses existing finite-distribution OPALX ELEMEDGE and OPAL 2022.1",
        "exit-monitor H5 files. Incomplete timestep directories are skipped.",
        "",
    ]
    if skipped:
        lines.extend(["Skipped:", ""])
        lines.extend(f"- {item}" for item in skipped)
        lines.append("")
    lines.extend(
        [
            "| DT [s] | n | OPALX step | OPAL step | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x rel diff [%] | centroid sigma rel diff [%] |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for _, row in summary.iterrows():
        lines.append(
            "| "
            f"{row['dt_s']:.9e} | "
            f"{int(row['n_particles'])} | "
            f"{row['opalx_step']} | "
            f"{row['opal_step']} | "
            f"{row['opalx_spos_m']:.12e} | "
            f"{row['opal_spos_m']:.12e} | "
            f"{row['dpx_slope_per_m']:.9e} | "
            f"{row['dpx_r2']:.6f} | "
            f"{row['sigma_x_rel_diff_percent']:.6e} | "
            f"{row['centroid_sigma_x_rel_diff_percent']:.6e} |"
        )
    lines.extend(
        [
            "",
            "## H5 Steps",
            "",
        ]
    )
    for _, row in summary.iterrows():
        lines.extend(
            [
                f"### DT {row['dt_s']:.9e} s",
                "",
                f"- OPALX: `{row['opalx_h5_steps']}`",
                f"- OPAL: `{row['opal_h5_steps']}`",
                "",
            ]
        )
    lines.extend(
        [
            "## Plots",
            "",
            "- `dt_monitor_slice_summary.png`",
            "- `dt_monitor_spos_difference.png`",
        ]
    )
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
