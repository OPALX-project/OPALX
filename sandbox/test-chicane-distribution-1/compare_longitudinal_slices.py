#!/usr/bin/env python3
"""Compare OPALX and old OPAL exit data by initial longitudinal slice."""

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

from compare_results import read_distribution, read_exit_monitor


def rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def centered_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return rms(values - np.mean(values))


def read_ordered_monitor(path: Path, n_particles: int) -> tuple[str, float, dict[str, np.ndarray]]:
    selected, spos, data = read_exit_monitor(path, n_particles)
    return selected, spos, {key: np.asarray(value, dtype=float) for key, value in data.items()}


def make_rows(
        initial: np.ndarray,
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray]) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    z0_values = np.unique(initial[:, 5])
    for z0 in z0_values:
        mask = np.isclose(initial[:, 5], z0, rtol=0.0, atol=1.0e-15)
        row: dict[str, float | int] = {
            "z0_m": float(z0),
            "n": int(np.count_nonzero(mask)),
            "delta_min": float(np.min(initial[mask, 7])),
            "delta_max": float(np.max(initial[mask, 7])),
        }
        for quantity in ("x", "px", "z", "pz", "y", "py"):
            vx = opalx[quantity][mask]
            vo = opal[quantity][mask]
            diff = vx - vo
            row[f"opalx_mean_{quantity}"] = float(np.mean(vx))
            row[f"opal_mean_{quantity}"] = float(np.mean(vo))
            row[f"diff_mean_{quantity}"] = float(np.mean(diff))
            row[f"opalx_rms_{quantity}"] = centered_rms(vx)
            row[f"opal_rms_{quantity}"] = centered_rms(vo)
            row[f"diff_rms_{quantity}"] = centered_rms(diff)
            row[f"diff_max_abs_{quantity}"] = float(np.max(np.abs(diff)))
        rows.append(row)
    return rows


def plot_two_panel(
        path: Path,
        frame: pd.DataFrame,
        top_columns: list[tuple[str, str, str]],
        bottom_column: str,
        ylabel_top: str,
        ylabel_bottom: str,
        title: str) -> None:
    z_mm = 1.0e3 * frame["z0_m"].to_numpy()
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(9.5, 6.3),
        sharex=True,
        gridspec_kw={"height_ratios": [1.35, 1.0]},
        constrained_layout=True,
    )
    for column, label, color in top_columns:
        ax_top.plot(z_mm, frame[column], marker="o", lw=1.7, ms=4.5, label=label, color=color)
    ax_top.set_ylabel(ylabel_top)
    ax_top.set_title(title)
    ax_top.grid(True, which="major", alpha=0.35)
    ax_top.grid(True, which="minor", alpha=0.18)
    ax_top.minorticks_on()
    ax_top.legend(frameon=False, loc="best")

    ax_bottom.plot(z_mm, frame[bottom_column], marker="o", lw=1.7, ms=4.5, color="#D55E00")
    ax_bottom.axhline(0.0, color="0.35", lw=0.9)
    ax_bottom.set_xlabel(r"initial $z_0$ [mm]")
    ax_bottom.set_ylabel(ylabel_bottom)
    ax_bottom.grid(True, which="major", alpha=0.35)
    ax_bottom.grid(True, which="minor", alpha=0.18)
    ax_bottom.minorticks_on()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def variance_decomposition(frame: pd.DataFrame, prefix: str) -> dict[str, float]:
    weights = frame["n"].to_numpy(dtype=float)
    weights /= np.sum(weights)
    means = frame[f"{prefix}_mean_x"].to_numpy(dtype=float)
    slice_variances = frame[f"{prefix}_rms_x"].to_numpy(dtype=float) ** 2
    within = float(np.sum(weights * slice_variances))
    centroid = float(np.sum(weights * (means - np.sum(weights * means)) ** 2))
    total = within + centroid
    return {
        "sigma_total": float(np.sqrt(total)),
        "sigma_within": float(np.sqrt(within)),
        "sigma_centroid": float(np.sqrt(centroid)),
        "within_var_fraction": within / total if total > 0.0 else float("nan"),
        "centroid_var_fraction": centroid / total if total > 0.0 else float("nan"),
    }


def linear_fit_r2(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    prediction = slope * x + intercept
    ss_res = float(np.sum((y - prediction) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return float(slope), float(intercept), 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")


def write_plots(out_dir: Path, frame: pd.DataFrame) -> None:
    plot_two_panel(
        out_dir / "slice_mean_x.png",
        frame,
        [
            ("opalx_mean_x", "OPALX", "#0072B2"),
            ("opal_mean_x", "OPAL 2022.1", "#009E73"),
        ],
        "diff_mean_x",
        r"slice centroid $\bar{x}$ [m]",
        r"mean $\Delta x$ [m]",
        r"Horizontal slice centroid by initial longitudinal slice",
    )
    plot_two_panel(
        out_dir / "slice_rms_x.png",
        frame,
        [
            ("opalx_rms_x", "OPALX", "#0072B2"),
            ("opal_rms_x", "OPAL 2022.1", "#009E73"),
        ],
        "diff_mean_x",
        r"slice $\sigma_x$ [m]",
        r"mean $\Delta x$ [m]",
        r"Horizontal size by initial longitudinal slice",
    )
    plot_two_panel(
        out_dir / "slice_rms_px.png",
        frame,
        [
            ("opalx_rms_px", "OPALX", "#0072B2"),
            ("opal_rms_px", "OPAL 2022.1", "#009E73"),
        ],
        "diff_mean_px",
        r"slice $\sigma_{p_x}$",
        r"mean $\Delta p_x$",
        r"Horizontal momentum by initial longitudinal slice",
    )
    plot_two_panel(
        out_dir / "slice_rms_z.png",
        frame,
        [
            ("opalx_rms_z", "OPALX", "#0072B2"),
            ("opal_rms_z", "OPAL 2022.1", "#009E73"),
        ],
        "diff_mean_z",
        r"slice $\sigma_z$ [m]",
        r"mean $\Delta z$ [m]",
        r"Longitudinal coordinate by initial longitudinal slice",
    )
    plot_two_panel(
        out_dir / "slice_max_abs_px_diff.png",
        frame,
        [
            ("diff_max_abs_px", r"max $|\Delta p_x|$", "#CC79A7"),
            ("diff_rms_px", r"rms $\Delta p_x$", "#0072B2"),
        ],
        "diff_mean_px",
        r"$p_x$ difference",
        r"mean $\Delta p_x$",
        r"OPALX - OPAL horizontal momentum difference",
    )


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--distribution",
        type=Path,
        default=default_base / "opalx_elemedge" / "test-chicane-distribution-1_distribution.txt",
    )
    parser.add_argument(
        "--opalx-monitor",
        type=Path,
        default=default_base / "opalx_elemedge" / "test-chicane-distribution-2_exit.h5",
    )
    parser.add_argument(
        "--opal-monitor",
        type=Path,
        default=default_base / "opal_2022" / "test-chicane-distribution-1_opal_exit.h5",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=default_base / "longitudinal_slice_compare",
    )
    args = parser.parse_args()

    initial = read_distribution(args.distribution)
    n_particles = len(initial)
    opalx_step, opalx_spos, opalx = read_ordered_monitor(args.opalx_monitor, n_particles)
    opal_step, opal_spos, opal = read_ordered_monitor(args.opal_monitor, n_particles)

    frame = pd.DataFrame(make_rows(initial, opalx, opal))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out_dir / "slice_summary.csv", index=False)
    write_plots(args.out_dir, frame)

    strongest = frame.iloc[np.argmax(np.abs(frame["diff_mean_px"].to_numpy()))]
    largest_rms_x = frame.iloc[np.argmax(np.abs((frame["opalx_rms_x"] - frame["opal_rms_x"]).to_numpy()))]
    opalx_decomp = variance_decomposition(frame, "opalx")
    opal_decomp = variance_decomposition(frame, "opal")
    dpx_slope, dpx_intercept, dpx_r2 = linear_fit_r2(
        frame["z0_m"].to_numpy(dtype=float),
        frame["diff_mean_px"].to_numpy(dtype=float),
    )
    lines = [
        "# Longitudinal Slice Comparison",
        "",
        f"Distribution: `{args.distribution}`",
        f"OPALX monitor: `{args.opalx_monitor}`",
        f"OPAL monitor: `{args.opal_monitor}`",
        f"OPALX selected step: `{opalx_step}`, SPOS `{opalx_spos:.12e} m`",
        f"OPAL selected step: `{opal_step}`, SPOS `{opal_spos:.12e} m`",
        "",
        "Rows are grouped by initial `z0`; each slice contains the full transverse and momentum grid.",
        "",
        "## Key Slice Findings",
        "",
        (
            "- Largest absolute mean horizontal kick difference: "
            f"`z0 = {strongest['z0_m']:.12e} m`, "
            f"`mean(OPALX-OPAL px) = {strongest['diff_mean_px']:.12e}`."
        ),
        (
            "- Largest absolute slice `sigma_x` difference: "
            f"`z0 = {largest_rms_x['z0_m']:.12e} m`, "
            f"`sigma_x(OPALX)-sigma_x(OPAL) = "
            f"{(largest_rms_x['opalx_rms_x'] - largest_rms_x['opal_rms_x']):.12e} m`."
        ),
        (
            "- Linear fit of slice `mean(OPALX-OPAL px)` versus `z0`: "
            f"slope `{dpx_slope:.12e} 1/m`, intercept `{dpx_intercept:.12e}`, "
            f"`R^2 = {dpx_r2:.6f}`."
        ),
        "",
        "## Horizontal Variance Decomposition",
        "",
        "| code | total sigma_x [m] | within-slice sigma [m] | slice-centroid sigma [m] | centroid variance fraction |",
        "| --- | --- | --- | --- | --- |",
        (
            "| OPALX | "
            f"{opalx_decomp['sigma_total']:.9e} | "
            f"{opalx_decomp['sigma_within']:.9e} | "
            f"{opalx_decomp['sigma_centroid']:.9e} | "
            f"{opalx_decomp['centroid_var_fraction']:.9e} |"
        ),
        (
            "| OPAL 2022.1 | "
            f"{opal_decomp['sigma_total']:.9e} | "
            f"{opal_decomp['sigma_within']:.9e} | "
            f"{opal_decomp['sigma_centroid']:.9e} | "
            f"{opal_decomp['centroid_var_fraction']:.9e} |"
        ),
        "",
        "## Compact Table",
        "",
        "| z0 [mm] | sigma_x OPALX [m] | sigma_x OPAL [m] | delta sigma_x [m] | mean delta px | rms delta px | mean delta z [m] |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in frame.iterrows():
        lines.append(
            "| "
            f"{1.0e3 * row['z0_m']:.9e} | "
            f"{row['opalx_rms_x']:.9e} | "
            f"{row['opal_rms_x']:.9e} | "
            f"{row['opalx_rms_x'] - row['opal_rms_x']:.9e} | "
            f"{row['diff_mean_px']:.9e} | "
            f"{row['diff_rms_px']:.9e} | "
            f"{row['diff_mean_z']:.9e} |"
        )
    lines.extend(
        [
            "",
            "## Plots",
            "",
            "- `slice_mean_x.png`",
            "- `slice_rms_x.png`",
            "- `slice_rms_px.png`",
            "- `slice_rms_z.png`",
            "- `slice_max_abs_px_diff.png`",
        ]
    )
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
