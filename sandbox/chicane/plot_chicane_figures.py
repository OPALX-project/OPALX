#!/usr/bin/env python3
"""Create all high-quality chicane comparison figures."""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-mplconfig")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
COMPARE_DIR = REPO / "sandbox" / "test-chicane-distribution-1"
sys.path.insert(0, str(COMPARE_DIR))

from compare_opal_opalx import (  # noqa: E402
    BEAM_COLUMNS,
    DIAGNOSTIC_COLUMNS,
    comparison_rows,
    field_interval_rows,
    markdown_table,
    read_sdds_stat,
    write_csv,
    write_plots,
)


P0_BETA_GAMMA = 1957.950928190185
EXPECTED_R56 = -4.333333333333334e-2


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


def read_distribution(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, skiprows=2)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    particle_id = np.arange(len(data), dtype=int)
    delta = data[:, 5] / P0_BETA_GAMMA - 1.0
    return particle_id, delta


def read_exit_z(path: Path, delta_by_id: np.ndarray) -> dict[str, object]:
    rows = []
    with h5py.File(path, "r") as h5:
        for name in sorted(h5.keys(), key=step_index):
            group = h5[name]
            if not {"id", "z", "pz"}.issubset(group.keys()):
                continue
            particle_id = np.asarray(group["id"], dtype=int)
            if np.any(particle_id < 0) or np.any(particle_id >= len(delta_by_id)):
                continue
            z = np.asarray(group["z"], dtype=float)
            pz = np.asarray(group["pz"], dtype=float)
            delta = delta_by_id[particle_id]
            slope, intercept = np.polyfit(delta, z, 1)
            rows.append(
                {
                    "step": name,
                    "spos": h5_attr_scalar(group, "SPOS"),
                    "particle_id": particle_id,
                    "delta": delta,
                    "z": z,
                    "pz": pz,
                    "slope": slope,
                    "intercept": intercept,
                    "r56": -slope,
                }
            )
    if not rows:
        raise RuntimeError(f"no usable Step# groups found in {path}")
    exit_rows = [row for row in rows if float(row["spos"]) > 0.5]
    return min(exit_rows or rows, key=lambda row: abs(float(row["intercept"])))


def plot_chicane_1(fig_dir: Path, summary_dir: Path) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    _, delta_by_id = read_distribution(ROOT / "runs" / "chicane-1-opalx" / "chicane-1-opalx_distribution.txt")
    opalx = read_exit_z(ROOT / "runs" / "chicane-1-opalx" / "chicane-1-opalx_exit.h5", delta_by_id)
    opal = read_exit_z(ROOT / "runs" / "chicane-1-opal" / "chicane-1-opal_exit.h5", delta_by_id)

    with plt.rc_context(
            {
                "font.family": "serif",
                "mathtext.fontset": "stix",
                "font.size": 10,
                "axes.labelsize": 10,
                "axes.titlesize": 11,
                "axes.linewidth": 0.8,
                "xtick.labelsize": 9,
                "ytick.labelsize": 9,
                "legend.fontsize": 8.5,
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "savefig.dpi": 300,
            }):
        fig, (ax_fit, ax_resid) = plt.subplots(
            2,
            1,
            figsize=(7.2, 5.25),
            sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1.0]},
        )
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.90, hspace=0.08)

        colors = {"OPALX": "#0072B2", "OPAL": "#D55E00"}
        dense_delta = np.linspace(float(np.min(delta_by_id)), float(np.max(delta_by_id)), 200)
        for label, result in (("OPALX", opalx), ("OPAL", opal)):
            delta = np.asarray(result["delta"], dtype=float)
            z = np.asarray(result["z"], dtype=float)
            slope = float(result["slope"])
            intercept = float(result["intercept"])
            fit = slope * dense_delta + intercept
            residual = 1.0e6 * (z - (slope * delta + intercept))
            ax_fit.plot(1.0e3 * delta, 1.0e3 * z, "o", color=colors[label], label=label)
            ax_fit.plot(1.0e3 * dense_delta, 1.0e3 * fit, color=colors[label], lw=1.5)
            ax_resid.plot(1.0e3 * delta, residual, "o-", color=colors[label], lw=1.1)

        ax_fit.axline(
            (1.0e3 * dense_delta[0], 1.0e3 * (-EXPECTED_R56 * dense_delta[0])),
            (1.0e3 * dense_delta[-1], 1.0e3 * (-EXPECTED_R56 * dense_delta[-1])),
            color="#555555",
            lw=1.0,
            ls=(0, (2.0, 2.0)),
            label="small-angle target",
        )
        ax_resid.axhline(0.0, color="#5f5f5f", lw=0.65, alpha=0.65)
        ax_fit.set_ylabel(r"monitor $z$ [mm]")
        ax_resid.set_ylabel(r"fit residual [$\mu$m]")
        ax_resid.set_xlabel(r"$\delta = (p_z-p_0)/p_0$ [$10^{-3}$]")
        for axis in (ax_fit, ax_resid):
            axis.grid(True, which="major", color="#d9d9d9", lw=0.6)
            axis.grid(True, which="minor", color="#eeeeee", lw=0.4)
            axis.minorticks_on()
            axis.spines["top"].set_visible(False)
        ax_fit.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.18))
        fig.savefig(fig_dir / "chicane-1-r56.png", bbox_inches="tight")
        fig.savefig(fig_dir / "chicane-1-r56.pdf", bbox_inches="tight")
        plt.close(fig)

    summary = {
        "expected_r56_m": EXPECTED_R56,
        "opalx": {
            "step": opalx["step"],
            "spos_m": opalx["spos"],
            "r56_m": opalx["r56"],
            "r56_difference_m": float(opalx["r56"]) - EXPECTED_R56,
            "intercept_m": opalx["intercept"],
        },
        "opal": {
            "step": opal["step"],
            "spos_m": opal["spos"],
            "r56_m": opal["r56"],
            "r56_difference_m": float(opal["r56"]) - EXPECTED_R56,
            "intercept_m": opal["intercept"],
        },
    }
    (summary_dir / "chicane-1-r56-summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def plot_chicane_2(fig_dir: Path, summary_dir: Path) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    opalx_stat = ROOT / "runs" / "chicane-2-opalx" / "chicane-2-opalx.stat"
    opal_stat = ROOT / "runs" / "chicane-2-opal" / "chicane-2-opal.stat"
    opalx = read_sdds_stat(opalx_stat)
    opal = read_sdds_stat(opal_stat)
    rows = comparison_rows(opalx, opal, BEAM_COLUMNS + DIAGNOSTIC_COLUMNS)
    intervals = field_interval_rows(opalx, opal)
    write_csv(summary_dir / "chicane-2-stat-summary.csv", rows)
    write_csv(summary_dir / "chicane-2-by-field-intervals.csv", intervals)
    write_plots(opalx, opal, fig_dir)

    summary = [
        "# Chicane 2 OPALX vs OPAL",
        "",
        f"OPALX stat: `{opalx_stat}`",
        f"OPAL stat: `{opal_stat}`",
        "",
        "Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.",
        "",
        markdown_table(rows),
    ]
    (summary_dir / "chicane-2-summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")


def main() -> int:
    figures = ROOT / "figures"
    summaries = ROOT / "summaries"
    plot_chicane_1(figures / "chicane-1", summaries / "chicane-1")
    plot_chicane_2(figures / "chicane-2", summaries / "chicane-2")
    print(figures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
