#!/usr/bin/env python3
"""Compare exact-timed track12 results from one and multiple MPI ranks."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DEFAULT_ONE_RANK = HERE / "a100_1rank_400k_1024x128x128_refpart_fix"
DEFAULT_FOUR_RANK = HERE / "a100_400k_1024x128x128_refpart_fix"
COLORS = {"electron": "#b51f2e", "positron": "#1f5fb5"}
LINE_STYLES = ["-", (0, (8, 6)), (0, (5, 5)), (0, (10, 8)), (0, (2, 3)), (0, (9, 4, 2, 4))]


def load_comparison(run_dir: Path) -> pd.DataFrame:
    return pd.read_csv(run_dir / "results" / "track12_pointwise_comparison.csv")


def configure_matplotlib(output_dir: Path) -> None:
    cache = output_dir / ".plot-cache-ranks"
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    cache.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def merge_runs(one_rank: pd.DataFrame, four_rank: pd.DataFrame) -> pd.DataFrame:
    keys = ["species", "pair", "reference_step"]
    merged = one_rank.merge(
        four_rank,
        on=keys,
        suffixes=("_1rank", "_4rank"),
        validate="one_to_one",
    )
    if len(merged) != len(one_rank) or len(merged) != len(four_rank):
        raise ValueError("rank-count runs do not contain the same trajectory samples")
    return merged.sort_values(keys).reset_index(drop=True)


def make_particle_summary(merged: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for (species, pair), group in merged.groupby(["species", "pair"], sort=True):
        dx = group["x_opalx_m_1rank"] - group["x_opalx_m_4rank"]
        dy = group["y_opalx_m_1rank"] - group["y_opalx_m_4rank"]
        ds = group["s_opalx_m_1rank"] - group["s_opalx_m_4rank"]
        dpx = group["px_opalx_1rank"] - group["px_opalx_4rank"]
        rows.append(
            {
                "species": species,
                "pair": int(pair),
                "samples": int(len(group)),
                "rmse_dx_um": float(np.sqrt(np.mean(np.square(dx))) * 1.0e6),
                "max_abs_dx_um": float(np.max(np.abs(dx)) * 1.0e6),
                "rmse_dy_um": float(np.sqrt(np.mean(np.square(dy))) * 1.0e6),
                "max_abs_dy_um": float(np.max(np.abs(dy)) * 1.0e6),
                "rmse_ds_um": float(np.sqrt(np.mean(np.square(ds))) * 1.0e6),
                "max_abs_ds_um": float(np.max(np.abs(ds)) * 1.0e6),
                "max_abs_dpx": float(np.max(np.abs(dpx))),
            }
        )
    return pd.DataFrame(rows)


def make_first_kick_summary(merged: pd.DataFrame) -> pd.DataFrame:
    first = merged.loc[merged["reference_step"].eq(1)].copy()
    first["one_rank_to_cain"] = np.abs(first["px_opalx_1rank"] / first["px_cain_1rank"])
    first["four_rank_to_cain"] = np.abs(first["px_opalx_4rank"] / first["px_cain_1rank"])
    first["four_rank_creation_rank"] = first["particle_id_4rank"].astype(np.int64) % 4
    return first[
        [
            "species",
            "pair",
            "global_step_1rank",
            "px_cain_1rank",
            "px_opalx_1rank",
            "px_opalx_4rank",
            "one_rank_to_cain",
            "four_rank_to_cain",
            "four_rank_creation_rank",
        ]
    ].rename(
        columns={
            "global_step_1rank": "global_step",
            "px_cain_1rank": "px_cain",
        }
    )


def charge_symmetry(first: pd.DataFrame, momentum_column: str) -> pd.Series:
    pivot = first.pivot(index="pair", columns="species", values=momentum_column)
    electron = np.abs(pivot["electron"])
    positron = np.abs(pivot["positron"])
    return 200.0 * np.abs(electron - positron) / (electron + positron)


def plot_first_kicks(first: pd.DataFrame, output: Path) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), dpi=180)

    markers = {"electron": "o", "positron": "s"}
    for species in ("electron", "positron"):
        group = first.loc[first["species"].eq(species)].sort_values("pair")
        axes[0].plot(
            group["pair"],
            group["one_rank_to_cain"],
            color=COLORS[species],
            marker=markers[species],
            linewidth=1.8,
            label=f"1 rank {species}",
        )
        axes[0].plot(
            group["pair"],
            group["four_rank_to_cain"],
            color=COLORS[species],
            marker=markers[species],
            markerfacecolor="white",
            linestyle="--",
            linewidth=1.4,
            label=f"4 ranks {species}",
        )
    axes[0].axhline(1.0, color="0.2", linewidth=0.9, linestyle=(0, (2, 4)))
    axes[0].set_yscale("log")
    axes[0].set_xlabel("CAIN pair / birth order")
    axes[0].set_ylabel(r"$|\Delta p_x^{\mathrm{OPALX}}/\Delta p_x^{\mathrm{CAIN}}|$")
    axes[0].set_title("First transverse kick")
    axes[0].set_xticks(range(1, 7))
    axes[0].grid(True, which="both", color="0.9", linewidth=0.6)
    axes[0].legend(frameon=False, fontsize=9, ncol=2)

    one_symmetry = charge_symmetry(first, "px_opalx_1rank")
    four_symmetry = charge_symmetry(first, "px_opalx_4rank")
    width = 0.34
    pair = np.arange(1, 7)
    axes[1].bar(pair - width / 2, one_symmetry, width, color="#5b8f29", label="1 rank")
    axes[1].bar(pair + width / 2, four_symmetry, width, color="#7a4eab", label="4 ranks")
    axes[1].set_yscale("log")
    axes[1].set_ylim(1.0e-9, 100.0)
    axes[1].set_xlabel("CAIN pair / birth order")
    axes[1].set_ylabel("electron/positron kick asymmetry [%]")
    axes[1].set_title("Charge-sign symmetry")
    axes[1].set_xticks(pair)
    axes[1].grid(True, axis="y", which="both", color="0.9", linewidth=0.6)
    axes[1].legend(frameon=False)

    fig.suptitle("Timed track12 MPI-rank diagnostic", fontsize=16)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(output)
    plt.close(fig)


def plot_rank_differences(merged: pd.DataFrame, output: Path) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 11, "axes.labelsize": 13})
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), dpi=180, sharex=True)
    for axis, (coordinate, label) in zip(axes, (("x", "x"), ("y", "y")), strict=True):
        for (species, pair), group in merged.groupby(["species", "pair"], sort=True):
            difference_um = (
                group[f"{coordinate}_opalx_m_4rank"]
                - group[f"{coordinate}_opalx_m_1rank"]
            ) * 1.0e6
            axis.plot(
                group["s_opalx_m_1rank"] * 1.0e3,
                difference_um,
                color=COLORS[species],
                linestyle=LINE_STYLES[int(pair) - 1],
                linewidth=1.0,
            )
        axis.axhline(0.0, color="black", linewidth=0.7, linestyle=(0, (2, 6)))
        axis.set_xlabel("S - IP [mm]")
        axis.set_ylabel(rf"$\Delta {label}_{{4r-1r}}$ [$\mu$m]")
        axis.grid(True, color="0.9", linewidth=0.5)
        axis.set_title(f"Rank-count difference in {label}(s)")
    fig.legend(
        handles=[
            Line2D([0], [0], color=COLORS["electron"], label="electron"),
            Line2D([0], [0], color=COLORS["positron"], label="positron"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Four MPI ranks minus one MPI rank", fontsize=16)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-rank", type=Path, default=DEFAULT_ONE_RANK)
    parser.add_argument("--four-rank", type=Path, default=DEFAULT_FOUR_RANK)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    output_dir = args.output_dir or args.one_rank / "results" / "rank_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    merged = merge_runs(load_comparison(args.one_rank), load_comparison(args.four_rank))
    particle_summary = make_particle_summary(merged)
    first_kicks = make_first_kick_summary(merged)

    particle_summary.to_csv(output_dir / "track12_rank_count_particle_summary.csv", index=False)
    first_kicks.to_csv(output_dir / "track12_rank_count_first_kicks.csv", index=False)
    report = {
        "matched_samples": int(len(merged)),
        "one_rank_charge_symmetry_error_pct": charge_symmetry(
            first_kicks, "px_opalx_1rank"
        ).to_dict(),
        "four_rank_charge_symmetry_error_pct": charge_symmetry(
            first_kicks, "px_opalx_4rank"
        ).to_dict(),
        "particles": particle_summary.to_dict(orient="records"),
    }
    (output_dir / "track12_rank_count_comparison.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    plot_first_kicks(first_kicks, output_dir / "track12_rank_count_first_kicks.png")
    plot_rank_differences(merged, output_dir / "track12_rank_count_trajectory_difference.png")
    print(first_kicks.to_string(index=False))
    print(particle_summary.to_string(index=False))
    print(f"results: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
