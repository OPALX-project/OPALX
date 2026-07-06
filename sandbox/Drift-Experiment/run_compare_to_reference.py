#!/usr/bin/env python3
"""Run drift convergence cases and compare metrics against saved reference data."""

from __future__ import annotations

import argparse
import importlib.util
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CONVERGENCE = HERE / "run_spacecharge_convergence.py"
DEFAULT_OPALX = REPO / "build_openmp" / "src" / "opalx"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=HERE / "reference" / "full27")
    parser.add_argument("--output-dir", type=Path, default=HERE / "post_merge_compare")
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--n-grid-values", default="16,32,64")
    parser.add_argument("--nppg-values", default="5,8,14")
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def load_convergence_module():
    spec = importlib.util.spec_from_file_location("run_spacecharge_convergence", CONVERGENCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {CONVERGENCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_matrix(args: argparse.Namespace) -> None:
    cmd = [
        str(Path.home() / ".venv-h6" / "bin" / "python"),
        str(CONVERGENCE),
        "--output-dir",
        str(args.output_dir / "current"),
        "--opalx",
        str(args.opalx),
        "--threads",
        str(args.threads),
        "--n-grid-values",
        args.n_grid_values,
        "--nppg-values",
        args.nppg_values,
        "--seeds",
        args.seeds,
    ]
    if args.force:
        cmd.append("--force")
    if args.resume:
        cmd.append("--resume")
    subprocess.run(cmd, cwd=REPO, check=True)


def compare_summary(reference_dir: Path, current_dir: Path, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref = pd.read_csv(reference_dir / "convergence_summary.csv")
    cur = pd.read_csv(current_dir / "convergence_summary.csv")
    keys = ["case", "n_grid", "nppg", "seed", "configured_particles"]
    merged = ref.merge(cur, on=keys, suffixes=("_ref", "_cur"), validate="one_to_one")
    metrics = [
        "relative_vector_L2_vs_analytic",
        "mean_pointwise_relative_vector_error",
        "p95_pointwise_relative_vector_error",
        "relative_magnitude_L2_vs_analytic",
        "mean_abs_pic_E",
        "mean_abs_analytic_E",
        "ratio_to_3_over_sqrt_N",
    ]
    for metric in metrics:
        merged[f"{metric}_diff"] = merged[f"{metric}_cur"] - merged[f"{metric}_ref"]
        denom = merged[f"{metric}_ref"].replace(0.0, np.nan)
        merged[f"{metric}_rel_diff"] = merged[f"{metric}_diff"] / denom
    merged.to_csv(output_dir / "case_metric_differences.csv", index=False)

    diff_cols = [f"{metric}_diff" for metric in metrics]
    rel_cols = [f"{metric}_rel_diff" for metric in metrics]
    grouped = merged.groupby(["n_grid", "nppg"], as_index=False).agg(
        runs=("case", "count"),
        relative_l2_diff_mean=("relative_vector_L2_vs_analytic_diff", "mean"),
        relative_l2_diff_std=("relative_vector_L2_vs_analytic_diff", "std"),
        relative_l2_rel_diff_mean=("relative_vector_L2_vs_analytic_rel_diff", "mean"),
        relative_l2_rel_diff_std=("relative_vector_L2_vs_analytic_rel_diff", "std"),
    )
    for col in diff_cols + rel_cols:
        grouped[f"{col}_max_abs"] = merged.groupby(["n_grid", "nppg"])[col].apply(
            lambda s: float(np.max(np.abs(s)))
        ).to_numpy()
    grouped.to_csv(output_dir / "grouped_metric_differences.csv", index=False)
    return merged, grouped


def save_difference_plots(diffs: pd.DataFrame, grouped: pd.DataFrame, output_dir: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output_dir / ".cache"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    for nppg, group in grouped.groupby("nppg"):
        ax.errorbar(
            group["n_grid"],
            group["relative_l2_diff_mean"],
            yerr=group["relative_l2_diff_std"].fillna(0.0),
            marker="o",
            capsize=3,
            label=f"NPPG={nppg}",
        )
    ax.axhline(0.0, color="0.2", lw=1.0, ls="--")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("N_GRID")
    ax.set_ylabel("current - reference relative L2")
    ax.grid(True, color="0.88", lw=0.6)
    ax.legend()
    fig.savefig(output_dir / "difference_relative_l2_vs_grid.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    for nppg, group in grouped.groupby("nppg"):
        ax.errorbar(
            group["n_grid"],
            100.0 * group["relative_l2_rel_diff_mean"],
            yerr=100.0 * group["relative_l2_rel_diff_std"].fillna(0.0),
            marker="o",
            capsize=3,
            label=f"NPPG={nppg}",
        )
    ax.axhline(0.0, color="0.2", lw=1.0, ls="--")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("N_GRID")
    ax.set_ylabel("relative L2 change [%]")
    ax.grid(True, color="0.88", lw=0.6)
    ax.legend()
    fig.savefig(output_dir / "difference_relative_l2_percent_vs_grid.png", dpi=220)
    plt.close(fig)

    heat = grouped.pivot(index="nppg", columns="n_grid", values="relative_l2_diff_mean")
    fig, ax = plt.subplots(figsize=(6.6, 4.8), constrained_layout=True)
    vmax = float(np.nanmax(np.abs(heat.to_numpy()))) if not heat.empty else 1.0
    vmax = vmax if vmax > 0 else 1.0
    im = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(heat.columns)), labels=[str(v) for v in heat.columns])
    ax.set_yticks(np.arange(len(heat.index)), labels=[str(v) for v in heat.index])
    ax.set_xlabel("N_GRID")
    ax.set_ylabel("NPPG")
    ax.set_title("Mean relative L2 difference")
    fig.colorbar(im, ax=ax, label="current - reference")
    fig.savefig(output_dir / "difference_relative_l2_heatmap.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    ordered = diffs.sort_values(["n_grid", "nppg", "seed"])
    ax.bar(np.arange(len(ordered)), ordered["relative_vector_L2_vs_analytic_diff"], width=0.8)
    ax.axhline(0.0, color="0.2", lw=1.0)
    ax.set_xlabel("case index")
    ax.set_ylabel("current - reference relative L2")
    ax.grid(True, axis="y", color="0.88", lw=0.6)
    fig.savefig(output_dir / "difference_relative_l2_by_case.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.output_dir.exists() and args.force:
        import shutil

        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_matrix(args)
    current_dir = args.output_dir / "current"
    diffs, grouped = compare_summary(args.reference_dir, current_dir, args.output_dir)
    save_difference_plots(diffs, grouped, args.output_dir)
    (args.output_dir / "README.md").write_text(
        "\n".join(
            [
                "# Drift convergence reference comparison",
                "",
                f"Reference: `{args.reference_dir}`",
                f"Current: `{current_dir}`",
                "",
                "Difference convention: `current - reference`.",
                "",
                "Files:",
                "- `case_metric_differences.csv`",
                "- `grouped_metric_differences.csv`",
                "- `difference_relative_l2_vs_grid.png`",
                "- `difference_relative_l2_percent_vs_grid.png`",
                "- `difference_relative_l2_heatmap.png`",
                "- `difference_relative_l2_by_case.png`",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
