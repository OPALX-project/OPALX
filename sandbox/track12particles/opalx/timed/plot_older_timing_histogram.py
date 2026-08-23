#!/usr/bin/env python3
"""Plot timing distributions from selected nonempty OPALX ``timing.dat`` files."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DEFAULT_OUTPUT = HERE / "timing_histogram"
MAIN_TIMER = re.compile(
    r"^mainTimer\.+\s+(?P<ranks>\d+)\s+(?P<wall>[0-9.eE+-]+)\s*$"
)
DETAIL_TIMER = re.compile(
    r"^(?P<timer>.*?)\.+\s+(?P<ranks>\d+)\s+"
    r"(?P<wall_max>[0-9.eE+-]+)\s+(?P<wall_min>[0-9.eE+-]+)\s+"
    r"(?P<wall_avg>[0-9.eE+-]+)\s*$"
)


def default_inputs() -> list[Path]:
    patterns = (
        "sandbox/Drift-Experiment/*/timing.dat",
        "sandbox/track12particles/opalx/timing_pairs/*/timing.dat",
        "sandbox/cain-opalx-reduced-order-model/opalx/timing.dat",
    )
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(ROOT.glob(pattern))
    return sorted(path for path in paths if path.stat().st_size > 0)


def run_label(path: Path) -> str:
    relative = path.relative_to(ROOT)
    if "timing_pairs" in relative.parts:
        return f"timing_pairs/{path.parent.name}"
    if "cain-opalx-reduced-order-model" in relative.parts:
        return "cain-opalx-reduced-order-model"
    return path.parent.name


def parse_timing(path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    main_wall_s: float | None = None
    ranks: int | None = None
    details: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if main_wall_s is None and (match := MAIN_TIMER.match(line)):
            ranks = int(match.group("ranks"))
            main_wall_s = float(match.group("wall"))
            continue
        if match := DETAIL_TIMER.match(line):
            details.append(
                {
                    "run": run_label(path),
                    "path": str(path.relative_to(ROOT)),
                    "timer": match.group("timer").strip(),
                    "ranks": int(match.group("ranks")),
                    "wall_max_s": float(match.group("wall_max")),
                    "wall_min_s": float(match.group("wall_min")),
                    "wall_avg_s": float(match.group("wall_avg")),
                }
            )
    if main_wall_s is None or ranks is None:
        raise ValueError(f"{path}: mainTimer total was not found")
    for row in details:
        row["fraction_of_main"] = float(row["wall_avg_s"]) / main_wall_s
    return {
        "run": run_label(path),
        "path": str(path.relative_to(ROOT)),
        "ranks": ranks,
        "main_wall_s": main_wall_s,
    }, details


def logarithmic_bins(values: pd.Series, count: int) -> np.ndarray:
    positive = values.loc[values > 0.0].to_numpy(dtype=float)
    if positive.size == 0:
        raise ValueError("histogram requires positive values")
    lower = 10.0 ** np.floor(np.log10(positive.min()))
    upper = 10.0 ** np.ceil(np.log10(positive.max()))
    if np.isclose(lower, upper):
        upper = lower * 10.0
    return np.geomspace(lower, upper, count + 1)


def plot_histograms(runs: pd.DataFrame, timers: pd.DataFrame, output: Path) -> None:
    cache = output.parent / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), constrained_layout=True)

    main_bins = logarithmic_bins(runs["main_wall_s"], 8)
    axes[0].hist(
        runs["main_wall_s"], bins=main_bins, color="C0", edgecolor="white", linewidth=0.8
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("mainTimer wall time [s]")
    axes[0].set_ylabel("number of runs")
    axes[0].set_title(f"Selected completed runs ($n={len(runs)}$)")

    nonzero = timers.loc[(timers["timer"] != "mainTimer") & (timers["fraction_of_main"] > 0.0)].copy()
    nonzero["percent_of_main"] = 100.0 * nonzero["fraction_of_main"]
    fraction_bins = logarithmic_bins(nonzero["percent_of_main"], 16)
    axes[1].hist(
        nonzero["percent_of_main"],
        bins=fraction_bins,
        color="C1",
        edgecolor="white",
        linewidth=0.8,
    )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("reported timer / mainTimer [%]")
    axes[1].set_ylabel("number of timer entries")
    axes[1].set_title(f"Nonzero named timers ($n={len(nonzero)}$; nested)")

    for axis in axes:
        axis.grid(True, which="both", axis="y", color="0.88", linewidth=0.6)
        axis.set_axisbelow(True)
    fig.suptitle("OPALX timing histograms")
    fig.savefig(output, dpi=240)
    plt.close(fig)


def plot_longest_run_top_timers(
    runs: pd.DataFrame, timers: pd.DataFrame, output: Path
) -> pd.DataFrame:
    import matplotlib.pyplot as plt

    slowest = runs.loc[runs["main_wall_s"].idxmax()]
    top = (
        timers.loc[
            (timers["run"] == slowest["run"]) & (timers["timer"] != "mainTimer")
        ]
        .nlargest(10, "wall_avg_s")
        .copy()
    )
    top["percent_of_main"] = 100.0 * top["fraction_of_main"]
    ordered = top.sort_values("wall_avg_s")

    fig, axis = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    bars = axis.barh(ordered["timer"], ordered["wall_avg_s"], color="C0")
    axis.bar_label(
        bars,
        labels=[
            f"{seconds:.1f} s  ({percent:.1f}%)"
            for seconds, percent in zip(
                ordered["wall_avg_s"], ordered["percent_of_main"]
            )
        ],
        padding=4,
        fontsize=8,
    )
    axis.set_xlabel("Wall avg [s]")
    axis.set_ylabel("reported timer (nested)")
    axis.set_title(
        "Ten longest timers in the longest selected run\n"
        f"{slowest['run']} — mainTimer {slowest['main_wall_s']:.0f} s",
        fontsize=11,
    )
    axis.set_xlim(0.0, 1.22 * ordered["wall_avg_s"].max())
    axis.grid(True, axis="x", color="0.88", linewidth=0.6)
    axis.set_axisbelow(True)
    fig.savefig(output, dpi=240)
    plt.close(fig)
    return top


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    inputs = args.inputs or default_inputs()
    if not inputs:
        raise FileNotFoundError("no nonempty timing.dat files found")
    run_rows: list[dict[str, object]] = []
    timer_rows: list[dict[str, object]] = []
    for path in inputs:
        run, details = parse_timing(path.resolve())
        run_rows.append(run)
        timer_rows.extend(details)

    runs = pd.DataFrame(run_rows).sort_values("main_wall_s").reset_index(drop=True)
    timers = pd.DataFrame(timer_rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    runs.to_csv(args.output_dir / "older_run_main_timers.csv", index=False)
    timers.to_csv(args.output_dir / "older_run_timer_rows.csv", index=False)
    plot_histograms(runs, timers, args.output_dir / "older_run_timing_histograms.png")
    top = plot_longest_run_top_timers(
        runs, timers, args.output_dir / "longest_run_top10_timers.png"
    )
    top.to_csv(args.output_dir / "longest_run_top10_timers.csv", index=False)
    print(runs.to_string(index=False))
    print("\nTen longest timers in the longest run:")
    print(top[["timer", "wall_avg_s", "percent_of_main"]].to_string(index=False))
    print(f"parsed {len(timers)} detailed timer rows from {len(runs)} runs")
    print(f"wrote {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
