#!/usr/bin/env python3
"""Visualize the reference trajectory from the ISIS ten-turn RING run."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm


CIRCUMFERENCE_M = 163.3630220240698


def read_design_path(path: Path) -> pd.DataFrame:
    rows: list[list[float]] = []
    with path.open() as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.split(maxsplit=15)
            rows.append([float(value) for value in fields[:15]])

    columns = [
        "s",
        "rx",
        "ry",
        "rz",
        "px",
        "py",
        "pz",
        "ex",
        "ey",
        "ez",
        "bx",
        "by",
        "bz",
        "energy",
        "time",
    ]
    return pd.DataFrame(rows, columns=columns)


def read_lattice_positions(path: Path) -> pd.DataFrame:
    rows: list[tuple[float, float]] = []
    with path.open() as stream:
        for line in stream:
            if not line.startswith(('"BEGIN:', '"END:', '"ENTRY EDGE:', '"EXIT EDGE:')):
                continue
            fields = line.rsplit('"', maxsplit=1)[-1].split()
            if len(fields) >= 2:
                rows.append((float(fields[0]), float(fields[1])))
    return pd.DataFrame(rows, columns=["x", "z"])


def read_statistics(path: Path) -> pd.DataFrame:
    rows: list[list[float]] = []
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 47:
                continue
            try:
                rows.append(
                    [
                        float(fields[0]),
                        float(fields[1]),
                        float(fields[4]),
                        float(fields[17]),
                        float(fields[18]),
                        float(fields[19]),
                    ]
                )
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["time", "s", "energy", "rx", "ry", "rz"])


def nearest_turn_samples(frame: pd.DataFrame, turns: int) -> pd.DataFrame:
    samples = []
    for turn in range(1, turns + 1):
        target = turn * CIRCUMFERENCE_M
        index = (frame["s"] - target).abs().idxmin()
        sample = frame.loc[index].copy()
        sample["turn"] = turn
        samples.append(sample)
    result = pd.DataFrame(samples)
    x0, z0 = frame.iloc[0][["rx", "rz"]]
    result["closure_error"] = np.hypot(result["rx"] - x0, result["rz"] - z0)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    base = Path(__file__).resolve().parent
    parser.add_argument(
        "--design-path",
        type=Path,
        default=base / "data" / "isis_sbend_ring_10t_DesignPath.dat",
    )
    parser.add_argument(
        "--element-positions",
        type=Path,
        default=base / "data" / "isis_sbend_ring_10t_ElementPositions.txt",
    )
    parser.add_argument(
        "--statistics",
        type=Path,
        default=base / "isis_sbend_ring_10t.stat",
    )
    parser.add_argument("--output", type=Path, default=base / "isis_sbend_ring_10t.png")
    args = parser.parse_args()

    frame = read_design_path(args.design_path)
    lattice = read_lattice_positions(args.element_positions)
    statistics = read_statistics(args.statistics)
    turns = 10
    samples = nearest_turn_samples(statistics, turns)

    # Downsample only for rendering; turn metrics use the statistics output.
    stride = max(1, len(frame) // 4_000)
    plotted = frame.iloc[::stride].copy()
    statistics["turn"] = np.minimum(
        (statistics["s"] / CIRCUMFERENCE_M).astype(int) + 1, turns
    )

    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(12.5, 6.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.65, 1.0])
    orbit_ax = fig.add_subplot(grid[:, 0])
    closure_ax = fig.add_subplot(grid[0, 1])
    energy_ax = fig.add_subplot(grid[1, 1])

    orbit_ax.scatter(
        lattice["x"], lattice["z"], s=5, color="0.75", alpha=0.6, label="lattice elements"
    )
    orbit_ax.plot(
        plotted["rx"],
        plotted["rz"],
        linewidth=1.5,
        color="0.35",
        label="one-turn threaded map",
    )
    points = statistics[["rx", "rz"]].to_numpy().reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = BoundaryNorm(np.arange(0.5, turns + 1.5), turns)
    collection = LineCollection(
        segments,
        cmap="viridis",
        norm=norm,
        linewidth=1.2,
        alpha=0.9,
    )
    collection.set_array(statistics["turn"].to_numpy()[:-1])
    orbit_ax.add_collection(collection)
    orbit_ax.autoscale()
    orbit_ax.set_aspect("equal", adjustable="datalim")
    orbit_ax.set_xlabel("Reference x [m]")
    orbit_ax.set_ylabel("Reference z [m]")
    orbit_ax.set_title("One-turn map and tracked reference samples")
    orbit_ax.legend(loc="best", frameon=False)
    colorbar = fig.colorbar(collection, ax=orbit_ax, ticks=np.arange(1, turns + 1))
    colorbar.set_label("Nominal turn")

    closure_ax.semilogy(
        samples["turn"], samples["closure_error"], marker="o", linewidth=1.8
    )
    closure_ax.set_xlabel("Nominal turn")
    closure_ax.set_ylabel("Distance from start [m]")
    closure_ax.set_title("Closure error at $s=nC$")
    closure_ax.set_xticks(np.arange(1, turns + 1))

    energy_ax.plot(samples["turn"], samples["energy"], marker="o", linewidth=1.8)
    energy_ax.set_xlabel("Nominal turn")
    energy_ax.set_ylabel("Reference kinetic energy [MeV]")
    energy_ax.set_title("Energy remains stable")
    energy_ax.set_xticks(np.arange(1, turns + 1))
    energy_ax.ticklabel_format(axis="y", style="plain", useOffset=False)

    fig.suptitle(
        "ISIS RING: 10-turn reference trajectory (without closed-orbit matching)",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(args.output, dpi=220)

    print(samples[["turn", "s", "rx", "rz", "closure_error", "energy"]].to_string(index=False))
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
