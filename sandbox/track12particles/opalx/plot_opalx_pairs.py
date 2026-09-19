#!/usr/bin/env python3
"""Plot OPALX CAIN e-/e+ witness trajectories in the note style."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def sorted_steps(h5file):
    return sorted(h5file.keys(), key=lambda key: int(key.split("#")[1]))


def collect_trajectory(path: Path):
    with h5py.File(path, "r") as h5:
        steps = sorted_steps(h5)
        first = h5[steps[0]]
        particle_ids = [int(v) for v in first["id"][:]]
        trajectories = {pid: {"z": [], "x": [], "y": []} for pid in particle_ids}
        for step in steps:
            group = h5[step]
            ids = group["id"][:]
            x = group["x"][:]
            y = group["y"][:]
            z = group["z"][:]
            for i, pid_raw in enumerate(ids):
                pid = int(pid_raw)
                if pid not in trajectories:
                    trajectories[pid] = {"z": [], "x": [], "y": []}
                trajectories[pid]["z"].append(float(z[i]))
                trajectories[pid]["x"].append(float(x[i]))
                trajectories[pid]["y"].append(float(y[i]))
    return trajectories


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--electron-h5", type=Path, default=Path("track12_primary_with_pairs_smoke_c1.h5"))
    parser.add_argument("--positron-h5", type=Path, default=Path("track12_primary_with_pairs_smoke_c2.h5"))
    parser.add_argument(
        "--cain-csv",
        type=Path,
        default=Path("../outputs/track12_trajectory.csv"),
        help="Parsed CAIN reference trajectory CSV to overlay.",
    )
    parser.add_argument("--output", type=Path, default=Path("track12_opalx_pairs_x_vs_z.png"))
    args = parser.parse_args()

    electron = collect_trajectory(args.electron_h5)
    positron = collect_trajectory(args.positron_h5)

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    if args.cain_csv.exists():
        cain = pd.read_csv(args.cain_csv)
        for name, group in cain.groupby("name", sort=True):
            color = "tab:blue" if group["species"].iloc[0] == "electron" else "tab:red"
            ax.plot(
                group["s"] * 1e3,
                group["x"] * 1e6,
                color=color,
                lw=1.0,
                alpha=0.25,
                label="CAIN e-" if name == "Te-1" else ("CAIN e+" if name == "Te+1" else None),
            )

    for pid, traj in sorted(electron.items()):
        ax.plot(
            [z * 1e3 for z in traj["z"]],
            [x * 1e6 for x in traj["x"]],
            lw=1.8,
            color="tab:blue",
            alpha=0.85,
            ls="--",
            label="OPALX e-" if pid == min(electron) else None,
        )
    for pid, traj in sorted(positron.items()):
        ax.plot(
            [z * 1e3 for z in traj["z"]],
            [x * 1e6 for x in traj["x"]],
            lw=1.8,
            color="tab:red",
            alpha=0.85,
            ls="--",
            label="OPALX e+" if pid == min(positron) else None,
        )

    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [um]")
    ax.set_title("CAIN reference and OPALX e-/e+ read-in particles, no BeamBeam")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(args.output, dpi=220)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
