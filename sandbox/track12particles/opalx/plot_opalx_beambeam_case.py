#!/usr/bin/env python3
"""Plot OPALX e-/e+ witness trajectories from BeamBeam H5 output."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


C_LIGHT = 299_792_458.0
WITNESS_BETA = 3.0**0.5 / 2.0

LINE_STYLES = {
    1: "-",
    2: (0, (8, 5)),
    3: (0, (5, 5)),
    4: (0, (9, 6)),
    5: (0, (2, 3)),
    6: (0, (10, 4, 2, 4)),
}

ELECTRON_COLOR = "#c00000"
POSITRON_COLOR = "#0017b8"


def sorted_steps(h5file: h5py.File) -> list[str]:
    return sorted(h5file.keys(), key=lambda key: int(key.split("#")[1]))


def reference_z(group: h5py.Group) -> float:
    if "RefPartR" in group.attrs:
        return float(group.attrs["RefPartR"][2])
    if "SPOS" in group.attrs:
        return float(group.attrs["SPOS"][0])
    return 0.0


def group_time(group: h5py.Group) -> float | None:
    if "TIME" in group.attrs:
        return float(group.attrs["TIME"][0])
    if "T" in group.attrs:
        return float(group.attrs["T"][0])
    return None


def new_trajectory() -> dict[str, list[float]]:
    return {
        "step": [],
        "time_s": [],
        "s": [],
        "x": [],
        "y": [],
        "z": [],
        "z_local": [],
        "px": [],
        "py": [],
        "pz": [],
    }


def collect_trajectories(
    path: Path,
    ip_s: float,
    track_by: str,
) -> dict[int, dict[str, list[float]]]:
    with h5py.File(path, "r") as h5:
        steps = sorted_steps(h5)
        first = h5[steps[0]]
        if track_by == "id":
            trajectories = {int(pid): new_trajectory() for pid in first["id"][:]}
        else:
            trajectories = {index: new_trajectory() for index in range(len(first["id"]))}

        for step_index, step in enumerate(steps):
            group = h5[step]
            ref_z = reference_z(group)
            time_s = group_time(group)
            ids = group["id"][:]
            if track_by == "id":
                particle_order = [(int(pid), index) for index, pid in enumerate(ids)]
            else:
                sorted_indices = sorted(range(len(ids)), key=lambda index: float(group["z"][index]))
                particle_order = list(enumerate(sorted_indices))

            for pid, index in particle_order:
                trajectories.setdefault(pid, new_trajectory())
                z_local = float(group["z"][index])
                z_absolute = ref_z + z_local
                trajectories[pid]["step"].append(step_index)
                trajectories[pid]["time_s"].append(time_s)
                trajectories[pid]["s"].append(z_absolute - ip_s)
                trajectories[pid]["z"].append(z_absolute)
                trajectories[pid]["z_local"].append(z_local)
                for name in ("x", "y", "z", "px", "py", "pz"):
                    if name != "z":
                        trajectories[pid][name].append(float(group[name][index]))

    return trajectories


def trajectories_to_frame(
    trajectories: dict[int, dict[str, list[float]]],
    species: str,
    kind: int,
    model: str,
) -> pd.DataFrame:
    rows = []
    for pair, (_pid, trajectory) in enumerate(sorted(trajectories.items()), start=1):
        first_time_s = next((time for time in trajectory["time_s"] if time is not None), None)
        first_s = trajectory["s"][0] if trajectory["s"] else 0.0
        initial_ct_mm = None
        if first_time_s is not None:
            initial_ct_mm = 1.0e3 * (first_s - WITNESS_BETA * C_LIGHT * first_time_s)
        for step, time_s, s, x, y in zip(
            trajectory["step"],
            trajectory["time_s"],
            trajectory["s"],
            trajectory["x"],
            trajectory["y"],
            strict=True,
        ):
            elapsed_time_fs = None
            ct_mm = None
            if time_s is not None and first_time_s is not None:
                elapsed_time_fs = 1.0e15 * (time_s - first_time_s)
            if time_s is not None and initial_ct_mm is not None:
                ct_mm = initial_ct_mm + 1.0e3 * C_LIGHT * time_s
            rows.append(
                {
                    "model": model,
                    "pair": pair,
                    "kind": kind,
                    "species": species,
                    "step": step,
                    "time_s": time_s,
                    "elapsed_time_fs": elapsed_time_fs,
                    "ct_mm": ct_mm,
                    "s_mm": 1.0e3 * s,
                    "x_um": 1.0e6 * x,
                    "y_um": 1.0e6 * y,
                }
            )
    return pd.DataFrame(rows)


def load_reference_frame(path: Path, model: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"pair", "kind", "species", "s_mm", "x_um", "y_um"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")
    columns = ["pair", "kind", "species", "s_mm", "x_um", "y_um"]
    if "step" in frame.columns:
        columns.append("step")
    if "time_s" in frame.columns:
        columns.append("time_s")
    elif "t" in frame.columns:
        columns.append("t")
    out = frame.loc[:, columns].copy()
    if "time_s" not in out.columns and "t" in out.columns:
        out["ct_mm"] = 1.0e3 * out["t"]
        out["time_s"] = out["t"] / 299_792_458.0
        out = out.drop(columns=["t"])
    elif "time_s" in out.columns:
        out["ct_mm"] = 1.0e3 * C_LIGHT * out["time_s"]
    if "time_s" in out.columns:
        out["elapsed_time_fs"] = 1.0e15 * (
            out["time_s"] - out.groupby(["species", "pair"])["time_s"].transform("min")
        )
    out["model"] = model
    return out


def filter_by_elapsed_time(frame: pd.DataFrame, max_elapsed_time_fs: float | None) -> pd.DataFrame:
    if max_elapsed_time_fs is None:
        return frame
    if "elapsed_time_fs" not in frame.columns:
        raise ValueError(
            f"{frame['model'].iloc[0]} frame does not contain elapsed_time_fs for time filtering"
        )
    return frame.loc[frame["elapsed_time_fs"] <= max_elapsed_time_fs].copy()


def summarize(label: str, trajectories: dict[int, dict[str, list[float]]]) -> None:
    for coordinate, scale, unit in (
        ("s", 1.0e3, "mm"),
        ("x", 1.0e6, "um"),
        ("y", 1.0e6, "um"),
        ("z", 1.0e3, "mm"),
        ("z_local", 1.0e3, "mm"),
        ("px", 1.0, ""),
        ("py", 1.0, ""),
        ("pz", 1.0, ""),
    ):
        values = [value for traj in trajectories.values() for value in traj[coordinate]]
        suffix = f" {unit}" if unit else ""
        print(
            f"{label} {coordinate}: "
            f"{min(values) * scale:.12g} to {max(values) * scale:.12g}{suffix}"
        )


def coordinate_axes(coordinate: str) -> tuple[str, str, str, str]:
    if coordinate == "x":
        return "s_mm", "x_um", "S - S_IP [mm]", "x [um]"
    if coordinate == "y":
        return "s_mm", "y_um", "S - S_IP [mm]", "y [um]"
    if coordinate == "z":
        return "ct_mm", "s_mm", "ct [mm]", "z = S - S_IP [mm]"
    raise ValueError(f"unsupported coordinate {coordinate!r}")


def plot_coordinate_vs_s(
    electrons: dict[int, dict[str, list[float]]],
    positrons: dict[int, dict[str, list[float]]],
    output: Path,
    title: str,
    coordinate: str,
) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 5.0), constrained_layout=True)
    _x_column, _y_column, xlabel, ylabel = coordinate_axes(coordinate)

    for pair, (_pid, trajectory) in enumerate(sorted(electrons.items()), start=1):
        ax.plot(
            [s * 1.0e3 for s in trajectory["s"]],
            [value * 1.0e6 for value in trajectory[coordinate]],
            color=ELECTRON_COLOR,
            ls=LINE_STYLES.get(pair, "-"),
            lw=1.7,
            label=f"e- {pair}",
        )

    for pair, (_pid, trajectory) in enumerate(sorted(positrons.items()), start=1):
        ax.plot(
            [s * 1.0e3 for s in trajectory["s"]],
            [value * 1.0e6 for value in trajectory[coordinate]],
            color=POSITRON_COLOR,
            ls=LINE_STYLES.get(pair, "-"),
            lw=1.7,
            label=f"e+ {pair}",
        )

    ax.axhline(0.0, color="black", lw=1.0, ls=(0, (1.5, 6)), alpha=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3, fontsize=8)
    fig.savefig(output, dpi=220)


def plot_model_overlay(
    frames: list[pd.DataFrame],
    output: Path,
    title: str,
    coordinate: str,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
) -> None:
    model_styles = {
        "OPALX": "-",
        "CAIN": "--",
        "manufactured anisotropic Gaussian": "-.",
    }
    model_linewidth = {
        "OPALX": 1.6,
        "CAIN": 1.2,
        "manufactured anisotropic Gaussian": 1.2,
    }
    species_colors = {
        "electron": ELECTRON_COLOR,
        "positron": POSITRON_COLOR,
    }

    data = pd.concat(frames, ignore_index=True)
    x_column, y_column, xlabel, ylabel = coordinate_axes(coordinate)
    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)

    for (model, species, pair), group in data.groupby(["model", "species", "pair"], sort=False):
        ax.plot(
            group[x_column],
            group[y_column],
            color=species_colors.get(str(species), "black"),
            ls=model_styles.get(str(model), "-"),
            lw=model_linewidth.get(str(model), 1.2),
            alpha=0.9,
        )

    ax.axhline(0.0, color="black", lw=1.0, ls=(0, (1.5, 6)), alpha=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    species_handles = [
        plt.Line2D([0], [0], color=ELECTRON_COLOR, lw=2.0, label="e-"),
        plt.Line2D([0], [0], color=POSITRON_COLOR, lw=2.0, label="e+"),
    ]
    model_handles = [
        plt.Line2D(
            [0],
            [0],
            color="0.2",
            lw=model_linewidth.get(model, 1.2),
            ls=model_styles.get(model, "-"),
            label=model,
        )
        for model in model_styles
        if model in set(data["model"])
    ]
    first_legend = ax.legend(handles=species_handles, loc="upper left", fontsize=8, title="species")
    ax.add_artist(first_legend)
    ax.legend(handles=model_handles, loc="upper right", fontsize=8, title="model")

    fig.savefig(output, dpi=220)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prefix",
        type=Path,
        default=Path("track12_beambeam_copied_source_1fs_tiny"),
        help="Output prefix before _c1.h5 and _c2.h5.",
    )
    parser.add_argument(
        "--ip-s",
        type=float,
        default=1.0e-3,
        help="BeamBeam IP longitudinal position in meters.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("track12_opalx_copied_source_1fs_tiny_x_vs_s.png"),
    )
    parser.add_argument(
        "--title",
        default="OPALX copied-source BeamBeam tracking, DT = 1 fs",
    )
    parser.add_argument(
        "--track-by",
        choices=("id", "z-order"),
        default="id",
        help="Track particles by H5 id or by longitudinal order at each step.",
    )
    parser.add_argument(
        "--coordinate",
        choices=("x", "y", "z"),
        default="x",
        help="Coordinate to plot. x/y are plotted versus S-S_IP; z is plotted versus ct.",
    )
    parser.add_argument(
        "--opalx-species-assignment",
        choices=("cain-branch", "h5-charge"),
        default="cain-branch",
        help=(
            "Assign OPALX c1/c2 to plotted species by CAIN trajectory branch "
            "(c2=e-, c1=e+) or by H5 charge sign (c1=e-, c2=e+)."
        ),
    )
    parser.add_argument(
        "--cain-csv",
        type=Path,
        help="Optional CAIN trajectory CSV with s_mm and x_um columns.",
    )
    parser.add_argument(
        "--manufactured-csv",
        type=Path,
        help="Optional manufactured anisotropic Gaussian trajectory CSV with s_mm and x_um columns.",
    )
    parser.add_argument(
        "--xlim",
        nargs=2,
        type=float,
        metavar=("MIN_MM", "MAX_MM"),
        help="Optional x-axis limits in mm for overlay plots.",
    )
    parser.add_argument(
        "--ylim",
        nargs=2,
        type=float,
        metavar=("MIN_UM", "MAX_UM"),
        help="Optional y-axis limits in um for overlay plots.",
    )
    parser.add_argument(
        "--max-elapsed-time-fs",
        type=float,
        help="Optional per-particle elapsed-time cutoff for overlay plots.",
    )
    args = parser.parse_args()

    c1 = collect_trajectories(args.prefix.with_name(args.prefix.name + "_c1.h5"), args.ip_s, args.track_by)
    c2 = collect_trajectories(args.prefix.with_name(args.prefix.name + "_c2.h5"), args.ip_s, args.track_by)

    if args.opalx_species_assignment == "h5-charge":
        electrons = c1
        positrons = c2
    else:
        # The H5 charges are c1=e- and c2=e+, but the copied-source BeamBeam
        # field sign makes c2 follow the CAIN/manufactured electron branch.
        electrons = c2
        positrons = c1
    if args.cain_csv or args.manufactured_csv:
        opalx_frames = [
            trajectories_to_frame(electrons, "electron", 2, "OPALX"),
            trajectories_to_frame(positrons, "positron", 3, "OPALX"),
        ]
        overlay_frames = [filter_by_elapsed_time(frame, args.max_elapsed_time_fs) for frame in opalx_frames]
        if args.cain_csv:
            overlay_frames.append(
                filter_by_elapsed_time(load_reference_frame(args.cain_csv, "CAIN"), args.max_elapsed_time_fs)
            )
        if args.manufactured_csv:
            overlay_frames.append(
                filter_by_elapsed_time(
                    load_reference_frame(
                        args.manufactured_csv,
                        "manufactured anisotropic Gaussian",
                    ),
                    args.max_elapsed_time_fs,
                )
            )
        xlim = tuple(args.xlim) if args.xlim is not None else None
        ylim = tuple(args.ylim) if args.ylim is not None else None
        plot_model_overlay(
            overlay_frames,
            args.output,
            args.title,
            args.coordinate,
            xlim=xlim,
            ylim=ylim,
        )
    else:
        plot_coordinate_vs_s(electrons, positrons, args.output, args.title, args.coordinate)

    print(args.output.resolve())
    summarize("e-", electrons)
    summarize("e+", positrons)


if __name__ == "__main__":
    main()
