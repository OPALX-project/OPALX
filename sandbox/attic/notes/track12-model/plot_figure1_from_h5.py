#!/usr/bin/env python3
"""Reproduce the note's Figure 1 layout using OPALX witness H5 output.

The left panel is parsed directly from TestParticleOrbit.dat.  The right panel
is read directly from OPALX H5 files; no CSV files are read or written.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


C_LIGHT = 299_792_458.0

SCRIPT_DIR = Path(__file__).resolve().parent
SANDBOX_DIR = SCRIPT_DIR.parents[1]
TRACK12_SCRIPT = SANDBOX_DIR / "track12particles" / "track12particles.py"
DEFAULT_REFERENCE = SANDBOX_DIR / "TestParticleOrbit.dat"
DEFAULT_RUN_DIR = (
    SANDBOX_DIR
    / "Drift-Experiment"
    / "withness_t0_4ps_active_1600fs_release_shifted_c0"
)

ELECTRON_COLOR = "#c00000"
POSITRON_COLOR = "#0017b8"
LINE_STYLES = [
    "-",
    (0, (8, 6)),
    (0, (5, 5)),
    (0, (10, 8)),
    (0, (2, 3)),
    (0, (9, 4, 2, 4)),
]


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure_matplotlib(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    mpl_dir = output.parent / ".matplotlib"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )


def sorted_steps(h5file: h5py.File) -> list[str]:
    return sorted(h5file.keys(), key=lambda key: int(key.split("#")[1]))


def reference_z(group: h5py.Group) -> float:
    if "RefPartR" in group.attrs:
        return float(group.attrs["RefPartR"][2])
    if "SPOS" in group.attrs:
        return float(group.attrs["SPOS"][0])
    return 0.0


def group_time_ps(group: h5py.Group, time_origin_s: float) -> float:
    if "TIME" not in group.attrs:
        return np.nan
    return 1.0e12 * (float(group.attrs["TIME"][0]) - time_origin_s)


def load_reference(path: Path):
    module = load_track12_module()
    frame = module.parse_reference_file(path).copy()
    if "s_mm" not in frame:
        frame["s_mm"] = 1.0e3 * frame["s"]
    if "x_um" not in frame:
        frame["x_um"] = 1.0e6 * frame["x"]
    if "t_ps" not in frame:
        frame["t_ps"] = 1.0e12 * frame["t"] / C_LIGHT
    return frame


def load_h5_trajectories(
    path: Path,
    *,
    kind: int,
    species: str,
    ip_s: float,
    time_origin_s: float,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    with h5py.File(path, "r") as h5:
        steps = sorted_steps(h5)
        for step_index, step_name in enumerate(steps):
            group = h5[step_name]
            ref_z = reference_z(group)
            time_ps = group_time_ps(group, time_origin_s)
            ids = group["id"][:] if "id" in group else np.arange(len(group["x"]))
            for index, particle_id in enumerate(ids):
                pair = int(particle_id) + 1
                rows.append(
                    {
                        "pair": pair,
                        "kind": kind,
                        "species": species,
                        "step": step_index,
                        "s_mm": 1.0e3 * (ref_z + float(group["z"][index]) - ip_s),
                        "t_ps": time_ps,
                        "x_um": 1.0e6 * float(group["x"][index]),
                    }
                )
    return rows


def x_column_for(axis: str) -> tuple[str, str]:
    if axis == "s":
        return "s_mm", "S (mm)"
    if axis == "t":
        return "t_ps", "t - T0 (ps)"
    raise ValueError(f"unsupported x axis {axis!r}")


def plot_groups(ax, rows, *, x_column: str) -> None:
    for pair in range(1, 7):
        style = LINE_STYLES[(pair - 1) % len(LINE_STYLES)]
        for kind, color in ((2, ELECTRON_COLOR), (3, POSITRON_COLOR)):
            group = [
                row for row in rows if int(row["pair"]) == pair and int(row["kind"]) == kind
            ]
            if not group:
                continue
            group = sorted(group, key=lambda row: float(row[x_column]))
            ax.plot(
                [float(row[x_column]) for row in group],
                [float(row["x_um"]) for row in group],
                color=color,
                linewidth=0.9,
                linestyle=style,
            )


def reference_rows(frame) -> list[dict[str, float | int | str]]:
    return [
        {
            "pair": int(row.pair),
            "kind": int(row.kind),
            "species": row.species,
            "s_mm": float(row.s_mm),
            "t_ps": float(row.t_ps),
            "x_um": float(row.x_um),
        }
        for row in frame.itertuples(index=False)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--electron-h5", type=Path)
    parser.add_argument("--positron-h5", type=Path)
    parser.add_argument("--ip-s", type=float, default=4.0e-3)
    parser.add_argument("--time-origin-s", type=float, default=4.0e-12)
    parser.add_argument("--x-axis", choices=("s", "t"), default="s")
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIR / "track12_figure1_opalx_h5.png",
    )
    parser.add_argument("--wide-x", action="store_true", help="Do not use Figure 1 S-axis limits.")
    parser.add_argument("--wide-y", action="store_true", help="Do not use Figure 1 y limits.")
    args = parser.parse_args()

    electron_h5 = args.electron_h5 or args.run_dir / "spacecharge_drift_withness_c1.h5"
    positron_h5 = args.positron_h5 or args.run_dir / "spacecharge_drift_withness_c2.h5"

    configure_matplotlib(args.output)
    x_column, xlabel = x_column_for(args.x_axis)
    reference = reference_rows(load_reference(args.reference))
    opalx = []
    opalx.extend(
        load_h5_trajectories(
            electron_h5,
            kind=2,
            species="electron",
            ip_s=args.ip_s,
            time_origin_s=args.time_origin_s,
        )
    )
    opalx.extend(
        load_h5_trajectories(
            positron_h5,
            kind=3,
            species="positron",
            ip_s=args.ip_s,
            time_origin_s=args.time_origin_s,
        )
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6), dpi=180, sharey=True)
    plot_groups(axes[0], reference, x_column=x_column)
    plot_groups(axes[1], opalx, x_column=x_column)

    axes[0].set_title("CAIN reference")
    axes[1].set_title("OPALX H5 BeamBeam")
    for ax in axes:
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle=(0, (2, 6)))
        ax.set_xlabel(xlabel)
        ax.tick_params(direction="out", top=True, right=True, which="both")
        ax.minorticks_on()
        if args.x_axis == "s" and not args.wide_x:
            ax.set_xlim(-0.75, 1.65)
        if not args.wide_y:
            ax.set_ylim(-9.0, 20.0)
    axes[0].set_ylabel(r"x ($\mu$m)")
    fig.suptitle("12 Test Particle Trajectories", fontsize=16)
    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)

    print(args.output.resolve())
    print(f"reference: {args.reference.resolve()}")
    print(f"electron H5: {electron_h5.resolve()}")
    print(f"positron H5: {positron_h5.resolve()}")


if __name__ == "__main__":
    main()
