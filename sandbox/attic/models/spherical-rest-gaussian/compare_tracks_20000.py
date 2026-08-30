#!/usr/bin/env python3
"""Compare the 20k-step OPALX witness tracks with manufactured tracking."""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_CAIN_CSV = ROOT / "sandbox" / "track12particles" / "outputs" / "track12_trajectory.csv"
DEFAULT_DATA_DIR = SCRIPT_DIR / "data"
DEFAULT_FIGURE_DIR = SCRIPT_DIR / "figures"

IP_S_M = 0.003
DT_MAX_S = 5.0e-15
SIGMA_Z_FITTED_M = 0.274e-3
CHARGE_SCALE_FITTED = 0.7750192814260582
MC_SOURCE_PARTICLES = 400_000
MC_SEED = 20260629
SOURCE_RETIRE_TIME_S = 10.0e-12

ELECTRON_COLOR = "#c00000"
POSITRON_COLOR = "#0017b8"
LINE_STYLES = {
    1: "-",
    2: (0, (8, 5)),
    3: (0, (5, 5)),
    4: (0, (9, 6)),
    5: (0, (2, 3)),
    6: (0, (10, 4, 2, 4)),
}


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def configure_matplotlib() -> None:
    cache_dir = SCRIPT_DIR / "atic" / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def sorted_step_names(h5: h5py.File) -> list[str]:
    return sorted(h5.keys(), key=lambda key: int(key.split("#")[1]))


def load_opalx_tracks(prefix: Path) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for container, species, kind in ((1, "electron", 2), (2, "positron", 3)):
        path = prefix.with_name(f"{prefix.name}_c{container}.h5")
        with h5py.File(path, "r") as h5:
            steps = sorted_step_names(h5)
            first = h5[steps[0]]
            z0 = float(first.attrs["RefPartR"][2]) + first["z"][:]
            ids0 = first["id"][:]
            pair_by_id = {int(ids0[index]): pair for pair, index in enumerate(np.argsort(z0), 1)}
            for step_name in steps:
                group = h5[step_name]
                ref = group.attrs["RefPartR"][:]
                time_s = float(group.attrs["TIME"][0])
                global_step = int(group.attrs["GlobalTrackStep"][0])
                for index, particle_id_raw in enumerate(group["id"][:]):
                    particle_id = int(particle_id_raw)
                    pair = pair_by_id[particle_id]
                    x_m = float(ref[0] + group["x"][index])
                    y_m = float(ref[1] + group["y"][index])
                    s_m = float(ref[2] + group["z"][index] - IP_S_M)
                    rows.append(
                        {
                            "model": "OPALX",
                            "container": container,
                            "id": particle_id,
                            "pair": pair,
                            "kind": kind,
                            "species": species,
                            "global_step": global_step,
                            "time_s": time_s,
                            "time_ps": 1.0e12 * time_s,
                            "x_m": x_m,
                            "y_m": y_m,
                            "s_m": s_m,
                            "x_um": 1.0e6 * x_m,
                            "s_mm": 1.0e3 * s_m,
                        }
                    )
    return pd.DataFrame(rows)


def initial_conditions(reference: pd.DataFrame) -> pd.DataFrame:
    return (
        reference.sort_values(["pair", "kind", "step"])
        .groupby(["pair", "species"], as_index=False)
        .first()
        .sort_values(["pair", "kind"])
    )


def make_sources(track12, source_selection: str, source_charge_sign: float):
    return track12.make_anisotropic_sources(
        track12.DEFAULT_SOURCE_KINETIC_MEV,
        source_charge_sign * abs(track12.DEFAULT_SOURCE_CHARGE_C) * CHARGE_SCALE_FITTED,
        track12.SIGMA_XY_M,
        track12.SIGMA_XY_M,
        SIGMA_Z_FITTED_M,
        source_selection,
        0.0,
        MC_SOURCE_PARTICLES,
        MC_SEED,
    )


def run_manufactured_on_opalx_grid(
    opalx: pd.DataFrame,
    reference: pd.DataFrame,
    track12,
    source_selection: str,
    source_charge_sign: float,
    clock_mode: str,
) -> pd.DataFrame:
    sources = make_sources(track12, source_selection, source_charge_sign)
    initials = initial_conditions(reference)
    rows: list[dict[str, float | int | str]] = []
    for _, first in initials.iterrows():
        pair = int(first["pair"])
        species = str(first["species"])
        charge_units = float(first["charge_units"])
        position = first[["x", "y", "s"]].to_numpy(dtype=float)
        momentum = first[["Px_beta_gamma", "Py_beta_gamma", "Ps_beta_gamma"]].to_numpy(dtype=float)
        initial_clock_s = float(first["t"]) / track12.C_LIGHT if clock_mode == "pair-clock" else 0.0
        group = opalx.loc[(opalx["pair"] == pair) & (opalx["species"] == species)].sort_values("time_s")
        current_elapsed_s = 0.0
        for _, row in group.iterrows():
            target_elapsed_s = float(row["time_s"])
            while current_elapsed_s < target_elapsed_s - 1.0e-21:
                dt_s = min(DT_MAX_S, target_elapsed_s - current_elapsed_s)
                field_time_s = initial_clock_s + current_elapsed_s
                if current_elapsed_s >= SOURCE_RETIRE_TIME_S:
                    position = position + track12.velocity_from_momentum(momentum) * dt_s
                else:
                    position, momentum, _e, _b = track12.advance_with_fields(
                        position,
                        momentum,
                        charge_units,
                        field_time_s,
                        dt_s,
                        track12.anisotropic_total_lab_fields,
                        sources,
                    )
                current_elapsed_s += dt_s
            rows.append(
                {
                    "model": f"manufactured {source_selection} {clock_mode} sign={source_charge_sign:g}",
                    "pair": pair,
                    "kind": int(first["kind"]),
                    "species": species,
                    "global_step": int(row["global_step"]),
                    "time_s": target_elapsed_s,
                    "time_ps": 1.0e12 * target_elapsed_s,
                    "x_m": position[0],
                    "y_m": position[1],
                    "s_m": position[2],
                    "x_um": 1.0e6 * position[0],
                    "s_mm": 1.0e3 * position[2],
                }
            )
    return pd.DataFrame(rows)


def trajectory_metrics(opalx: pd.DataFrame, manufactured: pd.DataFrame) -> dict[str, float]:
    merged = opalx.merge(
        manufactured,
        on=["pair", "species", "global_step"],
        suffixes=("_opalx", "_manufactured"),
        validate="one_to_one",
    )
    dx_um = merged["x_um_manufactured"] - merged["x_um_opalx"]
    ds_um = 1.0e3 * (merged["s_mm_manufactured"] - merged["s_mm_opalx"])
    return {
        "rmse_x_um": float(np.sqrt(np.mean(dx_um * dx_um))),
        "mae_x_um": float(np.mean(np.abs(dx_um))),
        "max_abs_x_um": float(np.max(np.abs(dx_um))),
        "rmse_s_um": float(np.sqrt(np.mean(ds_um * ds_um))),
    }


def plot_tracks(
    opalx: pd.DataFrame,
    manufactured: pd.DataFrame,
    output: Path,
    title: str,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
) -> None:
    configure_matplotlib()
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "savefig.dpi": 220,
        }
    )
    fig, ax = plt.subplots(figsize=(9.0, 5.8), constrained_layout=True)
    for frame, alpha, width, label_prefix in ((opalx, 0.85, 1.2, "OPALX"), (manufactured, 0.75, 0.9, "manufactured")):
        for (species, pair), group in frame.groupby(["species", "pair"], sort=True):
            color = ELECTRON_COLOR if species == "electron" else POSITRON_COLOR
            group = group.sort_values("s_mm")
            ax.plot(
                group["s_mm"],
                group["x_um"],
                color=color,
                ls=LINE_STYLES.get(int(pair), "-") if label_prefix == "OPALX" else ":",
                lw=width,
                alpha=alpha,
            )
    ax.axhline(0.0, color="black", lw=0.8, ls=(0, (2, 6)))
    ax.set_xlabel("S - IP [mm]")
    ax.set_ylabel(r"x [$\mu$m]")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.minorticks_on()
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    handles = [
        plt.Line2D([0], [0], color="0.2", lw=1.3, ls="-", label="OPALX"),
        plt.Line2D([0], [0], color="0.2", lw=1.3, ls=":", label="manufactured"),
        plt.Line2D([0], [0], color=ELECTRON_COLOR, lw=1.8, label="e-"),
        plt.Line2D([0], [0], color=POSITRON_COLOR, lw=1.8, label="e+"),
    ]
    ax.legend(handles=handles, loc="best")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx-prefix", type=Path, default=SCRIPT_DIR / "opalx-initial-condition")
    parser.add_argument("--cain-csv", type=Path, default=DEFAULT_CAIN_CSV)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--output-stem", default="track-20000")
    parser.add_argument("--title", default="20,000-step witness tracks: OPALX vs manufactured model")
    parser.add_argument("--zoom-title", default="Witness tracks through source-retirement time: OPALX vs manufactured")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    track12 = load_track12_module()
    reference = pd.read_csv(args.cain_csv)
    opalx = load_opalx_tracks(args.opalx_prefix)
    opalx.to_csv(args.data_dir / f"{args.output_stem}-opalx-witness-trajectory.csv", index=False)

    candidates = []
    for source_charge_sign in (-1.0, +1.0):
        manufactured = run_manufactured_on_opalx_grid(
            opalx,
            reference,
            track12,
            "copropagating",
            source_charge_sign,
            "global-clock",
        )
        metrics = trajectory_metrics(opalx, manufactured)
        candidates.append(
            {
                "source_selection": "copropagating",
                "source_charge_sign": source_charge_sign,
                "clock_mode": "global-clock",
                **metrics,
                "frame": manufactured,
            }
        )
    best = min(candidates, key=lambda row: row["rmse_x_um"])
    best_frame = best.pop("frame")
    best_frame.to_csv(args.data_dir / f"{args.output_stem}-manufactured-on-opalx-grid.csv", index=False)
    summary = pd.DataFrame([{k: v for k, v in row.items() if k != "frame"} for row in candidates])
    interaction_metrics = trajectory_metrics(
        opalx.loc[opalx["time_ps"] <= 10.0],
        best_frame.loc[best_frame["time_ps"] <= 10.0],
    )
    for key, value in interaction_metrics.items():
        summary.loc[summary["rmse_x_um"].idxmin(), f"interaction_{key}"] = value
    summary = summary.sort_values(["rmse_x_um", "mae_x_um"]).reset_index(drop=True)
    summary.to_csv(args.data_dir / f"{args.output_stem}-comparison-summary.csv", index=False)
    plot_tracks(
        opalx,
        best_frame,
        args.figure_dir / f"{args.output_stem}-opalx-vs-manufactured.png",
        args.title,
    )
    plot_tracks(
        opalx.loc[opalx["time_ps"] <= 10.0],
        best_frame.loc[best_frame["time_ps"] <= 10.0],
        args.figure_dir / f"{args.output_stem}-opalx-vs-manufactured-interaction-zoom.png",
        args.zoom_title,
        xlim=(-1.0, 4.5),
        ylim=(-500.0, 500.0),
    )
    print(summary.to_string(index=False))
    print(f"Wrote {args.figure_dir / f'{args.output_stem}-opalx-vs-manufactured.png'}")
    print(f"Wrote {args.figure_dir / f'{args.output_stem}-opalx-vs-manufactured-interaction-zoom.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
