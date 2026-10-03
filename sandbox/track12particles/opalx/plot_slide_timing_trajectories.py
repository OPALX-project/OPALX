#!/usr/bin/env python3
"""Compare slide-timed OPALX witness trajectories with the manufactured model."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TRACK12_MODULE = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_INITIAL = ROOT / "sandbox" / "track12particles" / "outputs" / "track12_initial_conditions.csv"
DEFAULT_CAIN = ROOT / "sandbox" / "track12particles" / "outputs" / "track12_trajectory.csv"
DEFAULT_PAIR_DIR = ROOT / "sandbox" / "track12particles" / "opalx" / "timing_pairs_100"
C_LIGHT = 299_792_458.0


def configure_matplotlib(output_dir: Path) -> None:
    cache_root = output_dir / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TRACK12_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_sources(module, row: pd.Series, source_selection: str):
    sources = module.make_anisotropic_sources(
        module.DEFAULT_SOURCE_KINETIC_MEV,
        module.DEFAULT_SOURCE_CHARGE_C,
        module.SIGMA_XY_M,
        module.SIGMA_XY_M,
        module.SIGMA_Z_M,
        source_selection,
        0.0,
        module.DEFAULT_MC_SOURCE_PARTICLES,
        module.DEFAULT_MC_SEED,
    )
    shifted = []
    for source in sources:
        center = source.center_t0_m.copy()
        center[2] -= float(row["t"])
        shifted.append(
            module.RigidAnisotropicGaussianSource(
                source.name,
                charge_C=source.charge_C,
                sigma_lab_m=(module.SIGMA_XY_M, module.SIGMA_XY_M, module.SIGMA_Z_M),
                beta_z=source.beta_z,
                center_t0_m=center,
                t0_s=0.0,
            )
        )
    return tuple(shifted)


def manufactured_trajectory(
    initial: pd.DataFrame,
    module,
    source_selection: str,
    steps: int,
    dt_s: float,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for _, first in initial.sort_values(["pair", "kind"]).iterrows():
        position = first[["x", "y", "s"]].to_numpy(dtype=float)
        momentum = first[["Px_beta_gamma", "Py_beta_gamma", "Ps_beta_gamma"]].to_numpy(dtype=float)
        charge_units = float(first["charge_units"])
        sources = make_sources(module, first, source_selection)
        for step in range(steps + 1):
            rows.append(
                {
                    "model": "manufactured",
                    "source_selection": source_selection,
                    "species": first["species"],
                    "pair": int(first["pair"]),
                    "step": step,
                    "time_s": step * dt_s,
                    "x_m": position[0],
                    "y_m": position[1],
                    "s_m": position[2],
                    "px": momentum[0],
                    "py": momentum[1],
                    "pz": momentum[2],
                }
            )
            if step == steps:
                break
            position, momentum, _e_field, _b_field = module.advance_with_fields(
                position,
                momentum,
                charge_units,
                step * dt_s,
                dt_s,
                module.anisotropic_total_lab_fields,
                sources,
            )
    return pd.DataFrame(rows)


def opalx_trajectory(pair_dir: Path, steps: int, species_assignment: str) -> pd.DataFrame:
    if species_assignment == "h5-charge":
        species_by_container = {1: "electron", 2: "positron"}
    elif species_assignment == "cain-branch":
        species_by_container = {1: "positron", 2: "electron"}
    else:
        raise ValueError(f"unsupported species assignment {species_assignment!r}")

    rows: list[dict[str, float | int | str]] = []
    for pair in range(1, 7):
        stem = f"track12_pair{pair}_slide_timed_one_source_1fs_400k_{steps}steps"
        path = pair_dir / f"pair{pair}" / f"{stem}_kicks.csv"
        kicks = pd.read_csv(path)
        # The diagnostic can include step == steps. Keep one pre-kick trajectory
        # sample per step, including step 0 through the requested final step.
        kicks = kicks.loc[kicks["step"].between(0, steps)].copy()
        for _, row in kicks.iterrows():
            container = int(row["container"])
            rows.append(
                {
                    "model": "opalx",
                    "species_assignment": species_assignment,
                    "container": container,
                    "species": species_by_container[container],
                    "pair": pair,
                    "step": int(row["step"]),
                    "time_s": float(row["time_s"]),
                    "x_m": float(row["x_m"]),
                    "y_m": float(row["y_m"]),
                    "s_m": float(row["s_minus_ip_m"]),
                    "px": float(row["Px_before"]),
                    "py": float(row["Py_before"]),
                    "pz": float(row["Pz_before"]),
                }
            )
    return pd.DataFrame(rows)


def compare(opalx: pd.DataFrame, manufactured: pd.DataFrame) -> pd.DataFrame:
    merged = opalx.merge(
        manufactured,
        on=["species", "pair", "step"],
        suffixes=("_opalx", "_manufactured"),
        validate="one_to_one",
    )
    for column in ("x_m", "y_m", "s_m", "px", "py", "pz"):
        merged[f"d{column}"] = merged[f"{column}_opalx"] - merged[f"{column}_manufactured"]
    return merged


def cain_trajectory(path: Path, steps: int, dt_s: float) -> pd.DataFrame:
    data = pd.read_csv(path)
    rows: list[pd.DataFrame] = []
    max_elapsed_s = steps * dt_s
    for (_species, _pair), group in data.groupby(["species", "pair"], sort=False):
        group = group.sort_values("step").copy()
        t0_m = float(group["t"].iloc[0])
        group["elapsed_time_s"] = (group["t"] - t0_m) / C_LIGHT
        group = group.loc[group["elapsed_time_s"] <= max_elapsed_s + 1.0e-18].copy()
        group["model"] = "cain"
        group["time_s"] = group["elapsed_time_s"]
        group["x_m"] = group["x"]
        group["y_m"] = group["y"]
        group["s_m"] = group["s"]
        group["px"] = group["Px_beta_gamma"]
        group["py"] = group["Py_beta_gamma"]
        group["pz"] = group["Ps_beta_gamma"]
        rows.append(group[["model", "species", "pair", "step", "time_s", "x_m", "y_m", "s_m", "px", "py", "pz"]])
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def plot_trajectories(
    comparison: pd.DataFrame,
    cain: pd.DataFrame,
    output: Path,
    title: str,
) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {"electron": "#b51f2e", "positron": "#1f5fb5"}
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for (species, pair), group in comparison.groupby(["species", "pair"], sort=True):
        group = group.sort_values("step")
        color = colors[species]
        alpha = 0.4 + 0.08 * pair
        ax.plot(
            group["s_m_opalx"] * 1.0e3,
            group["x_m_opalx"] * 1.0e6,
            color=color,
            lw=1.4,
            alpha=alpha,
        )
        ax.plot(
            group["s_m_manufactured"] * 1.0e3,
            group["x_m_manufactured"] * 1.0e6,
            color=color,
            lw=1.1,
            ls="--",
            alpha=alpha,
        )
    if not cain.empty:
        for (species, pair), group in cain.groupby(["species", "pair"], sort=True):
            group = group.sort_values("time_s")
            color = colors[species]
            alpha = 0.4 + 0.08 * pair
            ax.plot(
                group["s_m"] * 1.0e3,
                group["x_m"] * 1.0e6,
                color=color,
                lw=1.0,
                ls=":",
                alpha=alpha,
            )
    ax.set_xlabel(r"$S-S_{\rm IP}$ [mm]")
    ax.set_ylabel(r"$x$ [$\mu$m]")
    ax.set_title(title)
    ax.grid(True, color="0.88", lw=0.6)
    handles = [
        Line2D([0], [0], color="0.2", lw=1.4, label="OPALX"),
        Line2D([0], [0], color="0.2", lw=1.1, ls="--", label="manufactured"),
        Line2D([0], [0], color="0.2", lw=1.0, ls=":", label="CAIN"),
        Line2D([0], [0], color=colors["electron"], lw=1.4, label="electron"),
        Line2D([0], [0], color=colors["positron"], lw=1.4, label="positron"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def print_summary(comparison: pd.DataFrame) -> None:
    print("Trajectory OPALX minus manufactured:")
    for column in ("dx_m", "dy_m", "ds_m", "dpx", "dpy", "dpz"):
        values = comparison[column]
        scale = 1.0e6 if column.endswith("_m") else 1.0
        unit = " um" if column.endswith("_m") else ""
        print(f"  max |{column}| = {values.abs().max() * scale:.9g}{unit}")
    final = comparison.sort_values("step").groupby(["species", "pair"], as_index=False).tail(1)
    print("\nFinal-step OPALX minus manufactured:")
    for column in ("dx_m", "dpx", "ds_m"):
        values = final[column]
        scale = 1.0e6 if column.endswith("_m") else 1.0
        unit = " um" if column.endswith("_m") else ""
        print(f"  max |{column}| = {values.abs().max() * scale:.9g}{unit}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-csv", type=Path, default=DEFAULT_INITIAL)
    parser.add_argument("--cain-csv", type=Path, default=DEFAULT_CAIN)
    parser.add_argument("--no-cain", action="store_true")
    parser.add_argument("--pair-dir", type=Path, default=DEFAULT_PAIR_DIR)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--dt-s", type=float, default=1.0e-15)
    parser.add_argument("--source-selection", choices=("copropagating", "oncoming"), default="copropagating")
    parser.add_argument("--species-assignment", choices=("h5-charge", "cain-branch"), default="h5-charge")
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--output-plot", type=Path)
    args = parser.parse_args()

    if args.output_csv is None:
        args.output_csv = args.pair_dir / f"slide_timing_{args.steps}steps_{args.source_selection}_trajectory_comparison.csv"
    if args.output_plot is None:
        args.output_plot = args.pair_dir / f"slide_timing_{args.steps}steps_{args.source_selection}_x_vs_s.png"

    module = load_track12_module()
    initial = pd.read_csv(args.initial_csv)
    opalx = opalx_trajectory(args.pair_dir, args.steps, args.species_assignment)
    manufactured = manufactured_trajectory(initial, module, args.source_selection, args.steps, args.dt_s)
    cain = pd.DataFrame() if args.no_cain else cain_trajectory(args.cain_csv, args.steps, args.dt_s)
    comparison = compare(opalx, manufactured)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.output_csv, index=False)
    plot_trajectories(
        comparison,
        cain,
        args.output_plot,
        f"slide-timed 100-step OPALX, manufactured, CAIN ({args.source_selection})",
    )
    print_summary(comparison)
    print(f"\nWrote {args.output_csv}")
    print(f"Wrote {args.output_plot}")


if __name__ == "__main__":
    main()
