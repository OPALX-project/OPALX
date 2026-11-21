#!/usr/bin/env python3
"""Aggregate completed OPALX BeamBeam manufactured-field convergence cases."""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import sys
import warnings
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
PLOT_CACHE = MODEL_DIR / ".plot-cache"
PLOT_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(PLOT_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COMPARATOR_SCRIPT = SCRIPT_DIR / "compare_opalx_fields.py"
CASE_PATTERN = re.compile(
    r"N(?P<n_k>\d+)k_M(?P<nxy>\d+)_S(?P<seed>\d+)"
    r"(?P<fixed>_FIXED)?(?:_MPI(?P<ranks>\d+))?"
)
COLORS = {"E": "#0072B2", "B": "#D55E00", "rank": "#009E73"}


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
            "lines.linewidth": 1.8,
            "lines.markersize": 5.5,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.6,
            "savefig.bbox": "tight",
        }
    )


def save_figure(fig, path_without_suffix: Path) -> None:
    fig.savefig(path_without_suffix.with_suffix(".png"), dpi=300)
    fig.savefig(path_without_suffix.with_suffix(".pdf"))
    plt.close(fig)


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-root",
        type=Path,
        required=True,
        help="Directory containing N<k>k_M<nxy>_S<seed>/3sigma case directories.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260629)
    parser.add_argument("--fixed-mesh-nxy", type=int, default=64)
    parser.add_argument("--fixed-particles", type=int, default=400000)
    parser.add_argument("--occupancy-base-nxy", type=int, default=32)
    parser.add_argument("--occupancy-base-particles", type=int, default=400000)
    parser.add_argument("--probe-nx", type=int, default=51)
    parser.add_argument("--probe-nz", type=int, default=101)
    return parser.parse_args()


def read_compute_time(case_dir: Path) -> float | None:
    path = case_dir / "runtime_compute.txt"
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("real "):
            return float(line.split()[1])
    return None


def save_error_plot(
    frame: pd.DataFrame,
    x: str,
    xlabel: str,
    path_without_suffix: Path,
    *,
    log_x: bool,
    title: str,
) -> None:
    if frame.empty:
        return
    ordered = frame.sort_values(x)
    fig, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    axis.plot(
        ordered[x], 100.0 * ordered["relative_l2_E"], "o-",
        color=COLORS["E"], label=r"electric field $E$",
    )
    axis.plot(
        ordered[x], 100.0 * ordered["relative_l2_B"], "s-",
        color=COLORS["B"], label=r"magnetic field $B$",
    )
    if log_x:
        axis.set_xscale("log", base=2)
        tick_values = ordered[x].to_list()
        tick_labels = [
            f"{value / 1.0e6:g}M" if value >= 1.0e6 else f"{value / 1.0e3:g}k"
            for value in tick_values
        ]
        axis.set_xticks(tick_values, labels=tick_labels)
    else:
        axis.set_xticks(ordered[x].to_list())
    axis.set_xlabel(xlabel)
    axis.set_ylabel(r"relative $L_2$ error (\%)")
    axis.set_title(title)
    axis.legend(frameon=False)
    save_figure(fig, path_without_suffix)


def relative_field_difference(
    candidate: pd.DataFrame,
    baseline: pd.DataFrame,
    symbol: str,
    probe_nx: int,
    probe_nz: int,
) -> float:
    components = [f"{component}_{'opalx_V_per_m' if symbol == 'E' else 'opalx_T'}"
                  for component in (f"{symbol}x", f"{symbol}y", f"{symbol}z")]
    expected_count = probe_nx * probe_nz
    if len(candidate) != expected_count or len(baseline) != expected_count:
        raise RuntimeError(
            f"rank comparison expects a {probe_nx}x{probe_nz} probe grid "
            f"({expected_count} points), found {len(candidate)} and {len(baseline)}"
        )

    def add_grid_key(frame: pd.DataFrame) -> pd.DataFrame:
        keyed = frame.copy()
        x = keyed["x_m"].to_numpy()
        z = keyed["z_minus_ip_m"].to_numpy()
        ix = np.rint((x - np.min(x)) / np.ptp(x) * (probe_nx - 1)).astype(int)
        iz = np.rint((z - np.min(z)) / np.ptp(z) * (probe_nz - 1)).astype(int)
        keyed["probe_grid_key"] = iz * probe_nx + ix
        if keyed["probe_grid_key"].nunique() != expected_count:
            raise RuntimeError("could not recover unique logical probe-grid coordinates")
        return keyed

    candidate_keyed = add_grid_key(candidate)
    baseline_keyed = add_grid_key(baseline)
    position_columns = ["x_m", "y_m", "z_minus_ip_m"]
    merged = candidate_keyed[["probe_grid_key", *position_columns, *components]].merge(
        baseline_keyed[["probe_grid_key", *position_columns, *components]],
        on="probe_grid_key",
        suffixes=("_candidate", "_baseline"),
        validate="one_to_one",
    )
    candidate_values = merged[[f"{name}_candidate" for name in components]].to_numpy()
    baseline_values = merged[[f"{name}_baseline" for name in components]].to_numpy()
    denominator = np.linalg.norm(baseline_values)
    if denominator == 0.0:
        return float("nan")
    return float(np.linalg.norm(candidate_values - baseline_values) / denominator)


def save_rank_plots(rank_scan: pd.DataFrame, output_dir: Path) -> None:
    if rank_scan.empty:
        return
    ordered = rank_scan.sort_values("mpi_ranks")
    ranks = ordered["mpi_ranks"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), constrained_layout=True)
    axes[0].plot(
        ranks, 100.0 * ordered["relative_l2_E"], "o-",
        color=COLORS["E"], label=r"$E$ versus manufactured",
    )
    axes[0].plot(
        ranks, 100.0 * ordered["relative_l2_B"], "s-",
        color=COLORS["B"], label=r"$B$ versus manufactured",
    )
    axes[0].set_ylabel(r"relative $L_2$ error (\%)")
    axes[0].set_title("manufactured error versus MPI rank count")
    axes[0].legend(frameon=False)

    axes[1].plot(
        ranks, 100.0 * ordered["relative_l2_E_vs_rank1"], "o-",
        color=COLORS["E"], label=r"$E$ versus one rank",
    )
    axes[1].plot(
        ranks, 100.0 * ordered["relative_l2_B_vs_rank1"], "s-",
        color=COLORS["B"], label=r"$B$ versus one rank",
    )
    axes[1].set_ylabel(r"direct field difference (\%)")
    axes[1].set_title("fixed-source decomposition sensitivity")
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.set_xlabel("MPI ranks / A100 GPUs")
        axis.set_xticks(ranks)
    save_figure(fig, output_dir / "mpi_rank_sensitivity")

    runtime = ordered.dropna(subset=["gpu_compute_wall_s"]).copy()
    if runtime.empty or 1 not in set(runtime["mpi_ranks"]):
        return
    base_time = float(runtime.loc[runtime["mpi_ranks"] == 1, "gpu_compute_wall_s"].iloc[0])
    runtime["speedup"] = base_time / runtime["gpu_compute_wall_s"]
    runtime["parallel_efficiency"] = runtime["speedup"] / runtime["mpi_ranks"]
    runtime.to_csv(output_dir / "mpi_strong_scaling.csv", index=False)
    fig, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    axis.plot(
        runtime["mpi_ranks"], runtime["speedup"], "o-",
        color=COLORS["rank"], label="measured",
    )
    axis.plot(ranks, ranks, "--", color="0.45", label="ideal")
    axis.set_xlabel("MPI ranks / A100 GPUs")
    axis.set_ylabel("speedup of the two-step deck")
    axis.set_xticks(ranks)
    axis.set_title("strong scaling at fixed particles and mesh")
    axis.legend(frameon=False)
    save_figure(fig, output_dir / "mpi_strong_scaling")


def main() -> int:
    args = parse_args()
    configure_plots()
    comparator = load_script("beambeam_field_comparator", COMPARATOR_SCRIPT)
    field_model = comparator.load_field_model()
    parameters = field_model.load_parameters(field_model.DEFAULT_CONFIG)
    track12 = field_model.load_track12_module()

    rows: list[dict] = []
    comparisons: dict[str, pd.DataFrame] = {}
    for parent in sorted(args.case_root.glob("N*")):
        match = CASE_PATTERN.fullmatch(parent.name)
        case_dir = parent / "3sigma"
        if match is None or not (case_dir / "completed").exists():
            continue
        with warnings.catch_warnings():
            # A failed/empty field case deliberately produces NaN direction metrics;
            # retain it in invalid_cases.csv without noisy nanmedian warnings.
            warnings.filterwarnings("ignore", message=".*empty slice.*", category=RuntimeWarning)
            comparison, summary = comparator.analyze_case(
                case_dir,
                field_model,
                track12,
                parameters,
                step_index=1,
                source_model="physical-and-copied",
            )
        comparisons[parent.name] = comparison
        row = {
            "label": parent.name,
            "primary_macroparticles": 1000 * int(match.group("n_k")),
            "mesh_nxy": int(match.group("nxy")),
            "mesh_nz": 2 * int(match.group("nxy")),
            "seed": int(match.group("seed")),
            "mpi_ranks": int(match.group("ranks") or 1),
            "primary_sampling": (
                "fixed_fromfile" if match.group("fixed") else "gauss"
            ),
            "gpu_compute_wall_s": read_compute_time(case_dir),
        }
        row.update(summary)
        rows.append(row)

    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError(f"no completed convergence cases found below {args.case_root}")
    frame = frame.sort_values(
        ["primary_macroparticles", "mesh_nxy", "seed", "primary_sampling", "mpi_ranks"]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "a100_convergence_summary.csv", index=False)

    valid = frame[frame["uncovered_probe_fraction"] == 0.0].copy()
    invalid = frame[frame["uncovered_probe_fraction"] != 0.0].copy()
    invalid.to_csv(args.output_dir / "invalid_cases.csv", index=False)

    particle_scan = valid[
        (valid["mesh_nxy"] == args.fixed_mesh_nxy)
        & (valid["seed"] == args.seed)
        & (valid["mpi_ranks"] == 1)
        & (valid["primary_sampling"] == "gauss")
    ]
    mesh_scan = valid[
        (valid["primary_macroparticles"] == args.fixed_particles)
        & (valid["seed"] == args.seed)
        & (valid["mpi_ranks"] == 1)
        & (valid["primary_sampling"] == "gauss")
    ]
    seed_scan = valid[
        (valid["primary_macroparticles"] == args.fixed_particles)
        & (valid["mesh_nxy"] == args.fixed_mesh_nxy)
        & (valid["mpi_ranks"] == 1)
        & (valid["primary_sampling"] == "gauss")
    ]
    target_particles = args.occupancy_base_particles * (
        valid["mesh_nxy"] / args.occupancy_base_nxy
    ) ** 3
    occupancy_scan = valid[
        (valid["seed"] == args.seed)
        & (valid["mpi_ranks"] == 1)
        & (valid["primary_sampling"] == "gauss")
        & (valid["primary_macroparticles"] == target_particles.astype(int))
    ]
    rank_scan = valid[
        (valid["primary_macroparticles"] == args.fixed_particles)
        & (valid["mesh_nxy"] == args.fixed_mesh_nxy)
        & (valid["seed"] == args.seed)
        & (valid["primary_sampling"] == "fixed_fromfile")
    ].copy()
    rank1_rows = rank_scan[rank_scan["mpi_ranks"] == 1]
    if rank1_rows.empty:
        rank_scan["relative_l2_E_vs_rank1"] = np.nan
        rank_scan["relative_l2_B_vs_rank1"] = np.nan
    else:
        baseline_label = str(rank1_rows.iloc[0]["label"])
        baseline = comparisons[baseline_label]
        rank_scan["relative_l2_E_vs_rank1"] = [
            relative_field_difference(
                comparisons[str(label)], baseline, "E", args.probe_nx, args.probe_nz
            )
            for label in rank_scan["label"]
        ]
        rank_scan["relative_l2_B_vs_rank1"] = [
            relative_field_difference(
                comparisons[str(label)], baseline, "B", args.probe_nx, args.probe_nz
            )
            for label in rank_scan["label"]
        ]

    particle_scan.to_csv(args.output_dir / "particle_scan.csv", index=False)
    mesh_scan.to_csv(args.output_dir / "fixed_particle_mesh_scan.csv", index=False)
    occupancy_scan.to_csv(args.output_dir / "constant_occupancy_mesh_scan.csv", index=False)
    seed_scan.to_csv(args.output_dir / "seed_scan.csv", index=False)
    rank_scan.to_csv(args.output_dir / "mpi_rank_scan.csv", index=False)

    save_error_plot(
        particle_scan,
        "primary_macroparticles",
        "primary macroparticles",
        args.output_dir / "particle_convergence",
        log_x=True,
        title=rf"fixed ${args.fixed_mesh_nxy}\times{args.fixed_mesh_nxy}"
        rf"\times{2 * args.fixed_mesh_nxy}$ mesh",
    )
    save_error_plot(
        mesh_scan,
        "mesh_nxy",
        r"transverse mesh points $N_x=N_y$ ($N_z=2N_x$)",
        args.output_dir / "fixed_particle_mesh_convergence",
        log_x=False,
        title=f"fixed {args.fixed_particles:,} source macroparticles",
    )
    save_error_plot(
        occupancy_scan,
        "mesh_nxy",
        r"transverse mesh points $N_x=N_y$ ($N_z=2N_x$)",
        args.output_dir / "constant_occupancy_mesh_convergence",
        log_x=False,
        title="constant source macroparticles per mesh cell",
    )
    save_rank_plots(rank_scan, args.output_dir)

    columns = [
        "label",
        "primary_macroparticles",
        "mesh_nxy",
        "seed",
        "primary_sampling",
        "gpu_compute_wall_s",
        "relative_l2_E",
        "relative_l2_B",
        "median_E_magnitude_ratio",
        "median_B_magnitude_ratio",
        "median_E_direction_cosine",
        "median_B_direction_cosine",
        "uncovered_probe_fraction",
    ]
    print(frame[columns].to_string(index=False, float_format=lambda value: f"{value:.8e}"))
    if not invalid.empty:
        print("\nExcluded invalid copied-field cases:")
        print(invalid[["label", "uncovered_probe_fraction"]].to_string(index=False))
    print("\nBaseline seed spread:")
    print(
        seed_scan[["seed", "relative_l2_E", "relative_l2_B"]]
        .describe()
        .to_string(float_format=lambda value: f"{value:.8e}")
    )
    if not rank_scan.empty:
        print("\nFixed-primary MPI decomposition scan:")
        print(
            rank_scan[
                [
                    "mpi_ranks",
                    "gpu_compute_wall_s",
                    "relative_l2_E",
                    "relative_l2_B",
                    "relative_l2_E_vs_rank1",
                    "relative_l2_B_vs_rank1",
                ]
            ].to_string(index=False, float_format=lambda value: f"{value:.8e}")
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
