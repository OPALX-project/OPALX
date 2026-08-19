#!/usr/bin/env python3
"""Aggregate completed OPALX BeamBeam manufactured-field convergence cases."""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
COMPARATOR_SCRIPT = SCRIPT_DIR / "compare_opalx_fields.py"
CASE_PATTERN = re.compile(
    r"N(?P<n_k>\d+)k_M(?P<nxy>\d+)_S(?P<seed>\d+)(?:_MPI(?P<ranks>\d+))?"
)


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
    path: Path,
    *,
    log_x: bool,
    title: str,
) -> None:
    if frame.empty:
        return
    ordered = frame.sort_values(x)
    fig, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    axis.plot(ordered[x], 100.0 * ordered["relative_l2_E"], "o-", label=r"$E$")
    axis.plot(ordered[x], 100.0 * ordered["relative_l2_B"], "s-", label=r"$B$")
    if log_x:
        axis.set_xscale("log", base=2)
    else:
        axis.set_xticks(ordered[x].to_list())
    axis.set_xlabel(xlabel)
    axis.set_ylabel(r"relative $L_2$ error (\%)")
    axis.set_title(title)
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    comparator = load_script("beambeam_field_comparator", COMPARATOR_SCRIPT)
    field_model = comparator.load_field_model()
    parameters = field_model.load_parameters(field_model.DEFAULT_CONFIG)
    track12 = field_model.load_track12_module()

    rows: list[dict] = []
    for parent in sorted(args.case_root.glob("N*")):
        match = CASE_PATTERN.fullmatch(parent.name)
        case_dir = parent / "3sigma"
        if match is None or not (case_dir / "completed").exists():
            continue
        with warnings.catch_warnings():
            # A failed/empty field case deliberately produces NaN direction metrics;
            # retain it in invalid_cases.csv without noisy nanmedian warnings.
            warnings.filterwarnings("ignore", message=".*empty slice.*", category=RuntimeWarning)
            _, summary = comparator.analyze_case(
                case_dir,
                field_model,
                track12,
                parameters,
                step_index=1,
                source_model="physical-and-copied",
            )
        row = {
            "label": parent.name,
            "primary_macroparticles": 1000 * int(match.group("n_k")),
            "mesh_nxy": int(match.group("nxy")),
            "mesh_nz": 2 * int(match.group("nxy")),
            "seed": int(match.group("seed")),
            "mpi_ranks": int(match.group("ranks") or 1),
            "gpu_compute_wall_s": read_compute_time(case_dir),
        }
        row.update(summary)
        rows.append(row)

    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError(f"no completed convergence cases found below {args.case_root}")
    frame = frame.sort_values(["primary_macroparticles", "mesh_nxy", "seed"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "a100_convergence_summary.csv", index=False)

    valid = frame[frame["uncovered_probe_fraction"] == 0.0].copy()
    invalid = frame[frame["uncovered_probe_fraction"] != 0.0].copy()
    invalid.to_csv(args.output_dir / "invalid_cases.csv", index=False)

    particle_scan = valid[
        (valid["mesh_nxy"] == args.fixed_mesh_nxy)
        & (valid["seed"] == args.seed)
        & (valid["mpi_ranks"] == 1)
    ]
    mesh_scan = valid[
        (valid["primary_macroparticles"] == args.fixed_particles)
        & (valid["seed"] == args.seed)
        & (valid["mpi_ranks"] == 1)
    ]
    seed_scan = valid[
        (valid["primary_macroparticles"] == args.fixed_particles)
        & (valid["mesh_nxy"] == args.fixed_mesh_nxy)
        & (valid["mpi_ranks"] == 1)
    ]
    target_particles = args.occupancy_base_particles * (
        valid["mesh_nxy"] / args.occupancy_base_nxy
    ) ** 3
    occupancy_scan = valid[
        (valid["seed"] == args.seed)
        & (valid["mpi_ranks"] == 1)
        & (valid["primary_macroparticles"] == target_particles.astype(int))
    ]

    particle_scan.to_csv(args.output_dir / "particle_scan.csv", index=False)
    mesh_scan.to_csv(args.output_dir / "fixed_particle_mesh_scan.csv", index=False)
    occupancy_scan.to_csv(args.output_dir / "constant_occupancy_mesh_scan.csv", index=False)
    seed_scan.to_csv(args.output_dir / "seed_scan.csv", index=False)

    save_error_plot(
        particle_scan,
        "primary_macroparticles",
        "primary macroparticles",
        args.output_dir / "particle_convergence.png",
        log_x=True,
        title=rf"fixed ${args.fixed_mesh_nxy}\times{args.fixed_mesh_nxy}"
        rf"\times{2 * args.fixed_mesh_nxy}$ mesh",
    )
    save_error_plot(
        mesh_scan,
        "mesh_nxy",
        r"transverse mesh points $N_x=N_y$ ($N_z=2N_x$)",
        args.output_dir / "fixed_particle_mesh_convergence.png",
        log_x=False,
        title=f"fixed {args.fixed_particles:,} source macroparticles",
    )
    save_error_plot(
        occupancy_scan,
        "mesh_nxy",
        r"transverse mesh points $N_x=N_y$ ($N_z=2N_x$)",
        args.output_dir / "constant_occupancy_mesh_convergence.png",
        log_x=False,
        title="constant source macroparticles per mesh cell",
    )

    columns = [
        "label",
        "primary_macroparticles",
        "mesh_nxy",
        "seed",
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
