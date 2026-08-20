#!/usr/bin/env python3
"""Check timed BeamBeam witness fields for 1/2/4-rank invariance."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PREPARE_SCRIPT = HERE / "prepare_timed_track12.py"
TEMPLATE = HERE / "witness_gather_rank.in.template"
DEFAULT_OPALX = ROOT / "build_openmp" / "src" / "opalx"
DEFAULT_WORK_DIR = ROOT / "build_openmp" / "test" / "BeamBeam" / "witness-gather-mpi"
DEFAULT_RANKS = (1, 2, 4)
DEFAULT_RELATIVE_TOLERANCE = 5.0e-4
PRIMARY_MACROPARTICLES = 100_000
SIGMA_XY_M = 1.944325075701e-6
INITIAL_P = np.array([0.0, 0.0, np.sqrt(3.0)])
EXPECTED_XY_M = SIGMA_XY_M * np.array(
    [
        [+1.0, +0.25],
        [-1.0, -0.25],
        [+0.25, -1.0],
        [-0.25, +1.0],
        [+2.0, +0.50],
        [-2.0, -0.50],
        [+0.50, -2.0],
        [-0.50, +2.0],
    ],
    dtype=float,
)


def load_prepare_module():
    spec = importlib.util.spec_from_file_location("prepare_timed_track12", PREPARE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {PREPARE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_witness_file(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(EXPECTED_XY_M)}\n")
        stream.write("x y z px py pz birth_time\n")
        for x_m, y_m in EXPECTED_XY_M:
            stream.write(
                f"{x_m:.16e} {y_m:.16e} 0.0000000000000000e+00 "
                f"{INITIAL_P[0]:.16e} {INITIAL_P[1]:.16e} {INITIAL_P[2]:.16e} "
                "0.0000000000000000e+00\n"
            )


def prepare_inputs(work_dir: Path) -> dict[str, object]:
    input_dir = work_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    prepare = load_prepare_module()
    primary_path = input_dir / "primary_fixed.fromfile"
    metadata = prepare.write_fixed_primary(primary_path, PRIMARY_MACROPARTICLES)
    write_witness_file(input_dir / "witness_electrons.emittedfromfile")
    write_witness_file(input_dir / "witness_positrons.emittedfromfile")
    (input_dir / "primary_fixed_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata


def render_deck(run_dir: Path) -> Path:
    deck = TEMPLATE.read_text(encoding="utf-8")
    replacements = {
        "@PRIMARY_MACROPARTICLES@": str(PRIMARY_MACROPARTICLES),
        "@WITNESS_COUNT@": str(len(EXPECTED_XY_M)),
    }
    for token, value in replacements.items():
        deck = deck.replace(token, value)
    if "@" in deck:
        raise RuntimeError("unresolved token remains in witness-gather deck")
    path = run_dir / "witness_gather_rank.in"
    path.write_text(deck, encoding="utf-8")
    return path


def run_rank_case(
    opalx: Path,
    mpiexec: Path,
    rank_count: int,
    work_dir: Path,
    force: bool,
    openmpi_local: bool,
) -> tuple[Path, float]:
    run_dir = work_dir / f"rank{rank_count}"
    run_dir.mkdir(parents=True, exist_ok=True)
    input_path = render_deck(run_dir)
    electron_h5 = run_dir / "witness_gather_rank_c1.h5"
    positron_h5 = run_dir / "witness_gather_rank_c2.h5"
    if force:
        for generated in run_dir.glob("witness_gather_rank*"):
            if generated != input_path:
                generated.unlink()
    if electron_h5.exists() and positron_h5.exists():
        return run_dir, 0.0

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    command = [str(mpiexec)]
    if openmpi_local:
        command.extend(["--map-by", "slot:OVERSUBSCRIBE", "--bind-to", "none"])
    command.extend(["-n", str(rank_count), str(opalx), input_path.name])
    start = time.monotonic()
    completed = subprocess.run(
        command,
        cwd=run_dir,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    runtime_s = time.monotonic() - start
    (run_dir / "opalx.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"{rank_count}-rank OPALX failed with code {completed.returncode}; "
            f"see {run_dir / 'opalx.log'}"
        )
    if not electron_h5.exists() or not positron_h5.exists():
        raise RuntimeError(f"{rank_count}-rank OPALX did not produce both witness H5 files")
    return run_dir, runtime_s


def sorted_steps(h5_file: h5py.File) -> list[h5py.Group]:
    return sorted(h5_file.values(), key=lambda group: int(group.name.rsplit("#", 1)[1]))


def match_expected_xy(actual_xy: np.ndarray) -> np.ndarray:
    if len(actual_xy) != len(EXPECTED_XY_M):
        raise RuntimeError(
            f"expected {len(EXPECTED_XY_M)} witnesses, found {len(actual_xy)}"
        )
    remaining = set(range(len(actual_xy)))
    assignment: list[int] = []
    for expected in EXPECTED_XY_M:
        index = min(remaining, key=lambda candidate: np.linalg.norm(actual_xy[candidate] - expected))
        if np.linalg.norm(actual_xy[index] - expected) > 0.25 * SIGMA_XY_M:
            raise RuntimeError("could not match emitted witness to its input transverse position")
        assignment.append(index)
        remaining.remove(index)
    return np.asarray(assignment, dtype=int)


def read_first_complete_sample(path: Path, species: str, rank_count: int) -> pd.DataFrame:
    required = ("x", "y", "z", "px", "py", "pz", "Ex", "Ey", "Ez", "Bx", "By", "Bz")
    with h5py.File(path, "r") as h5_file:
        complete = [group for group in sorted_steps(h5_file) if len(group["x"]) == len(EXPECTED_XY_M)]
        if not complete:
            populations = [len(group["x"]) for group in sorted_steps(h5_file)]
            raise RuntimeError(f"{path}: no complete witness sample; populations={populations}")
        group = complete[0]
        missing = [name for name in required if name not in group]
        if missing:
            raise RuntimeError(f"{path}:{group.name} missing {missing}")
        ref_r = np.asarray(group.attrs["RefPartR"], dtype=float)
        actual_xy = np.column_stack((group["x"][:] + ref_r[0], group["y"][:] + ref_r[1]))
        order = match_expected_xy(actual_xy)
        frame = pd.DataFrame(
            {
                "witness": np.arange(len(order), dtype=int),
                "x_m": actual_xy[order, 0],
                "y_m": actual_xy[order, 1],
                "z_m": group["z"][:][order] + ref_r[2],
                "px": group["px"][:][order],
                "py": group["py"][:][order],
                "pz": group["pz"][:][order],
                "Ex": group["Ex"][:][order],
                "Ey": group["Ey"][:][order],
                "Ez": group["Ez"][:][order],
                "Bx": group["Bx"][:][order],
                "By": group["By"][:][order],
                "Bz": group["Bz"][:][order],
            }
        )
        frame.insert(0, "species", species)
        frame.insert(0, "rank_count", rank_count)
        frame["global_step"] = int(np.asarray(group.attrs["GlobalTrackStep"]).reshape(-1)[0])
        return frame


def relative_difference(candidate: np.ndarray, reference: np.ndarray) -> float:
    denominator = np.linalg.norm(reference)
    if denominator == 0.0:
        return 0.0 if np.linalg.norm(candidate) == 0.0 else float("inf")
    return float(np.linalg.norm(candidate - reference) / denominator)


def analyze(rank_frames: dict[int, pd.DataFrame], runtimes: dict[int, float]) -> pd.DataFrame:
    baseline = rank_frames[1].sort_values(["species", "witness"])
    rows: list[dict[str, float | int]] = []
    for rank_count, frame in sorted(rank_frames.items()):
        ordered = frame.sort_values(["species", "witness"])
        electron = ordered.loc[ordered["species"].eq("electron")]
        positron = ordered.loc[ordered["species"].eq("positron")]
        baseline_ordered = baseline
        momentum = ordered[["px", "py", "pz"]].to_numpy()
        baseline_momentum = baseline_ordered[["px", "py", "pz"]].to_numpy()
        momentum_kick = momentum - INITIAL_P
        baseline_momentum_kick = baseline_momentum - INITIAL_P
        row: dict[str, float | int] = {
            "rank_count": rank_count,
            "runtime_s": runtimes[rank_count],
            "witness_samples": len(ordered),
            "relative_E_vs_rank1": relative_difference(
                ordered[["Ex", "Ey", "Ez"]].to_numpy(),
                baseline_ordered[["Ex", "Ey", "Ez"]].to_numpy(),
            ),
            "relative_B_vs_rank1": relative_difference(
                ordered[["Bx", "By", "Bz"]].to_numpy(),
                baseline_ordered[["Bx", "By", "Bz"]].to_numpy(),
            ),
            "relative_momentum_vs_rank1": relative_difference(
                momentum,
                baseline_momentum,
            ),
            "relative_momentum_kick_vs_rank1": relative_difference(
                momentum_kick,
                baseline_momentum_kick,
            ),
            "relative_E_charge_symmetry": relative_difference(
                positron[["Ex", "Ey", "Ez"]].to_numpy(),
                electron[["Ex", "Ey", "Ez"]].to_numpy(),
            ),
            "relative_B_charge_symmetry": relative_difference(
                positron[["Bx", "By", "Bz"]].to_numpy(),
                electron[["Bx", "By", "Bz"]].to_numpy(),
            ),
        }
        electron_kick = electron["px"].to_numpy() - INITIAL_P[0]
        positron_kick = positron["px"].to_numpy() - INITIAL_P[0]
        kick_scale = max(np.linalg.norm(electron_kick), np.finfo(float).tiny)
        row["relative_px_kick_charge_symmetry"] = float(
            np.linalg.norm(electron_kick + positron_kick) / kick_scale
        )
        rows.append(row)
    return pd.DataFrame(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--mpiexec", type=Path, default=Path(shutil.which("mpiexec") or "mpiexec"))
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--ranks", type=int, nargs="+", default=list(DEFAULT_RANKS))
    parser.add_argument("--relative-tolerance", type=float, default=DEFAULT_RELATIVE_TOLERANCE)
    parser.add_argument(
        "--openmpi-local",
        action="store_true",
        help="Use unbound oversubscribed Open MPI mapping for a local workstation.",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-enforce", action="store_true")
    args = parser.parse_args()

    if 1 not in args.ranks:
        raise ValueError("--ranks must include the one-rank baseline")
    if any(rank <= 0 for rank in args.ranks):
        raise ValueError("all MPI rank counts must be positive")
    if args.relative_tolerance <= 0.0:
        raise ValueError("--relative-tolerance must be positive")
    if not args.opalx.is_file():
        raise FileNotFoundError(args.opalx)
    if not args.mpiexec.is_file():
        raise FileNotFoundError(args.mpiexec)
    args.opalx = args.opalx.resolve()
    args.mpiexec = args.mpiexec.resolve()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    primary_metadata = prepare_inputs(args.work_dir)
    rank_frames: dict[int, pd.DataFrame] = {}
    runtimes: dict[int, float] = {}
    for rank_count in args.ranks:
        run_dir, runtime_s = run_rank_case(
            args.opalx,
            args.mpiexec,
            rank_count,
            args.work_dir,
            args.force,
            args.openmpi_local,
        )
        runtimes[rank_count] = runtime_s
        rank_frames[rank_count] = pd.concat(
            [
                read_first_complete_sample(
                    run_dir / "witness_gather_rank_c1.h5", "electron", rank_count
                ),
                read_first_complete_sample(
                    run_dir / "witness_gather_rank_c2.h5", "positron", rank_count
                ),
            ],
            ignore_index=True,
        )

    samples = pd.concat(rank_frames.values(), ignore_index=True)
    summary = analyze(rank_frames, runtimes)
    samples.to_csv(args.work_dir / "witness_gather_samples.csv", index=False)
    summary.to_csv(args.work_dir / "witness_gather_rank_summary.csv", index=False)
    report = {
        "relative_tolerance": args.relative_tolerance,
        "ranks": args.ranks,
        "primary": primary_metadata,
        "opalx": {"path": str(args.opalx.resolve()), "sha256": sha256(args.opalx)},
        "summary": summary.to_dict(orient="records"),
    }
    (args.work_dir / "witness_gather_rank_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(f"results: {args.work_dir.resolve()}")

    if args.no_enforce:
        return 0

    baseline = rank_frames[1]
    baseline_norms = {
        "E": np.linalg.norm(baseline[["Ex", "Ey", "Ez"]].to_numpy()),
        "B": np.linalg.norm(baseline[["Bx", "By", "Bz"]].to_numpy()),
        "momentum kick": np.linalg.norm(
            baseline[["px", "py", "pz"]].to_numpy() - INITIAL_P
        ),
    }
    invalid_baselines = [
        f"{name}={value:.6g}"
        for name, value in baseline_norms.items()
        if not np.isfinite(value) or value <= 0.0
    ]
    if invalid_baselines:
        raise RuntimeError(
            "one-rank witness result does not contain finite nonzero fields and kick: "
            + "; ".join(invalid_baselines)
        )

    comparisons = summary.loc[summary["rank_count"].ne(1)]
    rank_comparison_columns = (
        "relative_E_vs_rank1",
        "relative_B_vs_rank1",
        "relative_momentum_vs_rank1",
        "relative_momentum_kick_vs_rank1",
    )
    symmetry_columns = (
        "relative_E_charge_symmetry",
        "relative_B_charge_symmetry",
        "relative_px_kick_charge_symmetry",
    )
    failures = [
        f"rank {int(row.rank_count)} {column}={getattr(row, column):.6g}"
        for row in comparisons.itertuples(index=False)
        for column in rank_comparison_columns
        if not np.isfinite(getattr(row, column))
        or getattr(row, column) > args.relative_tolerance
    ]
    failures.extend(
        f"rank {int(row.rank_count)} {column}={getattr(row, column):.6g}"
        for row in summary.itertuples(index=False)
        for column in symmetry_columns
        if not np.isfinite(getattr(row, column))
        or getattr(row, column) > args.relative_tolerance
    )
    if failures:
        raise RuntimeError(
            f"witness gather exceeds relative tolerance {args.relative_tolerance:.6g}: "
            + "; ".join(failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
