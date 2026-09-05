#!/usr/bin/env python3
"""Reproducible joint DBA DT/amplitude/Richardson sweep; never edits production code.

Run from the repository root with ~/.venv-h6/bin/python. Independent runs use
separate temporary directories and one OpenMP thread each. Existing validated
results are reused only with matching executable, source, input and settings.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import time

from convergence_dt import DBA, MAP2, dba, np, pd, replace_dt, run_case
from convergence_perturbations import configure_map

REPOSITORY = MAP2.parents[1]
DEFAULT_DT = [1e-10, 3e-11, 1e-11, 3e-12, 1e-12, 3e-13, 1e-13]
DEFAULT_EPSILON = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6]
# S = diag(1 m, 1, 1 m, 1, 1 m, 1): A = S^-1 M S is dimensionless.
# These numerical values equal the SI entries, but the scale choice is explicit.
COORDINATE_SCALES = np.ones(6)
SOURCE_PATHS = [
    *[f"src/Algorithms/{name}.{suffix}" for name in (
        "ExternalFieldRayTracker", "LinearTransferMapBuilder", "OrbitThreader"
    ) for suffix in ("h", "cpp")],
    "src/Structure/LinearTransferMap.h", "src/BasicActions/Option.cpp",
    "src/Algorithms/CompensatedSum.h",
    "src/Utilities/Options.h", "src/Utilities/Options.cpp",
    "sandbox/map-2/check_maps.py", "sandbox/map-2/dba/convergence_dt.py",
    "sandbox/map-2/dba/convergence_perturbations.py",
]


def digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def label(value: float) -> str:
    return f"{value:.8g}"


def matrix_path(directory: Path, dt: float) -> Path:
    return directory / f"map-2-dba-dt-{dt:.0e}-matrix.txt"


def identity(executable: Path, source: str) -> dict:
    return {
        "executable_sha256": digest(executable.read_bytes()),
        "input_template_sha256": digest(source.encode()),
        "distribution_sha256": digest((MAP2 / "reference-particle.txt").read_bytes()),
        "sources_sha256": {path: digest((REPOSITORY / path).read_bytes())
                           for path in SOURCE_PATHS},
    }


def run_one(executable: Path, source: str, output: Path, provenance: dict,
            dt: float, epsilon: float, levels: int, ranks: int, mpi_args: tuple[str, ...]) -> dict:
    directory = output / f"L{levels}" / f"dt-{label(dt)}" / f"eps-{label(epsilon)}" / f"ranks-{ranks}"
    directory.mkdir(parents=True, exist_ok=True)
    configured = replace_dt(configure_map(source, epsilon, levels), dt)
    settings = dict(dt_s=dt, epsilon=epsilon, richardson_levels=levels, ranks=ranks,
                    input_sha256=digest(configured.encode()), identity=provenance)
    signature = digest(json.dumps(settings, sort_keys=True).encode())
    result_file = directory / "result.json"
    if result_file.exists():
        previous = json.loads(result_file.read_text())
        if previous["signature"] != signature:
            raise RuntimeError(f"provenance mismatch: use a new output directory: {directory}")
        for name, sha in previous["artifacts_sha256"].items():
            if digest((directory / name).read_bytes()) != sha:
                raise RuntimeError(f"cached artifact changed: {directory / name}")
        return previous | {"cached": True}

    (directory / "map-2-dba.in").write_text(configured)
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="opalx-map-grid-") as temporary:
        result = run_case(executable, configured, dt, Path(temporary), directory, ranks, mpi_args)
    result.update(epsilon=epsilon, richardson_levels=levels, ranks=ranks,
                  finest_epsilon=epsilon / 2**levels, rays_per_segment=12 * (levels + 1),
                  elapsed_s=time.perf_counter() - started, integration_method="BORIS",
                  signature=signature, path=str(directory.relative_to(output)))
    stdout_file = directory / f"map-2-dba-dt-{dt:.0e}.out"
    stdout = stdout_file.read_text()
    if result["status"] == "OK":
        if f"Richardson levels: {levels};" not in stdout:
            raise RuntimeError(f"Richardson setting not confirmed: {stdout_file}")
        measured = np.loadtxt(matrix_path(directory, dt))
        if measured.shape != (6, 6) or not np.all(np.isfinite(measured)):
            raise RuntimeError(f"invalid matrix: {directory}")
        exact = dba()
        error = (measured - exact) * COORDINATE_SCALES[None, :] / COORDINATE_SCALES[:, None]
        scaled_exact = exact * COORDINATE_SCALES[None, :] / COORDINATE_SCALES[:, None]
        row, column = np.unravel_index(np.argmax(np.abs(error)), error.shape)
        result.update(scaled_max_error=float(np.max(np.abs(error))),
                      relative_frobenius_error=float(np.linalg.norm(error) / np.linalg.norm(scaled_exact)),
                      worst_entry=f"R{row + 1}{column + 1}")
        match = re.search(r"OrbThreader\.+ Wall max =\s*([\d.eE+\-]+)", stdout)
        result["orbit_threader_wall_s"] = float(match[1]) if match else None
        for i, j in itertools.product(range(6), repeat=2):
            result[f"R{i + 1}{j + 1}"] = float(measured[i, j])
    elif not ("time step is too long" in stdout and "element 'QACH'" in stdout):
        # Expected coarse-step rejection is data; other errors must not silently pass.
        raise RuntimeError(f"unexpected OPALX failure: inspect {stdout_file}")

    artifacts = [directory / "map-2-dba.in", stdout_file]
    if result["status"] == "OK":
        artifacts.append(matrix_path(directory, dt))
    result["artifacts_sha256"] = {file.name: digest(file.read_bytes()) for file in artifacts}
    pending = directory / "result.pending.json"
    pending.write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
    pending.replace(result_file)
    return result | {"cached": False}


def summarize(output: Path, plots: bool = True) -> pd.DataFrame:
    records = [json.loads(file.read_text()) for file in output.glob("L*/dt-*/eps-*/ranks-*/result.json")]
    if not records:
        raise RuntimeError("no recorded cases")
    frame = pd.DataFrame([{key: value for key, value in item.items()
                           if key not in ("signature", "artifacts_sha256")}
                          for item in records]).sort_values(
                              ["ranks", "richardson_levels", "dt_s", "epsilon"])
    for column in ("scaled_max_error", "relative_frobenius_error", "R16", "R26"):
        if column not in frame:
            frame[column] = np.nan
    valid = frame[frame.status == "OK"]
    for (_, _), group in valid.groupby(["richardson_levels", "epsilon"]):
        for _, row in group.iterrows():
            rank_group = group[group.ranks == row.ranks]
            finest = rank_group.loc[rank_group.dt_s.idxmin()]
            measured = np.loadtxt(matrix_path(output / row.path, row.dt_s))
            reference = np.loadtxt(matrix_path(output / finest.path, finest.dt_s))
            frame.loc[row.name, "difference_from_finest_dt"] = np.max(np.abs(measured - reference))
            frame.loc[row.name, "R16_difference_from_finest_dt"] = measured[0, 5] - reference[0, 5]
            frame.loc[row.name, "R26_difference_from_finest_dt"] = measured[1, 5] - reference[1, 5]
            peers = group[group.dt_s == row.dt_s]
            first_rank = peers.loc[peers.ranks.idxmin()]
            rank_reference = np.loadtxt(matrix_path(output / first_rank.path, first_rank.dt_s))
            frame.loc[row.name, "mpi_matrix_max_difference"] = np.max(np.abs(measured - rank_reference))
            frame.loc[row.name, "mpi_diagnostics_max_difference"] = max(
                abs(row.determinant_error - first_rank.determinant_error),
                abs(row.canonical_J_error - first_rank.canonical_J_error))

    frame.to_csv(output / "results.csv", index=False, float_format="%.12e")
    primary = frame[(frame.ranks == 1) & (frame.status == "OK")]
    best = primary.loc[primary.groupby("richardson_levels").scaled_max_error.idxmin()]
    best.to_csv(output / "best-by-level.csv", index=False, float_format="%.12e")
    best_by_dt = primary.loc[primary.groupby(["richardson_levels", "dt_s"]).scaled_max_error.idxmin()]
    best_by_dt.to_csv(output / "best-by-dt.csv", index=False, float_format="%.12e")
    orders = []
    for variable, fixed in [("dt_s", "epsilon"), ("epsilon", "dt_s")]:
        for (level, fixed_value), group in primary.groupby(["richardson_levels", fixed]):
            group = group.sort_values(variable, ascending=False)
            for coarse, fine in zip(group.iloc[:-1].to_dict("records"),
                                    group.iloc[1:].to_dict("records")):
                for metric in ("scaled_max_error", "abs_R16", "abs_R26", "difference_from_finest_dt"):
                    if coarse[metric] > 0 and fine[metric] > 0:
                        orders.append(dict(variable=variable, fixed_variable=fixed, fixed_value=fixed_value,
                                           richardson_levels=level, coarse=coarse[variable], fine=fine[variable],
                                           metric=metric, order=np.log(coarse[metric] / fine[metric]) /
                                           np.log(coarse[variable] / fine[variable])))
    pd.DataFrame(orders).to_csv(output / "observed-orders.csv", index=False, float_format="%.12e")
    scripts = [Path(__file__), Path(__file__).with_name("convergence_grid_plots.py"),
               Path(__file__).with_name("validate_convergence_grid.py"),
               Path(__file__).with_name("test_convergence_grid.py")]
    (output / "analysis-provenance.json").write_text(json.dumps(
        {"postprocessing_scripts_sha256": {path.name: digest(path.read_bytes()) for path in scripts},
         "coordinate_scales": "(1 m, 1, 1 m, 1, 1 m, 1)",
         "numpy_version": np.__version__, "pandas_version": pd.__version__}, indent=2) + "\n")
    if plots and not primary.empty:
        from convergence_grid_plots import make_plots
        make_plots(frame, output)
    print(best[["richardson_levels", "dt_s", "epsilon", "R16", "R26", "scaled_max_error",
                "relative_frobenius_error"]].to_string(index=False), flush=True)
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=REPOSITORY / "omp-build/src/opalx")
    parser.add_argument("--input", type=Path, default=DBA / "map-2-dba.in")
    parser.add_argument("--output-dir", type=Path, default=DBA / "convergence-grid")
    parser.add_argument("--dt", type=float, nargs="+", default=DEFAULT_DT)
    parser.add_argument("--epsilon", type=float, nargs="+", default=DEFAULT_EPSILON)
    parser.add_argument("--levels", type=int, nargs="+", default=list(range(5)))
    parser.add_argument("--ranks", type=int, nargs="+", default=[1])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--mpi-args", default="--map-by slot:OVERSUBSCRIBE --bind-to none")
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    if args.summarize_only:
        summarize(args.output_dir, not args.no_plots)
        return 0
    if args.workers < 1 or any(rank < 1 for rank in args.ranks):
        parser.error("workers and MPI ranks must be positive")
    if any(level not in range(5) for level in args.levels):
        parser.error("Richardson levels must be in 0..4")
    if any(not np.isfinite(value) or value <= 0 for value in args.dt + args.epsilon):
        parser.error("DT and epsilon must be finite and positive")
    if any(epsilon >= 1 for epsilon in args.epsilon):
        parser.error("epsilon must be less than one to retain forward momentum")
    # The shared DT runner stores/parses these values at one significant digit.
    if any(float(f"{dt:.0e}") != dt for dt in args.dt):
        parser.error("DT values must have one significant digit")
    for values in (args.dt, args.epsilon, args.levels, args.ranks):
        if len(set(values)) != len(values):
            parser.error("duplicate values are not allowed")
    executable = args.executable.resolve()
    source = args.input.read_text()
    if not re.search(r"\bMAXSTEPS\s*=\s*1\s*[,;]", source):
        parser.error("the DBA template must use MAXSTEPS=1")
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OMP_PROC_BIND"] = "false"
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    provenance = identity(executable, source)
    manifest_path = output / "provenance.json"
    if manifest_path.exists():
        if json.loads(manifest_path.read_text())["identity"] != provenance:
            raise RuntimeError("source/executable/input changed: use a new output directory")
    else:
        manifest = dict(identity=provenance, executable=str(executable), input_path=str(args.input.resolve()),
                        input_template=source, created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        git_revision=subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                             cwd=REPOSITORY, text=True).strip(),
                        omp_num_threads=1, workers=args.workers, mpi_args=args.mpi_args,
                        coordinate_scales="(1 m, 1, 1 m, 1, 1 m, 1)",
                        numpy_version=np.__version__, pandas_version=pd.__version__)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "source-changes.patch").write_bytes(subprocess.check_output(
            ["git", "diff", "HEAD", "--", "src"], cwd=REPOSITORY))
        # Preserve complete source snapshots as well: git diff omits new/untracked headers.
        for source_path in SOURCE_PATHS:
            archived = output / "source-snapshot" / source_path
            archived.parent.mkdir(parents=True, exist_ok=True)
            archived.write_bytes((REPOSITORY / source_path).read_bytes())
        (output / "reference-particle.txt").write_bytes((MAP2 / "reference-particle.txt").read_bytes())
        np.savetxt(output / "analytic-map.txt", dba(), fmt="% .16e")
    cases = list(itertools.product(args.dt, args.levels, args.epsilon, args.ranks))
    request = {"dt_s": args.dt, "epsilon": args.epsilon, "levels": args.levels,
               "ranks": args.ranks, "workers": args.workers, "cases": len(cases)}
    with (output / "invocations.jsonl").open("a") as log:
        log.write(json.dumps(request) + "\n")
    started = time.perf_counter()
    print(f"Requested {len(cases)} configurations; {args.workers} independent workers", flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_one, executable, source, output, provenance,
                               dt, epsilon, levels, ranks, tuple(shlex.split(args.mpi_args)))
                   for dt, levels, epsilon, ranks in cases]
        for index, future in enumerate(as_completed(futures), 1):
            row = future.result()
            print(f"{index:3d}/{len(cases)} L{row['richardson_levels']} "
                  f"dt={row['dt_s']:.0e} eps={row['epsilon']:.0e} np={row['ranks']} "
                  f"{row['status']} error={row.get('scaled_max_error', np.nan):.3e} "
                  f"{'cached' if row['cached'] else str(round(row['elapsed_s'], 2)) + ' s'}", flush=True)
    if identity(executable, source) != provenance or args.input.read_text() != source:
        raise RuntimeError("source/executable/template changed during study")
    summarize(output, not args.no_plots)
    print(f"Completed in {time.perf_counter() - started:.1f} s; results: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
