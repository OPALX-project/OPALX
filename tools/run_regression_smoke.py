#!/usr/bin/env python3
"""Run one OPALX regression smoke test without touching its source directory.

The selected regression-test directory is copied recursively below
``<build-dir>/smoke/<test-name>/<run-id>/work``. OPALX is launched only from
that private copy, and all generated output and logs remain next to it.

Unless ``--exit-only`` is used, assertions from ``<test-name>.rt`` are applied
to the generated ``<test-name>.stat`` using the same relative-error rule as
``.github/workflows/cpu-serial.yml``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Sequence, Tuple


EXIT_SUCCESS = 0
EXIT_TEST_FAILURE = 1
EXIT_RUNNER_ERROR = 2


class SmokeRunnerError(RuntimeError):
    """Report an invalid setup or an unreadable regression artifact."""


class SmokeHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    """Preserve CLI examples while also displaying argument defaults."""


def positive_integer(value: str) -> int:
    """Parse a strictly positive command-line integer."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def positive_float(value: str) -> float:
    """Parse a strictly positive command-line floating-point value."""
    parsed = float(value)
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def is_within(path: Path, directory: Path) -> bool:
    """Return whether ``path`` is ``directory`` or one of its descendants."""
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def utc_timestamp() -> str:
    """Create a sortable, collision-resistant run identifier."""
    now = dt.datetime.now(dt.timezone.utc)
    return now.strftime("%Y%m%dT%H%M%S.%fZ")


def write_json(path: Path, value: Dict[str, Any]) -> None:
    """Write a human-readable JSON artifact."""
    with path.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def parse_sdds(path: Path) -> Tuple[List[str], List[List[str]]]:
    """Read the SDDS columns and rows needed by the CI smoke checks."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise SmokeRunnerError(f"cannot read SDDS file {path}: {error}") from error

    columns: List[str] = []
    number_of_parameters = 0
    data_end_index = None
    section = None

    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("&") and stripped != "&end":
            section = stripped[1:]
            if section == "parameter":
                number_of_parameters += 1
        elif stripped == "&end":
            if section == "data":
                data_end_index = index + 1
            section = None
        elif section == "column" and stripped.startswith("name="):
            column = stripped.split("name=")[1].split(",")[0].strip().strip('"')
            columns.append(column)

    if data_end_index is None:
        raise SmokeRunnerError(f"no SDDS data section found in {path}")

    rows = [
        [value.strip() for value in line.strip().split("\t")]
        for line in lines[data_end_index + number_of_parameters :]
        if line.strip()
    ]
    if not rows:
        raise SmokeRunnerError(f"no SDDS data rows found in {path}")
    return columns, rows


def stat_value(
    columns: Sequence[str], rows: Sequence[Sequence[str]], column: str, aggregation: str
) -> float:
    """Return the CI-compatible last or average value for one SDDS column."""
    try:
        column_index = columns.index(column)
    except ValueError as error:
        raise SmokeRunnerError(f"SDDS output has no column named {column!r}") from error

    try:
        values = [float(row[column_index]) for row in rows]
    except (IndexError, ValueError) as error:
        raise SmokeRunnerError(f"invalid data in SDDS column {column!r}") from error

    if aggregation == "last":
        return values[-1]
    if aggregation == "avg":
        return sum(values) / len(values)
    raise SmokeRunnerError(f"unknown aggregation type: {aggregation}")


def parse_rt(path: Path) -> List[Tuple[str, str, float]]:
    """Parse ``stat \"column\" last|avg tolerance`` smoke assertions."""
    checks: List[Tuple[str, str, float]] = []
    try:
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, start=1):
                stripped = line.split("#")[0].strip()
                parts = stripped.split()
                if len(parts) >= 4 and parts[0] == "stat":
                    try:
                        tolerance = float(parts[3])
                    except ValueError as error:
                        raise SmokeRunnerError(
                            f"invalid tolerance in {path}:{line_number}: {parts[3]!r}"
                        ) from error
                    checks.append((parts[1].strip('"'), parts[2], tolerance))
    except OSError as error:
        raise SmokeRunnerError(f"cannot read assertion file {path}: {error}") from error
    return checks


def compare_stats(reference_path: Path, output_path: Path, rt_path: Path) -> Dict[str, Any]:
    """Apply all .rt assertions and return a serializable comparison report."""
    reference_columns, reference_rows = parse_sdds(reference_path)
    output_columns, output_rows = parse_sdds(output_path)
    checks = parse_rt(rt_path)

    results = []
    failed = False
    for column, aggregation, tolerance in checks:
        reference = stat_value(reference_columns, reference_rows, column, aggregation)
        output = stat_value(output_columns, output_rows, column, aggregation)
        relative_error = abs(output - reference) / (abs(reference) + 1.0e-30)
        passed = relative_error < tolerance
        failed = failed or not passed
        results.append(
            {
                "aggregation": aggregation,
                "column": column,
                "output": output,
                "passed": passed,
                "reference": reference,
                "relative_error": relative_error,
                "tolerance": tolerance,
            }
        )

    return {
        "checks": results,
        "number_of_output_rows": len(output_rows),
        "passed": not failed,
    }


def format_comparison(report: Dict[str, Any], test_name: str) -> str:
    """Format a comparison report like the repository CI output."""
    checks = report["checks"]
    lines = [
        f"Checking {len(checks)} assertion(s) from {test_name}.rt "
        f"({report['number_of_output_rows']} time steps)"
    ]
    for check in checks:
        status = "PASS" if check["passed"] else "FAIL"
        lines.append(
            f"  [{status}] {check['column']} ({check['aggregation']}): "
            f"ref={check['reference']:.15e}  got={check['output']:.15e}  "
            f"rel_err={check['relative_error']:.2e}  tol={check['tolerance']:.0e}"
        )
    return "\n".join(lines) + "\n"


def resolve_executable(value: str, description: str) -> Path:
    """Resolve an executable path or a command available through PATH."""
    candidate = Path(value).expanduser()
    if candidate.parent != Path(".") or candidate.is_absolute():
        resolved = candidate.resolve()
        if not resolved.is_file():
            raise SmokeRunnerError(f"{description} does not exist: {resolved}")
        if not os.access(resolved, os.X_OK):
            raise SmokeRunnerError(f"{description} is not executable: {resolved}")
        return resolved

    located = shutil.which(value)
    if located is None:
        raise SmokeRunnerError(f"{description} is not available through PATH: {value}")
    return Path(located).resolve()


def build_argument_parser() -> argparse.ArgumentParser:
    """Define the command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Copy one regression test to build-tree scratch space, run OPALX there, "
            "and optionally apply its .rt checks."
        ),
        epilog="""examples:
  python3 tools/run_regression_smoke.py \\
    --opalx ../build_refactor_take_two/src/opalx \\
    --test-dir ../regression-tests-x/RegressionTests/Drift-3-open-fromfile \\
    --build-dir ../build_refactor_take_two

  python3 tools/run_regression_smoke.py \\
    --opalx ../build_refactor_take_two/src/opalx \\
    --test-dir ../regression-tests-x/RegressionTests/SwissFEL-booster-SC-emittedfromfile \\
    --build-dir ../build_refactor_take_two --np 2 --exit-only
""",
        formatter_class=SmokeHelpFormatter,
    )
    parser.add_argument("--opalx", required=True, help="OPALX executable or command name")
    parser.add_argument(
        "--test-dir",
        required=True,
        type=Path,
        help="source directory for one existing RegressionTests test",
    )
    parser.add_argument(
        "--build-dir",
        required=True,
        type=Path,
        help="build tree below which smoke artifacts are stored",
    )
    parser.add_argument(
        "--test-name",
        help="test/input stem when it differs from the source directory name",
    )
    parser.add_argument("--np", type=positive_integer, default=1, help="number of MPI ranks")
    parser.add_argument(
        "--mpi-launcher", default="mpirun", help="MPI launcher executable or command name"
    )
    parser.add_argument(
        "--np-flag", default="-np", help="MPI launcher option used to set the rank count"
    )
    parser.add_argument(
        "--launcher-arg",
        action="append",
        default=[],
        help="extra MPI launcher argument; repeat as needed (use --launcher-arg=--flag)",
    )
    parser.add_argument("--info", type=int, default=2, help="value passed to OPALX --info")
    parser.add_argument(
        "--opalx-arg",
        action="append",
        default=[],
        help="extra OPALX argument; repeat as needed (use --opalx-arg=--flag)",
    )
    parser.add_argument(
        "--reference-stat",
        type=Path,
        help=(
            "alternate reference .stat file, for example a retained Stage 0 baseline; "
            "it is copied into the run artifacts before comparison"
        ),
    )
    parser.add_argument(
        "--exit-only",
        action="store_true",
        help="require only a zero OPALX exit code; do not parse .rt or compare .stat files",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        help="terminate the launcher if it exceeds this many seconds",
    )
    parser.add_argument(
        "--label",
        help="short artifact label prepended to the unique UTC run identifier",
    )
    return parser


def validate_label(label: str) -> str:
    """Keep user-provided artifact labels path-safe."""
    allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
    if not label or any(character not in allowed for character in label):
        raise SmokeRunnerError(
            "--label may contain only ASCII letters, digits, periods, underscores, and hyphens"
        )
    return label


def run(args: argparse.Namespace) -> int:
    """Copy, launch, compare, and retain one isolated smoke-test run."""
    test_dir = args.test_dir.expanduser().resolve()
    build_dir = args.build_dir.expanduser().resolve()
    if not test_dir.is_dir():
        raise SmokeRunnerError(f"test directory does not exist: {test_dir}")
    if is_within(build_dir, test_dir):
        raise SmokeRunnerError(
            f"build directory must not be inside the source test directory: {build_dir}"
        )

    test_name = args.test_name or test_dir.name
    if not test_name or Path(test_name).name != test_name:
        raise SmokeRunnerError(f"invalid test name: {test_name!r}")
    input_source = test_dir / f"{test_name}.in"
    if not input_source.is_file():
        raise SmokeRunnerError(f"test input does not exist: {input_source}")

    opalx = resolve_executable(args.opalx, "OPALX executable")
    mpi_launcher = resolve_executable(args.mpi_launcher, "MPI launcher")

    reference_source = None
    if args.reference_stat is not None:
        reference_source = args.reference_stat.expanduser().resolve()
        if not reference_source.is_file():
            raise SmokeRunnerError(f"reference stat file does not exist: {reference_source}")

    run_id = utc_timestamp()
    if args.label is not None:
        run_id = f"{validate_label(args.label)}-{run_id}"
    run_dir = build_dir / "smoke" / test_name / run_id
    work_dir = run_dir / "work"
    run_dir.mkdir(parents=True, exist_ok=False)

    manifest_path = run_dir / "run.json"
    manifest: Dict[str, Any] = {
        "build_directory": str(build_dir),
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "exit_only": args.exit_only,
        "mpi_ranks": args.np,
        "source_test_directory": str(test_dir),
        "status": "preparing",
        "test_name": test_name,
        "work_directory": str(work_dir),
    }
    write_json(manifest_path, manifest)
    print(f"Artifacts: {run_dir}")

    try:
        # Dereference symlinks so OPALX cannot write through a copied link into
        # the source checkout. All input and reference files become private.
        shutil.copytree(test_dir, work_dir, symlinks=False)
    except OSError as error:
        manifest.update({"status": "copy-error", "error": str(error)})
        write_json(manifest_path, manifest)
        raise SmokeRunnerError(f"cannot copy {test_dir} to {work_dir}: {error}") from error

    output_path = work_dir / f"{test_name}.stat"
    if output_path.exists():
        # A stale generated output may already be present in the source tree.
        # Remove only the private copy so this run cannot pass using old data.
        output_path.unlink()

    if reference_source is not None:
        reference_path = run_dir / "comparison-reference.stat"
        shutil.copy2(reference_source, reference_path)
    else:
        reference_path = work_dir / "reference" / f"{test_name}.stat"

    command = [
        str(mpi_launcher),
        *args.launcher_arg,
        args.np_flag,
        str(args.np),
        str(opalx),
        f"{test_name}.in",
        "--info",
        str(args.info),
        *args.opalx_arg,
    ]
    manifest.update(
        {
            "command": command,
            "opalx_executable": str(opalx),
            "reference_stat": None if args.exit_only else str(reference_path),
            "status": "running",
        }
    )
    write_json(manifest_path, manifest)
    print(f"Command: {shlex.join(command)}")

    stdout_path = run_dir / "opalx.stdout.log"
    stderr_path = run_dir / "opalx.stderr.log"
    try:
        with stdout_path.open("w", encoding="utf-8") as stdout_stream, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr_stream:
            completed = subprocess.run(
                command,
                cwd=work_dir,
                stdout=stdout_stream,
                stderr=stderr_stream,
                check=False,
                timeout=args.timeout,
            )
        return_code = completed.returncode
    except subprocess.TimeoutExpired:
        manifest.update(
            {
                "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "status": "timeout",
                "timeout_seconds": args.timeout,
            }
        )
        write_json(manifest_path, manifest)
        print(f"FAIL: OPALX exceeded the {args.timeout:g} second timeout", file=sys.stderr)
        return EXIT_TEST_FAILURE
    except OSError as error:
        manifest.update(
            {
                "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "error": str(error),
                "status": "launch-error",
            }
        )
        write_json(manifest_path, manifest)
        raise SmokeRunnerError(f"cannot launch smoke test: {error}") from error

    manifest["opalx_return_code"] = return_code
    if return_code != 0:
        manifest.update(
            {
                "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "status": "opalx-failed",
            }
        )
        write_json(manifest_path, manifest)
        print(f"FAIL: OPALX exited with status {return_code}", file=sys.stderr)
        return EXIT_TEST_FAILURE

    print("OPALX exited successfully")
    if args.exit_only:
        manifest.update(
            {
                "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "status": "passed-exit-only",
            }
        )
        write_json(manifest_path, manifest)
        return EXIT_SUCCESS

    try:
        rt_path = work_dir / f"{test_name}.rt"
        if not rt_path.is_file():
            raise SmokeRunnerError(f"assertion file does not exist in copied test: {rt_path}")
        if not reference_path.is_file():
            raise SmokeRunnerError(f"reference stat file does not exist: {reference_path}")
        if not output_path.is_file():
            raise SmokeRunnerError(f"OPALX did not produce expected stat output: {output_path}")
        comparison = compare_stats(reference_path, output_path, rt_path)
    except SmokeRunnerError as error:
        manifest.update(
            {
                "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "error": str(error),
                "status": "comparison-error",
            }
        )
        write_json(manifest_path, manifest)
        raise
    comparison_text = format_comparison(comparison, test_name)
    (run_dir / "comparison.log").write_text(comparison_text, encoding="utf-8")
    print(comparison_text, end="")

    manifest.update(
        {
            "comparison": comparison,
            "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "status": "passed" if comparison["passed"] else "comparison-failed",
        }
    )
    write_json(manifest_path, manifest)
    return EXIT_SUCCESS if comparison["passed"] else EXIT_TEST_FAILURE


def main() -> int:
    """Parse arguments and report expected runner errors without a traceback."""
    parser = build_argument_parser()
    args = parser.parse_args()
    try:
        return run(args)
    except SmokeRunnerError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_RUNNER_ERROR


if __name__ == "__main__":
    sys.exit(main())
