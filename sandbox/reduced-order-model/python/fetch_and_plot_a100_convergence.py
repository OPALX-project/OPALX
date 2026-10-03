#!/usr/bin/env python3
"""Fetch compact Merlin A100 convergence results and make local plots."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path, PurePosixPath


SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
ANALYZER = SCRIPT_DIR / "analyze_opalx_convergence.py"
DEFAULT_OUTPUT_ROOT = MODEL_DIR / "outputs" / "a100-convergence"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remote-host", default="merlin6")
    parser.add_argument(
        "--remote-run-root",
        required=True,
        help="Absolute Merlin study directory printed by submit_a100_convergence.sh.",
    )
    parser.add_argument(
        "--local-root",
        type=Path,
        default=None,
        help="Local destination (default: outputs/a100-convergence/<remote-name>).",
    )
    parser.add_argument("--skip-fetch", action="store_true")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Plot completed cases before the final four-rank job has finished.",
    )
    return parser.parse_args()


def fetch_results(host: str, remote_root: str, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    include_patterns = (
        "/study_manifest.csv",
        "/study_manifest.json",
        "/study_metadata.json",
        "/fixed_primary_metadata.json",
        "/opalx_source.patch",
        "/study-completed",
        "/rank-study-completed",
        "/cases/",
        "/cases/*/",
        "/cases/*/3sigma/",
        "/cases/*/3sigma/*.in",
        "/cases/*/3sigma/*.fromfile",
        "/cases/*/3sigma/*_c0.stat",
        "/cases/*/3sigma/*_c1.h5",
        "/cases/*/3sigma/runtime_compute.txt",
        "/cases/*/3sigma/opalx.log",
        "/cases/*/3sigma/completed",
    )
    # Copy content only. Preserving remote ownership/modes is unnecessary for
    # analysis and is rejected by some managed local workspaces.
    command = ["rsync", "-rz", "--prune-empty-dirs", "--partial"]
    for pattern in include_patterns:
        command.extend(("--include", pattern))
    command.extend(("--exclude", "*", f"{host}:{remote_root.rstrip('/')}/", f"{destination}/"))
    print("fetching compact study output with rsync", flush=True)
    subprocess.run(command, check=True)


def main() -> int:
    args = parse_args()
    remote_name = PurePosixPath(args.remote_run_root.rstrip("/")).name
    if not remote_name:
        raise ValueError("--remote-run-root must name a study directory")
    local_root = args.local_root or (DEFAULT_OUTPUT_ROOT / remote_name)
    raw_dir = local_root / "raw"
    plot_dir = local_root / "plots"

    if not args.skip_fetch:
        fetch_results(args.remote_host, args.remote_run_root, raw_dir)
    if not raw_dir.is_dir():
        raise FileNotFoundError(raw_dir)
    completed = (raw_dir / "study-completed").exists() or (
        raw_dir / "rank-study-completed"
    ).exists()
    if not completed and not args.allow_incomplete:
        raise RuntimeError(
            "the remote study is not complete; wait for the final four-rank job "
            "or pass --allow-incomplete"
        )

    command = [
        sys.executable,
        str(ANALYZER),
        "--case-root",
        str(raw_dir / "cases"),
        "--output-dir",
        str(plot_dir),
    ]
    subprocess.run(command, check=True)
    print(f"plots and CSV summaries: {plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
