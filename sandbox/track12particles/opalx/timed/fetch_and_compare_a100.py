#!/usr/bin/env python3
"""Fetch a completed Merlin track12 run and make the CAIN comparison locally."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_REMOTE = (
    "/psi/home/adelmann/opalx-dev/track12-a1cdd08-explicit-20260819/"
    "run/production"
)


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def fetch(host: str, remote_dir: str, destination: Path) -> None:
    run(["ssh", host, "test", "-f", f"{remote_dir}/completed"])
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "results").mkdir(exist_ok=True)
    (destination / "input").mkdir(exist_ok=True)

    root_files = [
        "completed",
        "runtime_compute.txt",
        "run_manifest.txt",
        "opalx.log",
        "track12_timed.in",
        "track12_timed_c1.h5",
        "track12_timed_c2.h5",
    ]
    run(
        ["scp"]
        + [f"{host}:{remote_dir}/{name}" for name in root_files]
        + [str(destination)]
    )
    run(
        [
            "scp",
            f"{host}:{remote_dir}/results/preparation_manifest.json",
            str(destination / "results"),
        ]
    )
    input_files = [
        "primary_fixed_metadata.json",
        "track12_electrons.emittedfromfile",
        "track12_positrons.emittedfromfile",
    ]
    run(
        ["scp"]
        + [f"{host}:{remote_dir}/input/{name}" for name in input_files]
        + [str(destination / "input")]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="merlin6")
    parser.add_argument("--remote-dir", default=DEFAULT_REMOTE)
    parser.add_argument(
        "--destination", type=Path, default=HERE / "a100_400k_1024x128x128"
    )
    parser.add_argument("--no-compare", action="store_true")
    args = parser.parse_args()

    fetch(args.host, args.remote_dir, args.destination)
    if not args.no_compare:
        run(
            [
                sys.executable,
                str(HERE / "compare_timed_track12.py"),
                "--run-dir",
                str(args.destination),
            ]
        )
    print(f"local run: {args.destination.resolve()}")


if __name__ == "__main__":
    main()
