#!/usr/bin/env python3
"""Fetch Merlin pair-4 field probes and analyze the mesh convergence locally."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_REMOTE = "/psi/home/adelmann/opalx-dev/track12-mesh-8da5b9e83-20260820"
DEFAULT_CASES = (
    "production_rect_1024x128_fixed400k",
    "production_rect_1536x192_fixed400k",
    "production_rect_2048x256_fixed400k",
)


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def fetch_case(host: str, remote_root: str, case: str, destination: Path) -> None:
    remote = f"{remote_root}/{case}"
    local = destination / case
    (local / "input").mkdir(parents=True, exist_ok=True)
    root_files = (
        "completed",
        "runtime_compute.txt",
        "run_manifest.txt",
        "opalx.log",
        "pair4_field_probe.in",
        "pair4_field_probe_c0.stat",
        "pair4_field_probe_c1.h5",
    )
    run(["ssh", host, "test", "-f", f"{remote}/completed"])
    run(
        ["scp"]
        + [f"{host}:{remote}/{filename}" for filename in root_files]
        + [str(local)]
    )
    run(
        [
            "scp",
            f"{host}:{remote}/input/pair4_probe.fromfile",
            str(local / "input"),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="merlin6")
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE)
    parser.add_argument("--case", action="append", dest="cases")
    parser.add_argument(
        "--destination", type=Path, default=HERE / "a100_first_kick_mesh"
    )
    parser.add_argument(
        "--fixed-primary",
        type=Path,
        default=HERE / "input" / "primary_fixed.fromfile",
    )
    args = parser.parse_args()
    cases = tuple(args.cases or DEFAULT_CASES)

    args.destination.mkdir(parents=True, exist_ok=True)
    for case in cases:
        fetch_case(args.host, args.remote_root.rstrip("/"), case, args.destination)

    command = [
        sys.executable,
        str(HERE / "scan_first_kick_fields.py"),
        "--analyze-only",
        "--output-dir",
        str(args.destination),
        "--fixed-primary",
        str(args.fixed_primary),
    ]
    for case in cases:
        command.extend(("--case", case))
    run(command)
    print(f"local results: {args.destination.resolve()}")


if __name__ == "__main__":
    main()
