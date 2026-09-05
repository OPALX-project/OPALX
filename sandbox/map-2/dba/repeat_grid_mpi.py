#!/usr/bin/env python3
"""Repeat exactly the archived MPI validation points with a new study executable.

Writes atomic per-case records, not summary CSVs, so this can run alongside the
primary sweep. Run convergence_grid.py --summarize-only after both finish.
"""
import argparse
import json
import os
from pathlib import Path
import shlex

from convergence_grid import identity, run_one


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    manifest = json.loads((args.output / "provenance.json").read_text())
    executable = Path(manifest["executable"])
    source = manifest["input_template"]
    provenance = identity(executable, source)
    if provenance != manifest["identity"]:
        raise RuntimeError("new-study source/executable mismatch")
    cases = [json.loads(path.read_text()) for path in args.baseline.glob("L*/dt-*/eps-*/ranks-*/result.json")]
    cases = sorted([case for case in cases if case["ranks"] > 1],
                   key=lambda row: (-row["dt_s"], row["epsilon"], row["richardson_levels"]))
    if not cases:
        raise RuntimeError("no baseline MPI validation cases")
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OMP_PROC_BIND"] = "false"
    for index, case in enumerate(cases, 1):
        row = run_one(executable, source, args.output, provenance, case["dt_s"], case["epsilon"],
                      case["richardson_levels"], case["ranks"], tuple(shlex.split(manifest["mpi_args"])))
        print(f"{index}/{len(cases)} {row['path']}: {row['status']}", flush=True)
    if identity(executable, source) != provenance:
        raise RuntimeError("numerical sources changed during MPI repeats")


if __name__ == "__main__":
    main()
