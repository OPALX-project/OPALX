#!/usr/bin/env python3
"""Run the three exact CI header-failure fixtures without editing their inputs.

Extract the selected RegressionTests directories from regression-tests-x
8fdc4502e2cac74f0bd13ba871d2a5023d80a59a into fixtures/ first. This runner
records completion, input hashes and finite SDDS physics output; it reports
reference differences without claiming to reproduce the remote CI comparator.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CASES = ("Dist-emittedfromfile", "AWAGun-1-emittedfromfile", "AWAGun-binned-p3m-emittedfromfile")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_stat(path):
    names, parameters, block, start = [], 0, None, None
    text = path.read_text().splitlines()
    for i, raw in enumerate(text):
        line = raw.strip()
        if line.startswith("&parameter"):
            parameters += 1
            block = "parameter"
        elif line.startswith("&column"):
            block = "column"
        elif line.startswith("&data"):
            block = "data"
        elif line.startswith("&end"):
            if block == "data":
                start = i+1
                break
            block = None
        elif block == "column":
            match = re.match(r"name\s*=\s*([^,\s]+)", line)
            if match:
                names.append(match[1])
    if start is None:
        raise ValueError(f"Missing SDDS data: {path}")
    rows = [line.split() for line in text[start:] if line.strip()][parameters:]
    return pd.DataFrame([[float(v) for v in row] for row in rows], columns=names)


def compare_runs(left, right):
    reports = []
    for stat in sorted(left.glob("n*/*/*.stat")):
        target = right / stat.relative_to(left)
        if not target.exists():
            continue
        a, b = read_stat(stat), read_stat(target)
        definition = stat.with_suffix(".rt")
        metrics = [line.split()[1].strip('"') for line in definition.read_text().splitlines()
                   if line.startswith("stat ")]
        clocks_match = len(a) == len(b) and np.array_equal(a.t, b.t)
        report = dict(case=stat.parent.name, ranks=stat.parent.parent.name,
                      left=str(stat), right=str(target), rows=[len(a), len(b)],
                      exact_clocks_match=clocks_match)
        if clocks_match:
            report["max_abs_difference"] = {m: float((a[m]-b[m]).abs().max()) for m in metrics}
        reports.append(report)
    if not reports:
        raise ValueError("No matching control statistics")
    (right / "old_binary_comparison.json").write_text(json.dumps(reports, indent=2)+"\n")
    print(json.dumps(reports, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx", type=Path)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--label", help="New output-directory name; never overwritten")
    parser.add_argument("--canonical-header", action="store_true",
                        help="Control: replace only the unit-annotated legacy header, never rows")
    parser.add_argument("--compare", type=Path, nargs=2, metavar=("OLD_RUN", "NEW_RUN"),
                        help="Compare saved physics statistics at identical output times; no tracking")
    args = parser.parse_args()
    if args.compare:
        compare_runs(*(p.resolve() for p in args.compare))
        return
    if args.opalx is None or args.label is None:
        parser.error("Runs require --opalx and a new --label")
    exe = args.opalx.resolve()
    root = HERE / "runs" / args.label
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_PROC_BIND="false", OPENBLAS_NUM_THREADS="1")
    reports = []
    for ranks in args.ranks:
        for case in args.cases:
            source = HERE / "fixtures" / case
            run = root / f"n{ranks}" / case
            shutil.copytree(source, run)
            inputs = {str(p.relative_to(source)): digest(p) for p in source.rglob("*") if p.is_file()}
            if args.canonical_header:
                # Generated control inputs only; the original fixture is never edited.
                data = run / ("awagun-flattop.dat" if case == "Dist-emittedfromfile" else "opal-AWAGun-1_DIST.dat")
                lines = data.read_bytes().splitlines(keepends=True)
                if not lines[1].lower().startswith(b"x [m]"):
                    raise ValueError("Unexpected legacy header in control fixture")
                lines[1] = b"x px y py t pz bin\n"
                data.write_bytes(b"".join(lines))
            command = ["mpiexec", "--oversubscribe", "--bind-to", "none", "-n", str(ranks),
                       str(exe), case+".in", "--info", "4"]
            print(f"Running n{ranks} {case}", flush=True)
            started = time.monotonic()
            with (run / "opalx.log").open("w") as log:
                result = subprocess.run(command, cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT)
            report = dict(case=case, ranks=ranks, exit_code=result.returncode,
                          wall_seconds=time.monotonic()-started, command=command,
                          test_revision="8fdc4502e2cac74f0bd13ba871d2a5023d80a59a",
                          executable_sha256=digest(exe), input_sha256=inputs,
                          canonical_header_control=args.canonical_header,
                          executed_input_sha256={p: digest(run / p) for p in inputs},
                          invocation_checkout_revision=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip())
            if result.returncode == 0:
                current = read_stat(run / (case+".stat"))
                reference = read_stat(run / "reference" / (case+".stat"))
                metrics = [line.split()[1].strip('"') for line in (run / (case+".rt")).read_text().splitlines() if line.startswith("stat ")]
                if not np.isfinite(current[metrics].to_numpy()).all():
                    raise ValueError(f"Nonfinite SDDS physics metrics: {run}")
                report.update(stat_rows=len(current), last_time_s=float(current.t.iloc[-1]),
                              finite_physics=True, final_particles=float(current.numParticles.iloc[-1]))
                final = pd.DataFrame([dict(metric=m, current=current[m].iloc[-1], reference=reference[m].iloc[-1],
                                            absolute_difference=abs(current[m].iloc[-1]-reference[m].iloc[-1])) for m in metrics])
                final.to_csv(run / "final_reference_differences.csv", index=False)
                report["reference_scope"] = "Final-value differences only, not the CI avg/last acceptance comparator; no tolerances changed"
            if any(digest(source / p) != expected for p, expected in inputs.items()):
                raise ValueError("Original fixture changed")
            (run / "manifest.json").write_text(json.dumps(report, indent=2)+"\n")
            reports.append(report)
            (root / "summary.json").write_text(json.dumps(reports, indent=2)+"\n")
            print(f"Finished n{ranks} {case}: exit={result.returncode}, {report['wall_seconds']:.2f}s", flush=True)
    if any(r["exit_code"] != 0 for r in reports):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
