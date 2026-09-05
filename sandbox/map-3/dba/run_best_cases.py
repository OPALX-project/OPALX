#!/usr/bin/env python3
"""Run exactly one sequential single-rank map-3 case per map integrator.

Select the minimum recorded full-map error from all valid map-2 observations.
This is a transfer of settings, not a map-3 convergence or matching claim.
Do not use the old hard-edge oracle to label map-3 dispersion as an error.
"""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import operator
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
MAP2 = REPO / "sandbox/map-2/dba"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-map3-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/opalx-map3-cache")
sys.path.insert(0, str(MAP2))
from convergence_perturbations import configure_map
from convergence_dt import parse_map, parse_reported_residual
from plot_field_history import read_element_spans
from analytic_fields import layout, load_parameters, np, pd

METHODS = ("BORIS", "RK4", "DOP853")
SELECTION = MAP2 / "integrator-comparison-stopped/all-completed.csv"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def numeric(expression, symbols):
    """Evaluate only arithmetic and known REAL names; never eval input text."""
    node = ast.parse(expression.strip(), mode="eval").body
    binary = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
              ast.Div: operator.truediv}
    def visit(item):
        if isinstance(item, ast.Constant) and type(item.value) in (int, float):
            return item.value
        if isinstance(item, ast.Name):
            return symbols[item.id]
        if isinstance(item, ast.UnaryOp) and isinstance(item.op, (ast.UAdd, ast.USub)):
            return visit(item.operand) * (-1 if isinstance(item.op, ast.USub) else 1)
        if isinstance(item, ast.BinOp) and type(item.op) in binary:
            return binary[type(item.op)](visit(item.left), visit(item.right))
        raise ValueError(f"Unsupported expression: {expression}")
    return float(visit(node))


def validate_physics(text, p):
    """Check every placed element against the already approved analytic model."""
    text = re.sub(r"//[^\n]*", "", text)
    symbols = {"PI": np.pi}
    for name, value in re.findall(r"\bREAL\s+(\w+)\s*=\s*([^;]+);", text):
        symbols[name] = numeric(value, symbols)
    np.testing.assert_allclose(symbols["P0"], p.momentum_GeV_c, rtol=0, atol=1e-15)
    expected = layout(p).set_index("name")
    matches = re.findall(r"\b(\w+)\s*:\s*(SBEND|DRIFT|QUADRUPOLE)\s*,([^;]+);", text)
    if len(matches) != 5 or {m[0] for m in matches} != set(expected.index):
        raise ValueError("Expected exactly B1, D1, QACH, D2, B2")
    for name, kind, attributes in matches:
        values = {key.strip(): numeric(value, symbols)
                  for key, value in (pair.split("=", 1) for pair in attributes.split(","))}
        row = expected.loc[name]
        wanted = {"L": row.body_end_m - row.body_start_m, "ELEMEDGE": row.body_start_m}
        if row.kind == "bend":
            assert kind == "SBEND"
            wanted |= {"ANGLE": np.deg2rad(p.body_angle_deg), "HGAP": p.hgap_m, "FINT": p.fint}
        elif row.kind == "quadrupole":
            assert kind == "QUADRUPOLE"
            wanted["K1"] = p.quadrupole_k1_seed_m2
        else:
            assert kind == "DRIFT"
        if values.keys() != wanted.keys():
            raise ValueError(f"Unexpected/missing physical attributes in {name}")
        np.testing.assert_allclose(list(values.values()), [wanted[k] for k in values], rtol=0, atol=2e-14)
    stop = numeric(re.search(r"\bZSTOP\s*=\s*([^,;]+)", text).group(1), symbols)
    np.testing.assert_allclose(stop, expected.support_end_m.max(), rtol=0, atol=2e-14)
    assert re.search(r"MAXSTEPS\s*=\s*1\s*[,;]", text)
    assert re.search(r"PARTICLE\s*=\s*ELECTRON", text)
    assert re.search(r"\bCHARGE\s*=\s*-1\s*[,;]", text)


def select_best(reselect=False):
    """Replay the published choices without requiring the incomplete map-2 archive.

    Explicit reselection retains the original archive audit. The default freezes
    the three settings as actually run; it makes no new optimality claim.
    """
    if not reselect:
        manifest = json.loads((ROOT / "opalx-best-settings/provenance.json").read_text())
        settings = manifest["settings"]
        if [item["integration_method"] for item in settings] != list(METHODS):
            raise ValueError("Expected one archived setting per integrator")
        for item in settings:
            result = json.loads((ROOT / "opalx-best-settings" /
                                 item["integration_method"].lower() / "result.json").read_text())
            for key in ("integration_method", "dt_s", "epsilon", "richardson_levels",
                        "previous_result_sha256"):
                if item[key] != result[key]:
                    raise ValueError(f"Archived setting mismatch: {key}")
            for name, sha in result["artifacts_sha256"].items():
                if digest(ROOT / "opalx-best-settings" / item["integration_method"].lower() / name) != sha:
                    raise ValueError(f"Archived artifact mismatch: {name}")
        return settings
    table = pd.read_csv(SELECTION, float_precision="round_trip")
    valid = table[(table.status == "OK") & (table.ranks == 1)]
    best = valid.loc[valid.groupby("integration_method").scaled_max_error.idxmin()].set_index("integration_method")
    selections = []
    for method in METHODS:
        row = best.loc[method]
        study = Path(row.study)
        # Historical tables retain the original machine's absolute checkout.
        case = MAP2 / study.name / row.path
        record = json.loads((case / "result.json").read_text())
        assert record["integration_method"] == method and record["status"] == "OK" and record["ranks"] == 1
        for key in ("dt_s", "epsilon", "richardson_levels", "scaled_max_error"):
            assert record[key] == row[key], (method, key)
        for filename, sha in record["artifacts_sha256"].items():
            assert digest(case / filename) == sha, case / filename
        selections.append(dict(integration_method=method, dt_s=record["dt_s"], epsilon=record["epsilon"],
                               richardson_levels=record["richardson_levels"],
                               previous_scaled_max_error=record["scaled_max_error"],
                               previous_case=str(case), previous_result_sha256=digest(case / "result.json")))
    return selections


def parse_result(case, selection, elapsed, returncode):
    stdout = (case / "map-3-dba.out").read_text()
    result = selection | dict(elapsed_s=elapsed, returncode=returncode, ranks=1,
                              status="FAILED", analytic_map_available=False)
    if returncode != 0:
        result["reason"] = "OPALX exited unsuccessfully; inspect map-3-dba.out"
        return result
    assert f'Linear-map integrator: {selection["integration_method"]};' in stdout
    assert f'Richardson levels: {selection["richardson_levels"]};' in stdout
    matrix = parse_map(stdout)
    if not np.isfinite(matrix).all():
        raise ValueError("Nonfinite map")
    np.savetxt(case / "matrix.txt", matrix, fmt="%.12e",
               header="OPALX map in (x,x',y,y',zeta,delta); NOT an analytic reference")
    result |= dict(status="OK", R16=float(matrix[0, 5]), R26=float(matrix[1, 5]),
                   determinant_error=parse_reported_residual(stdout, "|det(M) - 1|"),
                   canonical_J_error=parse_reported_residual(stdout, "max_ij |(M^T J M - J)_ij|"))
    for i in range(6):
        for j in range(6):
            result[f"R{i+1}{j+1}"] = float(matrix[i, j])
    timing = re.search(r"\bOrbThreader\s+([0-9.eE+\-]+)", stdout)
    if timing:
        result["orbit_threader_wall_s"] = float(timing.group(1))
    sdds = case / "data/map-3-dba_ElementPositions.sdds"
    data = re.search(r"&data\b.*?&end(.*)", sdds.read_text(), flags=re.S).group(1).strip().splitlines()
    assert int(data[0]) == 1, "SDDS must confirm a single rank"
    spans = read_element_spans(sdds)
    spans.to_csv(case / "element-spans.csv", index=False, float_format="%.12e")
    result["max_recorded_overlap_lanes"] = int(spans.lane.max() + 1)
    result["visited_elements"] = sorted(spans.element.unique())
    result["all_five_elements_visited"] = set(spans.element) == {"B1", "D1", "QACH", "D2", "B2"}
    result["recorded_support_end_m"] = float(spans.s_end_m.max())
    result["artifacts_sha256"] = {str(path.relative_to(case)): digest(path) for path in
                                 [case / "map-3-dba.in", case / "map-3-dba.out", case / "matrix.txt",
                                  sdds, case / "data/map-3-dba_DesignPath.dat", case / "timing.dat"]}
    return result


def report(output):
    results = [json.loads((output / method.lower() / "result.json").read_text()) for method in METHODS]
    pd.DataFrame(results).drop(columns=["artifacts_sha256"], errors="ignore").to_csv(
        output / "results.csv", index=False, float_format="%.12e")
    differences = []
    for i, left in enumerate(METHODS):
        for right in METHODS[i + 1:]:
            if any(row["status"] != "OK" for row in results if row["integration_method"] in (left, right)):
                continue
            a, b = (np.loadtxt(output / method.lower() / "matrix.txt") for method in (left, right))
            differences.append(dict(left=left, right=right, max_entry_difference=float(abs(a-b).max()),
                                    R16_difference=float(a[0, 5]-b[0, 5]), R26_difference=float(a[1, 5]-b[1, 5])))
    pd.DataFrame(differences).to_csv(output / "pairwise-differences.csv", index=False, float_format="%.12e")
    columns = ["integration_method", "status", "dt_s", "epsilon", "richardson_levels", "R16", "R26",
               "determinant_error", "canonical_J_error", "elapsed_s", "max_recorded_overlap_lanes"]
    frame = pd.DataFrame(results)
    print(frame[[c for c in columns if c in frame]].to_string(index=False), flush=True)
    return all(result["status"] == "OK" and result.get("all_five_elements_visited") for result in results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=REPO / "omp-build/src/opalx")
    parser.add_argument("--output", type=Path, default=ROOT / "opalx-best-settings")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--reselect-map2", action="store_true",
                        help="Reselect from the optional complete local map-2 campaign archive")
    args = parser.parse_args()
    output = args.output.resolve()
    if args.report_only:
        return 0 if report(output) else 1
    if output.exists():
        raise RuntimeError("Preserve existing runs: use --report-only or a new --output")
    source = ROOT / "map-3-dba.in"
    distribution = ROOT.parent / "reference-particle.txt"
    template = source.read_text()
    p = load_parameters()
    validate_physics(template, p)
    assert distribution.read_bytes() == (MAP2.parent / "reference-particle.txt").read_bytes()
    selections = select_best(args.reselect_map2)
    executable = args.executable.resolve()
    output.mkdir(parents=True)
    (output / "reference-particle.txt").write_bytes(distribution.read_bytes())
    sources = [source, distribution, ROOT / "parameters.json", ROOT / "analytic_fields.py", Path(__file__),
               *[REPO / f"src/Algorithms/{name}" for name in
                 ("OrbitThreader.cpp", "ExternalFieldRayTracker.cpp", "ExternalFieldRayTracker.h",
                  "RungeKuttaTableau.h", "LinearTransferMapBuilder.cpp")],
               REPO / "src/AbsBeamline/BendFieldModel.h", REPO / "src/AbsBeamline/SBend.cpp"]
    selection_source = SELECTION if args.reselect_map2 else ROOT / "opalx-best-settings/provenance.json"
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), ranks=1, omp_threads=1,
                      sequential=True, selection_rule=("minimum full-map error among saved map-2 cases"
                          if args.reselect_map2 else "frozen settings from published map-3 runs"),
                      selection_table=str(selection_source), selection_table_sha256=digest(selection_source),
                      settings=selections, executable=str(executable), executable_sha256=digest(executable),
                      source_hashes={str(path.relative_to(REPO)): digest(path) for path in sources})
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    env = os.environ | {"OMP_NUM_THREADS": "1", "OMP_PROC_BIND": "false"}
    for setting in selections:
        case = output / setting["integration_method"].lower()
        case.mkdir()
        text = configure_map(template, setting["epsilon"], setting["richardson_levels"], setting["integration_method"])
        text, count = re.subn(r"\bDT\s*=\s*[^,;]+", f'DT={setting["dt_s"]:.16e}', text)
        assert count == 1
        validate_physics(text, p)
        (case / "map-3-dba.in").write_text(text)
        command = [str(executable), "--info", "2", "map-3-dba.in"]
        print(f'Starting {setting["integration_method"]}: DT={setting["dt_s"]}, epsilon={setting["epsilon"]}, L={setting["richardson_levels"]}', flush=True)
        start = time.monotonic()
        with (case / "map-3-dba.out").open("w") as stream:
            completed = subprocess.run(command, cwd=case, env=env, stdout=stream, stderr=subprocess.STDOUT)
        result = parse_result(case, setting, time.monotonic() - start, completed.returncode)
        result["command"] = command
        (case / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f'{setting["integration_method"]}: {result["status"]}, {result["elapsed_s"]:.2f} s', flush=True)
    assert digest(executable) == provenance["executable_sha256"], "Executable changed during run"
    return 0 if report(output) else 1


if __name__ == "__main__":
    raise SystemExit(main())
