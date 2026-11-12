#!/usr/bin/env python3
"""Cyclotron PROBE integration/passivity checks using maintained cyclotron inputs.

Run with ~/.venv-h6/bin/python tools/check_cyclotron_probe.py --help.
Runs are separate from unit tests; all transport cases use one MPI rank.
No numeric regression references are generated or changed. The coasting and RF
fixtures supply transport; analytic crossing accuracy is checked by C++ units.
"""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import h5py
import numpy as np


def run(executable, fixture, target, probes, ascii_output=False):
    """Run a fresh copied fixture and return its final phase-space arrays."""
    shutil.copytree(fixture, target)
    source = target / (fixture.name + ".in")
    text = source.read_text()
    text = text.replace("PSDUMPFREQ=0", "PSDUMPFREQ=1000000000")
    if fixture.name == "Cyclotron-Coasting":
        # A single turn misses periodic field-index drift; check ten unchanged-dt turns.
        text = text.replace("MAXSTEPS=2000", "MAXSTEPS=20000").replace("TURNS=1", "TURNS=10")
    if ascii_output:
        text = text.replace("OPTION, VERSION=", "OPTION, ASCIIDUMP=TRUE;\nOPTION, VERSION=", 1)
    if probes:
        definitions = (
            'P1: PROBE, XSTART=-10000, XEND=10000, YSTART=1000, YEND=1000, STEP=100, OUTFN="hits";\n'
            'P2: PROBE, XSTART=10000, XEND=-10000, YSTART=1000, YEND=1000, STEP=100, OUTFN="reverse";\n'
            'P3: PROBE, XSTART=-10000, XEND=10000, YSTART=1000000, YEND=1000000, STEP=100, OUTFN="miss";\n'
        )
        text = text.replace("MYCYCL: RING=(", definitions + "MYCYCL: RING=(P1,P2,P3,", 1)
    source.write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS="2", IPPL_OPENMP_THREADS="2")
    command = ["mpirun", "-np", "1", str(executable), source.name, "--info", "2"]
    with (target / "run.log").open("w") as log:
        subprocess.run(command, cwd=target, env=env, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=300)
    files = [path for path in target.glob("*.h5") if path.stem == fixture.name]
    assert len(files) == 1, f"Expected phase-space file in {target}"
    with h5py.File(files[0]) as file:
        steps = sorted((key for key in file if key.startswith("Step#")),
                       key=lambda key: int(key.split("#")[1]))
        final = file[steps[-1]]
        arrays = {key: final[key][:] for key in ("id", "x", "y", "z", "px", "py", "pz")}
        for key in ("GlobalTrackStep", "TIME", "SPOS", "RefPartR", "RefPartP"):
            arrays[key] = final.attrs[key]
    if fixture.name == "Cyclotron-Coasting":
        separation = np.linalg.norm(np.column_stack([arrays[key] for key in ("x", "y", "z")]), axis=1)
        # 1 nm bounds accumulated roundoff for this zero-offset proton, far below
        # the 0.2 mm first missed-kick excursion of the old periodic candidate map.
        assert separation.max() < 1e-9, f"Coasting particle/reference separation: {separation} m"
        assert "completed directed turn 10" in (target / "run.log").read_text()
    return arrays


def check_outputs(directory, ascii_output):
    """Check hit clipping, suffixes, histogram units/counts and reversed geometry."""
    extension = ".loss" if ascii_output else ".h5"
    assert not (directory / ("miss" + extension)).exists()
    assert np.loadtxt(directory / "miss.hist").sum() == 0
    assert not [line for line in (directory / "miss.peaks").read_text().splitlines()
                if line and not line.startswith("#")]
    counts = []
    for name in ("hits", "reverse"):
        if ascii_output:
            data = np.atleast_2d(np.loadtxt(directory / (name + extension)))
            assert data.shape[1] == 10
            ids, times, x, z = data[:, 6], data[:, 9], data[:, 0], data[:, 2]
            # Existing ASCII loss writer uses six significant digits.
            np.testing.assert_allclose(z, 1, rtol=0, atol=1e-6)
        else:
            with h5py.File(directory / (name + extension)) as file:
                groups = [file[key] for key in file if key.startswith("Step#")]
                data = {key: np.concatenate([group[key][:] for group in groups])
                        for key in ("id", "time", "x", "z")}
            ids, times, x, z = (data[key] for key in ("id", "time", "x", "z"))
            np.testing.assert_allclose(z, 1, rtol=0, atol=1e-13)
        assert len(times) >= 1, directory
        assert np.all(times > 0) and len(np.unique(times)) == len(times)
        assert np.all(np.abs(x) <= 10)
        hist = np.loadtxt(directory / (name + ".hist"))
        assert hist.sum() == len(times)
        header = (directory / (name + ".hist")).read_text().splitlines()[0]
        assert "1000 mm" in header and "100 mm" in header
        peaks = np.atleast_1d(np.loadtxt(directory / (name + ".peaks")))
        radii = np.hypot(x, z) * 1000
        assert np.all(peaks >= radii.min() - .01) and np.all(peaks <= radii.max() + .01)
        counts.append((ids, times))
    np.testing.assert_array_equal(counts[0][0], counts[1][0])
    np.testing.assert_array_equal(counts[0][1], counts[1][1])
    return len(counts[0][0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--regression-root", type=Path, required=True,
                        help="regression-tests-x checkout")
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    executable = args.executable.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    report = {"executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(), "cases": []}
    for fixture_name in ("Cyclotron-Coasting", "Cyclotron-RF-Stop"):
        fixture = args.regression_root / "RegressionTests" / fixture_name
        off = run(executable, fixture, args.output / (fixture_name + "-off"), False)
        for ascii_output in (False, True):
            target = args.output / (fixture_name + ("-ascii" if ascii_output else "-hdf5"))
            on = run(executable, fixture, target, True, ascii_output)
            for key in off:
                np.testing.assert_array_equal(on[key], off[key], err_msg=f"PROBE changed {key}")
            hits = check_outputs(target, ascii_output)
            if fixture_name == "Cyclotron-Coasting":
                # The horizontal chord crosses the orbit twice each complete turn.
                assert hits == 20, f"Expected twenty coasting crossings, got {hits}"
            report["cases"].append({"fixture": fixture_name, "ascii": ascii_output,
                                    "hits": hits, "phase_space_identical": True, "time_steps_reference_identical": True,
                                    "input_sha256": hashlib.sha256(
                                        (target / (fixture_name + ".in")).read_bytes()).hexdigest()})
    # Explicit rejection contracts: these are integration cases, not unit tests.
    fixture = args.regression_root / "RegressionTests" / "Cyclotron-Coasting"
    text = (args.output / "Cyclotron-Coasting-hdf5" / "Cyclotron-Coasting.in").read_text()
    followup = text[text.index("TRACK,"):text.index("QUIT;")]
    noncyclotron = re.sub(r"MYCYCL: RING=.*?;",
                         "D0: DRIFT, L=1, ELEMEDGE=0; MYCYCL: LINE=(D0,P1);", text)
    noncyclotron = re.sub(r", TURNS=\d+", "", noncyclotron)
    rejections = {
        "step-zero": (text.replace("STEP=100", "STEP=0", 1), "positive STEP"),
        "degenerate": (text.replace("XEND=10000", "XEND=-10000", 1), "distinct endpoints"),
        "noncyclotron": (noncyclotron, "requires cyclotron-sector"),
        "reused-output": (text.replace("QUIT;", followup + "QUIT;"), "PROBE output already exists"),
    }
    for name, (input_text, expected) in rejections.items():
        target = args.output / ("reject-" + name)
        shutil.copytree(fixture, target)
        source = target / "Cyclotron-Coasting.in"
        source.write_text(input_text)
        result = subprocess.run(
            ["mpirun", "-np", "1", str(executable), source.name, "--info", "2"],
            cwd=target, env=dict(os.environ, OMP_NUM_THREADS="2"),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300)
        (target / "run.log").write_text(result.stdout)
        assert result.returncode != 0 and expected in result.stdout, (name, result.stdout[-2000:])
    report["rejections_checked"] = list(rejections)
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
