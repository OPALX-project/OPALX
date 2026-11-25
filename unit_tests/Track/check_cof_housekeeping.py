#!/usr/bin/env python3
"""Exercise COF file replacement, logging levels and collective failures.

Run with the project's Python environment and --executable <build>/bin/opalx.
Use --mpiexec <launcher> --ranks 2 to exercise the same checks on two ranks.
Each case runs in a private temporary directory; logs are printed on failure.
"""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--mpiexec")
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()
    command = [str(args.executable.resolve())]
    if args.mpiexec:
        command = ([args.mpiexec, args.numproc_flag, str(args.ranks)]
                   + args.mpi_preflag + command + args.mpi_postflag)
    elif args.ranks != 1:
        parser.error("--mpiexec is required for multiple ranks")
    template = Path(__file__).with_name("cof-housekeeping.in").read_text()
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_PROC_BIND="false")
    details = ("Nominal closure position_m=", "Integrator=", "iteration damping=",
               "Coordinates (x,px,y,py):", "Phase samples=",
               "lambda[0]", "lambda[1]", "lambda[2]", "lambda[3]")
    optics = ("COF matrix and tunes:", "One-turn-matrix (x, px, y, py):", "Full tunes")

    with tempfile.TemporaryDirectory(prefix="opalx-cof-") as temporary:
        root = Path(temporary)

        def run(name, *, level=3, output="orbit.json", time=0, steps=5000,
                failure=None, directory=None, track=False):
            work = directory or root / name
            work.mkdir(exist_ok=True)
            text = template.replace("@TIME@", str(time)).replace("@MAXSTEPS@", str(steps))
            text = text.replace("@OUTPUT@", "" if output is None else f', OUTPUT="{output}"')
            text = text.replace("@TRACK@", "" if not track else """
FS0: FIELDSOLVER, TYPE=NONE, NX=8, NY=8, NZ=8;
TRACK, LINE=FFA_RING, BEAM=B1, INITIALORBIT=FFA_ORBIT,
       DT=1e-9, MAXSTEPS=20000;
    RUN, METHOD="PARALLEL", FIELDSOLVER=FS0, TURNS=1;
ENDTRACK;
""")
            (work / "case.in").write_text(text)
            result = subprocess.run(command + ["--info", str(level), "case.in"], cwd=work,
                                    env=env, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=90)
            log = result.stdout
            try:
                if failure:
                    assert result.returncode != 0, "failure returned success"
                    assert failure in log, f"missing diagnostic: {failure}"
                else:
                    assert result.returncode == 0, f"exit status {result.returncode}"
                    for token in details:
                        assert (token in log) == (level >= 5), f"level {level}: {token}"
                    for token in optics:
                        assert (token in log) == (level >= 4), f"level {level}: {token}"
                    assert "Stability is a numerical matrix diagnostic" not in log
                    assert "Matrix eigenanalysis alone:" not in log
                    if level >= 4:
                        heading = next(line for line in log.splitlines()
                                       if "COF matrix and tunes:" in line)
                        assert "Near-integer warning: no" in heading
                        assert log.count("Near-integer warning:") == 1
                        tune_lines = log.split("Full tunes (uncoupled planar):", 1)[1].splitlines()
                        horizontal = tune_lines[1].split(">", 1)[-1]
                        vertical = tune_lines[2].split(">", 1)[-1]
                        assert horizontal.strip().startswith("horizontal/radial")
                        assert vertical.strip().startswith("vertical")
                        assert horizontal.index("nur=") == vertical.index("nuy=")
                        float(horizontal.split("nur=", 1)[1])
                        float(vertical.split("nuy=", 1)[1])
                        matrix_lines = log.split("One-turn-matrix (x, px, y, py):", 1)[1].splitlines()
                        for coordinate, line in zip(("x", "px", "y", "py"), matrix_lines[2:6]):
                            line = line.split(">", 1)[-1].strip()
                            assert line.startswith(coordinate + " [")
                            values = line.split("[", 1)[1].split("]", 1)[0].split()
                            assert len(values) == 4 and all("e" in value for value in values)
                    if level >= 4:
                        assert log.count("COF matrix and tunes:") == 1, "duplicate rank output"
                    if track:
                        for token in ("Creating new file for data:", "OrbitThreader work (",
                                      "* OT reference:", "* OT segmentation:",
                                      "* OT maps:", "* OT output:"):
                            assert (token in log) == (level >= 5), f"level {level}: {token}"
                        assert list(work.rglob("*_ElementPositions.sdds")), (
                            "logging level suppressed element-position file output")
                assert not list(work.glob("*.tmp.*")), "temporary output was not removed"
            except AssertionError:
                print(f"FAILED {name}:\n{log}")
                raise
            return work

        def read_result(work):
            def reject_constant(value):
                raise AssertionError(f"non-JSON number: {value}")
            data = json.loads((work / "orbit.json").read_text(), parse_constant=reject_constant)
            assert data["converged"] is True
            assert len(data["coordinates"]) == 4
            assert len(data["matrix"]) == 4 and all(len(row) == 4 for row in data["matrix"])
            assert data["stability"] == "STABLE"
            assert data["full_tunes"] is not None
            return data

        work = run("overwrite")
        first = read_result(work)
        # Holding the old inode open distinguishes replacement from in-place truncation.
        with (work / "orbit.json").open() as old:
            original = old.read()
            run("overwrite-again", directory=work, time=1e-9)
            old.seek(0)
            assert old.read() == original, "existing inode was modified"
        second = read_result(work)
        # The input parser can round a decimal literal by a few ulps.
        assert first["time_s"] == 0 and math.isclose(second["time_s"], 1e-9, rel_tol=1e-15), (
            first["time_s"], second["time_s"])
        saved = (work / "orbit.json").read_bytes()
        run("failed-solve", directory=work, steps=1, failure="No directed return within MAXSTEPS")
        assert (work / "orbit.json").read_bytes() == saved, "failed solve replaced result"
        for level in (4, 5):
            assert read_result(run(f"level-{level}", level=level, track=True)) == first, (
                "logging level changed the result")
        omitted = run("omitted", output=None, track=True)
        assert not list(omitted.glob("*.json")), "omitted OUTPUT wrote a result"
        run("empty", output="", failure="String value expected")
        run("missing-parent", output="missing/orbit.json", failure="Cannot create temporary COF OUTPUT")
        if os.geteuid() != 0:  # Root bypasses directory permission checks.
            readonly = work / "readonly"
            readonly.mkdir()
            previous = readonly / "orbit.json"
            previous.write_bytes(saved)
            readonly.chmod(0o555)
            try:
                run("unwritable", directory=work, output="readonly/orbit.json",
                    failure="Cannot create temporary COF OUTPUT")
                assert previous.read_bytes() == saved, "write failure replaced result"
                assert not list(readonly.glob("*.tmp.*"))
            finally:
                readonly.chmod(0o755)
        blocked = root / "blocked"
        blocked.mkdir()
        (blocked / "orbit.json").mkdir()
        marker = blocked / "orbit.json" / "keep"
        marker.write_text("previous contents")
        run("rename-failure", directory=blocked, failure="Cannot replace COF OUTPUT")
        assert marker.read_text() == "previous contents"
    print(f"PASS: COF housekeeping on {args.ranks} rank(s)")


if __name__ == "__main__":
    main()
