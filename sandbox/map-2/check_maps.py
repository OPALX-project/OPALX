#!/usr/bin/env python3
"""Run the map-2 inputs and compare their printed maps with analytic optics."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
P = 490.23677597553325
GAMMA2 = 1.0 + P * P


def drift(length: float) -> np.ndarray:
    result = np.eye(6)
    result[0, 1] = length
    result[2, 3] = length
    result[4, 5] = length / GAMMA2
    return result


def plane(k: float, length: float) -> np.ndarray:
    root = math.sqrt(abs(k))
    phase = root * length
    if k > 0.0:
        return np.array([[math.cos(phase), math.sin(phase) / root],
                         [-root * math.sin(phase), math.cos(phase)]])
    return np.array([[math.cosh(phase), math.sinh(phase) / root],
                     [root * math.sinh(phase), math.cosh(phase)]])


def quadrupole(k: float, length: float) -> np.ndarray:
    result = np.eye(6)
    # OPALX K1>0 defocuses horizontally and focuses vertically.
    result[0:2, 0:2] = plane(-k, length)
    result[2:4, 2:4] = plane(k, length)
    result[4, 5] = length / GAMMA2
    return result


def product(*maps: np.ndarray) -> np.ndarray:
    result = np.eye(6)
    for matrix in maps:
        result = matrix @ result
    return result


def sector_bend(rho: float, angle: float) -> np.ndarray:
    """Hard-edge sector-bend map in (x, x', y, y', zeta, delta)."""
    c = math.cos(angle)
    s = math.sin(angle)
    length = rho * angle
    result = np.eye(6)
    result[0, 0] = c
    result[0, 1] = rho * s
    result[0, 5] = rho * (1.0 - c)
    result[1, 0] = -s / rho
    result[1, 1] = c
    result[1, 5] = s
    result[2, 3] = length
    result[4, 0] = -s
    result[4, 1] = -rho * (1.0 - c)
    result[4, 5] = rho * (s - angle) + length / GAMMA2
    return result


def dba() -> np.ndarray:
    rho = 2.0
    angle = math.pi / 6.0
    bend = sector_bend(rho, angle)
    return product(bend, drift(1.0), quadrupole(-6.371966681365967, 0.2),
                   drift(1.0), bend)


def write_analytic_map(name: str, matrix: np.ndarray) -> None:
    output_file = ROOT / name / "analytic-map.txt"
    header = "Exact hard-edge map in (x, x', y, y', zeta, delta); rows follow."
    printable = matrix.copy()
    printable[np.abs(printable) < 1.0e-14] = 0.0
    np.savetxt(output_file, printable, fmt="% .12e", header=header)


def parse_map(output: str) -> np.ndarray:
    marker = "Combined linear transfer map"
    start = output.find(marker)
    if start < 0:
        raise RuntimeError("combined transfer map was not printed")
    rows = []
    for line in output[start:].splitlines()[1:]:
        values = re.findall(r"[-+]?\d\.\d+e[-+]\d+", line)
        if len(values) == 6:
            rows.append([float(value) for value in values])
        if len(rows) == 6:
            return np.asarray(rows)
    raise RuntimeError("could not parse six transfer-map rows")


def run(executable: Path, name: str) -> np.ndarray:
    input_file = ROOT / name / f"map-2-{name}.in"
    completed = subprocess.run(
        [str(executable), "--info", "2", input_file.name],
        cwd=input_file.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    output_file = input_file.with_suffix(".out")
    output_file.write_text(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{input_file.name} exited with status {completed.returncode}; see {output_file}"
        )
    return parse_map(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path, help="path to the OPALX executable")
    args = parser.parse_args()

    expected = {
        "drift": drift(1.0),
        "quadrupole": quadrupole(1.0, 0.4),
        "fodo": product(
            drift(0.4), quadrupole(-1.0, 0.2), drift(0.8),
            quadrupole(1.0, 0.2), drift(0.4)),
        "dba": dba(),
    }
    # The exact matrices are continuous-s models. The finite dt advances the reference endpoint
    # by at most about c*dt, so millimetre-scale absolute tolerances are used here.
    tolerances = {
        "drift": 1.0e-3,
        "quadrupole": 2.0e-3,
        "fodo": 5.0e-3,
        "dba": 2.0e-3,
    }

    failed = False
    for name, reference in expected.items():
        write_analytic_map(name, reference)
        measured = run(args.executable.resolve(), name)
        error = float(np.max(np.abs(measured - reference)))
        passed = error <= tolerances[name]
        detail = ""
        if name == "dba":
            achromat_error = max(abs(measured[0, 5]), abs(measured[1, 5]))
            detail = f"  max(|R16|,|R26|)={achromat_error:.6e}"
        print(
            f"{name:10s} max|M_OPALX-M_exact| = {error:.6e}  "
            f"tolerance={tolerances[name]:.1e}{detail}  {'PASS' if passed else 'FAIL'}"
        )
        failed |= not passed
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
