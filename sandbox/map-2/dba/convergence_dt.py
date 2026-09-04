#!/usr/bin/env python3
"""Convergence study of the full DBA transfer map with respect to TRACK::DT."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "opalx-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "opalx-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MAP2 = Path(__file__).resolve().parents[1]
DBA = Path(__file__).resolve().parent

import sys

sys.path.insert(0, str(MAP2))
from check_maps import dba, parse_map  # noqa: E402


DT_VALUES = (1.0e-10, 1.0e-11, 1.0e-12, 1.0e-13)


def replace_dt(input_text: str, dt: float) -> str:
    updated, count = re.subn(
        r"(\bDT\s*=\s*)[-+]?\d+(?:\.\d*)?[eE][-+]?\d+",
        rf"\g<1>{dt:.1e}",
        input_text,
        count=1,
    )
    if count != 1:
        raise RuntimeError("expected exactly one TRACK::DT assignment")
    return updated


def parse_reported_residual(output: str, label: str) -> float:
    match = re.search(rf"{re.escape(label)}\s*=\s*([-+]?\d\.\d+e[-+]\d+)", output)
    if not match:
        raise RuntimeError(f"could not parse reported residual: {label}")
    return float(match.group(1))


def run_case(
    executable: Path, source: str, dt: float, temporary_root: Path
) -> dict[str, float | str]:
    label = f"{dt:.0e}"
    run_directory = temporary_root / label / "dba"
    run_directory.mkdir(parents=True)
    shutil.copy2(MAP2 / "reference-particle.txt", run_directory.parent / "reference-particle.txt")
    input_file = run_directory / "map-2-dba.in"
    input_file.write_text(replace_dt(source, dt))

    completed = subprocess.run(
        [str(executable), "--info", "2", input_file.name],
        cwd=run_directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    output_file = DBA / f"map-2-dba-dt-{label}.out"
    output_file.write_text(completed.stdout)
    if completed.returncode != 0:
        reason = "OPALX rejected this time step"
        if "time step is too long" in completed.stdout:
            reason = "time step too long for QACH field extent"
        return {
            "dt_s": dt,
            "status": "INVALID",
            "reason": reason,
            "full_matrix_error": np.nan,
            "abs_R16": np.nan,
            "abs_R26": np.nan,
            "dispersion_error": np.nan,
            "determinant_error": np.nan,
            "canonical_J_error": np.nan,
        }

    measured = parse_map(completed.stdout)
    exact = dba()
    return {
        "dt_s": dt,
        "status": "OK",
        "reason": "",
        "full_matrix_error": float(np.max(np.abs(measured - exact))),
        "abs_R16": float(abs(measured[0, 5])),
        "abs_R26": float(abs(measured[1, 5])),
        "dispersion_error": float(max(abs(measured[0, 5]), abs(measured[1, 5]))),
        "determinant_error": parse_reported_residual(completed.stdout, "|det(M) - 1|"),
        "canonical_J_error": parse_reported_residual(
            completed.stdout, "max_ij |(M^T J M - J)_ij|"
        ),
    }


def add_observed_orders(frame: pd.DataFrame) -> pd.DataFrame:
    for column in ("full_matrix_error", "dispersion_error"):
        order = np.full(len(frame), np.nan)
        for index in range(1, len(frame)):
            order[index] = np.log(frame.loc[index - 1, column] / frame.loc[index, column]) / np.log(
                frame.loc[index - 1, "dt_s"] / frame.loc[index, "dt_s"]
            )
        frame[f"{column}_order"] = order
    return frame


def make_plot(frame: pd.DataFrame, output: Path) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    series = (
        ("full_matrix_error", r"$\max_{ij}|M_{ij}-M_{ij}^{\rm exact}|$", "o"),
        ("dispersion_error", r"$\max(|R_{16}|,|R_{26}|)$", "s"),
        ("determinant_error", r"$|\det M-1|$", "^"),
        ("canonical_J_error", r"$\max_{ij}|(M^TJM-J)_{ij}|$", "D"),
    )
    for column, label, marker in series:
        axis.loglog(frame["dt_s"], frame[column], marker=marker, linewidth=1.6, label=label)
    axis.set_xlabel(r"Time step $\Delta t$ [s]")
    axis.set_ylabel("Absolute residual")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(frameon=False)
    axis.invert_xaxis()
    figure.savefig(output, dpi=200)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path, help="path to the OPALX executable")
    args = parser.parse_args()
    executable = args.executable.resolve()
    source = (DBA / "map-2-dba.in").read_text()

    with tempfile.TemporaryDirectory(prefix="opalx-map-2-dba-") as directory:
        temporary_root = Path(directory)
        rows = [run_case(executable, source, dt, temporary_root) for dt in DT_VALUES]

    frame = add_observed_orders(pd.DataFrame(rows))
    frame.to_csv(DBA / "convergence-dt.csv", index=False, float_format="%.12e")
    make_plot(frame, DBA / "convergence-dt.png")
    print(frame.to_string(index=False, float_format=lambda value: f"{value:.6e}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
