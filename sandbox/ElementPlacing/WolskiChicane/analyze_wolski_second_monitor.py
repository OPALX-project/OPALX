#!/usr/bin/env python3
"""Compare the Wolski chicane exit plane between BMAD and OPALX.

For the no-RF negative-chirp benchmark, the BMAD final beam file is the bunch
at the downstream end of the lattice, which coincides with the ``chicane_out``
marker. On the OPALX side, the exit temporal monitor can record multiple rows.
The physically relevant row is chosen as the end-monitor sample with

.. math::

   s \approx 10\ \mathrm{m},

that is, the row with ``name == wolski_tesla_chicane_endmonitor`` and the
smallest ``|s - 10|`` among rows with ``s > 1``.

The OPALX monitor momenta are converted from local ``beta*gamma`` to the BMAD
canonical basis using

.. math::

   p_x^\mathrm{can} = \frac{\beta\gamma_x}{(\beta\gamma)_0},
   \qquad
   p_y^\mathrm{can} = \frac{\beta\gamma_y}{(\beta\gamma)_0},
   \qquad
   \delta \approx \frac{\beta\gamma_s}{(\beta\gamma)_0} - 1.

The longitudinal spread column is reported as ``rms_z`` for continuity with the
first-plane table, but note that BMAD writes phase-space ``z`` while OPALX
reports local ``s`` at the monitor plane.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parent
OPALX_DIR = ROOT / "opalx"
BMAD_EXIT = ROOT / "wolski_no_rf_neg.dat"
OPALX_MONITOR = OPALX_DIR / "wolski_tesla_chicane_opalx_negchirp_dt_1em11_monitor_both_Monitors.stat"
BMAD_INPUT = ROOT / "wolski_tesla_chicane_bmad_negchirp.beam"
OUT_MD = ROOT / "wolski_second_monitor_summary.md"
OUT_CSV = ROOT / "wolski_second_monitor_summary.csv"

ELECTRON_MASS_MEV = 0.51099895
END_MONITOR_NAME = "wolski_tesla_chicane_endmonitor"
EXPECTED_EXIT_S_M = 10.0


def pop_rms(values: list[float]) -> float:
    mean = fmean(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / len(values))


def read_bmad_reference_beta_gamma(path: Path) -> float:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# p0c ="):
            p0c_ev = float(line.split("=", 1)[1].strip())
            return p0c_ev / (ELECTRON_MASS_MEV * 1.0e6)
    raise RuntimeError(f"Could not read p0c from {path}")


def read_bmad_exit_rows(path: Path) -> list[list[float]]:
    """Read the BMAD final beam file rows as ``[x, px, y, py, z, delta]``."""
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        parts = line.split()
        rows.append(
            [
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
                float(parts[4]),
                float(parts[5]),
                float(parts[6]),
            ]
        )
    return rows


def summarize_phase_space(rows: list[list[float]]) -> dict[str, float]:
    columns = list(zip(*rows))
    return {
        "rms_x": pop_rms(list(columns[0])),
        "rms_px_can": pop_rms(list(columns[1])),
        "rms_y": pop_rms(list(columns[2])),
        "rms_py_can": pop_rms(list(columns[3])),
        "rms_z": pop_rms(list(columns[4])),
        "rms_delta": pop_rms(list(columns[5])),
    }


def read_opalx_monitor_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    column_names: list[str] = []
    rows: list[list[str]] = []
    in_column = False
    in_data = False

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("&column"):
            in_column = True
            continue
        if in_column and line.startswith("        name="):
            column_names.append(line.split("=", 1)[1].rstrip(","))
            continue
        if in_column and line.startswith("&end"):
            in_column = False
            continue
        if line.startswith("&data"):
            in_data = True
            continue
        if in_data and not line.startswith("&") and line.strip():
            rows.append(line.split())

    return column_names, rows[3:]


def select_exit_monitor_row(path: Path) -> dict[str, float]:
    column_names, rows = read_opalx_monitor_rows(path)
    candidates: list[dict[str, float]] = []

    for parts in rows:
        if parts[0] != END_MONITOR_NAME:
            continue
        row = {"name": parts[0]}
        for name, value in zip(column_names[1:], parts[1:]):
            row[name] = float(value)
        if row["s"] > 1.0:
            candidates.append(row)

    if not candidates:
        raise RuntimeError(f"No usable end monitor rows found in {path}")

    return min(candidates, key=lambda row: abs(row["s"] - EXPECTED_EXIT_S_M))


def format_value(value: float) -> str:
    return f"{value:.5f}"


def write_outputs(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "row",
        "source",
        "numParticles",
        "rms_x",
        "rms_y",
        "rms_z",
        "rms_px_can",
        "rms_py_can",
        "rms_delta",
    ]

    with OUT_CSV.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    header = "| " + " | ".join(fieldnames) + " |"
    rule = "| " + " | ".join("---" for _ in fieldnames) + " |"
    lines = [
        "# Wolski chicane second-monitor comparison",
        "",
        "BMAD and OPALX are compared at the downstream chicane plane in the no-RF negative-chirp case.",
        "For BMAD, the final beam file is the `chicane_out` plane.",
        "For OPALX, the selected row is the end-monitor sample with `s` closest to `10 m`.",
        "",
        header,
        rule,
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[name] for name in fieldnames) + " |")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    beta_gamma_ref = read_bmad_reference_beta_gamma(BMAD_INPUT)
    bmad_exit_rows = read_bmad_exit_rows(BMAD_EXIT)
    opalx_monitor = select_exit_monitor_row(OPALX_MONITOR)

    bmad_summary = summarize_phase_space(bmad_exit_rows)
    opalx_summary = {
        "rms_x": opalx_monitor["rms_x"],
        "rms_y": opalx_monitor["rms_y"],
        "rms_z": opalx_monitor["rms_s"],
        "rms_px_can": opalx_monitor["rms_px"] / beta_gamma_ref,
        "rms_py_can": opalx_monitor["rms_py"] / beta_gamma_ref,
        "rms_delta": opalx_monitor["rms_ps"] / beta_gamma_ref,
    }

    rows = [
        {
            "row": "SECOND",
            "source": "BMAD-exit",
            "numParticles": str(len(bmad_exit_rows)),
            **{key: format_value(value) for key, value in bmad_summary.items()},
        },
        {
            "row": "SECOND",
            "source": "OPALX-monitor",
            "numParticles": str(int(opalx_monitor["numParticles"])),
            **{key: format_value(value) for key, value in opalx_summary.items()},
        },
    ]
    write_outputs(rows)
    print(f"Wrote {OUT_MD}")
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
