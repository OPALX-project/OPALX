#!/usr/bin/env python3
"""Compare the Wolski chicane entrance plane between BMAD and OPALX.

The no-RF benchmark places the first comparison plane at the chicane entrance.
For BMAD, the stock ``beam_track_example`` program does not write an
intermediate beam snapshot at the ``chicane_in`` marker, so the exact entrance
phase space is the input beam file itself. For OPALX, the entrance temporal
monitor samples the same plane before any beamline transport occurs.

To make the momentum comparison apples-to-apples, the OPALX monitor RMS
momenta, stored in local ``beta*gamma`` coordinates, are converted to the BMAD
canonical basis using

.. math::

   p_x^\mathrm{can} = \frac{\beta\gamma_x}{(\beta\gamma)_0},
   \qquad
   p_y^\mathrm{can} = \frac{\beta\gamma_y}{(\beta\gamma)_0},
   \qquad
   \delta = \frac{\beta\gamma_s}{(\beta\gamma)_0} - 1.

At the entrance plane, where the reference axis is aligned with the beamline
and the injected distribution has zero transverse momentum, the RMS
longitudinal spread obeys

.. math::

   \sigma_\delta = \frac{\sigma_{\beta\gamma_s}}{(\beta\gamma)_0}.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parent
OPALX_DIR = ROOT / "opalx"
BMAD_INPUT = ROOT / "wolski_tesla_chicane_bmad_negchirp.beam"
OPALX_INPUT = ROOT / "wolski_tesla_chicane_distribution_negchirp.txt"
OPALX_MONITOR = OPALX_DIR / "wolski_tesla_chicane_opalx_negchirp_dt_1em11_monitor_entrance_Monitors.stat"
OUT_MD = ROOT / "wolski_first_monitor_summary.md"
OUT_CSV = ROOT / "wolski_first_monitor_summary.csv"

ELECTRON_MASS_MEV = 0.51099895


def pop_rms(values: list[float]) -> float:
    mean = fmean(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / len(values))


def read_bmad_beam(path: Path) -> tuple[float, list[list[float]]]:
    """Read a BMAD beam file and return ``(beta_gamma_ref, rows)``.

    Each data row is ``[x, px_can, y, py_can, z, delta]``.
    """
    p0c_ev = None
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith("# p0c ="):
            p0c_ev = float(line.split("=", 1)[1].strip())
            continue
        if line.startswith("#"):
            continue
        rows.append([float(value) for value in line.split()])

    if p0c_ev is None:
        raise RuntimeError(f"Could not read p0c from {path}")

    beta_gamma_ref = p0c_ev / (ELECTRON_MASS_MEV * 1.0e6)
    return beta_gamma_ref, rows


def read_opalx_distribution(path: Path, beta_gamma_ref: float) -> list[list[float]]:
    """Read an OPALX FROMFILE distribution and convert it to canonical momenta.

    The input rows are ``x px y py z pz`` with momenta in local ``beta*gamma``.
    The returned rows are ``[x, px_can, y, py_can, z, delta]``.
    """
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as stream:
        lines = [line.strip() for line in stream if line.strip()]

    for line in lines[2:]:
        x, bgx, y, bgy, z, bgs = (float(value) for value in line.split())
        rows.append([x, bgx / beta_gamma_ref, y, bgy / beta_gamma_ref, z, bgs / beta_gamma_ref - 1.0])

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


def read_opalx_monitor_row(path: Path, monitor_name: str) -> dict[str, float]:
    """Read one row from an OPALX ``*_Monitors.stat`` SDDS file."""
    column_names: list[str] = []
    data_lines: list[str] = []
    in_data = False
    in_column = False

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
        elif in_data and not line.startswith("&") and line.strip():
            data_lines.append(line.rstrip())

    if len(data_lines) < 4:
        raise RuntimeError(f"Unexpected monitor SDDS payload in {path}")

    row_lines = data_lines[3:]
    for row in row_lines:
        parts = row.split()
        if parts[0] != monitor_name:
            continue
        values = [parts[0], *[float(value) for value in parts[1:]]]
        row_dict: dict[str, float] = {}
        for name, value in zip(column_names, values):
            if name == "name":
                continue
            row_dict[name] = float(value)
        return row_dict

    raise RuntimeError(f"Monitor row '{monitor_name}' not found in {path}")


def format_value(value: float) -> str:
    return f"{value:.5f}"


def write_outputs(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "row",
        "source",
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
        "# Wolski chicane first-monitor comparison",
        "",
        "BMAD and OPALX are compared at the first plane in the no-RF negative-chirp case.",
        "For BMAD, the first plane is the input beam at `chicane_in`.",
        "For OPALX, the first plane is the `wolski_tesla_chicane_entrance_monitor` row.",
        "",
        header,
        rule,
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[name] for name in fieldnames) + " |")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    beta_gamma_ref, bmad_rows = read_bmad_beam(BMAD_INPUT)
    opalx_rows = read_opalx_distribution(OPALX_INPUT, beta_gamma_ref)
    opalx_monitor = read_opalx_monitor_row(
        OPALX_MONITOR, "wolski_tesla_chicane_entrance_monitor"
    )

    bmad_summary = summarize_phase_space(bmad_rows)
    opalx_input_summary = summarize_phase_space(opalx_rows)
    opalx_monitor_summary = {
        "rms_x": opalx_monitor["rms_x"],
        "rms_y": opalx_monitor["rms_y"],
        "rms_z": opalx_monitor["rms_s"],
        "rms_px_can": opalx_monitor["rms_px"] / beta_gamma_ref,
        "rms_py_can": opalx_monitor["rms_py"] / beta_gamma_ref,
        "rms_delta": opalx_monitor["rms_ps"] / beta_gamma_ref,
    }

    rows = [
        {
            "row": "FIRST",
            "source": "BMAD-input",
            **{key: format_value(value) for key, value in bmad_summary.items()},
        },
        {
            "row": "FIRST",
            "source": "OPALX-input",
            **{key: format_value(value) for key, value in opalx_input_summary.items()},
        },
        {
            "row": "FIRST",
            "source": "OPALX-monitor",
            **{key: format_value(value) for key, value in opalx_monitor_summary.items()},
        },
    ]
    write_outputs(rows)
    print(f"Wrote {OUT_MD}")
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
