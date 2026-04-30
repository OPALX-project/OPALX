#!/Users/adelmann/git/opalx/.venv-h6/bin/python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


WOLSKI_TARGET_SIGMA_Z_M = 0.3e-3
ELECTRON_MASS_EV = 0.510998950e6
BMAD_P0C_EV = 5.000510972840675e9
REFERENCE_BETA_GAMMA = BMAD_P0C_EV / ELECTRON_MASS_EV


@dataclass(frozen=True)
class BeamMoments:
    sigma_x_m: float
    sigma_y_m: float
    sigma_z_m: float
    sigma_px: float
    sigma_py: float
    sigma_delta: float


def _sample_std(values: pd.Series) -> float:
    return float(values.astype(float).std(ddof=1))


def load_bmad_beam(path: Path) -> BeamMoments:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            if len(tokens) < 6:
                continue
            try:
                float(tokens[0])
            except ValueError:
                continue
            if len(tokens) >= 16:
                rows.append([float(tok) for tok in tokens[1:7]])
            else:
                rows.append([float(tok) for tok in tokens[0:6]])

    df = pd.DataFrame(rows, columns=["x", "px", "y", "py", "z", "pz"])
    return BeamMoments(
        sigma_x_m=_sample_std(df["x"]),
        sigma_y_m=_sample_std(df["y"]),
        sigma_z_m=_sample_std(df["z"]),
        sigma_px=_sample_std(df["px"]),
        sigma_py=_sample_std(df["py"]),
        sigma_delta=_sample_std(df["pz"]),
    )


def load_opalx_distribution(path: Path) -> BeamMoments:
    df = pd.read_csv(path, sep=r"\s+", skiprows=2, names=["x", "px", "y", "py", "z", "pz"])
    delta = df["pz"] / REFERENCE_BETA_GAMMA - 1.0
    return BeamMoments(
        sigma_x_m=_sample_std(df["x"]),
        sigma_y_m=_sample_std(df["y"]),
        sigma_z_m=_sample_std(df["z"]),
        sigma_px=_sample_std(df["px"] / REFERENCE_BETA_GAMMA),
        sigma_py=_sample_std(df["py"] / REFERENCE_BETA_GAMMA),
        sigma_delta=_sample_std(delta),
    )


def load_opalx_stat(path: Path) -> BeamMoments:
    data_started = False
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("&data"):
                data_started = True
                continue
            if not data_started or not stripped or stripped.startswith("!"):
                continue
            try:
                rows.append([float(x) for x in stripped.split()])
            except ValueError:
                continue

    df = pd.DataFrame(rows)
    final = df.iloc[-1]
    return BeamMoments(
        sigma_x_m=float(final[5]),
        sigma_y_m=float(final[6]),
        sigma_z_m=float(final[7]),
        sigma_px=float(final[8]) / REFERENCE_BETA_GAMMA,
        sigma_py=float(final[9]) / REFERENCE_BETA_GAMMA,
        sigma_delta=float(final[10]) / REFERENCE_BETA_GAMMA,
    )


def fmt_mm(value: float | None) -> str:
    if value is None:
        return f"{'-':>10s}"
    return f"{value * 1e3:10.6f}"


def fmt_dimless(value: float | None) -> str:
    if value is None:
        return f"{'-':>10s}"
    return f"{value:10.6e}"


def print_row(label: str, moments: BeamMoments | None) -> None:
    if moments is None:
        print(
            f"{label:34s} {fmt_mm(None)} {fmt_mm(None)} {fmt_mm(WOLSKI_TARGET_SIGMA_Z_M)}"
            f" {fmt_dimless(None)} {fmt_dimless(None)} {fmt_dimless(None)}"
        )
        return

    print(
        f"{label:34s}"
        f" {fmt_mm(moments.sigma_x_m)}"
        f" {fmt_mm(moments.sigma_y_m)}"
        f" {fmt_mm(moments.sigma_z_m)}"
        f" {fmt_dimless(moments.sigma_px)}"
        f" {fmt_dimless(moments.sigma_py)}"
        f" {fmt_dimless(moments.sigma_delta)}"
    )


def main() -> None:
    root = Path(__file__).resolve().parent
    opalx_dir = root / "opalx"

    entrance_bmad_neg = load_bmad_beam(root / "wolski_tesla_chicane_bmad_negchirp.beam")
    entrance_opalx_neg = load_opalx_distribution(root / "wolski_tesla_chicane_distribution_negchirp.txt")

    exit_bmad_neg = load_bmad_beam(root / "wolski_no_rf_neg.dat")
    exit_bmad_e0_neg = load_bmad_beam(root / "wolski_no_rf_e0_neg.dat")
    exit_bmad_pos = load_bmad_beam(root / "wolski_no_rf_pos.dat")

    exit_opalx_neg = load_opalx_stat(opalx_dir / "wolski_tesla_chicane_opalx_negchirp.stat")
    exit_opalx_neg_dt = load_opalx_stat(opalx_dir / "wolski_tesla_chicane_opalx_negchirp_dt_1em11_cmp.stat")
    exit_opalx_zero_face = load_opalx_stat(opalx_dir / "wolski_tesla_chicane_opalx_signed_e0.stat")

    print("RMS beam-size and momentum comparison at entrance")
    print(
        f"{'Case':34s} {'sigma_x':>10s} {'sigma_y':>10s} {'sigma_z/s':>10s}"
        f" {'sigma_px':>10s} {'sigma_py':>10s} {'sigma_delta':>10s}"
    )
    print_row("Bmad no-RF input (-h)", entrance_bmad_neg)
    print_row("OPALX no-RF input (-h)", entrance_opalx_neg)
    print()

    print("RMS beam-size and momentum comparison at exit")
    print(
        f"{'Case':34s} {'sigma_x':>10s} {'sigma_y':>10s} {'sigma_z/s':>10s}"
        f" {'sigma_px':>10s} {'sigma_py':>10s} {'sigma_delta':>10s}"
    )
    print_row("Wolski target", None)
    print_row("Bmad no-RF (-h)", exit_bmad_neg)
    print_row("Bmad no-RF (-h) zero-face", exit_bmad_e0_neg)
    print_row("Bmad no-RF (+h)", exit_bmad_pos)
    print_row("OPALX no-RF (-h) [fixed-t]", exit_opalx_neg)
    print_row("OPALX no-RF (-h), dt=1e-11 [fixed-t]", exit_opalx_neg_dt)
    print_row("OPALX signed zero-face [fixed-t]", exit_opalx_zero_face)


if __name__ == "__main__":
    main()
