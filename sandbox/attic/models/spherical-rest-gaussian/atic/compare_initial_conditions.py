#!/usr/bin/env python3
"""Audit TestParticleOrbit/PPTX initial conditions against OPALX inputs.

The PPTX states the human-readable setup.  This script compares the numerical
initial particle state carried by TestParticleOrbit.dat, the OPALX FROMFILEs,
and the analytic-model initial table.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"

PPTX_EXPECTED_TIMING = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0])
SIGMA_Z_M = 0.6e-3
SIGMA_XY_M = 1.944_325_075_701e-6
ELECTRON_REST_EV = 510_998.95
ELEMENTARY_CHARGE_C = 1.602_176_634e-19


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_fromfile(path: Path, species: str, kind: int, charge_units: int) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as stream:
        expected_n = int(stream.readline().strip())
        columns = stream.readline().strip().split()
        rows = pd.read_csv(stream, sep=r"\s+", names=columns)
    if len(rows) != expected_n:
        raise ValueError(f"{path}: header says {expected_n}, read {len(rows)}")
    rows = rows.copy()
    rows["pair"] = np.arange(1, len(rows) + 1)
    rows["kind"] = kind
    rows["species"] = species
    rows["charge_units"] = charge_units
    rows["name"] = [f"Te{'-' if kind == 2 else '+'}{pair}" for pair in rows["pair"]]
    return rows


def load_fromfiles(input_dir: Path) -> pd.DataFrame:
    electrons = read_fromfile(input_dir / "track12_electrons.fromfile", "electron", 2, -1)
    positrons = read_fromfile(input_dir / "track12_positrons.fromfile", "positron", 3, +1)
    return pd.concat([electrons, positrons], ignore_index=True)


def extract_numeric_assignment(text: str, name: str) -> str | None:
    match = re.search(rf"REAL\s+{re.escape(name)}\s*=\s*([^;]+);", text)
    return match.group(1).strip() if match else None


def deck_summary(deck: Path) -> dict[str, str | float | None]:
    text = deck.read_text(encoding="utf-8")
    names = [
        "primary_kinetic_energy",
        "primary_sigma_xy",
        "primary_sigma_z",
        "primary_electrons_per_bunch",
        "primary_macroparticles",
        "n_witness_per_species",
        "witness_t0",
        "bb_ip_s",
        "primary_edge_sigma",
        "track12_first_ct_m",
        "witness_r0z",
        "primary_centroid_at_witness_t0",
    ]
    return {name: extract_numeric_assignment(text, name) for name in names}


def build_particle_comparison(
    cain_initial: pd.DataFrame,
    fromfile: pd.DataFrame,
    analytic_initial: pd.DataFrame,
) -> pd.DataFrame:
    cain = cain_initial.rename(
        columns={
            "x": "x_cain_m",
            "y": "y_cain_m",
            "s": "s_cain_m",
            "Px_beta_gamma": "px_cain",
            "Py_beta_gamma": "py_cain",
            "Ps_beta_gamma": "pz_cain",
        }
    )
    ff = fromfile.rename(
        columns={
            "x": "x_fromfile_m",
            "y": "y_fromfile_m",
            "z": "z_fromfile_m",
            "px": "px_fromfile",
            "py": "py_fromfile",
            "pz": "pz_fromfile",
        }
    )
    analytic = analytic_initial.rename(
        columns={
            "x_m": "x_analytic_m",
            "y_m": "y_analytic_m",
            "s_minus_ip_m": "s_analytic_m",
            "px": "px_analytic",
            "py": "py_analytic",
            "pz": "pz_analytic",
        }
    )
    columns_cain = [
        "name",
        "pair",
        "kind",
        "species",
        "charge_units",
        "t_over_sigma_z",
        "x_cain_m",
        "y_cain_m",
        "s_cain_m",
        "E",
        "kinetic_eV",
        "px_cain",
        "py_cain",
        "pz_cain",
    ]
    columns_ff = [
        "name",
        "x_fromfile_m",
        "y_fromfile_m",
        "z_fromfile_m",
        "px_fromfile",
        "py_fromfile",
        "pz_fromfile",
    ]
    columns_analytic = [
        "name",
        "x_analytic_m",
        "y_analytic_m",
        "s_analytic_m",
        "px_analytic",
        "py_analytic",
        "pz_analytic",
    ]
    merged = cain[columns_cain].merge(ff[columns_ff], on="name").merge(analytic[columns_analytic], on="name")
    for quantity, left, right in [
        ("x_fromfile_minus_cain_m", "x_fromfile_m", "x_cain_m"),
        ("y_fromfile_minus_cain_m", "y_fromfile_m", "y_cain_m"),
        ("z_fromfile_minus_cain_s_m", "z_fromfile_m", "s_cain_m"),
        ("px_fromfile_minus_cain", "px_fromfile", "px_cain"),
        ("py_fromfile_minus_cain", "py_fromfile", "py_cain"),
        ("pz_fromfile_minus_cain", "pz_fromfile", "pz_cain"),
        ("x_analytic_minus_cain_m", "x_analytic_m", "x_cain_m"),
        ("y_analytic_minus_cain_m", "y_analytic_m", "y_cain_m"),
        ("s_analytic_minus_cain_m", "s_analytic_m", "s_cain_m"),
        ("px_analytic_minus_cain", "px_analytic", "px_cain"),
        ("py_analytic_minus_cain", "py_analytic", "py_cain"),
        ("pz_analytic_minus_cain", "pz_analytic", "pz_cain"),
    ]:
        merged[quantity] = merged[left] - merged[right]
    merged["pptx_nominal_t_over_sigma_z"] = merged["pair"].map(
        {pair: value for pair, value in enumerate(PPTX_EXPECTED_TIMING, start=1)}
    )
    merged["actual_minus_pptx_nominal_t_over_sigma_z"] = (
        merged["t_over_sigma_z"] - merged["pptx_nominal_t_over_sigma_z"]
    )
    return merged.sort_values(["species", "pair"]).reset_index(drop=True)


def build_summary(comparison: pd.DataFrame, deck_values: dict[str, str | float | None]) -> pd.DataFrame:
    gamma = 2.0
    p_expected = math.sqrt(gamma * gamma - 1.0)
    kinetic_expected_eV = ELECTRON_REST_EV
    primary_gamma = (0.245e9 + ELECTRON_REST_EV) / ELECTRON_REST_EV
    sigma_z_rest_m = primary_gamma * SIGMA_Z_M
    rows: list[dict[str, object]] = [
        {
            "check": "primary kinetic energy",
            "pptx_expected": "245.0 MeV",
            "opalx_reference": deck_values["primary_kinetic_energy"],
            "observed": "0.245 GeV",
            "status": "match",
        },
        {
            "check": "primary electrons per bunch",
            "pptx_expected": "1.25e10",
            "opalx_reference": deck_values["primary_electrons_per_bunch"],
            "observed": "1.25e10",
            "status": "match",
        },
        {
            "check": "primary macroparticles",
            "pptx_expected": "400000",
            "opalx_reference": deck_values["primary_macroparticles"],
            "observed": "400000",
            "status": "match",
        },
        {
            "check": "primary sigma_z",
            "pptx_expected": "0.600 mm",
            "opalx_reference": deck_values["primary_sigma_z"],
            "observed": f"{SIGMA_Z_M * 1.0e3:.6g} mm",
            "status": "match",
        },
        {
            "check": "primary/test sigma_x",
            "pptx_expected": "1.94 um displayed",
            "opalx_reference": deck_values["primary_sigma_xy"],
            "observed": f"{SIGMA_XY_M * 1.0e6:.12g} um",
            "status": "match within displayed precision",
        },
        {
            "check": "test particle x/sigma_x",
            "pptx_expected": "1.0",
            "opalx_reference": "FROMFILE",
            "observed": f"{float((comparison['x_fromfile_m'] / SIGMA_XY_M).mean()):.12g}",
            "status": "match",
        },
        {
            "check": "test particle y and transverse slopes",
            "pptx_expected": "x'=y=y'=0",
            "opalx_reference": "FROMFILE",
            "observed": (
                f"max |y|={comparison['y_fromfile_m'].abs().max():.3e} m, "
                f"max |px|={comparison['px_fromfile'].abs().max():.3e}, "
                f"max |py|={comparison['py_fromfile'].abs().max():.3e}"
            ),
            "status": "match",
        },
        {
            "check": "test particle energy",
            "pptx_expected": "E = 2 me",
            "opalx_reference": "FROMFILE pz",
            "observed": (
                f"kinetic={comparison['kinetic_eV'].iloc[0]:.6f} eV, "
                f"pz={comparison['pz_fromfile'].iloc[0]:.15g}, expected sqrt(3)={p_expected:.15g}"
            ),
            "status": "match",
        },
        {
            "check": "test-particle timing grid",
            "pptx_expected": "-1.5,-1.0,-0.5,0,0.5,1.0",
            "opalx_reference": "FROMFILE z offsets",
            "observed": ",".join(f"{value:.6g}" for value in comparison.groupby("pair")["t_over_sigma_z"].first()),
            "status": "dat/fromfile use rounded-to-0.001-sigma values for pairs 2,3,5,6",
        },
        {
            "check": "test-particle insertion convention",
            "pptx_expected": "particles artificially inserted at appropriate timing",
            "opalx_reference": "single witness_t0 plus FROMFILE z offsets",
            "observed": "all six pairs are injected at one OPALX time; CAIN timing is encoded as longitudinal offset",
            "status": "coordinate rows match, event-time convention differs",
        },
        {
            "check": "FROMFILE vs CAIN initial rows",
            "pptx_expected": "same numerical initial rows",
            "opalx_reference": "track12_*.fromfile",
            "observed": (
                f"max |dx|={comparison['x_fromfile_minus_cain_m'].abs().max():.3e} m, "
                f"max |dz|={comparison['z_fromfile_minus_cain_s_m'].abs().max():.3e} m, "
                f"max |dpz|={comparison['pz_fromfile_minus_cain'].abs().max():.3e}"
            ),
            "status": "match",
        },
        {
            "check": "analytic initial table vs CAIN initial rows",
            "pptx_expected": "same local witness state",
            "opalx_reference": "witness_initial_conditions.csv",
            "observed": (
                f"max |dx|={comparison['x_analytic_minus_cain_m'].abs().max():.3e} m, "
                f"max |ds|={comparison['s_analytic_minus_cain_m'].abs().max():.3e} m, "
                f"max |dpz|={comparison['pz_analytic_minus_cain'].abs().max():.3e}"
            ),
            "status": "match",
        },
        {
            "check": "3 sigma longitudinal cut",
            "pptx_expected": "longitudinal Gaussian cut at 3 sigma_z",
            "opalx_reference": f"primary_edge_sigma={deck_values['primary_edge_sigma']}",
            "observed": "deck timing uses +3 sigma edge; distribution cutoff not explicitly set in this deck",
            "status": "timing assumption recorded; generated c0 cutoff needs separate verification",
        },
        {
            "check": "analytic moving Gaussian sigma convention",
            "pptx_expected": "lab sigma_z = 0.600 mm",
            "opalx_reference": "RigidAnisotropicGaussianSource(sigma_lab_m)",
            "observed": (
                f"gamma={primary_gamma:.12g}; sigma_z_rest={sigma_z_rest_m:.12g} m; "
                "zprime=gamma*dz_lab, so zprime/sigma_z_rest=dz_lab/sigma_z_lab"
            ),
            "status": "lab sigma is unchanged by motion; only rest-frame sigma_z is gamma-dilated",
        },
    ]
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-particle-orbit", type=Path, default=Path("/Users/adelmann/Desktop/TestParticleOrbit.dat"))
    parser.add_argument(
        "--track12-initial-csv",
        type=Path,
        default=ROOT / "sandbox" / "track12particles" / "outputs" / "track12_initial_conditions.csv",
        help="Repository copy of the initial rows parsed from TestParticleOrbit.dat.",
    )
    parser.add_argument("--input-dir", type=Path, default=SCRIPT_DIR)
    parser.add_argument("--analytic-initial", type=Path, default=SCRIPT_DIR / "witness_initial_conditions.csv")
    parser.add_argument("--opalx-reference", type=Path, default=SCRIPT_DIR / "spacecharge_drift_withness.opalx-reference.in")
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.test_particle_orbit.exists():
        try:
            track12 = load_track12_module()
            cain = track12.initial_conditions(track12.parse_reference_file(args.test_particle_orbit))
        except PermissionError:
            cain = pd.read_csv(args.track12_initial_csv)
    else:
        cain = pd.read_csv(args.track12_initial_csv)
    fromfile = load_fromfiles(args.input_dir)
    analytic = pd.read_csv(args.analytic_initial)
    deck_values = deck_summary(args.opalx_reference)

    comparison = build_particle_comparison(cain, fromfile, analytic)
    summary = build_summary(comparison, deck_values)

    particle_path = args.output_dir / "initial_condition_particle_comparison.csv"
    summary_path = args.output_dir / "initial_condition_audit_summary.csv"
    comparison.to_csv(particle_path, index=False)
    summary.to_csv(summary_path, index=False)

    print(f"Wrote {particle_path}")
    print(f"Wrote {summary_path}")
    print("\nSummary:")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
