#!/usr/bin/env python3
"""Validate the CAIN-to-OPALX conversion and three-container deck wiring."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from convert_cain_pairs import (
    ELECTRON_REST_ENERGY_EV,
    OUTPUT_COLUMNS,
    SPECIES,
    SPEED_OF_LIGHT_M_PER_S,
    read_cain,
)


def read_emitted(path: Path) -> pd.DataFrame:
    with path.open(encoding="utf-8") as stream:
        count_line = stream.readline().strip()
        header = stream.readline().split()
    try:
        declared_count = int(count_line)
    except ValueError as exc:
        raise ValueError(f"{path}: first line is not an integer particle count") from exc
    if header != OUTPUT_COLUMNS:
        raise ValueError(f"{path}: header {header!r} != {OUTPUT_COLUMNS!r}")

    data = pd.read_csv(path, sep=r"\s+", skiprows=2, names=header, engine="python")
    if len(data) != declared_count:
        raise ValueError(f"{path}: declared {declared_count} rows but read {len(data)}")
    if not np.isfinite(data.to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: non-finite output value")
    return data


def max_abs_difference(lhs: pd.DataFrame, rhs: pd.DataFrame, columns: list[str]) -> dict[str, float]:
    return {
        column: float(np.max(np.abs(lhs[column].to_numpy() - rhs[column].to_numpy())))
        for column in columns
    }


def validate_species(source: pd.DataFrame, species_id: int, output: pd.DataFrame) -> dict[str, object]:
    selected = source.loc[source["species"] == species_id].reset_index(drop=True)
    if len(selected) != len(output):
        raise ValueError(
            f"species {species_id}: CAIN has {len(selected)} rows, OPALX file has {len(output)}"
        )

    expected = selected.loc[:, ["x", "y", "z", "px", "py", "pz"]].copy()
    expected.loc[:, ["px", "py", "pz"]] /= ELECTRON_REST_ENERGY_EV
    expected["birth_time"] = selected["ct"] / SPEED_OF_LIGHT_M_PER_S
    expected = expected.loc[:, OUTPUT_COLUMNS]

    np.testing.assert_allclose(
        output.to_numpy(), expected.to_numpy(), rtol=5.0e-15, atol=5.0e-24
    )
    return {
        "count": int(len(output)),
        "max_abs_conversion_error": max_abs_difference(output, expected, OUTPUT_COLUMNS),
        "birth_time_min_s": float(output["birth_time"].min()),
        "birth_time_max_s": float(output["birth_time"].max()),
        "z_min_m": float(output["z"].min()),
        "z_max_m": float(output["z"].max()),
    }


def normalized_deck(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    return re.sub(r"\s+", "", text).upper()


def validate_deck(path: Path) -> dict[str, object]:
    deck = normalized_deck(path)
    required = {
        "container_order": "BEAMS={PRIMARYELECTRONS,GAMMAGAMMAELECTRONS,GAMMAGAMMAPOSITRONS}",
        "electron_species": "GAMMAGAMMAELECTRONS:BEAM,PARTICLE=ELECTRON",
        "positron_species": "GAMMAGAMMAPOSITRONS:BEAM,PARTICLE=POSITRON",
        "witness_indices": 'WITNESS_CONTAINERS="1,2"',
        "electron_file": 'FNAME="INPUT/CAIN_ELECTRONS.EMITTEDFROMFILE"',
        "positron_file": 'FNAME="INPUT/CAIN_POSITRONS.EMITTEDFROMFILE"',
    }
    checks = {name: fragment in deck for name, fragment in required.items()}
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError(f"{path}: failed deck wiring checks: {failed}")
    return {
        "checks": checks,
        "container_map": {
            "0": "primary electron source and field target",
            "1": "CAIN electron witness",
            "2": "CAIN positron witness",
        },
    }


def validate(input_path: Path, output_dir: Path, deck_path: Path) -> dict[str, object]:
    source = read_cain(input_path)
    electron = source.loc[source["species"] == 2].reset_index(drop=True)
    positron = source.loc[source["species"] == 3].reset_index(drop=True)
    pair_columns = ["ct", "x", "y", "z"]
    paired_spacetime = all(
        np.array_equal(electron[column].to_numpy(), positron[column].to_numpy())
        for column in pair_columns
    )
    if not paired_spacetime:
        raise ValueError("CAIN electron/positron rows do not have identical paired spacetime data")

    mass_shell = source["energy"].to_numpy() ** 2
    mass_shell -= np.square(source[["px", "py", "pz"]].to_numpy()).sum(axis=1)
    if np.any(mass_shell <= 0.0):
        raise ValueError("CAIN data contain a non-timelike four-momentum")
    inferred_mass = np.sqrt(mass_shell)

    report: dict[str, object] = {
        "source_rows": int(len(source)),
        "paired_spacetime_exact": paired_spacetime,
        "unique_birth_times": int(source["ct"].nunique()),
        "inferred_rest_energy_eV": {
            "min": float(inferred_mass.min()),
            "max": float(inferred_mass.max()),
        },
        "species": {},
        "deck": validate_deck(deck_path),
        "field_participation_contract": {
            "source_container": 0,
            "witness_containers": [1, 2],
            "witness_role": "gather-only; validated dynamically by the BeamBeam regression",
        },
    }
    species_report: dict[str, object] = {}
    for species_id, (name, filename) in SPECIES.items():
        species_report[name] = validate_species(
            source, species_id, read_emitted(output_dir / filename)
        )
    report["species"] = species_report
    return report


def parse_args() -> argparse.Namespace:
    workflow = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=workflow.parent / "track-e-p" / "fort98.txt"
    )
    parser.add_argument("--output-dir", type=Path, default=workflow / "opalx" / "input")
    parser.add_argument(
        "--deck", type=Path, default=workflow / "opalx" / "gamma_gamma_pairs_timed.in"
    )
    parser.add_argument(
        "--report", type=Path, default=workflow / "results" / "input_validation.json"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = validate(args.input, args.output_dir, args.deck)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report: {args.report}")


if __name__ == "__main__":
    main()
