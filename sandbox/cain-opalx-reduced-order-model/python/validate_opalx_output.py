#!/usr/bin/env python3
"""Validate timed witness emission and species assignment in OPALX H5 output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from validate_pair_input import read_emitted


PRIMARY_KINETIC_ENERGY_GEV = 0.245
ELECTRON_MASS_GEV = 0.000_510_998_95
SPEED_OF_LIGHT_M_PER_S = 299_792_458.0
IP_S_M = 0.035


def numeric_steps(h5: h5py.File) -> list[h5py.Group]:
    return sorted(h5.values(), key=lambda group: int(group.name.rsplit("#", 1)[1]))


def scalar_attribute(group: h5py.Group, name: str) -> float:
    return float(np.asarray(group.attrs[name]).reshape(-1)[0])


def expected_t0() -> float:
    gamma = (PRIMARY_KINETIC_ENERGY_GEV + ELECTRON_MASS_GEV) / ELECTRON_MASS_GEV
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return IP_S_M / (beta * SPEED_OF_LIGHT_M_PER_S)


def validate_witness(
    path: Path, emitted_path: Path, expected_species: int, expected_charge_sign: int
) -> dict[str, object]:
    emitted = read_emitted(emitted_path)
    births = np.sort(emitted["birth_time"].to_numpy() + expected_t0())

    with h5py.File(path, "r") as h5:
        steps = numeric_steps(h5)
        if not steps:
            raise ValueError(f"{path}: no H5 steps")
        observed_counts: list[int] = []
        expected_counts: list[int] = []
        times: list[float] = []
        for step in steps:
            time = scalar_attribute(step, "TIME")
            observed = int(step["q"].shape[0])
            expected = int(np.searchsorted(births, time + 1.0e-18, side="right"))
            times.append(time)
            observed_counts.append(observed)
            expected_counts.append(expected)

            species = np.unique(step["sp"][:])
            if not np.array_equal(species, np.array([expected_species], dtype=species.dtype)):
                raise ValueError(f"{path}:{step.name}: species {species}, expected {expected_species}")
            charges = step["q"][:]
            if not np.all(np.sign(charges) == expected_charge_sign):
                raise ValueError(f"{path}:{step.name}: wrong particle-charge sign")

        if observed_counts != expected_counts:
            raise ValueError(
                f"{path}: cumulative emitted counts differ from exact birth schedule: "
                f"observed={observed_counts}, expected={expected_counts}"
            )
        if observed_counts[-1] != len(emitted):
            raise ValueError(f"{path}: final count {observed_counts[-1]} != {len(emitted)}")

        final_q = steps[-1]["q"][:]
        return {
            "h5_steps": len(steps),
            "first_global_step": int(scalar_attribute(steps[0], "GlobalTrackStep")),
            "last_global_step": int(scalar_attribute(steps[-1], "GlobalTrackStep")),
            "first_dump_time_s": times[0],
            "last_dump_time_s": times[-1],
            "first_count": observed_counts[0],
            "final_count": observed_counts[-1],
            "species_code": expected_species,
            "particle_charge_C": float(final_q[0]),
            "birth_schedule_exact": True,
        }


def validate_primary_and_field(opalx_dir: Path) -> dict[str, object]:
    primary_path = opalx_dir / "gamma_gamma_pairs_timed_c0.h5"
    with h5py.File(primary_path, "r") as h5:
        steps = numeric_steps(h5)
        primary_q = steps[-1]["q"][:]
        primary_count = int(primary_q.size)
        primary_total_charge = float(primary_q.sum())
        if primary_count != 10_000 or not np.all(primary_q < 0.0):
            raise ValueError(f"{primary_path}: invalid primary container")

    rho_path = opalx_dir / "data" / "gamma_gamma_pairs_timed-RHO_scalar-beambeam_rho_pre.h5"
    with h5py.File(rho_path, "r") as h5:
        steps = numeric_steps(h5)
        source_charges = np.array(
            [scalar_attribute(step, "particle_total_charge") for step in steps]
        )
        all_container_counts = np.array(
            [int(scalar_attribute(step, "particle_total_num")) for step in steps]
        )
        np.testing.assert_allclose(source_charges, primary_total_charge, rtol=5.0e-12, atol=1.0e-24)
        if all_container_counts[-1] != 10_000 + 2 * 1_297:
            raise ValueError(f"{rho_path}: final all-container count is {all_container_counts[-1]}")

    return {
        "primary_count": primary_count,
        "primary_total_charge_C": primary_total_charge,
        "field_diagnostic_steps": len(source_charges),
        "all_container_count_initial": int(all_container_counts[0]),
        "all_container_count_final": int(all_container_counts[-1]),
        "source_charge_min_C": float(source_charges.min()),
        "source_charge_max_C": float(source_charges.max()),
        "source_charge_unchanged_as_witnesses_appear": True,
    }


def validate(opalx_dir: Path) -> dict[str, object]:
    input_dir = opalx_dir / "input"
    return {
        "t0_s": expected_t0(),
        "primary_and_field": validate_primary_and_field(opalx_dir),
        "electron_container": validate_witness(
            opalx_dir / "gamma_gamma_pairs_timed_c1.h5",
            input_dir / "cain_electrons.emittedfromfile",
            expected_species=1,
            expected_charge_sign=-1,
        ),
        "positron_container": validate_witness(
            opalx_dir / "gamma_gamma_pairs_timed_c2.h5",
            input_dir / "cain_positrons.emittedfromfile",
            expected_species=2,
            expected_charge_sign=1,
        ),
    }


def parse_args() -> argparse.Namespace:
    workflow = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx-dir", type=Path, default=workflow / "opalx")
    parser.add_argument(
        "--report", type=Path, default=workflow / "results" / "opalx_output_validation.json"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = validate(args.opalx_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report: {args.report}")


if __name__ == "__main__":
    main()
