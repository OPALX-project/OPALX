#!/usr/bin/env python3
"""Convert CAIN gamma-gamma pairs to OPALX explicit-birth input files.

CAIN ``fort.98`` rows are expected to contain

    species weight ct x y z energy px py pz

with positions and ``ct`` in metres and momenta/energy in eV/c and eV.  The
OPALX files contain named columns

    x y z px py pz birth_time

where positions are metres relative to ``EMISSIONSOURCE::R0``, momenta are
``p/(m_e c)``, and ``birth_time = ct/c`` is seconds relative to
``EMISSIONSOURCE::T0``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


COLUMNS = ["species", "weight", "ct", "x", "y", "z", "energy", "px", "py", "pz"]
OUTPUT_COLUMNS = ["x", "y", "z", "px", "py", "pz", "birth_time"]
SPECIES = {
    2: ("electron", "cain_electrons.emittedfromfile"),
    3: ("positron", "cain_positrons.emittedfromfile"),
}
ELECTRON_REST_ENERGY_EV = 510_998.95
SPEED_OF_LIGHT_M_PER_S = 299_792_458.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_cain(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep=r"\s+", names=COLUMNS, comment="#", engine="python")
    if data.empty:
        raise ValueError(f"No CAIN particles found in {path}")
    if not np.isfinite(data.to_numpy(dtype=float)).all():
        raise ValueError(f"CAIN input {path} contains a non-finite value")

    data["species"] = data["species"].astype(int)
    unknown = sorted(set(data["species"]) - set(SPECIES))
    if unknown:
        raise ValueError(f"Unsupported CAIN species IDs in {path}: {unknown}")
    return data


def convert_species(data: pd.DataFrame) -> pd.DataFrame:
    converted = data.loc[:, ["x", "y", "z", "px", "py", "pz"]].copy()
    converted.loc[:, ["px", "py", "pz"]] /= ELECTRON_REST_ENERGY_EV
    converted["birth_time"] = data["ct"] / SPEED_OF_LIGHT_M_PER_S
    return converted.loc[:, OUTPUT_COLUMNS]


def write_emitted(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(data)}\n")
        stream.write(" ".join(OUTPUT_COLUMNS) + "\n")
        data.to_csv(stream, sep=" ", header=False, index=False, float_format="%.16e")


def convert(input_path: Path, output_dir: Path, manifest_path: Path) -> dict[str, object]:
    source = read_cain(input_path)
    manifest: dict[str, object] = {
        "format": "OPALX explicit emitted spacetime v1",
        "source": str(input_path.resolve()),
        "source_sha256": sha256(input_path),
        "units": {
            "x_y_z": "m relative to EMISSIONSOURCE R0",
            "px_py_pz": "p/(m_e c)",
            "birth_time": "s relative to EMISSIONSOURCE T0; CAIN ct/c",
        },
        "c_m_per_s": SPEED_OF_LIGHT_M_PER_S,
        "electron_rest_energy_eV": ELECTRON_REST_ENERGY_EV,
        "cain_weight_policy": "recorded in source only; not written because witnesses do not deposit",
        "species": {},
    }

    species_manifest: dict[str, object] = {}
    for species_id, (name, filename) in SPECIES.items():
        selected = source.loc[source["species"] == species_id].copy()
        converted = convert_species(selected)
        output_path = output_dir / filename
        write_emitted(converted, output_path)
        species_manifest[name] = {
            "cain_species_id": species_id,
            "count": int(len(selected)),
            "output": str(output_path.resolve()),
            "output_sha256": sha256(output_path),
            "birth_time_min_s": float(converted["birth_time"].min()),
            "birth_time_max_s": float(converted["birth_time"].max()),
            "z_min_m": float(converted["z"].min()),
            "z_max_m": float(converted["z"].max()),
            "cain_weight_min": float(selected["weight"].min()),
            "cain_weight_max": float(selected["weight"].max()),
        }
    manifest["species"] = species_manifest

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def parse_args() -> argparse.Namespace:
    workflow = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=workflow.parent / "track-e-p" / "fort98.txt",
        help="CAIN fort.98 table (default: sandbox/track-e-p/fort98.txt)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workflow / "opalx" / "input",
        help="Directory for generated OPALX emitted files",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=workflow / "results" / "conversion_manifest.json",
        help="Conversion manifest path",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = convert(args.input, args.output_dir, args.manifest)
    for name, details in manifest["species"].items():
        print(f"{name}: {details['count']} -> {details['output']}")
    print(f"manifest: {args.manifest}")


if __name__ == "__main__":
    main()
