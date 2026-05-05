#!/usr/bin/env python3
"""Convert gamma-gamma pair output into OPALX FROMFILE distributions.

The source file is expected to contain the columns

    species weight ct x y z energy px py pz

where species 2 is an electron and species 3 is a positron.  The current
OPALX FROMFILE reader consumes only phase-space coordinates, so this converter
intentionally ignores the particle weight, creation time, and energy.  OPALX
stores particle-container momenta as beta*gamma = p / (m c).  The input momenta
are in eV/c, so the output momenta are divided by the electron rest energy in eV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SPECIES = {
    2: ("electron", "gamma_gamma_electrons.fromfile"),
    3: ("positron", "gamma_gamma_positrons.fromfile"),
}

COLUMNS = ["species", "weight", "ct", "x", "y", "z", "energy", "px", "py", "pz"]
FROMFILE_COLUMNS = ["x", "y", "z", "px", "py", "pz"]
ELECTRON_REST_ENERGY_EV = 510998.95000


def read_fort98(path: Path) -> pd.DataFrame:
    """Read the fort.98 gamma-gamma particle table into a dataframe."""
    data = pd.read_csv(path, sep=r"\s+", names=COLUMNS, comment="#", engine="python")
    if data.empty:
        raise ValueError(f"No particles found in {path}")

    data["species"] = data["species"].astype(int)
    unknown_species = sorted(set(data["species"]) - set(SPECIES))
    if unknown_species:
        raise ValueError(f"Unsupported species IDs in {path}: {unknown_species}")

    return data


def write_fromfile(data: pd.DataFrame, path: Path) -> None:
    """Write one OPALX FROMFILE distribution."""
    path.parent.mkdir(parents=True, exist_ok=True)

    phase_space = data.loc[:, FROMFILE_COLUMNS].copy()
    phase_space.loc[:, ["px", "py", "pz"]] /= ELECTRON_REST_ENERGY_EV

    with path.open("w", encoding="utf-8") as out:
        out.write(f"{len(phase_space)}\n")
        out.write(" ".join(FROMFILE_COLUMNS) + "\n")
        phase_space.to_csv(
            out,
            sep=" ",
            header=False,
            index=False,
            float_format="%.12e",
        )


def convert(input_path: Path, output_dir: Path) -> dict[str, int]:
    """Split the input by species and write the corresponding OPALX files."""
    data = read_fort98(input_path)
    counts: dict[str, int] = {}

    for species_id, (name, filename) in SPECIES.items():
        species_data = data[data["species"] == species_id]
        output_path = output_dir / filename
        write_fromfile(species_data, output_path)
        counts[str(output_path)] = len(species_data)

    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create electron and positron OPALX FROMFILE distributions from fort98.txt."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("/Users/adelmann/Desktop/fort98.txt"),
        help="Input gamma-gamma particle table.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for generated OPALX FROMFILE distributions.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts = convert(args.input, args.output_dir)

    for output_path, count in counts.items():
        print(f"{output_path}: {count} particles")


if __name__ == "__main__":
    main()
