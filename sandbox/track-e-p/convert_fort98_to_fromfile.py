#!/usr/bin/env python3
"""Convert gamma-gamma pair output into OPALX FROMFILE distributions.

The source file is expected to contain the columns

    species weight ct x y z energy px py pz

where species 2 is an electron and species 3 is a positron.  By default this
converter writes only the phase-space columns consumed by the OPALX FROMFILE
reader.  With ``--add-time`` it also writes the fort.98 creation time as an
extra ``t`` column, adds a constant ``bin_number`` column with value 1, and
appends ``-t`` to each output filename.  The fort.98 column is ``ct`` in metres,
so the written time is ``t = ct / c`` in seconds.
OPALX stores particle-container momenta as beta*gamma = p / (m c).  The input
momenta are in eV/c, so the output momenta are divided by the electron rest
energy in eV.
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
SPEED_OF_LIGHT_M_PER_S = 299792458.0


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


def time_filename(filename: str) -> str:
    """Return the time-carrying variant of a FROMFILE filename."""
    path = Path(filename)
    return f"{path.stem}-t{path.suffix}"


def write_fromfile(data: pd.DataFrame, path: Path, add_time: bool = False) -> None:
    """Write one OPALX FROMFILE distribution."""
    path.parent.mkdir(parents=True, exist_ok=True)

    phase_space = data.loc[:, FROMFILE_COLUMNS].copy()
    phase_space.loc[:, ["px", "py", "pz"]] /= ELECTRON_REST_ENERGY_EV
    if add_time:
        phase_space.loc[:, "t"] = data["ct"] / SPEED_OF_LIGHT_M_PER_S
        phase_space.loc[:, "bin_number"] = 1

    with path.open("w", encoding="utf-8") as out:
        out.write(f"{len(phase_space)}\n")
        out.write(" ".join(phase_space.columns) + "\n")
        phase_space.to_csv(
            out,
            sep=" ",
            header=False,
            index=False,
            float_format="%.12e",
        )


def convert(input_path: Path, output_dir: Path, add_time: bool = False) -> dict[str, int]:
    """Split the input by species and write the corresponding OPALX files."""
    data = read_fort98(input_path)
    counts: dict[str, int] = {}

    for species_id, (name, filename) in SPECIES.items():
        species_data = data[data["species"] == species_id]
        if add_time:
            filename = time_filename(filename)
        output_path = output_dir / filename
        write_fromfile(species_data, output_path, add_time=add_time)
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
    parser.add_argument(
        "--add-time",
        action="store_true",
        help=(
            "Also write the fort.98 creation time as t=ct/c in seconds and a "
            "constant bin_number=1 column. The output filenames receive a "
            "'-t' suffix before '.fromfile'."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts = convert(args.input, args.output_dir, add_time=args.add_time)

    for output_path, count in counts.items():
        print(f"{output_path}: {count} particles")


if __name__ == "__main__":
    main()
