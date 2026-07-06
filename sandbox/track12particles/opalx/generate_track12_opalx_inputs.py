#!/usr/bin/env python3
"""Generate OPALX FROMFILE witness inputs from CAIN TestParticleOrbit.dat.

The OPALX FROMFILE distribution uses columns

    x y z px py pz

where positions are in meters and momenta are normalized beta*gamma.  The CAIN
file uses x, y, s in meters and momenta in eV/c.  For electrons and positrons
we use the first row of each named trajectory as the injected witness state and
map s -> z, P/(m_e c) -> beta*gamma.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path

ELECTRON_REST_EV = 510_998.95
DEFAULT_CAIN_FILE = Path(__file__).resolve().parents[2] / "TestParticleOrbit.dat"
DEFAULT_INITIAL_CSV = Path(__file__).resolve().parents[1] / "outputs" / "track12_initial_conditions.csv"
TRACK12_SCRIPT = Path(__file__).resolve().parents[1] / "track12particles.py"


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def first_rows_by_particle(dataframe):
    return (
        dataframe.sort_values(["name", "t"])
        .groupby("name", sort=False, as_index=False)
        .first()
        .sort_values(["kind", "pair"])
    )


def load_initial_rows_from_csv(path: Path) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    numeric = {
        "pair",
        "kind",
        "t",
        "x",
        "y",
        "s",
        "E",
        "Px",
        "Py",
        "Ps",
    }
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        for raw in reader:
            row: dict[str, float | int | str] = dict(raw)
            for key in numeric:
                if key in row:
                    row[key] = float(row[key])
            row["pair"] = int(row["pair"])
            row["kind"] = int(row["kind"])
            rows.append(row)
    return sorted(rows, key=lambda row: (int(row["kind"]), int(row["pair"])))


def write_fromfile(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write(f"{len(rows)}\n")
        stream.write("x y z px py pz\n")
        for row in rows:
            stream.write(
                f"{row['x']:.16e} {row['y']:.16e} {row['s']:.16e} "
                f"{row['Px'] / ELECTRON_REST_EV:.16e} "
                f"{row['Py'] / ELECTRON_REST_EV:.16e} "
                f"{row['Ps'] / ELECTRON_REST_EV:.16e}\n"
            )


def write_metadata(path: Path, rows, source_file: Path) -> None:
    fieldnames = [
        "name",
        "pair",
        "kind",
        "species",
        "source_file",
        "source_mtime",
        "source_date",
        "t_m",
        "x_m",
        "y_m",
        "s_m",
        "E_eV",
        "Px_eVc",
        "Py_eVc",
        "Ps_eVc",
        "px_beta_gamma",
        "py_beta_gamma",
        "pz_beta_gamma",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "name": row["name"],
                    "pair": int(row["pair"]),
                    "kind": int(row["kind"]),
                    "species": row["species"],
                    "source_file": str(source_file.resolve()),
                    "source_mtime": row["source_mtime"],
                    "source_date": row["source_date"],
                    "t_m": f"{row['t']:.16e}",
                    "x_m": f"{row['x']:.16e}",
                    "y_m": f"{row['y']:.16e}",
                    "s_m": f"{row['s']:.16e}",
                    "E_eV": f"{row['E']:.16e}",
                    "Px_eVc": f"{row['Px']:.16e}",
                    "Py_eVc": f"{row['Py']:.16e}",
                    "Ps_eVc": f"{row['Ps']:.16e}",
                    "px_beta_gamma": f"{row['Px'] / ELECTRON_REST_EV:.16e}",
                    "py_beta_gamma": f"{row['Py'] / ELECTRON_REST_EV:.16e}",
                    "pz_beta_gamma": f"{row['Ps'] / ELECTRON_REST_EV:.16e}",
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_CAIN_FILE)
    parser.add_argument(
        "--initial-csv",
        type=Path,
        default=None,
        help="Use an existing track12_initial_conditions.csv instead of parsing the CAIN .dat file.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_file = args.input
    if args.initial_csv is not None:
        rows = load_initial_rows_from_csv(args.initial_csv)
    else:
        try:
            module = load_track12_module()
            dataframe = module.parse_reference_file(args.input)
            rows = first_rows_by_particle(dataframe).to_dict("records")
        except PermissionError:
            if not DEFAULT_INITIAL_CSV.exists():
                raise
            print(
                f"warning: could not read {args.input}; falling back to {DEFAULT_INITIAL_CSV}"
            )
            rows = load_initial_rows_from_csv(DEFAULT_INITIAL_CSV)

    electrons = [row for row in rows if row["species"] == "electron"]
    positrons = [row for row in rows if row["species"] == "positron"]

    write_fromfile(args.output_dir / "track12_electrons.fromfile", electrons)
    write_fromfile(args.output_dir / "track12_positrons.fromfile", positrons)
    write_metadata(
        args.output_dir / "track12_witness_metadata.csv",
        electrons + positrons,
        source_file,
    )

    print(f"wrote {len(electrons)} electron witnesses")
    print(f"wrote {len(positrons)} positron witnesses")
    print(f"source: {args.input.resolve()}")


if __name__ == "__main__":
    main()
