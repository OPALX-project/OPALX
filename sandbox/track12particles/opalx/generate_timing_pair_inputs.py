#!/usr/bin/env python3
"""Generate pair-wise OPALX inputs for the slide insertion-time convention.

The source slides define the witness ``ct`` values as insertion times, not as
one simultaneous spatial snapshot.  OPALX has one simulation clock per run, so
the clean way to keep the original witness coordinates and still evaluate the
primary bunch at each insertion time is to run one electron/positron pair per
input deck and shift only the primary source centroid.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ELECTRON_REST_EV = 510_998.95
DEFAULT_INITIAL_CSV = Path(__file__).resolve().parents[1] / "outputs" / "track12_initial_conditions.csv"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "timing_pairs"


INPUT_TEMPLATE = """OPTION, PSDUMPFREQ = 1001;
OPTION, STATDUMPFREQ = 1001;
OPTION, BOUNDPDESTROY = 10;
OPTION, VERSION = 10900;
OPTION, SEED = 20260629;

TITLE, STRING = "slide-timed one-source OPALX BeamBeam pair {pair} 400k/{steps}-step";

// Pair-wise timing diagnostic for TestParticleOrbitSimulation.pptx:
//   - witness coordinates are the numeric TestParticleOrbit initial values
//   - the primary source centroid is shifted to the oncoming-bunch position at
//     this pair's insertion time, ct = {ct_m:.16e} m
//   - COPY_TIME = 0 disables the mirrored copied-source field
REAL primary_kinetic_energy = 0.245;
REAL primary_gamma = (primary_kinetic_energy + EMASS) / EMASS;
REAL primary_beta = sqrt(1.0 - 1.0 / (primary_gamma^2));
REAL primary_p0 = primary_gamma * primary_beta * EMASS;
REAL primary_beta_gamma = primary_p0 / EMASS;

REAL primary_sigma_xy = 1.944325075701e-6;
REAL primary_sigma_z = 0.6e-3;
REAL primary_sigma_xyp = 1.0e-12;
REAL primary_sigma_pxy = primary_beta_gamma * primary_sigma_xyp;

REAL primary_electrons_per_bunch = 1.25e10;
REAL primary_macroparticles = 400000;
REAL elementary_charge = 1.602176634e-19;
REAL primary_bunch_charge = primary_electrons_per_bunch * elementary_charge;

REAL n_witness_per_species = 1;
REAL witness_bunch_charge = n_witness_per_species * elementary_charge;
REAL witness_p0 = sqrt(3.0) * EMASS;

REAL bb_length = 0.006;
REAL bb_ip_s = 0.5*bb_length;
REAL witness_ct_m = {ct_m:.16e};
REAL primary_source_r0z = bb_ip_s - witness_ct_m;
REAL dt_track = 1.0e-15;
REAL n_steps = {steps};
REAL NXY = 32;
REAL NZ = 64;

IP1: BEAMBEAM, L = bb_length, ELEMEDGE = 0.0,
    VISUALIZE = FALSE,
    COPY_TIME = 0.0,
    WITNESS_CONTAINERS = "1,2",
    RETIRE_TIME = 1000.0e-12,
    APERTURE = "CIRCLE(0.001)";

track12Line: LINE = (IP1);

SOURCEBINNING: BINNING,
    MAXBINS = 1,
    DESIREDWIDTH = 0.1,
    BINNINGALPHA = 1.0,
    BINNINGBETA = 1.0,
    PARAMETER = "VELOCITYZ",
    ADAPTIVEBINNING = FALSE,
    DUMPBINSFILE = "NONE",
    TABLEPRINTFREQ = 0;

FS1: FIELDSOLVER, TYPE = OPEN,
    BINS = "SOURCEBINNING",
    NX = NXY, NY = NXY, NZ = NZ,
    PARFFTX = FALSE, PARFFTY = FALSE, PARFFTZ = TRUE,
    BCFFTX = OPEN, BCFFTY = OPEN, BCFFTZ = OPEN,
    BBOXINCR = 1,
    GREENSF = STANDARD;

DistPrimaryElectrons: DISTRIBUTION, TYPE = GAUSS,
    SIGMAX = primary_sigma_xy, SIGMAY = primary_sigma_xy, SIGMAZ = primary_sigma_z,
    SIGMAPX = primary_sigma_pxy, SIGMAPY = primary_sigma_pxy, SIGMAPZ = 1.0e-12,
    NPARTDIST = primary_macroparticles;

DistWitnessElectrons: DISTRIBUTION, TYPE = FROMFILE,
    FNAME = "{electron_fromfile}",
    NPARTDIST = n_witness_per_species;

DistWitnessPositrons: DISTRIBUTION, TYPE = FROMFILE,
    FNAME = "{positron_fromfile}",
    NPARTDIST = n_witness_per_species;

SourcePrimaryElectrons: EMISSIONSOURCE, DISTRIBUTION = DistPrimaryElectrons, R0Z = primary_source_r0z;
SourceWitnessElectrons: EMISSIONSOURCE, DISTRIBUTION = DistWitnessElectrons, R0Z = bb_ip_s;
SourceWitnessPositrons: EMISSIONSOURCE, DISTRIBUTION = DistWitnessPositrons, R0Z = bb_ip_s;

PrimaryElectronSources: EMISSIONSOURCELIST = (SourcePrimaryElectrons);
WitnessElectronSources: EMISSIONSOURCELIST = (SourceWitnessElectrons);
WitnessPositronSources: EMISSIONSOURCELIST = (SourceWitnessPositrons);

PrimaryElectrons: BEAM,
    PARTICLE = ELECTRON,
    PC = primary_p0,
    NALLOC = primary_macroparticles,
    BCHARGE = primary_bunch_charge,
    SOURCES = PrimaryElectronSources;

WitnessElectrons: BEAM,
    PARTICLE = ELECTRON,
    PC = witness_p0,
    NALLOC = n_witness_per_species,
    BCHARGE = witness_bunch_charge,
    SOURCES = WitnessElectronSources;

WitnessPositrons: BEAM,
    PARTICLE = POSITRON,
    PC = witness_p0,
    NALLOC = n_witness_per_species,
    BCHARGE = witness_bunch_charge,
    SOURCES = WitnessPositronSources;

TRACK, LINE = track12Line,
    BEAMS = {{PrimaryElectrons, WitnessElectrons, WitnessPositrons}},
    MAXSTEPS = n_steps,
    ZSTOP = bb_length,
    DT = dt_track;
 RUN, METHOD = "PARALLEL", FIELDSOLVER = FS1;
ENDTRACK;

QUIT;
"""


def read_initial_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_fromfile(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("1\n")
        stream.write("x y z px py pz\n")
        stream.write(
            f"{float(row['x']):.16e} {float(row['y']):.16e} {float(row['s']):.16e} "
            f"{float(row['Px']) / ELECTRON_REST_EV:.16e} "
            f"{float(row['Py']) / ELECTRON_REST_EV:.16e} "
            f"{float(row['Ps']) / ELECTRON_REST_EV:.16e}\n"
        )


def write_metadata(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "pair",
        "ct_m",
        "ct_over_sigma_z",
        "primary_source_r0z_m",
        "witness_r0z_m",
        "witness_minus_source_m",
        "witness_minus_source_over_sigma_z",
        "electron_name",
        "positron_name",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-csv", type=Path, default=DEFAULT_INITIAL_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sigma-z-m", type=float, default=0.6e-3)
    parser.add_argument("--bb-ip-s-m", type=float, default=0.003)
    parser.add_argument("--steps", type=int, default=10)
    args = parser.parse_args()

    rows = read_initial_rows(args.initial_csv)
    by_pair: dict[int, dict[str, dict[str, str]]] = {}
    for row in rows:
        pair = int(row["pair"])
        by_pair.setdefault(pair, {})[row["species"]] = row

    metadata_rows: list[dict[str, str]] = []
    for pair in sorted(by_pair):
        species_rows = by_pair[pair]
        electron = species_rows["electron"]
        positron = species_rows["positron"]
        ct_m = float(electron["t"])
        if abs(ct_m - float(positron["t"])) > 1.0e-15:
            raise ValueError(f"pair {pair}: electron and positron insertion times differ")

        pair_dir = args.output_dir / f"pair{pair}"
        electron_fromfile = pair_dir / f"track12_pair{pair}_electron.fromfile"
        positron_fromfile = pair_dir / f"track12_pair{pair}_positron.fromfile"
        write_fromfile(electron_fromfile, electron)
        write_fromfile(positron_fromfile, positron)

        input_name = f"track12_pair{pair}_slide_timed_one_source_1fs_400k_{args.steps}steps.in"
        input_path = pair_dir / input_name
        input_path.write_text(
            INPUT_TEMPLATE.format(
                pair=pair,
                ct_m=ct_m,
                steps=args.steps,
                electron_fromfile=electron_fromfile.name,
                positron_fromfile=positron_fromfile.name,
            ),
            encoding="utf-8",
            newline="\n",
        )

        primary_r0z = args.bb_ip_s_m - ct_m
        witness_minus_source = 2.0 * ct_m
        metadata_rows.append(
            {
                "pair": str(pair),
                "ct_m": f"{ct_m:.16e}",
                "ct_over_sigma_z": f"{ct_m / args.sigma_z_m:.16e}",
                "primary_source_r0z_m": f"{primary_r0z:.16e}",
                "witness_r0z_m": f"{args.bb_ip_s_m:.16e}",
                "witness_minus_source_m": f"{witness_minus_source:.16e}",
                "witness_minus_source_over_sigma_z": f"{witness_minus_source / args.sigma_z_m:.16e}",
                "electron_name": electron["name"],
                "positron_name": positron["name"],
            }
        )

    write_metadata(args.output_dir / "slide_timing_metadata.csv", metadata_rows)
    print(f"wrote {len(metadata_rows)} pair inputs under {args.output_dir}")
    print(f"wrote {args.output_dir / 'slide_timing_metadata.csv'}")


if __name__ == "__main__":
    main()
