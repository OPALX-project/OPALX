#!/usr/bin/env python3
from __future__ import annotations

import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

sys.path.insert(0, "/private/tmp/opalx-py")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
BMAD_EXE = Path("/Users/adelmann/git/bmad-opalx/install/bmad-ecosystem/production/bin/beam_track_example")
OPALX_EXE = Path("/Users/adelmann/git/opalx/build/src/opalx")

ELECTRON_MASS_EV = 510_998.95
BEGINNING_E_TOT_EV = 1.00051099895e9
REFERENCE_BETA_GAMMA = math.sqrt((BEGINNING_E_TOT_EV / ELECTRON_MASS_EV) ** 2 - 1.0)
REFERENCE_P0C_EV = REFERENCE_BETA_GAMMA * ELECTRON_MASS_EV

N_PARTICLE = 1000
SIGMA = 1.0e-3
RNG_SEED = 12345
OPALX_DT = 1.0e-11

INPUT_CANONICAL_CSV = ROOT / "bend_bunch_rms_bmad_canonical_input.csv"
BMAD_INPUT_BEAM = ROOT / "bend_bunch_rms_bmad_input.dat"
OPALX_INPUT_DISTRIBUTION = ROOT / "bend_bunch_rms_opalx_distribution.txt"

SUMMARY_CSV = ROOT / "bend_bunch_rms_summary.csv"
SUMMARY_MD = ROOT / "bend_bunch_rms_summary.md"


@dataclass(frozen=True)
class Case:
    family: str
    angle: float
    angle_label: str
    suffix: str
    e1: float = 0.0

    @property
    def e1_label(self) -> str:
        return "0"

    @property
    def bmad_stem(self) -> str:
        return f"bend_bunch_rms_m_{self.suffix}"

    @property
    def opalx_stem(self) -> str:
        return f"bend_bunch_rms_o_{self.suffix}"


CASES = [
    Case("SBEND", math.pi / 4.0, "+pi/4", "sb45p"),
    Case("SBEND", -math.pi / 4.0, "-pi/4", "sb45m"),
    Case("SBEND", math.pi / 2.0, "+pi/2", "sb90p"),
    Case("SBEND", -math.pi / 2.0, "-pi/2", "sb90m"),
    Case("RBEND", math.pi / 4.0, "+pi/4", "rb45p"),
    Case("RBEND", -math.pi / 4.0, "-pi/4", "rb45m"),
    Case("RBEND", math.pi / 2.0, "+pi/2", "rb90p"),
    Case("RBEND", -math.pi / 2.0, "-pi/2", "rb90m"),
]


def angle_expr(angle: float) -> str:
    if math.isclose(angle, math.pi / 4.0):
        return "PI / 4"
    if math.isclose(angle, -math.pi / 4.0):
        return "-PI / 4"
    if math.isclose(angle, math.pi / 2.0):
        return "PI / 2"
    if math.isclose(angle, -math.pi / 2.0):
        return "-PI / 2"
    raise ValueError(f"unsupported angle {angle}")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def rows_to_markdown(df: pd.DataFrame) -> str:
    def fmt(value: object) -> str:
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.5f}"
        return str(value)

    headers = list(df.columns)
    body = [[fmt(value) for value in row] for row in df.itertuples(index=False, name=None)]
    widths = [len(header) for header in headers]
    for row in body:
        for i, value in enumerate(row):
            widths[i] = max(widths[i], len(value))

    def render(values: Iterable[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[i]) for i, value in enumerate(values)) + " |"

    lines = [render(headers), render(["-" * width for width in widths])]
    lines.extend(render(row) for row in body)
    return "\n".join(lines) + "\n"


def generate_canonical_bunch() -> pd.DataFrame:
    """
    Generate a Gaussian bunch in BMAD canonical coordinates.

    The sampled phase-space coordinates are

    .. math::

       (x, p_x, y, p_y, z, p_z),

    with zero mean and common sigma :math:`10^{-3}`.
    """

    rng = np.random.default_rng(RNG_SEED)
    data = rng.normal(loc=0.0, scale=SIGMA, size=(N_PARTICLE, 6))
    return pd.DataFrame(data, columns=["x", "px", "y", "py", "z", "pz"])


def canonical_to_opalx(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert BMAD canonical coordinates into OPALX beta-gamma momenta.

    For reference beta-gamma :math:`(\beta\gamma)_0`,

    .. math::

       (\beta\gamma)_x &= (\beta\gamma)_0 p_x, \\
       (\beta\gamma)_y &= (\beta\gamma)_0 p_y, \\
       (\beta\gamma)_s &= (\beta\gamma)_0
         \sqrt{(1 + p_z)^2 - p_x^2 - p_y^2}.
    """

    px = df["px"].to_numpy()
    py = df["py"].to_numpy()
    pz = df["pz"].to_numpy()
    ps = np.sqrt((1.0 + pz) ** 2 - px**2 - py**2)
    return pd.DataFrame(
        {
            "x": df["x"].to_numpy(),
            "px": REFERENCE_BETA_GAMMA * px,
            "y": df["y"].to_numpy(),
            "py": REFERENCE_BETA_GAMMA * py,
            "z": df["z"].to_numpy(),
            "pz": REFERENCE_BETA_GAMMA * ps,
        }
    )


def canonical_to_betagamma(df: pd.DataFrame, p0c_ev: float) -> pd.DataFrame:
    """Convert BMAD canonical coordinates to beta-gamma coordinates."""

    bg0 = p0c_ev / ELECTRON_MASS_EV
    px = df["px"].to_numpy()
    py = df["py"].to_numpy()
    pz = df["pz"].to_numpy()
    ps = np.sqrt((1.0 + pz) ** 2 - px**2 - py**2)
    return pd.DataFrame(
        {
            "x": df["x"].to_numpy(),
            "y": df["y"].to_numpy(),
            "z": df["z"].to_numpy(),
            "px": bg0 * px,
            "py": bg0 * py,
            "ps": bg0 * ps,
        }
    )


def write_bmad_beam_file(df: pd.DataFrame, path: Path) -> None:
    lines = [
        "# species = electron",
        "# charge_tot = 0.0",
        "# s_position = 0.0",
        "# time = 0.0",
        "#! index x px y py z pz state",
    ]
    for idx, row in enumerate(df.itertuples(index=False), start=1):
        lines.append(
            f"{idx} "
            f"{row.x:.16e} {row.px:.16e} {row.y:.16e} {row.py:.16e} "
            f"{row.z:.16e} {row.pz:.16e} Alive"
        )
    write_text(path, "\n".join(lines) + "\n")


def write_opalx_distribution(df: pd.DataFrame, path: Path) -> None:
    lines = [str(len(df)), "x px y py z pz"]
    for row in df.itertuples(index=False):
        lines.append(
            f"{row.x:.16e} {row.px:.16e} {row.y:.16e} {row.py:.16e} "
            f"{row.z:.16e} {row.pz:.16e}"
        )
    write_text(path, "\n".join(lines) + "\n")


def make_bmad_lattice(case: Case) -> str:
    if case.family == "SBEND":
        bend_def = (
            f"b1: sbend, l = 1.0, angle = {angle_expr(case.angle).lower()}, e1 = 0.0, e2 = 0.0, &\n"
            "    fringe_type = none, tracking_method = runge_kutta, num_steps = 800"
        )
    else:
        bend_def = (
            f"b1: rbend, l_rectangle = 1.0, angle = {angle_expr(case.angle).lower()}, e1 = 0.0, e2 = 0.0, &\n"
            "    fiducial_pt = entrance_end, fringe_type = none, tracking_method = runge_kutta, num_steps = 800"
        )
    return f"""parameter[particle] = electron
parameter[geometry] = open

beginning[e_tot] = 1.00051099895e9
beginning[beta_a] = 1.0
beginning[beta_b] = 1.0
beginning[alpha_a] = 0.0
beginning[alpha_b] = 0.0

{bend_def}

line1: line = (b1)
use, line1

no_digested
end_file
"""


def make_bmad_track_input(case: Case, output_beam_name: str) -> str:
    return f"""&beam_track_params
    lat_filename = '{case.bmad_stem}.bmad'
    outfile_name = '{output_beam_name}'
    beam_init%position_file = '{BMAD_INPUT_BEAM.name}'
    beam_init%n_particle = 0
    beam_init%n_bunch = 0
    ran_seed = {RNG_SEED}
/
"""


def make_opalx_deck(case: Case) -> str:
    return f"""OPTION, PSDUMPFREQ = 50000;
OPTION, STATDUMPFREQ = 1;
OPTION, BOUNDPDESTROY = 10;
OPTION, VERSION = 10900;

Title, string = "{case.opalx_stem} {case.family} ANGLE={case.angle_label} E1=0 PSI=0";

B1: {case.family},
    L = 1.0,
    ANGLE = {angle_expr(case.angle)},
    E1 = 0.0,
    PSI = 0.0,
    ELEMEDGE = 0.0;

Line1: Line = (B1);

Dist0: DISTRIBUTION,
    TYPE = FROMFILE,
    FNAME = "{OPALX_INPUT_DISTRIBUTION.name}",
    NPARTDIST = {N_PARTICLE};

ES0: EMISSIONSOURCE, DISTRIBUTION = Dist0;
Sources0: EMISSIONSOURCELIST = (ES0);

FS0: FIELDSOLVER,
    TYPE = NONE,
    NX = 16,
    NY = 16,
    NZ = 16,
    PARFFTX = false,
    PARFFTY = false,
    PARFFTZ = true,
    BCFFTX = open,
    BCFFTY = open,
    BCFFTZ = open,
    BBOXINCR = 1,
    GREENSF = INTEGRATED;

BEAM0: BEAM,
    PARTICLE = ELECTRON,
    NALLOC = {N_PARTICLE},
    BCHARGE = 1.6021766339999998e-19,
    SOURCES = Sources0,
    CHARGE = -1;

TRACK,
    LINE = Line1,
    BEAM = BEAM0,
    MAXSTEPS = 20000,
    DT = {OPALX_DT:.1e},
    ZSTOP = 1.5;

RUN,
    METHOD = "PARALLEL",
    FIELDSOLVER = FS0;

ENDTRACK;

QUIT;
"""


def run_logged(cmd: list[str], log_path: Path) -> None:
    proc = subprocess.run(
        cmd,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log_path.write_text(proc.stdout, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}\nSee {log_path}")


def parse_bmad_beam_file(path: Path) -> tuple[pd.DataFrame, float]:
    columns: list[str] = []
    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#!"):
            columns = line[2:].split()
            continue
        if not line or line.startswith("#"):
            continue
        rows.append(line.split())
    if not columns or not rows:
        raise RuntimeError(f"could not parse BMAD beam file: {path}")
    df = pd.DataFrame(rows, columns=columns)
    numeric_cols = [col for col in columns if col not in {"state", "location"}]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="raise")
    p0c = float(df["p0c"].iloc[0])
    return df, p0c


def read_opalx_stat_row(path: Path) -> dict[str, float]:
    columns: list[str] = []
    data_row: list[str] | None = None
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("&column"):
            name = ""
            i += 1
            while i < len(lines) and not lines[i].startswith("&end"):
                if "name=" in lines[i]:
                    name = lines[i].split("name=")[1].split(",")[0].strip()
                i += 1
            columns.append(name)
        elif line and not line.startswith(("SDDS1", "&", "!")) and columns:
            data_row = line.split()
        i += 1
    if data_row is None or len(data_row) < len(columns):
        raise RuntimeError(f"could not parse OPALX stat file: {path}")
    return {name: float(value) for name, value in zip(columns, data_row)}


def rms_summary_from_particles(df: pd.DataFrame) -> dict[str, float]:
    return {
        "rms_x": float(df["x"].std(ddof=0)),
        "rms_y": float(df["y"].std(ddof=0)),
        "rms_z": float(df["z"].std(ddof=0)),
        "rms_px": float(df["px"].std(ddof=0)),
        "rms_py": float(df["py"].std(ddof=0)),
        "rms_ps": float(df["ps"].std(ddof=0)),
    }


def entrance_summary(canonical_df: pd.DataFrame) -> dict[str, float]:
    return rms_summary_from_particles(canonical_to_betagamma(canonical_df, REFERENCE_P0C_EV))


def run_case(case: Case) -> dict[str, float | str]:
    bmad_lattice = ROOT / f"{case.bmad_stem}.bmad"
    bmad_input = ROOT / f"{case.bmad_stem}.in"
    bmad_output = ROOT / f"{case.bmad_stem}_output.dat"
    bmad_log = ROOT / f"{case.bmad_stem}.log"
    opalx_deck = ROOT / f"{case.opalx_stem}.in"
    opalx_log = ROOT / f"{case.opalx_stem}.log"
    opalx_stat = ROOT / f"{case.opalx_stem}.stat"

    write_text(bmad_lattice, make_bmad_lattice(case))
    write_text(bmad_input, make_bmad_track_input(case, bmad_output.name))
    write_text(opalx_deck, make_opalx_deck(case))

    run_logged([str(BMAD_EXE), bmad_input.name], bmad_log)
    run_logged([str(OPALX_EXE), "--info", "2", opalx_deck.name], opalx_log)

    bmad_df, p0c_exit = parse_bmad_beam_file(bmad_output)
    bmad_exit = rms_summary_from_particles(
        canonical_to_betagamma(bmad_df[["x", "px", "y", "py", "z", "pz"]], p0c_exit)
    )
    opalx_row = read_opalx_stat_row(opalx_stat)
    opalx_exit = {
        "rms_x": opalx_row["rms_x"],
        "rms_y": opalx_row["rms_y"],
        "rms_z": opalx_row["rms_s"],
        "rms_px": opalx_row["rms_px"],
        "rms_py": opalx_row["rms_py"],
        "rms_ps": opalx_row["rms_ps"],
    }

    return {
        "family": case.family,
        "ANGLE": case.angle_label,
        "E1": case.e1_label,
        "bmad_exit_rms_x": bmad_exit["rms_x"],
        "bmad_exit_rms_y": bmad_exit["rms_y"],
        "bmad_exit_rms_z": bmad_exit["rms_z"],
        "bmad_exit_rms_px": bmad_exit["rms_px"],
        "bmad_exit_rms_py": bmad_exit["rms_py"],
        "bmad_exit_rms_ps": bmad_exit["rms_ps"],
        "opalx_exit_rms_x": opalx_exit["rms_x"],
        "opalx_exit_rms_y": opalx_exit["rms_y"],
        "opalx_exit_rms_z": opalx_exit["rms_z"],
        "opalx_exit_rms_px": opalx_exit["rms_px"],
        "opalx_exit_rms_py": opalx_exit["rms_py"],
        "opalx_exit_rms_ps": opalx_exit["rms_ps"],
    }


def main() -> None:
    canonical_df = generate_canonical_bunch()
    canonical_df.to_csv(INPUT_CANONICAL_CSV, index=False)
    write_bmad_beam_file(canonical_df, BMAD_INPUT_BEAM)
    write_opalx_distribution(canonical_to_opalx(canonical_df), OPALX_INPUT_DISTRIBUTION)

    start = entrance_summary(canonical_df)
    rows = []
    for case in CASES:
        row = run_case(case)
        rows.append(
            {
                "family_code": f"{case.family}-O",
                "ANGLE": row["ANGLE"],
                "E1": row["E1"],
                "entrance_rms_x": start["rms_x"],
                "entrance_rms_y": start["rms_y"],
                "entrance_rms_z": start["rms_z"],
                "entrance_rms_px": start["rms_px"],
                "entrance_rms_py": start["rms_py"],
                "entrance_rms_ps": start["rms_ps"],
                "rms_x": row["opalx_exit_rms_x"],
                "rms_y": row["opalx_exit_rms_y"],
                "rms_z": row["opalx_exit_rms_z"],
                "rms_px": row["opalx_exit_rms_px"],
                "rms_py": row["opalx_exit_rms_py"],
                "rms_ps": row["opalx_exit_rms_ps"],
            }
        )
        rows.append(
            {
                "family_code": f"{case.family}-B",
                "ANGLE": row["ANGLE"],
                "E1": row["E1"],
                "entrance_rms_x": start["rms_x"],
                "entrance_rms_y": start["rms_y"],
                "entrance_rms_z": start["rms_z"],
                "entrance_rms_px": start["rms_px"],
                "entrance_rms_py": start["rms_py"],
                "entrance_rms_ps": start["rms_ps"],
                "rms_x": row["bmad_exit_rms_x"],
                "rms_y": row["bmad_exit_rms_y"],
                "rms_z": row["bmad_exit_rms_z"],
                "rms_px": row["bmad_exit_rms_px"],
                "rms_py": row["bmad_exit_rms_py"],
                "rms_ps": row["bmad_exit_rms_ps"],
            }
        )

    summary = pd.DataFrame(rows)
    summary = summary[
        [
            "family_code",
            "ANGLE",
            "E1",
            "rms_x",
            "rms_y",
            "rms_z",
            "rms_px",
            "rms_py",
            "rms_ps",
        ]
    ]
    numeric_cols = ["rms_x", "rms_y", "rms_z", "rms_px", "rms_py", "rms_ps"]
    summary[numeric_cols] = summary[numeric_cols].round(5)
    summary.to_csv(SUMMARY_CSV, index=False, float_format="%.5f")
    write_text(SUMMARY_MD, rows_to_markdown(summary))


if __name__ == "__main__":
    main()
