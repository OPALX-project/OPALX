#!/usr/bin/env python3
from __future__ import annotations

import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

sys.path.insert(0, "/private/tmp/opalx-py")

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
BMAD_EXE = Path(
    "/Users/adelmann/git/bmad-opalx/install/bmad-ecosystem/production/bin/beam_track_example"
)
OPALX_EXE = Path("/Users/adelmann/git/opalx/build/src/opalx")

# Electron rest energy in eV.
ELECTRON_MASS_EV = 510_998.95
BEGINNING_E_TOT_EV = 1.00051099895e9
REFERENCE_BETA_GAMMA = math.sqrt((BEGINNING_E_TOT_EV / ELECTRON_MASS_EV) ** 2 - 1.0)
SPEED_OF_LIGHT = 299_792_458.0

N_PARTICLE = 1000
SIGMA = 1.0e-3
RNG_SEED = 12345
OPALX_DT = 1.0e-11
MONITOR_LENGTH = 0.0
INPUT_CANONICAL_CSV = ROOT / "sbend_exit_monitor_input_canonical.csv"
SUMMARY_CSV = ROOT / "sbend_exit_monitor_summary.csv"
SUMMARY_MD = ROOT / "sbend_exit_monitor_summary.md"
DIAGNOSTICS_CSV = ROOT / "sbend_exit_monitor_diagnostics.csv"
DIAGNOSTICS_MD = ROOT / "sbend_exit_monitor_diagnostics.md"


@dataclass(frozen=True)
class E1Config:
    tag: str
    label: str
    bmad_expr: str
    opalx_expr: str


@dataclass(frozen=True)
class CaseConfig:
    tag: str
    angle_label: str
    bmad_angle_expr: str
    opalx_angle_expr: str
    angle_rad: float
    e1: E1Config

    @property
    def exit_x(self) -> float:
        curvature = self.angle_rad
        return (math.cos(self.angle_rad) - 1.0) / curvature

    @property
    def exit_z(self) -> float:
        curvature = self.angle_rad
        return math.sin(self.angle_rad) / curvature

    @property
    def exit_theta_opalx_expr(self) -> str:
        return {
            "sb45p": "-PI / 4",
            "sb45m": "PI / 4",
            "sb90p": "-PI / 2",
            "sb90m": "PI / 2",
        }[self.tag]

    @property
    def exit_theta_rad(self) -> float:
        return -self.angle_rad

    @property
    def case_prefix(self) -> str:
        return f"sbend_exit_monitor_{self.tag}_{self.e1.tag}"


E1_CASES = [
    E1Config("e10_0", "0", "0.0", "0.0"),
    E1Config("e10_p", "10deg", "10 * pi / 180", "10 * PI / 180"),
]


CASES = [
    CaseConfig(tag, angle_label, bmad_angle_expr, opalx_angle_expr, angle_rad, e1)
    for e1 in E1_CASES
    for tag, angle_label, bmad_angle_expr, opalx_angle_expr, angle_rad in [
        ("sb45p", "+pi/4", "pi / 4", "PI / 4", math.pi / 4.0),
        ("sb45m", "-pi/4", "-pi / 4", "-PI / 4", -math.pi / 4.0),
        ("sb90p", "+pi/2", "pi / 2", "PI / 2", math.pi / 2.0),
        ("sb90m", "-pi/2", "-pi / 2", "-PI / 2", -math.pi / 2.0),
    ]
]


def case_paths(case: CaseConfig) -> dict[str, Path]:
    monitor_base = f"{case.case_prefix}_monitor"
    return {
        "bmad_input_beam": ROOT / f"{case.case_prefix}_bmad_input.dat",
        "bmad_lattice": ROOT / f"{case.case_prefix}.bmad",
        "bmad_track_input": ROOT / f"{case.case_prefix}_beam_track.in",
        "bmad_track_log": ROOT / f"{case.case_prefix}_beam_track.log",
        "bmad_output_beam": ROOT / f"{case.case_prefix}_bmad_output.dat",
        "opalx_input_distribution": ROOT / f"{case.case_prefix}_opalx_distribution.txt",
        "opalx_deck": ROOT / f"{case.case_prefix}_opalx.in",
        "opalx_log": ROOT / f"{case.case_prefix}_opalx.log",
        "opalx_stat": ROOT / f"{case.case_prefix}_opalx.stat",
        "opalx_lbal": ROOT / f"{case.case_prefix}_opalx.lbal",
        "opalx_h5": ROOT / f"{case.case_prefix}_opalx.h5",
        "opalx_monitor_base": ROOT / monitor_base,
        "opalx_monitor_h5": ROOT / f"{monitor_base}.h5",
        "opalx_monitor_stat": ROOT / f"{case.case_prefix}_opalx_Monitors.stat",
        "element_positions_txt": DATA_DIR / f"{case.case_prefix}_opalx_ElementPositions.txt",
        "element_positions_py": DATA_DIR / f"{case.case_prefix}_opalx_ElementPositions.py",
        "element_positions_sdds": DATA_DIR / f"{case.case_prefix}_opalx_ElementPositions.sdds",
    }


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

    The input coordinates are

    .. math::

       (x, p_x, y, p_y, \zeta, \delta),

    where each component is sampled independently from a normal distribution
    with mean zero and sigma :math:`10^{-3}`.
    """

    rng = np.random.default_rng(RNG_SEED)
    data = rng.normal(loc=0.0, scale=SIGMA, size=(N_PARTICLE, 6))
    return pd.DataFrame(data, columns=["x", "px_can", "y", "py_can", "zeta", "delta"])


def canonical_to_opalx(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert BMAD canonical coordinates to OPALX local beta-gamma momentum components.

    OPALX stores the local momentum components as

    .. math::

       (\beta\gamma)_x,\; (\beta\gamma)_y,\; (\beta\gamma)_s,

    while BMAD uses the canonical variables

    .. math::

       p_x = P_x / P_0,\qquad
       p_y = P_y / P_0,\qquad
       \delta = (P - P_0) / P_0.

    For the reference beta-gamma :math:`(\beta\gamma)_0`, the conversion is

    .. math::

       (\beta\gamma)_x &= (\beta\gamma)_0 p_x, \\
       (\beta\gamma)_y &= (\beta\gamma)_0 p_y, \\
       (\beta\gamma)_s &= (\beta\gamma)_0
           \sqrt{(1 + \delta)^2 - p_x^2 - p_y^2}.
    """

    px_can = df["px_can"].to_numpy()
    py_can = df["py_can"].to_numpy()
    delta = df["delta"].to_numpy()
    ps = np.sqrt((1.0 + delta) ** 2 - px_can**2 - py_can**2)
    return pd.DataFrame(
        {
            "x": df["x"].to_numpy(),
            "px": REFERENCE_BETA_GAMMA * px_can,
            "y": df["y"].to_numpy(),
            "py": REFERENCE_BETA_GAMMA * py_can,
            "z": df["zeta"].to_numpy(),
            "pz": REFERENCE_BETA_GAMMA * ps,
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
            f"{row.x:.16e} {row.px_can:.16e} {row.y:.16e} {row.py_can:.16e} "
            f"{row.zeta:.16e} {row.delta:.16e} Alive"
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


def write_bmad_inputs(case: CaseConfig, paths: dict[str, Path]) -> None:
    write_text(
        paths["bmad_lattice"],
        f"""parameter[particle] = electron
parameter[geometry] = open

beginning[e_tot] = 1.00051099895e9
beginning[beta_a] = 1.0
beginning[beta_b] = 1.0
beginning[alpha_a] = 0.0
beginning[alpha_b] = 0.0

b1: sbend, l = 1.0, angle = {case.bmad_angle_expr}, e1 = {case.e1.bmad_expr}, e2 = 0.0, &
    fringe_type = none, tracking_method = runge_kutta, num_steps = 800

line1: line = (b1)
use, line1

no_digested
end_file
""",
    )
    write_text(
        paths["bmad_track_input"],
        f"""&beam_track_params
    lat_filename = '{paths["bmad_lattice"].name}'
    outfile_name = '{paths["bmad_output_beam"].name}'
    beam_init%position_file = '{paths["bmad_input_beam"].name}'
    beam_init%n_particle = 0
    beam_init%n_bunch = 0
    ran_seed = {RNG_SEED}
/
""",
    )


def write_opalx_deck(case: CaseConfig, paths: dict[str, Path]) -> None:
    write_text(
        paths["opalx_deck"],
        f"""OPTION, ENABLEHDF5 = TRUE;
OPTION, PSDUMPFREQ = 50000;
OPTION, STATDUMPFREQ = 1;
OPTION, BOUNDPDESTROY = 10;
OPTION, VERSION = 10900;

Title, string = "sbend_exit_monitor SBEND ANGLE={case.angle_label} E1={case.e1.label} PSI=0";

B1: SBEND,
    L = 1.0,
    ANGLE = {case.opalx_angle_expr},
    E1 = {case.e1.opalx_expr},
    PSI = 0.0,
    ELEMEDGE = 0.0;

MEXIT: MONITOR,
    L = {MONITOR_LENGTH:.16e},
    X = {case.exit_x:.16e},
    Y = 0.0,
    Z = {case.exit_z:.16e},
    THETA = {case.exit_theta_opalx_expr},
    PHI = 0.0,
    PSI = 0.0,
    TYPE = TEMPORAL,
    OUTFN = "{paths["opalx_monitor_base"].name}";

Line1: Line = (B1, MEXIT);

Dist0: DISTRIBUTION,
    TYPE = FROMFILE,
    FNAME = "{paths["opalx_input_distribution"].name}",
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
""",
    )


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


def cleanup_previous_outputs() -> None:
    patterns = [
        "sbend_exit_monitor_*.bmad",
        "sbend_exit_monitor_*_beam_track.in",
        "sbend_exit_monitor_*_beam_track.log",
        "sbend_exit_monitor_*_bmad_input.dat",
        "sbend_exit_monitor_*_bmad_output.dat",
        "sbend_exit_monitor_*_opalx.in",
        "sbend_exit_monitor_*_opalx.log",
        "sbend_exit_monitor_*_opalx.stat",
        "sbend_exit_monitor_*_opalx.lbal",
        "sbend_exit_monitor_*_opalx.h5",
        "sbend_exit_monitor_*_opalx_distribution.txt",
        "sbend_exit_monitor_*_monitor.h5",
        "sbend_exit_monitor_*_Monitors.stat",
        "sbend_exit_monitor_summary.csv",
        "sbend_exit_monitor_summary.md",
        "sbend_exit_monitor_diagnostics.csv",
        "sbend_exit_monitor_diagnostics.md",
    ]
    data_patterns = [
        "sbend_exit_monitor_*_opalx_ElementPositions.txt",
        "sbend_exit_monitor_*_opalx_ElementPositions.py",
        "sbend_exit_monitor_*_opalx_ElementPositions.sdds",
    ]
    for pattern in patterns:
        for path in ROOT.glob(pattern):
            if path.is_file():
                path.unlink()
    for pattern in data_patterns:
        for path in DATA_DIR.glob(pattern):
            if path.is_file():
                path.unlink()


def parse_bmad_beam_file(path: Path) -> pd.DataFrame:
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
    return df


def scalar_attr(value: object) -> float:
    arr = np.asarray(value, dtype=float).reshape(-1)
    return float(arr[0])


def read_opalx_stat_rows(path: Path) -> pd.DataFrame:
    columns: list[str] = []
    rows: list[list[float]] = []
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
        else:
            stripped = line.strip()
            if columns and stripped and not stripped.startswith(("SDDS1", "&", "!")):
                parts = stripped.split()
                if len(parts) == len(columns):
                    try:
                        rows.append([float(value) for value in parts])
                    except ValueError:
                        pass
        i += 1

    if not columns or not rows:
        raise RuntimeError(f"could not parse OPALX stat rows: {path}")
    return pd.DataFrame(rows, columns=columns)


def exit_plane_signed_distance(case: CaseConfig, stat_rows: pd.DataFrame) -> dict[str, float]:
    """
    Compute the signed reference-particle distance to the nominal SBEND exit plane.

    For explicitly placed straight elements the local longitudinal axis of the
    body chart follows

    .. math::

       \hat{\mathbf{s}} = (\sin\theta,\; 0,\; \cos\theta),

    in the laboratory :math:`x-z` plane. The signed plane coordinate is then

    .. math::

       s_\mathrm{plane} =
       (\mathbf{r}_\mathrm{ref} - \mathbf{r}_0)\cdot\hat{\mathbf{s}},

    with :math:`\mathbf{r}_0` the nominal exit-plane origin.
    """

    tangent_x = math.sin(case.exit_theta_rad)
    tangent_z = math.cos(case.exit_theta_rad)
    dx = stat_rows["ref_x"] - case.exit_x
    dz = stat_rows["ref_z"] - case.exit_z
    signed_distance = dx * tangent_x + dz * tangent_z
    final = stat_rows.iloc[-1]
    return {
        "ANGLE": case.angle_label,
        "E1": case.e1.label,
        "signed_distance_min": float(signed_distance.min()),
        "signed_distance_max": float(signed_distance.max()),
        "signed_distance_final": float(signed_distance.iloc[-1]),
        "crosses_nominal_exit_plane": float(
            (signed_distance.min() <= 0.0) and (signed_distance.max() >= 0.0)
        ),
        "final_ref_x": float(final["ref_x"]),
        "final_ref_z": float(final["ref_z"]),
        "final_t_ns": float(final["t"]),
        "final_s_path_m": float(final["s"]),
    }


def rms_from_canonical_dataframe(df: pd.DataFrame) -> dict[str, float]:
    return {
        "rms_x": float(df["x"].std(ddof=0)),
        "rms_y": float(df["y"].std(ddof=0)),
        "rms_zeta": float(df["zeta"].std(ddof=0)),
        "rms_px_can": float(df["px_can"].std(ddof=0)),
        "rms_py_can": float(df["py_can"].std(ddof=0)),
        "rms_delta": float(df["delta"].std(ddof=0)),
    }


def rms_from_bmad_output(df: pd.DataFrame) -> dict[str, float]:
    return {
        "rms_x": float(df["x"].std(ddof=0)),
        "rms_y": float(df["y"].std(ddof=0)),
        "rms_zeta": float(df["z"].std(ddof=0)),
        "rms_px_can": float(df["px"].std(ddof=0)),
        "rms_py_can": float(df["py"].std(ddof=0)),
        "rms_delta": float(df["pz"].std(ddof=0)),
    }


def rms_from_opalx_monitor(path: Path) -> tuple[dict[str, float], dict[str, float]]:
    with h5py.File(path, "r") as h5file:
        step_names = sorted(h5file.keys())
        if not step_names:
            raise RuntimeError(f"monitor file has no steps: {path}")

        step = h5file[step_names[-1]]
        x = step["x"][:]
        y = step["y"][:]
        px_bg = step["px"][:]
        py_bg = step["py"][:]
        pz_bg = step["pz"][:]
        time = step["time"][:]

        ref_time = scalar_attr(step.attrs["TIME"])
        ref_p = np.asarray(step.attrs["RefPartP"], dtype=float).reshape(3)
        bg0 = float(np.linalg.norm(ref_p))
        beta0 = bg0 / math.sqrt(1.0 + bg0 * bg0)

        zeta = -beta0 * SPEED_OF_LIGHT * (time - ref_time)
        px_can = px_bg / bg0
        py_can = py_bg / bg0
        delta = np.sqrt(px_bg**2 + py_bg**2 + pz_bg**2) / bg0 - 1.0

        summary = {
            "rms_x": float(np.std(x, ddof=0)),
            "rms_y": float(np.std(y, ddof=0)),
            "rms_zeta": float(np.std(zeta, ddof=0)),
            "rms_px_can": float(np.std(px_can, ddof=0)),
            "rms_py_can": float(np.std(py_can, ddof=0)),
            "rms_delta": float(np.std(delta, ddof=0)),
        }
        meta = {
            "monitor_s": scalar_attr(step.attrs["SPOS"]),
            "ref_time_s": ref_time,
            "ref_bg0": bg0,
            "particle_count": float(x.size),
        }

    return summary, meta


def build_markdown(summary: pd.DataFrame, diagnostics: pd.DataFrame) -> str:
    explanation = f"""# SBEND Exit Monitor Comparison

This study compares the same {N_PARTICLE}-particle Gaussian bunch at the entrance and at the SBEND exit plane for `ANGLE = ±pi/4, ±pi/2` and `E1 = 0, 10deg`.

Definitions:
- `rms_x`, `rms_y`: RMS transverse coordinates in the local reference frame at the comparison plane.
- `rms_zeta`: RMS longitudinal phase coordinate in metres, reconstructed as `zeta = -beta0 * c * (t - t_ref)`.
- `rms_px_can`, `rms_py_can`: RMS canonical transverse momenta `P_x / P_0` and `P_y / P_0`.
- `rms_delta`: RMS relative momentum error `(P - P_0) / P_0`.

For BMAD the exit quantities come directly from the canonical beam output file. For OPALX the exit quantities are reconstructed from the monitor H5 dump, which stores local particle coordinates, local beta-gamma momenta, and particle crossing times at the monitor plane.

## RMS Comparison
"""
    return explanation + rows_to_markdown(summary) + "\n## Exit-Plane Diagnostics\n" + rows_to_markdown(diagnostics)


def run_case(case: CaseConfig, canonical_df: pd.DataFrame) -> tuple[list[dict[str, object]], dict[str, object]]:
    paths = case_paths(case)
    write_bmad_beam_file(canonical_df, paths["bmad_input_beam"])
    write_opalx_distribution(canonical_to_opalx(canonical_df), paths["opalx_input_distribution"])
    write_bmad_inputs(case, paths)
    write_opalx_deck(case, paths)

    run_logged([str(BMAD_EXE), paths["bmad_track_input"].name], paths["bmad_track_log"])
    run_logged([str(OPALX_EXE), "--info", "2", paths["opalx_deck"].name], paths["opalx_log"])

    bmad_exit = rms_from_bmad_output(parse_bmad_beam_file(paths["bmad_output_beam"]))
    stat_rows = read_opalx_stat_rows(paths["opalx_stat"])
    plane_meta = exit_plane_signed_distance(case, stat_rows)

    monitor_meta: dict[str, float] | None = None
    if paths["opalx_monitor_h5"].exists():
        opalx_exit, monitor_meta = rms_from_opalx_monitor(paths["opalx_monitor_h5"])
        opalx_status = "OK"
    else:
        opalx_exit = {key: float("nan") for key in bmad_exit.keys()}
        opalx_status = "NO_CROSSING_NOMINAL_EXIT"

    if monitor_meta is None:
        monitor_meta = {
            "monitor_s": float("nan"),
            "ref_time_s": float("nan"),
            "ref_bg0": float("nan"),
            "particle_count": 0.0,
        }

    rows = [
        {
            "row": "EXIT",
            "code": "BMAD",
            "ANGLE": case.angle_label,
            "E1": case.e1.label,
            "status": "OK",
            **bmad_exit,
        },
        {
            "row": "EXIT",
            "code": "OPALX-MON",
            "ANGLE": case.angle_label,
            "E1": case.e1.label,
            "status": opalx_status,
            **opalx_exit,
        },
    ]

    diagnostic = {
        **plane_meta,
        "status": opalx_status,
        "monitor_s": monitor_meta["monitor_s"],
        "t_ref_ns": monitor_meta["ref_time_s"] * 1.0e9,
        "ref_bg0": monitor_meta["ref_bg0"],
        "n_particle": int(monitor_meta["particle_count"]),
    }
    return rows, diagnostic


def main() -> None:
    canonical_df = generate_canonical_bunch()
    canonical_df.to_csv(INPUT_CANONICAL_CSV, index=False)

    cleanup_previous_outputs()

    summary_rows = [
        {
            "row": "ENTRANCE",
            "code": "COMMON",
            "ANGLE": "ALL",
            "E1": "ALL",
            "status": "OK",
            **rms_from_canonical_dataframe(canonical_df),
        }
    ]
    diagnostics_rows: list[dict[str, object]] = []

    for case in CASES:
        case_rows, case_diag = run_case(case, canonical_df)
        summary_rows.extend(case_rows)
        diagnostics_rows.append(case_diag)

    summary = pd.DataFrame(summary_rows)
    diagnostics = pd.DataFrame(diagnostics_rows)
    summary.to_csv(SUMMARY_CSV, index=False)
    diagnostics.to_csv(DIAGNOSTICS_CSV, index=False)
    write_text(SUMMARY_MD, build_markdown(summary, diagnostics))
    write_text(DIAGNOSTICS_MD, rows_to_markdown(diagnostics))


if __name__ == "__main__":
    main()
