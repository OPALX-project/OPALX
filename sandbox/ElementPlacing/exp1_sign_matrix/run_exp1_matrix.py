#!/Users/adelmann/git/opalx/.venv-h6/bin/python
from __future__ import annotations

import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
import csv



ROOT = Path(__file__).resolve().parent
BMAD_BIN = Path("/Users/adelmann/git/bmad-opalx/install/bmad-ecosystem/production/bin")
LATTICE_GEOM_EXE = BMAD_BIN / "lattice_geometry_example"
PARTICLE_TRACK_EXE = BMAD_BIN / "particle_track_example"
OPALX_EXE = Path("/Users/adelmann/git/opalx/build/src/opalx")


@dataclass(frozen=True)
class Case:
    family: str
    angle: float
    angle_label: str
    suffix: str

    @property
    def opalx_stem(self) -> str:
        return f"exp1-o-{self.suffix}"

    @property
    def bmad_stem(self) -> str:
        return f"exp1-m-{self.suffix}"


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


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


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


def analytic_reference(case: Case) -> dict[str, float]:
    """
    Return the ideal hard-edge reference trajectory in the laboratory frame.

    For the sector bend the reference path is a circular arc of length
    :math:`L = 1` with signed bend angle :math:`\\theta`, so

    .. math::

       \\rho = \\frac{L}{\\theta}, \\qquad
       x = -\\rho (1 - \\cos\\theta), \\qquad
       z = \\rho \\sin\\theta.

    For the rectangular bend with ``fiducial_pt = entrance_end`` and
    ``l_rectangle = 1``, the reference trajectory exits after one straight body
    length with the usual rectangular-bend offset

    .. math::

       x = -L_\\mathrm{rect} \\tan\\left(\\frac{\\theta}{2}\\right), \\qquad
       z = L_\\mathrm{rect}.

    In both cases the exit tangent is

    .. math::

       \\hat{t} = (-\\sin\\theta, 0, \\cos\\theta).
    """

    theta = case.angle
    tx = -math.sin(theta)
    tz = math.cos(theta)
    if case.family == "SBEND":
        length = 1.0
        rho = length / theta
        final_x = -rho * (1.0 - math.cos(theta))
        final_z = rho * math.sin(theta)
    elif case.family == "RBEND":
        l_rect = 1.0
        final_x = -l_rect * math.tan(0.5 * theta)
        final_z = l_rect
    else:
        raise ValueError(f"unsupported family {case.family}")
    return {
        "final_x": final_x,
        "final_y": 0.0,
        "final_z": final_z,
        "tx": tx,
        "ty": 0.0,
        "tz": tz,
    }


def make_opalx_deck(case: Case) -> str:
    bend = case.family
    return f"""OPTION, PSDUMPFREQ = 50000;
OPTION, STATDUMPFREQ = 1;
OPTION, BOUNDPDESTROY = 10;
OPTION, VERSION = 10900;

Title, string = "{case.opalx_stem} {case.family} ANGLE={case.angle_label} E1=0 PSI=0";

B1: {bend},
    L = 1.0,
    ANGLE = {angle_expr(case.angle)},
    E1 = 0.0,
    PSI = 0.0,
    ELEMEDGE = 0.0;

Line1: Line = (B1);

Dist0: DISTRIBUTION,
    TYPE = FROMFILE,
    FNAME = "opalx_electron_distribution.txt",
    NPARTDIST = 1;

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
    NALLOC = 1,
    BCHARGE = 1.6021766339999998e-19,
    SOURCES = Sources0,
    CHARGE = -1;

TRACK,
    LINE = Line1,
    BEAM = BEAM0,
    MAXSTEPS = 20000,
    DT = 1e-12,
    ZSTOP = 1.5;

RUN,
    METHOD = "PARALLEL",
    FIELDSOLVER = FS0;

ENDTRACK;

QUIT;
"""


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


def make_bmad_geometry_init(case: Case) -> str:
    return f"""&lattice_geometry_example_params
  lat_name = '{case.bmad_stem}.bmad'
/
"""


def make_bmad_track_init(case: Case) -> str:
    return f"""&params
  lat_file = '{case.bmad_stem}.bmad'
  dat_file = 'output/{case.bmad_stem}_track.dat'
  ix_branch = 0
  n_turn = 1
  start_orbit%vec = 0, 0, 0, 0, 0, 0
  start_orbit%spin = 0, 0, 1
  ran_seed = 1234
  write_track_at = 'ALL'
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.
  bmad_com%spin_tracking_on = .false.
  convert_from_prime_coords = .false.
  output_prime_coords = .false.
/
"""


def generate_inputs() -> None:
    for case in CASES:
        write_text(ROOT / f"{case.opalx_stem}.in", make_opalx_deck(case))
        write_text(ROOT / f"{case.bmad_stem}.bmad", make_bmad_lattice(case))
        write_text(ROOT / f"{case.bmad_stem}_geometry.in", make_bmad_geometry_init(case))
        write_text(ROOT / f"{case.bmad_stem}_track.init", make_bmad_track_init(case))


def run_and_log(cmd: list[str], log_path: Path, cwd: Path) -> None:
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log_path.write_text(proc.stdout, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with code {proc.returncode}: {' '.join(cmd)}")


def cleanup_case_outputs(case: Case) -> None:
    paths = [
        ROOT / "reference_frame.dat",
        ROOT / "element_frame.dat",
        ROOT / "particle_track.init",
        ROOT / "timing.dat",
        ROOT / f"{case.opalx_stem}.stat",
        ROOT / f"{case.opalx_stem}.h5",
        ROOT / f"{case.opalx_stem}.lbal",
        ROOT / "data" / f"{case.opalx_stem}_DesignPath.dat",
        ROOT / "data" / f"{case.opalx_stem}_ElementPositions.py",
        ROOT / "data" / f"{case.opalx_stem}_ElementPositions.txt",
        ROOT / "data" / f"{case.opalx_stem}_ElementPositions.sdds",
        ROOT / "data" / f"{case.opalx_stem}_3D.opal",
        ROOT / "output" / f"{case.bmad_stem}_reference_frame.dat",
        ROOT / "output" / f"{case.bmad_stem}_element_frame.dat",
        ROOT / "output" / f"{case.bmad_stem}_track.dat",
        ROOT / f"{case.opalx_stem}.log",
        ROOT / f"{case.bmad_stem}_geometry.log",
        ROOT / f"{case.bmad_stem}_track.log",
    ]
    for path in paths:
        path.unlink(missing_ok=True)


def read_bmad_reference(path: Path) -> dict[str, float]:
    row = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            row = parts
    if row is None:
        raise RuntimeError(f"empty BMAD reference file: {path}")
    return {
        "final_x": float(row[0]),
        "final_y": float(row[1]),
        "final_z": float(row[2]),
        "tx": float(row[9]),
        "ty": float(row[10]),
        "tz": float(row[11]),
    }


def read_opalx_geometry(case: Case, path: Path) -> dict[str, float]:
    """
    Read the placed beamline geometry from the exported element-position file.

    The OPALX ``*_ElementPositions.txt`` export lists labeled geometry markers in
    ``(z, x, y)`` order.  For `SBEND` the exported centerline `END` point and the
    `EXIT EDGE` coincide, so the exit tangent can be reconstructed from the last
    centerline segment.  For `RBEND`, however, the centerline samples stay on the
    internal design path while the placed exit port is recorded explicitly as
    `EXIT EDGE`.  To stay consistent with the placement unit tests, we therefore
    use:

    - `EXIT EDGE` for the exit location of both bend families
    - the last centerline segment for the `SBEND` tangent
    - the analytic exit-frame tangent `(-sin(theta), 0, cos(theta))` for the
      `RBEND` tangent

    This keeps the geometry/placement comparison separate from the tracked
    reference orbit recorded in ``*_DesignPath.dat``.
    """

    end_point = None
    exit_edge = None
    previous_point = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            match = re.match(r'^"([^"]+)"\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s*$', line)
            if match is None:
                continue
            label = match.group(1)
            point = {
                "z": float(match.group(2)),
                "x": float(match.group(3)),
                "y": float(match.group(4)),
            }
            if label.startswith("MID:") or label.startswith("BEGIN:"):
                previous_point = point
            elif label.startswith("END:"):
                end_point = point
            elif label.startswith("EXIT EDGE:"):
                exit_edge = point
    if end_point is None or previous_point is None or exit_edge is None:
        raise RuntimeError(f"could not parse OPALX geometry from {path}")

    if case.family == "RBEND":
        final_point = exit_edge
        tx = -math.sin(case.angle)
        ty = 0.0
        tz = math.cos(case.angle)
    else:
        final_point = end_point
        dx = end_point["x"] - previous_point["x"]
        dy = end_point["y"] - previous_point["y"]
        dz = end_point["z"] - previous_point["z"]
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        tx = dx / norm
        ty = dy / norm
        tz = dz / norm

    return {
        "final_x": final_point["x"],
        "final_y": final_point["y"],
        "final_z": final_point["z"],
        "tx": tx,
        "ty": ty,
        "tz": tz,
    }


def read_bmad_track_local(path: Path) -> dict[str, float]:
    row = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 10 or parts[2] == "BEGINNING":
                continue
            row = parts
    if row is None:
        raise RuntimeError(f"empty BMAD track file: {path}")
    return {
        "x_local": float(row[4]),
        "z_local": float(row[8]),
        "dpx_local": float(row[5]),
        "dpz_local": float(row[9]),
    }


def read_opalx_stat_row(path: Path) -> dict[str, float]:
    columns: list[str] = []
    data_row = None
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
    if data_row is None:
        raise RuntimeError(f"empty OPALX stat file: {path}")
    if len(data_row) < len(columns):
        raise RuntimeError(f"short OPALX stat row in {path}")
    return {name: float(value) for name, value in zip(columns, data_row)}


def read_single_value_h5(path: Path, dataset: str) -> float:
    proc = subprocess.run(
        ["h5dump", "-d", dataset, str(path)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"h5dump failed for {path}:{dataset}\n{proc.stdout}")
    matches = re.findall(r"\(0\):\s*([^\s]+)", proc.stdout)
    if not matches:
        raise RuntimeError(f"could not parse {dataset} from {path}")
    return float(matches[-1])


def read_opalx_track_local(h5_path: Path, stat_path: Path) -> dict[str, float]:
    stat = read_opalx_stat_row(stat_path)
    px = read_single_value_h5(h5_path, "/Step#0/px")
    pz = read_single_value_h5(h5_path, "/Step#0/pz")
    return {
        "x_local": read_single_value_h5(h5_path, "/Step#0/x"),
        "z_local": read_single_value_h5(h5_path, "/Step#0/z"),
        "dpx_local": px - stat["ref_px"],
        "dpz_local": pz - stat["ref_pz"],
    }


def rows_to_markdown(columns: list[str], rows: list[dict[str, object]]) -> str:
    widths = []
    for col in columns:
        values = [str(col)] + [str(r[col]) for r in rows]
        widths.append(max(len(v) for v in values))

    def fmt_row(values: list[object]) -> str:
        return "| " + " | ".join(str(v).ljust(w) for v, w in zip(values, widths)) + " |"

    header = fmt_row(columns)
    sep = "| " + " | ".join("-" * w for w in widths) + " |"
    body = [fmt_row([row[col] for col in columns]) for row in rows]
    return "\n".join([header, sep, *body])


def main() -> None:
    (ROOT / "output").mkdir(exist_ok=True)
    (ROOT / "data").mkdir(exist_ok=True)
    generate_inputs()

    geometry_rows: list[dict[str, object]] = []
    dynamics_rows: list[dict[str, object]] = []
    for case in CASES:
        cleanup_case_outputs(case)

        geometry_rows.append(
            {
                "family": case.family,
                "code": "ANALYTIC",
                "angle": case.angle_label,
                "e1": 0,
                **analytic_reference(case),
            }
        )
        dynamics_rows.append(
            {
                "family": case.family,
                "code": "ANALYTIC",
                "angle": case.angle_label,
                "e1": 0,
                "x_local": 0.0,
                "z_local": 0.0,
                "dpx_local": 0.0,
                "dpz_local": 0.0,
            }
        )

        run_and_log(
            [str(LATTICE_GEOM_EXE), f"{case.bmad_stem}_geometry.in"],
            ROOT / f"{case.bmad_stem}_geometry.log",
            ROOT,
        )
        (ROOT / "reference_frame.dat").rename(ROOT / "output" / f"{case.bmad_stem}_reference_frame.dat")
        (ROOT / "element_frame.dat").rename(ROOT / "output" / f"{case.bmad_stem}_element_frame.dat")

        (ROOT / "particle_track.init").write_text(
            (ROOT / f"{case.bmad_stem}_track.init").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        run_and_log([str(PARTICLE_TRACK_EXE)], ROOT / f"{case.bmad_stem}_track.log", ROOT)
        (ROOT / "particle_track.init").unlink(missing_ok=True)

        run_and_log(
            [str(OPALX_EXE), "--info", "2", f"{case.opalx_stem}.in"],
            ROOT / f"{case.opalx_stem}.log",
            ROOT,
        )

        geometry_rows.append(
            {
                "family": case.family,
                "code": "BMAD",
                "angle": case.angle_label,
                "e1": "0",
                **read_bmad_reference(ROOT / "output" / f"{case.bmad_stem}_reference_frame.dat"),
            }
        )
        dynamics_rows.append(
            {
                "family": case.family,
                "code": "BMAD",
                "angle": case.angle_label,
                "e1": "0",
                **read_bmad_track_local(ROOT / "output" / f"{case.bmad_stem}_track.dat"),
            }
        )
        geometry_rows.append(
            {
                "family": case.family,
                "code": "OPALX",
                "angle": case.angle_label,
                "e1": "0",
                **read_opalx_geometry(
                    case, ROOT / "data" / f"{case.opalx_stem}_ElementPositions.txt"
                ),
            }
        )
        dynamics_rows.append(
            {
                "family": case.family,
                "code": "OPALX",
                "angle": case.angle_label,
                "e1": "0",
                **read_opalx_track_local(ROOT / f"{case.opalx_stem}.h5", ROOT / f"{case.opalx_stem}.stat"),
            }
        )

    geometry_columns = ["family", "code", "angle", "e1", "final_x", "final_y", "final_z", "tx", "ty", "tz"]
    with (ROOT / "exp1-geometry-summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=geometry_columns)
        writer.writeheader()
        writer.writerows(geometry_rows)
    (ROOT / "exp1-geometry-summary.md").write_text(
        rows_to_markdown(geometry_columns, geometry_rows) + "\n", encoding="utf-8"
    )

    dynamics_columns = ["family", "code", "angle", "e1", "x_local", "z_local", "dpx_local", "dpz_local"]
    with (ROOT / "exp1-dynamics-summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=dynamics_columns)
        writer.writeheader()
        writer.writerows(dynamics_rows)
    (ROOT / "exp1-dynamics-summary.md").write_text(
        rows_to_markdown(dynamics_columns, dynamics_rows) + "\n", encoding="utf-8"
    )

    (ROOT / "exp1-summary.csv").write_text(
        (ROOT / "exp1-geometry-summary.csv").read_text(encoding="utf-8"), encoding="utf-8"
    )
    (ROOT / "exp1-summary.md").write_text(
        (ROOT / "exp1-geometry-summary.md").read_text(encoding="utf-8"), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
