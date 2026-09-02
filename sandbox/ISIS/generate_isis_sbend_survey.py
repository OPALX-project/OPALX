#!/usr/bin/env python3
"""Generate an explicit-placement SBEND survey of the ISIS ring.

The source deck uses ELEMEDGE and curved MULTIPOLET elements.  For this
geometry-only model, all zero-angle components become DRIFT envelopes and the
ten curved dipoles become native SBENDs.  Every element receives an explicit
global X/Y/Z/THETA/PHI/PSI placement, following the convention exercised by
sandbox/test-sbend-ring-1.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "isisring_with_elemedge.in"
OUTPUT = ROOT / "isis_sbend_survey.in"


@dataclass(frozen=True)
class SourceElement:
    name: str
    kind: str
    length: float
    elemedge: float
    angle: float


def source_elements() -> list[SourceElement]:
    text = SOURCE.read_text()
    definitions: dict[str, SourceElement] = {}
    for match in re.finditer(
        r"^\s*([A-Za-z][\w]*)\s*:\s*(\w+)\s*,\s*L\s*=\s*"
        r"([0-9.eE+-]+)\s*,\s*ELEMEDGE\s*=\s*([0-9.eE+-]+)(.*?);",
        text,
        re.MULTILINE,
    ):
        name, kind, length, elemedge, tail = match.groups()
        angle_match = re.search(r"\bANGLE\s*=\s*([0-9.eE+-]+)", tail)
        definitions[name.lower()] = SourceElement(
            name=name,
            kind=kind.upper(),
            length=float(length),
            elemedge=float(elemedge),
            angle=float(angle_match.group(1)) if angle_match else 0.0,
        )

    line_match = re.search(
        r"^line0\s*:\s*LINE\s*=\s*\((.*?)\);", text, re.MULTILINE | re.DOTALL
    )
    if line_match is None:
        raise RuntimeError("could not find line0 in source deck")
    names = [name.strip() for name in line_match.group(1).split(",")]
    elements = [definitions[name.lower()] for name in names]

    s = 0.0
    for element in elements:
        if not math.isclose(element.elemedge, s, rel_tol=0.0, abs_tol=1.0e-7):
            raise RuntimeError(
                f"{element.name}: ELEMEDGE={element.elemedge:.12g}, expected {s:.12g}"
            )
        s += element.length

    bends = [element for element in elements if element.angle != 0.0]
    if len(elements) != 192 or len(bends) != 10:
        raise RuntimeError(
            f"unexpected ISIS topology: {len(elements)} elements, {len(bends)} bends"
        )
    return elements


def number(value: float) -> str:
    if abs(value) < 5.0e-14:
        value = 0.0
    return f"{value:.12g}"


def main() -> None:
    elements = source_elements()
    bend_angle = 2.0 * math.pi / 10.0
    source_angles = {element.angle for element in elements if element.angle != 0.0}
    if len(source_angles) != 1 or not math.isclose(
        source_angles.pop(), bend_angle, rel_tol=0.0, abs_tol=5.0e-7
    ):
        raise RuntimeError("source dipole angles are not ten approximately equal sectors")

    x = 0.0
    z = 0.0
    tangent = 0.0
    definitions: list[str] = []
    line_names = ["MSTART"]
    bend_count = 0

    for element in elements:
        if element.angle != 0.0:
            bend_count += 1
            radius = element.length / bend_angle
            local_full_x = radius * (math.cos(bend_angle) - 1.0)
            local_full_z = radius * math.sin(bend_angle)
            definitions.append(
                f"{element.name}: SBEND, L = {number(element.length)}, "
                f"ANGLE = BEND_ANGLE, HGAP = 2.0, HAPERT = 2.0, FINT = 0.5,\n"
                f"    X = {number(x)}, Y = 0.0, Z = {number(z)}, "
                f"THETA = {number(tangent)}, PHI = 0.0, PSI = 0.0;"
            )
            x += math.cos(tangent) * local_full_x + math.sin(tangent) * local_full_z
            z += -math.sin(tangent) * local_full_x + math.cos(tangent) * local_full_z
            tangent -= bend_angle
        else:
            definitions.append(
                f"// Source type: {element.kind}\n"
                f"{element.name}: DRIFT, L = {number(element.length)}, "
                f"X = {number(x)}, Y = 0.0, Z = {number(z)}, "
                f"THETA = {number(tangent)}, PHI = 0.0, PSI = 0.0;"
            )
            x += element.length * math.sin(tangent)
            z += element.length * math.cos(tangent)
        line_names.append(element.name)

    line_names.append("MCLOSE")
    closure = math.hypot(x, z)
    if closure > 1.0e-10 or abs(tangent + 2.0 * math.pi) > 1.0e-12:
        raise RuntimeError(
            f"survey does not close: x={x:.16g}, z={z:.16g}, theta={tangent:.16g}"
        )

    wrapped_line = []
    for offset in range(0, len(line_names), 8):
        chunk = ", ".join(line_names[offset : offset + 8])
        if offset + 8 < len(line_names):
            chunk += ","
        wrapped_line.append("    " + chunk)

    total_length = sum(element.length for element in elements)
    straight_length = total_length - sum(
        element.length for element in elements if element.angle != 0.0
    )
    deck = f"""OPTION, PSDUMPFREQ = 50000;
OPTION, STATDUMPFREQ = 1;
OPTION, BOUNDPDESTROY = 10;
OPTION, VERSION = 10900;
OPTION, ASCIIDUMP = TRUE;

Title, string = "ISIS native-SBEND explicit-placement survey";

// Generated by generate_isis_sbend_survey.py from isisring_with_elemedge.in.
// This is a geometry survey, not a physics-equivalent ISIS lattice:
// straight devices are represented by DRIFT envelopes with the same names and L.
// SBEND L is arc length in OPALX, as it is for the source curved MULTIPOLET.
// Source ELEMEDGE length = {number(total_length)} m
// Straight length       = {number(straight_length)} m
// Ten bend arc lengths  = {number(total_length - straight_length)} m
// Exact bend sum        = 2*pi rad
// Analytic closure      = {closure:.3e} m

REAL BEND_ANGLE = 2.0 * PI / 10.0;

MSTART: MONITOR, L = 0.0, X = 0.0, Y = 0.0, Z = 0.0,
    THETA = 0.0, PHI = 0.0, PSI = 0.0, TYPE = TEMPORAL,
    OUTFN = "isis_sbend_survey_start";

{chr(10).join(definitions)}

MCLOSE: MONITOR, L = 0.0, X = {number(x)}, Y = 0.0, Z = {number(z)},
    THETA = {number(tangent)}, PHI = 0.0, PSI = 0.0, TYPE = TEMPORAL,
    OUTFN = "isis_sbend_survey_close";

ISIS_SURVEY: LINE = (\n{chr(10).join(wrapped_line)}\n);

Dist0: DISTRIBUTION, TYPE = FROMFILE,
    FNAME = "isis_sbend_survey_distribution.txt", NPARTDIST = 1;
ES0: EMISSIONSOURCE, DISTRIBUTION = Dist0;
Sources0: EMISSIONSOURCELIST = (ES0);
FS0: FIELDSOLVER, TYPE = NONE, NX = 16, NY = 16, NZ = 16,
    PARFFTX = true, PARFFTY = true, PARFFTZ = true,
    BCFFTX = open, BCFFTY = open, BCFFTZ = open,
    BBOXINCR = 1, GREENSF = INTEGRATED;
BEAM0: BEAM, PARTICLE = PROTON, NALLOC = 1,
    BCHARGE = 1.602176634e-19, SOURCES = Sources0, CHARGE = 1;

// One step is sufficient to construct and export the complete placement survey.
TRACK, LINE = ISIS_SURVEY, BEAM = BEAM0, MAXSTEPS = 1,
    DT = 1.0e-10, ZSTOP = 1.0e6;
RUN, METHOD = "PARALLEL", FIELDSOLVER = FS0;
ENDTRACK;

QUIT;
"""
    OUTPUT.write_text(deck)
    print(f"wrote {OUTPUT}")
    print(f"elements={len(elements)}, bends={bend_count}, source_length_m={total_length:.12g}")
    print(f"closure_x_m={x:.16g}, closure_z_m={z:.16g}, closure_norm_m={closure:.3e}")


if __name__ == "__main__":
    main()
