#!/usr/bin/env python3
"""Generate the compact ELEMEDGE/RING ISIS fixture from the golden deck."""

from __future__ import annotations

import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "orig" / "isis_sbend_survey.in"
OUTPUT = HERE / "isis_sbend_ring.in"


def main() -> None:
    source = REFERENCE.read_text()
    definitions = {
        match.group(1).upper(): (match.group(2).upper(), match.group(3))
        for match in re.finditer(
            r"^([A-Za-z][A-Za-z0-9_]*):\s*(DRIFT|SBEND),\s*(.*?);",
            source,
            re.MULTILINE | re.DOTALL,
        )
    }
    line_match = re.search(
        r"^ISIS_SURVEY:\s*LINE\s*=\s*\((.*?)\);",
        source,
        re.MULTILINE | re.DOTALL,
    )
    if line_match is None:
        raise RuntimeError("ISIS_SURVEY LINE was not found in the golden deck")

    members = [name.strip() for name in line_match.group(1).split(",")]
    lattice_members = members[1:-1]
    if members[0].upper() != "MSTART" or members[-1].upper() != "MCLOSE":
        raise RuntimeError("golden line does not have the expected closure monitors")

    element_lines: list[str] = []
    elemedge = 0.0
    for name in lattice_members:
        kind, attributes = definitions[name.upper()]
        length_match = re.search(r"\bL\s*=\s*([0-9.eE+-]+)", attributes)
        if length_match is None:
            raise RuntimeError(f"{name} has no numeric length")
        length = float(length_match.group(1))
        if kind == "SBEND":
            element_lines.append(
                f"{name}: SBEND, L = {length:.12g}, ANGLE = BEND_ANGLE, "
                "HGAP = 2.0, HAPERT = 2.0, FINT = 0.5, "
                f"ELEMEDGE = {elemedge:.12g};"
            )
        else:
            element_lines.append(
                f"{name}: DRIFT, L = {length:.12g}, ELEMEDGE = {elemedge:.12g};"
            )
        elemedge += length

    wrapped_members = []
    for offset in range(0, len(members), 8):
        suffix = "," if offset + 8 < len(members) else ""
        wrapped_members.append("    " + ", ".join(members[offset : offset + 8]) + suffix)

    deck = f"""OPTION, PSDUMPFREQ = 50000;
OPTION, STATDUMPFREQ = 1;
OPTION, BOUNDPDESTROY = 10;
OPTION, VERSION = 10900;
OPTION, ASCIIDUMP = TRUE;

Title, string = "ISIS native-SBEND ELEMEDGE RING survey";

// Compact counterpart of orig/isis_sbend_survey.in.
// It differs only in placement representation and the closed-topology keyword.
REAL BEND_ANGLE = 2.0 * PI / 10.0;

MSTART: MONITOR, L = 0.0, ELEMEDGE = 0.0, TYPE = TEMPORAL,
    OUTFN = "isis_sbend_ring_start";

{chr(10).join(element_lines)}

MCLOSE: MONITOR, L = 0.0, ELEMEDGE = {elemedge:.12g}, TYPE = TEMPORAL,
    OUTFN = "isis_sbend_ring_close";

ISIS_RING: RING = (
{chr(10).join(wrapped_members)}
);

Dist0: DISTRIBUTION, TYPE = FROMFILE,
    FNAME = "orig/isis_sbend_survey_distribution.txt", NPARTDIST = 1;
ES0: EMISSIONSOURCE, DISTRIBUTION = Dist0;
Sources0: EMISSIONSOURCELIST = (ES0);
FS0: FIELDSOLVER, TYPE = NONE, NX = 16, NY = 16, NZ = 16,
    PARFFTX = true, PARFFTY = true, PARFFTZ = true,
    BCFFTX = open, BCFFTY = open, BCFFTZ = open,
    BBOXINCR = 1, GREENSF = INTEGRATED;
BEAM0: BEAM, PARTICLE = PROTON, NALLOC = 1,
    BCHARGE = 1.602176634e-19, SOURCES = Sources0, CHARGE = 1;

TRACK, LINE = ISIS_RING, BEAM = BEAM0, MAXSTEPS = 1,
    DT = 1.0e-10, ZSTOP = 1.0e6;
RUN, METHOD = "PARALLEL", FIELDSOLVER = FS0;
ENDTRACK;

QUIT;
"""
    OUTPUT.write_text(deck)
    print(f"wrote {OUTPUT}")
    print(f"members={len(lattice_members)}, circumference_m={elemedge:.12g}")


if __name__ == "__main__":
    main()
