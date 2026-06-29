#!/usr/bin/env python3
"""Run truncated chicane decks with one terminal temporal monitor per section."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


SECTIONS = {
    "after_b1": {
        "edge": "D_PRE + LB",
        "line": "Line1: LINE = (DPRE, B1, MSEC);",
    },
    "after_b2": {
        "edge": "D_PRE + 2 * LB + D_CH",
        "line": "Line1: LINE = (DPRE, B1, D12, B2, MSEC);",
    },
    "after_b3": {
        "edge": "D_PRE + 3 * LB + D_CH + C_CH",
        "line": "Line1: LINE = (DPRE, B1, D12, B2, D23, B3, MSEC);",
    },
    "after_b4": {
        "edge": "S_EXIT",
        "line": "Line1: LINE = (DPRE, B1, D12, B2, D23, B3, D34, B4, MSEC);",
    },
}


def copy_run_dir(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(
        src,
        dst,
        ignore=shutil.ignore_patterns("*.h5", "*.stat", "*.loss", "timing.dat", "data"),
    )


def patch_deck(path: Path, section: str, edge: str, line: str, out_prefix: str) -> None:
    text = path.read_text(encoding="utf-8")
    monitor = (
        f"MSEC: MONITOR,\n"
        f"    L = 0.01,\n"
        f"    ELEMEDGE = {edge},\n"
        f"    TYPE = TEMPORAL,\n"
        f"    OUTFN = \"{out_prefix}_{section}\";\n"
    )
    text = text.replace("MEXIT: MONITOR,", monitor + "\nMEXIT: MONITOR,")
    text = text.replace("Line1: LINE = (DPRE, B1, D12, B2, D23, B3, D34, B4, MEXIT);", line)
    text = text.replace("ZSTOP = S_EXIT + 0.2;", f"ZSTOP = {edge} + 0.2;")
    path.write_text(text, encoding="utf-8")


def run_command(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    with (cwd / "section_terminal_monitor_run.log").open("w", encoding="utf-8") as log:
        subprocess.run(cmd, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=default_base)
    parser.add_argument("--opalx-exe", type=Path, default=Path("/Users/adelmann/git/opalx/build/src/opalx"))
    parser.add_argument("--opal-exe", type=Path, default=Path("/Users/adelmann/OPAL-2022.1/bin/opal"))
    parser.add_argument("--opal-prefix", type=Path, default=Path("/Users/adelmann/OPAL-2022.1"))
    parser.add_argument("--only", choices=tuple(SECTIONS), action="append")
    args = parser.parse_args()

    selected = args.only or list(SECTIONS)
    env = os.environ.copy()
    env["OPAL_PREFIX"] = str(args.opal_prefix)

    for section in selected:
        spec = SECTIONS[section]
        opalx_dst = args.base_dir / f"section_monitor_{section}" / "opalx_elemedge"
        opal_dst = args.base_dir / f"section_monitor_{section}" / "opal_2022"
        copy_run_dir(args.base_dir / "opalx_elemedge", opalx_dst)
        copy_run_dir(args.base_dir / "opal_2022", opal_dst)
        patch_deck(opalx_dst / "test-chicane-distribution-2.in", section, spec["edge"], spec["line"], "test-chicane-distribution-2")
        patch_deck(opal_dst / "test-chicane-distribution-1_opal.in", section, spec["edge"], spec["line"], "test-chicane-distribution-1_opal")
        run_command([str(args.opalx_exe), "test-chicane-distribution-2.in"], opalx_dst)
        run_command([str(args.opal_exe), "test-chicane-distribution-1_opal.in"], opal_dst, env=env)
        print(section)
        print(opalx_dst)
        print(opal_dst)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
