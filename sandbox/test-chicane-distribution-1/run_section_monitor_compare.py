#!/usr/bin/env python3
"""Run OPALX/old-OPAL chicane variants with section temporal monitors."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


MONITORS = [
    ("MB1", "D_PRE + LB"),
    ("MD12", "D_PRE + LB + D_CH"),
    ("MB2", "D_PRE + 2 * LB + D_CH"),
    ("MD23", "D_PRE + 2 * LB + D_CH + C_CH"),
    ("MB3", "D_PRE + 3 * LB + D_CH + C_CH"),
    ("MD34", "D_PRE + 3 * LB + 2 * D_CH + C_CH"),
    ("MB4", "S_EXIT"),
]


def copy_run_dir(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(
        src,
        dst,
        ignore=shutil.ignore_patterns("*.h5", "*.stat", "*.loss", "timing.dat", "data"),
    )


def monitor_block(name: str, edge: str, out_prefix: str) -> str:
    return (
        f"{name}: MONITOR,\n"
        f"    L = 0.01,\n"
        f"    ELEMEDGE = {edge},\n"
        f"    TYPE = TEMPORAL,\n"
        f"    OUTFN = \"{out_prefix}_{name.lower()}\";\n"
    )


def patch_deck(path: Path, out_prefix: str) -> None:
    text = path.read_text(encoding="utf-8")
    blocks = "\n".join(monitor_block(name, edge, out_prefix) for name, edge in MONITORS)
    text = text.replace("MEXIT: MONITOR,", blocks + "\nMEXIT: MONITOR,")
    old_line = "Line1: LINE = (DPRE, B1, D12, B2, D23, B3, D34, B4, MEXIT);"
    new_line = (
        "Line1: LINE = (DPRE, B1, MB1, D12, MD12, B2, MB2, "
        "D23, MD23, B3, MB3, D34, MD34, B4, MB4, MEXIT);"
    )
    text = text.replace(old_line, new_line)
    path.write_text(text, encoding="utf-8")


def run_command(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    with (cwd / "section_monitor_run.log").open("w", encoding="utf-8") as log:
        subprocess.run(cmd, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=default_base)
    parser.add_argument("--opalx-exe", type=Path, default=Path("/Users/adelmann/git/opalx/build/src/opalx"))
    parser.add_argument("--opal-exe", type=Path, default=Path("/Users/adelmann/OPAL-2022.1/bin/opal"))
    parser.add_argument("--opal-prefix", type=Path, default=Path("/Users/adelmann/OPAL-2022.1"))
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()

    opalx_src = args.base_dir / "opalx_elemedge"
    opal_src = args.base_dir / "opal_2022"
    opalx_dst = args.base_dir / "opalx_elemedge_section_monitors"
    opal_dst = args.base_dir / "opal_2022_section_monitors"

    copy_run_dir(opalx_src, opalx_dst)
    copy_run_dir(opal_src, opal_dst)
    patch_deck(opalx_dst / "test-chicane-distribution-2.in", "test-chicane-distribution-2")
    patch_deck(opal_dst / "test-chicane-distribution-1_opal.in", "test-chicane-distribution-1_opal")

    if not args.skip_run:
        run_command([str(args.opalx_exe), "test-chicane-distribution-2.in"], opalx_dst)
        env = os.environ.copy()
        env["OPAL_PREFIX"] = str(args.opal_prefix)
        run_command([str(args.opal_exe), "test-chicane-distribution-1_opal.in"], opal_dst, env=env)

    print(opalx_dst)
    print(opal_dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
