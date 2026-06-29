#!/usr/bin/env python3
"""Run unmodified chicane decks truncated at selected longitudinal stops."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


STOPS = {
    "after_b1": "D_PRE + LB + 0.05",
    "after_b2": "D_PRE + 2 * LB + D_CH + 0.05",
    "after_b3": "D_PRE + 3 * LB + D_CH + C_CH + 0.05",
    "after_b4": "S_EXIT + 0.05",
}


def copy_run_dir(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(
        src,
        dst,
        ignore=shutil.ignore_patterns("*.h5", "*.stat", "*.loss", "timing.dat", "data"),
    )


def patch_zstop(path: Path, expression: str) -> None:
    text = path.read_text(encoding="utf-8")
    text = text.replace("ZSTOP = S_EXIT + 0.2;", f"ZSTOP = {expression};")
    path.write_text(text, encoding="utf-8")


def run_command(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    with (cwd / "section_stop_run.log").open("w", encoding="utf-8") as log:
        subprocess.run(cmd, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=default_base)
    parser.add_argument("--opalx-exe", type=Path, default=Path("/Users/adelmann/git/opalx/build/src/opalx"))
    parser.add_argument("--opal-exe", type=Path, default=Path("/Users/adelmann/OPAL-2022.1/bin/opal"))
    parser.add_argument("--opal-prefix", type=Path, default=Path("/Users/adelmann/OPAL-2022.1"))
    parser.add_argument("--only", choices=tuple(STOPS), action="append")
    args = parser.parse_args()

    selected = args.only or list(STOPS)
    env = os.environ.copy()
    env["OPAL_PREFIX"] = str(args.opal_prefix)

    for label in selected:
        expression = STOPS[label]
        opalx_dst = args.base_dir / f"section_stop_{label}" / "opalx_elemedge"
        opal_dst = args.base_dir / f"section_stop_{label}" / "opal_2022"
        copy_run_dir(args.base_dir / "opalx_elemedge", opalx_dst)
        copy_run_dir(args.base_dir / "opal_2022", opal_dst)
        patch_zstop(opalx_dst / "test-chicane-distribution-2.in", expression)
        patch_zstop(opal_dst / "test-chicane-distribution-1_opal.in", expression)
        run_command([str(args.opalx_exe), "test-chicane-distribution-2.in"], opalx_dst)
        run_command([str(args.opal_exe), "test-chicane-distribution-1_opal.in"], opal_dst, env=env)
        print(label)
        print(opalx_dst)
        print(opal_dst)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
