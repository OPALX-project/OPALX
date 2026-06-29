#!/usr/bin/env python3
"""Run a configurable-GAP ELEMEDGE diagnostic for the 10000-particle chicane comparison."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
SOURCE_BASE = (
    ROOT
    / "comparison"
    / "dt_scan"
    / "finite_10000"
    / "dt_5p0em12"
)
OPALX_EXE = Path("/Users/adelmann/git/opalx/build/src/opalx")
OPAL_EXE = Path("/Users/adelmann/OPAL-2022.1/bin/opal")
OPAL_PREFIX = Path("/Users/adelmann/OPAL-2022.1")


def format_tag(value: float) -> str:
    return f"{value:.1e}".replace(".", "p").replace("-", "m").replace("+", "")


def patch_input(text: str, gap: float, statdumpfreq: int) -> str:
    gap_text = f"GAP = {gap:.16e}"
    text = re.sub(r"GAP\s*=\s*[-+0-9.eEdD]+", gap_text, text)
    text = re.sub(
        r"(STATDUMPFREQ\s*=\s*)[-+0-9]+",
        rf"\g<1>{statdumpfreq}",
        text,
        flags=re.IGNORECASE,
    )
    return text


def run_command(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    print(f"$ {' '.join(command)}", flush=True)
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (cwd / "run.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout, flush=True)
        raise RuntimeError(f"command failed with exit code {result.returncode}: {' '.join(command)}")


def prepare_run_dir(out_dir: Path, gap: float, statdumpfreq: int, force: bool) -> tuple[Path, Path, Path]:
    opalx_dir = out_dir / "opalx_elemedge"
    opal_dir = out_dir / "opal_2022"
    compare_dir = out_dir / "compare_elemedge_beam_only"
    marker_dir = out_dir / "marker_compare"
    for directory in (opalx_dir, opal_dir, compare_dir, marker_dir):
        directory.mkdir(parents=True, exist_ok=True)

    opalx_input = opalx_dir / "test-chicane-distribution-2.in"
    opal_input = opal_dir / "test-chicane-distribution-1_opal.in"
    if force or not opalx_input.exists():
        template = SOURCE_BASE / "opalx_elemedge" / opalx_input.name
        opalx_input.write_text(
            patch_input(template.read_text(encoding="utf-8"), gap, statdumpfreq),
            encoding="utf-8",
        )
    if force or not opal_input.exists():
        template = SOURCE_BASE / "opal_2022" / opal_input.name
        opal_input.write_text(
            patch_input(template.read_text(encoding="utf-8"), gap, statdumpfreq),
            encoding="utf-8",
        )

    for name in ("test-chicane-distribution-1_distribution.txt", "test-chicane-distribution-1_opal_distribution.txt"):
        src = SOURCE_BASE / "opalx_elemedge" / name
        dst = opalx_dir / name
        if force or not dst.exists():
            shutil.copy2(src, dst)

    for name in ("test-chicane-distribution-1_distribution.txt", "test-chicane-distribution-1_opal_distribution.txt"):
        src = SOURCE_BASE / "opal_2022" / name
        dst = opal_dir / name
        if force or not dst.exists():
            shutil.copy2(src, dst)

    return opalx_dir, opal_dir, compare_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gap", type=float, default=0.0)
    parser.add_argument("--statdumpfreq", type=int, default=20)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-runs", action="store_true")
    args = parser.parse_args()

    out_dir = args.out_dir
    if out_dir is None:
        out_dir = SOURCE_BASE.parent / f"dt_5p0em12_gap{format_tag(args.gap)}_sdf{args.statdumpfreq}"

    opalx_dir, opal_dir, compare_dir = prepare_run_dir(out_dir, args.gap, args.statdumpfreq, args.force)
    if not args.skip_runs:
        if args.force or not (opalx_dir / "test-chicane-distribution-2.stat").exists():
            run_command([str(OPALX_EXE), "test-chicane-distribution-2.in"], opalx_dir)
        if args.force or not (opal_dir / "test-chicane-distribution-1_opal.stat").exists():
            env = dict(os.environ, OPAL_PREFIX=str(OPAL_PREFIX))
            run_command([str(OPAL_EXE), "test-chicane-distribution-1_opal.in"], opal_dir, env=env)

    run_command(
        [
            "/Users/adelmann/.venv-h6/bin/python",
            str(ROOT / "compare_opal_opalx.py"),
            "--beam-only",
            "--opalx-stat",
            str(opalx_dir / "test-chicane-distribution-2.stat"),
            "--opal-stat",
            str(opal_dir / "test-chicane-distribution-1_opal.stat"),
            "--out-dir",
            str(compare_dir),
            "--plots",
        ],
        out_dir,
    )
    run_command(
        [
            "/Users/adelmann/.venv-h6/bin/python",
            str(ROOT / "compare_element_markers.py"),
            "--opalx-markers",
            str(opalx_dir / "data" / "test-chicane-distribution-2_ElementPositions.txt"),
            "--opal-markers",
            str(opal_dir / "data" / "test-chicane-distribution-1_opal_ElementPositions.txt"),
            "--out-dir",
            str(out_dir / "marker_compare"),
        ],
        out_dir,
    )
    print(f"Wrote {compare_dir / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
