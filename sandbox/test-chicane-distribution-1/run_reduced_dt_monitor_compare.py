#!/usr/bin/env python3
"""Run reduced-grid or full-grid exit-monitor comparisons for selected timesteps."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE_BASE = ROOT / "comparison" / "dt_scan" / "finite_10000" / "dt_5p0em12"
OPALX_EXE = Path("/Users/adelmann/git/opalx/build/src/opalx")
OPAL_EXE = Path("/Users/adelmann/OPAL-2022.1/bin/opal")
OPAL_PREFIX = Path("/Users/adelmann/OPAL-2022.1")


def dt_tag(dt: float) -> str:
    return f"dt_{dt:.1e}".replace("+", "").replace("-", "m").replace(".", "p")


def selected_distribution(source: Path, transverse_mode: str) -> np.ndarray:
    data = np.loadtxt(source, skiprows=2)
    if transverse_mode == "all":
        return data
    x_values = np.unique(data[:, 0])[::2]
    y_values = np.unique(data[:, 2])[::2]
    if transverse_mode == "symmetric4":
        x_all = np.unique(data[:, 0])
        y_all = np.unique(data[:, 2])
        x_values = x_all[[0, 3, 6, 9]]
        y_values = y_all[[0, 3, 6, 9]]
    elif transverse_mode != "biased5":
        raise ValueError(f"unknown transverse mode: {transverse_mode}")
    keep = np.isin(data[:, 0], x_values) & np.isin(data[:, 2], y_values)
    return data[keep]


def write_distribution(path: Path, data: np.ndarray, with_header: bool) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(data)}\n")
        if with_header:
            stream.write("x px y py z pz\n")
        np.savetxt(stream, data, fmt="%.16e")


def patch_deck(path: Path, dt: float, n_particles: int, statdumpfreq: int) -> None:
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"(STATDUMPFREQ\s*=\s*)[-+0-9]+", rf"\g<1>{statdumpfreq}", text)
    text = re.sub(r"(DT\s*=\s*)[-+0-9.eE]+", rf"\g<1>{dt:.16e}", text)
    text = re.sub(r"(REAL TIME_STEP\s*=\s*)[-+0-9.eE]+", rf"\g<1>{dt:.16e}", text)
    text = re.sub(r"(NPARTDIST\s*=\s*)[-+0-9]+", rf"\g<1>{n_particles}", text)
    text = re.sub(r"(NALLOC\s*=\s*)[-+0-9]+", rf"\g<1>{n_particles}", text)
    text = re.sub(
        r"(REAL NUMBER_OF_PARTICLES\s*=\s*)[-+0-9]+",
        rf"\g<1>{n_particles}",
        text,
    )
    path.write_text(text, encoding="utf-8")


def prepare_run_dir(
        out_dir: Path,
        dt: float,
        statdumpfreq: int,
        transverse_mode: str) -> tuple[Path, Path]:
    run_dir = out_dir / dt_tag(dt)
    opalx_dir = run_dir / "opalx_elemedge"
    opal_dir = run_dir / "opal_2022"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(
        SOURCE_BASE / "opalx_elemedge",
        opalx_dir,
        ignore=shutil.ignore_patterns("*.h5", "*.stat", "*.loss", "timing.dat", "data", "run.log"),
    )
    shutil.copytree(
        SOURCE_BASE / "opal_2022",
        opal_dir,
        ignore=shutil.ignore_patterns("*.h5", "*.stat", "*.loss", "timing.dat", "data", "run.log"),
    )

    data = selected_distribution(
        SOURCE_BASE / "opalx_elemedge" / "test-chicane-distribution-1_distribution.txt",
        transverse_mode,
    )
    write_distribution(opalx_dir / "test-chicane-distribution-1_distribution.txt", data, with_header=True)
    write_distribution(opalx_dir / "test-chicane-distribution-1_opal_distribution.txt", data, with_header=False)
    write_distribution(opal_dir / "test-chicane-distribution-1_distribution.txt", data, with_header=True)
    write_distribution(opal_dir / "test-chicane-distribution-1_opal_distribution.txt", data, with_header=False)

    patch_deck(opalx_dir / "test-chicane-distribution-2.in", dt, len(data), statdumpfreq)
    patch_deck(opal_dir / "test-chicane-distribution-1_opal.in", dt, len(data), statdumpfreq)
    return opalx_dir, opal_dir


def run_command(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    print(f"{cwd}: $ {' '.join(command)}", flush=True)
    with (cwd / "run.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "comparison" / "dt_scan" / "finite_2500",
    )
    parser.add_argument("--dt", type=float, nargs="+", default=[5.0e-12, 2.0e-12])
    parser.add_argument("--statdumpfreq", type=int, default=20)
    parser.add_argument("--transverse-mode", choices=["biased5", "symmetric4", "all"], default="biased5")
    parser.add_argument("--skip-opal", action="store_true")
    parser.add_argument("--opalx-exe", type=Path, default=OPALX_EXE)
    parser.add_argument("--opal-exe", type=Path, default=OPAL_EXE)
    args = parser.parse_args()

    env = os.environ.copy()
    env["OPAL_PREFIX"] = str(OPAL_PREFIX)

    for dt in args.dt:
        opalx_dir, opal_dir = prepare_run_dir(args.out_dir, dt, args.statdumpfreq, args.transverse_mode)
        run_command([str(args.opalx_exe), "test-chicane-distribution-2.in"], opalx_dir)
        if not args.skip_opal:
            run_command([str(args.opal_exe), "test-chicane-distribution-1_opal.in"], opal_dir, env=env)
        print(f"finished {dt_tag(dt)}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
