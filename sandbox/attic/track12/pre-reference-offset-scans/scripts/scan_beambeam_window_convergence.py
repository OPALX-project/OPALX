#!/usr/bin/env python3
"""Scan BeamBeam window aperture/grid convergence for the track12 pair-4 kick.

The study intentionally uses the slide-timed pair-4 one-step witness kick because
that is the smallest reproducible case that exercises the BeamBeam frozen-window
mesh and passive witness gather path.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TRACK12_MODULE = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_TEMPLATE_DIR = (
    ROOT / "sandbox" / "track12particles" / "opalx" / "timing_pairs_100" / "pair4"
)
DEFAULT_TEMPLATE = (
    DEFAULT_TEMPLATE_DIR / "track12_pair4_slide_timed_one_source_1fs_400k_100steps.in"
)
DEFAULT_OUTPUT_DIR = (
    ROOT / "sandbox" / "track12particles" / "opalx" / "beambeam_window_scan"
)
DEFAULT_OPALX = ROOT / "build_openmp" / "src" / "opalx"


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TRACK12_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manufactured_pair4_field(source_selection: str) -> tuple[np.ndarray, np.ndarray]:
    module = load_track12_module()
    initial = pd.read_csv(
        ROOT / "sandbox" / "track12particles" / "outputs" / "track12_initial_conditions.csv"
    )
    row = initial.loc[initial["pair"].eq(4) & initial["species"].eq("electron")].iloc[0]
    sources = module.make_anisotropic_sources(
        module.DEFAULT_SOURCE_KINETIC_MEV,
        module.DEFAULT_SOURCE_CHARGE_C,
        module.SIGMA_XY_M,
        module.SIGMA_XY_M,
        module.SIGMA_Z_M,
        source_selection,
        0.0,
        module.DEFAULT_MC_SOURCE_PARTICLES,
        module.DEFAULT_MC_SEED,
    )
    shifted = []
    for source in sources:
        center = source.center_t0_m.copy()
        center[2] -= float(row["t"])
        shifted.append(
            module.RigidAnisotropicGaussianSource(
                source.name,
                charge_C=source.charge_C,
                sigma_lab_m=(module.SIGMA_XY_M, module.SIGMA_XY_M, module.SIGMA_Z_M),
                beta_z=source.beta_z,
                center_t0_m=center,
                t0_s=0.0,
            )
        )
    position = row[["x", "y", "s"]].to_numpy(dtype=float)
    return module.anisotropic_total_lab_fields(position, 0.0, tuple(shifted))


def patch_input(template: str, aperture_m: float, nxy: int, nz: int, particles: int | None) -> str:
    text = template
    text = re.sub(r"OPTION, PSDUMPFREQ = \d+;", "OPTION, PSDUMPFREQ = 2;", text)
    text = re.sub(r"OPTION, STATDUMPFREQ = \d+;", "OPTION, STATDUMPFREQ = 2;", text)
    text = re.sub(r"REAL n_steps = \d+;", "REAL n_steps = 1;", text)
    text = re.sub(r"REAL NXY = \d+;", f"REAL NXY = {nxy};", text)
    text = re.sub(r"REAL NZ = \d+;", f"REAL NZ = {nz};", text)
    text = re.sub(
        r'APERTURE = "CIRCLE\([^)]*\)"',
        f'APERTURE = "CIRCLE({aperture_m:.12g})"',
        text,
    )
    if particles is not None:
        text = re.sub(
            r"REAL primary_macroparticles = \d+;",
            f"REAL primary_macroparticles = {particles};",
            text,
        )
    return text


def prepare_case(
    output_dir: Path,
    template_path: Path,
    aperture_m: float,
    nxy: int,
    nz: int,
    particles: int | None,
) -> tuple[Path, str]:
    particle_tag = "" if particles is None else f"_np{particles // 1000}k"
    case_name = f"ap{aperture_m * 1.0e6:.0f}um_nxy{nxy}_nz{nz}{particle_tag}"
    case_dir = output_dir / case_name
    case_dir.mkdir(parents=True, exist_ok=True)

    for filename in ("track12_pair4_electron.fromfile", "track12_pair4_positron.fromfile"):
        shutil.copy2(template_path.parent / filename, case_dir / filename)

    input_text = patch_input(template_path.read_text(encoding="utf-8"), aperture_m, nxy, nz, particles)
    input_path = case_dir / f"{case_name}.in"
    input_path.write_text(input_text, encoding="utf-8")
    return case_dir, input_path.name


def run_case(
    opalx: Path,
    case_dir: Path,
    input_name: str,
    omp_threads: int,
    force: bool,
) -> dict[str, object]:
    stem = Path(input_name).stem
    kick_csv = case_dir / f"{stem}_kicks.csv"
    if kick_csv.exists() and not force:
        return {"status": "cached", "returncode": 0, "runtime_s": 0.0, "kick_csv": kick_csv}

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["OPALX_BB_WITNESS_KICK_CSV"] = kick_csv.name
    env["OPALX_BB_WITNESS_KICK_STEPS"] = "1"

    start = time.monotonic()
    completed = subprocess.run(
        [str(opalx), input_name],
        cwd=case_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    runtime_s = time.monotonic() - start
    (case_dir / f"{stem}.log").write_text(completed.stdout, encoding="utf-8")
    return {
        "status": "ok" if completed.returncode == 0 else "failed",
        "returncode": completed.returncode,
        "runtime_s": runtime_s,
        "kick_csv": kick_csv,
    }


def read_case_result(kick_csv: Path) -> dict[str, float]:
    kicks = pd.read_csv(kick_csv)
    row = kicks.loc[kicks["step"].eq(0) & kicks["container"].eq(1)].iloc[0]
    return {
        "Ex_opalx": float(row["Ex_V_per_m"]),
        "Ey_opalx": float(row["Ey_V_per_m"]),
        "Ez_opalx": float(row["Ez_V_per_m"]),
        "Bx_opalx": float(row["Bx_T"]),
        "By_opalx": float(row["By_T"]),
        "Bz_opalx": float(row["Bz_T"]),
        "E_abs_opalx": float(row["E_abs_V_per_m"]),
        "B_abs_opalx": float(row["B_abs_T"]),
        "dPx_electron": float(row["dPx"]),
        "dPy_electron": float(row["dPy"]),
        "dPz_electron": float(row["dPz"]),
    }


def write_plots(results: pd.DataFrame, output_dir: Path) -> None:
    cache_root = output_dir / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    ok = results.loc[results["status"].isin(["ok", "cached"])].copy()
    if ok.empty:
        return

    for quantity, ylabel, filename in [
        ("Ex_opalx", r"$E_x$ [V/m]", "pair4_ex_vs_cell_size.png"),
        ("By_opalx", r"$B_y$ [T]", "pair4_by_vs_cell_size.png"),
        ("E_abs_ratio", r"$|E|_\mathrm{OPALX}/|E|_\mathrm{manuf.}$", "pair4_e_ratio_vs_cell_size.png"),
    ]:
        fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
        for nxy, group in ok.sort_values("cell_size_um").groupby("nxy"):
            ax.plot(
                group["cell_size_um"],
                group[quantity],
                marker="o",
                lw=1.2,
                label=f"NXY={nxy}",
            )
        if quantity == "Ex_opalx":
            ax.axhline(ok["Ex_manufactured"].iloc[0], color="0.2", ls="--", lw=1.1, label="manufactured")
        if quantity == "By_opalx":
            ax.axhline(ok["By_manufactured"].iloc[0], color="0.2", ls="--", lw=1.1, label="manufactured")
        if quantity == "E_abs_ratio":
            ax.axhline(1.0, color="0.2", ls="--", lw=1.1, label="manufactured")
        ax.set_xscale("log")
        ax.set_xlabel(r"BeamBeam transverse cell size [$\mu$m]")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", color="0.88", lw=0.6)
        ax.legend(fontsize=8)
        fig.savefig(output_dir / filename, dpi=220)
        plt.close(fig)

    if ok["primary_macroparticles"].nunique() > 1:
        for quantity, ylabel, filename in [
            ("Ex_opalx", r"$E_x$ [V/m]", "pair4_ex_vs_particles.png"),
            ("By_opalx", r"$B_y$ [T]", "pair4_by_vs_particles.png"),
            (
                "E_abs_ratio",
                r"$|E|_\mathrm{OPALX}/|E|_\mathrm{manuf.}$",
                "pair4_e_ratio_vs_particles.png",
            ),
        ]:
            fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
            for cell_size, group in ok.sort_values("primary_macroparticles").groupby("cell_size_um"):
                ax.plot(
                    group["primary_macroparticles"],
                    group[quantity],
                    marker="o",
                    lw=1.2,
                    label=f"cell={cell_size:g} um",
                )
            if quantity == "Ex_opalx":
                ax.axhline(
                    ok["Ex_manufactured"].iloc[0],
                    color="0.2",
                    ls="--",
                    lw=1.1,
                    label="manufactured",
                )
            if quantity == "By_opalx":
                ax.axhline(
                    ok["By_manufactured"].iloc[0],
                    color="0.2",
                    ls="--",
                    lw=1.1,
                    label="manufactured",
                )
            if quantity == "E_abs_ratio":
                ax.axhline(1.0, color="0.2", ls="--", lw=1.1, label="manufactured")
            ax.set_xscale("log")
            ax.set_xlabel("Primary macroparticles")
            ax.set_ylabel(ylabel)
            ax.grid(True, which="both", color="0.88", lw=0.6)
            ax.legend(fontsize=8)
            fig.savefig(output_dir / filename, dpi=220)
            plt.close(fig)


def parse_float_list(text: str) -> list[float]:
    return [float(item) for item in text.split(",") if item.strip()]


def parse_int_list(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--apertures-m", default="0.001,0.0002,0.0001,0.00005,0.00002")
    parser.add_argument("--nxy-values", default="32,64")
    parser.add_argument("--nz", type=int, default=64)
    parser.add_argument("--particles", type=int)
    parser.add_argument(
        "--particle-values",
        help="Comma-separated primary macroparticle counts. Overrides --particles when set.",
    )
    parser.add_argument("--omp-threads", type=int, default=8)
    parser.add_argument("--source-selection", choices=("copropagating", "oncoming"), default="copropagating")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    e_manuf, b_manuf = manufactured_pair4_field(args.source_selection)
    apertures = parse_float_list(args.apertures_m)
    nxy_values = parse_int_list(args.nxy_values)
    particle_values = (
        parse_int_list(args.particle_values)
        if args.particle_values
        else [args.particles if args.particles is not None else 400000]
    )

    rows: list[dict[str, object]] = []
    for particles in particle_values:
        particle_arg = None if particles == 400000 and args.particles is None and not args.particle_values else particles
        for aperture in apertures:
            for nxy in nxy_values:
                case_dir, input_name = prepare_case(
                    args.output_dir, args.template, aperture, nxy, args.nz, particle_arg
                )
                run = run_case(args.opalx, case_dir, input_name, args.omp_threads, args.force)
                row: dict[str, object] = {
                    "case": case_dir.name,
                    "aperture_m": aperture,
                    "aperture_um": aperture * 1.0e6,
                    "nxy": nxy,
                    "nz": args.nz,
                    "primary_macroparticles": particles,
                    "cell_size_um": 2.0 * aperture / nxy * 1.0e6,
                    "source_selection": args.source_selection,
                    "status": run["status"],
                    "returncode": run["returncode"],
                    "runtime_s": run["runtime_s"],
                    "Ex_manufactured": e_manuf[0],
                    "Ey_manufactured": e_manuf[1],
                    "Ez_manufactured": e_manuf[2],
                    "Bx_manufactured": b_manuf[0],
                    "By_manufactured": b_manuf[1],
                    "Bz_manufactured": b_manuf[2],
                    "E_abs_manufactured": float(np.linalg.norm(e_manuf)),
                    "B_abs_manufactured": float(np.linalg.norm(b_manuf)),
                }
                kick_csv = run["kick_csv"]
                if row["status"] in ("ok", "cached") and Path(kick_csv).exists():
                    row.update(read_case_result(Path(kick_csv)))
                    row["E_abs_ratio"] = row["E_abs_opalx"] / row["E_abs_manufactured"]
                    row["B_abs_ratio"] = row["B_abs_opalx"] / row["B_abs_manufactured"]
                    row["Ex_error_V_per_m"] = row["Ex_opalx"] - row["Ex_manufactured"]
                    row["By_error_T"] = row["By_opalx"] - row["By_manufactured"]
                rows.append(row)
                print(
                    f"{case_dir.name}: {row['status']} return={row['returncode']} "
                    f"cell={row['cell_size_um']:.3g} um particles={particles}",
                    flush=True,
                )

    results = pd.DataFrame(rows)
    csv_path = args.output_dir / "pair4_beambeam_window_scan.csv"
    results.to_csv(csv_path, index=False)
    write_plots(results, args.output_dir)
    print(f"wrote {csv_path}")
    print(f"wrote plots in {args.output_dir}")


if __name__ == "__main__":
    main()
