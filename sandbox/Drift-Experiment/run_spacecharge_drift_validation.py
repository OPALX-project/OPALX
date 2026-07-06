#!/usr/bin/env python3
"""Run the c0-only drift space-charge field validation."""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DEFAULT_OPALX = REPO / "build_openmp" / "src" / "opalx"
DEFAULT_PYTHON = Path.home() / ".venv-h6" / "bin" / "python"
COMPARATOR = REPO / "sandbox" / "compare-e-fields" / "compare_gaussian_pic_fields.py"

SIGMA_M = 1.0e-3
TOTAL_CHARGE_C = -1.0e-9


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=HERE / "spacecharge_drift_30cm.in")
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    parser.add_argument("--force", action="store_true", help="Replace an existing output directory")
    parser.add_argument("--skip-run", action="store_true", help="Reuse existing raw field CSV")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--success-sigma-factor",
        type=float,
        default=3.0,
        help="Pass threshold multiplier for the nominal 1/sqrt(N) statistical scale",
    )
    return parser.parse_args()


def prepare_output(output_dir: Path, force: bool) -> tuple[Path, Path, Path]:
    if output_dir.exists() and force:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "raw"
    inputs_dir = output_dir / "inputs"
    comparison_dir = output_dir / "comparison"
    raw_dir.mkdir(exist_ok=True)
    inputs_dir.mkdir(exist_ok=True)
    comparison_dir.mkdir(exist_ok=True)
    return raw_dir, inputs_dir, comparison_dir


def run_opalx(args: argparse.Namespace, raw_csv: Path, work_dir: Path) -> None:
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.threads)
    env["OPALX_SC_FIELD_CSV"] = str(raw_csv)
    env["OPALX_SC_FIELD_STEPS"] = "0"
    subprocess.run(
        [str(args.opalx), str(args.input)],
        cwd=work_dir,
        env=env,
        check=True,
        stdout=(work_dir / "opalx.out").open("w", encoding="utf-8"),
        stderr=subprocess.STDOUT,
    )


def read_raw_field_csv(raw_csv: Path) -> pd.DataFrame:
    candidates = sorted(raw_csv.parent.glob(raw_csv.name + ".rank*"))
    if raw_csv.exists():
        candidates.insert(0, raw_csv)
    if not candidates:
        raise FileNotFoundError(f"No OPALX space-charge field CSV found for {raw_csv}")
    return pd.concat((pd.read_csv(path) for path in candidates), ignore_index=True)


def export_pic_fields(raw: pd.DataFrame, pic_csv: Path) -> pd.DataFrame:
    step0 = raw.loc[raw["step"].eq(0) & raw["container"].eq(0)].copy()
    if step0.empty:
        raise RuntimeError("The OPALX diagnostic did not contain any c0 step-0 field rows")
    step0 = step0.sort_values(["rank", "local_index"]).reset_index(drop=True)
    pic = pd.DataFrame(
        {
            "x_m": step0["x_m"],
            "y_m": step0["y_m"],
            "z_m": step0["z_m"],
            "s": step0["id"],
            "Ex": step0["Ex_V_per_m"],
            "Ey": step0["Ey_V_per_m"],
            "Ez": step0["Ez_V_per_m"],
            "E_abs": step0["E_abs_V_per_m"],
        }
    )
    pic.to_csv(pic_csv, index=False)
    return pic


def run_comparator(
        args: argparse.Namespace, pic_csv: Path, comparison_dir: Path,
        center: tuple[float, float, float]) -> None:
    env = dict(os.environ)
    env.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    env.setdefault("XDG_CACHE_HOME", str(args.output_dir / ".cache"))
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(args.python),
            str(COMPARATOR),
            "--pic",
            str(pic_csv),
            "--sigma",
            str(SIGMA_M),
            f"--charge={TOTAL_CHARGE_C}",
            "--center",
            str(center[0]),
            str(center[1]),
            str(center[2]),
            "--outdir",
            str(comparison_dir),
            "--max-plot-points",
            "50000",
        ],
        env=env,
        check=True,
    )


def write_summary(
        args: argparse.Namespace, pic: pd.DataFrame, comparison_dir: Path,
        center: tuple[float, float, float]) -> None:
    metrics_path = comparison_dir / "comparison_metrics.csv"
    metrics = pd.read_csv(metrics_path)
    n_particles = int(metrics.loc[0, "n_particles"])
    relative_l2 = float(metrics.loc[0, "relative_vector_L2_vs_analytic"])
    nominal_statistical_scale = 1.0 / math.sqrt(n_particles)
    threshold = args.success_sigma_factor * nominal_statistical_scale
    passed = relative_l2 <= threshold

    summary = pd.DataFrame(
        [
            {
                "n_particles": n_particles,
                "sigma_m": SIGMA_M,
                "analytic_total_charge_C": TOTAL_CHARGE_C,
                "relative_vector_L2_vs_analytic": relative_l2,
                "nominal_statistical_relative_scale_1_over_sqrt_N": nominal_statistical_scale,
                "success_sigma_factor": args.success_sigma_factor,
                "success_threshold": threshold,
                "passed": passed,
                "mean_x_m": float(pic["x_m"].mean()),
                "mean_y_m": float(pic["y_m"].mean()),
                "mean_z_m": float(pic["z_m"].mean()),
                "analytic_center_x_m": center[0],
                "analytic_center_y_m": center[1],
                "analytic_center_z_m": center[2],
                "rms_x_m": float(pic["x_m"].std(ddof=0)),
                "rms_y_m": float(pic["y_m"].std(ddof=0)),
                "rms_z_m": float(pic["z_m"].std(ddof=0)),
            }
        ]
    )
    summary.to_csv(args.output_dir / "validation_summary.csv", index=False)
    (args.output_dir / "README.md").write_text(
        "\n".join(
            [
                "# c0-only drift space-charge validation",
                "",
                f"Input deck: `{args.input.relative_to(REPO)}`",
                "",
                "The OPALX run contains only container 0, a 30 cm drift, one integration step,",
                "and an isotropic Gaussian electron bunch. `OPALX_SC_FIELD_CSV` dumps the",
                "normal space-charge field after `computeSelfFields()` and transformation back",
                "to the reference frame. The dumped particle fields are compared with",
                "`sandbox/compare-e-fields/compare_gaussian_pic_fields.py`. The analytic",
                "Gaussian is centered on the dumped c0 centroid in the same reference frame.",
                "",
                "Files:",
                "",
                "- `raw/spacecharge_fields.csv`: OPALX diagnostic field dump.",
                "- `inputs/pic_fields.csv`: comparator input adapted from the OPALX dump.",
                "- `comparison/`: output from `compare_gaussian_pic_fields.py`.",
                "- `validation_summary.csv`: pass/fail summary.",
                "",
                "Success criterion used by this runner:",
                "",
                f"`relative_vector_L2_vs_analytic <= {args.success_sigma_factor:g} / sqrt(N)`.",
                "",
                f"Observed relative L2: `{relative_l2:.12e}`",
                f"Threshold: `{threshold:.12e}`",
                f"Passed: `{passed}`",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


def main() -> None:
    args = parse_args()
    raw_dir, inputs_dir, comparison_dir = prepare_output(args.output_dir, args.force)
    work_dir = args.output_dir / "run"
    work_dir.mkdir(exist_ok=True)
    run_input = work_dir / args.input.name
    shutil.copy2(args.input, run_input)

    raw_csv = raw_dir / "spacecharge_fields.csv"
    if not args.skip_run:
        run_opalx(args, raw_csv, work_dir)
    raw = read_raw_field_csv(raw_csv)
    pic = export_pic_fields(raw, inputs_dir / "pic_fields.csv")
    center = (
        float(pic["x_m"].mean()),
        float(pic["y_m"].mean()),
        float(pic["z_m"].mean()),
    )
    run_comparator(args, inputs_dir / "pic_fields.csv", comparison_dir, center)
    write_summary(args, pic, comparison_dir, center)


if __name__ == "__main__":
    main()
