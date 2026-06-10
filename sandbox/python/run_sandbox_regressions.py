#!/usr/bin/env python3
"""Run lightweight regression checks for the BeamBeam sandbox.

The checks deliberately consume the same sandbox artifacts used in the manual
workflow:

* manufactured BeamBeam field dumps in ``data/sandbox``;
* gamma-gamma e-/e+ FROMFILE inputs plus OPALX ``gamma_gamma_pairs-2`` outputs;
* the OPALX-vs-IMPACT drift comparison CSV/stat pair.

Use ``--update-baseline`` only after intentionally accepting new physics output.
The default baseline and current-metrics CSV live in ``sandbox/note`` so this
script remains usable after removing the old ``sandbox/regression`` directory.

Typical use from the repository root::

    source .venv-h6/bin/activate
    python sandbox/python/run_sandbox_regressions.py
    python sandbox/python/run_sandbox_regressions.py --check fields
    python sandbox/python/run_sandbox_regressions.py --run-opalx --opalx-exe build_openmp/src/opalx
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
NOTE_DIR = ROOT / "sandbox/note"
DEFAULT_BASELINE = NOTE_DIR / "sandbox_regression_baseline.json"
DEFAULT_CSV = NOTE_DIR / "current_metrics.csv"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_numpy():
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        raise RuntimeError("numpy is required; activate .venv-h6 first") from exc
    return np


def flatten(prefix: str, value: Any) -> dict[str, float | str]:
    rows: dict[str, float | str] = {}
    if isinstance(value, dict):
        for key, item in value.items():
            next_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.update(flatten(next_prefix, item))
    elif isinstance(value, (int, float)):
        rows[prefix] = float(value)
    elif isinstance(value, str):
        rows[prefix] = value
    else:
        rows[prefix] = str(value)
    return rows


def read_sdds_stat(path: Path) -> tuple[pd.DataFrame, list[str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    lines = path.read_text(errors="replace").splitlines()

    columns: list[str] = []
    data_start: int | None = None
    idx = 0
    while idx < len(lines):
        tag = lines[idx].strip()
        if tag == "&column":
            idx += 1
            name: str | None = None
            while idx < len(lines) and lines[idx].strip() != "&end":
                line = lines[idx].strip().rstrip(",")
                if line.startswith("name="):
                    name = line.split("=", 1)[1].strip().strip('"')
                idx += 1
            columns.append(name or f"col{len(columns) + 1}")
        elif tag == "&data":
            while idx < len(lines) and lines[idx].strip() != "&end":
                idx += 1
            data_start = idx + 1
            break
        idx += 1

    if data_start is None or not columns:
        raise ValueError(f"{path}: not an ASCII SDDS/stat file")

    rows: list[list[float]] = []
    for raw in lines[data_start:]:
        parts = raw.split()
        if len(parts) != len(columns):
            continue
        try:
            rows.append([float(part) for part in parts])
        except ValueError:
            continue

    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    return pd.DataFrame(rows, columns=columns), columns


def field_metrics() -> dict[str, Any]:
    np = load_numpy()
    manufactured = load_module(
        "beam_beam_manufactured_solution",
        ROOT / "sandbox/python/beam-beam-manufactured-solution.py",
    )

    rho_path = ROOT / "data/sandbox/BeamBeam-static-1V-RHO_scalar-beambeam_rho_pre-000003.dat"
    phi_path = ROOT / "data/sandbox/BeamBeam-static-1V-PHI_scalar-beambeam_phi-000004.dat"
    e_path = ROOT / "data/sandbox/BeamBeam-static-1V-EF_vector-beambeam_e-000004.dat"
    for path in (rho_path, phi_path, e_path):
        if not path.exists():
            raise FileNotFoundError(
                f"missing BeamBeam field dump {path}; run sandbox/BeamBeam-static-1V.in first"
            )

    setup = manufactured.manufactured_setup_from_ascii(np, rho_path, phi_path, e_path)
    origin = setup["origin"]
    spacing = setup["spacing"]
    shape = setup["shape"]
    x, y, z = manufactured.build_grid(np, origin, spacing, shape)
    analytic = manufactured.gaussian_pair_fields(
        np, x, y, z, 1.0e-3, setup["charge"], setup["centers"]
    )

    metrics: dict[str, Any] = {
        "snapshot_kind": setup["snapshot_kind"],
        "charge_c": setup["charge"],
        "n_centers": len(setup["centers"]),
        "rho_integral_c": manufactured.compute_volume_integral(
            np, setup["opalx"]["rho"], spacing
        ),
    }
    for name in ("rho", "phi", "Ex", "Ey", "Ez"):
        item = manufactured.compute_error_metrics(np, analytic[name], setup["opalx"][name])
        metrics[name] = {
            "max_abs": item["max_abs"],
            "l2": item["l2"],
            "rel_l2": item["rel_l2"],
        }
    return metrics


def pair_metrics() -> dict[str, Any]:
    pair = load_module(
        "analyze_pair_momenta",
        ROOT / "sandbox/track-e-p/analyze_pair_momenta.py",
    )
    track_dir = ROOT / "sandbox/track-e-p"
    frames = {
        "electron": pair.add_momentum_columns(
            pair.read_fromfile(track_dir / "gamma_gamma_electrons.fromfile")
        ),
        "positron": pair.add_momentum_columns(
            pair.read_fromfile(track_dir / "gamma_gamma_positrons.fromfile")
        ),
    }
    species = pair.summarize_species(frames)
    species_wide = species.pivot(index="quantity", columns="species", values="value")
    pair_summary = pair.summarize_pairs(frames["electron"], frames["positron"])

    selected_species_quantities = [
        "N",
        "<px> beta_gamma",
        "<py> beta_gamma",
        "<pz> beta_gamma",
        "<theta_z> deg",
        "median theta_z deg",
        "<p_abs> beta_gamma",
        "<K> keV",
    ]
    selected_pair_quantities = [
        "<u_e dot u_p>",
        "median opening angle deg",
        "fraction opening > 150 deg",
        "fraction opening < 30 deg",
        "median |p_e+p_p| beta_gamma",
        "median |p_e+p_p|/sum|p|",
    ]

    metrics: dict[str, Any] = {
        "fromfile": {
            species_name: {
                quantity: float(species_wide.loc[quantity, species_name])
                for quantity in selected_species_quantities
            }
            for species_name in ("electron", "positron")
        },
        "pair_summary": {
            quantity: float(pair_summary[quantity]) for quantity in selected_pair_quantities
        },
        "track_stat": {},
    }

    for stem in ("gamma_gamma_pairs-2",):
        metrics["track_stat"][stem] = {}
        for container in ("c0", "c1", "c2"):
            stat_path = track_dir / f"{stem}_{container}.stat"
            df, _ = read_sdds_stat(stat_path)
            last = df.iloc[-1]
            metrics["track_stat"][stem][container] = {
                "rows": float(len(df)),
                "final_t_ns": float(last["t"]),
                "final_s_m": float(last["s"]),
                "final_numParticles": float(last["numParticles"]),
                "final_energy_MeV": float(last["energy"]),
                "final_rms_x_m": float(last["rms_x"]),
                "final_rms_y_m": float(last["rms_y"]),
                "final_mean_s_m": float(last["mean_s"]),
            }
    return metrics


def impact_metrics() -> dict[str, Any]:
    np = load_numpy()
    impact = load_module(
        "plot_extracted_graph",
        ROOT / "sandbox/OPALX-IMPACT/plot_extracted_graph.py",
    )
    impact_dir = ROOT / "sandbox/OPALX-IMPACT"
    digitized = pd.read_csv(impact_dir / "extracted_graph_data.csv")
    impact.validate_columns(digitized)
    opalx, _, _, _ = impact.read_opalx_stat(impact_dir / "drift-1.stat")

    t = opalx["t_ns"].to_numpy(float)
    impact_x = np.interp(t, digitized["t_ns"].to_numpy(float), digitized["ImpactT_x_m"].to_numpy(float))
    impact_eps = np.interp(
        t,
        digitized["t_ns"].to_numpy(float),
        digitized["ImpactT_epsx_mm_mr"].to_numpy(float),
    )
    dx = opalx["OPALX_x_m"].to_numpy(float) - impact_x
    deps = opalx["OPALX_epsx_mm_mr"].to_numpy(float) - impact_eps

    return {
        "digitized_rows": float(len(digitized)),
        "opalx_rows": float(len(opalx)),
        "rms_x_diff_m": float(np.sqrt(np.mean(dx * dx))),
        "max_abs_x_diff_m": float(np.max(np.abs(dx))),
        "rms_epsx_diff_mm_mr": float(np.sqrt(np.mean(deps * deps))),
        "max_abs_epsx_diff_mm_mr": float(np.max(np.abs(deps))),
        "final": {
            "t_ns": float(opalx["t_ns"].iloc[-1]),
            "OPALX_x_m": float(opalx["OPALX_x_m"].iloc[-1]),
            "OPALX_epsx_mm_mr": float(opalx["OPALX_epsx_mm_mr"].iloc[-1]),
        },
    }


CHECKS = {
    "fields": field_metrics,
    "pairs": pair_metrics,
    "impact": impact_metrics,
}


def git_metadata() -> dict[str, Any]:
    def git(*args: str) -> str:
        try:
            result = subprocess.run(
                ["git", *args],
                cwd=ROOT,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError as exc:
            return f"unavailable: {exc.stderr.strip() or exc}"

    status = git("status", "--short")
    return {
        "branch": git("branch", "--show-current"),
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status_short": status,
    }


def run_opalx_cases(opalx_exe: Path) -> None:
    if not opalx_exe.is_absolute():
        opalx_exe = ROOT / opalx_exe
    opalx_exe = opalx_exe.resolve()
    if not opalx_exe.exists():
        raise FileNotFoundError(opalx_exe)
    cases = [
        (ROOT, ["sandbox/BeamBeam-static-1V.in"]),
        (ROOT / "sandbox/OPALX-IMPACT", ["drift-1.in"]),
        (ROOT / "sandbox/track-e-p", ["gamma_gamma_pairs-2.in"]),
    ]
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    for cwd, input_args in cases:
        cmd = [str(opalx_exe), *input_args]
        print(f"running: {' '.join(cmd)}  (cwd={cwd})")
        for attempt in (1, 2):
            try:
                subprocess.run(cmd, cwd=cwd, env=env, check=True)
                break
            except subprocess.CalledProcessError:
                if attempt == 2:
                    raise
                print("  first OPALX launch failed; retrying once after MPI cleanup delay")
                time.sleep(2.0)


def compare_metrics(
        actual: dict[str, Any],
        baseline: dict[str, Any],
        rtol: float,
        atol: float,
) -> list[str]:
    failures: list[str] = []
    actual_flat = flatten("", actual)
    baseline_flat = flatten("", baseline)

    for key in sorted(set(actual_flat) | set(baseline_flat)):
        if key not in actual_flat:
            failures.append(f"{key}: missing from actual metrics")
            continue
        if key not in baseline_flat:
            failures.append(f"{key}: missing from baseline")
            continue

        actual_value = actual_flat[key]
        baseline_value = baseline_flat[key]
        if isinstance(actual_value, str) or isinstance(baseline_value, str):
            if actual_value != baseline_value:
                failures.append(f"{key}: {actual_value!r} != {baseline_value!r}")
            continue

        if not math.isclose(float(actual_value), float(baseline_value), rel_tol=rtol, abs_tol=atol):
            failures.append(
                f"{key}: actual={actual_value:.17g}, baseline={baseline_value:.17g}, "
                f"rtol={rtol:g}, atol={atol:g}"
            )
    return failures


def metrics_table(metrics: dict[str, Any]) -> pd.DataFrame:
    rows = [{"metric": key, "value": value} for key, value in flatten("", metrics).items()]
    return pd.DataFrame(rows).sort_values("metric")


def baseline_payload(metrics: dict[str, Any]) -> dict[str, Any]:
    return {
        "_metadata": {
            "git": git_metadata(),
            "checks": sorted(metrics),
        },
        **metrics,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        default=DEFAULT_BASELINE,
        help="JSON baseline file",
    )
    parser.add_argument(
        "--check",
        action="append",
        choices=sorted(CHECKS),
        help="Run only one check; may be repeated. Default: all checks.",
    )
    parser.add_argument(
        "--update-baseline",
        action="store_true",
        help="Write current metrics to the baseline instead of comparing.",
    )
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="CSV path for the current flattened metrics table. Default: %(default)s",
    )
    parser.add_argument(
        "--run-opalx",
        action="store_true",
        help="Run the sandbox OPALX inputs before checking metrics.",
    )
    parser.add_argument(
        "--opalx-exe",
        type=Path,
        default=ROOT / "build_openmp/src/opalx",
        help="OPALX executable used with --run-opalx.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "sandbox/data/.plot-cache/matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "sandbox/data/.plot-cache/xdg"))
    selected = args.check or sorted(CHECKS)

    if args.run_opalx:
        run_opalx_cases(args.opalx_exe)

    metrics = {name: CHECKS[name]() for name in selected}
    table = metrics_table(metrics)
    print(table.to_string(index=False))

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.csv, index=False)
    print(f"wrote current metrics: {args.csv}")

    if args.update_baseline:
        args.baseline.parent.mkdir(parents=True, exist_ok=True)
        args.baseline.write_text(json.dumps(baseline_payload(metrics), indent=2, sort_keys=True) + "\n")
        print(f"updated baseline: {args.baseline}")
        return 0

    if not args.baseline.exists():
        raise FileNotFoundError(f"missing baseline {args.baseline}; run with --update-baseline")

    baseline = json.loads(args.baseline.read_text())
    baseline = {name: baseline[name] for name in selected}
    failures = compare_metrics(metrics, baseline, args.rtol, args.atol)
    if failures:
        print("\nRegression mismatches:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1

    print("\nAll selected sandbox regression checks match the baseline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
