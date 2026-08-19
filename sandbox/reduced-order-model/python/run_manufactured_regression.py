#!/usr/bin/env python3
"""Run and check the persistent 3-sigma BeamBeam field regression."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import shutil
import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
ROOT = MODEL_DIR.parents[1]
GENERATOR_SCRIPT = SCRIPT_DIR / "make_opalx_field_cases.py"
COMPARATOR_SCRIPT = SCRIPT_DIR / "compare_opalx_fields.py"
DEFAULT_CONFIG = MODEL_DIR / "regression.json"
DEFAULT_OPALX = ROOT / "build_openmp" / "src" / "opalx"


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_regression_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        config = json.load(stream)
    case = config["case"]
    if case["probe_nx"] % 2 == 0 or case["probe_nz"] % 2 == 0:
        raise ValueError("regression probe dimensions must be odd to sample the IP")
    return config


def check_summary(summary: dict, config: dict) -> list[str]:
    acceptance = config["acceptance"]
    failures = []

    upper_bounds = (
        ("relative_l2_E", "relative_l2_E_max"),
        ("relative_l2_B", "relative_l2_B_max"),
        ("uncovered_probe_fraction", "uncovered_probe_fraction_max"),
    )
    for metric, limit_name in upper_bounds:
        value = summary[metric]
        if not math.isfinite(value):
            failures.append(f"{metric} is not finite: {value}")
        elif value > acceptance[limit_name]:
            failures.append(
                f"{metric}={value:.8g} exceeds "
                f"{acceptance[limit_name]:.8g}"
            )

    bounded_ratios = (
        (
            "median_E_magnitude_ratio",
            "median_E_magnitude_ratio_min",
            "median_E_magnitude_ratio_max",
        ),
        (
            "median_B_magnitude_ratio",
            "median_B_magnitude_ratio_min",
            "median_B_magnitude_ratio_max",
        ),
    )
    for metric, minimum_name, maximum_name in bounded_ratios:
        value = summary[metric]
        if not math.isfinite(value):
            failures.append(f"{metric} is not finite: {value}")
        elif not acceptance[minimum_name] < value < acceptance[maximum_name]:
            failures.append(
                f"{metric}={value:.8g} is outside "
                f"({acceptance[minimum_name]:.8g}, {acceptance[maximum_name]:.8g})"
            )

    lower_bounds = (
        ("median_E_direction_cosine", "median_E_direction_cosine_min"),
        ("median_B_direction_cosine", "median_B_direction_cosine_min"),
    )
    for metric, limit_name in lower_bounds:
        value = summary[metric]
        if not math.isfinite(value):
            failures.append(f"{metric} is not finite: {value}")
        elif value <= acceptance[limit_name]:
            failures.append(
                f"{metric}={value:.8g} does not exceed "
                f"{acceptance[limit_name]:.8g}"
            )
    return failures


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--regression-config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    opalx = args.opalx.resolve()
    if not opalx.is_file():
        raise FileNotFoundError(opalx)

    config = load_regression_config(args.regression_config)
    case_config = config["case"]
    generator = load_script("beambeam_field_generator", GENERATOR_SCRIPT)
    comparator = load_script("beambeam_field_comparator", COMPARATOR_SCRIPT)
    field_model = generator.load_field_model()
    parameters = field_model.load_parameters(field_model.DEFAULT_CONFIG)
    parameters = replace(
        parameters,
        nx=int(case_config["probe_nx"]),
        nz=int(case_config["probe_nz"]),
    )
    field_model.validate_parameters(parameters)
    track12 = field_model.load_track12_module()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    case_dir = args.work_dir / "3sigma"
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir()

    template = generator.DEFAULT_TEMPLATE.read_text(encoding="utf-8")
    input_path = generator.render_case(
        template,
        parameters,
        track12,
        float(case_config["separation_sigma_z"]),
        case_dir,
        int(case_config["primary_macroparticles"]),
        int(case_config["mesh_nxy"]),
        int(case_config["mesh_nz"]),
    )
    generator.run_case(
        opalx,
        input_path,
        int(case_config["omp_threads"]),
        force=True,
    )

    source_h5 = list(case_dir.glob("*_c0.h5"))
    if len(source_h5) > 1:
        raise RuntimeError(f"expected at most one source H5 shell, found {source_h5}")
    if source_h5:
        with comparator.h5py.File(source_h5[0], "r") as source_file:
            source_steps = list(source_file.keys())
            if source_steps:
                raise RuntimeError(
                    f"source particle H5 output must be disabled, found steps "
                    f"{source_steps}"
                )

    _, summary = comparator.analyze_case(
        case_dir,
        field_model,
        track12,
        parameters,
        step_index=1,
        source_model="physical-and-copied",
    )
    expected_probes = int(case_config["probe_nx"]) * int(case_config["probe_nz"])
    if summary["field_probe_count"] != expected_probes:
        raise AssertionError(
            f"expected {expected_probes} probes, found {summary['field_probe_count']}"
        )

    metric_names = (
        "relative_l2_E",
        "relative_l2_B",
        "median_E_magnitude_ratio",
        "median_B_magnitude_ratio",
        "median_E_direction_cosine",
        "median_B_direction_cosine",
        "uncovered_probe_fraction",
    )
    metric_frame = pd.DataFrame(
        {"metric": metric_names, "value": [summary[name] for name in metric_names]}
    )
    metric_frame.to_csv(args.work_dir / "manufactured_regression_summary.csv", index=False)
    print(metric_frame.to_string(index=False, float_format=lambda value: f"{value:.8e}"))

    failures = check_summary(summary, config)
    if failures:
        print("BeamBeam manufactured-field regression failed:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print("BeamBeam manufactured-field regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
