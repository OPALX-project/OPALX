#!/usr/bin/env python3
"""Generate the CAIN-only manufactured analytic solution.

This script does not read OPALX output.  It tracks the CAIN/TestParticleOrbit
test particles through the analytic boosted-Gaussian source model and writes
CSV data products.  Plot generation is intentionally left to make_plots.py.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_CAIN_CSV = ROOT / "sandbox" / "track12particles" / "outputs" / "track12_trajectory.csv"
DEFAULT_DATA_DIR = SCRIPT_DIR / "data"
DEFAULT_CONFIG = SCRIPT_DIR / "manufactured_sol.json"
DEFAULT_SOURCE_SIGMA_Z_M = 0.6e-3


def resolve_config_path(config_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return config_dir / path


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config_dir = path.resolve().parent
    return {
        "config_path": path.resolve(),
        "cain_csv": resolve_config_path(config_dir, config.get("cain_csv")) or DEFAULT_CAIN_CSV,
        "data_dir": resolve_config_path(config_dir, config.get("data_dir")) or DEFAULT_DATA_DIR,
        "trajectory_csv": resolve_config_path(config_dir, config.get("trajectory_csv")),
        "summary_csv": resolve_config_path(config_dir, config.get("summary_csv")),
        "model": config.get("model", "anisotropic"),
        "charge_scale": float(config.get("charge_scale", 0.5)),
        "source_selection": config.get("source_selection", "oncoming"),
        "source_time_offset_ps": float(config.get("source_time_offset_ps", 0.0)),
        "source_sigma_z_m": float(config.get("source_sigma_z_m", DEFAULT_SOURCE_SIGMA_Z_M)),
        "mc_source_particles": int(config.get("mc_source_particles", 400000)),
        "max_substep_s": float(config.get("max_substep_s", 5.0e-15)),
    }


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    config["data_dir"].mkdir(parents=True, exist_ok=True)

    trajectory_csv = config["trajectory_csv"] or config["data_dir"] / "cain-fit-half-charge-analytic-trajectory.csv"
    summary_csv = config["summary_csv"] or config["data_dir"] / "cain-fit-half-charge-summary.csv"

    track12 = load_track12_module()
    reference = pd.read_csv(config["cain_csv"])
    source_charge_C = track12.DEFAULT_SOURCE_CHARGE_C * config["charge_scale"]
    manufactured = track12.track_manufactured_on_reference_grid(
        reference,
        config["model"],
        track12.DEFAULT_SOURCE_KINETIC_MEV,
        source_charge_C,
        config["source_sigma_z_m"],
        track12.SIGMA_XY_M,
        track12.SIGMA_XY_M,
        config["source_sigma_z_m"],
        config["max_substep_s"],
        config["source_selection"],
        config["source_time_offset_ps"] * 1.0e-12,
        config["mc_source_particles"],
        track12.DEFAULT_MC_SEED,
    )
    metrics = track12.trajectory_error_metrics(reference, manufactured)
    summary = pd.DataFrame(
        [
            {
                "config": str(config["config_path"]),
                "model": config["model"],
                "source_selection": config["source_selection"],
                "source_time_offset_ps": config["source_time_offset_ps"],
                "charge_scale": config["charge_scale"],
                "source_charge_C": source_charge_C,
                "source_sigma_z_m": config["source_sigma_z_m"],
                "mc_source_particles": config["mc_source_particles"],
                **metrics,
            }
        ]
    )

    trajectory_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    manufactured.to_csv(trajectory_csv, index=False)
    summary.to_csv(summary_csv, index=False)

    print(summary.to_string(index=False))
    print(f"Wrote {summary_csv}")
    print(f"Wrote {trajectory_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
