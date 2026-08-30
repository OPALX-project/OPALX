#!/usr/bin/env python3
"""Prepare the reproducible BeamBeam A100 convergence case matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
ROOT = MODEL_DIR.parents[1]
GENERATOR_SCRIPT = SCRIPT_DIR / "make_opalx_field_cases.py"
DEFAULT_STUDY = MODEL_DIR / "a100_convergence_study.json"


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--study-config", type=Path, default=DEFAULT_STUDY)
    parser.add_argument("--opalx", type=Path, required=True)
    return parser.parse_args()


def add_case(
    cases: dict[tuple[int, int, int, int, str], set[str]],
    particles: int,
    nxy: int,
    seed: int,
    ranks: int,
    scan: str,
    primary_sampling: str = "gauss",
) -> None:
    if particles <= 0 or particles % 1000 != 0:
        raise ValueError(f"particle count must be a positive multiple of 1000: {particles}")
    if nxy <= 0 or ranks <= 0:
        raise ValueError(f"mesh and rank counts must be positive: nxy={nxy}, ranks={ranks}")
    if primary_sampling not in {"gauss", "fixed_fromfile"}:
        raise ValueError(f"unknown primary sampling mode: {primary_sampling}")
    cases.setdefault((particles, nxy, seed, ranks, primary_sampling), set()).add(scan)


def case_label(
    particles: int, nxy: int, seed: int, ranks: int, primary_sampling: str = "gauss"
) -> str:
    label = f"N{particles // 1000}k_M{nxy}_S{seed}"
    if primary_sampling == "fixed_fromfile":
        return f"{label}_FIXED_MPI{ranks}"
    return label if ranks == 1 else f"{label}_MPI{ranks}"


def git_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    return completed.stdout.strip()


def sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    with args.study_config.open("r", encoding="utf-8") as stream:
        study = json.load(stream)

    base_seed = int(study["base_seed"])
    fixed_nxy = int(study["fixed_mesh_nxy"])
    fixed_particles = int(study["fixed_particles"])
    cases: dict[tuple[int, int, int, int, str], set[str]] = {}

    for particles in study["particle_scan"]:
        add_case(cases, int(particles), fixed_nxy, base_seed, 1, "particle")
    for nxy in study["fixed_particle_mesh_scan"]:
        add_case(cases, fixed_particles, int(nxy), base_seed, 1, "fixed_particle_mesh")
    for entry in study["constant_occupancy_scan"]:
        add_case(
            cases,
            int(entry["primary_macroparticles"]),
            int(entry["mesh_nxy"]),
            base_seed,
            1,
            "constant_occupancy",
        )
    for seed in study["seed_scan"]:
        add_case(cases, fixed_particles, fixed_nxy, int(seed), 1, "seed")
    for ranks in study["mpi_rank_scan"]:
        add_case(
            cases,
            fixed_particles,
            fixed_nxy,
            base_seed,
            int(ranks),
            "mpi_rank",
            "fixed_fromfile",
        )

    generator = load_script("beambeam_field_generator", GENERATOR_SCRIPT)
    field_model = generator.load_field_model()
    parameters = field_model.load_parameters(field_model.DEFAULT_CONFIG)
    parameters = replace(
        parameters,
        nx=int(study["probe_nx"]),
        nz=int(study["probe_nz"]),
    )
    field_model.validate_parameters(parameters)
    track12 = field_model.load_track12_module()
    template = generator.DEFAULT_TEMPLATE.read_text(encoding="utf-8")

    fixed_source_config = study["fixed_rank_source"]
    fixed_source_path = args.run_root / (
        f"primary_N{fixed_particles}_S{base_seed}_fixed.fromfile"
    )
    fixed_source_metadata = generator.write_deterministic_primary_file(
        fixed_source_path,
        parameters,
        track12,
        fixed_particles,
        base_seed,
        cutoff_sigma=float(fixed_source_config["position_cutoff_sigma"]),
        sigma_xyp=float(fixed_source_config["sigma_xyp"]),
        sigma_pz=float(fixed_source_config["sigma_pz"]),
    )
    fixed_source_metadata["path"] = fixed_source_path.name
    (args.run_root / "fixed_primary_metadata.json").write_text(
        json.dumps(fixed_source_metadata, indent=2) + "\n", encoding="utf-8"
    )

    args.run_root.mkdir(parents=True, exist_ok=True)
    case_root = args.run_root / "cases"
    case_root.mkdir(exist_ok=True)
    rows = []
    for (particles, nxy, seed, ranks, primary_sampling), scans in sorted(cases.items()):
        label = case_label(particles, nxy, seed, ranks, primary_sampling)
        case_dir = case_root / label / "3sigma"
        case_dir.mkdir(parents=True, exist_ok=True)
        primary_particle_filename = None
        primary_file_sha256 = None
        if primary_sampling == "fixed_fromfile":
            primary_particle_filename = os.path.relpath(fixed_source_path, case_dir)
            primary_file_sha256 = fixed_source_metadata["sha256"]
        input_path = generator.render_case(
            template,
            parameters,
            track12,
            float(study["separation_sigma_z"]),
            case_dir,
            particles,
            nxy,
            2 * nxy,
            seed,
            primary_particle_filename=primary_particle_filename,
        )
        rows.append(
            {
                "label": label,
                "primary_macroparticles": particles,
                "mesh_nxy": nxy,
                "mesh_nz": 2 * nxy,
                "seed": seed,
                "mpi_ranks": ranks,
                "primary_sampling": primary_sampling,
                "primary_file_sha256": primary_file_sha256,
                "scans": ";".join(sorted(scans)),
                "case_directory": str(case_dir.relative_to(args.run_root)),
                "input_file": input_path.name,
            }
        )

    columns = list(rows[0])
    with (args.run_root / "study_manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    (args.run_root / "study_manifest.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    for ranks in (1, 2, 4):
        labels = [row["label"] for row in rows if row["mpi_ranks"] == ranks]
        (args.run_root / f"rank{ranks}_labels.txt").write_text(
            "".join(f"{label}\n" for label in labels), encoding="utf-8"
        )
        mpi_labels = [
            row["label"]
            for row in rows
            if row["mpi_ranks"] == ranks and "mpi_rank" in row["scans"].split(";")
        ]
        (args.run_root / f"mpi_rank{ranks}_labels.txt").write_text(
            "".join(f"{label}\n" for label in mpi_labels), encoding="utf-8"
        )
    baseline_label = case_label(
        fixed_particles, fixed_nxy, base_seed, 1, "fixed_fromfile"
    )
    (args.run_root / "rank_baseline_label.txt").write_text(
        f"{baseline_label}\n", encoding="utf-8"
    )

    source_patch = subprocess.run(
        ["git", "diff", "--no-ext-diff", "--binary", "HEAD", "--", "src"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    (args.run_root / "opalx_source.patch").write_bytes(source_patch)
    source_status = subprocess.run(
        ["git", "status", "--short", "--untracked-files=all"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.splitlines()

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository": str(ROOT),
        "git_revision": git_revision(),
        "source_worktree_status": source_status,
        "opalx_source_patch_sha256": hashlib.sha256(source_patch).hexdigest(),
        "opalx_executable": str(args.opalx),
        "opalx_executable_sha256": sha256_file(args.opalx),
        "study_config": study,
        "fixed_primary_source": fixed_source_metadata,
        "case_count": len(rows),
        "physics": {
            "field_step": 1,
            "source_model": "physical-and-copied",
            "tracking_steps": 2,
            "probe_grid": [parameters.nx, parameters.nz],
            "pair_fields_included": False,
            "mpi_rank_scan_primary_sampling": "fixed_fromfile",
        },
    }
    (args.run_root / "study_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"prepared {len(rows)} cases below {case_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
