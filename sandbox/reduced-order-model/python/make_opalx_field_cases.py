#!/usr/bin/env python3
"""Generate and optionally run OPALX rigid-source field-snapshot cases."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
ROOT = MODEL_DIR.parents[1]
FIELD_MODEL_SCRIPT = SCRIPT_DIR / "rigid_two_gaussian_fields.py"
DEFAULT_TEMPLATE = MODEL_DIR / "opalx" / "rigid_two_gaussian_field_snapshot.in.template"
DEFAULT_OUTPUT_DIR = MODEL_DIR / "outputs" / "opalx-fields"
DEFAULT_OPALX = ROOT / "build_openmp" / "src" / "opalx"


def load_field_model():
    spec = importlib.util.spec_from_file_location("rigid_two_gaussian_fields", FIELD_MODEL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {FIELD_MODEL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_probe_file(path: Path, parameters) -> int:
    x_values = np.linspace(
        -parameters.x_half_width_sigma_x * parameters.sigma_x_m,
        +parameters.x_half_width_sigma_x * parameters.sigma_x_m,
        parameters.nx,
    )
    z_values = np.linspace(
        -parameters.z_half_width_sigma_z * parameters.sigma_z_m,
        +parameters.z_half_width_sigma_z * parameters.sigma_z_m,
        parameters.nz,
    )
    z_grid, x_grid = np.meshgrid(z_values, x_values, indexing="ij")
    count = x_grid.size
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{count}\n")
        stream.write("x y z px py pz\n")
        for x_value, z_value in zip(x_grid.ravel(), z_grid.ravel()):
            stream.write(
                f"{x_value:.16e} {parameters.plane_y_m:.16e} {z_value:.16e} "
                "0.0000000000000000e+00 0.0000000000000000e+00 "
                "1.7320508075688772e+00\n"
            )
    return count


def _truncated_standard_normal(
    rng: np.random.Generator, shape: tuple[int, ...], cutoff_sigma: float
) -> np.ndarray:
    """Draw independent standard normals conditioned on ``abs(x) <= cutoff``."""
    if cutoff_sigma <= 0.0:
        raise ValueError("cutoff_sigma must be positive")
    values = np.empty(int(np.prod(shape)), dtype=np.float64)
    remaining = np.arange(values.size)
    while remaining.size:
        candidates = rng.standard_normal(remaining.size)
        accepted = np.abs(candidates) <= cutoff_sigma
        values[remaining[accepted]] = candidates[accepted]
        remaining = remaining[~accepted]
    return values.reshape(shape)


def write_deterministic_primary_file(
    path: Path,
    parameters,
    track12,
    particle_count: int,
    seed: int,
    *,
    cutoff_sigma: float = 3.0,
    sigma_xyp: float = 1.0e-12,
    sigma_pz: float = 1.0e-12,
) -> dict[str, object]:
    """Write one fixed primary Gaussian sample for MPI-decomposition tests.

    The position construction matches the active OPALX ``GAUSS`` semantics:
    independent normals truncated at +/- ``cutoff_sigma`` and then translated
    so that the finite sample has exactly zero centroid. File momenta are
    absolute normalized momenta (beta*gamma), as required by ``FROMFILE``.
    """
    if particle_count <= 0:
        raise ValueError("particle_count must be positive")
    if sigma_xyp < 0.0 or sigma_pz < 0.0:
        raise ValueError("momentum widths must be non-negative")

    rng = np.random.Generator(np.random.PCG64(seed))
    normalized_positions = _truncated_standard_normal(
        rng, (particle_count, 3), cutoff_sigma
    )
    position_sigmas = np.array(
        [parameters.sigma_x_m, parameters.sigma_y_m, parameters.sigma_z_m],
        dtype=np.float64,
    )
    positions = normalized_positions * position_sigmas
    positions -= np.mean(positions, axis=0, dtype=np.float64)

    gamma = 1.0 + parameters.kinetic_energy_MeV * 1.0e6 / track12.ELECTRON_REST_EV
    beta_gamma = np.sqrt(gamma * gamma - 1.0)
    sigma_pxy = beta_gamma * sigma_xyp
    momenta = np.empty_like(positions)
    momenta[:, 0] = rng.normal(0.0, sigma_pxy, particle_count)
    momenta[:, 1] = rng.normal(0.0, sigma_pxy, particle_count)
    momenta[:, 2] = rng.normal(beta_gamma, sigma_pz, particle_count)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{particle_count}\n")
        stream.write("x y z px py pz\n")
        np.savetxt(stream, np.column_stack((positions, momenta)), fmt="%.16e")

    digest_state = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest_state.update(chunk)
    return {
        "path": str(path),
        "sha256": digest_state.hexdigest(),
        "particle_count": particle_count,
        "seed": seed,
        "rng": "NumPy PCG64",
        "numpy_version": np.__version__,
        "position_cutoff_sigma_before_recentering": cutoff_sigma,
        "position_centroid_m": np.mean(positions, axis=0).tolist(),
        "sigma_xyp": sigma_xyp,
        "sigma_pz": sigma_pz,
        "beta_gamma": float(beta_gamma),
    }


def render_case(
    template: str,
    parameters,
    track12,
    separation: float,
    case_dir: Path,
    primary_macroparticles: int,
    nxy: int,
    nz: int,
    seed: int = 20260629,
    primary_particle_filename: str | None = None,
) -> Path:
    label = f"{separation:g}".replace(".", "p") + "sigma"
    case_name = f"rigid_fields_{label}"
    probe_name = f"{case_name}_probes.fromfile"
    probe_count = write_probe_file(case_dir / probe_name, parameters)

    beta = track12.beta_from_kinetic_energy(parameters.kinetic_energy_MeV)
    dt_s = 1.0e-15
    ip_s_m = 0.008
    # The copied source first participates in the second solve. Its physical
    # source is advanced by approximately 1.5*beta*c*dt before that midpoint
    # field evaluation; compensate this negligible offset explicitly.
    source_r0z_m = (
        ip_s_m
        - 0.5 * separation * parameters.sigma_z_m
        - 1.5 * beta * track12.C_LIGHT * dt_s
    )
    if primary_particle_filename is None:
        primary_distribution = """DIST_PRIMARY_ELECTRONS: DISTRIBUTION,
    TYPE = GAUSS,
    SIGMAX = primary_sigma_xy,
    SIGMAY = primary_sigma_xy,
    SIGMAZ = primary_sigma_z,
    SIGMAPX = primary_sigma_pxy,
    SIGMAPY = primary_sigma_pxy,
    SIGMAPZ = 1.0e-12,
    CUTOFFX = 6.0,
    CUTOFFY = 6.0,
    CUTOFFLONG = 6.0,
    NPARTDIST = primary_macroparticles;"""
        primary_beam_energy = "    PC = primary_p0,"
    else:
        primary_distribution = f"""DIST_PRIMARY_ELECTRONS: DISTRIBUTION,
    TYPE = FROMFILE,
    FNAME = \"{primary_particle_filename}\",
    NPARTDIST = primary_macroparticles;"""
        # FROMFILE contains absolute beta*gamma and rejects PC/ENERGY/GAMMA.
        primary_beam_energy = ""

    replacements = {
        "@CASE_TITLE@": f"Rigid two-Gaussian field snapshot, d={separation:g} sigma_z",
        "@SEPARATION_SIGMA_Z@": f"{separation:.16e}",
        "@N_FIELD_PROBES@": str(probe_count),
        "@PRIMARY_SOURCE_R0Z_M@": f"{source_r0z_m:.16e}",
        "@FIELD_PROBE_FILENAME@": probe_name,
        "@PRIMARY_MACROPARTICLES@": str(primary_macroparticles),
        "@NXY@": str(nxy),
        "@NZ@": str(nz),
        "@SEED@": str(seed),
        "@PRIMARY_DISTRIBUTION@": primary_distribution,
        "@PRIMARY_BEAM_ENERGY@": primary_beam_energy,
    }
    rendered = template
    for token, value in replacements.items():
        rendered = rendered.replace(token, value)
    unresolved = [token for token in replacements if token in rendered]
    if unresolved:
        raise RuntimeError(f"unresolved template tokens: {unresolved}")

    input_path = case_dir / f"{case_name}.in"
    input_path.write_text(rendered, encoding="utf-8")
    return input_path


def run_case(opalx: Path, input_path: Path, omp_threads: int, force: bool) -> None:
    log_path = input_path.with_suffix(".log")
    witness_h5 = input_path.with_name(input_path.stem + "_c1.h5")
    if witness_h5.exists() and not force:
        print(f"cached {witness_h5}")
        return
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(omp_threads)
    completed = subprocess.run(
        [str(opalx), input_path.name],
        cwd=input_path.parent,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"OPALX failed for {input_path.name} with code {completed.returncode}; see {log_path}"
        )
    print(f"ran {input_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--omp-threads", type=int, default=4)
    parser.add_argument("--primary-macroparticles", type=int, default=400000)
    parser.add_argument("--nxy", type=int, default=64)
    parser.add_argument("--nz", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260629)
    parser.add_argument(
        "--primary-particle-file",
        type=str,
        default=None,
        help="Use a FROMFILE primary sample instead of OPALX GAUSS sampling.",
    )
    parser.add_argument(
        "--probe-nx",
        type=int,
        default=None,
        help="Override the configured number of passive probes along x.",
    )
    parser.add_argument(
        "--probe-nz",
        type=int,
        default=None,
        help="Override the configured number of passive probes along z.",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--separations-sigma-z",
        type=float,
        nargs="+",
        default=None,
        help="Optional subset/override of configured centroid separations.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.primary_macroparticles <= 0:
        raise ValueError("--primary-macroparticles must be positive")
    if args.nxy <= 0 or args.nz <= 0:
        raise ValueError("--nxy and --nz must be positive")
    field_model = load_field_model()
    config = args.config or field_model.DEFAULT_CONFIG
    parameters = field_model.load_parameters(config)
    if args.probe_nx is not None or args.probe_nz is not None:
        parameters = replace(
            parameters,
            nx=args.probe_nx if args.probe_nx is not None else parameters.nx,
            nz=args.probe_nz if args.probe_nz is not None else parameters.nz,
        )
        field_model.validate_parameters(parameters)
    track12 = field_model.load_track12_module()
    template = args.template.read_text(encoding="utf-8")
    separations = args.separations_sigma_z or parameters.centroid_separations_sigma_z

    if args.run and not args.opalx.exists():
        raise FileNotFoundError(args.opalx)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for separation in separations:
        label = f"{separation:g}".replace(".", "p") + "sigma"
        case_dir = args.output_dir / label
        case_dir.mkdir(parents=True, exist_ok=True)
        input_path = render_case(
            template,
            parameters,
            track12,
            separation,
            case_dir,
            args.primary_macroparticles,
            args.nxy,
            args.nz,
            args.seed,
            args.primary_particle_file,
        )
        print(f"wrote {input_path}")
        if args.run:
            run_case(args.opalx, input_path, args.omp_threads, args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
