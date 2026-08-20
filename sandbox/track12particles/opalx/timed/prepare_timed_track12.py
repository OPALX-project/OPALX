#!/usr/bin/env python3
"""Prepare exact-timed OPALX inputs for the 12 CAIN test particles."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REFERENCE = ROOT / "sandbox" / "TestParticleOrbit.dat"
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
TEMPLATE = HERE / "track12_timed.in.template"
ELECTRON_REST_EV = 510_998.95
C_LIGHT = 299_792_458.0
CAIN_CT_STEP_M = 1.8e-6
SEED = 20260629


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def truncated_standard_normal(
    rng: np.random.Generator, shape: tuple[int, ...], cutoff_sigma: float
) -> np.ndarray:
    values = np.empty(int(np.prod(shape)), dtype=np.float64)
    remaining = np.arange(values.size)
    while remaining.size:
        candidates = rng.standard_normal(remaining.size)
        accepted = np.abs(candidates) <= cutoff_sigma
        values[remaining[accepted]] = candidates[accepted]
        remaining = remaining[~accepted]
    return values.reshape(shape)


def write_fixed_primary(path: Path, particle_count: int) -> dict[str, object]:
    rng = np.random.Generator(np.random.PCG64(SEED))
    sigmas = np.array([1.944325075701e-6, 1.944325075701e-6, 0.6e-3])
    positions = truncated_standard_normal(rng, (particle_count, 3), 3.0) * sigmas
    positions -= np.mean(positions, axis=0, dtype=np.float64)

    gamma = 1.0 + 245.0e6 / ELECTRON_REST_EV
    beta_gamma = np.sqrt(gamma * gamma - 1.0)
    momenta = np.empty_like(positions)
    momenta[:, 0] = rng.normal(0.0, beta_gamma * 1.0e-12, particle_count)
    momenta[:, 1] = rng.normal(0.0, beta_gamma * 1.0e-12, particle_count)
    momenta[:, 2] = rng.normal(beta_gamma, 1.0e-12, particle_count)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{particle_count}\n")
        stream.write("x y z px py pz\n")
        np.savetxt(stream, np.column_stack((positions, momenta)), fmt="%.16e")

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": str(path.resolve()),
        "sha256": digest,
        "particle_count": particle_count,
        "seed": SEED,
        "rng": "NumPy PCG64",
        "numpy_version": np.__version__,
        "position_cutoff_sigma_before_recentering": 3.0,
        "position_centroid_m": np.mean(positions, axis=0).tolist(),
        "beta_gamma": float(beta_gamma),
    }


def write_witness(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(rows)}\n")
        stream.write("x y z px py pz birth_time\n")
        for row in rows.itertuples(index=False):
            stream.write(
                f"{row.x:.16e} {row.y:.16e} {row.s:.16e} "
                f"{row.Px / ELECTRON_REST_EV:.16e} "
                f"{row.Py / ELECTRON_REST_EV:.16e} "
                f"{row.Ps / ELECTRON_REST_EV:.16e} "
                f"{row.t / C_LIGHT:.16e}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--primary-macroparticles", type=int, default=400_000)
    parser.add_argument("--nx", type=int, default=1024)
    parser.add_argument("--ny", type=int, default=128)
    parser.add_argument("--nz", type=int, default=128)
    args = parser.parse_args()
    if args.primary_macroparticles <= 0:
        raise ValueError("--primary-macroparticles must be positive")

    module = load_track12_module()
    reference = module.parse_reference_file(args.reference)
    initial = (
        reference.sort_values(["pair", "kind", "step"])
        .groupby("name", sort=False, as_index=False)
        .first()
        .sort_values(["kind", "pair"])
    )
    electrons = initial.loc[initial["kind"] == 2].sort_values("pair")
    positrons = initial.loc[initial["kind"] == 3].sort_values("pair")
    if len(electrons) != 6 or len(positrons) != 6:
        raise ValueError("expected exactly six electron and six positron initial rows")
    np.testing.assert_allclose(electrons["t"], positrons["t"], rtol=0.0, atol=0.0)

    input_dir = HERE / "input"
    write_witness(input_dir / "track12_electrons.emittedfromfile", electrons)
    write_witness(input_dir / "track12_positrons.emittedfromfile", positrons)
    primary_metadata = write_fixed_primary(
        input_dir / "primary_fixed.fromfile", args.primary_macroparticles
    )
    (input_dir / "primary_fixed_metadata.json").write_text(
        json.dumps(primary_metadata, indent=2) + "\n", encoding="utf-8"
    )

    deck = TEMPLATE.read_text(encoding="utf-8")
    replacements = {
        "@PRIMARY_MACROPARTICLES@": str(args.primary_macroparticles),
        "@NX@": str(args.nx),
        "@NY@": str(args.ny),
        "@NZ@": str(args.nz),
        "@SEED@": str(SEED),
    }
    for token, value in replacements.items():
        deck = deck.replace(token, value)
    if "@" in deck:
        raise ValueError("unreplaced token remains in generated deck")
    (HERE / "track12_timed.in").write_text(deck, encoding="utf-8")

    birth_steps = np.rint(
        1.0 + (initial.sort_values(["kind", "pair"])["t"].to_numpy() - initial["t"].min())
        / CAIN_CT_STEP_M
    ).astype(int)
    manifest = {
        "reference": str(args.reference.resolve()),
        "reference_rows": int(len(reference)),
        "witnesses_per_species": 6,
        "cain_ct_step_m": CAIN_CT_STEP_M,
        "track_dt_s": CAIN_CT_STEP_M / C_LIGHT,
        "pair_t0_s": float(-initial["t"].min() / C_LIGHT + CAIN_CT_STEP_M / C_LIGHT),
        "birth_global_steps_by_species_then_pair": birth_steps.tolist(),
        "final_global_step": 1500,
        "final_cain_ct_m": 1.8e-3,
        "primary": primary_metadata,
        "mesh": [args.nx, args.ny, args.nz],
    }
    results_dir = HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "preparation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    print(f"deck: {HERE / 'track12_timed.in'}")


if __name__ == "__main__":
    main()
