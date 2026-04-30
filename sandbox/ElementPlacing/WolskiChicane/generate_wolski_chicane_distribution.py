#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


N_PARTICLES = 5000
RNG_SEED = 1234
BUNCH_CHARGE_C = 8.010883169999999e-16

ELECTRON_MASS_MEV = 0.51099895
KINETIC_ENERGY_MEV = 5000.0
SIGMA_Z_M = 6.0e-3
SIGMA_DELTA_UNCORRELATED = 1.3e-3
R56_M = -0.230480589
H_CHIRP_PER_M = 4.327913277


def build_sample(rng: np.random.Generator) -> dict[str, np.ndarray]:
    gamma = (KINETIC_ENERGY_MEV + ELECTRON_MASS_MEV) / ELECTRON_MASS_MEV
    beta_gamma = math.sqrt(gamma * gamma - 1.0)

    z = rng.normal(loc=0.0, scale=SIGMA_Z_M, size=N_PARTICLES)
    delta_uncorrelated = rng.normal(
        loc=0.0, scale=SIGMA_DELTA_UNCORRELATED, size=N_PARTICLES
    )
    return {
        "x": np.zeros(N_PARTICLES),
        "px": np.zeros(N_PARTICLES),
        "y": np.zeros(N_PARTICLES),
        "py": np.zeros(N_PARTICLES),
        "z": z,
        "delta_uncorrelated": delta_uncorrelated,
        "beta_gamma_ref": np.full(N_PARTICLES, beta_gamma),
    }


def build_opalx_distribution(
    sample: dict[str, np.ndarray], chirp_sign: float
) -> dict[str, np.ndarray]:
    delta = chirp_sign * H_CHIRP_PER_M * sample["z"] + sample["delta_uncorrelated"]
    return {
        "x": sample["x"],
        "px": sample["px"],
        "y": sample["y"],
        "py": sample["py"],
        "z": sample["z"],
        "pz": sample["beta_gamma_ref"] * (1.0 + delta),
    }


def build_bmad_distribution(
    sample: dict[str, np.ndarray], chirp_sign: float | None
) -> dict[str, np.ndarray]:
    """Build a BMAD particle distribution in canonical coordinates.

    The BMAD beam files use ``pz = delta = (P - P0) / P0``. For the RF-driven
    benchmark, the input beam is deliberately left unchirped and the cavity
    generates the chirp internally. For the no-RF benchmarks, the chirp is
    embedded directly in the input beam via

    .. math::

        \delta(z) = h z + \delta_u.
    """
    if chirp_sign is None:
        delta = sample["delta_uncorrelated"]
    else:
        delta = chirp_sign * H_CHIRP_PER_M * sample["z"] + sample["delta_uncorrelated"]
    return {
        "x": sample["x"],
        "px": sample["px"],
        "y": sample["y"],
        "py": sample["py"],
        "z": sample["z"],
        "pz": delta,
    }


def column_stack_distribution(df: dict[str, np.ndarray]) -> np.ndarray:
    return np.column_stack((df["x"], df["px"], df["y"], df["py"], df["z"], df["pz"]))


def write_opalx_distribution(path: Path, df: dict[str, np.ndarray]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"{N_PARTICLES}\n")
        fh.write("x px y py z pz\n")
        np.savetxt(fh, column_stack_distribution(df), fmt="%.16e")


def write_bmad_distribution(
    path: Path, df: dict[str, np.ndarray], reference_p0c_ev: float
) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("# species = electron\n")
        fh.write(f"# charge_tot = {BUNCH_CHARGE_C:.16e}\n")
        fh.write(f"# p0c = {reference_p0c_ev:.16e}\n")
        fh.write("# state = alive\n")
        fh.write("# location = upstream_end\n")
        fh.write("# direction = 1\n")
        fh.write("# time_dir = 1\n")
        fh.write("#! x px y py z pz\n")
        np.savetxt(fh, column_stack_distribution(df), fmt="%.16e")


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    sample = build_sample(rng)

    opalx_pos = build_opalx_distribution(sample, chirp_sign=+1.0)
    opalx_neg = build_opalx_distribution(sample, chirp_sign=-1.0)
    bmad_unchirped = build_bmad_distribution(sample, chirp_sign=None)
    bmad_pos = build_bmad_distribution(sample, chirp_sign=+1.0)
    bmad_neg = build_bmad_distribution(sample, chirp_sign=-1.0)
    reference_p0c_ev = float(sample["beta_gamma_ref"][0] * ELECTRON_MASS_MEV * 1.0e6)

    out_opalx_pos = Path(__file__).with_name("wolski_tesla_chicane_distribution.txt")
    out_opalx_neg = Path(__file__).with_name("wolski_tesla_chicane_distribution_negchirp.txt")
    out_bmad_unchirped = Path(__file__).with_name("wolski_tesla_chicane_bmad_unchirped.beam")
    out_bmad_pos = Path(__file__).with_name("wolski_tesla_chicane_bmad_poschirp.beam")
    out_bmad_neg = Path(__file__).with_name("wolski_tesla_chicane_bmad_negchirp.beam")

    write_opalx_distribution(out_opalx_pos, opalx_pos)
    write_opalx_distribution(out_opalx_neg, opalx_neg)
    write_bmad_distribution(out_bmad_unchirped, bmad_unchirped, reference_p0c_ev)
    write_bmad_distribution(out_bmad_pos, bmad_pos, reference_p0c_ev)
    write_bmad_distribution(out_bmad_neg, bmad_neg, reference_p0c_ev)

    sigma_z_final_linear = abs(R56_M) * SIGMA_DELTA_UNCORRELATED
    beta_gamma = float(sample["beta_gamma_ref"][0])
    print(f"Wrote {out_opalx_pos}")
    print(f"Wrote {out_opalx_neg}")
    print(f"Wrote {out_bmad_unchirped}")
    print(f"Wrote {out_bmad_pos}")
    print(f"Wrote {out_bmad_neg}")
    print(f"beta*gamma = {beta_gamma:.12f}")
    print(f"Linear target sigma_z,f = {sigma_z_final_linear * 1e3:.6f} mm")


if __name__ == "__main__":
    main()
