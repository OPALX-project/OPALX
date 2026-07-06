#!/usr/bin/env python3
"""Generate one-particle transverse witness FROMFILE distributions.

The OPALX FROMFILE columns use positions in metres and normalized momenta
``p/(m c)``.  This generator keeps the witness physics explicit and avoids
opaque hand-edited numbers in the deck.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


EMASS_GEV = 0.00051099895


def normalized_momentum_from_kinetic_energy(kinetic_energy_gev: float) -> float:
    gamma = 1.0 + kinetic_energy_gev / EMASS_GEV
    return math.sqrt(gamma * gamma - 1.0)


def write_fromfile(path: Path, x: float, y: float, z: float, px: float, py: float, pz: float) -> None:
    path.write_text(
        "\n".join(
            [
                "1",
                "x y z px py pz",
                (
                    f"{x:.16e} {y:.16e} {z:.16e} "
                    f"{px:.16e} {py:.16e} {pz:.16e}"
                ),
                "",
            ]
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("sandbox/Drift-Experiment"))
    parser.add_argument("--kinetic-energy-kev", type=float, default=313.0)
    parser.add_argument("--x", type=float, default=1.944325075701e-6)
    parser.add_argument("--y", type=float, default=0.0)
    parser.add_argument("--z", type=float, default=0.0)
    parser.add_argument(
        "--direction",
        choices=("y", "x"),
        default="y",
        help="Transverse axis for the electron witness. The positron uses the opposite sign.",
    )
    args = parser.parse_args()

    kinetic_energy_gev = args.kinetic_energy_kev * 1.0e-6
    p_norm = normalized_momentum_from_kinetic_energy(kinetic_energy_gev)
    electron = {"px": 0.0, "py": 0.0, "pz": 0.0}
    positron = {"px": 0.0, "py": 0.0, "pz": 0.0}
    if args.direction == "y":
        electron["py"] = p_norm
        positron["py"] = -p_norm
    else:
        electron["px"] = p_norm
        positron["px"] = -p_norm

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_fromfile(
        args.output_dir / "track12_pair4_electron.fromfile",
        args.x,
        args.y,
        args.z,
        electron["px"],
        electron["py"],
        electron["pz"],
    )
    write_fromfile(
        args.output_dir / "track12_pair4_positron.fromfile",
        args.x,
        args.y,
        args.z,
        positron["px"],
        positron["py"],
        positron["pz"],
    )

    gamma = 1.0 + kinetic_energy_gev / EMASS_GEV
    beta = math.sqrt(1.0 - 1.0 / (gamma * gamma))
    print(f"kinetic_energy_keV={args.kinetic_energy_kev:.9g}")
    print(f"gamma={gamma:.16e}")
    print(f"beta={beta:.16e}")
    print(f"p_norm={p_norm:.16e}")
    print(f"direction={args.direction}")


if __name__ == "__main__":
    main()
