#!/usr/bin/env python3
"""Manufactured data model for test-chicane-distribution-1."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np


P0_BETA_GAMMA = 1957.950928190185
LB = 1.0
TH = 0.1
D_CH = 1.5
C_CH = 2.0
R_CH = LB / (2.0 * math.sin(TH / 2.0))
ARC_LENGTH = R_CH * TH
PATH_LENGTH = 4.0 * ARC_LENGTH + 2.0 * D_CH + C_CH

# Same small-angle target as test-chicane-1/check-r56.py.
R56_REF = -(TH**2) * (4.0 * LB / 3.0 + 2.0 * D_CH)


@dataclass(frozen=True)
class ManufacturedPaths:
    distribution: Path
    opal_distribution: Path
    solution_csv: Path
    solution_json: Path


def build_distribution() -> np.ndarray:
    """Return rows with columns id, x, px, y, py, z, pz, delta.

    The grid is deterministic, centered, and deliberately small enough that the
    linear manufactured longitudinal model is the dominant observable.  It has
    finite spatial extent in x, y, and z so OPALX computes non-degenerate bunch
    bounds before tracking.
    """

    deltas = np.array([-2.0e-3, -1.0e-3, 0.0, 1.0e-3, 2.0e-3])
    z_values = np.array([-1.0e-3, -5.0e-4, 0.0, 5.0e-4, 1.0e-3])
    transverse = np.array(
        [
            [-1.0e-4, -8.0e-5],
            [-5.0e-5, 4.0e-5],
            [0.0, 0.0],
            [5.0e-5, -4.0e-5],
            [1.0e-4, 8.0e-5],
        ]
    )

    rows = []
    particle_id = 0
    for delta in deltas:
        for z in z_values:
            for x, y in transverse:
                rows.append(
                    [
                        particle_id,
                        x,
                        0.0,
                        y,
                        0.0,
                        z,
                        P0_BETA_GAMMA * (1.0 + delta),
                        delta,
                    ]
                )
                particle_id += 1
    return np.asarray(rows, dtype=float)


def manufactured_solution(distribution: np.ndarray) -> np.ndarray:
    """Return expected final rows for the independent linear target."""

    result = distribution.copy()
    z_initial = distribution[:, 5]
    delta = distribution[:, 7]
    result[:, 5] = z_initial - R56_REF * delta
    return result


def write_distribution(path: Path, distribution: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(distribution)}\n")
        stream.write("x px y py z pz\n")
        np.savetxt(stream, distribution[:, 1:7], fmt="%.16e")


def write_opal_distribution(path: Path, distribution: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(distribution)}\n")
        np.savetxt(stream, distribution[:, 1:7], fmt="%.16e")


def write_solution_csv(path: Path, distribution: np.ndarray, solution: np.ndarray) -> None:
    header = (
        "id,x_initial,px_initial,y_initial,py_initial,z_initial,pz_initial,"
        "delta,z_expected,y_expected,pz_expected"
    )
    data = np.column_stack(
        (
            distribution[:, 0],
            distribution[:, 1],
            distribution[:, 2],
            distribution[:, 3],
            distribution[:, 4],
            distribution[:, 5],
            distribution[:, 6],
            distribution[:, 7],
            solution[:, 5],
            solution[:, 3],
            solution[:, 6],
        )
    )
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.16e")


def write_solution_json(path: Path, distribution: np.ndarray, solution: np.ndarray) -> None:
    delta = distribution[:, 7]
    z_expected = solution[:, 5]
    payload = {
        "p0_beta_gamma": P0_BETA_GAMMA,
        "lb_m": LB,
        "theta_rad": TH,
        "outer_drift_m": D_CH,
        "central_drift_m": C_CH,
        "rho_m": R_CH,
        "arc_length_m": ARC_LENGTH,
        "path_length_m": PATH_LENGTH,
        "r56_ref_m": R56_REF,
        "n_particles": int(len(distribution)),
        "delta_values": sorted(float(value) for value in np.unique(delta)),
        "z_expected_mean_m": float(np.mean(z_expected)),
        "z_expected_rms_m": float(np.std(z_expected)),
        "y_expected_rms_m": float(np.std(solution[:, 3])),
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_all(directory: Path) -> ManufacturedPaths:
    distribution = build_distribution()
    solution = manufactured_solution(distribution)

    paths = ManufacturedPaths(
        distribution=directory / "test-chicane-distribution-1_distribution.txt",
        opal_distribution=directory / "test-chicane-distribution-1_opal_distribution.txt",
        solution_csv=directory / "manufactured_solution.csv",
        solution_json=directory / "manufactured_solution.json",
    )
    write_distribution(paths.distribution, distribution)
    write_opal_distribution(paths.opal_distribution, distribution)
    write_solution_csv(paths.solution_csv, distribution, solution)
    write_solution_json(paths.solution_json, distribution, solution)
    return paths
