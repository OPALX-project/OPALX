#!/usr/bin/env python3
"""Compare the RING OrbitThreader outputs with the explicit golden run."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def element_positions(path: Path) -> pd.DataFrame:
    pattern = re.compile(
        r'^"([^"]+)"\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)'
    )
    rows = []
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            rows.append((match.group(1), *(float(value) for value in match.groups()[1:])))
    frame = pd.DataFrame(rows, columns=["label", "x", "z", "y"])
    frame["occurrence"] = frame.groupby("label").cumcount()
    return frame.sort_values(["label", "occurrence"]).reset_index(drop=True)


def design_path(path: Path) -> tuple[pd.DataFrame, list[tuple[str, ...]]]:
    numeric_rows = []
    active = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=15)
        numeric_rows.append([float(value) for value in fields[:15]])
        names = fields[15].strip() if len(fields) == 16 else ""
        active.append(tuple(sorted(name.strip() for name in names.split(",") if name.strip())))
    return pd.DataFrame(numeric_rows), active


def statistics_row_near_s(path: Path, target_s: float) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 23:
            continue
        try:
            rows.append(np.asarray([float(value) for value in fields[:23]]))
        except ValueError:
            continue
    if not rows:
        raise AssertionError(f"no statistics rows found in {path}")
    return min(rows, key=lambda row: abs(row[1] - target_s))


def assert_numeric_equal(
    label: str, reference: pd.DataFrame, candidate: pd.DataFrame, atol: float
) -> None:
    if reference.shape != candidate.shape:
        raise AssertionError(
            f"{label}: shape differs: reference={reference.shape}, candidate={candidate.shape}"
        )
    if not np.allclose(reference.to_numpy(), candidate.to_numpy(), rtol=0.0, atol=atol):
        difference = np.abs(reference.to_numpy() - candidate.to_numpy())
        raise AssertionError(f"{label}: maximum absolute difference is {difference.max():.6g}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path, help="candidate data directory")
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path(__file__).resolve().parent / "orig" / "data",
    )
    # ElementPositions.txt prints only about ten significant decimal digits;
    # 2e-8 m covers its final-place rounding while remaining far below the
    # geometry scales in this lattice.
    parser.add_argument("--atol", type=float, default=2.0e-8)
    parser.add_argument("--circumference", type=float, default=163.3630220240698)
    parser.add_argument("--closure-atol", type=float, default=2.0e-2)
    parser.add_argument("--statistics", type=Path)
    args = parser.parse_args()

    reference_positions = element_positions(
        args.reference / "isis_sbend_survey_ElementPositions.txt"
    )
    candidate_positions = element_positions(
        args.candidate / "isis_sbend_ring_ElementPositions.txt"
    )
    if reference_positions[["label", "occurrence"]].to_records(index=False).tolist() != (
        candidate_positions[["label", "occurrence"]].to_records(index=False).tolist()
    ):
        raise AssertionError("ElementPositions labels or multiplicities differ")
    assert_numeric_equal(
        "ElementPositions",
        reference_positions[["x", "z", "y"]],
        candidate_positions[["x", "z", "y"]],
        args.atol,
    )

    candidate_path, _ = design_path(
        args.candidate / "isis_sbend_ring_DesignPath.dat"
    )
    if candidate_path.empty:
        raise AssertionError("DesignPath is empty")

    final = candidate_path.iloc[-1].to_numpy()
    final_s = final[0]
    closure_error = np.linalg.norm(final[1:4])
    if abs(final_s - args.circumference) > args.closure_atol:
        raise AssertionError(
            f"DesignPath ends at s={final_s:.9g} m, expected {args.circumference:.9g} m"
        )
    if closure_error > args.closure_atol:
        raise AssertionError(
            f"DesignPath one-turn position closure error is {closure_error:.9g} m"
        )
    if np.linalg.norm(final[4:6]) > 1.0e-3 or final[6] <= 0.0:
        raise AssertionError("DesignPath momentum does not return to the entrance direction")

    statistics_path = args.statistics or args.candidate.parent / "isis_sbend_ring.stat"
    statistics = statistics_row_near_s(statistics_path, args.circumference)
    particle_s = statistics[1]
    particle_closure_error = np.linalg.norm(statistics[17:20])
    if abs(particle_s - args.circumference) > args.closure_atol:
        raise AssertionError(
            f"particle ends at s={particle_s:.9g} m, expected {args.circumference:.9g} m"
        )
    if particle_closure_error > args.closure_atol:
        raise AssertionError(
            f"particle one-turn position closure error is {particle_closure_error:.9g} m"
        )
    if np.linalg.norm(statistics[20:22]) > 1.0e-3 or statistics[22] <= 0.0:
        raise AssertionError("particle momentum does not return to the entrance direction")

    print(
        f"PASS: ElementPositions ({len(reference_positions)} rows) and "
        f"one-turn DesignPath ({len(candidate_path)} rows, "
        f"closure={closure_error:.6g} m), particle closure={particle_closure_error:.6g} m, "
        f"placement atol={args.atol:g}"
    )


if __name__ == "__main__":
    main()
