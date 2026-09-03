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

    reference_path, reference_active = design_path(
        args.reference / "isis_sbend_survey_DesignPath.dat"
    )
    candidate_path, candidate_active = design_path(
        args.candidate / "isis_sbend_ring_DesignPath.dat"
    )
    assert_numeric_equal("DesignPath", reference_path, candidate_path, args.atol)
    if reference_active != candidate_active:
        raise AssertionError("DesignPath active-element names differ")

    print(
        f"PASS: ElementPositions ({len(reference_positions)} rows) and "
        f"DesignPath ({len(reference_path)} rows), atol={args.atol:g}"
    )


if __name__ == "__main__":
    main()
