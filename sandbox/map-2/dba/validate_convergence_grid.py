#!/usr/bin/env python3
"""Audit archived DBA grid results independently of the sweep's cached summaries."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from convergence_dt import dba, np, parse_map, pd
from convergence_grid import DEFAULT_DT, DEFAULT_EPSILON


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check-legacy", action="store_true",
                        help="Require exact reproduction of the earlier pre-bookkeeping matrices")
    args = parser.parse_args()
    root = args.output
    manifest = json.loads((root / "provenance.json").read_text())
    if (root / "source-snapshot").exists():
        for name, sha in manifest["identity"]["sources_sha256"].items():
            assert hashlib.sha256((root / "source-snapshot" / name).read_bytes()).hexdigest() == sha, name
    frame = pd.read_csv(root / "results.csv", float_precision="round_trip")
    primary = frame[frame.ranks == 1]
    expected_cases = {(dt, level, epsilon) for dt in DEFAULT_DT for level in range(5)
                      for epsilon in DEFAULT_EPSILON}
    actual_cases = set(zip(primary.dt_s, primary.richardson_levels, primary.epsilon))
    assert expected_cases == actual_cases, "incomplete or unexpected primary grid"
    assert len(primary) == len(expected_cases) == 385, "duplicate cases"
    rejected = frame[frame.status != "OK"]
    assert (rejected.dt_s == 1e-10).all(), "unexpected rejection"
    assert len(primary[primary.status != "OK"]) == 55
    for _, row in frame.iterrows():
        case = root / row.path
        record = json.loads((case / "result.json").read_text())
        for name, sha in record["artifacts_sha256"].items():
            assert hashlib.sha256((case / name).read_bytes()).hexdigest() == sha, case / name
        stdout = (case / f"map-2-dba-dt-{row.dt_s:.0e}.out").read_text()
        if row.status != "OK":
            assert "time step is too long" in stdout and "element 'QACH'" in stdout
            continue
        matrix = parse_map(stdout)
        saved = np.loadtxt(case / f"map-2-dba-dt-{row.dt_s:.0e}-matrix.txt")
        assert np.array_equal(matrix, saved), case
        assert np.isfinite(matrix).all()
        assert f"Richardson levels: {row.richardson_levels};" in stdout
        error = matrix - dba()
        np.testing.assert_allclose(row.scaled_max_error, np.max(np.abs(error)), rtol=1e-11, atol=0)
        np.testing.assert_allclose(row.relative_frobenius_error,
                                   np.linalg.norm(error) / np.linalg.norm(dba()), rtol=1e-11, atol=0)
        for i in range(6):
            for j in range(6):
                np.testing.assert_allclose(row[f"R{i + 1}{j + 1}"], matrix[i, j], rtol=1e-11, atol=0)
        if row.ranks > 1:
            peer = primary[(primary.dt_s == row.dt_s) & (primary.epsilon == row.epsilon) &
                           (primary.richardson_levels == row.richardson_levels)].iloc[0]
            reference = np.loadtxt(root / peer.path / f"map-2-dba-dt-{row.dt_s:.0e}-matrix.txt")
            assert np.array_equal(matrix, reference), f"MPI discrepancy: {case}"
            assert row.determinant_error == peer.determinant_error
            assert row.canonical_J_error == peer.canonical_J_error
    legacy_checks = 0
    for level in (0, 1, 2):
        legacy = root.parent / "richardson" / f"levels-{level}" / "eps-1e-03" / "ranks-1" / \
            "map-2-dba-dt-1e-13-matrix.txt"
        if args.check_legacy and legacy.exists():
            measured = root / f"L{level}" / "dt-1e-13" / "eps-0.001" / "ranks-1" / \
                "map-2-dba-dt-1e-13-matrix.txt"
            assert np.array_equal(np.loadtxt(measured), np.loadtxt(legacy))
            legacy_checks += 1
    legacy_summary = (f"{legacy_checks} earlier-study matrices reproduced exactly"
                      if args.check_legacy else "legacy cross-version equality not requested")
    print(f"PASS: {len(primary)} primary cases, {len(rejected)} expected rejections; "
          f"{len(frame[frame.ranks > 1])} MPI repeats; artifact hashes, 36 entries and norms verified; "
          f"{legacy_summary}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
