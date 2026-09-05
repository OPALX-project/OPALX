#!/usr/bin/env python3
"""Fast tests of sweep provenance and failure handling; no OPALX launches."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from convergence_grid import dba, np, run_one, summarize


class GridRecordingTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="opalx-grid-test-")
        self.root = Path(self.temporary.name)
        self.kwargs = dict(executable=Path("unused"), source="TRACK, DT=1e-12, MAXSTEPS=1;\n",
                           output=self.root, provenance={"binary": "test"}, dt=1e-12,
                           epsilon=1e-3, levels=1, ranks=1, mpi_args=())

    def tearDown(self):
        self.temporary.cleanup()

    def fake_run(self, executable, source, dt, temporary, output, ranks, mpi_args):
        matrix = dba()
        stdout = "Richardson levels: 1;\nOrbThreader......... Wall max = 0.25\n"
        (output / "map-2-dba-dt-1e-12.out").write_text(stdout)
        np.savetxt(output / "map-2-dba-dt-1e-12-matrix.txt", matrix, fmt="%.16e")
        return dict(status="OK", dt_s=dt, R16=matrix[0, 5], R26=matrix[1, 5],
                    determinant_error=0.0, canonical_J_error=0.0)

    def record(self):
        with patch("convergence_grid.run_case", side_effect=self.fake_run):
            return run_one(**self.kwargs)

    def test_records_full_matrix_and_norm(self):
        result = self.record()
        self.assertEqual(result["scaled_max_error"], 0.0)
        self.assertEqual(result["relative_frobenius_error"], 0.0)
        self.assertEqual(result["rays_per_segment"], 24)
        self.assertEqual(result["orbit_threader_wall_s"], 0.25)
        self.assertEqual(len([name for name in result if len(name) == 3 and name.startswith("R")]), 36)

    def test_matching_cache_does_not_launch(self):
        self.record()
        with patch("convergence_grid.run_case", side_effect=AssertionError("must use cache")):
            result = run_one(**self.kwargs)
        self.assertTrue(result["cached"])

    def test_changed_provenance_rejected(self):
        self.record()
        with self.assertRaisesRegex(RuntimeError, "provenance mismatch"):
            run_one(**(self.kwargs | {"provenance": {"binary": "different"}}))

    def test_corrupt_artifact_rejected(self):
        result = self.record()
        (self.root / result["path"] / "map-2-dba.in").write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "cached artifact changed"):
            run_one(**self.kwargs)

    def test_unconfirmed_level_rejected(self):
        with patch("convergence_grid.run_case", side_effect=self.fake_run):
            with self.assertRaisesRegex(RuntimeError, "Richardson setting not confirmed"):
                run_one(**(self.kwargs | {"levels": 2}))

    def rejected(self, reason):
        def run(executable, source, dt, temporary, output, ranks, mpi_args):
            (output / "map-2-dba-dt-1e-12.out").write_text(reason)
            return dict(status="INVALID", dt_s=dt)
        return run

    def test_expected_guard_is_recorded(self):
        with patch("convergence_grid.run_case", side_effect=self.rejected(
                "The time step is too long compared to the length of the\nelement 'QACH'")):
            result = run_one(**self.kwargs)
        self.assertEqual(result["status"], "INVALID")

    def test_unexpected_failure_is_not_hidden(self):
        with patch("convergence_grid.run_case", side_effect=self.rejected("unknown input keyword")):
            with self.assertRaisesRegex(RuntimeError, "unexpected OPALX failure"):
                run_one(**self.kwargs)

    def test_rejected_only_sweep_can_be_summarized(self):
        with patch("convergence_grid.run_case", side_effect=self.rejected(
                "The time step is too long compared to the length of the\nelement 'QACH'")):
            run_one(**self.kwargs)
        with patch("builtins.print"):
            frame = summarize(self.root, plots=False)
        self.assertEqual(len(frame), 1)
        self.assertEqual(frame.iloc[0].status, "INVALID")


if __name__ == "__main__":
    unittest.main()
