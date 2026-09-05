"""Fast guards for matched numerical-study comparisons."""
import json
from pathlib import Path
import tempfile
import unittest

from compare_convergence_grids import load_paired, np, pd


class ComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.before = Path(self.temporary.name) / "before"
        self.after = Path(self.temporary.name) / "after"
        for root, error in [(self.before, 4.e-8), (self.after, 1.e-8)]:
            (root / "case").mkdir(parents=True)
            (root / "provenance.json").write_text(json.dumps({"identity": {
                "input_template_sha256": "same", "distribution_sha256": "same"}}))
            np.savetxt(root / "analytic-map.txt", np.eye(6))
            (root / "case/result.json").write_text(json.dumps({"artifacts_sha256": {
                "map-2-dba.in": "same-configured-input"}}))
            pd.DataFrame([dict(dt_s=1.e-13, epsilon=.003, richardson_levels=1, ranks=1,
                               status="OK", path="case", scaled_max_error=error,
                               abs_R16=error, abs_R26=error)]).to_csv(root / "results.csv", index=False)

    def test_gain_direction(self):
        paired = load_paired(self.before, self.after)
        self.assertEqual(paired.scaled_max_error_gain.iloc[0], 4.0)

    def test_rejects_different_grids(self):
        frame = pd.read_csv(self.after / "results.csv")
        frame["epsilon"] = .004
        frame.to_csv(self.after / "results.csv", index=False)
        with self.assertRaisesRegex(ValueError, "primary grids differ"):
            load_paired(self.before, self.after)

    def test_rejects_changed_acceptance(self):
        frame = pd.read_csv(self.after / "results.csv")
        frame["status"] = "INVALID"
        frame.to_csv(self.after / "results.csv", index=False)
        with self.assertRaisesRegex(ValueError, "acceptance/rejection"):
            load_paired(self.before, self.after)

    def test_rejects_changed_configured_input(self):
        (self.after / "case/result.json").write_text(json.dumps({"artifacts_sha256": {
            "map-2-dba.in": "different-input"}}))
        with self.assertRaisesRegex(ValueError, "configured inputs differ"):
            load_paired(self.before, self.after)


if __name__ == "__main__":
    unittest.main()
