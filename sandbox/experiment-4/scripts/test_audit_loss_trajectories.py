"""Synthetic controls for the immutable-ID loss audit (no OPALX run needed)."""
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from audit_loss_trajectories import IP_TIME, audit


class LossAuditTest(unittest.TestCase):
    def case(self, mode):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        path = Path(temporary.name) / "witness.h5"
        with h5py.File(path, "w") as h5:
            for step, ids in enumerate(([10, 11, 12], [12, 10], [10, 12])):
                group = h5.create_group(f"Step#{step}")
                group.attrs.update(TIME=[IP_TIME + step * 1e-12], GlobalTrackStep=[step],
                                   RefPartR=[0., 0., .17], TaitBryantAngles=[0., 0., 0.])
                group["id"] = [13 if mode == "renumber" and step > 0 and i == 12 else i for i in ids]
                for field in ("x", "y", "z", "px", "py", "pz"):
                    values = np.zeros(len(ids))
                    if field == "x" and mode == "jump" and step == 1:
                        values[0] = .01
                    if field == "z" and mode in ("outside", "outside_kick"):
                        values[:] = .02  # Lab z=.19 m, outside the fixed [.16,.18] m window.
                    if field == "px" and mode == "outside_kick" and step > 0:
                        values[0] = .5
                    group[field] = values
        return audit(path)[0]

    def test_deletion_and_slot_permutation_are_valid(self):
        result = self.case("valid")
        self.assertEqual(result["lost_ids"], 1)
        for key in ("duplicate_ids", "reappearing_ids", "new_ids_after_initial_injection",
                    "causal_displacement_violations"):
            self.assertEqual(result[key], 0)

    def test_same_id_with_unphysical_jump_is_detected(self):
        self.assertGreater(self.case("jump")["causal_displacement_violations"], 0)

    def test_renumbered_survivor_is_detected(self):
        self.assertGreater(self.case("renumber")["new_ids_after_initial_injection"], 0)

    def test_outside_window_ballistic_intervals_are_checked(self):
        result = self.case("outside")
        self.assertGreater(result["outside_same_side_intervals"], 0)
        self.assertEqual(result["max_delta_p_outside_window_mc"], 0.)
        self.assertEqual(result["max_drift_residual_outside_window_m"], 0.)

    def test_nonzero_outside_window_kick_is_detected(self):
        self.assertGreater(self.case("outside_kick")["max_delta_p_outside_window_mc"], 0.)


if __name__ == "__main__":
    unittest.main()
