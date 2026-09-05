"""Validate input reproduction, map-2 selection and independent field poses."""
from pathlib import Path
import unittest

from analytic_fields import lattice_field, layout, load_parameters, np
from run_best_cases import ROOT, numeric, select_best, validate_physics
from validate_opalx_fields import arc_coordinates, global_field, placed_elements, rotation_y


class OpalxCaseTests(unittest.TestCase):
    def test_input_matches_analytic_model(self):
        source = (ROOT / "map-3-dba.in").read_text()
        validate_physics(source, load_parameters())
        with self.assertRaises(AssertionError):
            validate_physics(source.replace("REAL HALF_GAP = 0.1", "REAL HALF_GAP = 0.01"), load_parameters())

    def test_numeric_parser_rejects_code(self):
        self.assertEqual(numeric("2*W+1", {"W": 1}), 3)
        with self.assertRaises(ValueError):
            numeric("__import__('os').getcwd()", {})

    def test_archived_best_settings_and_hashes(self):
        selected = select_best()
        expected = [("BORIS", 1e-13, 1e-5, 1), ("RK4", 3e-11, 3e-5, 0), ("DOP853", 1e-11, .003, 1)]
        self.assertEqual([(s["integration_method"], s["dt_s"], s["epsilon"], s["richardson_levels"])
                          for s in selected], expected)

    def test_global_fields_match_local_analytic_design_line(self):
        p = load_parameters()
        for row, origin, rotation in placed_elements(p):
            length = row.body_end_m - row.body_start_m
            for local_z in (-0.2, 0.1, length / 2, length + 0.2) if row.kind == "bend" else (length / 2,):
                local = np.array([0.0, p.probe_y_m, local_z])
                tangent = rotation
                if row.kind == "bend" and local_z > 0:
                    phi = min(local_z, length) / p.radius_m
                    local[0] = p.radius_m * (np.cos(phi) - 1)
                    local[2] = p.radius_m * np.sin(phi)
                    if local_z > length:
                        local += rotation_y(-phi) @ np.array([0, 0, local_z - length])
                    tangent = rotation @ rotation_y(-phi)
                position = origin + rotation @ local
                measured = global_field(position, p)[0]
                expected = tangent @ lattice_field(row.body_start_m + local_z, 0, p.probe_y_m, p)
                np.testing.assert_allclose(measured, expected, rtol=2e-12, atol=2e-14)

    def test_arc_chart_round_trip(self):
        p = load_parameters()
        phi = 0.2
        expected = np.array([0.003, 0.001, phi * p.radius_m])
        local = [[(p.radius_m + expected[0]) * np.cos(phi) - p.radius_m,
                  expected[1], (p.radius_m + expected[0]) * np.sin(phi)]]
        np.testing.assert_allclose(arc_coordinates(np.array(local), p)[0], expected, atol=1e-15)


if __name__ == "__main__":
    unittest.main()
