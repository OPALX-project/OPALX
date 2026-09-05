"""Regression/sanity checks for the source-matched map-3 analytic field model."""

from dataclasses import replace
import unittest

import numpy as np
from scipy.integrate import quad

from analytic_fields import (bend_field, enge, envelope, lattice_field, layout,
                             load_parameters, sample_fields, summary, ROOT)


class AnalyticFieldTests(unittest.TestCase):
    def setUp(self):
        self.p = load_parameters()

    def test_layout_reserves_full_fringe_supports(self):
        p = self.p
        spans = layout(p)
        self.assertEqual(p.gap, 0.2)
        self.assertEqual(p.half_width, 1.0)
        self.assertEqual(list(spans.name), ["B1", "D1", "QACH", "D2", "B2"])
        np.testing.assert_array_equal(spans.support_end_m.iloc[:-1], spans.support_start_m.iloc[1:])
        self.assertAlmostEqual(spans.support_end_m.iloc[-1], 8.294395102393195)
        for row in spans.itertuples():
            if row.kind == "bend":
                self.assertAlmostEqual(row.body_start_m - row.support_start_m, 1.0)
                self.assertAlmostEqual(row.support_end_m - row.body_end_m, 1.0)

    def test_drifts_field_free_even_off_axis(self):
        for row in layout(self.p).itertuples():
            if row.kind == "drift":
                s = np.linspace(row.support_start_m + 1e-10, row.support_end_m - 1e-10, 101)
                np.testing.assert_array_equal(lattice_field(s, 0.001, 0.001, self.p), 0)

    def test_support_masks_and_hard_edge_limit(self):
        p = self.p
        np.testing.assert_array_equal(bend_field([-1.01, p.body_length + 1.01], 0, 0.001, p), 0)
        p = replace(p, hgap_m=0)
        z = np.array([-0.01, 0, p.body_length / 2, p.body_length, p.body_length + 0.01])
        expected = np.zeros((5, 3))
        expected[1:3, 1] = p.b0
        np.testing.assert_array_equal(bend_field(z, 0, 0.001, p), expected)

    def test_enge_face_and_saturation(self):
        f, fp, fpp = enge(np.array([-2, 0, 2]), self.p.gap)
        self.assertAlmostEqual(f[1], 1 / (1 + np.exp(0.478959)), places=15)
        np.testing.assert_array_equal(f[[0, 2]], [1, 0])
        np.testing.assert_array_equal(fp[[0, 2]], [0, 0])
        np.testing.assert_array_equal(fpp[[0, 2]], [0, 0])

    def test_enge_derivatives_and_active_face_chain_rule(self):
        # Central differences have O(h²) truncation and subtraction roundoff.
        # Tolerances below cover these derivative-check errors, not map errors.
        d = np.linspace(-0.25, 0.35, 17)
        step = 2e-6
        f, fp, fpp = enge(d, self.p.gap)
        plus, minus = enge(d + step, self.p.gap), enge(d - step, self.p.gap)
        np.testing.assert_allclose(fp, (plus[0] - minus[0]) / (2 * step), rtol=2e-8, atol=1e-9)
        np.testing.assert_allclose(fpp, (plus[1] - minus[1]) / (2 * step), rtol=2e-8, atol=2e-8)
        length = self.p.body_length
        z = np.array([-0.1, 0.1, length - 0.1, length + 0.1])
        a, ap, app, face = envelope(z, length, self.p.gap)
        plus = envelope(z + step, length, self.p.gap)
        minus = envelope(z - step, length, self.p.gap)
        np.testing.assert_array_equal(face, [0, 0, 1, 1])
        np.testing.assert_allclose(ap, (plus[0] - minus[0]) / (2 * step), rtol=2e-8, atol=1e-9)
        np.testing.assert_allclose(app, (plus[1] - minus[1]) / (2 * step), rtol=2e-8, atol=2e-8)

    def test_fint_invisible_on_axis_and_only_changes_bx(self):
        p = self.p
        s = np.linspace(0, layout(p).support_end_m.iloc[-1], 1001)
        np.testing.assert_array_equal(lattice_field(s, 0, 0, p), lattice_field(s, 0, 0, p, fint=0))
        with_fint = lattice_field(s, 0, p.probe_y_m, p)
        without = lattice_field(s, 0, p.probe_y_m, p, fint=0)
        np.testing.assert_array_equal(with_fint[:, 1:], without[:, 1:])
        self.assertGreater(np.max(abs(with_fint[:, 0] - without[:, 0])), 1e-6)

    def test_entry_exit_symmetry_and_edge_integral(self):
        p = self.p
        z = np.linspace(-0.3, 0.3, 71)
        left = bend_field(z, 0, p.probe_y_m, p)
        right = bend_field(p.body_length - z, 0, p.probe_y_m, p)
        np.testing.assert_allclose(left[:, :2], right[:, :2], rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(left[:, 2], -right[:, 2], rtol=1e-12, atol=1e-15)
        self.assertTrue(np.all(left[:, 0] < 0))
        h = 1 / p.radius_m
        expected = p.b0 / h * h * np.tan(h * p.hgap_m * p.fint)
        # At these dimensions the envelope reaches 1 before the active-face switch.
        # Each half-integral then recovers the distributed edge coefficient.
        for begin, end in ((-p.half_width, p.body_length / 2),
                           (p.body_length / 2, p.body_length + p.half_width)):
            measured = quad(lambda zz: float(bend_field(zz, 0, 1e-3, p)[0] / 1e-3),
                            begin, end, epsabs=1e-12, epsrel=1e-12)[0]
            self.assertAlmostEqual(measured, expected, delta=2e-12)

    def test_quadrupole_sign_matches_current_standalone_multipole(self):
        row = layout(self.p).query("kind == 'quadrupole'").iloc[0]
        s = (row.body_start_m + row.body_end_m) / 2
        field = lattice_field(s, 0.001, 0.002, self.p)
        self.assertLess(self.p.quad_gradient, 0)
        np.testing.assert_allclose(field, [self.p.quad_gradient * 0.002,
                                          self.p.quad_gradient * 0.001, 0])

    def test_summary_is_field_integral_not_orbit_or_map(self):
        report = summary(self.p, ROOT / "parameters.json")
        self.assertFalse(report["simulation_run"])
        self.assertFalse(report["matched_achromat"])
        self.assertLess(report["quadrature_error_estimate_m"], 1e-10)
        # The default Enge coefficients nearly balance the added/lost field;
        # the measured difference is only about 2.94 micrometres here.
        self.assertNotAlmostEqual(report["envelope_integral_m"], self.p.body_length, places=9)
        table = sample_fields(self.p)
        self.assertTrue(np.isfinite(table.to_numpy()).all())
        self.assertTrue(np.all(np.diff(table.s_m) > 0))
        np.testing.assert_array_equal(table[["Bx_axis_T", "Bs_axis_T"]], 0)


if __name__ == "__main__":
    unittest.main()
