"""Numerical-analysis regressions; run with Python's unittest module."""
import unittest
from pathlib import Path
import json

import numpy as np

import school_reference as reference
from school_linear_optics import optics, field_and_gradients
from school_tunes import peak


class TuneAnalysisTest(unittest.TestCase):
    def test_complex_spectrum_resolves_both_frequency_orientations(self):
        n = np.arange(768)
        for frequency in [.2131662782, .7880554921]:
            signal = 2e-4 * np.exp(2j * np.pi * frequency * n)
            # Add DC and a weaker counter-rotating component to challenge the
            # real-signal alias which would be indistinguishable from 1-f.
            signal += 4e-5 * np.exp(-2j * np.pi * frequency * n) + .02
            self.assertAlmostEqual(peak(signal), frequency, delta=1e-8)

    def test_bad_spectral_input_is_rejected(self):
        for signal in [np.ones(100), np.arange(20), np.r_[np.ones(99), np.nan]]:
            with self.assertRaises(ValueError):
                peak(signal)

    def test_profile_gradients_match_independent_median_plane_field(self):
        for radius in [4.39, 4.42, 4.48]:
            for phi in np.radians([0., 6., 9., 13., 15., 19.]):
                value, dr, dphi = field_and_gradients(radius, phi)
                self.assertAlmostEqual(value, reference.field(radius, phi), delta=1e-13)
                h = 1e-6
                radial = (reference.field(radius + h, phi) - reference.field(radius - h, phi)) / (2 * h)
                angle_step = 1e-7
                angular = (reference.field(radius, phi + angle_step) - reference.field(radius, phi - angle_step)) / (2 * angle_step)
                self.assertAlmostEqual(dr, radial, delta=2e-7)
                self.assertAlmostEqual(dphi, angular, delta=2e-7)

    def test_cell_maps_reproduce_full_ring_cof_jacobian(self):
        linear = optics()
        matrix = np.array(json.loads((Path(__file__).parent / 'school-closed-orbit.json').read_text())['matrix'])
        for plane, indices in [('horizontal', [0, 1]), ('vertical', [2, 3])]:
            predicted = np.linalg.matrix_power(linear[plane]['matrix'], 16)
            # The stored COF matrix uses central FDSTEP=1e-5 m; the vertical
            # block has ~8e-7 truncation error. This bounds map comparison,
            # not nonlinear fixed-point closure or tune-fit convergence.
            self.assertLess(np.max(abs(predicted - matrix[np.ix_(indices, indices)])), 2e-6)
            self.assertAlmostEqual(linear[plane]['determinant'], 1., delta=1e-9)
            self.assertGreater(linear[plane]['beta_m_per_mc'], 0.)


if __name__ == '__main__':
    unittest.main()
