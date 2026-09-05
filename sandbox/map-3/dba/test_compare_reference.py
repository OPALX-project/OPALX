"""Sanity checks of the separate CAS thin-edge analytic comparison."""
from dataclasses import replace
import unittest
import numpy as np

from analytic_fields import load_parameters
from compare_reference import analytic_thin_edge
from test_variational_reference import exact


class ThinEdgeTests(unittest.TestCase):
    def test_hard_edge_limit_recovers_map2_analytic_matrix(self):
        p = replace(load_parameters(), hgap_m=0)
        np.testing.assert_allclose(analytic_thin_edge(p), exact.dba(), rtol=0, atol=3e-14)

    def test_thin_edge_map_is_symplectic_to_roundoff(self):
        matrix = analytic_thin_edge(load_parameters())
        j = np.kron(np.eye(3), np.array([[0., 1.], [-1., 0.]]))
        np.testing.assert_allclose(matrix.T @ j @ matrix, j, rtol=0, atol=1e-13)
        self.assertLess(abs(np.linalg.det(matrix)-1), 1e-13)

    def test_fint_only_changes_vertical_block_for_zero_face_rotation(self):
        p = load_parameters()
        delta = analytic_thin_edge(p) - analytic_thin_edge(replace(p, fint=0))
        self.assertGreater(abs(delta[2:4, 2:4]).max(), 0.1)
        delta[2:4, 2:4] = 0
        np.testing.assert_array_equal(delta, 0)


if __name__ == "__main__":
    unittest.main()
