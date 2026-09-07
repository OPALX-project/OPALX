"""Sanity checks for the independent hard-edge benchmark optics."""
import unittest
import numpy as np
from analytic_scan import cell, describe, OPTICS


class TestAnalyticScan(unittest.TestCase):
    def test_original_achromat(self):
        np.testing.assert_allclose(cell(-6.371966681365967), OPTICS['dba'](), atol=1e-14)
        self.assertFalse(describe(-6.371966681365967)['stable'])

    def test_zero_strength_is_drift(self):
        self.assertTrue(np.isfinite(cell(0)).all())
        np.testing.assert_allclose(cell(2.84, 0), cell(2.84), atol=1e-14)

    def test_selected_ring(self):
        m = np.linalg.matrix_power(cell(2.84), 6)
        self.assertTrue(describe(2.84)['stable'])
        self.assertGreater(abs(m[0, 5]), 1.)
        np.testing.assert_allclose(abs(np.linalg.eigvals(m[:4, :4])), 1., atol=1e-12)
        j = np.array([[0, 1], [-1, 0]])
        for block in (m[:2, :2], m[2:4, 2:4]):
            np.testing.assert_allclose(block.T @ j @ block, j, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
