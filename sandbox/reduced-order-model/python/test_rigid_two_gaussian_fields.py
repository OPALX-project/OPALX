#!/usr/bin/env python3
"""Regression checks for the rigid two-Gaussian field snapshots."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parent / "rigid_two_gaussian_fields.py"
SPEC = importlib.util.spec_from_file_location("rigid_two_gaussian_fields", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"could not load {SCRIPT}")
MODEL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODEL
SPEC.loader.exec_module(MODEL)


class RigidTwoGaussianFieldsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.parameters = MODEL.load_parameters(MODEL.DEFAULT_CONFIG)
        cls.track12 = MODEL.load_track12_module()

    def test_centroids_have_requested_separation_and_approach_ip(self) -> None:
        sources = MODEL.make_sources(self.track12, self.parameters, 3.0)
        expected = 3.0 * self.parameters.sigma_z_m
        self.assertAlmostEqual(sources[1].center(0.0)[2] - sources[0].center(0.0)[2], expected)
        self.assertGreater(sources[0].beta_z, 0.0)
        self.assertLess(sources[1].beta_z, 0.0)

    def test_vectorized_fields_match_established_scalar_evaluator(self) -> None:
        sources = MODEL.make_sources(self.track12, self.parameters, 2.0)
        positions = np.array(
            [
                [0.73 * self.parameters.sigma_x_m, 0.0, -1.1 * self.parameters.sigma_z_m],
                [-1.2 * self.parameters.sigma_x_m, 0.4 * self.parameters.sigma_y_m, 0.2 * self.parameters.sigma_z_m],
                [2.0 * self.parameters.sigma_x_m, -0.8 * self.parameters.sigma_y_m, 1.7 * self.parameters.sigma_z_m],
            ]
        )
        e_batch, b_batch = MODEL.total_lab_fields_batch(
            positions,
            sources,
            self.track12,
            self.parameters.quadrature_chunk_size,
        )
        for index, position in enumerate(positions):
            e_scalar, b_scalar = self.track12.anisotropic_total_lab_fields(
                position, 0.0, sources
            )
            np.testing.assert_allclose(e_batch[index], e_scalar, rtol=2.0e-14, atol=1.0e-3)
            np.testing.assert_allclose(b_batch[index], b_scalar, rtol=2.0e-14, atol=1.0e-12)

    def test_fields_vanish_at_ip_by_pair_symmetry(self) -> None:
        for separation in (3.0, 2.0, 1.0, 0.0):
            sources = MODEL.make_sources(self.track12, self.parameters, separation)
            e_field, b_field = MODEL.total_lab_fields_batch(
                np.zeros((1, 3)),
                sources,
                self.track12,
                self.parameters.quadrature_chunk_size,
            )
            np.testing.assert_allclose(e_field, 0.0, atol=1.0e-7)
            np.testing.assert_allclose(b_field, 0.0, atol=1.0e-14)

    def test_magnetic_fields_cancel_everywhere_at_exact_overlap(self) -> None:
        sources = MODEL.make_sources(self.track12, self.parameters, 0.0)
        positions = np.array(
            [
                [self.parameters.sigma_x_m, 0.0, 0.0],
                [0.5 * self.parameters.sigma_x_m, 0.0, self.parameters.sigma_z_m],
                [-2.0 * self.parameters.sigma_x_m, self.parameters.sigma_y_m, -0.7 * self.parameters.sigma_z_m],
            ]
        )
        _e_field, b_field = MODEL.total_lab_fields_batch(
            positions,
            sources,
            self.track12,
            self.parameters.quadrature_chunk_size,
        )
        np.testing.assert_allclose(b_field, 0.0, atol=1.0e-14)


if __name__ == "__main__":
    unittest.main()
