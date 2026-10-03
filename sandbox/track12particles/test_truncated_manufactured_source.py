#!/usr/bin/env python3
"""Focused checks for the component-wise truncated manufactured source."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().with_name("track12particles.py")


def load_model():
    spec = importlib.util.spec_from_file_location("track12_truncated_test_model", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TruncatedManufacturedSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = load_model()
        beta = cls.model.beta_from_kinetic_energy(
            cls.model.DEFAULT_SOURCE_KINETIC_MEV
        )
        gamma = 1.0 / np.sqrt(1.0 - beta * beta)
        cls.sigma_rest = np.array(
            [
                cls.model.SIGMA_XY_M,
                cls.model.SIGMA_XY_M,
                gamma * cls.model.SIGMA_Z_M,
            ]
        )

    def test_large_cutoff_recovers_untruncated_field(self) -> None:
        position = np.array(
            [0.7 * self.model.SIGMA_XY_M, -0.4 * self.model.SIGMA_XY_M, 0.2e-3]
        )
        expected = self.model.anisotropic_gaussian_rest_field(
            position, self.model.DEFAULT_SOURCE_CHARGE_C, self.sigma_rest
        )
        actual = self.model.anisotropic_truncated_gaussian_rest_field(
            position,
            self.model.DEFAULT_SOURCE_CHARGE_C,
            self.sigma_rest,
            10.0,
        )
        np.testing.assert_allclose(actual, expected, rtol=2.0e-10, atol=1.0e-6)

    def test_field_is_odd_about_source_center(self) -> None:
        position = np.array(
            [self.model.SIGMA_XY_M, 0.3 * self.model.SIGMA_XY_M, 0.1e-3]
        )
        positive = self.model.anisotropic_truncated_gaussian_rest_field(
            position,
            self.model.DEFAULT_SOURCE_CHARGE_C,
            self.sigma_rest,
            3.0,
        )
        negative = self.model.anisotropic_truncated_gaussian_rest_field(
            -position,
            self.model.DEFAULT_SOURCE_CHARGE_C,
            self.sigma_rest,
            3.0,
        )
        np.testing.assert_allclose(negative, -positive, rtol=2.0e-13, atol=1.0e-9)

    def test_three_sigma_cutoff_changes_pair4_reference_by_expected_amount(self) -> None:
        position = np.array([self.model.SIGMA_XY_M, 0.0, 0.0])
        untruncated = self.model.anisotropic_gaussian_rest_field(
            position, self.model.DEFAULT_SOURCE_CHARGE_C, self.sigma_rest
        )
        truncated = self.model.anisotropic_truncated_gaussian_rest_field(
            position,
            self.model.DEFAULT_SOURCE_CHARGE_C,
            self.sigma_rest,
            3.0,
        )
        self.assertAlmostEqual(
            np.linalg.norm(truncated) / np.linalg.norm(untruncated),
            1.0082068311,
            places=8,
        )


if __name__ == "__main__":
    unittest.main()
