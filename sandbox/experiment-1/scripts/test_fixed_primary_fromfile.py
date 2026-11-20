#!/usr/bin/env python3
"""Tests for the deterministic primary FROMFILE rank-scan source."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
GENERATOR_PATH = SCRIPT_DIR / "make_opalx_field_cases.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("fixed_primary_generator", GENERATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {GENERATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FixedPrimaryFromFileTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.generator = load_generator()
        field_model = cls.generator.load_field_model()
        cls.parameters = field_model.load_parameters(field_model.DEFAULT_CONFIG)
        cls.track12 = field_model.load_track12_module()
        cls.template = cls.generator.DEFAULT_TEMPLATE.read_text(encoding="utf-8")

    def test_fixed_file_is_reproducible_recentred_and_absolute_momentum(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "first.fromfile"
            second = root / "second.fromfile"
            metadata_first = self.generator.write_deterministic_primary_file(
                first, self.parameters, self.track12, 256, 20260629
            )
            metadata_second = self.generator.write_deterministic_primary_file(
                second, self.parameters, self.track12, 256, 20260629
            )

            self.assertEqual(metadata_first["sha256"], metadata_second["sha256"])
            self.assertEqual(first.read_bytes(), second.read_bytes())
            data = np.loadtxt(first, skiprows=2)
            self.assertEqual(data.shape, (256, 6))
            np.testing.assert_allclose(np.mean(data[:, :3], axis=0), 0.0, atol=1.0e-18)

            position_sigmas = np.array(
                [
                    self.parameters.sigma_x_m,
                    self.parameters.sigma_y_m,
                    self.parameters.sigma_z_m,
                ]
            )
            self.assertLess(np.max(np.abs(data[:, :3] / position_sigmas)), 3.1)
            np.testing.assert_allclose(
                np.mean(data[:, 5]), metadata_first["beta_gamma"], rtol=0.0, atol=1.0e-11
            )

    def test_fixed_deck_uses_fromfile_and_omits_primary_pc(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            input_path = self.generator.render_case(
                self.template,
                self.parameters,
                self.track12,
                3.0,
                case_dir,
                256,
                16,
                32,
                primary_particle_filename="../../../primary_fixed.fromfile",
            )
            rendered = input_path.read_text(encoding="utf-8")
            primary_distribution = rendered.split(
                "DIST_PRIMARY_ELECTRONS: DISTRIBUTION,", maxsplit=1
            )[1].split("DIST_FIELD_PROBES:", maxsplit=1)[0]
            primary_beam = rendered.split("PRIMARY_ELECTRONS: BEAM,", maxsplit=1)[1].split(
                "FIELD_PROBES: BEAM,", maxsplit=1
            )[0]
            self.assertIn("TYPE = FROMFILE", primary_distribution)
            self.assertIn('FNAME = "../../../primary_fixed.fromfile"', primary_distribution)
            self.assertNotIn("PC =", primary_beam)


if __name__ == "__main__":
    unittest.main()
