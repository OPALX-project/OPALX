#!/usr/bin/env python3
"""Input/coordinate regressions; these tests never launch OPALX."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("compare_timed_track12", HERE / "compare_timed_track12.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def deck_text(length: str = "32.0e-3", edge: str = "0.0") -> str:
    return (
        f"REAL bb_length = {length};\n"
        f"REAL bb_edge = {edge};\n"
        "IP1: BEAMBEAM, L = bb_length, ELEMEDGE = bb_edge,\n"
        'APERTURE = "RECTANGLE(2.4e-3, 2.4e-4)";\n'
        "FS1: FIELDSOLVER, BCFFTX = OPEN, BCFFTY = OPEN, BCFFTZ = OPEN;\n"
    )


class CoordinateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="track12-analysis-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.deck = self.directory / "track12_timed.in"
        self.deck.write_text(deck_text(), encoding="utf-8")

    def write_samples(self, xs=(0.4e-6, -0.4e-6)) -> Path:
        path = self.directory / "track12_timed_c1.h5"
        with h5py.File(path, "w") as h5:
            for step, x in enumerate(xs, start=1):
                group = h5.create_group(f"Step#{step}")
                group.attrs["TIME"] = [step * comparison.CAIN_CT_STEP_M / comparison.C_LIGHT]
                group.attrs["GlobalTrackStep"] = [step]
                group.attrs["RefPartR"] = [0.0, 0.0, 0.016]
                for name, value in {
                    "id": 17, "x": x, "y": 0.0, "z": step * 1.0e-6,
                    "px": 0.1, "py": 0.0, "pz": 1.0,
                }.items():
                    group.create_dataset(name, data=[value])
        return path

    def test_deck_midpoint_and_open_boundary(self):
        self.assertEqual(comparison.interaction_point(self.deck), 0.016)
        self.assertIsNone(comparison.particle_periods(self.deck, [256, 256, 128]))
        self.deck.write_text(deck_text("2 * 4e-3", "1e-3"), encoding="utf-8")
        self.assertEqual(comparison.interaction_point(self.deck), 0.005)

    def test_unsupported_geometry_fails_instead_of_assuming_ip(self):
        self.deck.write_text(deck_text("unknown_length"), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "use --ip-s-m"):
            comparison.interaction_point(self.deck)

    def test_open_coordinates_are_not_unwrapped(self):
        samples = comparison.load_opalx(
            self.write_samples(), "electron", 0.0, {1: 0.0}, None
        )
        np.testing.assert_array_equal(samples["x_opalx_m"], [0.4e-6, -0.4e-6])
        np.testing.assert_allclose(samples["s_opalx_m"], [1e-6, 2e-6], rtol=0, atol=2e-18)
        self.assertIsNone(
            comparison.make_identity_summary(samples, None)["electron"]["transverse_wraps_by_pair"]
        )

    def test_explicit_legacy_periodic_unwrap(self):
        samples = comparison.load_opalx(
            self.write_samples(), "electron", 0.0, {1: 0.0}, (1e-6, 1e-6)
        )
        np.testing.assert_allclose(samples["x_opalx_m"], [0.4e-6, 0.6e-6], rtol=0, atol=1e-21)
        self.assertEqual(
            comparison.make_identity_summary(samples, (1e-6, 1e-6))["electron"]
            ["transverse_wraps_by_pair"]["1"]["x"], 1
        )
        np.testing.assert_array_equal(
            comparison.particle_periods(self.deck, [256, 256, 128], "legacy-periodic"),
            comparison.transverse_periods(self.deck, [256, 256, 128]),
        )

    def test_nonfinite_phase_space_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-finite phase space"):
            comparison.load_opalx(
                self.write_samples((np.nan,)), "electron", 0.0, {1: 0.0}, None
            )


class PreparationTests(unittest.TestCase):
    def test_isolated_deterministic_inputs_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="track12-prepare-") as temporary:
            run_dirs = [Path(temporary) / name for name in ("first", "second")]
            command = [
                sys.executable, str(HERE / "prepare_timed_track12.py"),
                "--primary-macroparticles", "64", "--nx", "16", "--ny", "16", "--nz", "32",
            ]
            for run_dir in run_dirs:
                subprocess.run(command + ["--output-dir", str(run_dir)], check=True, capture_output=True)
            first, second = run_dirs
            for name in ("primary_fixed.fromfile", "track12_electrons.emittedfromfile",
                         "track12_positrons.emittedfromfile"):
                self.assertEqual((first / "input" / name).read_bytes(), (second / "input" / name).read_bytes())
            manifest = json.loads((first / "results" / "preparation_manifest.json").read_text())
            self.assertEqual(manifest["reference_rows"], 13012)
            self.assertEqual(manifest["birth_global_steps_by_species_then_pair"], [1, 168, 334, 501, 668, 834] * 2)
            before = {path: path.read_bytes() for path in first.rglob("*") if path.is_file()}
            rerun = subprocess.run(command + ["--output-dir", str(first)], capture_output=True, text=True)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn("Refusing to overwrite", rerun.stderr)
            self.assertEqual(before, {path: path.read_bytes() for path in first.rglob("*") if path.is_file()})


class HistoricalH5Tests(unittest.TestCase):
    """Optional read-only checks using locally retained scientific run artifacts."""

    def test_legacy_samples_match_retained_comparison(self):
        run = HERE / "a100_4rank_400k_4096x256x128_twofield_816d11ff8"
        h5_path = run / "track12_timed_c1.h5"
        csv_path = run / "results" / "track12_pointwise_comparison.csv"
        if not h5_path.is_file() or not csv_path.is_file():
            self.skipTest("historical H5/CSV artifacts are not installed")
        manifest = json.loads((run / "results" / "preparation_manifest.json").read_text())
        periods = comparison.particle_periods(run / "track12_timed.in", manifest["mesh"], "legacy-periodic")
        # A prefix exercises real H5 metadata and identity tracking without
        # rerunning the expensive manufactured trajectory evaluator.
        with tempfile.TemporaryDirectory(prefix="track12-legacy-h5-") as temporary:
            prefix = Path(temporary) / "prefix.h5"
            with h5py.File(h5_path, "r") as source, h5py.File(prefix, "w") as destination:
                names = sorted(source, key=lambda name: int(name.split("#")[-1]))[:10]
                for name in names:
                    source.copy(name, destination)
            samples = comparison.load_opalx(
                prefix, "electron", manifest["pair_t0_s"], {1: -0.9e-3}, periods,
                comparison.interaction_point(run / "track12_timed.in"),
            )
        expected = pd.read_csv(csv_path)
        expected = expected.loc[(expected["species"] == "electron") & (expected["pair"] == 1)]
        merged = samples.merge(expected, on="reference_step", suffixes=("_actual", "_stored"), validate="one_to_one")
        self.assertEqual(len(merged), len(samples))
        for coordinate in ("x", "y", "s"):
            np.testing.assert_allclose(
                merged[f"{coordinate}_opalx_m_actual"], merged[f"{coordinate}_opalx_m_stored"],
                rtol=0, atol=1e-15,
            )

    def test_current_open_32mm_h5_preserves_absolute_coordinates(self):
        run = HERE / "a100_1rank_400k_256x256x128_open_efedfcc31_partial_step5231"
        h5_path = run / "track12_timed_c1.h5"
        if not h5_path.is_file():
            self.skipTest("OPEN-BC H5 artifacts are not installed")
        manifest = json.loads((run / "results" / "preparation_manifest.json").read_text())
        deck = run / "track12_timed.in"
        ip_s_m = comparison.interaction_point(deck)
        self.assertEqual(ip_s_m, 0.016)
        self.assertIsNone(comparison.particle_periods(deck, manifest["mesh"]))
        expected = []
        with tempfile.TemporaryDirectory(prefix="track12-open-h5-") as temporary:
            prefix = Path(temporary) / "prefix.h5"
            with h5py.File(h5_path, "r") as source, h5py.File(prefix, "w") as destination:
                names = sorted(source, key=lambda name: int(name.split("#")[-1]))[:10]
                for name in names:
                    source.copy(name, destination)
                    group = source[name]
                    ref = comparison.reference_position(group)
                    expected.extend([
                        [ref[0] + group["x"][i], ref[1] + group["y"][i], ref[2] + group["z"][i] - ip_s_m]
                        for i in range(len(group["id"]))
                    ])
            samples = comparison.load_opalx(
                prefix, "electron", manifest["pair_t0_s"], {1: -0.9e-3}, None, ip_s_m
            )
        np.testing.assert_array_equal(samples[["x_opalx_m", "y_opalx_m", "s_opalx_m"]], expected)


if __name__ == "__main__":
    unittest.main()
