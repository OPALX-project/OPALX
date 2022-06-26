"""Check that beam can be set up okay"""

import unittest
import pyopal.objects.beam

class TestBeam(unittest.TestCase):
    """Very light beam test class"""
    def test_init(self):
        my_beam = pyopal.objects.beam.Beam()
        my_beam.particle = "proton"
        self.assertEqual(my_beam.particle, "PROTON")
        my_beam.mass = 1
        self.assertAlmostEqual(my_beam.mass, 1)
        my_beam.charge = 2
        self.assertAlmostEqual(my_beam.charge, 2)
        my_beam.energy = 3
        self.assertAlmostEqual(my_beam.energy, 3)
        my_beam.momentum = 4
        self.assertAlmostEqual(my_beam.momentum, 4)
        my_beam.gamma = 5
        self.assertAlmostEqual(my_beam.gamma, 5)
        my_beam.beam_current = 6
        self.assertAlmostEqual(my_beam.beam_current, 6)
        my_beam.beam_frequency = 7
        self.assertAlmostEqual(my_beam.beam_frequency, 7)
        my_beam.number_of_particles = 8
        self.assertAlmostEqual(my_beam.number_of_particles, 8)

    def test_bad_particle_type(self):
        """See what happens if we set a bad particle type"""
        my_beam = pyopal.objects.beam.Beam()
        my_beam.particle = "bad"
        # would not be terrible if this threw exception

    def test_register(self):
        """Check we can register the beam okay"""
        my_beam = pyopal.objects.beam.Beam()
        my_beam.register()


if __name__ == "__main__":
    unittest.main()
