"""Test that distribution parses okay"""
import unittest
import pyopal.objects.distribution

class TestDistribution(unittest.TestCase):
    """Test that distribution parses okay"""
    def test_init(self):
        """Check we can initialise variables without a problem"""
        my_distribution = pyopal.objects.distribution.Distribution()
        my_distribution.type = "cheese"
        self.assertEqual(my_distribution.type, "CHEESE")
        my_distribution.fname = "disttest.dat"
        self.assertEqual(my_distribution.fname, "disttest.dat")
        my_distribution.momentum_units = "asdasd"
        self.assertEqual(my_distribution.momentum_units, "ASDASD")

    def test_register(self):
        """Check we can call register function"""
        my_distribution = pyopal.objects.distribution.Distribution()
        my_distribution.register()

if __name__ == "__main__":
    unittest.main()
