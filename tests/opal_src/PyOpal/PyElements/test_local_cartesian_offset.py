"""Test LocalCartesianOffset python implementation"""
import unittest
import pyopal.elements.local_cartesian_offset

class LocalCartesianOffsetTest(unittest.TestCase):
    """Tests for LocalCartesianOffset"""
    def test_init(self):
        """Check that we can build an offset okay"""
        offset = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        offset.end_position_x = 1.0
        self.assertEqual(offset.end_position_x, 1.0)
        offset.end_position_y = 2.0
        self.assertEqual(offset.end_position_y, 2.0)
        offset.end_normal_x = 3.0
        self.assertEqual(offset.end_normal_x, 3.0)
        offset.end_normal_y = 4.0
        self.assertEqual(offset.end_normal_y, 4.0)

if __name__ == "__main__":
    unittest.main()
