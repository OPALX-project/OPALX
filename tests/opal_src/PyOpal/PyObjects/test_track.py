import unittest
import pyopal.objects.track

class TestTrack(unittest.TestCase):
    def setUp(self):
        """Set up some data"""
        self.track = pyopal.objects.track.Track()

    def test_init(self):
        """Check I didn't make any typos"""
        pass

    def test_execute(self):
        """Check we can register the beam"""
        self.track.execute()


if __name__ == "__main__":
    unittest.main()
