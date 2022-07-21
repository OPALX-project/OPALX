"""Test the track module"""
import unittest
import pyopal.objects.line
import pyopal.objects.beam
import pyopal.objects.track

class TestTrack(unittest.TestCase):
    """Quick check that we can run track okay"""

    def setUp(self):
        """Set up some data"""
        self.track = pyopal.objects.track.Track()

    def test_init(self):
        """Check I didn't make any typos"""
        self.track.line = "test_line"
        self.track.beam = "test_beam"
        self.track.dt = 1.0
        self.track.dt_space_charge = 2.0
        self.track.dtau = 3.0
        self.track.t0 = 4.0
        self.track.max_steps = [5]
        self.track.steps_per_turn = 6.0
        self.track.z_start = 7.0
        self.track.z_stop = [8.0]
        self.track.time_integrator = "integrator"
        self.track.map_order = 9

    def test_execute(self):
        """Check we can register the beam"""
        beam = pyopal.objects.beam.Beam()
        beam.set_opal_name("test_track::beam")
        beam.register()

        line = pyopal.objects.line.Line()
        line.register()
        #line.set_opal_name("test_track::line") does not work

        self.track.beam = "test_track::beam"
        self.track.line = "LINE" # all lines are called LINE
        self.track.execute()


if __name__ == "__main__":
    unittest.main()
