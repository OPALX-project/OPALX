"""Module to test track run (and by extension the entire PyOpal workflow)"""

import os
import sys
import tempfile
import subprocess
import unittest

import pyopal.elements.vertical_ffa_magnet
import pyopal.objects.track_run
import pyopal.objects.beam
import pyopal.objects.distribution
import pyopal.objects.line
import pyopal.elements.ring_definition
import pyopal.elements.local_cartesian_offset
import pyopal.objects.field_solver
import pyopal.objects.track
import pyopal.objects.parser

class TrackRunExecute():
    """Module to test track run (and by extension the entire PyOpal workflow)"""
    def __init__(self):
        """Set up the test"""
        self.here = os.getcwd()
        self.tmp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.tmp_dir.name)
        self.field_solver = None
        self.line = None
        self.ring = None
        self.offset = None
        self.distribution = None
        self.distribution_file = tempfile.NamedTemporaryFile("w+")
        self.track_run = pyopal.objects.track_run.TrackRun()


    def make_field_solver(self):
        """Make an empty fieldsolver"""
        self.field_solver = pyopal.objects.field_solver.FieldSolver()
        self.field_solver.type = "NONE"
        self.field_solver.mesh_size_x = 5
        self.field_solver.mesh_size_y = 5
        self.field_solver.mesh_size_t = 5
        self.field_solver.parallelize_x = False
        self.field_solver.parallelize_y = False
        self.field_solver.parallelize_t = False
        self.field_solver.boundary_x = "open"
        self.field_solver.boundary_y = "open"
        self.field_solver.boundary_t = "open"
        self.field_solver.bounding_box_increase = 2
        self.field_solver.register()

    @classmethod
    def make_drift(cls):
        """Returns a drift of length 0"""
        drift = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        drift.end_position_x=0.0
        drift.end_position_y=0.0
        drift.end_normal_x=1.0
        drift.end_normal_y=0.0
        return drift

    def make_line(self):
        """Make a beamline, just consisting of a drift section"""
        drift = self.make_drift()
        self.line = pyopal.objects.line.Line()
        self.ring = pyopal.elements.ring_definition.RingDefinition()
        self.ring.lattice_initial_r = 4.0
        self.ring.beam_initial_r = 0.0
        self.ring.minimum_r = 0.5
        self.ring.maximum_r = 10.0
        self.ring.is_closed = False
        self.offset = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        self.offset.end_position_x = 0.0
        self.offset.end_position_y = 1.0
        self.offset.normal_x = 1.0

        self.line.append(self.ring)
        self.line.append(drift)
        self.line.register()

    def make_distribution(self):
        """Make a distribution, from tempfile data"""
        self.distribution_file.write(self.distribution_str)
        self.distribution_file.flush()
        self.distribution = pyopal.objects.distribution.Distribution()
        self.distribution.set_opal_name("SuperDist")
        self.distribution.type = "FROMFILE"
        self.distribution.fname = self.distribution_file.name
        self.distribution.register()

        return self.distribution

    def run_one(self):
        """Set up and run a simulation"""
        self.make_line()
        self.make_distribution()
        self.make_field_solver()

        beam = pyopal.objects.beam.Beam()
        beam.set_opal_name("SuperBeam")
        beam.mass = 0.938272
        beam.charge = 1.0
        beam.momentum = 0.1
        beam.beam_frequency = 1.0
        beam.number_of_slices = 10
        beam.number_of_particles = 1
        beam.register()

        track = pyopal.objects.track.Track()
        track.line = "LINE"
        track.beam = "SuperBeam"
        run = pyopal.objects.track_run.TrackRun()
        run.method = "CYCLOTRON-T"
        run.keep_alive = True
        run.beam_name = "SuperBeam"
        run.distribution = ["SuperDist"]
        run.field_solver = "FIELDSOLVER"
        track.execute()
        run.execute()
        track.execute()
        run.execute()

    def __del__(self):
        """move back to the old cwd"""
        os.chdir(self.here)
        self.tmp_dir.cleanup()
        self.distribution_file.close()


    distribution_str = """1
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
"""


class TestTrackRun(unittest.TestCase):
    """Test class for track_run"""
    def test_run_one(self):
        """
        Test that we can run okay without an exception.

        If running from the command line, will spit out the OPAL log to screen
        """
        log = tempfile.TemporaryFile("w+")
        proc = subprocess.run(["python", __file__, "run_test_track_run"],
            stdout=log, stderr=subprocess.STDOUT, check=False)
        log.seek(0)
        error_message=""
        if self.verbose:
            error_message = log.read()
        log.close()
        self.assertEqual(proc.returncode, 0, msg=error_message)

    verbose = False

if __name__ == "__main__":
    if "run_test_track_run" in sys.argv:
        TrackRunExecute().run_one()
    else:
        TestTrackRun.verbose = True
        unittest.main()
