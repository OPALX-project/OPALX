import os
import unittest
import tempfile
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

DISTRIBUTION = """10
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
"""

class TestTrackRun(unittest.TestCase):
    def make_field_solver(self):
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
        drift = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        drift.end_position_x=0.0
        drift.end_position_y=0.0
        drift.end_normal_x=1.0
        drift.end_normal_y=0.0
        return drift

    def make_line(self):
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
        global DISTRIBUTION
        self.distribution_file = tempfile.NamedTemporaryFile("w+")
        self.distribution_file.write(DISTRIBUTION)
        self.distribution_file.flush()
        self.distribution = pyopal.objects.distribution.Distribution()
        self.distribution.set_name("SuperDist")
        self.distribution.type = "FROMFILE"
        self.distribution.fname = self.distribution_file.name
        self.distribution.register()

        return self.distribution

    def run_one(self):
        beam = pyopal.objects.beam.Beam()
        beam.mass = 0.938272
        beam.charge = 1.0
        beam.momentum = 0.1
        beam.beam_frequency = 1.0
        beam.number_of_slices = 10
        beam.number_of_particles = 10
        beam.register()


        line = self.make_line()
        field_solver = self.make_field_solver()
        distribution = self.make_distribution()

        track = pyopal.objects.track.Track()
        track.line = "LINE"
        track.beam = "BEAM"
        run = pyopal.objects.track_run.TrackRun()
        run.method = "CYCLOTRON-T"
        run.keep_alive = True
        run.beam_name = "BEAM"
        run.distribution = ["SuperDist"]
        run.field_solver = "FIELDSOLVER"
        print(pyopal.objects.parser.list_objects())
        track.execute()
        run.execute()

    def setUp(self):
        """Define a few default variables"""
        self.track_run = pyopal.objects.track_run.TrackRun()

    def test_initialisation(self):
        """Test that we can set member data okay"""
        self.run_one()

if __name__ == "__main__":
    unittest.main()