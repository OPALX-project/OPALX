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

DISTRIBUTION = """4
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
"""

class TestTrackRun(unittest.TestCase):
    def make_field_solver(self):
        self.field_solver = pyopal.objects.field_solver.FieldSolver()
        self.field_solver.field_solver_type = "NONE"
        self.field_solver.mesh_size_x = 1
        self.field_solver.mesh_size_y = 1
        self.field_solver.mesh_size_t = 1
        self.field_solver.parallelize_x = False
        self.field_solver.parallelize_y = False
        self.field_solver.parallelize_t = False
        self.field_solver.boundary_x = "open"
        self.field_solver.boundary_y = "open"
        self.field_solver.boundary_t = "open"
        self.field_solver.bounding_box_increase = 2

        self.field_solver.register()

    @classmethod
    def make_magnet(cls):
        magnet1 = pyopal.elements.vertical_ffa_magnet.VerticalFFAMagnet()
        magnet1.b0 = 1.0
        magnet1.field_index = 1.31
        magnet1.max_horizontal_power = 12
        magnet1.centre_length = 0.5
        magnet1.end_length = 0.15
        magnet1.width = 0.5
        magnet1.height_neg_extent = 1.0
        magnet1.height_pos_extent = 1.0
        magnet1.bb_length = 12.0
        return magnet1

    def make_line(self):
        magnet1 = self.make_magnet()
        magnet2 = self.make_magnet()
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
        self.line.append(magnet1)
        self.line.append(self.offset)
        self.line.append(magnet2)
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
        beam.number_of_particles = 4
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
        field0 = line[0].get_field_value(0.0, 1.0, 0.0, 0.0)
        field1 = pyopal.objects.field.get_field_value(0.0, 1.0, 0.0, 0.0)
        field2 = line[1].get_field_value(0.0, 0.0, 0.0, 0.0)

    def setUp(self):
        """Define a few default variables"""
        self.track_run = pyopal.objects.track_run.TrackRun()

    def test_initialisation(self):
        """Test that we can set member data okay"""
        self.run_one()

if __name__ == "__main__":
    unittest.main()