import unittest
import tempfile
import subprocess

class ParserTest(unittest.TestCase):
    def make_temp(self, a_string):
        """Dump string to a temporary file for use in testing"""
        my_temp = tempfile.NamedTemporaryFile(mode='w+')
        my_temp.write(a_string)
        my_temp.flush()
        return my_temp

    def encapsulate_parser_in_subprocess(self):
        """OPAL can kill python execution so we hide the test in a subprocess"""
        temp_file = self.make_temp(self.good_lattice)
        temp_stdout = tempfile.TemporaryFile()
        proc = subprocess.run(["python3",
            "-c", self.command+"'"+temp_file.name+"')"],
            stdout=temp_stdout, stderr=subprocess.STDOUT, check=False)
        if proc.returncode != 0:
            temp_stdout.seek(0)
            for line in temp_stdout:
                print(line[:-1])
        temp_stdout.close()
        temp_file.close()
        return proc.returncode

    def test_parser_initialise(self):
        """Test that we can initialise some dummy lattice"""
        is_error = self.encapsulate_parser_in_subprocess()
        self.assertFalse(is_error)

    command = """
import pyopal.objects.parser
pyopal.objects.parser.initialise_from_opal_file(
"""
    good_lattice = """
///////////////////////////////////////
Title,string="Dummy lattice for testing";
///

OPTION, VERSION=20210100;

null: LOCAL_CARTESIAN_OFFSET,
                end_position_x=0., end_position_y=0.,
                end_normal_x=1.0, end_normal_y=0.;


ringdef: RINGDEFINITION, HARMONIC_NUMBER=1, LAT_RINIT=1, LAT_PHIINIT=0,
         LAT_THETAINIT=0.0, BEAM_PHIINIT=0, BEAM_PRINIT=0,
         BEAM_RINIT=1, SYMMETRY=1, RFFREQ=1, IS_CLOSED=false,
         MIN_R=0.1, MAX_R=2;

l1: Line = (ringdef, null);

Dist1: DISTRIBUTION, TYPE=gauss;

//Dist1: DISTRIBUTION, TYPE=fromfile, FNAME="disttest.dat", INPUTMOUNITS=NONE, EMITTED=TRUE, EMISSIONSTEPS=1, NBIN=1;

Fs1:FIELDSOLVER, FSTYPE=None, MX=5, MY=5, MT=5,
                 PARFFTX=true, PARFFTY=false, PARFFTT=false,
                 BCFFTX=open, BCFFTY=open, BCFFTT=open,BBOXINCR=2;

beam1: BEAM, PARTICLE=PROTON, pc=1, NPART=1, BCURRENT=1.6e-19, CHARGE=1.0, BFREQ=1;

TRACK, LINE=l1, BEAM=beam1, MAXSTEPS=1, STEPSPERTURN=1;
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=Fs1, DISTRIBUTION=Dist1;
ENDTRACK;
STOP;
"""

if __name__ == "__main__":
    unittest.main()
