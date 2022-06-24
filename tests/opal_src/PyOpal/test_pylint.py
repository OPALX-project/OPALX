import tempfile
import os
import unittest
import subprocess

class PyLintTest(unittest.TestCase):
    """
    Because there is no actual python code in PyOpal, we use a pylint run 
    against the tests as a proxy for checking that code is reasonably pythonic

    Would be better not to - sometimes test code is (deliberately) not pythonic
    """
    def setUp(self):
        """set up the test"""
        self.log = tempfile.NamedTemporaryFile("w+")
        # need an environment variable to define OPAL test location - does this
        # exist?
        self.opal_test_path = "tests/opal_src/PyOpal/"
        self.pass_score = 8
        self.pylint_cmd = ["pylint"]


    def run_pylint(self, dirpath, fname):
        """Run pylint and fill the log tempfile"""
        full_path = os.path.join(dirpath, fname)
        subprocess.run(self.pylint_cmd+[full_path],
            stdout=self.log, stderr=subprocess.STDOUT)

    def get_score_from_line(self, line):
        if "Your code has been rated at" not in line:
            return -1
        score = line.split("rated at ")[1:]
        score = score[0].split("/10")[0]
        score = float(score)
        return score

    def check_logfile(self):
        self.log.flush()
        self.log.seek(0)
        buffer = ""
        for line in self.log.readlines():
            buffer += line
            score = self.get_score_from_line(line)
            if score < 0:
                continue
            msg="""

Failed static code analysis with output

{0}

Score must be greater than {1} to pass. Try running pylint to check manually.
            """.format(buffer, self.pass_score)
            self.assertGreater(score, self.pass_score, msg)
            buffer = ""

    def test_check_pylint(self):
        for dirpath, dirnames, filenames in os.walk(top=self.opal_test_path):
            for fname in filenames:
                if fname[-3:] != ".py":
                    continue
                self.log.truncate(0)
                self.run_pylint(dirpath, fname)
                self.check_logfile()

    def test_TODO(Self):
        self.assertTrue(False, msg="""
        * PyOpal module initialisation which wraps BOOST_PYTHON_MODULE and does
        ** initialise globals
        ** initialise exception handling
        ** handle exceptions in module initialisation (which is not handled by 
          Boost exception handling)
        * Also a test module for pyobject.h
        * Some way to share the gmsg object rather than initialise once for each module
        e.g. a "core" module that is always imported during pyopal import and 
        binding to pull the object across during Global initialisation
        * Get predefined string to check against the list of predefines
        """

    verbose = False

if __name__ == "__main__":
    unittest.main()