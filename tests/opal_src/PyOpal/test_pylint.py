"""Module to drive pylint tests"""

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
        self.pass_score = 8.0
        self.disables = ["missing-module-docstring", "missing-class-docstring"]
        self.pylint_cmd = ["pylint"]


    def run_pylint(self, dirpath, fname):
        """Run pylint and fill the log tempfile"""
        dis_list = []
        for disable_item in self.disables:
            dis_list += ["--disable", disable_item]
        full_path = os.path.join(dirpath, fname)
        subprocess.run(self.pylint_cmd+dis_list+[full_path],
            stdout=self.log, stderr=subprocess.STDOUT)

    def get_score_from_line(self, line):
        """Extract the pylint score by scraping the text output"""
        if "Your code has been rated at" not in line:
            return -1
        score = line.split("rated at ")[1:]
        score = score[0].split("/10")[0]
        score = float(score)
        return score

    def check_logfile(self):
        """Read the log and look for errors"""
        self.log.flush()
        self.log.seek(0)
        buffer = ""
        for line in self.log.readlines():
            buffer += line
            score = self.get_score_from_line(line)
            if score < 0:
                continue
            msg=f"""

Failed static code analysis with output

{buffer}

Score must be greater than {self.pass_score} to pass. Try running pylint to check manually.
            """
            self.assertGreater(score, self.pass_score, msg)
            buffer = ""

    def test_check_pylint(self):
        """Recurse up the directory structure and run pyline on each .py file"""
        for dirpath, dirnames, filenames in os.walk(top=self.opal_test_path):
            for fname in filenames:
                if fname[-3:] != ".py":
                    continue
                self.log.truncate(0)
                self.run_pylint(dirpath, fname)
                self.check_logfile()

    verbose = False

if __name__ == "__main__":
    unittest.main()
