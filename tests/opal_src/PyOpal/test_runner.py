"""Module to run the unit tests"""

import sys
import unittest

def main():
    """A very simple test runner script"""
    runner = unittest.TextTestRunner()
    suite = unittest.defaultTestLoader.discover(
                start_dir = "tests/opal_src/PyOpal/",
                pattern = "test*"
    )
    result = runner.run(suite)
    n_errors = len(result.errors)
    n_failures = len(result.failures)
    print(f"Ran tests with {n_errors} errors and {n_failures} failures")
    if result.wasSuccessful():
        print("Tests passed")
        sys.exit(0)
    else:
        print("Tests failed (don't forget to make install...)")
        sys.exit(len(result.errors)+len(result.failures))

if __name__ == "__main__":
    main()
