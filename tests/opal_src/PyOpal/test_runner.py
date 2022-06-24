import sys
import os
import unittest

def main():
    runner = unittest.TextTestRunner()
    suite = unittest.defaultTestLoader.discover(
                start_dir = "tests/opal_src/PyOpal/",
                pattern = "test*"
    )
    result = runner.run(suite)
    print("Ran tests with {0} errors and {1} failures".format(len(result.errors), 
                                                            len(result.failures)))
    if result.wasSuccessful():
        print("Tests passed")
        sys.exit(0)
    else:
        print("Tests failed (don't forget to make install...)")
        sys.exit(len(result.errors)+len(result.failures))


if __name__ == "__main__":
    main()