#!/usr/bin/env python3
"""Patch OpalRegressionTests/stattest.py to tolerate a missing gnuplot.

The regression runner calls gnuplot via subprocess.Popen(['gnuplot'], ...)
to generate comparison PNGs. In CSCS CI environments gnuplot is not part of
the uenv image, so the subprocess raises FileNotFoundError and breaks the
test. This patch wraps the gnuplot call in a try/except and skips plot
generation when gnuplot is unavailable; the numerical pass/fail check is
still performed.
"""

import sys
from pathlib import Path


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path/to/stattest.py>", file=sys.stderr)
        sys.exit(1)

    p = Path(sys.argv[1])
    s = p.read_text()

    old = """        plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
        plot.communicate(bytes(plotcmd, "UTF-8"))
        os.remove(stat_plot_file)
        os.remove(reference_plot_file)
"""

    new = """        try:
            plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
            plot.communicate(bytes(plotcmd, "UTF-8"))
        except FileNotFoundError:
            print("WARNING: gnuplot not found, skipping plot generation for " + self.name)
            os.remove(stat_plot_file)
            os.remove(reference_plot_file)
            return ""
        except Exception as e:
            print("WARNING: gnuplot failed for " + self.name + ": " + str(e))
            os.remove(stat_plot_file)
            os.remove(reference_plot_file)
            return ""
        os.remove(stat_plot_file)
        os.remove(reference_plot_file)
"""

    if old not in s:
        print("ERROR: could not find gnuplot Popen block in stattest.py", file=sys.stderr)
        sys.exit(1)

    p.write_text(s.replace(old, new, 1))
    print("Patched stattest.py to tolerate missing gnuplot")


if __name__ == "__main__":
    main()
