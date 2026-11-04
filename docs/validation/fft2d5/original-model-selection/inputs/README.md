# 100-step, 10 pC inputs

Exact archived inputs for the plotted comparison:

- `opalx-10pC-100steps.in`: local ORIGINAL benchmark (executable SHA a2c87276106782554324753c889db203629f1b7f645d35e2f1bd30b00e32fa9e).
- `pyorbit-10pC-100steps.py`: native PyORBIT driver from the run's immutable provenance snapshot. Imports `check_bend_kinematics.py` alongside it.

Both track the same 5 m bunch, 10 mm transverse RMS, 10 pC, through a uniform vertical magnetic field (R=2 m) over 10 m in 100 steps. The charged grid is 512³.

These are unmodified run records, including the original absolute paths. They require the external particle distribution/NPZ, circular reference path and OPALX magnetic field map; those large source assets are not bundled here. `source.json` records the source paths and checksums. `provenance.json` records the checksums of these archived files. The OPALX input's opening comment is inherited from its committed-658 template; the executable above identifies the actual local ORIGINAL run.

`pyorbit-command.json` records the invocation. For an adapted local run, point `--source-json` to the source metadata, `--out` to a new directory, and `--periodic-check` to `integer-roll-check.json`; omit `--zero-sc`. PyORBIT and its native libraries must be installed. Source paths in the metadata must also resolve.
