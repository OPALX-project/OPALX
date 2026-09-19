# Isolated timed Track-12 runs

Prepare a new case without changing existing input or result files:

```bash
~/.venv-h6/bin/python sandbox/track12particles/opalx/timed/prepare_timed_track12.py \
  --output-dir /tmp/track12-new-case \
  --primary-macroparticles 400000 --nx 256 --ny 256 --nz 128
```

The generator refuses to overwrite a prepared or executed case. Its manifest
records the source, reference, and input hashes. Copy the exact input files
between CPU and GPU runs when comparing architectures. The local wrapper
`run_timed_track12.sh [new-output-directory]` prepares, runs, and compares a
case; without an argument it creates a unique directory under `timed/runs/`.

Analyze a completed case with:

```bash
~/.venv-h6/bin/python sandbox/track12particles/opalx/timed/compare_timed_track12.py \
  --run-dir /tmp/track12-new-case
```

The analyzer reconstructs positions using all three H5 reference-position
components and derives the IP from the input BeamBeam midpoint. For geometry
expressions it cannot resolve, supply the absolute IP in metres with
`--ip-s-m`. Default boundary handling follows `BCFFTX/Y`: OPEN coordinates are
never periodically unwrapped. Historical executables wrapped particles even
for OPEN fields; use `--particle-boundary legacy-periodic` explicitly to
reproduce that old fixed-aperture analysis. This mode does not describe a
dynamic field domain. The archived one-off 32 mm plotting script predates this
IP handling; use the current analyzer instead of its manual IP translation.

The comparison requires all 13,012 CAIN samples and rejects non-finite H5
phase space. Short smoke runs therefore require a separate finiteness/first-kick
check. Position units remain metres, H5 momenta are `p/(m_e c)`, and CAIN time
is `ct` in metres. See [the parent README](../README.md) for source assumptions,
historical results, and Merlin submission.

Input and coordinate tests do not launch simulations:

```bash
~/.venv-h6/bin/python -m unittest discover \
  -s sandbox/track12particles/opalx/timed -p test_timed_workflow.py -v
```

The two historical-H5 checks skip when the ignored run artifacts are absent;
synthetic OPEN/legacy coordinate tests and the isolated-input checks still run.
