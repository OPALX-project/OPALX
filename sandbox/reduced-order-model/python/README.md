# Python reference model

`rigid_two_gaussian_fields.py` evaluates the lab-frame electric and magnetic
fields of the two rigid counter-propagating primary Gaussian beams. It produces
snapshots when their centroids are 3, 2, 1, and 0 longitudinal rms sizes apart.

Run from the repository root:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/rigid_two_gaussian_fields.py
```

The script writes one CSV field table per separation, a numerical summary, a
JSON record of assumptions and derived parameters, and a combined x-z field
figure under `sandbox/reduced-order-model/outputs/fields/`. Every CSV contains
all three components of E and B in SI units; the plot shows the nonzero
components Ex, Ez, and By in the symmetry plane y=0.

Run the focused checks with:

```bash
~/.venv-h6/bin/python -m unittest \
  sandbox/reduced-order-model/python/test_rigid_two_gaussian_fields.py \
  sandbox/reduced-order-model/python/test_fixed_primary_fromfile.py
```

The historical manufactured-solution material remains in
`sandbox/attic/models/spherical-rest-gaussian/`; it is retained only as a
historical record. It is not an active reference model or acceptance test.

For the A100 convergence study, `prepare_a100_convergence.py` renders the
deduplicated case matrix, `analyze_opalx_convergence.py` creates CSV summaries
and PNG/PDF convergence figures, and `fetch_and_plot_a100_convergence.py` copies
only compact results from Merlin before invoking the analyzer locally. The
one-command Slurm submission entry point is
`../merlin/submit_a100_convergence.sh`.

`make_opalx_field_cases.py` also creates the deterministic primary source used
by the MPI scan. It writes positions in metres and absolute normalized momenta
in beta-gamma units, matching OPALX `FROMFILE`. The sampled positions reproduce
the active `GAUSS` behavior: independent normals truncated at three sigma and
then translated to an exactly zero finite-sample centroid. Preparation records
the generated file's SHA-256 and RNG metadata; the same file is referenced by
all rank counts.
