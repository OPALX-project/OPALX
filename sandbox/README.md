# BeamBeam sandbox

The authoritative model, implementation, and validation summary is
[`BEAMBEAM_PHYSICS_AND_VALIDATION.md`](BEAMBEAM_PHYSICS_AND_VALIDATION.md).
This README is only the operational index; detailed physics conclusions should
not be duplicated here.

## Active workspaces

| path | purpose |
|---|---|
| [`reduced-order-model/`](reduced-order-model/README.md) | independent anisotropic rigid two-Gaussian manufactured fields, OPALX field comparison, and A100 convergence workflow |
| [`cain-opalx-reduced-order-model/`](cain-opalx-reduced-order-model/README.md) | conversion and validation of the 1,297 + 1,297 exact-timed CAIN input |
| [`track12particles/opalx/`](track12particles/opalx/README.md) | exact-timed 6 e- + 6 e+ CAIN/OPALX/manufactured trajectory benchmark |
| [`track-e-p/`](track-e-p/INPUT_INVENTORY.md) | 32 cm, 15 cm-radius production inputs for the complete CAIN population |
| [`OPALX-IMPACT/`](OPALX-IMPACT/) | independent drift-with-space-charge field-solver sanity comparison |
| [`regression/`](regression/README.md) | compact broad-regression baseline and generated diagnostic figures |
| [`python/`](python/) | maintained analysis and visualization utilities |
| [`attic/`](attic/README.md) | superseded notes, models, scans, and one-off scripts retained for provenance |

The chronological design and run ledger is
[`../BEAMBEAM_REDESIGN_STATE.md`](../BEAMBEAM_REDESIGN_STATE.md). It records
historical intermediate results; the consolidated note above determines which
results are current.

## Authoritative cases

The short validation case and the production case are intentionally different:

| case | element and field domain | witnesses |
|---|---|---:|
| Track-12 validation | 8 mm; `RECTANGLE(2.4e-3,2.4e-4)` | 6 e- + 6 e+ |
| Production | 32 cm; `CIRCLE(0.30)` = 15 cm radius | 1,297 e- + 1,297 e+ |

Both use the element midpoint as the sole IP. The manufactured source is the
anisotropic lab-frame Gaussian; the archived spherical model is not an
acceptance reference.

## Python environment

Use the shared OPALX Python environment:

```bash
source ~/.venv-h6/bin/activate
python --version
```

Assume Python 3.11 with NumPy, pandas, matplotlib, Pillow, h5py, SciPy,
PyVista/VTK, and imageio. See [`../VENV.md`](../VENV.md) for the full
environment recipe. Set a run-local or temporary `MPLCONFIGDIR` on machines
where the home cache is unavailable.

## Core checks

Run the available static-field component of the broad sandbox regression from
the repository root:

```bash
python sandbox/python/run_sandbox_regressions.py --check fields
```

The complete broad wrapper also expects historical OPALX--IMPACT and pair-stat
outputs that are not present in this checkout; see
[`regression/README.md`](regression/README.md). The maintained focused checks
below are the current BeamBeam acceptance path.

Run the focused manufactured Python checks with:

```bash
~/.venv-h6/bin/python -m unittest \
  sandbox/reduced-order-model/python/test_rigid_two_gaussian_fields.py \
  sandbox/reduced-order-model/python/test_fixed_primary_fromfile.py
```

The end-to-end manufactured OPALX comparison is an opt-in CTest:

```bash
OPALX_RUN_MANUFACTURED_REGRESSION=1 \
ctest --test-dir build_openmp -R TestBeamBeamManufacturedFields --output-on-failure
```

The exact-timed full-CAIN input pipeline is:

```bash
sandbox/cain-opalx-reduced-order-model/run_pipeline.sh
```

The persistent one-, two-, and four-rank timed witness-gather regression is:

```bash
~/.venv-h6/bin/python \
  sandbox/track12particles/opalx/timed/run_witness_gather_mpi_regression.py \
  --opalx build_openmp/src/opalx \
  --mpiexec /opt/homebrew/bin/mpiexec \
  --openmpi-local --force
```

Detailed input generation, Merlin submission, data fetching, and plotting
commands live with the corresponding active workspace. Every Merlin run must
keep its submitted Slurm script, scheduler logs, application logs, timing
files, manifests, and exact inputs below that run's local directory.

## Output and units

Active broad-regression output belongs under `sandbox/regression/`; scientific
case output belongs under its case-specific workspace. Active scripts must not
write into `sandbox/attic/`.

OPALX electric-field data are interpreted as V/m and magnetic fields as T. The
current H5 electric-field unit metadata incorrectly says MV/m; this known issue
is documented in the consolidated note and must be fixed separately.
