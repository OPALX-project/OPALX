# Primary-only field diagnostic

This case samples the physical primary field on a passive `101 x 201` x-z probe sheet. The probes
do not deposit charge and are limited to `+/-3 sigma_x`, so they do not intentionally enlarge the
three-sigma primary transverse support. `EBDUMP=TRUE` stores the gathered field components in the
probe H5 file. The BeamBeam copied source is disabled.

Generate and run the `256 x 256 x 128` case from the repository root:

```bash
~/.venv-h6/bin/python sandbox/reduced-order-model/python/make_opalx_field_cases.py \
  --config sandbox/track-e-p-new/field-diagnostics/primary_field_parameters.json \
  --template sandbox/track-e-p-new/field-diagnostics/primary_field_snapshot.in.template \
  --output-dir sandbox/track-e-p-new/field-diagnostics/outputs \
  --run --opalx "$PWD/build_openmp/src/opalx" --omp-threads 4 \
  --primary-macroparticles 400000 \
  --primary-particle-file \
    "$PWD/sandbox/track12particles/opalx/timed/input/primary_fixed.fromfile" \
  --nxy 256 --nz 128 --separations-sigma-z 0 --force
```

The generator invokes OPALX with `--info 4`. Analyze and plot the probe fields against the
component-wise three-sigma-truncated manufactured Gaussian with:

```bash
MPLBACKEND=Agg ~/.venv-h6/bin/python \
  sandbox/track-e-p-new/field-diagnostics/plot_primary_field_diagnostic.py
```

Generated H5, STAT, log, CSV, JSON, PNG, and PDF artifacts remain below the ignored `outputs/`
directory.

With the temporary device-side stage counters enabled, search the log with:

```bash
rg "field finiteness|witness finiteness" \
  sandbox/track-e-p-new/field-diagnostics/outputs/0sigma/rigid_fields_0sigma.log
```

The integrated case reports every rest-frame E mesh component as non-finite immediately after
`runSolver(true)`. The otherwise identical standard-Green control reports zero non-finite mesh and
probe components. This does not identify the integrated solver as defective: the old validated
Track-12 case also used it successfully. The diagnostic demonstrates that the redesigned
particle-bounds/fixed-z field geometry is outside the previously validated numerical regime.

## Transverse-domain control matrix

With the 20 mm longitudinal window fixed, particle-tight transverse domains produced non-finite
integrated-Green fields for all three tested meshes: `64 x 64 x 128`, `128 x 128 x 128`, and
`256 x 256 x 128`. A second matrix kept the transverse lab-frame cell width near `0.394 um` by
using domains of `25 um`, `50 um`, and `100 um` with 64, 128, and 256 transverse cells,
respectively. All three cases were finite. Relative L2 errors against the three-sigma-truncated
manufactured Gaussian were, respectively, `1.9662e-2`, `1.9562e-2`, and `1.9515e-2` for E, with
the same values to reported precision for B.

The results are under `outputs/padded-n64-w25um`, `outputs/padded-n128-w50um`, and
`outputs/padded-n256-w100um`. These directories are ignored generated artifacts. The current
uncommitted code used the physical `APERTURE` as a convenient way to impose the diagnostic floors;
that is deliberately not the intended interface. Physical acceptance and numerical field-domain
sizing must remain separate.

The off-diagonal cases are under `outputs/ratio-matrix-*`. With `Delta z' = 75.66 mm` in the
primary rest frame, the `12 um / 32` and `25 um / 64` cases were finite at aspect ratios
`1.95e5` and `1.91e5`. Holding the 25 um domain fixed and increasing to 128 transverse cells raised
the aspect ratio to `3.84e5` and made every probe E and B component non-finite. The `50 um / 64`
case was finite at `9.53e4`, although its coarser `0.794 um` transverse cells increased the relative
L2 field error to 4.38%. This matrix identifies transformed cell aspect ratio as the cause of the
non-finite solve, while also showing that a separate resolution criterion is required for accuracy.
