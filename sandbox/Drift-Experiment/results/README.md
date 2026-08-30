# c0-only drift space-charge validation

Input deck: `sandbox/Drift-Experiment/spacecharge_drift_30cm.in`

The OPALX run contains only container 0, a 30 cm drift, one integration step,
and an isotropic Gaussian electron bunch. `OPALX_SC_FIELD_CSV` dumps the
normal space-charge field after `computeSelfFields()` and transformation back
to the reference frame. The dumped particle fields are compared with
`sandbox/compare-e-fields/compare_gaussian_pic_fields.py`. The analytic
Gaussian is centered on the dumped c0 centroid in the same reference frame.

Files:

- `raw/spacecharge_fields.csv`: OPALX diagnostic field dump.
- `inputs/pic_fields.csv`: comparator input adapted from the OPALX dump.
- `comparison/`: output from `compare_gaussian_pic_fields.py`.
- `validation_summary.csv`: pass/fail summary.

Success criterion used by this runner:

`relative_vector_L2_vs_analytic <= 3 / sqrt(N)`.

Observed relative L2: `1.054309982766e-02`
Threshold: `4.743416490253e-03`
Passed: `False`
