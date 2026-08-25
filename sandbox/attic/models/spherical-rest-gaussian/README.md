# Analytic Gaussian Witness Model

This sandbox reproduces the timing and witness setup of the clean OPALX
BeamBeam diagnostic run:

```text
sandbox/Drift-Experiment/one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs_extended_bb_nobound
```

The analytic model uses the rigid boosted triaxial Gaussian field evaluator in
`sandbox/track12particles/track12particles.py`.  It tracks the six electron and
six positron FROMFILE witnesses with the same initial positions, momentum,
species charge sign, witness bunch charge magnitude, c0 bunch charge, c0
sigmas, c0 timing, IP location, and c0 retirement time as the OPALX input copy
stored here as `spacecharge_drift_withness.opalx-reference.in`.

Run:

```bash
/Users/adelmann/.venv-h6/bin/python sandbox/analytic-model/run_analytic_witness.py
```

Main products:

- `setup_summary.json`
- `witness_initial_conditions.csv`
- `analytic_witness_trajectory.csv`
- `opalx_witness_trajectory_sampled.csv`
- `analytic_vs_opalx_samples.csv`
- `analytic_vs_opalx_summary.csv`
- `analytic_vs_opalx_side_by_side_x_vs_s.png`
- `analytic_vs_opalx_side_by_side_x_vs_t.png`
- `analytic_vs_opalx_overlay_x_vs_s.png`
- `analytic_vs_opalx_overlay_x_vs_t.png`
- `analytic_field_timing.png`
- `analytic_minus_opalx_differences.png`

By default the analytic source charge is the physical c0 H5 charge,
`-1.25e10 e`.  For sign-convention diagnostics only, pass
`--source-charge-scale -1` to flip the analytic source field.
