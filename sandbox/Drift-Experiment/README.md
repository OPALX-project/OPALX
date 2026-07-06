# Drift Space-Charge Experiment

This directory contains a c0-only OPALX drift validation for the normal
space-charge path.

Run from the repository root:

```sh
~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_drift_validation.py --force
```

The OPALX deck is `spacecharge_drift_30cm.in`: one electron bunch, one 30 cm
drift, one integration step, and no BeamBeam element. The runner enables the
opt-in `OPALX_SC_FIELD_CSV` diagnostic, converts the dumped c0 field rows to
the `sandbox/compare-e-fields` PIC format, and compares against the analytic
isotropic Gaussian field.

Current single-rank result with 400000 macroparticles:

- `relative_vector_L2_vs_analytic = 1.0543e-2`
- strict total-particle scale `3/sqrt(N) = 4.743e-3`
- `passed = False` under that strict criterion

The analytic charge sign is negative for the electron bunch. A positive-charge
sanity comparison gives relative L2 near 2, confirming the sign convention.

Convergence study:

```sh
~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_convergence.py --force --threads 8
```

The 27-run matrix uses `N_GRID = 16,32,64`, `NPPG = 5,8,14`, and seeds
`42,43,44`. Results are collected in `convergence_results/`.
