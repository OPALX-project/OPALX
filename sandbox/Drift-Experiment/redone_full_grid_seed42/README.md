# Drift Space-Charge Convergence Study

N_GRID values: `[16, 32, 64]`
NPPG values: `[5, 8, 14]`
Seeds: `[42]`

Files:

- `convergence_summary.csv`: one row per OPALX run.
- `convergence_grouped_summary.csv`: grouped mean/std by `(N_GRID, NPPG)`.
- `relative_l2_vs_grid.png`, `relative_l2_vs_nppg.png`,
  `relative_l2_vs_particles.png`: aggregate convergence plots.
- `radial_field_vs_analytic_all_cases.png`,
  `signed_radial_field_vs_analytic_all_cases.png`, and
  `radial_field_relative_error_all_cases.png`: radial profile grids.
- Per-case directories contain `opalx.out`, the rendered input deck,
  `metrics.csv`, `radial_profile.csv`, radial PNGs, and a deterministic
  sampled comparison CSV.

Raw H5 particle dumps are deleted after each successful case unless
`--keep-raw` is used.
