# Drift Space-Charge Convergence Study

N_GRID values: `[16]`
NPPG values: `[5]`
Seeds: `[42]`

Files:

- `convergence_summary.csv`: one row per OPALX run.
- `convergence_grouped_summary.csv`: grouped mean/std by `(N_GRID, NPPG)`.
- `relative_l2_vs_grid.png`, `relative_l2_vs_nppg.png`,
  `relative_l2_vs_particles.png`: aggregate convergence plots.
- Per-case directories contain `opalx.out`, the rendered input deck,
  `metrics.csv`, and a deterministic sampled comparison CSV.

Raw particle field data are read from the OPALX H5 dump.
