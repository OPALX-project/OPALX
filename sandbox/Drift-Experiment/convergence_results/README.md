# Drift Space-Charge Convergence Study

N_GRID values: `[16, 32, 64]`
NPPG values: `[5, 8, 14]`
Seeds: `[42, 43, 44]`

Files:

- `convergence_summary.csv`: one row per OPALX run.
- `convergence_grouped_summary.csv`: grouped mean/std by `(N_GRID, NPPG)`.
- `relative_l2_vs_grid.png`, `relative_l2_vs_nppg.png`,
  `relative_l2_vs_particles.png`: aggregate convergence plots.
- Per-case directories contain `opalx.out`, the rendered input deck,
  `metrics.csv`, and a deterministic sampled comparison CSV.

Raw field CSVs are deleted after each successful case unless `--keep-raw` is used.
