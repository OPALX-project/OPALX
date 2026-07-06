Field comparison output
=======================

Files:
- comparison_metrics.csv: one-row global norm and error summary.
- particle_field_comparison.csv: per-particle PIC, analytic, and error columns.
- radial_bin_statistics.csv: radial means/stds, if radial bins were possible.
- field_projection_xy.png: x-y scatter projections for PIC, analytic, and errors.
- component_scatter_pic_vs_analytic.png: Ex, Ey, Ez, and |E| PIC-vs-analytic scatter plots.
- error_histograms.png: histograms for field magnitudes and errors.
- radial_profile_Eabs.png and radial_profile_errors.png: binned radial profiles.
- quiver_3d_pic_and_analytic.png: sampled 3D vector visualization, if enabled.

Run settings:
- pic: sandbox/Drift-Experiment/results/inputs/pic_fields.csv
- analytic: None
- match: row
- sigma: 0.000985
- charge: -1e-09
- center: [2.723169e-19, 6.487866e-20, 0.15]
- eps0: 8.8541878128e-12
- coord_tol: 1e-12
- rel_floor: 1e-30
- outdir: sandbox/Drift-Experiment/results/comparison_sigma_sample
- max_plot_points: 50000
- radial_bins: 80
- make_quiver: False
- quiver_points: 1000
- quiver_scale: None
- seed: 12345