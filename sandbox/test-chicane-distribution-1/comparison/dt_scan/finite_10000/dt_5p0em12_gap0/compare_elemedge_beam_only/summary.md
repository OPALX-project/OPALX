# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12_gap0/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12_gap0/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

Comparison mode: stat beam quantities only. Temporal monitors, reference-particle fields, and R56 fits are not used.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx  | final_opal   | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------ | ------------ | ------------- |
| energy   | 630     | 1.830358e-10 | 9.967616e-11 | 1.000000e+03 | 1.000000e+03 | -1.830358e-10 |
| dE       | 630     | 2.177503e-11 | 1.000344e-11 | 1.277221e+00 | 1.277221e+00 | -2.144906e-11 |
| rms_x    | 630     | 5.087035e-04 | 2.227305e-04 | 5.686176e-04 | 5.991417e-05 | 5.087035e-04  |
| rms_y    | 630     | 1.418530e-07 | 6.572536e-08 | 5.106278e-05 | 5.120463e-05 | -1.418530e-07 |
| rms_s    | 630     | 9.312326e-07 | 5.507249e-07 | 6.410776e-04 | 6.407976e-04 | 2.800274e-07  |
| rms_px   | 630     | 2.064869e-01 | 8.748455e-02 | 2.191512e-01 | 1.266426e-02 | 2.064869e-01  |
| rms_py   | 630     | 5.936786e-05 | 3.326040e-05 | 0.000000e+00 | 5.936785e-05 | -5.936785e-05 |
| rms_ps   | 630     | 2.254586e-05 | 1.028082e-05 | 2.499460e+00 | 2.499460e+00 | -2.016532e-07 |
| emit_x   | 630     | 1.434750e-04 | 6.915888e-05 | 5.687091e-05 | 7.516880e-07 | 5.611922e-05  |
| emit_y   | 630     | 2.642141e-10 | 1.505977e-10 | 0.000000e+00 | 2.601400e-10 | -2.601400e-10 |
| emit_s   | 630     | 2.419143e-06 | 1.412462e-06 | 1.596176e-03 | 1.595368e-03 | 8.087206e-07  |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_rms_relative_opalx_opal.png`
- `stat_rms_relative_opalx_opal.pdf`
- `stat_rms_momentum_difference_opalx_opal.png`
- `stat_rms_momentum_difference_opalx_opal.pdf`
- `stat_emittance_difference_opalx_opal.png`
- `stat_emittance_difference_opalx_opal.pdf`
- `stat_energy_opalx_opal.png`
- `stat_dE_opalx_opal.png`

