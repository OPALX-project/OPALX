# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

Comparison mode: stat beam quantities only. Temporal monitors, reference-particle fields, and R56 fits are not used.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx  | final_opal   | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------ | ------------ | ------------- |
| energy   | 630     | 3.469722e-10 | 2.178407e-10 | 1.000000e+03 | 1.000000e+03 | -3.420837e-10 |
| dE       | 630     | 2.888090e-11 | 1.747829e-11 | 1.277221e+00 | 1.277221e+00 | -2.133183e-11 |
| rms_x    | 630     | 1.747598e-05 | 4.647197e-06 | 9.166109e-05 | 7.418512e-05 | 1.747598e-05  |
| rms_y    | 630     | 4.463877e-09 | 2.095704e-09 | 5.106768e-05 | 5.107214e-05 | -4.463877e-09 |
| rms_s    | 630     | 3.347392e-07 | 1.640814e-07 | 6.408123e-04 | 6.407592e-04 | 5.311967e-08  |
| rms_px   | 630     | 1.349527e-02 | 3.133896e-03 | 3.003737e-02 | 1.654210e-02 | 1.349527e-02  |
| rms_py   | 630     | 7.775251e-07 | 4.837915e-07 | 6.303362e-06 | 7.080887e-06 | -7.775251e-07 |
| rms_ps   | 630     | 4.541884e-06 | 1.249647e-06 | 2.499461e+00 | 2.499460e+00 | 9.413754e-08  |
| emit_x   | 630     | 2.583172e-06 | 1.378382e-06 | 1.981628e-06 | 1.208974e-06 | 7.726538e-07  |
| emit_y   | 630     | 1.343635e-11 | 5.201054e-12 | 3.041756e-10 | 3.059840e-10 | -1.808427e-12 |
| emit_s   | 630     | 8.458359e-07 | 4.182756e-07 | 1.595372e-03 | 1.595252e-03 | 1.206216e-07  |

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

