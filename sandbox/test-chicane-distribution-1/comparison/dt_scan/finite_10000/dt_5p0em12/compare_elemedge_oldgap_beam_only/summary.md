# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge_oldgap/test-chicane-distribution-2.stat`
OPAL stat: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

Comparison mode: stat beam quantities only. Temporal monitors, reference-particle fields, and R56 fits are not used.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx  | final_opal   | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------ | ------------ | ------------- |
| energy   | 630     | 4.343974e-10 | 2.578939e-10 | 1.000000e+03 | 1.000000e+03 | -4.343974e-10 |
| dE       | 630     | 2.711986e-11 | 1.612392e-11 | 1.277221e+00 | 1.277221e+00 | -1.914890e-11 |
| rms_x    | 630     | 2.103572e-05 | 8.724355e-06 | 9.394852e-05 | 7.418512e-05 | 1.976340e-05  |
| rms_y    | 630     | 4.534955e-09 | 2.129540e-09 | 5.106761e-05 | 5.107214e-05 | -4.534955e-09 |
| rms_s    | 630     | 1.257058e-07 | 6.106081e-08 | 6.407719e-04 | 6.407592e-04 | 1.269754e-08  |
| rms_px   | 630     | 1.046534e-02 | 4.866209e-03 | 1.874510e-02 | 1.654210e-02 | 2.203006e-03  |
| rms_py   | 630     | 7.999741e-07 | 4.972654e-07 | 6.280913e-06 | 7.080887e-06 | -7.999741e-07 |
| rms_ps   | 630     | 1.500554e-06 | 6.837485e-07 | 2.499460e+00 | 2.499460e+00 | 3.444199e-08  |
| emit_x   | 630     | 6.117939e-06 | 3.249049e-06 | 1.750237e-06 | 1.208974e-06 | 5.412630e-07  |
| emit_y   | 630     | 1.343635e-11 | 4.865328e-12 | 3.051728e-10 | 3.059840e-10 | -8.112302e-13 |
| emit_s   | 630     | 3.151630e-07 | 1.535020e-07 | 1.595289e-03 | 1.595252e-03 | 3.750589e-08  |

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

