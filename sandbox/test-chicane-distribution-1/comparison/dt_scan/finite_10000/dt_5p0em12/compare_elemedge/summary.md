# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 630     | 3.469722e-10 | 2.176399e-10 | 1.000000e+03  | 1.000000e+03  | -3.420837e-10 |
| rms_x    | 630     | 1.747598e-05 | 4.647197e-06 | 9.166109e-05  | 7.418512e-05  | 1.747598e-05  |
| rms_y    | 630     | 4.463877e-09 | 2.095704e-09 | 5.106768e-05  | 5.107214e-05  | -4.463877e-09 |
| rms_s    | 630     | 3.347392e-07 | 1.640814e-07 | 6.408123e-04  | 6.407592e-04  | 5.311967e-08  |
| rms_px   | 630     | 1.349527e-02 | 3.133896e-03 | 3.003737e-02  | 1.654210e-02  | 1.349527e-02  |
| rms_py   | 630     | 7.775251e-07 | 4.837915e-07 | 6.303362e-06  | 7.080887e-06  | -7.775251e-07 |
| rms_ps   | 630     | 4.541884e-06 | 1.249647e-06 | 2.499461e+00  | 2.499460e+00  | 9.413745e-08  |
| emit_x   | 630     | 2.583172e-06 | 1.378382e-06 | 1.981628e-06  | 1.208974e-06  | 7.726538e-07  |
| emit_y   | 630     | 1.343635e-11 | 5.201054e-12 | 3.041756e-10  | 3.059840e-10  | -1.808427e-12 |
| emit_s   | 630     | 8.458359e-07 | 4.182756e-07 | 1.595372e-03  | 1.595252e-03  | 1.206216e-07  |
| mean_x   | 630     | 1.318672e-05 | 8.452783e-06 | -3.889248e-05 | -2.570576e-05 | -1.318672e-05 |
| mean_y   | 630     | 7.647329e-18 | 4.478532e-18 | 3.317659e-21  | -7.641038e-18 | 7.644356e-18  |
| mean_s   | 630     | 1.677794e-06 | 1.117151e-06 | 4.444794e-06  | 3.174517e-06  | 1.270277e-06  |
| ref_x    | 630     | 4.371114e-04 | 2.278567e-04 | 2.195652e-05  | 1.349320e-05  | 8.463316e-06  |
| ref_y    | 630     | 1.146076e-16 | 5.970998e-17 | 0.000000e+00  | 1.146076e-16  | -1.146076e-16 |
| ref_z    | 630     | 1.704346e-05 | 6.869742e-06 | 9.427804e+00  | 9.427787e+00  | 1.704345e-05  |
| ref_px   | 630     | 8.487936e-01 | 3.417715e-01 | 8.967867e-03  | 2.745284e-03  | 6.222583e-03  |
| ref_py   | 630     | 4.788343e-14 | 2.521073e-14 | 0.000000e+00  | 4.788343e-14  | -4.788343e-14 |
| ref_pz   | 630     | 8.413723e-02 | 1.973397e-02 | 1.957951e+03  | 1.957951e+03  | -1.879312e-08 |
| By_ref   | 630     | 3.239108e-08 | 1.290603e-09 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.548236e-01  | 2.548236e-01 | 4.518608e-11  | 1.244139e+00 | 1.244139e+00 | -1.199041e-11 | 67            | 67           |
| 2        | 2.758090e+00  | 2.758090e+00 | -4.473755e-09 | 3.747405e+00 | 3.747405e+00 | -4.666045e-09 | 67            | 67           |
| 3        | 5.756014e+00  | 5.756014e+00 | -1.002043e-08 | 6.745329e+00 | 6.745329e+00 | -1.030198e-08 | 67            | 67           |
| 4        | 8.259281e+00  | 8.259281e+00 | -1.662848e-08 | 9.248596e+00 | 9.248596e+00 | -1.740474e-08 | 67            | 67           |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_rms_relative_opalx_opal.png`
- `stat_rms_relative_opalx_opal.pdf`
- `stat_rms_momentum_difference_opalx_opal.png`
- `stat_rms_momentum_difference_opalx_opal.pdf`
- `stat_emittance_difference_opalx_opal.png`
- `stat_emittance_difference_opalx_opal.pdf`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.pdf`

## R56 Fits

- OPALX: $R_{56}=-4.45195920e-02$ m, $\Delta R_{56}=-1.19e-03$ m
- OPAL: $R_{56}=-4.44786273e-02$ m, $\Delta R_{56}=-1.15e-03$ m

