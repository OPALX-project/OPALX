# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_2p0em12/opalx_body/test-chicane-distribution-1.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_2p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 15728   | 3.519744e-10 | 1.607786e-10 | 1.000000e+03  | 1.000000e+03  | -3.229843e-10 |
| rms_x    | 15728   | 3.654190e-04 | 1.410487e-04 | 3.978384e-04  | 3.241947e-05  | 3.654190e-04  |
| rms_y    | 15728   | 1.618415e-19 | 8.403225e-20 | 0.000000e+00  | 1.618415e-19  | -1.618415e-19 |
| rms_s    | 15728   | 1.993188e-05 | 1.017847e-05 | 7.128572e-04  | 7.099132e-04  | 2.944043e-06  |
| mean_x   | 15728   | 4.715749e-04 | 2.089535e-04 | 4.721277e-04  | 5.528290e-07  | 4.715749e-04  |
| mean_y   | 15728   | 2.605475e-22 | 1.373042e-22 | 0.000000e+00  | 2.605475e-22  | -2.605475e-22 |
| mean_s   | 15728   | 2.601920e-05 | 1.287347e-05 | 1.406985e-05  | -1.321834e-07 | 1.420204e-05  |
| ref_x    | 15728   | 5.587980e-04 | 2.612048e-04 | -5.596141e-04 | -8.161141e-07 | -5.587980e-04 |
| ref_y    | 15728   | 1.141550e-16 | 5.933382e-17 | 0.000000e+00  | 1.141550e-16  | -1.141550e-16 |
| ref_z    | 15728   | 1.546464e-05 | 7.120402e-06 | 9.408601e+00  | 9.408616e+00  | -1.546464e-05 |
| ref_px   | 15728   | 2.346065e-01 | 1.319009e-01 | -2.346937e-01 | -8.719816e-05 | -2.346065e-01 |
| ref_py   | 15728   | 4.787607e-14 | 2.506375e-14 | 0.000000e+00  | 4.787607e-14  | -4.787607e-14 |
| ref_pz   | 15728   | 2.305811e-02 | 8.998075e-03 | 1.957951e+03  | 1.957951e+03  | -1.406602e-05 |
| By_ref   | 15728   | 2.019104e-01 | 4.159428e-03 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.500269e-01  | 2.476285e-01 | 2.398339e+00  | 1.250134e+00 | 1.252533e+00 | -2.398339e+00 | 1669          | 1677         |
| 2        | 2.750895e+00  | 2.748497e+00 | 2.398339e+00  | 3.750403e+00 | 3.753401e+00 | -2.997924e+00 | 1668          | 1677         |
| 3        | 5.751218e+00  | 5.748819e+00 | 2.398339e+00  | 6.750726e+00 | 6.753724e+00 | -2.997924e+00 | 1668          | 1677         |
| 4        | 8.251487e+00  | 8.249088e+00 | 2.398339e+00  | 9.251594e+00 | 9.253992e+00 | -2.398339e+00 | 1669          | 1677         |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`

