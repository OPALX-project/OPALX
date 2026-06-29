# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_2p0em12/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_2p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 15761   | 4.610001e-10 | 2.912974e-10 | 1.000000e+03  | 1.000000e+03  | 3.969944e-10  |
| rms_x    | 15761   | 8.337766e-06 | 1.059157e-06 | 3.412711e-05  | 3.255103e-05  | 1.576081e-06  |
| rms_y    | 15761   | 1.625265e-19 | 8.427167e-20 | 0.000000e+00  | 1.625265e-19  | -1.625265e-19 |
| rms_s    | 15761   | 8.273099e-07 | 4.157305e-07 | 7.098940e-04  | 7.099132e-04  | -1.913257e-08 |
| mean_x   | 15761   | 8.503836e-07 | 5.477620e-07 | 7.506535e-07  | 5.530953e-07  | 1.975582e-07  |
| mean_y   | 15761   | 2.615066e-22 | 1.376795e-22 | 0.000000e+00  | 2.615066e-22  | -2.615066e-22 |
| mean_s   | 15761   | 1.218156e-07 | 7.204336e-08 | -2.212387e-07 | -1.321838e-07 | -8.905488e-08 |
| ref_x    | 15761   | 1.454483e-04 | 7.608342e-05 | -1.081442e-06 | -8.169952e-07 | -2.644471e-07 |
| ref_y    | 15761   | 1.146388e-16 | 5.950239e-17 | 0.000000e+00  | 1.146388e-16  | -1.146388e-16 |
| ref_z    | 15761   | 4.938664e-06 | 2.118983e-06 | 9.428407e+00  | 9.428402e+00  | 4.938664e-06  |
| ref_px   | 15761   | 2.844710e-01 | 1.150331e-01 | -2.963946e-04 | -8.719816e-05 | -2.091964e-04 |
| ref_py   | 15761   | 4.787607e-14 | 2.513315e-14 | 0.000000e+00  | 4.787607e-14  | -4.787607e-14 |
| ref_pz   | 15761   | 2.834759e-02 | 6.644501e-03 | 1.957951e+03  | 1.957951e+03  | -1.887202e-11 |
| By_ref   | 15761   | 2.612373e-01 | 5.518510e-03 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.476285e-01  | 2.476285e-01 | 8.326673e-14  | 1.252533e+00 | 1.252533e+00 | 1.887379e-11  | 1677          | 1677         |
| 2        | 2.747897e+00  | 2.748497e+00 | -5.995848e-01 | 3.752801e+00 | 3.753401e+00 | -5.995848e-01 | 1677          | 1677         |
| 3        | 5.747620e+00  | 5.748819e+00 | -1.199170e+00 | 6.752524e+00 | 6.753724e+00 | -1.199170e+00 | 1677          | 1677         |
| 4        | 8.247889e+00  | 8.249088e+00 | -1.199170e+00 | 9.252793e+00 | 9.253992e+00 | -1.199170e+00 | 1677          | 1677         |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`

