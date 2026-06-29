# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_1p0em11/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_1p0em11/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 3151    | 1.920171e-10 | 1.368542e-10 | 1.000000e+03  | 1.000000e+03  | 1.840590e-10  |
| rms_x    | 3151    | 9.040743e-05 | 4.291229e-05 | 4.715497e-04  | 3.814291e-04  | 9.012067e-05  |
| rms_y    | 3151    | 1.698081e-19 | 8.907469e-20 | 0.000000e+00  | 1.698081e-19  | -1.698081e-19 |
| rms_s    | 3151    | 8.467323e-06 | 4.620693e-06 | 7.092943e-04  | 7.119672e-04  | -2.672928e-06 |
| mean_x   | 3151    | 4.756502e-04 | 2.770853e-04 | 4.962627e-04  | 2.061246e-05  | 4.756502e-04  |
| mean_y   | 3151    | 3.384308e-20 | 2.128767e-20 | 0.000000e+00  | -3.384308e-20 | 3.384308e-20  |
| mean_s   | 3151    | 4.421954e-05 | 2.248776e-05 | 1.498626e-05  | 1.422430e-05  | 7.619643e-07  |
| ref_x    | 3151    | 1.336724e-03 | 7.778377e-04 | -5.641787e-04 | -1.315506e-04 | -4.326281e-04 |
| ref_y    | 3151    | 1.145743e-16 | 5.953720e-17 | 0.000000e+00  | 1.145743e-16  | -1.145743e-16 |
| ref_z    | 3151    | 5.542753e-05 | 2.922418e-05 | 9.424789e+00  | 9.424734e+00  | 5.542698e-05  |
| ref_px   | 3151    | 1.725948e+00 | 7.314865e-01 | -2.193414e-01 | -1.603470e-01 | -5.899446e-02 |
| ref_py   | 3151    | 4.791051e-14 | 2.512306e-14 | 0.000000e+00  | 4.791051e-14  | -4.791051e-14 |
| ref_pz   | 3151    | 1.716580e-01 | 4.284297e-02 | 1.957951e+03  | 1.957951e+03  | -5.720138e-06 |
| By_ref   | 3151    | 3.335999e-01 | 1.834222e-02 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.488277e-01  | 2.488277e-01 | -9.159340e-13 | 1.250134e+00 | 1.253132e+00 | -2.997924e+00 | 335           | 336          |
| 2        | 2.749096e+00  | 2.752094e+00 | -2.997924e+00 | 3.750403e+00 | 3.753401e+00 | -2.997924e+00 | 335           | 335          |
| 3        | 5.750019e+00  | 5.753016e+00 | -2.997924e+00 | 6.751325e+00 | 6.757321e+00 | -5.995848e+00 | 335           | 336          |
| 4        | 8.250287e+00  | 8.256283e+00 | -5.995848e+00 | 9.251594e+00 | 9.260588e+00 | -8.993773e+00 | 335           | 336          |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`

