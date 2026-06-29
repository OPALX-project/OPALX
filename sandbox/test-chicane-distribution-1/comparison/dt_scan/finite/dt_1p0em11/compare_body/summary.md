# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_1p0em11/opalx_body/test-chicane-distribution-1.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_1p0em11/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 3146    | 2.150955e-10 | 1.206858e-10 | 1.000000e+03  | 1.000000e+03  | 2.019078e-10  |
| rms_x    | 3146    | 1.245657e-04 | 8.654894e-05 | 3.068620e-04  | 3.889055e-04  | -8.204350e-05 |
| rms_y    | 3146    | 9.939924e-09 | 4.786872e-09 | 5.656854e-05  | 5.657848e-05  | -9.939924e-09 |
| rms_s    | 3146    | 4.397300e-05 | 2.057194e-05 | 6.679930e-04  | 7.119657e-04  | -4.397271e-05 |
| mean_x   | 3146    | 2.110394e-04 | 1.579431e-04 | -1.493017e-04 | 1.786160e-05  | -1.671633e-04 |
| mean_y   | 3146    | 7.828796e-12 | 4.299878e-12 | 0.000000e+00  | -6.529509e-12 | 6.529509e-12  |
| mean_s   | 3146    | 3.427259e-05 | 2.078011e-05 | 3.576276e-05  | 1.421089e-05  | 2.155186e-05  |
| ref_x    | 3146    | 8.102670e-04 | 4.107310e-04 | -2.858750e-14 | -1.303231e-04 | 1.303231e-04  |
| ref_y    | 3146    | 1.142075e-16 | 5.940961e-17 | 0.000000e+00  | 1.142075e-16  | -1.142075e-16 |
| ref_z    | 3146    | 3.744990e-05 | 1.380997e-05 | 9.409782e+00  | 9.409744e+00  | 3.744990e-05  |
| ref_px   | 3146    | 1.709095e+00 | 6.850273e-01 | -1.154632e-14 | -1.603470e-01 | 1.603470e-01  |
| ref_py   | 3146    | 4.791051e-14 | 2.507037e-14 | 0.000000e+00  | 4.791051e-14  | -4.791051e-14 |
| ref_pz   | 3146    | 1.699891e-01 | 3.964428e-02 | 1.957951e+03  | 1.957951e+03  | 6.565838e-06  |
| By_ref   | 3146    | 3.335999e-01 | 1.680054e-02 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.518256e-01  | 2.488277e-01 | 2.997924e+00  | 1.250134e+00 | 1.253132e+00 | -2.997924e+00 | 334           | 336          |
| 2        | 2.752094e+00  | 2.752094e+00 | -3.954392e-08 | 3.750403e+00 | 3.753401e+00 | -2.997924e+00 | 334           | 335          |
| 3        | 5.753016e+00  | 5.753016e+00 | -3.929124e-08 | 6.748327e+00 | 6.757321e+00 | -8.993773e+00 | 333           | 336          |
| 4        | 8.253285e+00  | 8.256283e+00 | -2.997924e+00 | 9.251594e+00 | 9.260588e+00 | -8.993773e+00 | 334           | 336          |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`

