# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 6304    | 4.461072e-10 | 2.564781e-10 | 1.000000e+03  | 1.000000e+03  | -4.239382e-10 |
| rms_x    | 6304    | 2.677172e-06 | 7.269126e-07 | 4.002661e-05  | 3.823993e-05  | 1.786681e-06  |
| rms_y    | 6304    | 1.624486e-19 | 8.429054e-20 | 0.000000e+00  | 1.624486e-19  | -1.624486e-19 |
| rms_s    | 6304    | 6.191313e-07 | 3.369450e-07 | 7.098854e-04  | 7.094118e-04  | 4.736228e-07  |
| mean_x   | 6304    | 1.294850e-05 | 8.288589e-06 | -2.878679e-05 | -1.853241e-05 | -1.025438e-05 |
| mean_y   | 6304    | 5.145727e-21 | 3.098146e-21 | 0.000000e+00  | -4.589797e-21 | 4.589797e-21  |
| mean_s   | 6304    | 1.735285e-06 | 1.134099e-06 | 4.488828e-06  | 3.031978e-06  | 1.456850e-06  |
| ref_x    | 6304    | 4.371580e-04 | 2.277853e-04 | 2.195652e-05  | 1.349320e-05  | 8.463316e-06  |
| ref_y    | 6304    | 1.146076e-16 | 5.950753e-17 | 0.000000e+00  | 1.146076e-16  | -1.146076e-16 |
| ref_z    | 6304    | 1.704346e-05 | 6.832291e-06 | 9.427804e+00  | 9.427787e+00  | 1.704345e-05  |
| ref_px   | 6304    | 8.487936e-01 | 3.407999e-01 | 8.967867e-03  | 2.745284e-03  | 6.222583e-03  |
| ref_py   | 6304    | 4.788343e-14 | 2.512690e-14 | 0.000000e+00  | 4.788343e-14  | -4.788343e-14 |
| ref_pz   | 6304    | 8.450848e-02 | 1.969400e-02 | 1.957951e+03  | 1.957951e+03  | -1.860417e-08 |
| By_ref   | 6304    | 3.336817e-01 | 1.158390e-02 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.488277e-01  | 2.488277e-01 | -2.803313e-12 | 1.251633e+00 | 1.251633e+00 | 5.195844e-11  | 670           | 670          |
| 2        | 2.747598e+00  | 2.749096e+00 | -1.498962e+00 | 3.751902e+00 | 3.753401e+00 | -1.498962e+00 | 671           | 671          |
| 3        | 5.748520e+00  | 5.751518e+00 | -2.997924e+00 | 6.752824e+00 | 6.754323e+00 | -1.498962e+00 | 671           | 670          |
| 4        | 8.248788e+00  | 8.251786e+00 | -2.997924e+00 | 9.251594e+00 | 9.256091e+00 | -4.496886e+00 | 670           | 671          |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`

