# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_5p0em12/opalx_body/test-chicane-distribution-1.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/zero/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 6291    | 3.603873e-11 | 2.009717e-11 | 1.000000e+03  | 1.000000e+03  | -2.899014e-11 |
| rms_x    | 6291    | 3.538087e-04 | 1.614418e-04 | 3.919300e-04  | 3.812136e-05  | 3.538087e-04  |
| rms_y    | 6291    | 1.617745e-19 | 8.405506e-20 | 0.000000e+00  | 1.617745e-19  | -1.617745e-19 |
| rms_s    | 6291    | 9.575353e-06 | 4.687743e-06 | 7.098579e-04  | 7.094118e-04  | 4.461579e-07  |
| mean_x   | 6291    | 5.673612e-04 | 4.052224e-04 | 5.488713e-04  | -1.848987e-05 | 5.673612e-04  |
| mean_y   | 6291    | 5.145727e-21 | 3.094353e-21 | 0.000000e+00  | -4.565892e-21 | 4.565892e-21  |
| mean_s   | 6291    | 8.024248e-05 | 4.377940e-05 | -2.409014e-05 | 3.031979e-06  | -2.712212e-05 |
| ref_x    | 6291    | 4.388390e-04 | 2.690461e-04 | -6.513696e-14 | 1.346588e-05  | -1.346588e-05 |
| ref_y    | 6291    | 1.141311e-16 | 5.934160e-17 | 0.000000e+00  | 1.141311e-16  | -1.141311e-16 |
| ref_z    | 6291    | 4.255750e-05 | 2.114046e-05 | 9.408343e+00  | 9.408300e+00  | 4.255750e-05  |
| ref_px   | 6291    | 6.177809e-01 | 2.529348e-01 | -1.498801e-14 | 2.745284e-03  | -2.745284e-03 |
| ref_py   | 6291    | 4.788343e-14 | 2.505849e-14 | 0.000000e+00  | 4.788343e-14  | -4.788343e-14 |
| ref_pz   | 6291    | 6.154213e-02 | 1.611391e-02 | 1.957951e+03  | 1.957951e+03  | 1.931994e-09  |
| By_ref   | 6291    | 3.336817e-01 | 1.076316e-02 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.503267e-01  | 2.488277e-01 | 1.498962e+00  | 1.250134e+00 | 1.251633e+00 | -1.498962e+00 | 668           | 670          |
| 2        | 2.750595e+00  | 2.749096e+00 | 1.498962e+00  | 3.750403e+00 | 3.753401e+00 | -2.997924e+00 | 668           | 671          |
| 3        | 5.751518e+00  | 5.751518e+00 | -1.625367e-09 | 6.749826e+00 | 6.754323e+00 | -4.496886e+00 | 667           | 670          |
| 4        | 8.251786e+00  | 8.251786e+00 | -1.742606e-09 | 9.251594e+00 | 9.256091e+00 | -4.496886e+00 | 668           | 671          |

## Plots

- `stat_rms_opalx_opal.png`
- `stat_energy_opalx_opal.png`
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`

