# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 6304    | 3.469722e-10 | 2.171243e-10 | 1.000000e+03  | 1.000000e+03  | -3.420837e-10 |
| rms_x    | 6304    | 1.747598e-05 | 4.590712e-06 | 9.166109e-05  | 7.418512e-05  | 1.747598e-05  |
| rms_y    | 6304    | 4.463877e-09 | 2.087104e-09 | 5.106768e-05  | 5.107214e-05  | -4.463877e-09 |
| rms_s    | 6304    | 3.360315e-07 | 1.640138e-07 | 6.408123e-04  | 6.407592e-04  | 5.311967e-08  |
| rms_px   | 6304    | 1.032608e-01 | 4.096185e-03 | 3.003737e-02  | 1.654210e-02  | 1.349527e-02  |
| rms_py   | 6304    | 8.437728e-07 | 4.821854e-07 | 6.303362e-06  | 7.080887e-06  | -7.775251e-07 |
| rms_ps   | 6304    | 6.571014e-06 | 1.254045e-06 | 2.499461e+00  | 2.499460e+00  | 9.413760e-08  |
| emit_x   | 6304    | 3.317922e-05 | 1.533026e-06 | 1.981628e-06  | 1.208974e-06  | 7.726538e-07  |
| emit_y   | 6304    | 7.410045e-11 | 5.504476e-12 | 3.041756e-10  | 3.059840e-10  | -1.808427e-12 |
| emit_s   | 6304    | 8.485023e-07 | 4.181158e-07 | 1.595372e-03  | 1.595252e-03  | 1.206216e-07  |
| mean_x   | 6304    | 1.318672e-05 | 8.433038e-06 | -3.889248e-05 | -2.570576e-05 | -1.318672e-05 |
| mean_y   | 6304    | 7.648397e-18 | 4.465139e-18 | 1.708974e-21  | -7.641038e-18 | 7.642747e-18  |
| mean_s   | 6304    | 1.677830e-06 | 1.115590e-06 | 4.444794e-06  | 3.174517e-06  | 1.270277e-06  |
| ref_x    | 6304    | 4.371580e-04 | 2.277853e-04 | 2.195652e-05  | 1.349320e-05  | 8.463316e-06  |
| ref_y    | 6304    | 1.146076e-16 | 5.950753e-17 | 0.000000e+00  | 1.146076e-16  | -1.146076e-16 |
| ref_z    | 6304    | 1.704346e-05 | 6.832291e-06 | 9.427804e+00  | 9.427787e+00  | 1.704345e-05  |
| ref_px   | 6304    | 8.487936e-01 | 3.407999e-01 | 8.967867e-03  | 2.745284e-03  | 6.222583e-03  |
| ref_py   | 6304    | 4.788343e-14 | 2.512690e-14 | 0.000000e+00  | 4.788343e-14  | -4.788343e-14 |
| ref_pz   | 6304    | 8.450848e-02 | 1.969400e-02 | 1.957951e+03  | 1.957951e+03  | -1.879312e-08 |
| By_ref   | 6304    | 3.336817e-01 | 1.158390e-02 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.488277e-01  | 2.488277e-01 | -2.803313e-12 | 1.251633e+00 | 1.251633e+00 | 4.884981e-11  | 670           | 670          |
| 2        | 2.747598e+00  | 2.749096e+00 | -1.498962e+00 | 3.751902e+00 | 3.753401e+00 | -1.498962e+00 | 671           | 671          |
| 3        | 5.748520e+00  | 5.751518e+00 | -2.997924e+00 | 6.752824e+00 | 6.754323e+00 | -1.498962e+00 | 671           | 670          |
| 4        | 8.248788e+00  | 8.251786e+00 | -2.997924e+00 | 9.251594e+00 | 9.256091e+00 | -4.496886e+00 | 670           | 671          |

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

