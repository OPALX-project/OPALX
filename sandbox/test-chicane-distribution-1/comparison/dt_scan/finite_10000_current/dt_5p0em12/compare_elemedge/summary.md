# test-chicane-distribution-1 OPALX vs OPAL 2022.1

OPALX stat: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000_current/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat`
OPAL stat: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat`

Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.

| quantity | samples | max_abs_diff | rms_diff     | final_opalx   | final_opal    | final_diff    |
| -------- | ------- | ------------ | ------------ | ------------- | ------------- | ------------- |
| energy   | 316     | 1.189164e-10 | 5.766318e-11 | 1.000000e+03  | 1.000000e+03  | 1.500666e-11  |
| dE       | 316     | 2.617284e-11 | 1.562986e-11 | 1.277221e+00  | 1.277221e+00  | -1.847300e-11 |
| rms_x    | 316     | 4.237542e-06 | 1.838001e-06 | 7.652320e-05  | 7.418512e-05  | 2.338081e-06  |
| rms_y    | 316     | 4.664462e-09 | 2.169550e-09 | 5.106748e-05  | 5.107214e-05  | -4.664462e-09 |
| rms_s    | 316     | 3.707864e-07 | 1.744628e-07 | 6.409059e-04  | 6.407592e-04  | 1.467259e-07  |
| rms_px   | 316     | 6.316123e-03 | 2.020383e-03 | 2.285822e-02  | 1.654210e-02  | 6.316123e-03  |
| rms_py   | 316     | 1.274452e-06 | 7.343085e-07 | 5.806435e-06  | 7.080887e-06  | -1.274452e-06 |
| rms_ps   | 316     | 5.422460e-06 | 1.423778e-06 | 2.499460e+00  | 2.499460e+00  | 6.203274e-08  |
| emit_x   | 316     | 2.057667e-06 | 9.507319e-07 | 1.445149e-06  | 1.208974e-06  | 2.361755e-07  |
| emit_y   | 316     | 2.709006e-11 | 1.554698e-11 | 2.806749e-10  | 3.059840e-10  | -2.530908e-11 |
| emit_s   | 316     | 8.937009e-07 | 4.217367e-07 | 1.595636e-03  | 1.595252e-03  | 3.840578e-07  |
| mean_x   | 316     | 2.011814e-05 | 9.709903e-06 | -2.367641e-05 | -2.570576e-05 | 2.029348e-06  |
| mean_y   | 316     | 7.647424e-18 | 4.487100e-18 | 3.451829e-21  | -7.641038e-18 | 7.644490e-18  |
| mean_s   | 316     | 1.333839e-06 | 7.091847e-07 | 2.373145e-06  | 3.174517e-06  | -8.013721e-07 |
| ref_x    | 316     | 3.132790e-04 | 1.679236e-04 | 1.393298e-05  | 1.349320e-05  | 4.397780e-07  |
| ref_y    | 316     | 1.146076e-16 | 5.987637e-17 | 0.000000e+00  | 1.146076e-16  | -1.146076e-16 |
| ref_z    | 316     | 1.468031e-05 | 5.805016e-06 | 9.427801e+00  | 9.427787e+00  | 1.468030e-05  |
| ref_px   | 316     | 5.854791e-01 | 2.330698e-01 | 1.025613e-02  | 2.745284e-03  | 7.510850e-03  |
| ref_py   | 316     | 4.788343e-14 | 2.527866e-14 | 0.000000e+00  | 4.788343e-14  | -4.788343e-14 |
| ref_pz   | 316     | 5.720230e-02 | 1.345837e-02 | 1.957951e+03  | 1.957951e+03  | -2.512115e-08 |
| By_ref   | 316     | 3.246396e-08 | 1.826267e-09 | 0.000000e+00  | 0.000000e+00  | 0.000000e+00  |

## By Field Intervals

Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`.

| interval | opalx_start_m | opal_start_m | start_diff_mm | opalx_end_m  | opal_end_m   | end_diff_mm   | opalx_samples | opal_samples |
| -------- | ------------- | ------------ | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ |
| 1        | 2.698132e-01  | 2.548236e-01 | 1.498962e+01  | 1.229149e+00 | 1.244139e+00 | -1.498962e+01 | 33            | 67           |
| 2        | 2.758090e+00  | 2.758090e+00 | -4.123812e-09 | 3.747405e+00 | 3.747405e+00 | -4.300116e-09 | 34            | 67           |
| 3        | 5.756014e+00  | 5.756014e+00 | -6.745715e-09 | 6.745329e+00 | 6.745329e+00 | -7.085887e-09 | 34            | 67           |
| 4        | 8.274271e+00  | 8.259281e+00 | 1.498962e+01  | 9.233606e+00 | 9.248596e+00 | -1.498962e+01 | 33            | 67           |

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
- `stat_reference_xz_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.png`
- `stat_bfield_ref_opalx_opal.pdf`

