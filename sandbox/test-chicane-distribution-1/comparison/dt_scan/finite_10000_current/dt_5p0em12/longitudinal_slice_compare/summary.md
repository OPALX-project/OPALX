# Longitudinal Slice Comparison

Distribution: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000_current/dt_5p0em12/opalx_elemedge/test-chicane-distribution-1_distribution.txt`
OPALX monitor: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000_current/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2_exit.h5`
OPAL monitor: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal_exit.h5`
OPALX selected step: `Step#0`, SPOS `9.251669234960e+00 m`
OPAL selected step: `Step#0`, SPOS `9.255537067310e+00 m`

Rows are grouped by initial `z0`; each slice contains the full transverse and momentum grid.

## Key Slice Findings

- Largest absolute mean horizontal kick difference: `z0 = 1.000000000000e-03 m`, `mean(OPALX-OPAL px) = 2.802699476773e-01`.
- Largest absolute slice `sigma_x` difference: `z0 = 3.333333333333e-04 m`, `sigma_x(OPALX)-sigma_x(OPAL) = -5.298283720719e-06 m`.
- Linear fit of slice `mean(OPALX-OPAL px)` versus `z0`: slope `7.952424721976e+01 1/m`, intercept `2.116600675612e-01`, `R^2 = 0.824036`.

## Horizontal Variance Decomposition

| code | total sigma_x [m] | within-slice sigma [m] | slice-centroid sigma [m] | centroid variance fraction |
| --- | --- | --- | --- | --- |
| OPALX | 7.525249562e-05 | 6.004780032e-05 | 4.535636421e-05 | 3.632742825e-01 |
| OPAL 2022.1 | 7.392105877e-05 | 6.031931605e-05 | 4.273058671e-05 | 3.341499147e-01 |

## Compact Table

| z0 [mm] | sigma_x OPALX [m] | sigma_x OPAL [m] | delta sigma_x [m] | mean delta px | rms delta px | mean delta z [m] |
| --- | --- | --- | --- | --- | --- | --- |
| -1.000000000e+00 | 5.919442812e-05 | 6.272356183e-05 | -3.529133704e-06 | 9.475476163e-02 | 5.129989042e-03 | -2.331096334e-06 |
| -7.777777778e-01 | 6.178165385e-05 | 5.904007005e-05 | 2.741583805e-06 | 1.295599887e-01 | 7.112801894e-03 | -2.266787854e-06 |
| -5.555555556e-01 | 5.984219231e-05 | 5.663342540e-05 | 3.208766915e-06 | 1.824381613e-01 | 1.291339904e-02 | 1.373403037e-07 |
| -3.333333333e-01 | 5.700335868e-05 | 5.764521124e-05 | -6.418525618e-07 | 2.158426200e-01 | 7.814278476e-03 | 1.179733544e-06 |
| -1.111111111e-01 | 6.007946879e-05 | 6.107606942e-05 | -9.966006255e-07 | 2.258856935e-01 | 6.510092938e-03 | 4.339947238e-07 |
| 1.111111111e-01 | 6.251824142e-05 | 6.414479108e-05 | -1.626549654e-06 | 2.474626102e-01 | 8.086535350e-03 | -1.639254648e-07 |
| 3.333333333e-01 | 5.923312500e-05 | 6.453140872e-05 | -5.298283721e-06 | 2.530767505e-01 | 2.249713901e-03 | -1.346024164e-06 |
| 5.555555556e-01 | 6.000259770e-05 | 6.203257430e-05 | -2.029976598e-06 | 2.395435931e-01 | 6.383784744e-03 | -2.497226587e-06 |
| 7.777777778e-01 | 6.181847046e-05 | 5.787070107e-05 | 3.947769393e-06 | 2.477665491e-01 | 3.253169730e-03 | -1.807396774e-06 |
| 1.000000000e+00 | 5.880040879e-05 | 5.681291270e-05 | 1.987496090e-06 | 2.802699477e-01 | 8.099340764e-03 | 6.795998005e-07 |

## Plots

- `slice_mean_x.png`
- `slice_rms_x.png`
- `slice_rms_px.png`
- `slice_rms_z.png`
- `slice_max_abs_px_diff.png`
