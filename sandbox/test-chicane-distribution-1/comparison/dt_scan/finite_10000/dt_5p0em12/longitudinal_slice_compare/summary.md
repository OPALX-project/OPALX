# Longitudinal Slice Comparison

Distribution: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-1_distribution.txt`
OPALX monitor: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2_exit.h5`
OPAL monitor: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal_exit.h5`
OPALX selected step: `Step#0`, SPOS `9.250753364361e+00 m`
OPAL selected step: `Step#0`, SPOS `9.255537067310e+00 m`

Rows are grouped by initial `z0`; each slice contains the full transverse and momentum grid.

## Key Slice Findings

- Largest absolute mean horizontal kick difference: `z0 = -1.000000000000e-03 m`, `mean(OPALX-OPAL px) = -2.286316100580e-01`.
- Largest absolute slice `sigma_x` difference: `z0 = -5.555555555556e-04 m`, `sigma_x(OPALX)-sigma_x(OPAL) = 4.838562498123e-06 m`.
- Linear fit of slice `mean(OPALX-OPAL px)` versus `z0`: slope `1.857330071435e+02 1/m`, intercept `-4.398443558713e-02`, `R^2 = 0.974542`.

## Horizontal Variance Decomposition

| code | total sigma_x [m] | within-slice sigma [m] | slice-centroid sigma [m] | centroid variance fraction |
| --- | --- | --- | --- | --- |
| OPALX | 8.954720529e-05 | 6.069065463e-05 | 6.584334755e-05 | 5.406543888e-01 |
| OPAL 2022.1 | 7.392105877e-05 | 6.031931605e-05 | 4.273058671e-05 | 3.341499147e-01 |

## Compact Table

| z0 [mm] | sigma_x OPALX [m] | sigma_x OPAL [m] | delta sigma_x [m] | mean delta px | rms delta px | mean delta z [m] |
| --- | --- | --- | --- | --- | --- | --- |
| -1.000000000e+00 | 5.886796023e-05 | 6.272356183e-05 | -3.855601593e-06 | -2.286316101e-01 | 1.977038070e-03 | -7.493151534e-04 |
| -7.777777778e-01 | 6.169227309e-05 | 5.904007005e-05 | 2.652203038e-06 | -2.148770443e-01 | 1.148182409e-02 | -7.495055905e-04 |
| -5.555555556e-01 | 6.147198789e-05 | 5.663342540e-05 | 4.838562498e-06 | -1.485790116e-01 | 2.344184948e-02 | -7.467957410e-04 |
| -3.333333333e-01 | 5.846707809e-05 | 5.764521124e-05 | 8.218668474e-07 | -6.715866091e-02 | 1.467628065e-02 | -7.446181175e-04 |
| -1.111111111e-01 | 5.849000824e-05 | 6.107606942e-05 | -2.586061175e-06 | -5.085778106e-02 | 3.084539852e-03 | -7.475743829e-04 |
| 1.111111111e-01 | 6.341484924e-05 | 6.414479108e-05 | -7.299418389e-07 | -3.447704695e-02 | 1.093861194e-02 | -7.503371964e-04 |
| 3.333333333e-01 | 6.229422210e-05 | 6.453140872e-05 | -2.237186622e-06 | 1.975013997e-02 | 9.410759539e-03 | -7.494119796e-04 |
| 5.555555556e-01 | 5.897214460e-05 | 6.203257430e-05 | -3.060429700e-06 | 4.554964777e-02 | 4.258315818e-03 | -7.494661223e-04 |
| 7.777777778e-01 | 6.231273838e-05 | 5.787070107e-05 | 4.442037308e-06 | 7.714853708e-02 | 1.753405028e-02 | -7.490966005e-04 |
| 1.000000000e+00 | 6.067365409e-05 | 5.681291270e-05 | 3.860741396e-06 | 1.622884742e-01 | 2.477685575e-02 | -7.459233685e-04 |

## Plots

- `slice_mean_x.png`
- `slice_rms_x.png`
- `slice_rms_px.png`
- `slice_rms_z.png`
- `slice_max_abs_px_diff.png`
