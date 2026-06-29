# Section Temporal Monitor Longitudinal Slice Comparison

Each run uses the original upstream lattice truncated at one section and adds
a single terminal temporal monitor at that section.  The reported coordinates
therefore use the same monitor-local convention as the final-exit diagnostic.

| section | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x OPALX [m] | sigma_x OPAL [m] | centroid sigma OPALX [m] | centroid sigma OPAL [m] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| after_b1 | 1.250416594e+00 | 1.251382143e+00 | 9.626202694e+01 | 0.848348 | 9.051549502e-05 | 9.061533863e-05 | 9.427077541e-06 | 9.512189147e-06 |
| after_b2 | 3.750835182e+00 | 3.752768812e+00 | -1.226961835e+02 | 0.816073 | 3.237417463e-04 | 3.243036311e-04 | 2.187291203e-05 | 2.733236025e-05 |
| after_b3 | 6.751254186e+00 | 6.754156797e+00 | -5.074452371e+01 | 0.933698 | 2.511173597e-04 | 2.535000978e-04 | 3.996537796e-05 | 4.462316082e-05 |
| after_b4 | 9.251669235e+00 | 9.255537067e+00 | 7.952424722e+01 | 0.824036 | 7.525249562e-05 | 7.392105877e-05 | 4.535636421e-05 | 4.273058671e-05 |

## Plots

- `monitor_centroid_sigma_x.png`
- `monitor_dpx_slope.png`
- `slice_mean_x_by_section.png`
- `slice_diff_mean_x_by_section.png`
- `slice_mean_px_by_section.png`
- `slice_diff_mean_px_by_section.png`
- `slice_rms_x_by_section.png`
- `slice_diff_rms_x_by_section.png`
