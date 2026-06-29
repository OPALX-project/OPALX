# Exit Monitor Longitudinal-Slice DT Summary

This uses existing finite-distribution OPALX ELEMEDGE and OPAL 2022.1
exit-monitor H5 files. Incomplete timestep directories are skipped.

| DT [s] | n | OPALX step | OPAL step | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x rel diff [%] | centroid sigma rel diff [%] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.000000000e-12 | 2500 | Step#0 | Step#0 | 9.251667575285e+00 | 9.251933018620e+00 | 3.642769548e+01 | 0.916993 | 1.805761e+00 | 1.241079e+03 |
| 5.000000000e-12 | 2500 | Step#0 | Step#0 | 9.251669234959e+00 | 9.255537067310e+00 | 7.955888371e+01 | 0.822915 | 1.765055e+00 | 6.018098e+00 |

## H5 Steps

### DT 2.000000000e-12 s

- OPALX: `Step#0:n=2500:s=9.251667575285e+00`
- OPAL: `Step#0:n=2500:s=9.251933018620e+00`

### DT 5.000000000e-12 s

- OPALX: `Step#0:n=2500:s=9.251669234959e+00`
- OPAL: `Step#0:n=2500:s=9.255537067310e+00`

## Plots

- `dt_monitor_slice_summary.png`
- `dt_monitor_spos_difference.png`
