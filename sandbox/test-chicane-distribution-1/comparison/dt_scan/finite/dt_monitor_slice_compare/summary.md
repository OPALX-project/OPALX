# Exit Monitor Longitudinal-Slice DT Summary

This uses existing finite-distribution OPALX ELEMEDGE and OPAL 2022.1
exit-monitor H5 files. Incomplete timestep directories are skipped.

Skipped:

- dt_1p0em12: missing None

| DT [s] | n | OPALX step | OPAL step | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x rel diff [%] | centroid sigma rel diff [%] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.000000000e-12 | 125 | Step#1 | Step#0 | 9.250299677828e+00 | 9.251933018620e+00 | 8.495468905e+01 | 0.977239 | 1.318362e+00 | 1.238300e+03 |
| 5.000000000e-12 | 125 | Step#0 | Step#0 | 9.250753364361e+00 | 9.255537067310e+00 | 1.898318964e+02 | 0.994398 | 5.917372e-01 | 2.807532e+01 |
| 1.000000000e-11 | 125 | Step#1 | Step#0 | 9.251519195800e+00 | 9.261572943508e+00 | 2.062170763e+02 | 0.879509 | 2.349991e+01 | 2.407417e+01 |

## H5 Steps

### DT 2.000000000e-12 s

- OPALX: `Step#0:n=125:s=9.249999885412e+00; Step#1:n=125:s=9.250299677828e+00`
- OPAL: `Step#0:n=125:s=9.251933018620e+00`

### DT 5.000000000e-12 s

- OPALX: `Step#0:n=125:s=9.250753364361e+00; Step#1:n=125:s=9.250003883314e+00`
- OPAL: `Step#0:n=125:s=9.255537067310e+00`

### DT 1.000000000e-11 s

- OPALX: `Step#0:n=125:s=9.250020234125e+00; Step#1:n=125:s=9.251519195800e+00`
- OPAL: `Step#0:n=125:s=9.261572943508e+00`

## Plots

- `dt_monitor_slice_summary.png`
- `dt_monitor_spos_difference.png`
