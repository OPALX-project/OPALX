# Section Stop Longitudinal Slice Comparison

The original lattice is unchanged. Each run only changes `ZSTOP` and compares
the normal H5 bunch output at the truncated stop.

| stop | OPALX SPOS [m] | OPAL SPOS [m] | dpx slope [1/m] | R^2 | sigma_x OPALX [m] | sigma_x OPAL [m] | centroid sigma OPALX [m] | centroid sigma OPAL [m] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| after_b1 | 1.299600133e+00 | 1.300000000e+00 | 2.235180504e-01 | 0.000346 | 9.503931872e-05 | 9.507439573e-05 | 9.689259527e-06 | 9.609129347e-06 |
| after_b2 | 3.799868904e+00 | 3.800000000e+00 | 4.243602259e+00 | 0.012030 | 3.235657696e-04 | 3.241703734e-04 | 2.204552150e-05 | 2.757601163e-05 |
| after_b3 | 6.799292052e+00 | 6.800000000e+00 | 8.970090267e-01 | 0.020921 | 2.449422505e-04 | 2.476194053e-04 | 4.002003746e-05 | 4.472901757e-05 |
| after_b4 | 9.299560823e+00 | 9.300000000e+00 | 2.626572051e+00 | 0.017000 | 7.555034088e-05 | 7.397850546e-05 | 4.577934313e-05 | 4.274002274e-05 |

## Plots

- `section_centroid_sigma_x.png`
- `section_dpx_slope.png`
