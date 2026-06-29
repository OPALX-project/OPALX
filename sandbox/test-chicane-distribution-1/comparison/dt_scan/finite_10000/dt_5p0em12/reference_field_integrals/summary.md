# Reference Field Integral Comparison

Field column: `By_ref`
Threshold for interval detection: `1e-06 T`
Samples per bend: `4001`

| bend | s_start_m       | s_end_m         | opalx_integral_Tm | opal_integral_Tm | integral_diff_Tm | integral_rel_diff | max_abs_cumulative_diff_Tm | max_abs_field_diff_T |
| ---- | --------------- | --------------- | ----------------- | ---------------- | ---------------- | ----------------- | -------------------------- | -------------------- |
| B1   | 2.548235561e-01 | 1.244138536e+00 | 3.300309896e-01   | 3.300309896e-01  | 3.315125952e-13  | 1.004489292e-12   | 3.315681063e-13            | 4.352046501e-11      |
| B2   | 2.758090251e+00 | 3.747405231e+00 | -3.300309896e-01  | -3.300309896e-01 | 9.388601008e-13  | -2.844763463e-12  | 9.388601008e-13            | 2.289589074e-10      |
| B3   | 5.756014437e+00 | 6.745329417e+00 | -3.300309896e-01  | -3.300309898e-01 | 2.444530134e-10  | -7.406971494e-10  | 2.444530134e-10            | 3.239107893e-08      |
| B4   | 8.259281132e+00 | 9.248596112e+00 | 3.300309896e-01   | 3.300309896e-01  | -2.164934898e-15 | -6.559792766e-15  | 4.576894419e-14            | 3.873443233e-10      |

## Total Over Four Bend Intervals

- OPALX: `7.555067682574e-13 T m`
- OPAL: `-2.449657143799e-10 T m`
- Difference: `2.457212211482e-10 T m`
- Relative difference: `-1.003084132692e+00`

This diagnostic compares stat-file reference fields only. It does not yet
force both codes to evaluate fields on an externally prescribed 3D trajectory.
