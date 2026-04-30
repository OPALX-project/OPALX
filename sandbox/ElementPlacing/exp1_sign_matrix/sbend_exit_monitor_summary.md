# SBEND Exit Monitor Comparison

This study compares the same 1000-particle Gaussian bunch at the entrance and at the SBEND exit plane for `ANGLE = ±pi/4, ±pi/2` and `E1 = 0, 10deg`.

Definitions:
- `rms_x`, `rms_y`: RMS transverse coordinates in the local reference frame at the comparison plane.
- `rms_zeta`: RMS longitudinal phase coordinate in metres, reconstructed as `zeta = -beta0 * c * (t - t_ref)`.
- `rms_px_can`, `rms_py_can`: RMS canonical transverse momenta `P_x / P_0` and `P_y / P_0`.
- `rms_delta`: RMS relative momentum error `(P - P_0) / P_0`.

For BMAD the exit quantities come directly from the canonical beam output file. For OPALX the exit quantities are reconstructed from the monitor H5 dump, which stores local particle coordinates, local beta-gamma momenta, and particle crossing times at the monitor plane.

## RMS Comparison
| row      | code      | ANGLE | E1    | status | rms_x   | rms_y   | rms_zeta | rms_px_can | rms_py_can | rms_delta |
| -------- | --------- | ----- | ----- | ------ | ------- | ------- | -------- | ---------- | ---------- | --------- |
| ENTRANCE | COMMON    | ALL   | ALL   | OK     | 0.00099 | 0.00099 | 0.00101  | 0.00105    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | +pi/4 | 0     | OK     | 0.00126 | 0.00145 | 0.00130  | 0.00115    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | +pi/4 | 0     | OK     | 0.00141 | 0.00145 | 0.00000  | 0.00151    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | -pi/4 | 0     | OK     | 0.00120 | 0.00145 | 0.00131  | 0.00118    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | -pi/4 | 0     | OK     | 0.00133 | 0.00145 | 0.00000  | 0.00156    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | +pi/2 | 0     | OK     | 0.00093 | 0.00145 | 0.00163  | 0.00178    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | +pi/2 | 0     | OK     | 0.00128 | 0.00145 | 0.00000  | 0.00231    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | -pi/2 | 0     | OK     | 0.00091 | 0.00145 | 0.00162  | 0.00192    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | -pi/2 | 0     | OK     | 0.00123 | 0.00145 | 0.00000  | 0.00252    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | +pi/4 | 10deg | OK     | 0.00126 | 0.00145 | 0.00130  | 0.00115    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | +pi/4 | 10deg | OK     | 0.00145 | 0.00145 | 0.00000  | 0.00149    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | -pi/4 | 10deg | OK     | 0.00120 | 0.00145 | 0.00131  | 0.00118    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | -pi/4 | 10deg | OK     | 0.00129 | 0.00145 | 0.00000  | 0.00148    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | +pi/2 | 10deg | OK     | 0.00093 | 0.00145 | 0.00163  | 0.00178    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | +pi/2 | 10deg | OK     | 0.00129 | 0.00145 | 0.00000  | 0.00234    | 0.00102    | 0.00100   |
| EXIT     | BMAD      | -pi/2 | 10deg | OK     | 0.00091 | 0.00145 | 0.00162  | 0.00192    | 0.00102    | 0.00100   |
| EXIT     | OPALX-MON | -pi/2 | 10deg | OK     | 0.00124 | 0.00145 | 0.00000  | 0.00239    | 0.00102    | 0.00100   |

## Exit-Plane Diagnostics
| ANGLE | E1    | signed_distance_min | signed_distance_max | signed_distance_final | crosses_nominal_exit_plane | final_ref_x | final_ref_z | final_t_ns | final_s_path_m | status | monitor_s | t_ref_ns | ref_bg0    | n_particle |
| ----- | ----- | ------------------- | ------------------- | --------------------- | -------------------------- | ----------- | ----------- | ---------- | -------------- | ------ | --------- | -------- | ---------- | ---------- |
| +pi/4 | 0     | -0.89819            | 0.49895             | 0.49895               | 1.00000                    | -0.72606    | 1.25280     | 5.00000    | 1.49896        | OK     | 1.00001   | 3.33569  | 1957.97462 | 1000       |
| -pi/4 | 0     | -0.89819            | 0.49897             | 0.49897               | 1.00000                    | 0.72615     | 1.25275     | 5.00000    | 1.49896        | OK     | 0.99999   | 3.33560  | 1957.97462 | 1000       |
| +pi/2 | 0     | -0.63661            | 0.49894             | 0.49894               | 1.00000                    | -1.13556    | 0.63564     | 5.00000    | 1.49896        | OK     | 1.00002   | 3.33573  | 1957.97462 | 1000       |
| -pi/2 | 0     | -0.63661            | 0.49898             | 0.49898               | 1.00000                    | 1.13560     | 0.63556     | 5.00000    | 1.49896        | OK     | 0.99998   | 3.33557  | 1957.97462 | 1000       |
| +pi/4 | 10deg | -0.89819            | 0.49809             | 0.49809               | 1.00000                    | -0.74529    | 1.23235     | 5.00000    | 1.49896        | OK     | 1.00001   | 3.33569  | 1957.97462 | 1000       |
| -pi/4 | 10deg | -0.89819            | 0.49852             | 0.49852               | 1.00000                    | 0.71011     | 1.26814     | 5.00000    | 1.49896        | OK     | 1.00002   | 3.33570  | 1957.97462 | 1000       |
| +pi/2 | 10deg | -0.63661            | 0.49138             | 0.49138               | 1.00000                    | -1.12800    | 0.55366     | 5.00000    | 1.49896        | OK     | 1.00002   | 3.33573  | 1957.97462 | 1000       |
| -pi/2 | 10deg | -0.63661            | 0.49193             | 0.49193               | 1.00000                    | 1.12855     | 0.72302     | 5.00000    | 1.49896        | OK     | 1.00082   | 3.33839  | 1957.97462 | 1000       |
