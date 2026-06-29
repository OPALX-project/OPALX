# Element marker comparison

OPALX markers: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/data/test-chicane-distribution-2_ElementPositions.txt`
OPAL markers: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/data/test-chicane-distribution-1_opal_ElementPositions.txt`

Differences are `OPALX - OPAL`.  Coordinates are the three columns written by the generated element-position files.

| marker                                  | opalx_c0_m      | opalx_c1_m       | opal_c0_m       | opal_c1_m        | diff_c0_mm       | diff_c1_mm       | diff_norm_mm    |
| --------------------------------------- | --------------- | ---------------- | --------------- | ---------------- | ---------------- | ---------------- | --------------- |
| OPALX FIELD BEGIN: B1 vs OPAL BEGIN: B1 | 2.450000000e-01 | 0.000000000e+00  | 2.450000000e-01 | 0.000000000e+00  | 0.000000000e+00  | 0.000000000e+00  | 0.000000000e+00 |
| OPALX FIELD END: B1 vs OPAL END: B1     | 1.255416788e+00 | 0.000000000e+00  | 1.255101141e+00 | 5.061357466e-02  | 3.156470000e-01  | -5.061357466e+01 | 5.061455890e+01 |
| OPALX FIELD BEGIN: B2 vs OPAL BEGIN: B2 | 2.735866781e+00 | 1.991885178e-01  | 2.737242363e+00 | 1.993265363e-01  | -1.375582000e+00 | -1.380185000e-01 | 1.382488678e+00 |
| OPALX FIELD END: B2 vs OPAL END: B2     | 3.741235694e+00 | 3.000618780e-01  | 3.747350132e+00 | 2.498076668e-01  | -6.114438000e+00 | 5.025421120e+01  | 5.062481699e+01 |
| OPALX FIELD BEGIN: B3 vs OPAL BEGIN: B3 | 5.734175274e+00 | 2.496668541e-01  | 5.736933345e+00 | 2.498048727e-01  | -2.758071000e+00 | -1.380186000e-01 | 2.761522184e+00 |
| OPALX FIELD END: B3 vs OPAL END: B3     | 6.744592063e+00 | 2.496668541e-01  | 6.747034486e+00 | 1.991912980e-01  | -2.442423000e+00 | 5.047555610e+01  | 5.053461382e+01 |
| OPALX FIELD BEGIN: B4 vs OPAL BEGIN: B4 | 8.225042056e+00 | 5.047833635e-02  | 8.229175708e+00 | 5.047833635e-02  | -4.133652000e+00 | 0.000000000e+00  | 4.133652000e+00 |
| OPALX FIELD END: B4 vs OPAL END: B4     | 9.230410969e+00 | -5.039502385e-02 | 9.239283477e+00 | -2.794139361e-06 | -8.872508000e+00 | -5.039222971e+01 | 5.116735496e+01 |
