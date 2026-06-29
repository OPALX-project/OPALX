# Element marker comparison

OPALX markers: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/data/test-chicane-distribution-2_ElementPositions.txt`
OPAL markers: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/data/test-chicane-distribution-1_opal_ElementPositions.txt`

Differences are `OPALX - OPAL`.  Coordinates are the three columns written by the generated element-position files.

| marker                                  | opalx_c0_m      | opalx_c1_m       | opal_c0_m       | opal_c1_m        | diff_c0_mm       | diff_c1_mm       | diff_norm_mm    |
| --------------------------------------- | --------------- | ---------------- | --------------- | ---------------- | ---------------- | ---------------- | --------------- |
| OPALX FIELD BEGIN: B1 vs OPAL BEGIN: B1 | 2.450000000e-01 | 1.734723476e-18  | 2.450000000e-01 | 0.000000000e+00  | 0.000000000e+00  | 1.734723476e-15  | 1.734723476e-15 |
| OPALX FIELD END: B1 vs OPAL END: B1     | 1.253725281e+00 | 5.047833635e-02  | 1.255101141e+00 | 5.061357466e-02  | -1.375860000e+00 | -1.352383100e-01 | 1.382490557e+00 |
| OPALX FIELD BEGIN: B2 vs OPAL BEGIN: B2 | 2.735866781e+00 | 1.991885178e-01  | 2.737242363e+00 | 1.993265363e-01  | -1.375582000e+00 | -1.380185000e-01 | 1.382488678e+00 |
| OPALX FIELD END: B2 vs OPAL END: B2     | 3.744592063e+00 | 2.496668541e-01  | 3.747350132e+00 | 2.498076668e-01  | -2.758069000e+00 | -1.408127000e-01 | 2.761661244e+00 |
| OPALX FIELD BEGIN: B3 vs OPAL BEGIN: B3 | 5.734175274e+00 | 2.496668541e-01  | 5.736933345e+00 | 2.498048727e-01  | -2.758071000e+00 | -1.380186000e-01 | 2.761522184e+00 |
| OPALX FIELD END: B3 vs OPAL END: B3     | 6.742900556e+00 | 1.991885178e-01  | 6.747034486e+00 | 1.991912980e-01  | -4.133930000e+00 | -2.780200000e-03 | 4.133930935e+00 |
| OPALX FIELD BEGIN: B4 vs OPAL BEGIN: B4 | 8.225042056e+00 | 5.047833635e-02  | 8.229175708e+00 | 5.047833635e-02  | -4.133652000e+00 | 0.000000000e+00  | 4.133652000e+00 |
| OPALX FIELD END: B4 vs OPAL END: B4     | 9.233767337e+00 | -7.042977312e-16 | 9.239283477e+00 | -2.794139361e-06 | -5.516140000e+00 | 2.794139360e-03  | 5.516140708e+00 |
