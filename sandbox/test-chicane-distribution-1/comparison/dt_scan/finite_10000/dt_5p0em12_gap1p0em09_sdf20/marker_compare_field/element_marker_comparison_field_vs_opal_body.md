# Element marker comparison

OPALX markers: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12_gap1p0em09_sdf20/opalx_elemedge/data/test-chicane-distribution-2_ElementPositions.txt`
OPAL markers: `sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12_gap1p0em09_sdf20/opal_2022/data/test-chicane-distribution-1_opal_ElementPositions.txt`

Differences are `OPALX - OPAL`.  Coordinates are the three columns written by the generated element-position files.

| marker                                  | opalx_c0_m      | opalx_c1_m       | opal_c0_m       | opal_c1_m        | diff_c0_mm       | diff_c1_mm       | diff_norm_mm    |
| --------------------------------------- | --------------- | ---------------- | --------------- | ---------------- | ---------------- | ---------------- | --------------- |
| OPALX FIELD BEGIN: B1 vs OPAL BEGIN: B1 | 2.499999950e-01 | 0.000000000e+00  | 2.499999950e-01 | 0.000000000e+00  | 0.000000000e+00  | 0.000000000e+00  | 0.000000000e+00 |
| OPALX FIELD END: B1 vs OPAL END: B1     | 1.250416793e+00 | 0.000000000e+00  | 1.249635706e+00 | 5.006793758e-02  | 7.810870000e-01  | -5.006793758e+01 | 5.007402990e+01 |
| OPALX FIELD BEGIN: B2 vs OPAL BEGIN: B2 | 2.740841797e+00 | 1.996876843e-01  | 2.741727231e+00 | 1.997765240e-01  | -8.854340000e-01 | -8.883970000e-02 | 8.898796889e-01 |
| OPALX FIELD END: B2 vs OPAL END: B2     | 3.736260679e+00 | 2.995627114e-01  | 3.741367380e+00 | 2.497557660e-01  | -5.106701000e+00 | 4.980694540e+01  | 5.006805574e+01 |
| OPALX FIELD BEGIN: B3 vs OPAL BEGIN: B3 | 5.739175269e+00 | 2.496668541e-01  | 5.740950582e+00 | 2.497556938e-01  | -1.775313000e+00 | -8.883970000e-02 | 1.777534455e+00 |
| OPALX FIELD END: B3 vs OPAL END: B3     | 6.739592068e+00 | 2.496668541e-01  | 6.740586293e+00 | 1.996877562e-01  | -9.942250000e-01 | 4.997909790e+01  | 4.998898589e+01 |
| OPALX FIELD BEGIN: B4 vs OPAL BEGIN: B4 | 8.230017072e+00 | 4.997916977e-02  | 8.232677818e+00 | 4.997916977e-02  | -2.660746000e+00 | 0.000000000e+00  | 2.660746000e+00 |
| OPALX FIELD END: B4 vs OPAL END: B4     | 9.225435953e+00 | -4.989585727e-02 | 9.232317967e+00 | -7.221976868e-08 | -6.882014000e+00 | -4.989578505e+01 | 5.036815941e+01 |
