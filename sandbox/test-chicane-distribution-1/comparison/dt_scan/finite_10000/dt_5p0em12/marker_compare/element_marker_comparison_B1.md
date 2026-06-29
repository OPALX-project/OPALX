# Element marker comparison

OPALX markers: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/data/test-chicane-distribution-2_ElementPositions.txt`
OPAL markers: `/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/data/test-chicane-distribution-1_opal_ElementPositions.txt`

Differences are `OPALX - OPAL`.  Coordinates are the three columns written by the generated element-position files.

| marker         | opalx_c0_m      | opalx_c1_m      | opal_c0_m       | opal_c1_m       | diff_c0_mm       | diff_c1_mm       | diff_norm_mm    |
| -------------- | --------------- | --------------- | --------------- | --------------- | ---------------- | ---------------- | --------------- |
| ENTRY EDGE: B1 | 2.500000000e-01 | 0.000000000e+00 | 2.500000000e-01 | 0.000000000e+00 | 0.000000000e+00  | 0.000000000e+00  | 0.000000000e+00 |
| BEGIN: B1      | 2.500000000e-01 | 0.000000000e+00 | 2.450000000e-01 | 0.000000000e+00 | 5.000000000e+00  | 0.000000000e+00  | 5.000000000e+00 |
| END: B1        | 1.248750260e+00 | 4.997916927e-02 | 1.255101141e+00 | 5.061357466e-02 | -6.350881000e+00 | -6.344053900e-01 | 6.382488517e+00 |
| EXIT EDGE: B1  | 1.248750260e+00 | 4.997916927e-02 | 1.248750260e+00 | 4.997916927e-02 | 0.000000000e+00  | 0.000000000e+00  | 6.120682967e-15 |
