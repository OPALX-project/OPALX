# test-chicane-distribution-1 comparison

Selected monitor step: `Step#0`
Monitor `SPOS`: `9.250753364361e+00` m
Particles at exit: `10000`

## Manufactured longitudinal check

- Reference `R56`: `-4.333333333333e-02` m
- Fitted transport `R56`: `-4.451959196508e-02` m
- Difference: `-1.186258631743e-03` m
- Fitted monitor offset: `-7.450292687103e-04` m
- Raw `z` residual RMS: `7.450336473000e-04` m
- Raw `z` residual max abs: `7.515662654433e-04` m
- Intercept-corrected `z` residual RMS: `2.554285459839e-06` m
- Intercept-corrected `z` residual max abs: `6.536996733051e-06` m

## Transverse closure diagnostics

- `y` residual RMS: `1.455387964380e-08` m
- `y` residual max abs: `5.613584278150e-08` m
- final `x` mean: `-3.732412769944e-05` m
- final `x` RMS: `8.954720529147e-05` m

## Checks

- `all_particles_at_exit`: `True`
- `r56_within_tolerance`: `False`
- `z_rms_within_tolerance`: `True`
- `z_max_within_tolerance`: `True`
- `y_rms_within_tolerance`: `False`

Overall: `FAIL`
