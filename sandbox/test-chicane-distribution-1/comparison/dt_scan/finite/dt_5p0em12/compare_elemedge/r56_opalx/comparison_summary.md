# test-chicane-distribution-1 comparison

Selected monitor step: `Step#0`
Monitor `SPOS`: `9.250753364361e+00` m
Particles at exit: `125`

## Manufactured longitudinal check

- Reference `R56`: `-4.333333333333e-02` m
- Fitted transport `R56`: `-4.445817259791e-02` m
- Difference: `-1.124839264580e-03` m
- Fitted monitor offset: `-7.449643967315e-04` m
- Raw `z` residual RMS: `7.449701888950e-04` m
- Raw `z` residual max abs: `7.517791304347e-04` m
- Intercept-corrected `z` residual RMS: `2.937676766969e-06` m
- Intercept-corrected `z` residual max abs: `6.814733703166e-06` m

## Transverse closure diagnostics

- `y` residual RMS: `1.724756495082e-08` m
- `y` residual max abs: `5.251586918013e-08` m
- final `x` mean: `-2.753103777550e-05` m
- final `x` RMS: `6.904089791971e-05` m

## Checks

- `all_particles_at_exit`: `True`
- `r56_within_tolerance`: `True`
- `z_rms_within_tolerance`: `True`
- `z_max_within_tolerance`: `True`
- `y_rms_within_tolerance`: `True`

Overall: `PASS`
