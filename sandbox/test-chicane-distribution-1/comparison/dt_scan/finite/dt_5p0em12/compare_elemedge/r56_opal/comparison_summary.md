# test-chicane-distribution-1 comparison

Selected monitor step: `Step#0`
Monitor `SPOS`: `9.255537067310e+00` m
Particles at exit: `125`

## Manufactured longitudinal check

- Reference `R56`: `-4.333333333333e-02` m
- Fitted transport `R56`: `-4.444156476076e-02` m
- Difference: `-1.108231427431e-03` m
- Fitted monitor offset: `3.045530522560e-06` m
- Raw `z` residual RMS: `4.183581535324e-06` m
- Raw `z` residual max abs: `9.772692469332e-06` m
- Intercept-corrected `z` residual RMS: `2.868291878254e-06` m
- Intercept-corrected `z` residual max abs: `6.727161946773e-06` m

## Transverse closure diagnostics

- `y` residual RMS: `2.012220802610e-08` m
- `y` residual max abs: `6.682376340005e-08` m
- final `x` mean: `-1.774147563692e-05` m
- final `x` RMS: `6.863476051286e-05` m

## Checks

- `all_particles_at_exit`: `True`
- `r56_within_tolerance`: `True`
- `z_rms_within_tolerance`: `True`
- `z_max_within_tolerance`: `True`
- `y_rms_within_tolerance`: `True`

Overall: `PASS`
