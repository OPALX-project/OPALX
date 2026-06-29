# test-chicane-distribution-1 comparison

Selected monitor step: `Step#0`
Monitor `SPOS`: `9.255537067310e+00` m
Particles at exit: `10000`

## Manufactured longitudinal check

- Reference `R56`: `-4.333333333333e-02` m
- Fitted transport `R56`: `-4.447862731188e-02` m
- Difference: `-1.145293978543e-03` m
- Fitted monitor offset: `3.175156557922e-06` m
- Raw `z` residual RMS: `4.324453248464e-06` m
- Raw `z` residual max abs: `9.907862746753e-06` m
- Intercept-corrected `z` residual RMS: `2.935860475369e-06` m
- Intercept-corrected `z` residual max abs: `6.827202302503e-06` m

## Transverse closure diagnostics

- `y` residual RMS: `1.653817787715e-08` m
- `y` residual max abs: `6.787245499365e-08` m
- final `x` mean: `-2.529901666722e-05` m
- final `x` RMS: `7.392105876699e-05` m

## Checks

- `all_particles_at_exit`: `True`
- `r56_within_tolerance`: `False`
- `z_rms_within_tolerance`: `True`
- `z_max_within_tolerance`: `True`
- `y_rms_within_tolerance`: `False`

Overall: `FAIL`
