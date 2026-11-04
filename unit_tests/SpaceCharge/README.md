# FFT2D5 unit checks

`TestFFT2D5Algorithm` checks that a fixed laboratory source gives the same
physical fields when the moving reference frame changes. It also covers
multiple active containers, repeated solves, inactive containers and restoration
of coordinates and time steps after an open-path coverage error.

The complete solve prepares path geometry and particle projections once. A
focused comparison against the uncached diagnostic stages checks deposited
density and all electric/magnetic field components on nonuniform paths, including
zero-length segments, corner ties, periodic-seam probes, invalid particles and
changes of particle count, ordering and momentum. Standalone diagnostic stages
recompute projections; cached entries are refreshed for every complete solve.

`TestFFT2D5Numerics` checks the transverse slice solver against independent
rectangular-cell logarithmic-potential quadrature and analytic Gaussian fields. It covers
anisotropic cells, changes of length unit, boundary-cell self-force, potential
halos, mesh convergence and reuse of the slice workspace. Empty slices must
clear their fields; neutral slices with nonzero charge density must still be
solved.

The potential reference is the logarithmic Green function averaged over a
rectangular source cell, including its finite self term. The force is the
centered difference of that potential. Older snapshots of a point-sampled
kernel with a zero self term are replaced by this independent reference;
passing the tests does not establish convergence of a tracking simulation.

`TestSpaceChargeConfig` checks the supported solver settings, including the
explicit `FFT2D5MODEL="ORIGINAL"` selection. The original slice model does not
include a finite-source 3D operator or a curvature-correction operator.

Build and run the normal targets, for example:

```sh
cmake --build build-omp --target TestFFT2D5Algorithm TestFFT2D5Numerics TestSpaceChargeConfig --parallel 4
OMP_NUM_THREADS=1 OMP_PROC_BIND=false ctest --test-dir build-omp -R '^(TestSpaceChargeConfig|TestFFT2D5Algorithm|TestFFT2D5Numerics)$' --output-on-failure
OMP_NUM_THREADS=4 OMP_PROC_BIND=false ctest --test-dir build-omp -R '^(TestSpaceChargeConfig|TestFFT2D5Algorithm|TestFFT2D5Numerics)$' --output-on-failure
```

FFT2D5 currently supports one MPI rank. Its rejection of multiple ranks can be
checked separately:

```sh
OMP_NUM_THREADS=1 OMP_PROC_BIND=false mpiexec -n 2 build-omp/unit_tests/SpaceCharge/TestFFT2D5Numerics --gtest_filter=TestSolve2d5.Configuration_MultipleRanksRejected
```

Tracking inputs, published plots and the analytic ZeroSC comparison are stored
in the separate `OPALX-project/opalx-documents` repository under
`examples/2026/2d5-1`.
