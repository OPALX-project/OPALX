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

The nearest-segment hierarchy is checked directly against an independent exhaustive
scan: segment indices and all three cached coordinates must match exactly, with
no numerical tolerance. Probes include self-intersections, repeated segments and
vertices, corner/seam ties, adjacent floating-point coordinates, stale seeds after
movement/reordering, allocation changes, large translations and very small/large
coordinate scales (including subnormal squared lengths). Bounding-box pruning is
conservative; neither the previous segment nor path-index proximity restricts the
search. The existing field/density checks also compare the complete cached solve
with the uncached conversion.

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

FFT2D5 supports MPI decomposition along the reference-path coordinate with
`PARFFTX=FALSE, PARFFTY=FALSE, PARFFTZ=TRUE`. Each rank owns complete transverse
planes and uses private heFFTe plans. Particle migration and field halo exchange
reuse the same IPPL infrastructure as OPEN; IPPL is unchanged. The small global
line-density profile is reduced separately for the longitudinal derivative.

`TestFFT2D5Parallel` compares the distributed fields with direct discrete Coulomb
convolution using the cell-integrated logarithmic Green function, independently
of FFT execution. It checks open/closed paths, the periodic seam, slab boundaries,
both longitudinal deposition shapes, empty ranks and repeated migration. It also
checks all particle electric/magnetic components and restoration of coordinates,
momenta and IDs. Run it on 1, 2 and 4 ranks:

```sh
OMP_NUM_THREADS=1 mpiexec -n 1 build-omp/unit_tests/SpaceCharge/TestFFT2D5Parallel
OMP_NUM_THREADS=1 mpiexec -n 2 build-omp/unit_tests/SpaceCharge/TestFFT2D5Parallel
OMP_NUM_THREADS=1 mpiexec -n 4 build-omp/unit_tests/SpaceCharge/TestFFT2D5Parallel
```

Tracking inputs, published plots and the analytic ZeroSC comparison are stored
in the separate `OPALX-project/opalx-documents` repository under
`examples/2026/2d5-1`.

### Optional OPALX-side slice batching

`opalx case.in --2D5BatchSize 8` groups eight independent transverse slice
transforms using existing heFFTe R2C and complex row executors. One R2C
operation handles every x row across the batch; a transpose makes the y columns
contiguous for a second batched transform. The inverse reverses these steps.
Only x and y are transformed. Plans are cached for the full batch and current
remainder; source normalization is folded into the Green multiplication.
The option requires a positive integer; without it the serial slice path remains active. It is a
CPU and CUDA optimization and does not transform the longitudinal axis.
With z decomposition it batches each rank's local planes. If no batch size is
specified, multiple ranks use a local batch size of one. No IPPL source changes
are required. Inspection policies retain per-slice potential callbacks; with MPI,
callbacks receive local slice views and global longitudinal indices.

The batch tests compare potential including its transverse halo, all field
components, unchanged input charge and cleared empty/vacated slices with the
serial solver for sizes 1, 8, 32 and 64. They cover neutral sources, incomplete
batches, anisotropic rectangular grids, grid/spacing changes and the complete
solve under both open and closed longitudinal boundary modes. CLI parsing,
normalization and the same 100-step tracking pilot are checked separately.

### FFT2D5 timing and warmup

Like OPEN, FFT2D5 prepares its backend on a zero source before physical solves.
`2D5 Initialize` measures storage/path initialization; `2D5 Warmup` measures
backend preparation, including the separately reported `2D5 Plan` and
`2D5 Green Init` subphases. Warmup does not move particles, advance tracking,
or increment physical backend-solve counts. FFT execution during warmup is
excluded from runtime phase counters without resetting other instances' samples.

`2D5 First Solve` measures the first physical call, including first-use particle
and projection storage. `2D5 Solve` measures subsequent complete space-charge
calls, from frame entry through frame restoration and moments. A 100-step run
therefore has one first-solve sample and 99 steady-state samples. No physical
step is discarded. The old `mainTimer - OrbThreader` difference also includes
other tracking and output work and is not a solver-only measurement.

Phase timers cover geometry/projections, migration, deposition, density/field
halos, line density, transverse Poisson work, packing, forward/inverse FFTs,
Green multiplication, field extraction and gathering. They use the existing
IPPL timing output (`timing.dat` and the log), with local Kokkos fences so CUDA
measurements include completed device work. No extra MPI barrier is added.
Timers are inclusive: do not add parent and child totals. Phase totals include
all physical calls; geometry also includes post-migration projection refresh.
With `PIPEMODE=NONE`, the line-density timer measures only clearing disabled
longitudinal-field buffers. The batch path has no separate `2D5 Extract` work;
its field extraction is included in `2D5 Field`.

For comparison with OPEN, `2D5 Poisson` is the backend timer; `2D5 Solve` has the
broader complete-space-charge scope described above. The warmup/timer regression
also checks that warming a second solver does not erase earlier measurements.

### Optional heFFTe slice-batching microbenchmark

`BenchmarkFFT2D5Batch` is an explicit, non-CTest experiment using the existing
heFFTe public plan and executor factories. It compares transverse potential and
fields with the actual `FFT2D5Poisson::solveBatched()` before timing. It does not
modify the production solver or transform the longitudinal coordinate.

Build only this target in the existing backend build directory, then run it from
an empty output directory. `--validate-only` checks rectangular, partial, sparse,
neutral and empty batches; `--quick` adds one small timing sweep. The default
measures 64×64×256/512 with batch sizes 1, 8, 32 and 64 and writes
`validation.csv` and `samples.csv`. Setup and warmup are excluded. The pinned
heFFTe no-reorder option is tested only at batch size one because multi-batch
numerical equivalence fails. This executable requires one MPI rank; it measures
local Poisson work, not MPI exchange or full tracking.
