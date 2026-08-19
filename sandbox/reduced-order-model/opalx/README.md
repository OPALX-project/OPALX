# OPALX reduced-order inputs

`rigid_two_gaussian_field_snapshot.in.template` is the tracked source for the
four matching OPALX field-sampling cases. Generate them from the repository root:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/make_opalx_field_cases.py
```

Add `--run` to execute the generated cases with `build_openmp/src/opalx`, then
compare the first two-source (`Step#1`) witness fields with:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/compare_opalx_fields.py --cases 3sigma
```

The baseline generator defaults are 400000 source macroparticles and a
64x64x128 mesh. Its diagnostic/plotting grid remains 101x201. The persistent
regression overrides only the passive-probe grid to 51x101; both dimensions are
odd so that the IP is sampled exactly. For a fast parser/runtime smoke test,
use a separate output directory and reduced numerical parameters, for example
`--primary-macroparticles 20000 --nxy 16 --nz 32 --probe-nx 21 --probe-nz 41`.

The template deliberately selects `GREENSF=INTEGRATED`. The relativistic binned
solver transforms the longitudinal mesh spacing to the primary-beam rest frame
by multiplying it by the bin Lorentz factor. For this 245 MeV case the resulting
cell aspect ratio is about 190000:1. `GREENSF=STANDARD` treats each deposited
cell as a point source and overestimates the sampled field by several orders of
magnitude on that mesh; the cell-integrated kernel is required for the
reduced-order comparison.

Each case uses a passive electron probe container on the same x-z grid as the
Python model. `EBDUMP=TRUE` writes all E and B components to its `_c1.h5` file.
The probes do not deposit charge. `BBRIGID=TRUE` removes the source response but
does not remove the field gathered to these probes.

The otherwise identical self-consistent control deck should set
`BBRIGID=FALSE` (or omit the option).

`Step#0` contains only the physical primary. Reproduce that single-source check
with:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/compare_opalx_fields.py \
  --cases 3sigma --step 0 --source-model physical
```

After `COPY_TIME`, `Step#1` contains the physical and mirrored
counter-propagating primaries. Reproduce the two-source check with:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/compare_opalx_fields.py \
  --cases 3sigma --step 1 --source-model physical-and-copied
```

Keep generated `.h5`, `.stat`, timing, and visualization output out of the input
set. Each validation case should record its OPALX revision, MPI rank count, mesh,
time step, source macro-particle count, and random seed.

The persistent, plot-free 3-sigma comparison is registered as an opt-in
integration test. Configure it with the OPALX Python environment and run it via
CTest:

```bash
cmake -S . -B build_openmp \
  -DOPALX_ENABLE_TESTS=ON \
  -DPython3_EXECUTABLE="$HOME/.venv-h6/bin/python"
ctest --test-dir build_openmp \
  -R '^TestBeamBeamManufacturedFields$' --output-on-failure
```

The test uses `regression.json` for its 100000-particle source, 51x101 probe grid,
and acceptance thresholds. The particle count intentionally covers the
short-tracking-horizon case that previously clipped the OrbitThreader range and
ended the BeamBeam window after the first solve on CUDA. It writes its generated
inputs, small H5/STAT outputs, and compact summary CSV below
`build_openmp/test/BeamBeam/`; it never creates plots.

On Merlin, run a multi-rank GPU case in its own batch allocation. For example,
request one node, two tasks, and two GPUs in the job header, then use the
cluster's Slurm-aware Open MPI launcher:

```bash
mpiexec -n 2 /path/to/opalx rigid_fields_3sigma.in \
  --kokkos-map-device-id-by=mpi_rank
```

Do not append this launch to the four-GPU batch that runs many exclusive
single-rank `srun` steps. In job 353949 both OPALX ranks completed, but the
parent PRRTE launcher remained alive after rank exit in that mixed workflow.
Also do not use plain multi-task `srun` on this cluster: its PMI2 plugin is not
connected to this Open MPI build, so it starts independent MPI singletons. The
dedicated two-task `mpiexec` job returned normally for the full 400000-particle
case and reported `MPI_Comm_size=2`.

The current A100 coverage should not be overstated. Job 353948 used four
connected MPI ranks on four A100s, but only for the distributed 8x8x8
mirror-field unit test; it did not execute the OPALX tracking loop. The full
manufactured deck has been validated with at most two ranks on two A100s. It
uses `MAXSTEPS=2`: Step#0 contains the physical-source solve and Step#1 the
physical-plus-copied-source solve. A four-rank/four-A100 full-deck run remains
future validation work.

For convergence studies, keep each completed case in a directory named
`N<particles-in-thousands>k_M<transverse-mesh>_S<seed>/3sigma` and add an empty
`completed` marker only after OPALX exits successfully. Aggregate the Step#1
physical-plus-copied comparisons with:

```bash
MPLBACKEND=Agg ~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/analyze_opalx_convergence.py \
  --case-root /path/to/convergence-cases \
  --output-dir /path/to/convergence-summary
```

The analyzer reports fixed-mesh particle convergence, fixed-particle mesh
behavior, constant-particles-per-cell mesh convergence, and seed spread. It
keeps cases with missing probe coverage in `invalid_cases.csv` but excludes them
from convergence plots. Generated H5, CSV, logs, and PNG files remain build/run
artifacts and must not be committed.
