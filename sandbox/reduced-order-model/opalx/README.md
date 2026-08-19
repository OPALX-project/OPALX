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

The `CUTOFFX`, `CUTOFFY`, and `CUTOFFLONG` values in the Gaussian form of the
template are currently accepted but ignored by OPALX `GAUSS`; its active
position sampler uses fixed three-sigma truncation and then recentres the
finite sample. The deterministic MPI source reproduces those semantics.

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

Job 353965 completed the full manufactured deck with four connected MPI ranks
on four A100s. It uses `MAXSTEPS=2`: Step#0 contains the physical-source solve
and Step#1 the physical-plus-copied-source solve. The run took 9.41 seconds and
gave relative L2(E)=1.3643%, relative L2(B)=1.2356%, and zero uncovered probes.
This required every MPI rank to participate in the distributed witness gather,
including ranks that locally own zero witness particles.

## Reproducible Merlin A100 convergence study

The authoritative matrix is `../a100_convergence_study.json`. From a Merlin6
login shell, one command prepares all inputs and submits a dependency chain:

```bash
cd /psi/home/adelmann/opalx-dev/opalx-a100-5a354101f
sandbox/reduced-order-model/merlin/submit_a100_convergence.sh
```

The command prints the absolute study directory and three Slurm job IDs. The
first allocation times the one-rank baseline in isolation, then runs at most
four other independent one-rank cases concurrently on four A100s. Dedicated
two-rank/two-A100 and four-rank/four-A100 allocations then run the baseline full
deck using the cluster's Slurm-aware `mpiexec`.
Every deck has exactly two tracking steps. A `study-completed` marker is written
only when all cases have returned successfully. The run directory can be fixed
explicitly with `--run-root`; `--opalx`, `--cluster`, `--partition`, `--account`,
and `--omp-threads` are also configurable.

To run only the MPI-decomposition check without repeating the convergence
matrix, add `--rank-only`. This requests one, two, and four A100s in three
dependent allocations and writes `rank-study-completed` after all three fixed-
source cases succeed:

```bash
sandbox/reduced-order-model/merlin/submit_a100_convergence.sh \
  --rank-only \
  --run-root /psi/home/adelmann/opalx-dev/a100-fixed-rank-<revision>
```

The matrix uses the fixed 51x101 passive-probe grid at 3 sigma separation:

| scan | source macroparticles | mesh `Nx x Ny x Nz` | seeds/ranks |
| --- | --- | --- | --- |
| particle sampling | 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4 million | 64x64x128 | seed 20260629, one rank |
| fixed-particle mesh | 0.4 million | 32x32x64 through 128x128x256 | seed 20260629, one rank |
| constant occupancy | 0.4, 3.2, 10.8, 25.6 million | 32x32x64 through 128x128x256 | seed 20260629, one rank |
| seed spread | 0.4 million | 64x64x128 | seeds 20260629--20260631, one rank |
| MPI decomposition | one fixed 0.4-million-particle `FROMFILE` sample | 64x64x128 | identical coordinates and momenta on 1/2/4 ranks and A100s |

Preparation writes the fixed primary once at the study root. Its positions are
metres; its momenta are absolute beta-gamma, so the corresponding primary
`BEAM` deliberately omits `PC`. Each rank reads the same file, and OPALX assigns
contiguous non-overlapping record ranges across ranks. `fixed_primary_metadata.json`
records the PCG64 seed, NumPy version, distribution parameters, centroid, and
file SHA-256. The MPI plots therefore measure deposition, solve, and gather
decomposition sensitivity without changing the primary Monte Carlo sample.

Fetch only the compact inputs, witness H5 fields, STAT files, logs, and timings,
then produce publication-quality PNG and PDF figures locally with:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/fetch_and_plot_a100_convergence.py \
  --remote-run-root /psi/home/adelmann/opalx-dev/a100-convergence-<revision>-<timestamp>
```

The local result is written below the ignored
`sandbox/reduced-order-model/outputs/a100-convergence/` directory. Replot an
already fetched study with the same command plus `--skip-fetch`. The analyzer
can also be invoked directly:

```bash
MPLBACKEND=Agg ~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/analyze_opalx_convergence.py \
  --case-root /path/to/convergence-cases \
  --output-dir /path/to/convergence-summary
```

Each completed case is named
`N<particles-in-thousands>k_M<transverse-mesh>_S<seed>[_MPI<ranks>]/3sigma`;
fixed-source cases use the explicit suffix `_FIXED_MPI<ranks>`.
`runtime_compute.txt` contains `/usr/bin/time -p` wall time, `opalx.log` contains
the complete application output, and `completed` is created only after a zero
exit status. `study_metadata.json` records the executable SHA-256, source
revision, source-patch SHA-256, and worktree status; `opalx_source.patch` stores
the actual tracked source delta. The analyzer reports the three convergence
sequences, seed spread, fixed-source MPI sensitivity, and strong scaling. Cases
with missing probe coverage remain in `invalid_cases.csv` but are excluded from
plots. Generated H5, CSV, logs, PNG, and PDF files remain run artifacts and must
not be committed.
