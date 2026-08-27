# Exact-timed OPALX Track-12 benchmark

This workspace compares six artificial CAIN electron--positron pairs with
OPALX and the anisotropic rigid two-Gaussian manufactured solution. The current
one-clock workflow is under [`timed/`](timed/); older pair-by-pair and
pre-reference-offset studies are retained under
[`../../attic/track12/`](../../attic/track12/).

See [`../../BEAMBEAM_PHYSICS_AND_VALIDATION.md`](../../BEAMBEAM_PHYSICS_AND_VALIDATION.md)
for the physics interpretation and current validation summary.

## Definition of the validation case

```text
BeamBeam length            8 mm
interaction point          element midpoint = 4 mm
longitudinal field window  20 mm, centered on the IP
transverse field domain    non-shrinking primary+witness envelope with aspect floor
physical aperture          2.4 mm x 0.24 mm rectangle
primary macroparticles     400000 deterministic particles
primary rms sizes          1.944325075701 um, 1.944325075701 um, 0.6 mm
primary support            component-wise 3 sigma truncation
witnesses                  6 electrons + 6 positrons
timestep                   1.8 um/c = 6.004153713566736 fs
production steps           1501
fine mesh                  4096 x 256 x 128
```

This is a validation geometry, not the 32 cm, 15 cm-radius production
BeamBeam element used for the 1,297 + 1,297 CAIN population.

The CAIN births occur at

```text
ct = -0.9000, -0.5994, -0.3006, 0.0000, 0.3006, 0.5994 mm.
```

The tracker begins one interval before the first birth. Every birth and output
sample lies on the CAIN time grid, so all 13,012 trajectory samples are matched
without interpolation. Electron and positron emitted files deliberately carry
the same numeric initial phase space; their separate OPALX beam containers
provide the species and charge sign.

The primary and its reflection about the element midpoint provide the two
source fields. `BBRIGID=TRUE` is used here to keep the manufactured comparison
controlled. Witness containers 1 and 2 gather the source field but never
deposit charge.

## Prepare and run locally

From this directory:

```bash
./timed/run_timed_track12.sh
```

The preparation step regenerates the deterministic primary, exact emitted
witness files, and manifest. The comparison reconstructs absolute coordinates
using all three H5 `RefPartR` components. Pair identity is recovered from birth
order and phase-space continuity because raw H5 IDs can change after MPI
redistribution.

## Run on Merlin

Stage the case in a persistent run directory and submit through the tracked
launcher. The efficient single-A100 form is:

```bash
./timed/merlin/submit_track12_a100.sh submit \
  --run-root /path/to/run \
  --opalx /path/to/cuda-build/src/opalx \
  --mpi-ranks 1 \
  --production-time 04:00:00
```

Multi-rank runs use `mpiexec`, one rank per A100, with
`--kokkos-map-device-id-by=mpi_rank`. Plain multi-task `srun` is not used on
Merlin because it creates MPI singleton processes. Launchers pass `--info 4`
and require a nonempty `timing.dat` before marking a case complete.

Every run directory must retain:

- the exact submitted Slurm script and `%j` scheduler logs under `slurm/`;
- input decks, emitted files, fixed primary, and hashes;
- OPALX logs, `runtime_compute.txt`, and nonempty `timing.dat`;
- H5 results and the preparation/run manifests.

Fetch completed H5 and compact metadata, then analyze locally:

```bash
~/.venv-h6/bin/python timed/fetch_and_compare_a100.py \
  --remote-dir /path/to/run/production
```

## Persistent field and MPI checks

The pair-4 field/kick refinement driver is
[`timed/scan_first_kick_fields.py`](timed/scan_first_kick_fields.py). It compares
both the historical untruncated Gaussian and the authoritative component-wise
three-sigma-truncated anisotropic source.

The timed MPI witness-gather regression runs the same fixed source and newborn
witnesses on one, two, and four ranks:

```bash
~/.venv-h6/bin/python timed/run_witness_gather_mpi_regression.py \
  --opalx ../../../../build_openmp/src/opalx \
  --mpiexec /opt/homebrew/bin/mpiexec \
  --openmpi-local --force
```

It requires E, B, and first-kick agreement with the one-rank result within
`5e-4`; the current deterministic result agrees below `1.5e-15`.

## Current fine-grid result

The authoritative complete legacy-domain baseline is Merlin job 354018 on four A100-SXM4-40GB
GPUs. It used the `4096x256x128` mesh for all 1501 steps, completed normally in
25,278.83 seconds, produced a nonempty 68-line `timing.dat`, matched all 13,012
CAIN-grid samples, and had no transverse wraps.

| comparison | x relative L2 | x RMSE | longitudinal relative L2 | longitudinal RMSE |
|---|---:|---:|---:|---:|
| OPALX vs three-sigma manufactured | 0.6004% | 1.996 um | 0.0971% | 0.807 um |
| OPALX vs CAIN | 9.545% | 30.648 um | 2.9774% | 24.477 um |

The twelve OPALX/manufactured first-x-kick ratios are 0.92947--0.97490,
whereas OPALX/CAIN spans 1.505--110.39. Mesh refinement therefore validates
OPALX against the like-for-like manufactured source but does not make it
converge to CAIN's different collective-field model.

The optimized one-solve/two-field implementation ran this identical deck 1.992
times faster than the preceding two-solve version. A direct H5 comparison gives
maximum accumulated coordinate changes of 1.87 nm in x, 21.62 nm in y, and
0.257 nm longitudinally.

The local compact result is under
`timed/a100_4rank_400k_4096x256x128_twofield_816d11ff8/`. Generated plots and
large H5 files are intentionally ignored by Git.

## Coordinate and field conventions

```text
x_OPALX  = x_CAIN
y_OPALX  = y_CAIN
z_OPALX  = s_CAIN
p_OPALX  = p_CAIN / (m_e c)
```

Emission-source `R0Z` places the CAIN longitudinal coordinate around the
BeamBeam midpoint. A witness born inside a step receives only the remaining
fractional drift and Boris kick. OPALX electric-field values in H5 are V/m even
though the current metadata incorrectly labels them MV/m.

Historical generated decks at this directory level remain useful for tracing
the early investigation, but they are not current validation evidence. Use the
single `timed/` deck and its manifest for new comparisons.
