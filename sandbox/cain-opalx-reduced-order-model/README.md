# CAIN–OPALX reduced-order model

This is the reproducible workspace for timed CAIN electron/positron witnesses
in the OPALX `BEAMBEAM` element. The authoritative raw table remains
`../track-e-p/fort98.txt`; the generated inputs and JSON reports are intentionally
not committed.

The distinction between this complete 1,297 + 1,297 input dataset, the
12-particle trajectory benchmark, and the production geometry is consolidated
in [`../BEAMBEAM_PHYSICS_AND_VALIDATION.md`](../BEAMBEAM_PHYSICS_AND_VALIDATION.md).

## Explicit emitted-file contract

The converter writes

```text
particle_count
x y z px py pz birth_time
...
```

with the following semantics:

- `x`, `y`, and `z` are metre offsets from `EMISSIONSOURCE::R0`;
- `px`, `py`, and `pz` are normalized momenta `p/(m_e c)`;
- `birth_time = ct/c` is a second offset from `EMISSIONSOURCE::T0`;
- the tracker birth event is therefore exactly
  `t_birth = T0 + birth_time`, `R_birth = R0 + (x,y,z)`.

The named header distinguishes this format from the retained legacy old-OPAL
positional format `x px y py t pz [bin]`. For compatibility with the first
prototype converter, `t` is accepted as an alias for `birth_time` only when a
named `z` column is also present.

CAIN species 2 is written to the electron file and species 3 to the positron
file. The source contains 1,297 paired rows of each species. Paired electrons
and positrons have identical `(ct,x,y,z)` values. The CAIN weight is recorded in
the conversion manifest but is not put into the particle file because these
containers are passive witnesses.

## Tracker and field contract

The OPALX deck fixes the container order:

| Container | Beam | Role |
|---:|---|---|
| 0 | primary electrons | deposits charge and receives the collective field |
| 1 | CAIN electrons | passive witness; gathers the source field |
| 2 | CAIN positrons | passive witness; gathers the source field |

`PartBunch::computeSelfFields()` and `BinnedFieldSolver` scatter only container
0. `BeamBeamInteraction::AfterEmission` then gathers that already-solved mesh
field onto containers 1 and 2. Thus a particle born inside a tracker step sees
the primary field in its creation step, while neither witness species
contributes to the solve.

The post-run validator checks that the BeamBeam field diagnostic's source
charge remains equal to container 0 while its all-container particle count
grows from 10,000 to 12,594. It also checks every witness H5 dump against the
exact cumulative CAIN birth schedule, including species code and charge sign.

Each newborn receives `dt_particle = t_step_end - t_birth`. Its initial
leapfrog half drift and Boris kick therefore use only the remaining fraction of
the birth step. The focused `TestEmittedFromFile` regression checks this and the
preservation of the signed CAIN `z` coordinate on one and multiple MPI ranks.

## Run

Build the focused test once:

```bash
cmake --build build_openmp -j8 --target TestEmittedFromFile
```

Then convert, validate, and run the focused single-rank regression:

```bash
sandbox/cain-opalx-reduced-order-model/run_pipeline.sh
```

To also run the full reduced-order OPALX deck locally:

```bash
RUN_OPALX=1 sandbox/cain-opalx-reduced-order-model/run_pipeline.sh
```

The deck sets `T0 = s_IP/(beta*c)`, so CAIN time zero coincides with the primary
centroid at the element midpoint. It sets `BBRIGID=TRUE`, zero primary momentum
spread/divergence, and `GREENSF=INTEGRATED`. This stage validates input timing,
container wiring, and witness trajectories; it does not yet establish physical
convergence of the mirroring onset or field mesh.
