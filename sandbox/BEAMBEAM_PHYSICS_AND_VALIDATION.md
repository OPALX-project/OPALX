# BeamBeam physics model and validation

This is the authoritative description of the gamma--gamma BeamBeam model and
its validation. Historical derivations, superseded scans, and generated
publication assets are retained under [`attic/`](attic/README.md). The
chronological engineering record remains in
[`../BEAMBEAM_REDESIGN_STATE.md`](../BEAMBEAM_REDESIGN_STATE.md).

## 1. Problem statement

OPALX models a gamma--gamma interaction region in which two high-energy
electron primary beams pass through one another. Their photons collide near the
interaction point (IP) and create low-energy electron--positron pairs. The
primary beams are present when these pairs are born, so the primary collective
electric and magnetic fields can strongly alter the pair trajectories.

The present model does not generate pairs. It reads CAIN electron and positron
phase space, including each particle's individual birth time, and transports
the particles through the primary-beam fields. These low-statistics pair
particles are passive **witnesses**: they sample the fields but do not deposit
charge and do not contribute to a field solve.

Two geometries serve different purposes and must not be mixed:

| case | BeamBeam element | witness population | purpose |
|---|---:|---:|---|
| Track-12 validation | 8 mm long; `2.4 mm x 0.24 mm` rectangular field domain | 6 electrons + 6 positrons | exact-timed trajectory, field, mesh, MPI, and performance validation |
| Production | 32 cm long; 15 cm radius (`CIRCLE(0.30)`) | 1,297 electrons + 1,297 positrons | complete CAIN pair transport |

The IP is always the BeamBeam element midpoint. There is no independent `IP_S`
input. In the production element the IP is therefore 16 cm from either element
edge; its global path coordinate is determined by placement of the element.

The primary source is an anisotropic bunch in the laboratory frame. The
reference parameters used by the current reduced-order validation are

```text
primary kinetic energy       245 MeV
electrons per primary bunch  1.25e10
sigma_x = sigma_y            1.944325075701 um
sigma_z                      0.6 mm
source support               component-wise |x_i| <= 3 sigma_i
```

The old spherical rest-frame Gaussian is not the physical reference model. It
is archived only to retain provenance.

## 2. OPALX implementation

### Architecture and ownership

`BeamBeam` follows the same element-level ownership model as other OPALX
elements: it contains immutable input and placement configuration. Its mutable
per-run state belongs to `BeamBeamInteraction`, created through the generic
`ElementBase::createInteraction()` hook and managed by
`ElementInteractionManager`.

`ParallelTracker` dispatches generic `SelfField`, `AfterEmission`, and
`Diagnostics` phases through `ElementInteractionContext`. It does not contain a
BeamBeam-specific branch. Currently only `BeamBeamInteraction` uses this
context; every other element inherits the null `createInteraction()` behavior
and continues through the ordinary visitor/`apply()` mechanism.

The runtime interaction owns the window lifecycle, fixed field-domain state,
primary-copy state, witness gathering, diagnostics, and timers. Particle
containers remain owned by `PartBunch`.

### Active-step data flow

The maintained production and Track-12 configurations use exactly one source
bin and no cathode/image-charge model. They therefore take the guarded
two-field path:

```text
tracked primary particles
    -> scatter container 0 charge to rho
    -> one integrated-Green open-boundary field solve
    -> physical-primary E_p and B_p
    -> reflect the solved mesh about the IP to construct E_v and B_v
    -> accumulate E_p + E_v in E_m and B_p + B_v in B_m
    -> gather the combined field to the tracked primary when BBRIGID=FALSE
    -> emit witnesses due during this step
    -> gather the same combined field to witness containers 1 and 2
    -> fractional Boris update from each exact birth time to the step end
```

Only `E_m` and `B_m` persist on this fast path. Reflection and accumulation are
device kernels; the full fields are not staged through host memory. Generic
multi-bin and image-charge configurations retain their separate temporary-field
path.

The virtual primary is a reflected field/source representation rather than a
second independently tracked particle bunch. With the production default
`BBRIGID=FALSE`, the tracked physical primary receives the combined collective
field. With `BBRIGID=TRUE`, it still deposits charge and drifts, but the
BeamBeam collective field is removed from its momentum update; this is the
controlled rigid-source validation mode. External fields are unaffected.

Mirroring is activated by the configured `COPY_TIME`. It is not assumed to be
physically converged at a prescribed number of longitudinal sigmas. Likewise,
the source is not retired at `3 sigma_z`; retirement remains an explicit
`RETIRE_TIME` choice. The adequacy of both times must be checked by convergence.

### Timed witnesses and field participation

The explicit emitted-file format is

```text
particle_count
x y z px py pz birth_time
...
```

Positions are metres relative to the emission-source `R0`; momenta are
dimensionless `p/(m_e c)`; and `birth_time` is seconds relative to the source
`T0`. A particle born inside a step is advanced only over

```text
dt_particle = t_step_end - t_birth.
```

The container contract is fixed:

| container | particles | field role |
|---:|---|---|
| 0 | primary electrons | deposits source charge and, unless rigid, receives the collective field |
| 1 | CAIN electrons | passive witness; gathers fields only |
| 2 | CAIN positrons | passive witness; gathers fields only |

Witness gathering uses each container's full reference-particle offset before
sampling the source-local mesh and performs an MPI owner-rank reduction. This
is required because newly emitted witnesses need not reside on the rank that
owns the sampled field cell.

### Numerics, units, and assumptions

- Field solves use the cell-integrated Green function and open boundaries.
- Source coordinates are transformed to the instantaneous beam frame for the
  solve and fields are transformed back to the reference-path frame.
- `APERTURE` arguments are full widths/diameters; `CIRCLE(0.30)` gives a 15 cm
  radius.
- Electric and magnetic fields have SI units V/m and T internally.
- The current H5 electric-field values are V/m, but the writer metadata labels
  them MV/m. This known metadata defect must be corrected separately with a
  downstream compatibility audit; analyses in this sandbox interpret the
  stored values as V/m.
- Mesh refinement, MPI reduction order, finite-source sampling, and reflected
  mesh reconstruction can change floating-point results. Deterministic fixed
  source files are used when decomposition invariance is the quantity tested.

### Current implementation evidence

| test | particles | mesh | step / steps | ranks / A100s | result | wall time |
|---|---:|---:|---:|---:|---|---:|
| Persistent 3-sigma manufactured CTest | 100,000 source; 5,151 probes | `64x64x128` | 1 fs / 2 | 1 / 0 | relative L2 `(E,B)=(2.0976%,2.1984%)` | about 5 s |
| Optimized fixed-source rank scan | 400,000 source; 5,151 probes | `64x64x128` | 1 fs / 2 | 1,2,4 / 1,2,4 | four-rank fields agree with one rank to `2.4e-16`; L2 errors about 1.22% | 4.94, 8.29, 10.50 s |
| Timed witness-gather rank regression | 100,000 source; 8 e- + 8 e+ | `32x32x64` | 6.004154 fs / 5 | 1,2,4 / 1,2,4 | rank differences in E, B, kick below `1.5e-15` | 5.85, 8.30, 9.82 s |
| One-solve timing probe | 400,000 source; 1 e- + 1 e+ | `1024x128x128` | 6.004154 fs / 10 | 1 / 1 | solves 20 -> 11; BeamBeam self-field time -45.2% | 6.62 s |
| Optimized fine Track-12, job 354018 | 400,000 source; 6 e- + 6 e+ | `4096x256x128` | 6.004154 fs / 1501 | 4 / 4 | all 13,012 samples; no wraps; normal completion | 25,278.83 s |

The extra solve in the ten-step count is solver warm-up: ten physical-source
solves plus one warm-up. The optimized full trajectory was 1.992 times faster
than the identical two-solve case and changed the final stored trajectory only
at the nanometre-scale accumulated numerical level.

## 3. Manufactured solution

The authoritative independent reference is
[`reduced-order-model/`](reduced-order-model/README.md). It represents two
unchanging, uniformly moving anisotropic Gaussian primary bunches. The IP is at
the origin, their centroids move in opposite longitudinal directions, and the
witnesses do not generate fields.

For a lab-frame density

```text
rho(x,y,z) proportional to
exp[-0.5 ((x/sigma_x)^2 + (y/sigma_y)^2 + (z/sigma_z)^2)],
```

the independent Python evaluator computes the anisotropic rest-frame Coulomb
integral and applies the exact uniform-motion Lorentz transformation to obtain
lab-frame E and B. For like-for-like trajectory comparisons the density is
truncated and renormalized independently at three sigma in x, y, and z,
matching the deterministic OPALX source-generation rule without finite-sample
noise.

The field snapshot study compares centroid separations of 3, 2, 1, and 0
`sigma_z`. Its diagnostic grid is `101x201`: odd dimensions place a sample
exactly at the IP, 101 points resolve the transverse line, and 201 points cover
the longer longitudinal interval. The persistent test uses `51x101` probes for
lower cost while retaining the exact central sample.

At three-sigma separation, the corrected physical-source comparison gave:

| probes | OPALX max E [V/m] | manufactured max E [V/m] | relative L2 E | median magnitude ratio | median direction cosine | relative L2 B |
|---:|---:|---:|---:|---:|---:|---:|
| 20,301 | `5.545751e9` | `5.553866e9` | 1.1098% | 0.997923 | 0.999970 | 1.1098% |

The Track-12 transverse refinement uses the same fixed `2.4 mm x 0.24 mm`
domain and `Nz=128`:

| transverse mesh | ranks/A100s | OPALX/truncated E | OPALX/truncated kick | runtime [s] |
|---:|---:|---:|---:|---:|
| `1024x128` | 1/1 | 0.672577 | 0.672591 | 5.69 |
| `1536x192` | 1/1 | 0.835404 | 0.835416 | 6.93 |
| `2048x256` | 1/1 | 0.912683 | 0.912698 | 7.59 |
| `2304x288` | 3/3 | 0.916806 | 0.916822 | 58.11 |
| `2560x320` | 3/3 | 0.930909 | 0.930931 | 69.51 |
| `2688x336` | 3/3 | 0.940164 | 0.940189 | 76.38 |
| `3072x384` | 4/4 | 0.971438 | 0.971466 | 121.56 |

The direction cosine is at least 0.999993. Refinement corrects field magnitude,
not direction, and the observed asymptotic behavior approaches second order.
The separate `4096x256x128` point gives an E ratio of 0.960704 because it
refines x but coarsens y relative to `3072x384x128`; the error is genuinely
anisotropic.

For the full `4096x256x128`, 1501-step trajectory, OPALX agrees with the
three-sigma-truncated manufactured trajectory to 0.6004% relative L2 in x and
0.0971% longitudinally. The twelve first-x-kick ratios span 0.92947--0.97490
with mean 0.94860. This is the principal current physics validation of OPALX.

Reproduction instructions for local analysis and remote A100 data generation
are in [`reduced-order-model/opalx/README.md`](reduced-order-model/opalx/README.md)
and [`track12particles/opalx/README.md`](track12particles/opalx/README.md).

## 4. CAIN solution

CAIN provides two distinct datasets with different roles:

| dataset | contents | role |
|---|---:|---|
| `track-e-p/fort98.txt` | 1,297 electrons + 1,297 positrons | production timed pair population and input-format validation |
| `TestParticleOrbit.dat` | 6 artificial electrons + 6 artificial positrons | exact-timed trajectory comparison on the CAIN output grid |

The conversion and emission contract for the full dataset is maintained in
[`cain-opalx-reduced-order-model/`](cain-opalx-reduced-order-model/README.md).
CAIN species 2 maps to electrons and species 3 to positrons. Paired particles
have matching creation coordinates and times. Their statistical weights are
retained in the conversion manifest but are deliberately not used for witness
deposition.

The Track-12 births occur at

```text
ct = -0.9000, -0.5994, -0.3006, 0.0000, 0.3006, 0.5994 mm
```

and CAIN samples every `1.8 um/c = 6.004153713566736 fs` through `ct=1.8 mm`.
The OPALX run begins one interval before the first birth and uses the same time
grid, so all 13,012 electron/positron trajectory samples match without time
interpolation.

CAIN and OPALX do not use the same collective-field model. CAIN is therefore a
physics comparison, not the numerical oracle for OPALX's rigid two-Gaussian
model. This distinction is borne out by refinement:

| comparison at `4096x256x128` | x relative L2 | x RMSE | longitudinal relative L2 | longitudinal RMSE |
|---|---:|---:|---:|---:|
| OPALX vs manufactured | 0.6004% | 1.996 um | 0.0971% | 0.807 um |
| OPALX vs CAIN | 9.545% | 30.648 um | 2.9774% | 24.477 um |

As the OPALX transverse mesh is refined, OPALX converges strongly toward the
like-for-like manufactured solution and not toward CAIN. The OPALX/CAIN
first-kick ratios span 1.505--110.39, while charge-conjugate OPALX kicks agree
to a maximum relative residual of `7.57e-10`. Thus timed input, species
assignment, charge sign, MPI gathering, and integration are validated; the
remaining CAIN discrepancy is attributed to the difference in physical source
model or conventions and should be investigated separately rather than tuned
away in OPALX.

The next production step is to run the 32 cm, 15 cm-radius element with the
full 1,297 + 1,297 timed CAIN witness population, then repeat the necessary
time-, aperture-, mesh-, and primary-sampling convergence checks before drawing
physics conclusions from the pair losses or trajectories.
