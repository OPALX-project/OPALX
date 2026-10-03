# Fixed-z BeamBeam redesign status

## Goal

Validate the first incremental redesign before running all 1,297 CAIN electrons and 1,297 CAIN
positrons:

- solve collective fields only while the complete primary is inside the IP +/- 10 mm window;
- recompute transverse mesh bounds from the primary and active witnesses each field step;
- keep witnesses passive in deposition;
- keep primary container 0 alive after the window, but disable subsequent self-field solves;
- remove `RETIRE_TIME`.

## Implementation decisions

- The 32 cm BeamBeam element remains the physical tracking/aperture region.
- The 20 mm field window is a fixed internal model constant, not a new input attribute.
- Existing `PartBunch::computeBoundsForFieldSolve()` supplies the all-container bounds and box
  increment. Witness positions are temporarily mapped to the source-field frame for this operation.
- The interaction becomes active after the complete incoming primary is inside the upstream field
  boundary and becomes completed after its tail passes the downstream boundary.
- Completion suppresses default self fields but does not delete or deactivate any container.
- The existing `sandbox/track-e-p` inputs and results remain unchanged as the baseline.
- `sandbox/note/bb-note.tex` is intentionally not updated until physics comparisons pass.

## Checks completed locally

- `TestBeamBeam`: 15/15 passed on one rank.
- `TestBeamBeam`: 15/15 passed on two MPI ranks.
- `TestBinnedFieldSolver`: 10 passed, one existing four-rank-only test skipped.
- Primary-only lifecycle smoke: parsed and completed 1,607 steps. Diagnostics showed
  `Inactive -> Active -> Completed`, while `source_active=TRUE` remained true.

## Merlin6 one-A100 validation (commit `af1a1bbd6`)

All run directories are below
`/psi/home/adelmann/opalx-dev/runs/beambeam-af1a1bbd6/`; each contains its input, Slurm script,
logs, manifest, timing file, and H5 output.

- Job 354057 used the staged 1 fs input and was intentionally cancelled after about 15 minutes.
  OPALX was still serially orbit-threading the low-energy witness reference orbits on one CPU core;
  the A100 was idle. The preserved run directory is `full-witness-lifecycle-1a100`.
- Job 354058 used a coarse staged schedule. It loaded all 1,297 electrons and 1,297 positrons and
  activated BeamBeam, but stopped after the first 567-step segment. A witness orbit-map
  out-of-bounds condition leaves `globalEOL` set, which terminates tracking at a staged-segment
  boundary. The run directory is `full-witness-lifecycle-coarse-1a100`.
- Job 354059 used one scalar 608-step, 1 ps segment. It completed with exit code zero, retained both
  1,297-particle witness containers, and executed 67 BeamBeam self-field steps. It did not reach the
  BeamBeam `Completed` state. The run directory is `full-witness-lifecycle-scalar-1a100`.
- Job 354060 extended the scalar run to 630 steps. It completed with exit code zero, but BeamBeam
  remained `Active`. At H5 global step 589 the primary bounds were approximately
  `z=[-1.573, +1.730]` mm in the source frame. By global step 599, particles crossing the fixed
  downstream field-layout boundary had been remapped and the bounds jumped to approximately
  `z=[-20.240, +1.091]` mm. By step 609 all remaining source particles occupied approximately
  `z=[-22.600, -19.669]` mm. The corrupted tail can never satisfy the downstream completion test;
  995 of 1,000 primary macroparticles remained in the final H5 dump. The run directory is
  `full-witness-lifecycle-scalar630-1a100`.

The final scalar job used 211.00 s wall time, of which 200.6 s was `OrbThreader`. The BeamBeam
window itself used 0.7581 s, including 0.6144 s in self fields. These timings are smoke-test values
for a 16^3 mesh and must not be used as production performance estimates.

The same 630-step scalar case was reproduced on the Mac with the Release/OpenMP build, one MPI
singleton, four OpenMP threads, and `--info 4`. Total wall time was 12.05 s and `OrbThreader` used
8.7419 s. A matched primary-only case used only 0.001542 s in `OrbThreader`. The excess therefore
comes from constructing one orbit threader for each un-emitted witness container. Their initialized
reference momenta were `(0,0,1.95695e-6)` in beta-gamma units, so these nearly stationary reference
orbits consume the configured step budget without producing a meaningful witness path. This also
explains why their maps contain only the upstream drift and later generate out-of-bounds queries.
All subsequent OPALX validation runs must use `--info 4`.

A temporary `BEAM, ORBITTHREADER=FALSE` control now lets a secondary container reuse container
0's design-orbit element map without constructing an independent threader. With the flag on both
witness beams, the matched full-witness Mac case loaded both 1,297-particle files, completed all
630 requested steps, and reduced `OrbThreader` to 0.001553 s. Total wall time fell to 3.10 s. The
default is `TRUE`, and container 0 is required to keep an independent design threader. The focused
`TestBeam` and `TestBeamBeam` suites pass (18/18 and 15/15 respectively).

A coarse/fine/coarse Mac integration check (`567 x 1 ps`, `1000 x 1 fs`, `40 x 1 ps`) also
completed all 1,607 steps with `--info 4`; it no longer terminated at the first segment boundary.
`OrbThreader` used 0.290736 s and total wall time was 16.44 s. This confirms the shared primary map
works with staged time steps locally and keeps the primary-only minimum-step threading below one
second.

## Merlin6 one-A100 validation (commit `7b2a64523`)

Job 354062 ran the same reduced staged schedule (`567 x 1 ps`, `1000 x 1 fs`, `40 x 1 ps`) on one
A100 and one MPI rank with four CPU threads and `--info 4`. The 16^3 smoke mesh used 1,000 primary
macroparticles and both complete 1,297-particle witness files. The job completed all 1,607
integration steps with exit code zero in 44.88 s application wall time (46 s Slurm elapsed). The
timing file is nonempty. `TIntegration1`, `TIntegration2`, and `External field eval` each have 1,607
calls; `BB self fields` has 1,066 calls. Both witness gathers continued through the final step.

The reported `OrbThreader` interval fell from 177 s in the earlier 630-step full-witness A100 run
to 4.607 s and its count remained one. This demonstrates that the two nearly stationary witness
reference orbits are no longer threaded, but the number is not a pure primary-threading time:
`ParallelTracker` currently starts the timer before the initial `writePhaseSpace()` call and stops
it after `fieldSolver->orbitThreadersReady()`. Thus it also contains initial H5 output and
first-use field-solver/GPU setup. The timer boundaries need to be narrowed before using the
nominal `< 1 s` criterion for the primary `OrbitThreader` itself.

The self-contained run directory is
`/psi/home/adelmann/opalx-dev/runs/beambeam-7b2a64523/orbitthreader-staged-1a100`.

## Track-12 old/new comparison attempt (commit `7b2a64523`)

The preserved one-A100 job 353986 at commit `8da5b9e83` was selected as the old trajectory
baseline: 400,000 deterministic primary macroparticles, six timed electrons, six timed positrons,
a `1024 x 128 x 128` mesh, and 1,501 steps of 6.0041537 fs. Its deterministic input hashes were
reused byte-for-byte. The new deck differed only by removing the obsolete `RETIRE_TIME` attribute
and setting `ORBITTHREADER=FALSE` on the two witness beams. All jobs used `--info 4`.

The three-step smoke job 354063 failed on the A100 with `cudaErrorIllegalAddress` in
`ippl::ParticleAttrib::gather()` during the second witness-field gather. Production job 354064 was
dependency-blocked and cancelled without running. Isolation job 354065 restored independent
witness OrbitThreaders and failed at the identical call, proving that shared OrbitThreader maps are
not the cause.

Two-step first-kick job 354066 completed and wrote H5 output. Direct comparison with old job
353986 showed that all stored witness phase-space components (`x`, `y`, `z`, `px`, `py`, `pz`) are
finite in the old result but `NaN` in the new result immediately after the first kick, for both the
electron and positron. The second-step CUDA fault is therefore downstream damage: the next gather
receives NaN witness coordinates. The failure is in the new fixed-z/dynamic-x-y binned BeamBeam
field path and must be localized before a full trajectory or production-population run.

The failure was localized to the dynamic transverse mesh shrinking to the particle bounds. In the
Track-12 case this produced approximately 11.6 nm x 93.7 nm transverse cells, while the Lorentz
transform stretched the longitudinal cell to approximately 75.7 mm in the source rest frame. The
extreme aspect ratio makes the integrated Green-function expression evaluate `log(a+r)` with
`a+r == 0` after floating-point cancellation for some cells, producing the first-kick NaNs. Using
the standard Green function avoided the NaNs but changed the first electron kick from
`1.71198e-3` to approximately `-4.05`, so it is not a physically acceptable workaround.

The local fix treats the configured BeamBeam `APERTURE` as the minimum transverse field domain.
The domain still expands whenever the aggregate primary-plus-witness bounds exceed the aperture,
but it no longer shrinks inside it. With the Track-12 `RECTANGLE(2.4e-3,2.4e-4)` aperture this
restores the 2.4 mm x 0.24 mm transverse domain and keeps the integrated Green function finite.
A focused unit test verifies both the aperture floor and expansion for an outlying witness.

A three-step local Release/OpenMP smoke run with the exact 400,000-particle and 12-witness inputs
completed two witness gathers without NaNs or a crash. The first electron `px` after the kick was
`1.73676176e-3`, versus `1.71197961e-3` in old job 353986 (a 1.45% change). The remaining difference
is not yet an error measure: the old run used the 8 mm element-length longitudinal field domain,
whereas the redesigned model deliberately uses the fixed IP +/-10 mm domain. A full A100 trajectory
comparison is required after this fix is committed and pulled to Merlin6.

The aperture-floor implementation must not be committed in this form: `APERTURE` is the physical
acceptance/loss boundary, not a field-mesh sizing parameter. The physical and numerical roles need
separate semantics.

### Particle-bounds mesh experiment

A separate three-step local run removed `APERTURE` so that the uncommitted aperture-floor code could
not affect field sizing, and changed only the field mesh from `1024 x 128 x 128` to
`256 x 256 x 128`. The same deterministic 400,000-particle source and integrated Green function
were used with `--info 4`. The particle-derived transverse domain was approximately
`11.90 um x 11.90 um`, giving lab-frame cells of `46.66 nm x 46.66 nm x 157.48 um`. After the
gamma=480.45 rest-frame stretch, the longitudinal cell was approximately 75.66 mm, or about
1.62 million times either transverse cell width.

The program completed all three steps and both witness gathers, but every stored witness component
(`x`, `y`, `z`, `px`, `py`, `pz`) was NaN immediately after the first kick. Thus balancing the x/y
mesh at `256 x 256` does not cure the integrated-Green failure. Successful process completion is
not evidence of valid physics in this case; output finiteness must be checked explicitly. Local
`BB self fields` time was 131.1 s for the three solves. The temporary run is
`/tmp/track12-dynamic-256x256x128-rerun`.

### Primary-only field dump

The reproducible diagnostic under `sandbox/track-e-p-new/field-diagnostics` uses the same fixed
400,000-particle primary, a `256 x 256 x 128` field mesh, and a passive `101 x 201` x-z probe sheet.
The probes cover `+/-3 sigma_x` and `+/-4 sigma_z`, do not deposit charge, and the mirrored source
is disabled. `EBDUMP=TRUE` stores the gathered `Ex`, `Ey`, `Ez`, `Bx`, `By`, and `Bz` components in
the probe H5 file. The run and plotting commands use `--info 4` and the existing reduced-order
manufactured Gaussian evaluator.

With `GREENSF=INTEGRATED`, zero of 20,301 probe locations had finite E or B: all six dumped field
components were NaN at every probe immediately after the first solve. An otherwise identical
`GREENSF=STANDARD` control produced finite E and B at every probe. This verifies that source
deposition, passive-probe gathering, and H5 field output can complete with finite values, and
localizes where the invalid values first become visible rather than blaming the integrated solver
itself. The old validated Track-12 run also used `GREENSF=INTEGRATED` successfully.
The standard kernel is not a physical substitute on this mesh: relative to the component-wise
three-sigma-truncated manufactured Gaussian, both relative L2(E) and relative L2(B) were
approximately `1.6847e4`.

The integrated and standard diagnostic artifacts, including PNG/PDF field visualizations, CSV
samples, JSON summaries, full logs, and H5 files, are in the ignored directories
`field-diagnostics/outputs/0sigma` and `field-diagnostics/outputs/standard-control` respectively.

Temporary device-side reductions then checked the complete distributed field immediately after
`runSolver(true)`, after the two-field lab-frame finalization, and after passive-probe gathering.
For the integrated kernel, all `256 x 256 x 128 x 3 = 25,165,824` rest-frame E components were
already non-finite at solver return. Finalization produced the same number of non-finite components
in each of lab E and lab B. Gathering then transferred `20,301 x 3 = 60,903` non-finite E and B
components to the probe container. In the otherwise identical standard-kernel control, every
counter was zero at all three stages. This proves that the first NaNs are produced inside the
integrated-Green field solve; neither BeamBeam finalization, witness gathering, nor the Boris kick
creates them. The causal change is the redesigned field-domain geometry passed to the unchanged
integrated solver. The old Track-12 domain was `2.4 mm x 0.24 mm x 8 mm` on a
`1024 x 128 x 128` mesh. The new particle-bounds case is approximately
`11.9 um x 11.9 um x 20 mm` on a `256 x 256 x 128` mesh. After the gamma=480.45 longitudinal
stretch, the largest old cell aspect ratio was approximately `1.6e4`, whereas the new case is
approximately `1.62e6`. The previously tested temporary restoration of the old transverse domain,
even with the new 20 mm longitudinal window, produced finite integrated-Green fields. Therefore the
fix belongs in the new BeamBeam field-domain policy; modifying IPPL is neither justified nor the
next step.

The full-field reductions are currently enabled only when the active `Inform` output level is at
least 4. They are useful during this diagnosis but traverse the complete mesh and should be removed
or placed behind an explicit diagnostic switch before production benchmarking.

### Transverse-domain control matrix

A controlled matrix kept the 20 mm longitudinal window, the integrated Green function, the
deterministic 400,000-particle primary, and the passive `101 x 201` probe sheet unchanged. First,
the approximately `11.9 um x 11.9 um` particle-tight transverse domain was tested with
`64 x 64 x 128`, `128 x 128 x 128`, and `256 x 256 x 128` meshes. Every rest-frame E component was
already non-finite at solver return in all three cases. Merely reducing the transverse mesh count
therefore does not repair the redesigned geometry.

For diagnosis only, the uncommitted aperture-floor mechanism was then used to impose successively
larger transverse domains while keeping the lab-frame transverse cell width approximately
`0.394 um`. This use of `APERTURE` is not the proposed design and must not be committed. All three
controlled cases were finite and agreed closely with the component-wise three-sigma-truncated
manufactured Gaussian:

| Mesh | Transverse domain | Transverse cell | Relative L2(E) | Relative L2(B) | Local BB solve |
|---|---:|---:|---:|---:|---:|
| `64 x 64 x 128` | `25 um x 25 um` | `0.3968 um` | `1.9662e-2` | `1.9662e-2` | `1.125 s` |
| `128 x 128 x 128` | `50 um x 50 um` | `0.3937 um` | `1.9562e-2` | `1.9562e-2` | `4.403 s` |
| `256 x 256 x 128` | `100 um x 100 um` | `0.3922 um` | `1.9515e-2` | `1.9515e-2` | `22.050 s` |

After the `gamma=480.45` rest-frame stretch, these cases have a largest cell aspect ratio of about
`1.9e5`, within the previously exercised manufactured-solution regime. The results isolate the
failure to the new field-domain geometry policy: a particle-tight transverse box is incompatible
with the fixed 20 mm longitudinal window for this source and mesh. They do not justify changing
the integrated-Green implementation. The production design needs a dedicated numerical field
domain (primary envelope plus an explicit numerical margin, expanded by the active witness
envelope) that remains distinct from the physical aperture.

An off-diagonal matrix then varied domain width and mesh count independently. In the primary rest
frame, the 20 mm lab window gives `Delta z' = 75.66 mm`. The results establish that field
finiteness follows the transformed cell aspect ratio rather than the distance to the transverse
boundary:

| Transverse domain | Mesh | `Delta x` | `Delta z'/Delta x` | Relative L2(E/B) | Result |
|---:|---:|---:|---:|---:|---|
| `12 um x 12 um` | `32 x 32 x 128` | `0.3871 um` | `1.95e5` | `1.9352e-2` | finite |
| `25 um x 25 um` | `64 x 64 x 128` | `0.3968 um` | `1.91e5` | `1.9662e-2` | finite |
| `25 um x 25 um` | `128 x 128 x 128` | `0.1969 um` | `3.84e5` | unavailable | all NaN |
| `50 um x 50 um` | `64 x 64 x 128` | `0.7937 um` | `9.53e4` | `4.3792e-2` | finite |

The `12 um` case shows that a tight boundary is not by itself responsible: it is finite at an
acceptable cell ratio and has the lowest measured error in this matrix. Conversely, the same
`25 um` boundary changes from finite at 64 cells to entirely non-finite at 128 cells. The current
evidence places the transition between aspect ratios `1.95e5` and `3.84e5`. The larger `50 um`
case remains finite but its coarser transverse cells increase the relative L2 error to 4.38%, so
finiteness and field accuracy must be treated as separate constraints.

### Dynamic aspect-limited domain implementation

The replacement field-domain policy no longer reads `APERTURE`. Each active step transforms the
emitted witness containers into the primary-source frame and reuses the existing distributed
all-container bounds. `BBOXINCR` supplies the ordinary deposition/gather margin. In each transverse
dimension, the interaction then enforces

`gamma * Delta z / Delta x <= 1.5e5`

by symmetrically expanding the numerical domain when necessary. The value `1.5e5` is conservative
relative to the measured finite/failing bracket. The retained transverse bounds are unioned with
the previous step, so the field domain can expand with rapidly diverging witnesses but cannot
shrink during an active BeamBeam passage. The state is reset on entry and exit. `--info 4` reports
the aggregate and per-container bounds, retained domain, rest-frame longitudinal cell, and both
aspect ratios.

A primary-only `256 x 256 x 128` run without an aperture floor selected a
`128.63 um x 128.63 um` transverse domain, was completely finite, and had relative L2(E/B)
`2.4073e-2` against the three-sigma-truncated manufactured source. The full-mesh finiteness
reductions used during diagnosis were subsequently removed because they would add two complete
field traversals to every `--info 4` production step. The temporary witness-component reduction
was also removed after CUDA rejected its extended device lambda inside the private interaction
method; H5 trajectory validation remains the authoritative finiteness check.

Two three-step local Track-12 runs used `1024 x 128 x 128` and `512 x 64 x 128`. Both selected
`Delta x = Delta y = 0.5044 um`, produced finite first and second witness kicks, and agreed in the
first electron `px` kick to `8.9e-10` relative despite the factor-two difference in transverse
domain extent. The new first kick is `2.5593662e-3`, 3.95% above the old authoritative
`4096 x 256 x 128` result (`2.4619981e-3`) but closer to the manufactured kick: the new ratio is
96.63%, versus 92.95% for the old fine result. The stable same-resolution comparison shows that the
difference is not a transverse boundary artifact. A complete A100 trajectory is still needed to
compare the redesigned 20 mm model over all six pair birth times.

Remote artifacts are preserved under:

- `/psi/home/adelmann/opalx-dev/runs/track12-compare-8da5b9e83-vs-7b2a64523-1a100-1024x128x128`
- `/psi/home/adelmann/opalx-dev/runs/track12-isolate-independent-oth-7b2a64523`
- `/psi/home/adelmann/opalx-dev/runs/track12-first-kick-compare-7b2a64523`

## Blocking design issue

The fixed 20 mm field layout cannot own the physical primary particles unchanged over the requested
lifecycle. Activation waits until the complete source is inside the upstream boundary, whereas
completion waits until its tail exits the downstream boundary. Between the head exiting and the
tail exiting, some physical primary particles are outside the fixed field layout. Updating the
physical container against that layout remaps or removes those particles and corrupts the physical
extent used by the completion test.

Before production, separate the physical primary particle lifecycle from the population deposited
onto the fixed field mesh. The deposited population must be clipped or copied without attaching the
physical primary container to a layout that excludes some of its particles. The field and witness
gather semantics near the fixed-window boundaries then require a focused regression test.

## Remaining validation

1. Run the complete 1,501-step `4096 x 256 x 128` Track-12 case on A100 and compare the trajectory,
   first kicks, and final momenta with the legacy fine result and manufactured solution.
2. Use the `--info 4` geometry history to confirm that every emitted witness remains inside the
   retained source-frame domain and measure when witness expansion overtakes the aspect floor.
3. Resolve the fixed-window/physical-primary layout ownership issue above.
4. Check one-rank and multi-rank agreement for the new non-shrinking domain history.
5. Run the approximately 1,297-electron plus 1,297-positron input and assess the resulting
   transverse resolution before treating it as a production configuration.
6. Update `sandbox/note/bb-note.tex` after the complete Track-12 result is accepted.
