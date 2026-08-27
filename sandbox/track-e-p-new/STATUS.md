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

1. Resolve the fixed-window/physical-primary layout ownership issue above.
2. Verify `ORBITTHREADER=FALSE` gives `OrbThreader < 1 s` on Merlin6/A100.
3. Confirm per-step transverse mesh bounds contain both witness species in the source frame.
4. Compare witness trajectories and final momenta with the previous validated OPALX result,
   manufactured solution, and CAIN comparison.
5. Check one-rank and multi-rank agreement.
6. Update `sandbox/note/bb-note.tex` only after these results are accepted.
