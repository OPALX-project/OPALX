# Cyclotron acceleration: 72 to 590 MeV

Date: 2026-09-06. Status: approved; first-turn old-OPAL baseline established.
The first OPALX accelerating turn is implemented and compared; energy-stop and
longer-run milestones remain pending.

## OPALX implementation checkpoint (2026-09-06)

- Added Fields/CyclotronRFProfile.h: strict profile reader, cubic Hermite values
  and derivatives, shared host/device kick formulas with explicit SI contracts.
- RFCavity SINGLEGAP now loads radial profiles through a separate branch; static
  parameters only, SI gap bounds/width, standard Cartesian element poses.
- ParallelTracker has a deliberately restricted host event path for one
  unpolarized median-plane proton, one container, one rank, without restart.
  Both particle and reference use the same spatial-sector Boris integration,
  40-bisection plane localization, chronological kicks and remainder substeps.
  Particle coordinates are transformed to/from its moving frame. Host mirrors
  are intentional for this first one-particle milestone, not a GPU performance path.
- Only sectors and gaps are accepted; leaving magnetic field support is an error.
  Continuous field accumulation does not receive the RF impulse. The nominal
  coasting threader's path index is bypassed for this specialized field path.
- Initial end-to-end run exposed and corrected a moving-frame/lab-frame mismatch.
  Final 2880-step case: K=73.742318682136 MeV; old LF2=73.742310921580 MeV;
  maximum matched-time position error 0.758156 micrometres, momentum 2.648e-7;
  particle-reference offset <2.94e-12 m; five expected kicks.
- 5760-step case: K=73.742318288388 MeV; old LF2=73.742321726258 MeV;
  maximum position error 0.507069 micrometres, momentum 1.751e-7;
  particle-reference offset <6.20e-12 m. OPALX energy changes by 0.394 eV.
- First-turn regression guards: 2 micrometres maximum position error, 5e-7
  normalized momentum error, 20 eV final energy difference, 10 nm particle offset.
  These cover observed localization differences and ASCII precision, not a claim
  of asymptotic convergence. Energy gain is not constrained to be positive at every
  gap: the third harmonic decelerates in the supplied reference.
- Added six TestCyclotronRF cases (interpolation/endpoints, malformed files,
  energy/phase/transit factor, focusing/device parity, zero/invalid kicks, cloning
  and separation from continuous fields). TestCyclotronRF, TestRFCavity,
  TestCyclotronSector, TestOrbitThreader pass. Build uses -j 10.
- Full single-rank coasting DT study passes unchanged after RF changes:
  matched LF2 maximum position error remains 8.5256 nm.
- Doxygen kernel check passes. Physics/user manual pages in opalx-manual updated,
  validation and both HTML renders pass. GPU execution not checked.
- Reproducible generator/comparison: compare_acceleration_turn.py. Results under
  acceleration-opalx/dt2880-final and dt5760-final; DT/4 run in progress under
  dt11520-final. Original legacy inputs are preserved.

DT/4 review passed: 11520 nominal steps, DT=10.282987825 ps; directed return
at step 11523, t=118.490868707 ns. OPALX K=73.742315297581 MeV versus old
LF2 K=73.742317224502 MeV (difference -1.927 eV). Maximum matched-time position
error is 0.186150 micrometres and momentum error 6.0442e-8; particle/reference
offset is below 1.21e-11 m. Five expected gap events. The finer timestep changes
the sampled return endpoint; no exact-return interpolation is claimed for TURNS.
The generated comparison.png was visually reviewed. First accelerating-turn
milestone is complete; the final field-solver guard passed the DT smoke test
(acceleration-opalx/dt2880-verified), with unchanged metrics. TestCyclotronRF
and TestRFCavity were rerun and passed after that guard was added.
Particle E/B diagnostic views are not populated by this event path; future field
diagnostics must sample explicitly. The event timestep must resolve individual
crossings (no unobserved recrossing within a step).

Next: extend to 15 turns. EKINSTOP and 590 MeV completion
are not implemented yet. Independently verify historical trim conversion before
using the outer-radius comparison. No source/manual commits or pushes in this turn.

## First-turn baseline (2026-09-06)

User refined the sequence: establish one accelerating turn before expanding.
Added run_old_acceleration.py, which creates isolated inputs and runs one MPI rank.
The original cyclotron2.in remains unchanged. Parser adaptations add REAL, rename
the reserved Ring object, update distribution/output syntax, and replace obsolete
inline trim attributes with the named mirrored coil from supplied cyclotron1.in.
That historical inline-to-named trim equivalence still needs independent checking
before high-energy validation; its effect at this launch radius is negligible.
Initial attempts failed on legacy parser syntax; corrected generated inputs run.

- LF2, 2880 steps/nominal turn: DT=41.1319513 ps, directed return at step 2881,
  t=118.5011517 ns, K=73.74231092 MeV, r=2.057997398 m.
- LF2, 5760: DT=20.56597565 ps, return at step 5762 at the same sampled time,
  K=73.74232173 MeV, r=2.057993778 m.
- RK4, 2880: K=73.74231053 MeV, r=2.057994287 m.
- Five first-turn events: RF1, RF3, RF2 (third harmonic, decelerating), RF4, RF0.
  Approximate gains: +0.460019, +0.465394, -0.124854, +0.472531, +0.469220 MeV.
- The bounded reference run has a 10% margin to bracket the directed crossing;
  first-turn.csv is truncated at that crossing. Full logs include a sixth kick
  belonging to the next turn. Energy estimates use rounded old-OPAL orbit data.
- Refinement changes return energy by 10.8 eV and radius by 3.62 micrometres.
  This is a baseline check, not a demonstrated asymptotic convergence order.

Results are under acceleration-reference/lf2-2880-v5, lf2-5760 and rk4-2880;
subsequent validated script runs also record/assert the five first-turn kicks.
Next: implement and compare OPALX's first accelerating turn before extending
to 15 turns and ultimately 590 MeV. Resolve historical trim conversion before
relying on the high-radius comparison.

## Entry condition and scope

The first milestone is complete: one 72 MeV proton completes a directed turn,
with single-rank old-OPAL LF2 agreement of 8.53 nm maximum position difference.
Relevant tests and timestep comparisons passed; see CYCLOTRON_ONE_PARTICLE_PLAN.md.
Uncommitted source-documentation changes still need delivery. GPU execution is
unvalidated and is not a prerequisite for this single-rank CPU milestone.

Next goal: reproduce the supplied accelerating orbit from 72 MeV until the first
RF kick reaching or exceeding 590 MeV, including the magnetic map and trim coil.
One proton, one MPI rank, no space charge, no tune calculation, no continuous
3D RF maps. Keep development on 545-add-cyclotron-coasting-beam; builds use -j 10.

## Inspected reference

Source input: sandbox/cyclotron/opal/cyclotron2.in, unmodified.

- Launch: 72 MeV, r=2.037 m, pr=-0.0164 beta-gamma, PHIINIT=110 degrees.
  Do not reuse the coasting launch state or confuse launch azimuth with RF phase.
- Four fundamental SINGLEGAP cavities: 50.65 MHz, 0.847 MV, angles
  35/125/215/305 degrees, alternating phases 204/384 degrees.
- One third-harmonic cavity: 151.95 MHz, voltage 0.847*4*0.112*1.465 MV,
  angle 260 degrees, phase 149.40 degrees.
- Profiles rffield1.dat and rffield2.dat contain 21 support points with profile
  values and derivative data. These are voltage profiles, not 3D field maps.
- Offset gap planes: PDIS=416 mm for the fundamental and 452 mm for the harmonic;
  finite transit-time widths are 300 and 250 mm respectively.
- Magnetic map bfield.dat; legacy trim parameters TCR1=4350 mm, TCR2=4470 mm,
  MBTC=14e-3 and SLPTC=6.0. Verify their legacy conversion in the baseline.
- Existing MAXSTEPS=720*15 covers only 15 nominal turns and has no energy stop.

Legacy source inspected under /tmp/opal-2022.1-src:
RFCavity::getMomentaKick and spline in src/Classic/AbsBeamline/RFCavity.cpp;
checkGapCross, RFkick and gapCrossKick_m in
src/Algorithms/ParallelCyclotronTracker.cpp. Old tracking detects an oriented
offset-plane crossing, restores the previous state, advances to an estimated
crossing time, kicks, and advances the remaining step. Its crossing-time estimate
uses distance divided by speed, not a general root solve. Compatibility differences
must be measured rather than silently treating these algorithms as identical.

## Implementation sequence

1. Establish the old-OPAL accelerated baseline.
   Create isolated reproducible inputs without changing cyclotron2.in. First run
   its 15-turn case, explicitly selecting LF2, and also RK4 for comparison. Verify
   launch coordinates, actual trim model, RF time origin, phase signs, energy gains,
   and cavity order. Extend a bounded run sufficiently to bracket 590 MeV; do not
   assume the supplied settings reach that energy within the magnetic aperture.
   Record the first kick crossing 590 MeV and its actual post-kick energy.

2. Add the legacy radial-profile RF model.
   Keep RFCAVITY as the user-facing element with an explicit SINGLEGAP model,
   separate from its existing continuous-field path. Port profile loading,
   interpolation and derivative conventions, voltage scaling, transit-time factor,
   and focusing rotation. Internally use SI units and beta-gamma momentum; retain
   existing public VOLT [MV], FREQ [MHz], PHI0 [degrees] conventions. New geometric
   inputs use metres and X/Y/Z/THETA/PHI/PSI poses, not legacy millimetre ANGLE/PDIS
   placement attributes. Document the exact conversion of offset gap geometry and
   the profile coordinate along the gap. Preserve legacy proton behavior; flag
   questionable charge, off-plane or focusing conventions rather than generalizing
   them without tests. No API implementation before review of this plan.

3. Integrate discrete gap events with tracking.
   Detect oriented finite-support crossings, order events by crossing time, split
   the magnetic integration step, apply one complete kick, then integrate the
   remainder. Handle endpoint crossings without double kicks and multiple events
   within a step deterministically. Use a shared kick kernel for reference and
   particle tracking so a zero-offset particle remains on the reference orbit.
   Verify chronological RF time, including the remainder substep. Resolve whether
   a legacy-compatible time estimate suffices or a localized plane intersection
   is needed using measured convergence, with differences documented explicitly.
   Keep RF impulses out of continuous Boris field accumulation. Audit OrbitThreader,
   IndexMap and field-element selection for an expanding, non-closed trajectory:
   a threaded coasting orbit or fixed circumference must not exclude outer sectors
   or cavities. Prefer conservative spatial support where necessary.

4. Add kinetic-energy termination to TRACK.
   Proposed syntax: TRACK, ..., EKINSTOP=0.590, MAXSTEPS=...; (GeV).
   The name distinguishes kinetic energy from total energy and momentum; 0.590 GeV
   follows the BEAM energy-unit convention, while diagnostics print MeV.
   Evaluate reference kinetic energy K=(sqrt(1+|p|^2)-1)*m*c^2 immediately after
   each complete RF kick. Stop at the first K>=target, preserving the full physical
   kick: no fractional kick or artificial energy clamp. Write that post-kick state,
   time, cavity, actual energy, and directed-turn count before any remainder drift.
   Initially support this contract for the one-particle/single-container RING case;
   reject unsupported multi-container and restart combinations explicitly.
   Keep MAXSTEPS as a hard safety budget; budget exhaustion or loss before target
   is an unmet target, not successful completion. Validate finite positive target
   greater than launch energy. Initially reject simultaneous explicit TURNS and
   EKINSTOP to avoid competing success criteria. Current TURNS lives on RUN, so
   moving it to TRACK is a separate compatibility decision, not implicit in this work.

5. Validate from isolated kernels to the full orbit.
   Test profile values/derivatives, support edges, malformed maps, phase/voltage
   signs, zero voltage, transit-time limit, focusing and coordinate transforms.
   Test directed crossings, offset planes, no double counting, multiple events,
   reference/particle agreement and energy-stop output/failure semantics.
   Compare old and new states immediately before/after each cavity during the first
   turn, then reproduce 15 nominal turns, then extend to 590 MeV. Compare energy,
   RF phase, position, momentum, turn count and crossing time, not just final energy.
   Run DT, DT/2 and DT/4 convergence on one rank only. Establish tolerances from
   profile precision, event localization and accumulated phase error; do not reuse
   the coasting 8.53 nm result as an acceleration tolerance. Retain the RF-off
   coasting regression and relevant placement/threader tests.

6. Document and deliver.
   Add Doxygen contracts, user syntax/example and physics equations to the existing
   manuals. Save reproducible inputs, scripts, compact diagnostics and plots of
   energy versus turn, accelerating orbit and old/new error. Review diffs and record
   tested backend, tolerances, legacy discrepancies and remaining restrictions.

## Acceptance criteria

- Old OPAL and OPALX both reach >=590 MeV with the same supplied physical settings,
  or an input/aperture limitation is demonstrated before changing those settings.
- OPALX stops at the first qualifying complete kick and writes the actual state.
- Single-rank comparisons show understood timestep convergence and bounded phase,
  energy and orbit differences; event counts/order agree or discrepancies are explained.
- The first coasting milestone remains passing and no existing RFCAVITY mode changes.

## Current next step

The user approved the plan with a one-turn-first sequence. The old-OPAL baseline
is now available; proceed to profile/kick tests and one-turn OPALX RF support.
