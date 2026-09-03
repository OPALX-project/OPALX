# RING implementation plan

## Current state

- Branch: `543-ring-modelling-capabilities`.
- PR #540 was integrated as commits `afcb81755` and `f95c0c8ca`.
- The unused `ElementType::RING` was removed in commit `e8c43e150`.
- Focused OpenMP tests pass after integration:
  - `TestSolenoid`: 10 tests
  - `TestIndexMap`: 4 tests
  - `TestOrbitThreader`: 2 tests
- Issue #543 defines `RING` as a closed-topology beam sequence with per-element ring membership.
- The first-stage golden reference is `sandbox/ring-1/orig`.
- The golden input is an explicitly positioned, expanded ISIS survey. It runs one tracking step;
  therefore its strongest oracle is lattice placement and the initial OrbitThreader trajectory,
  not multi-turn particle tracking.

## First-stage goal

Accept an OPAL input containing a `RING` sequence built from the corresponding ISIS elements
placed with `ELEMEDGE`. Build the same physical one-turn lattice and initial reference orbit as the
expanded explicit-placement reference in `sandbox/ring-1/orig`.

Out of scope for this stage:

- multi-turn particle tracking semantics;
- automatic geometric correction of a non-closing ring;
- nested RINGs or a RING embedded in a LINE;
- GPU/device access to topology metadata;
- changes to field models, integrators, MPI reductions, or numerical tolerances.

## Design decisions

1. Introduce a topology type independent of `ElementType`. The unused `ElementType::RING` enum
   value has been removed; sequence topology must not be encoded as a physical element type.
2. Store host-side membership on `ElementBase`: topology defaults to `LINE`; `RING` membership
   carries the canonical outer ring name. The copy constructor must preserve it.
3. Store the same topology and name on the root `FlaggedBeamline`, so setup code can identify the
   selected sequence without inferring topology from leaf elements.
4. Add a distinct `Ring` beam-sequence parser registered as `RING`. Share LINE parsing helpers,
   but enforce issue #543's policy: no `N*member` repetition and no `-member` reflection.
5. Preserve nested LINE containers inside a RING. Recursively clone and tag each user-visible leaf
   occurrence so registered prototypes and ordinary LINE occurrences are not contaminated.
   Synthetic `#S` and `#E` markers remain untagged.
6. Reject nested RINGs and RING-as-a-LINE-member with a parse error in this stage.
7. Materialize RING occurrences immediately before tracking, after `OpalData::update()` has
   populated the registered element prototypes. Parser-time cloning is too early and would retain
   default geometry. LINE preparation remains a no-op.

## Implementation sequence

### 1. Lock down the golden comparison

- Create a new `sandbox/ring-1` input using the compact `ELEMEDGE` ISIS lattice and `RING` keyword.
- Keep beam, distribution, time step, `MAXSTEPS=1`, element apertures, and SBEND model identical to
  the reference wherever possible.
- Add a reproducible comparison script using `~/.venv-h6` and pandas for tabular outputs.
- Compare semantic/numeric content rather than timestamps, revision strings, or file ordering.

Primary comparisons:

- `_ElementPositions.txt`: element names plus entry/body/exit coordinates;
- `_ElementPositions.sdds`: first-pass longitudinal intervals and flags;
- `_DesignPath.dat`: numeric orbit columns and active-element names over the common path range;
- VTK mesh vertices/topology after canonical ordering;
- final PNG/PDF as visual confirmation, not the sole numerical oracle.

Start with tight absolute geometry tolerances (approximately `1e-10` m and `1e-10` rad) and only
relax them if a measured floating-point-order difference is justified.

### 2. Add topology and membership metadata

- Add `BeamlineTopology { LINE, RING }` and a small membership value type near `ElementBase`.
- Add const queries and validated mutation methods.
- Initialize all elements as LINE members with no ring name.
- Copy membership in `ElementBase(const ElementBase&)`.
- Unit-test default state, validation, and clone preservation.

Expected files:

- `src/AbsBeamline/ElementBase.h`
- `src/AbsBeamline/ElementBase.cpp`
- focused `unit_tests/AbsBeamline` test

### 3. Add the RING sequence

- Add `Ring.h/.cpp` and register `new Ring()` in `Configure::makeElements()`.
- Extract the reusable definition/length/attribute parts of LINE parsing without changing LINE
  behavior.
- Give RING a policy-driven member parser that rejects repetition and reflection.
- Preserve named nested LINEs. Decide anonymous-list representation in code before enabling it;
  if preservation cannot be cleanly represented with correct dynamic length expressions, reject
  anonymous nested lists in this first stage rather than silently flattening them.
- Ensure `print()` emits `name:RING=(...)`.

Expected files:

- `src/Lines/Line.h`, `src/Lines/Line.cpp`
- new `src/Lines/Ring.h`, `src/Lines/Ring.cpp`
- `src/Lines/CMakeLists.txt`
- `src/OpalConfigure/Configure.cpp`
- new parser/structure unit tests

### 4. Isolate and tag RING occurrences

- Recursively traverse the parsed ring structure.
- Clone every physical leaf occurrence before assigning membership.
- Preserve nested LINE nodes and tag their physical descendants with the outer ring name.
- Verify that the source element prototype, a normal LINE using it, and a second RING remain
  independent.
- Verify explicit duplicate members such as `(d1,d1)` become distinct tagged occurrences.

### 5. Make OrbitThreader topology-aware only where required

- First run the end-to-end RING case with the existing OrbitThreader. PR #541 already made
  `IndexMap` tolerate non-contiguous re-entry into the same element, so avoid duplicating that fix.
- Pass root topology/name (and, if needed, circumference) explicitly from the selected Beamline
  through `ParallelTracker`/`OpalBeamline`; do not rediscover it from element names.
- If the ring run fails to terminate or threads beyond the intended first pass, add an explicit
  first-turn completion rule based on ring topology and accumulated path length/circumference.
  Do not use the global bounding-box exit criterion for a closed lattice.
- Keep LINE termination byte-for-byte behaviorally unchanged.
- Add an OrbitThreader test with a small closed lattice that reaches the first element again and
  terminates deterministically while retaining first-crossing ranges.

### 6. End-to-end validation

- Run the compact RING input in a clean output directory.
- Compare it with `sandbox/ring-1/orig` using the committed comparison script.
- Inspect the generated lattice image after PR #540's curved-SBEND mesh fix.
- Run parser, membership, placement, IndexMap, and OrbitThreader unit tests.
- Run the relevant single-rank example, then repeat the end-to-end case with two MPI ranks to
  confirm identical geometry/reference-orbit results.
- Review the diff and document any justified tolerance or golden-reference exclusions.

## Acceptance criteria

- `RING` parses and can be selected by `TRACK, LINE=<ring-name>`.
- Every runtime physical element in the selected ring reports topology `RING` and the correct
  canonical ring name.
- Normal LINEs and registered element prototypes retain default LINE membership.
- Forbidden repetition/reflection and mixed/nested RING composition fail clearly.
- OrbitThreader completes deterministically for the ISIS ring and produces a usable IndexMap.
- The compact RING case matches the semantic/numeric outputs in `sandbox/ring-1/orig` within
  justified tolerances.
- Relevant focused tests pass in single-rank and multi-rank configurations where applicable.

## Next step

Add closed-orbit finding before treating turn-to-turn particle stability or geometric closure as a
RING acceptance criterion.

## Completed first-pass results

- Added and registered the `RING` sequence parser while preserving LINE behavior.
- RING rejects reflection, repetition, nested RINGs, and RING membership in a LINE.
- Named nested LINEs are accepted and remain nested; anonymous lists retain LINE's established
  flattening behavior.
- Removing the old `RING` string-constant registration preserves unquoted `TYPE=RING`: the string
  parser treats the registered RING sequence object name as the literal string `RING`.
- Added `sandbox/ring-1/isis_sbend_ring.in`, its reproducible generator, and a pandas comparison
  script. This remains a sandbox end-to-end regression and is intentionally not a unit test.
- One-rank and two-rank OpenMP runs complete successfully. Their 768 ElementPositions records and
  45 DesignPath records match `orig` numerically within `2e-8`; this tolerance covers the final
  printed digit of the text placement output.
- Added host-side `BeamlineMembership` metadata to `ElementBase`. The default is LINE with no
  owner; RING requires a non-empty owner name, and cloning preserves the value.
- RING roots carry their canonical name. Immediately before tracking, every non-synthetic member
  occurrence is recursively cloned and tagged, leaving registered prototypes and ordinary LINEs
  unchanged. This timing preserves element geometry populated during the global update.
- Corrected LINE member construction to share the registered element's existing `shared_ptr`
  control block; this prevents occurrence replacement from deleting a registered prototype.
- Added in-memory ownership tests for validation, clone preservation, duplicate occurrences,
  recursive nested-LINE tagging, and independence between two rings. ISIS remains intentionally
  outside the unit-test suite.
- Final clean one-rank and two-rank OpenMP runs complete successfully. Each produces 768
  ElementPositions records and 45 DesignPath records matching `orig` within `2e-8`.
- Focused tests pass: TestRing (6), TestOrbitThreader (2), TestIndexMap (4), TestSolenoid (10).
- Explicit `RUN,TURNS=N` on a RING now derives the path target from `N` times the circumference and
  supplies a conservative integration-step budget. The default/omitted TURNS path and LINE
  schedules are unchanged; the first implementation accepts one DT/MAXSTEPS/ZSTOP segment.
- `isis_sbend_ring_10t.in` completed ten nominal circumferences (1633.630220 m) with one particle.
  The reference trajectory is not expected to close or remain stable until a closed-orbit finder
  is available. `plot_ring_10t.py` provides a reproducible visualization of this result.
- RING OrbitThreader construction now stops after exactly one circumference, independently of the
  particle-tracking `MAXSTEPS` budget. The resulting IndexMap is reused for all later turns by
  folding absolute path-length queries into that one-turn interval, including queries that cross
  the turn seam. LINE map construction and lookup remain non-periodic.
- IndexMap exposes the signed turn quotient and a topological phase in `[0, 2*pi)`, both derived
  from accumulated path length and circumference. This is intentionally not a geometric azimuth:
  it remains well-defined for the non-closing reference trajectory used before closed-orbit finding.
- The one-rank `TURNS=10` example again reaches 1633.630220 m with one particle while its
  DesignPath output contains only the first 163.36 m turn. A two-rank one-turn run produces the
  same 768 placements and matches the 45-row explicit golden trajectory prefix within `2e-8`.
- Focused tests pass: TestRing, TestOrbitThreader, TestIndexMap, and TestSolenoid.
