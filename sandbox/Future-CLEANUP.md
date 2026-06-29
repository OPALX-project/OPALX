# Future Cleanup Notes

## Tracker-side query policy

Two recent tracker call sites exposed the same structural issue: the tracker is compensating for query-model gaps with local special cases.

### 1. `ParallelTracker::computeExternalFields(...)`

Previously the tracker:

- queried the traced occupancy model via `OrbitThreader::query(...)`, and then
- manually rescanned the action-range registration segments to re-add `SBEND` / `RBEND`.

That worked as a compatibility patch, but the layering was poor:

- `ParallelTracker` should not know how to reconstruct bend action-range semantics from raw segments.
- The distinction between traced occupancy and registered action range should be expressed by a named query API.

Current direction:

- `OrbitThreader::query(...)` remains the traced-occupancy query.
- `OrbitThreader::queryActionRangeElements(...)` now exposes the registration-model query explicitly.

Later cleanup:

- revisit whether the tracker should union these two models at all, or whether there should be a single higher-level "active-for-bunch-field-application" query owned outside the tracker.
- check whether the action-range model should carry richer semantics than a flat interval union.

### 2. `ParallelTracker::updateReferenceParticles(...)`

Previously the tracker:

- queried `OpalBeamline::getElements(refPosition)`, then
- rescanned all beamline elements to re-add online `MONITOR` elements.

Reason for the old behavior:

- `getElements(position)` is a geometric containment query.
- online monitors are diagnostics with thin or zero body length and do not fit cleanly into "ordinary geometry only" activation logic.

This means the re-add was compensating for a mismatch between:

- geometric containment, and
- "acts on the reference particle now" semantics.

Current direction:

- `OpalBeamline::getReferenceElements(...)` now hides that policy from `ParallelTracker`.

Later cleanup:

- revisit monitor activation in a more general framework, not as a monitor-only exception.
- define a general query vocabulary for:
  - geometric containment,
  - traced occupancy,
  - action-range overlap,
  - reference-particle-active elements,
  - bunch-field-active elements.
- decide whether monitors, bends, probes, and similar thin/support-driven elements should all participate through one common activation model instead of element-type-specific exceptions.

## Design question to revisit

The real issue is not "why re-add monitors?" or "why re-add bends?" in isolation.

The real issue is that different subsystems currently need different notions of "active element":

- geometry-active,
- reference-path-active,
- registered-action-range-active,
- reference-particle-active,
- bunch-application-active.

Those should probably become explicit first-class concepts with named query interfaces and tests, instead of being reconstructed ad hoc inside tracker loops.
