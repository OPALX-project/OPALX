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

## Remaining validation

1. Run the full-witness staged input on Merlin6/A100.
2. Confirm per-step transverse mesh bounds contain both witness species in the source frame.
3. Compare witness trajectories and final momenta with the previous validated OPALX result,
   manufactured solution, and CAIN comparison.
4. Check one-rank and multi-rank agreement.
5. Update `sandbox/note/bb-note.tex` only after these results are accepted.
