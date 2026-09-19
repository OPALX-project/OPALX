# BeamBeam sandbox attic

This directory preserves superseded studies and publication assets for
provenance. Nothing here defines the current BeamBeam physics model, validation
acceptance criteria, or production workflow. Do not use an attic path as the
default output of an active script.

The authoritative description is
[`../BEAMBEAM_PHYSICS_AND_VALIDATION.md`](../BEAMBEAM_PHYSICS_AND_VALIDATION.md).

## Contents

| path | archived material | reason archived | replacement |
|---|---|---|---|
| `notes/boosted-gaussian-witness/` | LaTeX note, PDF, talk, generated tables, and figures | built around the superseded spherical rest-frame Gaussian and mixed historical validation stages | `../BEAMBEAM_PHYSICS_AND_VALIDATION.md` and `../reduced-order-model/` |
| `notes/track12-model/` | older Track-12 LaTeX derivation and rendered figures | predates the exact one-clock emitted-particle workflow and later field/gather corrections | `../track12particles/opalx/timed/` and `../track12particles/opalx/README.md` |
| `notes/beam-beam-summary-legacy.md` | former repository-root BeamBeam summary | contains the obsolete isotropic source description and implementation gaps that have since been fixed | `../BEAMBEAM_PHYSICS_AND_VALIDATION.md` |
| `notes/master-integration-2026-08-18.md` | completed master-merge state report | useful merge provenance but no longer an active task state | `../../BEAMBEAM_REDESIGN_STATE.md` |
| `models/spherical-rest-gaussian/` | former `sandbox/analytic-model` calculations and generated plots | isotropic/spherical source is not the agreed production or manufactured model | `../reduced-order-model/python/rigid_two_gaussian_fields.py` |
| `track12/pre-reference-offset-scans/` | grid-free, particle-count, window-size, and c0-field scans plus their dedicated scripts | generated before the witness reference-coordinate and MPI owner-rank fixes; not acceptance evidence for the current implementation | `../track12particles/opalx/timed/scan_first_kick_fields.py` and the timed A100 result directories |
| `legacy-python/oldpy/` | older one-off visualization scripts | superseded by maintained scripts under `../python/` | `../python/` |

## Interpretation rule

Archived outputs may be useful for debugging regressions or reconstructing the
development history, but their numbers must be labeled with the historical
model and code state. In particular:

- the spherical Gaussian must not be presented as the current manufactured
  source;
- pre-reference-offset Track-12 kicks must not be compared as current OPALX
  results;
- old window and particle scans must not be used to set production mesh or
  aperture parameters;
- old CAIN comparisons do not supersede the exact-timed 13,012-sample analysis.

Large tracked result directories were moved rather than regenerated so Git can
retain their provenance without introducing new binary data.
