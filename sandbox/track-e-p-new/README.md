# BeamBeam fixed-z, dynamic-transverse validation

This directory isolates inputs and results for the redesigned BeamBeam field domain from the
established `sandbox/track-e-p` baseline.

The model under test is:

- physical BeamBeam element: 32 cm long with a 15 cm radius aperture;
- collective-field window: fixed at IP +/- 10 mm;
- transverse field domain: recomputed every active step from the physical primary and all emitted
  witness particles, expressed in the source-field frame;
- charge deposition: physical primary only, with the mirrored-primary field constructed by the
  BeamBeam solver;
- completion: after the complete physical primary crosses the downstream field-window boundary;
- post-completion tracking: all particle containers remain active, but no further self-field solve
  or BeamBeam witness gather is performed.

`gamma_gamma_pairs-20mm-dynamic-xy-staged.in` is the first production-population validation deck.
It reuses the immutable CAIN-derived witness files from `../track-e-p` and deliberately contains no
`RETIRE_TIME` setting.

`gamma_gamma_pairs-20mm-dynamic-xy-smoke.in` uses only 32 primary macroparticles and a
coarse/fine/coarse schedule that crosses both field-window boundaries. Use it for local lifecycle
and parser checks before submitting the production-population deck. Witness-domain coverage is
kept in the focused unit test because constructing the full witness reference orbits dominates a
small local smoke run.

Run artifacts and conclusions belong here until the physics and trajectory comparisons pass. Only
then should `sandbox/note/bb-note.tex` be updated.
