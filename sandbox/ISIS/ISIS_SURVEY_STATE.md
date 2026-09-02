# ISIS SBEND survey task state

- Goal: replace the ten legacy curved `MULTIPOLET` dipoles with native
  `SBEND`s and survey the explicitly placed ISIS ring.
- Source: `isisring_with_elemedge.in` (192 elements, ten bends, 163.36282 m
  reference-path length).
- Design decision: preserve `L=4.4 m` as bend arc length; use the exact angle
  `2*pi/10`; model every zero-angle component as a same-length `DRIFT` for this
  geometry-only deck.
- Generated files: `isis_sbend_survey.in` and
  `isis_sbend_survey_distribution.txt`.
- Generator: `generate_isis_sbend_survey.py` validates source `ELEMEDGE`
  continuity and analytic ring closure before writing the deck.
- Placement convention verified from the OPALX export: `SBEND.X/Y/Z` is the
  arc entrance and `SBEND.THETA` is the entrance tangent.
- Validation:
  - generator analytic closure: `1.986e-14 m`;
  - OPALX one-rank run: passed;
  - OPALX two-rank run: passed;
  - exported MSTART-to-MCLOSE closure: `0.0 m` at printed precision;
  - maximum interface gap: `1.000e-08 m` (text-export rounding);
  - final SBEND exit-to-MCLOSE gap: `5.235e-12 m`;
  - generated element-position Python passes `py_compile`;
  - all ten SBENDs export as curved meshes with 76 vertices and 148 triangles
    each;
  - unbounded-aperture markers have empty meshes instead of `1e6 m` radii;
    the exported model is bounded by 52.75 m and the VTK totals are 26,968
    points and 27,688 triangle cells.
- OPALX exporter change: `MeshGenerator` now builds native `SBEND` meshes and
  safely serializes empty triangle arrays. It also rejects the `1e6 m`
  unbounded-aperture sentinel as a physical mesh or drift-size reference.
  Focused `TestSolenoid` and `TestOpalBeamlinePlacement` tests pass.
- Reports: `isis_sbend_survey_report.md`, `isis_sbend_survey_summary.csv`, and
  `data/isis_sbend_survey.{png,pdf}`.
- Remaining scope limitation: straight magnets and RF devices are drift
  geometry proxies; optics and full-turn tracking are intentionally not tested.
