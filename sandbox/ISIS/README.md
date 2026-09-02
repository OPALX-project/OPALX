# ISIS SBEND placement survey

`isis_sbend_survey.in` is a geometry-only replacement for the legacy
`isisring_with_elemedge.in` deck. It keeps all 192 source element names and
lengths, replaces the ten curved `MULTIPOLET` dipoles with native `SBEND`
elements, and assigns explicit global placement coordinates in the style of
`sandbox/test-sbend-ring-1`.

All non-bending devices are represented as `DRIFT` envelopes. Consequently,
this input checks geometry and ring closure only; it is not a replacement for
the ISIS optics, RF, fringe fields, or tracking model.

The native OPALX `SBEND` uses `L` as reference-arc length. The old curved
`MULTIPOLET` uses the same convention (`rho = L / ANGLE`), so the ten dipoles
retain `L = 4.4 m`. Their source angle `0.628319 rad` is the rounded value of
36 degrees; the survey uses exactly `2*pi/10` so ten sectors close exactly.
For explicit placement, `SBEND.X/Y/Z` is the reference-arc entrance and
`SBEND.THETA` is the entrance-tangent orientation. This is verified against
the generated OPALX element-position export, including every interface.

Regenerate and run from this directory with:

```sh
/Users/adelmann/.venv-h6/bin/python generate_isis_sbend_survey.py
/Users/adelmann/git/opalx/omp-build/src/opalx isis_sbend_survey.in
/Users/adelmann/.venv-h6/bin/python -m py_compile data/isis_sbend_survey_ElementPositions.py
/Users/adelmann/.venv-h6/bin/python check_isis_sbend_survey.py
```

The placement export is
`data/isis_sbend_survey_ElementPositions.txt`. Compare `BEGIN: MSTART` and
`BEGIN: MCLOSE` for the numerical closure check. The checker writes a
per-element CSV, a Markdown report, and PNG/PDF survey plots.

The element-position Python/VTK exporter now represents each native `SBEND`
as a triangulated curved body following the reference arc. Empty element
meshes are also emitted as valid empty Python lists. OPALX's `1e6 m`
unbounded-aperture sentinel is ignored for mesh sizing and rendering, so
markers without a physical aperture no longer dominate the viewer scale.
