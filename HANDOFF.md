# OPALX BeamBeam Handoff

Last updated: 2026-06-10

## Active Task: Refresh BeamBeam Diagnostic Timeline Simulation

Goal:

- Rerun the sandbox simulation workflow that feeds the BeamBeam Diagnostic
  Timeline in `sandbox/regression/sandbox_regression_overview.{tex,pdf}`.
- Regenerate `sandbox/regression/current_metrics.csv`,
  `sandbox/regression/sandbox_regression_overview.tex`, and
  `sandbox/regression/sandbox_regression_overview.pdf`.

Planned commands from the overview How To Reproduce section:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/run_sandbox_regressions.py \
  --run-opalx \
  --opalx-exe build_openmp/src/opalx
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py
```

Assumptions:

- Use the existing `build_openmp/src/opalx` executable.
- Preserve the existing uncommitted `sandbox/track-e-p/e-e+.gpl` deletion and
  untracked `sandbox/track-e-p/attic/` archive.

Status:

- Completed by Codex on 2026-06-10.
- The simulation harness completed all three OPALX cases and wrote
  `sandbox/regression/current_metrics.csv`, but exited with regression
  mismatches against the accepted baseline.
- The mismatches are in `pairs.track_stat.gamma_gamma_pairs-2.*`.
  Container `c0` now ends retired with zero particles/energy, while witness
  containers `c1` and `c2` remain populated with small final metric shifts.
- The overview was regenerated successfully:
  - `sandbox/regression/opalx_impact_drift_comparison.png`
  - `sandbox/regression/sandbox_regression_overview.tex`
  - `sandbox/regression/sandbox_regression_overview.pdf`
- Ghostscript rendered page 4 successfully to
  `/private/tmp/sandbox_regression_overview_page4.png`; visual inspection
  showed the BeamBeam Diagnostic Timeline page is coherent.

Verification commands run:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/run_sandbox_regressions.py \
  --run-opalx \
  --opalx-exe build_openmp/src/opalx
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py
git diff --check -- HANDOFF.md sandbox/regression/current_metrics.csv sandbox/regression/sandbox_regression_overview.tex
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_overview_pdf.py sandbox/regression/run_sandbox_regressions.py
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=4 -dLastPage=4 \
  -sOutputFile=/private/tmp/sandbox_regression_overview_page4.png \
  sandbox/regression/sandbox_regression_overview.pdf
```

## Active Task: Cylinder BeamBeam Aperture And Figure 4

Goal:

- Change `sandbox/track-e-p/gamma_gamma_pairs-2.in` from a rectangular BeamBeam
  aperture to a cylindrical transverse aperture.
- Rerun the simulation and confirm the Figure 3 BeamBeam diagnostic timeline
  still appears.
- Track beyond 121 ps and inspect whether witness particles are deleted/lost as
  they reach the cylinder boundary.
- Add Figure 4 to the regression overview from the available loss/count data,
  matching the uploaded sketch: lost `e-`/`e+` counts versus `z` above a
  cylinder schematic.

Implementation decisions:

- OPALX aperture syntax uses `CIRCLE(diameter)` for a cylindrical transverse
  aperture. The input now uses `APERTURE="CIRCLE(0.02)"`, preserving the
  previous 1 cm half-aperture from `RECTANGLE(0.02, 0.02, 0.05)`.
- `MAXSTEPS` was increased from 130 to 220 so the run extends well beyond the
  `RETIRE_TIME = 121e-12`.
- The existing output path provides per-container `.stat` files with
  `numParticles` and `partsOutside`; no BeamBeam aperture-specific `.loss` file
  was found for this element path. Figure 4 is therefore generated from
  positive step-to-step drops in `numParticles` for containers `c1` and `c2`.

Changed files so far:

- `sandbox/track-e-p/gamma_gamma_pairs-2.in`
- `sandbox/regression/make_overview_pdf.py`
- `HANDOFF.md`

Next step:

- Completed. `gamma_gamma_pairs-2.in` ran successfully with
  `APERTURE="CIRCLE(0.02)"` and `MAXSTEPS = 220`.
- Figure 3 behavior is preserved. The direct run emitted the expected
  `BB-DIAG` sequence: inactive, active, witness injection, source overlap,
  completed with source retirement pending, then completed with
  `retired_bunches=1` and `source_active=FALSE`.
- The run reached 220 ps. Final witness populations remained
  `c1 = 1297`, `c2 = 1297`, and `partsOutside = 0` in both stat files.
  Therefore, the current OPALX BeamBeam path records witness trajectories but
  does not delete witness particles at the cylindrical BeamBeam aperture
  boundary.
- H5 trajectory inspection showed particles physically crossed the 1 cm
  cylinder radius:
  - `c1`: 1059 particles outside 1 cm at the last H5 step, max radius
    `3.153819082413395e-2 m`.
  - `c2`: 1053 particles outside 1 cm at the last H5 step, max radius
    `3.229996618876009e-2 m`.
- Figure 4 was added as a cylinder-edge crossing histogram from the first H5
  particle position with radius at least 1 cm for each witness particle. It is
  intentionally captioned as crossing data, not deleted-loss data.

Generated/updated artifacts:

- `sandbox/regression/current_metrics.csv`
- `sandbox/regression/gamma_gamma_cylinder_losses.png`
- `sandbox/regression/sandbox_regression_overview.tex`
- `sandbox/regression/sandbox_regression_overview.pdf`

Verification commands run for this task:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_overview_pdf.py
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx gamma_gamma_pairs-2.in
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/run_sandbox_regressions.py
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_overview_pdf.py sandbox/regression/run_sandbox_regressions.py
git diff --check -- HANDOFF.md sandbox/track-e-p/gamma_gamma_pairs-2.in sandbox/regression/make_overview_pdf.py sandbox/regression/sandbox_regression_overview.tex sandbox/regression/current_metrics.csv
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=5 -dLastPage=5 \
  -sOutputFile=/private/tmp/sandbox_regression_overview_page5.png \
  sandbox/regression/sandbox_regression_overview.pdf
```

## Active Task: Large Cylinder BeamBeam Setup And Figures 5-6

Goal:

- Keep `sandbox/track-e-p/gamma_gamma_pairs-2.in` as the proof-of-principle
  cylinder setup.
- Add a new setup with BeamBeam cylinder radius 15 cm and BeamBeam window
  length 32 cm.
- Adapt the witness injection and source retirement timing consistently with
  the new interaction-point position.
- Record the new setup in Figures 5 and 6 of
  `sandbox/regression/sandbox_regression_overview.{tex,pdf}`.

Implementation decisions:

- Added `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in`.
- The new input sets:
  - `bb_length = 0.32`
  - `bb_ip_s = drift_length + 0.5*bb_length = 0.17 m`
  - `APERTURE="CIRCLE(0.30)"`, because OPALX `CIRCLE` takes diameter, so this
    is a 15 cm radius cylinder.
  - `ZSTOP = 2*drift_length + bb_length`
  - `MAXSTEPS = 1600`
- The timing preserves the proof-of-principle offsets relative to the source IP
  crossing:
  - witness pair injected at `primary_ip_time - 16.747e-12`, about `550.3 ps`
  - source reaches IP at about `567.1 ps`
  - source retired at `primary_ip_time + 4.253e-12`, about `571.3 ps`

Status:

- Completed by Codex on 2026-06-10.
- The large-cylinder input ran successfully with
  `/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx`.
- The run emitted the expected `BB-DIAG` state sequence: inactive, active,
  witness injection, source overlap, completed with retirement pending, and
  completed with `retired_bunches=1 source_active=FALSE`.
- H5 trajectory inspection showed first crossings of the 15 cm radius:
  - `c1`: 701 first crossings; max final radius `0.27683411836950467 m`
  - `c2`: 692 first crossings; max final radius `0.28372368226326133 m`
- As in the proof-of-principle setup, OPALX records witness trajectories but
  does not delete witness particles at the BeamBeam cylinder boundary in this
  path. Figure 6 is therefore a first-crossing histogram, not a deleted-loss
  histogram.
- The generated large H5/stat outputs are intentionally ignored by `.gitignore`
  (`*.h5`, `*.stat`). The generated PNGs are also ignored (`*.png`), while the
  regenerated PDF embeds the plots.

Generated/updated artifacts:

- `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in`
- `sandbox/regression/gamma_gamma_large_cylinder_crossings.png`
- `sandbox/regression/sandbox_regression_overview.tex`
- `sandbox/regression/sandbox_regression_overview.pdf`

Verification commands run for this task:

```sh
cd sandbox/track-e-p
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx gamma_gamma_pairs-large-cylinder.in
cd ../..
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_overview_pdf.py sandbox/regression/run_sandbox_regressions.py
git diff --check -- HANDOFF.md sandbox/track-e-p/gamma_gamma_pairs-2.in sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in sandbox/regression/make_overview_pdf.py sandbox/regression/sandbox_regression_overview.tex sandbox/regression/current_metrics.csv
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=6 -dLastPage=6 \
  -sOutputFile=/private/tmp/sandbox_regression_overview_page6.png \
  sandbox/regression/sandbox_regression_overview.pdf
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=7 -dLastPage=7 \
  -sOutputFile=/private/tmp/sandbox_regression_overview_page7.png \
  sandbox/regression/sandbox_regression_overview.pdf
```

## Active Task: Reduced Primary Charge Large-Cylinder Comparison

Goal:

- Rerun the 15 cm radius, 32 cm length large-cylinder setup with the primary
  beam charge reduced by a factor of `100000`.
- Compare the nominal and reduced-charge cylinder-edge crossing histograms.

Implementation decisions:

- Added `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in`.
- The new input keeps the nominal large-cylinder geometry, timing, field solver,
  witness pair input, and tracking length unchanged.
- The only intended physics change is:
  - `primary_charge_scale = 1.0e-5`
  - `primary_electrons_per_bunch = 1.25e10 * primary_charge_scale`
- Added `sandbox/regression/gamma_gamma_large_cylinder_charge_compare.png`
  generation in `sandbox/regression/make_overview_pdf.py`.
- Added Figure 7 to the overview PDF. It overlays nominal and reduced-charge
  first-crossing histograms for `e-` and `e+` and shows the reduced-minus-
  nominal bin difference.

Status:

- Completed by Codex on 2026-06-10.
- The reduced-charge OPALX run completed successfully and emitted the expected
  `BB-DIAG` sequence through witness injection, overlap, completion, and source
  retirement.
- Histogram comparison result:
  - Nominal large-cylinder run: `c1 = 701`, `c2 = 692` first crossings.
  - Reduced primary charge run: `c1 = 701`, `c2 = 692` first crossings.
  - With the current 18-bin histogram, the reduced-minus-nominal bin-count
    difference is zero for both species.
  - Matching particle IDs shift in first-crossing `z` by at most about
    `141 um` for `c1` and `142 um` for `c2`.
- Figure 7 is a first-crossing comparison, not a deleted-loss comparison,
  because this BeamBeam path still does not delete witness particles at the
  cylinder boundary.

Generated/updated artifacts:

- `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in`
- `sandbox/regression/gamma_gamma_large_cylinder_charge_compare.png`
- `sandbox/regression/sandbox_regression_overview.tex`
- `sandbox/regression/sandbox_regression_overview.pdf`

Verification commands run for this task:

```sh
cd sandbox/track-e-p
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx gamma_gamma_pairs-large-cylinder-q1em5.in
cd ../..
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_overview_pdf.py sandbox/regression/run_sandbox_regressions.py
git diff --check -- HANDOFF.md sandbox/track-e-p/gamma_gamma_pairs-2.in sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in sandbox/regression/make_overview_pdf.py sandbox/regression/sandbox_regression_overview.tex sandbox/regression/current_metrics.csv
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=7 -dLastPage=7 \
  -sOutputFile=/private/tmp/sandbox_regression_overview_page7_qcompare_clean.png \
  sandbox/regression/sandbox_regression_overview.pdf
```

## Active Task: Large-Cylinder Field/Histogram Movie

Goal:

- Build an animation based on Figure 6 that shows:
  - the 15 cm radius, 32 cm long BeamBeam cylinder,
  - the electric field in the cylinder,
  - the running simulation time,
  - actual OPALX `e-`/`e+` witness particles,
  - the cumulative first-crossing histogram.

Implementation decisions:

- Added `sandbox/regression/make_large_cylinder_field_movie.py`.
- The movie uses actual OPALX H5 particle trajectories from
  `gamma_gamma_pairs-large-cylinder_c1.h5` and `_c2.h5`.
- The field overlay is diagnostic/analytic, not an OPALX field dump:
  - It is a vectorized version of the spherical boosted Gaussian source-field
    model from `sandbox/python/boosted_gaussian_witness.py`.
  - It uses the nominal primary charge `1.25e10 e`, primary kinetic energy
    `245 MeV`, and a spherical source width `sigma = 0.6 mm`.
  - The overlay is set to zero after `RETIRE_TIME = 571.35 ps`, matching the
    input's source retirement behavior.
- The lower panel is a cumulative version of the Figure 6 histogram. It records
  first radial crossings of 15 cm as a function of crossing `z`.

Status:

- Completed by Codex on 2026-06-10.
- The MP4 was generated:
  - `sandbox/data/gamma_gamma_large_cylinder_field_histogram.mp4`
- A preview frame was generated:
  - `sandbox/data/gamma_gamma_large_cylinder_field_histogram_preview.png`
- An active-field verification frame was extracted to:
  - `/private/tmp/gamma_gamma_field_movie_active_frame.png`
- Key physics observation from the movie/data:
  - The analytic primary field is visible during the source-active window near
    `t = 567 ps`.
  - The source is retired at about `571.4 ps`, so the overlay is zero afterward.
  - The first 15 cm crossings occur much later, starting around `1106 ps` for
    `c2` and `1120 ps` for `c1`, with medians around `1322-1329 ps`.
  - This timing supports the hypothesis that the Coulomb interaction seen by
    the low-energy pair is either very weak or effectively absent after source
    retirement in the current BeamBeam path.

Verification commands run for this task:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_large_cylinder_field_movie.py
git diff --check -- sandbox/regression/make_large_cylinder_field_movie.py
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_large_cylinder_field_movie.py
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -c "import imageio.v2 as imageio; from pathlib import Path; movie=Path('sandbox/data/gamma_gamma_large_cylinder_field_histogram.mp4'); reader=imageio.get_reader(movie); frame=reader.get_data(8); imageio.imwrite('/private/tmp/gamma_gamma_field_movie_active_frame.png', frame); reader.close(); print(frame.shape)"
```

## Active Task: Large-Cylinder Run With 1000 ps Primary Retirement

Goal:

- Redo the large-cylinder OPALX simulation with primary source retirement at
  `1000 ps`.
- Regenerate the field/histogram movie from the new OPALX H5 trajectories.

Implementation decisions:

- Added `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps.in`.
- The setup preserves the previous large-cylinder geometry, source/witness
  charge, witness injection time, field solver, and tracking length.
- The only intended BeamBeam timing change is:
  - `primary_retire_time = 1000e-12`
- The movie was regenerated with:
  - `--stem gamma_gamma_pairs-large-cylinder-retire1000ps`
  - `--retire-time-s 1000e-12`
  - a short display label to avoid title clipping.

Status:

- Completed by Codex on 2026-06-10.
- The OPALX run completed successfully.
- `BB-DIAG` showed:
  - inactive,
  - active before witness injection,
  - witness injection,
  - source/witness overlap,
  - overlap ended while source remained active,
  - completed with source retirement pending,
  - completed with `retired_bunches=1 source_active=FALSE`.
- Corrected first-crossing counts at 15 cm:
  - Previous large-cylinder run: `c1 = 701`, `c2 = 692`.
  - `1000 ps` retirement run: `c1 = 701`, `c2 = 691`.
- First-crossing time ranges in the `1000 ps` retirement run:
  - `c1`: min/median/max `1120 / 1329 / 1599 ps`
  - `c2`: min/median/max `1106 / 1322 / 1600 ps`
- New generated movie:
  - `sandbox/data/gamma_gamma_large_cylinder_retire1000ps_field_histogram.mp4`
- New generated preview:
  - `sandbox/data/gamma_gamma_large_cylinder_retire1000ps_field_histogram_preview.png`

Verification commands run for this task:

```sh
cd sandbox/track-e-p
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx gamma_gamma_pairs-large-cylinder-retire1000ps.in
cd ../..
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_large_cylinder_field_movie.py \
  --stem gamma_gamma_pairs-large-cylinder-retire1000ps \
  --label "large cylinder, retire 1000 ps" \
  --retire-time-s 1000e-12 \
  --output sandbox/data/gamma_gamma_large_cylinder_retire1000ps_field_histogram.mp4 \
  --preview sandbox/data/gamma_gamma_large_cylinder_retire1000ps_field_histogram_preview.png \
  --preview-time-ps 1120
```

## Active Task: 1000 ps Retirement With 1e-5 Primary Charge Comparison

Goal:

- Run the 1000 ps retirement large-cylinder setup again with the primary beam
  charge reduced by `1e-5`.
- Make a comparison histogram against the nominal 1000 ps retirement run.

Implementation decisions:

- Added `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5.in`.
- The input preserves the 15 cm radius, 32 cm BeamBeam length, 1000 ps source
  retirement, witness injection timing, pair charge, and tracking length.
- The only intended charge change is:
  - `primary_charge_scale = 1.0e-5`
  - `primary_electrons_per_bunch = 1.25e10 * primary_charge_scale`
- Added reusable comparison script:
  - `sandbox/regression/compare_cylinder_crossing_histograms.py`
- Generated comparison histogram:
  - `sandbox/regression/gamma_gamma_large_cylinder_retire1000ps_q1em5_compare.png`

Status:

- Completed by Codex on 2026-06-10.
- The reduced-charge 1000 ps OPALX run completed successfully and emitted the
  expected delayed-retirement `BB-DIAG` sequence.
- Histogram comparison against nominal 1000 ps retirement:
  - `e-`: nominal/reduced `701/701`, absolute bin difference sum `2`,
    max matched-particle first-crossing shift `262.2 um`.
  - `e+`: nominal/reduced `691/692`, absolute bin difference sum `1`,
    max matched-particle first-crossing shift `243.6 um`.
- The reduced-charge histogram is nearly coincident with the nominal 1000 ps
  retirement histogram.

Verification commands run for this task:

```sh
cd sandbox/track-e-p
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5.in
cd ../..
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/compare_cylinder_crossing_histograms.py
git diff --check -- sandbox/regression/compare_cylinder_crossing_histograms.py sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5.in
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/compare_cylinder_crossing_histograms.py
```

## Active Task: Manufactured-Solution Note Updated For Large BeamBeam Window

Goal:

- Update the manufactured boosted-Gaussian witness note and Python generator to
  account for the current large BeamBeam diagnostic geometry:
  - cylinder radius `15 cm`
  - BeamBeam length `32 cm`
  - IP at `s = 0.17 m`
- Regenerate the note tables, figures, and PDF with the new data.

Implementation decisions:

- Updated `sandbox/python/boosted_gaussian_witness.py` defaults:
  - `--sigma-m` default changed to `0.6e-3`
  - `--pair-duration-ps` default changed to `1050`
  - `--pair-dt-ps` default changed to `0.1`
  - added `--beambeam-radius-m` default `0.15`
  - added `--beambeam-length-m` default `0.32`
- Updated the pair-path plot to mark the 15 cm radial aperture scale and the
  32 cm BeamBeam longitudinal window.
- Updated `sandbox/note/boosted_gaussian_witness.tex` to document the large
  BeamBeam geometry, timing relation to the OPALX large-cylinder input, and new
  manufactured propagation data.

Status:

- Completed by Codex on 2026-06-10.
- Regenerated:
  - `sandbox/note/boosted_gaussian_witness_initial_cases_table.tex`
  - `sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex`
  - `sandbox/note/figs/boosted_gaussian_witness_physical_field_t0.png`
  - `sandbox/note/figs/boosted_gaussian_witness_physical_paths.png`
  - `sandbox/note/boosted_gaussian_witness.pdf`
- Key regenerated data:
  - one-kick off-axis symmetric field is now `8.682 GV/m` for
    `sigma'=0.6 mm`;
  - manufactured pair propagation uses `T=1050 ps` and reaches about
    `246.9 mm` transverse path length;
  - symmetric collision final field-on minus free displacement is about
    `0.986 um` in `x`;
  - asymmetric collision adds about `0.126 um` longitudinal displacement.
- Rendered PDF pages 7-10 to `/private/tmp/boosted_gaussian_witness_page*.png`
  and visually checked the updated tables and figures.

Verification commands run for this task:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/python/boosted_gaussian_witness.py --latex-table --write-latex-table sandbox/note/boosted_gaussian_witness_initial_cases_table.tex
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/python/boosted_gaussian_witness.py --pair-demo --write-pair-latex-table sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex --output-prefix sandbox/data/boosted_gaussian_witness
cp sandbox/data/boosted_gaussian_witness_field_t0.png sandbox/note/figs/boosted_gaussian_witness_physical_field_t0.png
cp sandbox/data/boosted_gaussian_witness_paths.png sandbox/note/figs/boosted_gaussian_witness_physical_paths.png
latexmk -pdf -interaction=nonstopmode -halt-on-error boosted_gaussian_witness.tex
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/python/boosted_gaussian_witness.py sandbox/python/beam-beam-manufactured-solution.py
git diff --check -- sandbox/python/boosted_gaussian_witness.py sandbox/note/boosted_gaussian_witness.tex sandbox/note/boosted_gaussian_witness_initial_cases_table.tex sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=7 -dLastPage=10 \
  -sOutputFile=/private/tmp/boosted_gaussian_witness_page%02d.png \
  sandbox/note/boosted_gaussian_witness.pdf
```

## Active Task: Manufactured Pair Integration Mode

Goal:

- Update the manufactured analysis mode so the two low-energy particles can be
  integrated directly in the manufactured boosted-Gaussian field from the
  combined analysis front-end.

Implementation decisions:

- Added a `pair` subcommand to `sandbox/python/beambeam_analysis.py`.
- The mode delegates to `sandbox/python/boosted_gaussian_witness.py` instead of
  duplicating the field or Boris-kick implementation.
- The new mode uses the same large-cylinder defaults as the regenerated note:
  - `sigma_m = 0.6e-3`
  - `pair_duration_ps = 1050`
  - `pair_dt_ps = 0.1`
  - `beambeam_radius_m = 0.15`
  - `beambeam_length_m = 0.32`
- Fixed the dynamic module loader in `beambeam_analysis.py` to register loaded
  modules in `sys.modules` before execution. This is required for dataclasses in
  dynamically loaded modules.

Status:

- Completed by Codex on 2026-06-10.
- Verified command:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/python/beambeam_analysis.py pair \
  --output-prefix sandbox/data/boosted_gaussian_witness_pair_mode_check \
  --write-pair-latex-table /private/tmp/boosted_gaussian_witness_pair_mode_check_table.tex
```

- The mode wrote:
  - `sandbox/data/boosted_gaussian_witness_pair_mode_check_trajectories.csv`
  - `sandbox/data/boosted_gaussian_witness_pair_mode_check_field_t0.png`
  - `sandbox/data/boosted_gaussian_witness_pair_mode_check_paths.png`
  - `/private/tmp/boosted_gaussian_witness_pair_mode_check_table.tex`
- The numerical output matches the regenerated note data:
  - free path about `246.9 mm` over `1050 ps`
  - symmetric field-on minus free displacement about `0.986 um`
  - asymmetric longitudinal displacement about `0.126 um`

## Repository State

- Repository: `/Users/adelmann/git/opalx-beambeam`
- Current branch: `271-implement-interation-point-element`
- Branch is in sync with `origin/271-implement-interation-point-element`.
- Current HEAD: `cb9c3186b Document BeamBeam handoff state and timeline updates`
- The BeamBeam retirement/timeline work, generated FROMFILE inputs, `gamma_gamma_pairs-3.in`, and this handoff file were committed and pushed in `cb9c3186b`.
- The tree still has local changes. Do not assume a clean checkout and do not revert unrelated files.

Current `git status --short` after this handoff refresh should include:

```text
 M HANDOFF.md
 M sandbox/regression/current_metrics.csv
 M sandbox/regression/make_large_cylinder_field_movie.py
 M sandbox/regression/make_overview_pdf.py
 M sandbox/regression/sandbox_regression_overview.pdf
 M sandbox/regression/sandbox_regression_overview.tex
 D sandbox/track-e-p/e-e+.gpl
 M sandbox/track-e-p/gamma_gamma_pairs-2.in
?? sandbox/track-e-p/attic/
 ?? sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in
 ?? sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in
 ?? sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps.in
 ?? sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5.in
```

Local uncommitted/untracked items:

- `sandbox/track-e-p/e-e+.gpl` is deleted in the worktree.
- `sandbox/track-e-p/attic/` is untracked and contains 273 archived sandbox files:
  - 265 files under `sandbox/track-e-p/attic/data/`
  - archived `gamma_gamma_pairs.in`
  - archived `e-e+.gpl`
  - archived `gamma_gamma_pairs_c{0,1,2}.h5`
  - archived `gamma_gamma_pairs_c{0,1,2}.stat`

Tracked additions now in `cb9c3186b` include:

- `HANDOFF.md`
- `sandbox/track-e-p/gamma_gamma_electrons-t.fromfile`
- `sandbox/track-e-p/gamma_gamma_positrons-t.fromfile`
- `sandbox/track-e-p/gamma_gamma_pairs-3.in`

## Important Working Rules

- Preserve user changes. Do not use `git reset --hard` or `git checkout --` unless explicitly requested.
- Use `apply_patch` for manual edits.
- For Python in this repo, use `/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python`.
- For generated regression overview content, edit `sandbox/regression/make_overview_pdf.py`; then regenerate `sandbox/regression/sandbox_regression_overview.tex` and `.pdf`.
- For C++ physics-facing changes, update tests and document the expected numerical/physics behavior.

## Implemented BeamBeam Changes

### BeamBeam source retirement

The primary BeamBeam source can now be deleted by time using `RETIRE_TIME` on the `BEAMBEAM` element.

Key files:

- `src/Elements/OpalBeamBeam.h`
- `src/Elements/OpalBeamBeam.cpp`
- `src/AbsBeamline/BeamBeamDefinitions.h`
- `src/Algorithms/ParallelTracker.h`
- `src/Algorithms/ParallelTracker.cpp`
- `unit_tests/PartBunch/TestBeamBeam.cpp`

Current input usage in `sandbox/track-e-p/gamma_gamma_pairs-2.in`:

```opal
REAL primary_retire_time = 121e-12;

IP1: BEAMBEAM, L = bb_length,
    ELEMEDGE = drift_length,
    VISUALIZE = FALSE,
    COPY = TRUE,
    WITNESS_CONTAINERS = "1,2",
    RETIRE_TIME = primary_retire_time,
    APERTURE="RECTANGLE(0.02, 0.02, 0.05)";
```

Behavior:

- `RETIRE_TIME` is specified in seconds.
- `RETIRE_TIME = 0` disables retirement.
- Negative `RETIRE_TIME` is rejected in `OpalBeamBeam::update`.
- When `RETIRE_TIME` is reached, the primary source container is retired.
- If retirement fires while BeamBeam is active, `ParallelTracker` first leaves the BeamBeam window, then retires the source container.
- Retirement marks all source particles invalid, deletes invalid particles, and sets container 0 inactive for BeamBeam source use.
- The source container object still exists, but its particle population is removed and it no longer contributes charge/fields.
- Witness containers remain active and continue integrating.

### Witness space-charge switch

The current working input uses:

```opal
WITNESS_CONTAINERS = "1,2"
```

The intended "off" form is:

```opal
WITNESS_CONTAINERS = "NONE"
```

The witness kick implementation is still a source-frame field sampling approximation. It rotates/translates the geometry but does not perform a full Lorentz transform of the source fields into the lab/reference frame seen by low-energy witnesses. This remains an important physics caveat.

### Grepable diagnostics

BeamBeam state diagnostics use the key:

```text
BB-DIAG
```

Current diagnostics include BeamBeam state/counts plus changed booleans such as:

- `source_active`
- `source_retirement_pending`
- `source_bunches_overlap`
- `retired_bunches`
- witness active/inactive state

The earlier noisy diagnostics `stage=after_sc step=...` and `stage=after_emission step=...` were removed.

Expected transition pattern from the current gamma-gamma sandbox run includes:

```text
BB-DIAG ... source_bunches_overlap=TRUE
BB-DIAG BB-state=Completed ... source_retirement_pending=TRUE source_bunches_overlap=FALSE
BB-DIAG BB-state=Completed active_bunches=2 retired_bunches=1 ... source_active=FALSE source_retirement_pending=FALSE
```

## Regression Overview / Timeline Figure

The BeamBeam diagnostic timeline is generated by:

```text
sandbox/regression/make_overview_pdf.py
```

Generated files:

```text
sandbox/regression/sandbox_regression_overview.tex
sandbox/regression/sandbox_regression_overview.pdf
```

Recent timeline changes:

- Added `time [ps]` at the top of the time arrow.
- Added a `116 ps` frame between `111 ps` and `121 ps`.
- In the `116 ps` frame, the source bunches have passed the IP:
  - green bunch on the left of the IP
  - black bunch on the right of the IP
- Added the primary deletion state at `121 ps`:
  - retired primary source shown as dashed grey crossed-out bunches
  - witness pair remains active
- The lower four frames are event-spaced, not to scale.
- The lower four frame spacings were increased by 3x in drawing coordinates.
- Because of the larger spacing, the overview PDF is currently 5 pages; "How To Reproduce" moved to page 5.

To regenerate:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py --skip-impact-plot
```

To render timeline pages for inspection:

```sh
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 \
  -dFirstPage=4 -dLastPage=4 \
  -sOutputFile=/private/tmp/sandbox_regression_overview_page4.png \
  sandbox/regression/sandbox_regression_overview.pdf
```

At handoff time, the preview image was:

```text
/private/tmp/sandbox_regression_overview_page4.png
```

## Older Static Gaussian BeamBeam Validation / Manual Context

This older chat also established and documented a solver-frame validation for
the static electrostatic Gaussian BeamBeam pair. It is not part of the current
BeamBeam retirement/timeline commit, but it is useful context for future
physics/manual work.

Analytic case:

- Solver-frame electrostatic comparison, not a lab-frame Lorentz-boosted field
  comparison.
- Two identical spherical Gaussian bunches, mirrored about the local IP.
- Reference parameters:
  - `sigma = 1 mm`
  - nominal `d = 10 mm`
  - `|Q| = 1.112650055e-14 C`
- Exact field at the mathematical IP is zero by symmetry.
- The scalar potential at the IP is nonzero because the two scalar potential
  contributions add while their gradients cancel.

Local OPALX artifacts still present in this checkout:

```text
sandbox/BeamBeam-static-1V.in
data/sandbox/BeamBeam-static-1V_on_axis_compare.csv
data/sandbox/BeamBeam-static-1V_Ez_on_axis_compare.png
```

The active OPALX diagnostic snapshot used in the write-up had actual
solver-frame separation:

```text
d_OPALX = 9.89313053174 mm
```

Key comparison values from that run:

```text
rho relL2 = 5.221859e-02
phi relL2 = 3.042248e-03
Ex  relL2 = 1.503024e-02
Ey  relL2 = 1.582187e-02
Ez  relL2 = 2.427369e-02
Ez nearest on-axis grid sample to IP:
  analytic = 1.038298 V/m
  OPALX    = 1.072345 V/m
Interpolated Ez at mathematical IP:
  analytic ~= -5.88e-08 V/m
  OPALX    ~=  3.77e-06 V/m
```

Physics manual work was written outside this repository in:

```text
/Users/adelmann/git/physics-manual-opalx/sections/beam-beam/index.qmd
/Users/adelmann/git/physics-manual-opalx/sections/beam-beam/figures/beambeam-static-1v-ez-on-axis.png
/Users/adelmann/git/physics-manual-opalx/sections/beam-beam/reproducers/
```

The manual reproducer bundle contains:

```text
BeamBeam-static-1V.in
beam-beam-manufactured-solution.py
README.md
```

The rendered local preview used in that chat was:

```text
file:///tmp/beam-beam-render/sections/beam-beam/index.html
```

## Verification History

Most recent documentation/timeline checks:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile sandbox/regression/make_overview_pdf.py
git diff --check -- sandbox/regression/make_overview_pdf.py sandbox/regression/sandbox_regression_overview.tex
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/regression/make_overview_pdf.py --skip-impact-plot
```

All passed. LaTeX emits only pre-existing table-width underfull/overfull warnings.

Earlier C++ verification after implementing `RETIRE_TIME`:

```sh
cmake --build build_openmp --target opalx_exe TestBeamBeam -j 8
ctest --test-dir build_openmp -R '^TestBeamBeam$' --output-on-failure
```

The direct and one-rank MPI gamma-gamma runs were also checked after the retirement implementation:

```sh
cd sandbox/track-e-p
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx gamma_gamma_pairs-2.in

/opt/homebrew/Cellar/open-mpi/5.0.8/bin/mpiexec \
  --map-by slot --oversubscribe -np 1 \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx \
  gamma_gamma_pairs-2.in
```

No `Error`, `Exception`, `Abort`, `MPI_ABORT`, `nan`, `Segmentation`, or `Quaternion` failure was seen in those passing runs.

Known caveat from previous investigation:

- A broader two-rank `TestBeamBeam`/PartBunch path had shown a hang before this handoff.
- A two-rank full gamma input had aborted before writing steps in earlier testing.
- These were not counted as passing.
- If multi-rank behavior matters next, rerun from a clean terminal and keep the exact command/output.

## Useful Search Anchors

```sh
rg -n "RETIRE_TIME|BB-DIAG|source_bunches_overlap|retired_bunches|source_active" src sandbox unit_tests
rg -n "collision_timeline_figure|source bunches passed IP|not to scale" sandbox/regression/make_overview_pdf.py
```

## Likely Next Steps

1. Decide whether the 5-page overview PDF is acceptable after the 3x timeline spacing. If not, keep the lower-frame spacing and move/compress surrounding content rather than reducing the requested spacing.
2. Rebuild and rerun the BeamBeam tests if C++ files are touched again.
3. Decide what to do with the local `sandbox/track-e-p/e-e+.gpl` deletion and untracked `sandbox/track-e-p/attic/` archive.
4. If preparing another commit, inspect `git diff` carefully because the remaining local changes are not part of the pushed BeamBeam retirement/timeline commit.

## 2026-06-10 Manufactured Solution Timestep Note Update

User asked to update the manufactured-solution PDF with the frame/magnetic-field/timestep information from the pair-in-field check.

Updated:

```text
sandbox/note/boosted_gaussian_witness.tex
sandbox/note/boosted_gaussian_witness.pdf
```

Content added after the pair kinematics table:

- The manufactured pair calculation is in the lab/reference frame.
- Each source is evaluated by transforming the witness event into the source rest frame, evaluating the rest-frame Gaussian Coulomb field, then boosting `E` and `B` back to the lab frame.
- The Boris pusher includes the magnetic field.
- In the default symmetric transverse launch, net `B` cancels by symmetry along the `z=0` trajectory.
- The asymmetric case samples nonzero magnetic field.
- The current `dt = 0.1 ps` pair table under-resolves the bunch-overlap transient.
- With `gamma_s = 480.453` and `sigma' = 0.6 mm`, `sigma_z/c ~= 0.00416 ps`.
- A 5 ps convergence check for the symmetric electron gives:
  - `dt=0.1 ps`: `dK ~= 5.20 eV`
  - `dt=0.01 ps`: `dK ~= 48.64 eV`
  - `dt=0.005 ps`: `dK ~= 58.68 eV`
  - `dt=0.002 ps`: `dK ~= 60.86 eV`
  - `dt=0.001 ps`: `dK ~= 61.16 eV`
- The note now recommends `dt ~= 1e-3 ps` or adaptive integration near the overlap event for production manufactured-solution comparisons.

Verification:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error boosted_gaussian_witness.tex
git diff --check -- sandbox/note/boosted_gaussian_witness.tex sandbox/note/boosted_gaussian_witness.pdf
```

Both passed. `pdftotext` was not installed, so PDF text extraction was not used for verification.

## 2026-06-10 Manufactured Solution Figure 2 Data Update

User asked to update Figure 2 data in the manufactured-solution PDF.

Updated the pair generator to support multirate integration:

```text
sandbox/python/boosted_gaussian_witness.py
sandbox/python/beambeam_analysis.py
```

New pair controls:

```text
--pair-fine-dt-ps
--pair-fine-duration-ps
--pair-output-dt-ps
```

Regenerated pair data/figures with:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python \
  sandbox/python/boosted_gaussian_witness.py \
  --pair-demo \
  --pair-fine-dt-ps 0.001 \
  --pair-fine-duration-ps 5.0 \
  --pair-dt-ps 0.1 \
  --pair-output-dt-ps 0.1 \
  --write-pair-latex-table sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex \
  --output-prefix sandbox/data/boosted_gaussian_witness
```

The regenerated Figure 2 asset is:

```text
sandbox/note/figs/boosted_gaussian_witness_physical_paths.png
```

and the regenerated table is:

```text
sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex
```

Resolved pair results now shown in the table/Figure 2:

- symmetric `e-`: `Delta r = (11.464, 0, 0) um`, `Delta K = +61.224 eV`
- symmetric `e+`: `Delta r = (11.466, 0, 0) um`, `Delta K = -61.220 eV`
- asymmetric `e-`: `Delta r = (10.317, 0, -1.461) um`, `Delta K = +55.102 eV`
- asymmetric `e+`: `Delta r = (10.319, 0, 1.462) um`, `Delta K = -55.098 eV`

The note text now describes Figure 2 as using `0.001 ps` steps for the first `5 ps`, then `0.1 ps` tail steps with `0.1 ps` output.

Verification:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile \
  sandbox/python/boosted_gaussian_witness.py sandbox/python/beambeam_analysis.py
git diff --check -- sandbox/python/boosted_gaussian_witness.py \
  sandbox/python/beambeam_analysis.py \
  sandbox/note/boosted_gaussian_witness.tex \
  sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error boosted_gaussian_witness.tex
```

All passed. The regenerated Figure 2 PNG was visually checked.

## 2026-06-10 OPALX-IMPACT Drift Staged-DT Check

User asked to test OPALX staged timesteps on the OPALX-IMPACT drift case with a mild factor-of-two timestep change and compare results.

Added:

```text
sandbox/OPALX-IMPACT/drift-1-staged-dt.in
```

The staged input keeps the same drift setup as `drift-1.in`, but changes:

```text
MAXSTEPS = {500, 1000}
DT       = {1.0e-10, 5.0e-11}
```

This is the same nominal integration time as the original uniform case:

```text
1000 * 1.0e-10 = 500 * 1.0e-10 + 1000 * 5.0e-11 = 1.0e-7 s
```

Commands run:

```sh
cd sandbox/OPALX-IMPACT
../../build_openmp/src/opalx drift-1.in
../../build_openmp/src/opalx drift-1-staged-dt.in
```

The full baseline completed and refreshed:

```text
sandbox/OPALX-IMPACT/drift-1.stat
sandbox/OPALX-IMPACT/drift-1.h5
```

The staged run completed and produced:

```text
sandbox/OPALX-IMPACT/drift-1-staged-dt.stat
sandbox/OPALX-IMPACT/drift-1-staged-dt.h5
```

Comparison artifacts:

```text
sandbox/OPALX-IMPACT/drift_1_staged_dt_comparison_summary.csv
sandbox/OPALX-IMPACT/drift_1_staged_dt_comparison.png
```

Switch verification from the staged stat file:

```text
base dt unique:   0.1 ns
staged dt unique: 0.1 ns, 0.05 ns
```

Result summary:

```text
base rows:   73
staged rows: 95
base final t:   72.29999999999963 ns
staged final t: 72.30000000000067 ns

base final rms_x:   0.03479139645306627 m
staged final rms_x: 0.03479139704389219 m
final Delta rms_x:  5.908259229081558e-10 m

base final eps_x:   1.929820568321231 mm mrad
staged final eps_x: 1.9298205441440939 mm mrad
final Delta eps_x: -2.4177137181169428e-08 mm mrad

common-time RMS Delta rms_x:       2.1377084664500679e-10 m
common-time max abs Delta rms_x:   5.908253886133252e-10 m
common-time RMS Delta eps_x:       8.868710224850992e-09 mm mrad
common-time max abs Delta eps_x:   2.4177137181169428e-08 mm mrad
```

Interpretation: the staged timestep switch happens and the drift result is numerically unchanged at the expected tiny level for this factor-of-two timestep change.

## 2026-06-10 BeamBeam Staged-DT Fix/Understanding

User asked to go back to BeamBeam after the OPALX-IMPACT drift check and first fix/understand different timesteps.

Findings from code inspection:

- `TRACK, DT={...}, MAXSTEPS={...}` creates sequential tracking segments.
- Segment switching is by accumulated step counts in `ParallelTracker::execute`, not by BeamBeam element boundaries.
- `ZSTOP` can remain scalar; `TrackCmd` pads scalar arrays to the number of tracks.
- `FROMFILE` emission does not override the global staged `DT`. It emits once when `[t, t+dt]` crosses `T0` and assigns newly emitted particles the remaining fractional substep `tEnd - t0`.
- Therefore, for pair birth to happen at fine resolution, the fine `DT` segment must begin before `tinj`.

Production staged BeamBeam input added:

```text
sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in
```

It is copied from `gamma_gamma_pairs-large-cylinder-retire1000ps.in` with this staged schedule:

```text
dt_coarse = 1.0e-12
dt_fine   = 1.0e-15

coarse_steps_before_overlap = 550
fine_overlap_steps          = 25001
coarse_steps_after_overlap  = 1025

MAXSTEPS = {550, 25001, 1025}
DT       = {1.0e-12, 1.0e-15, 1.0e-12}
```

Timing covered by the fine segment:

```text
coarse segment end: 550.000 ps
tinj:               550.312 ps
primary IP time:    567.059 ps
fine segment end:   575.001 ps
tail end:          1600.001 ps
```

Small BeamBeam staged-DT smoke input:

```text
sandbox/track-e-p/gamma_gamma_pairs-staged-dt-smoke.in
```

It uses 32 primary particles, 4 pairs per species, `N=4`, and:

```text
MAXSTEPS = {567, 1000}
DT       = {1.0e-12, 1.0e-15}
tinj     = primary_ip_time
```

Smoke validation run:

```sh
cd sandbox/track-e-p
../../build_openmp/src/opalx gamma_gamma_pairs-staged-dt-smoke.in
../../build_openmp/src/opalx --info 2 gamma_gamma_pairs-staged-dt-smoke.in
```

The run completed. Diagnostics showed source overlap and witness emission:

```text
BB-DIAG ... source_bunches_overlap=TRUE
BB-DIAG ... witness_states=c1:active:n=4,c2:active:n=4 source_bunches_overlap=TRUE
```

The final smoke stat samples reported:

```text
gamma_gamma_pairs-staged-dt-smoke_c0.stat dt unique ns: [1e-06], particles 32
gamma_gamma_pairs-staged-dt-smoke_c1.stat dt unique ns: [1e-06], particles 4
gamma_gamma_pairs-staged-dt-smoke_c2.stat dt unique ns: [1e-06], particles 4
```

`1e-06 ns` is `1 fs`, confirming the fine segment was active when witness particles existed.

Earlier full staged production run attempt:

- Started
  `../../build_openmp/src/opalx gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in`
  from `sandbox/track-e-p`.
- The original staged production input used `PSDUMPFREQ = 1` and
  `STATDUMPFREQ = 1`.
- The run confirmed the intended sequencing:
  - BeamBeam active before witness injection.
  - `c0.stat` showed `dt = 1.0e-06 ns`, i.e. `1 fs`, at about
    `550.17 ps`.
  - `BB-DIAG` reported `witness_states=c1:active:n=1297,c2:active:n=1297`
    in the fine segment.
- The run was stopped intentionally before completion because only about
  `3 GB` of disk space remained and dump-every-step output projected to tens
  of GB. Partial staged outputs from this aborted attempt were deleted.
- Updated
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in`
  to keep the same 1 fs integration but reduce output cadence:

```text
OPTION, PSDUMPFREQ = 100;
OPTION, STATDUMPFREQ = 100;
```

Additional note/PDF integration work during the staged run:

- Merged the regression overview content into
  `sandbox/note/boosted_gaussian_witness.tex` as a native
  "Sandbox Regression Overview" section.
- Copied the regression overview figures into `sandbox/note/figs` so the note
  no longer depends on `sandbox/regression` figure paths at compile time:
  - `opalx_impact_drift_comparison.png`
  - `gamma_gamma_cylinder_losses.png`
  - `gamma_gamma_large_cylinder_crossings.png`
  - `gamma_gamma_large_cylinder_charge_compare.png`
- Added a note subsection on staged timestep rationale:
  - OPALX treats `MAXSTEPS` as per-segment counts, not cumulative endpoints.
  - The fine 1 fs window brackets `tinj = primary_ip_time - 16.747 ps` and the
    primary IP crossing.
  - The histogram-focused run uses `PSDUMPFREQ = 10` and `STATDUMPFREQ = 0`.
  - Field debug `.dat` files are not used for the histogram and can be deleted.
- Rebuilt `sandbox/note/boosted_gaussian_witness.pdf`; build succeeded with
  only layout warnings from long paths/options.
- Checked and updated `sandbox/README.md` after the regression scripts were
  relocated.  The README now points to `sandbox/python/run_sandbox_regressions.py`,
  `sandbox/python/make_sandbox_regression_overview.py`, note-local
  `sandbox_regression_baseline.json`, note-local `current_metrics.csv`, and
  the integrated `sandbox/note/boosted_gaussian_witness.tex` PDF workflow.
- Rebuilt `sandbox/note/boosted_gaussian_witness.pdf` again after replacing
  stale `sandbox/regression` reproduction commands with `sandbox/python`
  commands.
- Validation after relocation:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python -m py_compile \
  sandbox/python/make_sandbox_regression_overview.py \
  sandbox/python/run_sandbox_regressions.py \
  sandbox/python/compare_cylinder_crossing_histograms.py \
  sandbox/python/make_large_cylinder_field_movie.py \
  sandbox/python/plot_cylinder_crossings.py
git diff --check -- sandbox/README.md sandbox/note/boosted_gaussian_witness.tex \
  sandbox/python/make_sandbox_regression_overview.py \
  sandbox/python/run_sandbox_regressions.py \
  sandbox/python/compare_cylinder_crossing_histograms.py \
  sandbox/python/make_large_cylinder_field_movie.py \
  sandbox/python/plot_cylinder_crossings.py
cd sandbox/note
latexmk -pdf -interaction=nonstopmode -halt-on-error boosted_gaussian_witness.tex
```

All three checks passed.  The LaTeX build still reports only layout warnings
from long paths/options.

Completed staged histogram run:

- Single-rank OPALX run completed:
  `../../build_openmp/src/opalx gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in`
- Current histogram-focused input settings are:

```text
OPTION, PSDUMPFREQ = 10;
OPTION, STATDUMPFREQ = 0;
MAXSTEPS = {550, 25001, 1025}
DT       = {1.0e-12, 1.0e-15, 1.0e-12}
```

- MPI two-rank attempt was tried to avoid field debug `.dat` dumps, but local
  OpenMPI mapping failed before OPALX started.
- During the run, periodically deleted only staged
  `sandbox/track-e-p/data/gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt-*`
  `.dat` byproducts.
- H5 witness trajectories remain the deliverable for the crossing histogram.
- OPALX reported completion and source retirement:

```text
BB-DIAG BB-state=Completed active_bunches=3 retired_bunches=0 witness_states=c1:active:n=1297,c2:active:n=1297 BB-active=FALSE source_retirement_pending=TRUE
BB-DIAG BB-state=Completed active_bunches=2 retired_bunches=1 witness_states=c1:active:n=1297,c2:active:n=1297 source_active=FALSE source_retirement_pending=FALSE
real 4850.86
```

- Final H5 summaries:

```text
c1 2572 samples, final time 1600.001 ps, final particle count 1297
c2 2572 samples, final time 1600.001 ps, final particle count 1297
```

- Final crossing histogram and low-charge overlay:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python sandbox/python/plot_cylinder_crossings.py \
  --compare-stem gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5 \
  --label "nominal staged" \
  --compare-label "1e-5 charge" \
  --title "1000 ps retirement: nominal staged vs 1e-5 charge" \
  --output sandbox/note/figs/gamma_gamma_large_cylinder_staged_dt_crossings.png
```

Output:

```text
sandbox/note/figs/gamma_gamma_large_cylinder_staged_dt_crossings.png
e-: N=701/701, |bin diff|=10, max |dz|=2026.5 um
e+: N=692/692, |bin diff|=12, max |dz|=1847.3 um
```

- A partial preview was also generated before the run completed:
  `sandbox/note/figs/gamma_gamma_large_cylinder_staged_dt_crossings_preview.png`;
  at that time both counts were zero because no particles had yet reached the
  15 cm cylinder radius.
- Added the final staged-DT nominal histogram to
  `sandbox/note/boosted_gaussian_witness.tex` and rebuilt
  `sandbox/note/boosted_gaussian_witness.pdf` successfully.  The PDF is now
  20 pages.  LaTeX still reports only layout warnings from long paths/options.
- The overlaid low-charge curve uses the existing
  `gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5` H5 files, not a new
  staged-DT rerun.  Those H5 files have 1050 stored samples per witness and
  also end at 1600 ps.

## BeamBeam mirrored E-field validation

Goal: validate the current BeamBeam mirror model visually.  The model computes
the self field from primary container/bunch 0 and mirrors that field/source
around the interaction point.  For this validation we returned to the existing
uniform coarse `DT = 1 ps` large-cylinder data rather than the staged-DT run.

Input/data used:

- Existing coarse debug dumps under `sandbox/track-e-p/data`:
  `gamma_gamma_pairs-large-cylinder-EF_vector-beambeam_e-*.dat`
- These dumps cover active BeamBeam snapshots from about `73 ps` to `571 ps`.
- The corresponding input is the uniform coarse large-cylinder setup:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in`
  with `MAXSTEPS = 1600` and `DT = 1.0e-12`.

Added reusable converter/visualizer:

```text
sandbox/python/convert_efield_dumps_to_h5.py
```

Command run:

```sh
/Users/adelmann/git/opalx-beambeam/.venv-h6/bin/python \
  sandbox/python/convert_efield_dumps_to_h5.py --gallery-frames 12
```

Outputs:

```text
sandbox/data/gamma_gamma_large_cylinder_efield_debug.h5
sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug.png
sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug_ez.png
sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug_lab_timeline.png
sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug_z_profile.png
```

The H5 contains one group per EF dump with datasets `x`, `y`, `z`,
`E/Ex`, `E/Ey`, `E/Ez`, and `E/Eabs`, plus the original dump metadata as H5
attributes.  It also records explicit marker positions in both the
BeamBeam-local frame (`primary_local_z_m`, `mirror_local_z_m`) and the
lab/global frame (`primary_lab_z_m`, `mirror_lab_z_m`,
`interaction_point_lab_z_m`).  The PNG galleries show a central-y `x-z` slice
on the BeamBeam-local dump grid.  The vertical markers are:

- white: local IP plane;
- cyan: bunch[0] centroid from `particle_mean_r[2]`;
- magenta dashed: inferred mirror centroid `2*interaction_point_local_z -
  particle_mean_r[2]`.

Important coordinate note: the field galleries use local dump coordinates, not
lab coordinates.  In this local frame, `bunch[0]` remains close to the local
origin while the IP plane sweeps through the mesh.  The lab-frame timeline
confirms the physical motion: `bunch[0]` moves in increasing global `z/s`, and
the mirrored source moves in decreasing global `z/s`, toward the fixed IP at
`s = 0.169981955 m`.

The `x-z` heatmaps of `Ez` or `|E|` can look like transverse `x` motion because
the bunch field is dominated by off-axis transverse lobes in a central-y slice.
The longitudinal profile plot integrates `|E|` over transverse `x-y`; it shows
the two field peaks moving toward each other along local `z`.

Representative metadata from the generated H5:

```text
step=75  time=73.000 ps  local: ip_z=0.148097 primary_z=-0.000120068 mirror_z=0.296314; lab: primary_z=0.0217647 mirror_z=0.318199
step=324 time=322.000 ps local: ip_z=0.073449 primary_z=-0.00104215  mirror_z=0.14794;  lab: primary_z=0.0954908 mirror_z=0.244473
step=573 time=571.000 ps local: ip_z=-0.00119917 primary_z=-0.00196425 mirror_z=-0.000434081; lab: primary_z=0.169217 mirror_z=0.170747
```

The current field dumps stop essentially at overlap, so they verify the
incoming z-direction geometry but do not yet show a post-overlap left/right
exchange.

Follow-up 1000 ps retirement run:

- Reran the uniform coarse `DT = 1 ps` input
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps.in`.
- This input keeps the primary source active until `RETIRE_TIME = 1000e-12`
  while preserving the large-cylinder witness injection timing.
- The run completed successfully.  OPALX reported the expected diagnostics:
  active BeamBeam window, witness injection, source overlap becoming true, then
  false again after overlap, completed with source-retirement pending, and final
  source retirement.
- The regenerated EF debug conversion read `927` dumps for stem
  `gamma_gamma_pairs-large-cylinder-retire1000ps`, ending at about `999 ps`.
- New debug products:
  - `sandbox/data/gamma_gamma_large_cylinder_retire1000ps_efield_debug.h5`
  - `sandbox/note/figs/gamma_gamma_large_cylinder_retire1000ps_efield_debug.png`
  - `sandbox/note/figs/gamma_gamma_large_cylinder_retire1000ps_efield_debug_ez.png`
  - `sandbox/note/figs/gamma_gamma_large_cylinder_retire1000ps_efield_debug_lab_timeline.png`
  - `sandbox/note/figs/gamma_gamma_large_cylinder_retire1000ps_efield_debug_z_profile.png`
- Representative metadata:
  - `step=75`, `time=73 ps`: lab `primary_z=0.0217647`,
    `mirror_z=0.318199`
  - `step=538`, `time=536 ps`: lab `primary_z=0.158854`,
    `mirror_z=0.18111`
  - `step=1001`, `time=999 ps`: lab `primary_z=0.295943`,
    `mirror_z=0.0440212`
- The 1000 ps longitudinal profile plot shows the pass-through: after overlap,
  the local marker ordering reverses, with `bunch[0]` on the downstream/right
  side and the mirrored source on the upstream/left side.

Input inventory/update:

- Added `sandbox/track-e-p/INPUT_INVENTORY.md`, which documents each active and
  legacy BeamBeam `.in` file and its purpose.
- Updated the default nominal and reduced-charge large-cylinder inputs to match
  the current field-debug configuration:
  - `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in`
  - `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in`
- These now use the same large-cylinder timing as the explicit
  `retire1000ps` variants: `primary_retire_time = 1000e-12`, preserving
  `tinj = primary_ip_time - 16.747e-12`.
- Left `gamma_gamma_pairs-2.in` unchanged as the small proof-of-principle
  cylinder case with `RETIRE_TIME = 121 ps`.
- Left `gamma_gamma_pairs-3.in` and `attic/gamma_gamma_pairs.in` unchanged as
  legacy/historical inputs.

OpenMP 2-thread/4-thread stability and timing check:

- Created benchmark input copies:
  - `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-omp2.in`
  - `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-omp4.in`
- Both are content copies of the nominal `1000 ps` coarse input; the separate
  filenames isolate OPALX output stems.
- Commands:

```sh
/usr/bin/time -p env OMP_NUM_THREADS=2 OMP_PROC_BIND=spread OMP_PLACES=threads \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx \
  gamma_gamma_pairs-large-cylinder-retire1000ps-omp2.in
/usr/bin/time -p env OMP_NUM_THREADS=4 OMP_PROC_BIND=spread OMP_PLACES=threads \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx \
  gamma_gamma_pairs-large-cylinder-retire1000ps-omp4.in
```

- Timings:
  - 2 threads: `real 236.60 s`, `user 292.74 s`, `sys 24.38 s`
  - 4 threads: `real 205.75 s`, `user 297.51 s`, `sys 38.03 s`
  - 4-vs-2 thread wall-time speedup: about `1.15x`
- Both runs completed with the expected BeamBeam diagnostics and final source
  retirement.
- Final stats were stable at the population level:
  - c0 retired, `numParticles = 0`
  - c1/c2 active, `numParticles = 1297`, `partsOutside = 0`
- Final witness H5 states differ slightly between thread counts after sorting
  by particle id:
  - c1 max absolute differences `[x,y,z,px,py,pz]` =
    `[4.02e-6 m, 4.24e-6 m, 2.40e-6 m, 3.07e-5, 3.54e-5, 2.20e-5]`
  - c2 max absolute differences `[x,y,z,px,py,pz]` =
    `[4.20e-6 m, 4.58e-6 m, 2.26e-6 m, 3.33e-5, 3.81e-5, 2.15e-5]`
- Primary c0 particle-by-particle states are not stable across thread counts,
  likely from OpenMP-sensitive source sampling/ordering; the primary is retired
  by the final stat row, so this does not affect the final witness crossing
  histogram directly.
- Cylinder crossing histogram comparison:
  - Output:
    `sandbox/note/figs/gamma_gamma_large_cylinder_omp2_vs_omp4_crossings.png`
  - e-: `N=701/701`, `|bin diff|=0`, max first-crossing `|dz|=180.3 um`
  - e+: `N=691/691`, `|bin diff|=0`, max first-crossing `|dz|=195.5 um`

Implementation inspection in `/Users/adelmann/git/opalx`:

- `src/PartBunch/PartBunch.cpp` expands the mesh bounds using
  `mirroredMinZ = 2.0 * planeZ - upper[2]` and
  `mirroredMaxZ = 2.0 * planeZ - lower[2]`.
- `src/PartBunch/ImageChargeScatterController.tpp` applies
  `rView(i)[2] = 2.0 * planeZ - rView(i)[2]`.
- The same image path then calls `flipChargeSignAll`, so the mirrored deposit
  has opposite charge sign.

Current interpretation: the geometry of the left/right mirror is consistent
with `z' = 2*z_IP - z` and the visual debug confirms the mirror marker crosses
the source marker around the IP.  However, the current implementation path is
an image-charge path and flips the sign of the mirrored charge.  For two
colliding electron primary beams, this sign flip is suspicious and should be
checked against the intended BeamBeam physics model before trusting witness
kicks.

2026-06-14 COPY_TIME sandbox rerun status:

- Working tree: `/Users/adelmann/git/opalx-beambeam`, branch
  `271-implement-interation-point-element`.
- BeamBeam `COPY` boolean has been replaced by `COPY_TIME`; current sandbox
  inputs use `COPY_TIME = 100e-12`.
- Built `build_openmp/src/opalx` after the `COPY_TIME` changes and verified the
  smoke input accepts the new attribute.
- Completed coarse nominal run:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in`.
- Completed coarse reduced-charge run:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in`.
- Regenerated Figure 8 preview only, without copying into notes:
  `/tmp/gamma_gamma_large_cylinder_charge_compare_preview.png`.
  Result: e- and e+ crossing histograms are bin-identical between nominal and
  reduced charge for this coarse rerun; max first-crossing shifts are about
  340.4 um (e-) and 324.0 um (e+).
- Active run at handoff time:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in`
  writing `gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt_c*.h5`.
- The staged-DT run completed cleanly, with BeamBeam diagnostics showing
  `source_bunches_overlap=TRUE`, then `FALSE`, followed by primary retirement.
- Regenerated Figure 9 preview only, without copying into notes:
  `/tmp/gamma_gamma_large_cylinder_staged_dt_crossings_preview.png`.
  Result: first cylinder-edge crossings count `N=172` for e- and `N=142` for e+.
- Both Figure 8 and Figure 9 previews were shown from `/tmp`; no files were
  copied to `sandbox/note/figs` yet.

2026-06-14 COPY_TIME 50 ps rerun status:

- Updated all active BeamBeam sandbox `.in` files under `sandbox/track-e-p` from
  `COPY_TIME = 100e-12` to `COPY_TIME = 50e-12`.
- Smoke run completed:
  `sandbox/track-e-p/gamma_gamma_pairs-staged-dt-smoke.in`.
- Completed coarse nominal run:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder.in`.
- Completed coarse reduced-charge run:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-q1em5.in`.
- Regenerated Figure 8 preview only:
  `/tmp/gamma_gamma_large_cylinder_charge_compare_copy50_preview.png`.
  Result: e- and e+ crossing histograms remain bin-identical between nominal and
  reduced charge; max first-crossing shifts are about 340.4 um (e-) and
  324.0 um (e+).
- Completed staged-DT run:
  `sandbox/track-e-p/gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in`
  with `COPY_TIME = 50e-12`, writing
  `gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt_c*.h5`.
- BeamBeam diagnostics showed `copy_active=TRUE` before BB activation,
  `source_bunches_overlap=TRUE`, then `FALSE`, followed by primary retirement.
- Regenerated Figure 9 preview only:
  `/tmp/gamma_gamma_large_cylinder_staged_dt_crossings_copy50_preview.png`.
  Result: first cylinder-edge crossings count `N=172` for e- and `N=142` for e+.
- Both COPY_TIME 50 ps preview figures were shown from `/tmp`; no files were
  copied to `sandbox/note/figs` yet.
- Plotting correction after comparing with note figures: the H5 crossing `z`
  values are local BeamBeam coordinates for this diagnostic.  The note figures
  use the sandbox coordinate, so `sandbox/python/compare_cylinder_crossing_histograms.py`
  and `sandbox/python/plot_cylinder_crossings.py` now add a default
  `--z-offset-m 0.33`.  Regenerated previews:
  - `/tmp/gamma_gamma_large_cylinder_charge_compare_copy50_preview.png`
  - `/tmp/gamma_gamma_large_cylinder_staged_dt_crossings_copy50_preview.png`
  Full counts are restored: e- `N=701`, e+ `N=692`.
