# Codex Handoff: OPALX Element Placement and Chicane Regression

Date: 2026-06-09

## 1. Original Goal

Finish the OPALX element-placement work for analytic `SBEND`/`RBEND` support,
with special focus on body-pose placement, old-OPAL compatibility placement,
and regression tests that can show whether OPALX and OPAL 2022.1 agree for a
C-type chicane at the moment hard-edge. 

The immediate sandbox goal was to build and document a distribution-level
chicane regression:

- keep the bend model analytic and hard-edge;
- do not introduce field maps yet;
- compare OPALX body-pose placement, OPALX `ELEMEDGE` compatibility placement,
  and old OPAL 2022.1 reference-line placement;
- make publication-quality comparison plots for RMS sizes, RMS momenta,
  emittances, and reference `By` field;
- record the Python tools and current findings so the work can restart cleanly.

## 2. Current Branch and Repo State

Repository: `/Users/adelmann/git/opalx`

Current branch:

```text
391-placement-of-sbend-and-rbend-analytic-fields-including-fringe-fields
```

Current HEAD:

```text
eee613d Save some tests
```

The branch is currently aligned with its upstream:

```text
HEAD -> 391-placement-of-sbend-and-rbend-analytic-fields-including-fringe-fields
origin/391-placement-of-sbend-and-rbend-analytic-fields-including-fringe-fields -> 44121ac2b
```

The worktree is dirty.  `git status --short` currently shows tracked
modifications in:

```text
cmake/Dependencies.cmake
sandbox/test-rbend-1/test-rbend-1.stat
sandbox/test-rbend-1/test-rbend-1_exit_monitor.loss
sandbox/test-rbend-2/test-rbend-2.stat
sandbox/test-rbend-2/test-rbend-2_exit_monitor.loss
src/PartBunch/FieldMirror.hpp
```

There are also many untracked sandbox and documentation files, including:

```text
AGENT.md
KS.md
hooks/post-index-change
sandbox/ElementPlacing/
sandbox/LaTeX/
sandbox/OPALX-IMPACT/
sandbox/pedro/
```

The most recent Codex work in this session is under:

```text
sandbox/test-chicane-distribution-1/
```

## 3. Files Changed and Why

### Chicane sandbox files

`sandbox/test-chicane-distribution-1/README.md`

- Newly expanded handoff-style documentation for the chicane test.
- Documents the analytic model, manufactured `R56`, input decks, Python tools,
  run commands, DT scan, plot locations, and current fitted values.

`sandbox/test-chicane-distribution-1/compare_opal_opalx.py`

- Extended to generate the same plot set for small and 10000-particle
  OPALX-vs-OPAL comparisons.
- Added Matplotlib publication-style plots for:
  - RMS sizes with `|\Delta| [%]`;
  - RMS momenta with `|\Delta| [%]`;
  - emittances with `|\Delta| [%]`;
  - `By_ref` with `|\Delta| [%]`.
- Relative-difference panels now show unsigned percent difference:
  `100 * abs(OPALX - OPAL) / abs(OPAL)`.
- Plot labels intentionally use only `|\Delta| [%]`, not the full formula.
- The `By_ref` relative curve masks near-zero field denominators so hard-edge
  zero crossings do not dominate the scale.  Edge timing shifts are captured in
  `by_field_intervals.csv`.

Generated comparison outputs:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_5p0em12/compare_elemedge/
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/compare_elemedge/
```

Each contains:

```text
summary.md
stat_summary.csv
by_field_intervals.csv
stat_rms_relative_opalx_opal.png/.pdf
stat_rms_momentum_difference_opalx_opal.png/.pdf
stat_emittance_difference_opalx_opal.png/.pdf
stat_bfield_ref_opalx_opal.png/.pdf
stat_energy_opalx_opal.png
stat_reference_xz_opalx_opal.png
```

### Other dirty tracked files

`cmake/Dependencies.cmake`

- Current diff unsets stale `Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK` before
  IPPL/Kokkos configuration.  This was from a separate build issue and is not
  part of the chicane sandbox.

`src/PartBunch/FieldMirror.hpp`

- Current diff changes `solver_recv` call arguments from a `Vector<bool, Dim>`
  to explicit `flipX`, `flipY`, `flipZ` booleans.  This is also from separate
  work and is not part of the chicane sandbox.

## 4. Important Design Decisions

- The primary OPALX placement model is body-pose placement using `X/Y/Z/THETA`.
  This is the target style for OPALX.
- `ELEMEDGE` is retained only as a compatibility/control path.  It is useful
  for comparing OPALX against old OPAL reference-line placement, but it should
  not become the new OPALX design direction.
- The chicane regression stays analytic and hard-edge:
  - no field maps;
  - no fringe fields;
  - no pole-face focusing;
  - `E1 = E2 = HGAP = FINT = 0`.
- The manufactured longitudinal target uses the small-angle formula:

```text
R56_ref = -TH^2 * (4 LB / 3 + 2 D_CH)
```

For the current chicane constants this is:

```text
R56_ref = -4.333333333333334e-2 m
```

- The OPALX manufactured regression compares the fitted slope after removing a
  monitor-frame intercept.  The raw `z` residual includes a fitted monitor
  offset and should not be interpreted directly.
- The OPALX-vs-OPAL comparison is a code-to-code diagnostic, not an exact
  truth model.  OPALX and OPAL 2022.1 use different placement/model concepts.
- For plots, the requested convention is `|\Delta| [%]` everywhere.  Do not
  relabel these panels with the full formula or with signed differences unless
  the user asks.

## 5. Constraints and Invariants That Must Not Be Broken

- Do not revert unrelated dirty files.  The worktree already contains user and
  previous-session changes.
- Keep the chicane test hard-edge and analytic until field maps/fringes are
  explicitly requested.
- Preserve the distinction between:
  - OPALX body-pose deck: `test-chicane-distribution-1.in`;
  - OPALX `ELEMEDGE` compatibility deck: `test-chicane-distribution-2.in`;
  - OPAL 2022.1 deck: `test-chicane-distribution-1_opal.in`.
- Preserve the deterministic particle distribution unless intentionally
  creating a new test variant.
- Keep `FIELDSOLVER = NONE` and `BCHARGE = 1e-9 C` for the distribution test.
  Earlier attempts without charge could not integrate the particles correctly.
- Do not use `ELEMEDGE` in new OPALX body-pose examples unless the task is
  explicitly about OPAL compatibility.
- For the comparison plots:
  - OPAL is interpolated onto OPALX `s` samples;
  - relative curves are unsigned percent differences;
  - labels are `|\Delta| [%]`;
  - the `By_ref` relative plot masks near-zero denominators.
- The 10000-particle input files now contain `STATDUMPFREQ=10`, but the current
  10000-particle `.stat` files still have 6304 samples from the dense run made
  before that change.  Rerun before relying on the lower dump frequency.

## 6. Commands/Tests Already Run and Results

### Manufactured OPALX regression

Commands used from `sandbox/test-chicane-distribution-1`:

```sh
/Users/adelmann/.venv-h6/bin/python generate_distribution.py
/Users/adelmann/git/opalx/build/src/opalx --info 5 test-chicane-distribution-1.in
/Users/adelmann/git/opalx/build/src/opalx test-chicane-distribution-1.in
/Users/adelmann/.venv-h6/bin/python compare_results.py --plots
```

Result from `comparison/comparison_summary.md`:

```text
Selected monitor step: Step#1
Particles at exit: 125
Reference R56: -4.333333333333e-02 m
Fitted transport R56: -4.327467367520e-02 m
Difference: 5.865965813612e-05 m
Intercept-corrected z residual RMS: 4.785727528281e-05 m
Intercept-corrected z residual max abs: 8.407208624630e-05 m
y residual RMS: 0
Overall: PASS
```

The `--info 5` run confirmed finite initial distribution bounds:

```text
MIN R ~= (-0.10000, -0.08000, -1.00000) mm
MAX R ~= ( 0.10000,  0.08000,  1.00000) mm
```

### DT scan

The scan generated results under:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/
```

Representative selected RMS-difference summary:

```text
finite/body     DT=2e-12: By_ref rms=4.159428e-03 T, rms_s rms=1.017762e-05 m
finite/body     DT=5e-12: By_ref rms=1.076316e-02 T, rms_s rms=4.698987e-06 m
finite/elemedge DT=2e-12: By_ref rms=5.518510e-03 T, rms_s rms=4.903372e-07 m
finite/elemedge DT=5e-12: By_ref rms=1.158390e-02 T, rms_s rms=3.531604e-07 m
```

The `ELEMEDGE` compatibility path strongly reduces the longitudinal RMS
discrepancy versus OPAL compared with OPALX body-pose placement, which is
expected because it mimics the old OPAL reference-line model more closely.

### Small 125-particle OPALX-vs-OPAL comparison

Command used to regenerate plots:

```sh
PYTHONPYCACHEPREFIX=/private/tmp/opalx-pycache \
MPLCONFIGDIR=/private/tmp/opalx-matplotlib \
/Users/adelmann/.venv-h6/bin/python \
/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/compare_opal_opalx.py \
  --opalx-stat /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat \
  --opal-stat /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat \
  --out-dir /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_5p0em12/compare_elemedge \
  --r56-summary OPALX=/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_5p0em12/compare_elemedge/r56_opalx/comparison_summary.json \
  --r56-summary OPAL=/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_5p0em12/compare_elemedge/r56_opal/comparison_summary.json \
  --plots
```

Result:

```text
OPALX fitted R56: -4.44581726e-02 m
OPAL fitted R56:  -4.44415648e-02 m
energy max abs diff: 3.449259e-10 MeV
rms_x final diff: 4.640679e-07 m
rms_y final diff: -5.313809e-09 m
rms_s final diff: 4.814756e-07 m
```

`By_ref` interval shifts for the small `ELEMEDGE` comparison:

```text
interval 1 start diff ~= 0 mm, end diff ~= 0 mm
interval 2 start/end diff ~= -1.498962 mm
interval 3 start diff ~= -2.997924 mm, end diff ~= -1.498962 mm
interval 4 start diff ~= -2.997924 mm, end diff ~= -4.496886 mm
```

### 10000-particle OPALX-vs-OPAL comparison

Command used to regenerate plots:

```sh
PYTHONPYCACHEPREFIX=/private/tmp/opalx-pycache \
MPLCONFIGDIR=/private/tmp/opalx-matplotlib \
/Users/adelmann/.venv-h6/bin/python \
/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/compare_opal_opalx.py \
  --opalx-stat /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat \
  --opal-stat /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat \
  --out-dir /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/compare_elemedge \
  --r56-summary OPALX=/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/compare_elemedge/r56_opalx/comparison_summary.json \
  --r56-summary OPAL=/Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/compare_elemedge/r56_opal/comparison_summary.json \
  --plots
```

Result:

```text
OPALX fitted R56: -4.45195920e-02 m
OPAL fitted R56:  -4.44786273e-02 m
energy max abs diff: 3.469722e-10 MeV
rms_x final diff: 1.747598e-05 m
rms_y final diff: -4.463877e-09 m
rms_s final diff: 5.311967e-08 m
```

Plots were visually checked after regeneration.  The right-axis label issue was
fixed; labels now read `|\Delta| [%]`.

### Other tests/builds from earlier branch work

Earlier in the branch, CUDA/HIP build issues were addressed around bend tests:

- CUDA error in `unit_tests/BeamlineCore/TestBendRep.cpp` from extended
  `__host__ __device__` lambdas inside private/protected test fixture methods.
- HIP error in `src/AbsBeamline/BendBase.cpp` from host-only helper functions
  called inside a device lambda.

Those fixes were already committed before the current sandbox work.  The
current branch HEAD is pushed to origin.

## 7. Current Failures or Suspicious Behavior

- The manufactured OPALX regression passes, but the fitted `R56` from the
  OPALX/OPAL `ELEMEDGE` code-to-code comparisons remains about
  `1.1e-3 m` more negative than the small-angle manufactured reference.
- The `By_ref` interval table still shows millimeter-scale hard-edge shifts
  between OPALX and OPAL:
  - roughly one half-step to one and a half time steps depending on edge;
  - one time step at `DT=1e-11 s` is about `2.99792458 mm`.
- The first/entrance field region agrees, while later bend intervals show
  systematic shifts.  This suggests remaining reference-line/path-length,
  monitor, or hard-edge boundary convention differences.
- The 10000-particle outputs are generated from dense stat dumps despite the
  current 10000-particle input files now having `STATDUMPFREQ=10`.
- `sandbox/test-chicane-distribution-1/compare_opal_opalx.py` uses Matplotlib
  for the newer publication plots but still has older simple PIL plotting code
  for some legacy plots.  This is functional but stylistically mixed.
- The branch worktree contains unrelated dirty changes in CMake/FieldMirror and
  generated sandbox outputs.  Future agents must not assume all dirty files are
  part of the chicane work.
- The 10000 runs perfrom much worth than 125 this needs to investigated further

## 8. Exact Next Steps, in Priority Order

1. Re-run the 10000-particle `ELEMEDGE` comparison with the current
   `STATDUMPFREQ=10` inputs if smaller stat files are desired:

```sh
cd /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge
/Users/adelmann/git/opalx/build/src/opalx test-chicane-distribution-2.in

cd /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022
OPAL_PREFIX=/Users/adelmann/OPAL-2022.1 /Users/adelmann/OPAL-2022.1/bin/opal test-chicane-distribution-1_opal.in
```

Then regenerate plots with the `compare_opal_opalx.py` command shown above.

2. Investigate the remaining `By_ref` interval shifts:
   - compare exact bend edge positions in OPALX and OPAL generated placement
     files;
   - check whether shifts match `DT * beta * c` sampling, inclusive/exclusive
     edge handling, or different field boundary conventions;
   - inspect whether OPALX `ELEMEDGE` uses exactly the same hard-edge start/end
     definition as OPAL 2022.1.

3. Investigate why fitted `R56` in the code-to-code comparison differs from the
   small-angle manufactured value by about `1.1e-3 m`:
   - verify formula applicability for chord length vs arc length;
   - compare using exact chicane transfer formula rather than small-angle
     approximation;
   - confirm whether the monitor/exit plane and `z`/`s` definitions match the
     manufactured model.

4. Clean up `compare_opal_opalx.py` if it will be committed:
   - keep Matplotlib output for all comparison plots;
   - remove unused `difference` tuple entries and legacy PIL code if no longer
     needed;
   - preserve existing `|\Delta| [%]` labeling.

5. Add a lightweight regression test only after the physical comparison target
   is agreed:
   - manufactured OPALX-only check is currently the safest first regression;
   - OPALX-vs-OPAL plots are diagnostics and should not become strict pass/fail
     without agreed tolerances.

6. Before final commit, run at least:

```sh
cmake --build /Users/adelmann/git/opalx/build --target opalx_exe -j 8
/Users/adelmann/git/opalx/build/src/opalx /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/test-chicane-distribution-1.in
/Users/adelmann/.venv-h6/bin/python /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/compare_results.py --plots
```

Also run the relevant bend unit/regression tests if committing code outside the
sandbox.

## 9. Things Still Uncertain

- Whether the `R56` comparison should use the current small-angle formula or a
  more exact hard-edge chicane expression for the chosen geometry.
- Whether OPALX `ELEMEDGE` compatibility should exactly reproduce OPAL
  reference-line hard-edge edge timing, or only be close enough to isolate
  body-pose placement differences.
- Whether the millimeter-level `By_ref` interval shifts are purely time-step
  sampling artifacts or encode a real placement convention mismatch.
- Whether the 10000-particle comparison should be kept in the repository or
  regenerated on demand because the stat/H5/plot files are large.
- Whether the dirty tracked changes in `cmake/Dependencies.cmake` and
  `src/PartBunch/FieldMirror.hpp` are intended to be part of the next commit.
- Whether generated sandbox outputs in `sandbox/test-rbend-*`,
  `sandbox/test-sbend-*`, and `sandbox/test-chicane-*` should be removed,
  ignored, or committed as regression artifacts.

## 10. 2026-06-10 Current Session Update

User asked to proceed first with:

- direct audit of `sandbox/test-chicane-distribution-1`;
- reproduction of the OPALX-only manufactured chicane regression.

Files inspected directly:

```text
sandbox/test-chicane-distribution-1/README.md
sandbox/test-chicane-distribution-1/test-chicane-distribution-1.in
sandbox/test-chicane-distribution-1/manufactured_chicane.py
sandbox/test-chicane-distribution-1/compare_results.py
sandbox/test-chicane-distribution-1/comparison/comparison_summary.md
```

Findings from audit:

- The sandbox still matches the handoff description: analytic hard-edge
  four-`SBEND` C-type chicane, deterministic 125-particle distribution, finite
  transverse extent, `FIELDSOLVER = NONE`, `BCHARGE = 1e-9 C`, and no fringe or
  pole-face terms.
- `test-chicane-distribution-1.in` is the OPALX body-pose deck using
  `X/Y/Z/THETA`; `test-chicane-distribution-2.in` remains the `ELEMEDGE`
  compatibility/control deck.
- The OPALX-only manufactured check still uses the small-angle target
  `R56_REF = -4.333333333333334e-2 m` and compares the fitted slope after
  removing the monitor-frame intercept.

Commands run from `sandbox/test-chicane-distribution-1`:

```sh
/Users/adelmann/.venv-h6/bin/python generate_distribution.py
/Users/adelmann/git/opalx/build/src/opalx --info 5 test-chicane-distribution-1.in
/Users/adelmann/git/opalx/build/src/opalx test-chicane-distribution-1.in
/Users/adelmann/.venv-h6/bin/python compare_results.py --plots
```

`--info 5` confirmed the expected finite loaded distribution bounds:

```text
MIN R = (-0.10000, -0.08000, -1.00000) mm
MAX R = ( 0.10000,  0.08000,  1.00000) mm
```

Regenerated manufactured comparison result:

```text
Selected monitor step: Step#1
Monitor SPOS: 9.253200162550e+00 m
Particles at exit: 125
Reference R56: -4.333333333333e-02 m
Fitted transport R56: -4.327467367520e-02 m
Difference: 5.865965813020e-05 m
Fitted monitor offset: -1.463199339425e-03 m
Raw z residual RMS: 1.463981770956e-03 m
Intercept-corrected z residual RMS: 4.785727528130e-05 m
Intercept-corrected z residual max abs: 8.407208623454e-05 m
y residual RMS: 0.000000000000e+00 m
Overall: PASS
```

This reproduces the previous PASS baseline to last-digit floating-point
differences only.  Next useful step is the handoff item 8.1/8.4: investigate
the 10000-particle `ELEMEDGE` comparison and determine whether the worse
behavior comes from stale dense outputs, finite-statistics behavior, or a real
model/placement difference.

### 10000-particle `ELEMEDGE` investigation

The 10000-particle `ELEMEDGE` input decks already contained
`STATDUMPFREQ = 10`, but the `.stat` files were stale:

```text
opalx_elemedge/test-chicane-distribution-2.in   newer than .stat
opal_2022/test-chicane-distribution-1_opal.in   newer than .stat
old comparison samples: 6304
```

Fresh reruns were made from:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022
```

Commands:

```sh
/Users/adelmann/git/opalx/build/src/opalx test-chicane-distribution-2.in
OPAL_PREFIX=/Users/adelmann/OPAL-2022.1 \
  /Users/adelmann/OPAL-2022.1/bin/opal test-chicane-distribution-1_opal.in
```

Then the R56 summaries and OPALX-vs-OPAL plots were regenerated with:

```sh
/Users/adelmann/.venv-h6/bin/python compare_results.py \
  --distribution .../opalx_elemedge/test-chicane-distribution-1_distribution.txt \
  --monitor .../opalx_elemedge/test-chicane-distribution-2_exit.h5 \
  --out-dir .../compare_elemedge/r56_opalx --plots

/Users/adelmann/.venv-h6/bin/python compare_results.py \
  --distribution .../opal_2022/test-chicane-distribution-1_distribution.txt \
  --monitor .../opal_2022/test-chicane-distribution-1_opal_exit.h5 \
  --out-dir .../compare_elemedge/r56_opal --plots

/Users/adelmann/.venv-h6/bin/python compare_opal_opalx.py \
  --opalx-stat .../opalx_elemedge/test-chicane-distribution-2.stat \
  --opal-stat .../opal_2022/test-chicane-distribution-1_opal.stat \
  --out-dir .../compare_elemedge \
  --r56-summary OPALX=.../compare_elemedge/r56_opalx/comparison_summary.json \
  --r56-summary OPAL=.../compare_elemedge/r56_opal/comparison_summary.json \
  --plots
```

Fresh stat comparison now has 630 shared samples, not 6304.  The sampled
`By_ref` interval shifts essentially vanish:

```text
By_ref max_abs_diff = 3.239107893016069e-08 T
By_ref rms_diff     = 1.290602775731538e-09 T
interval start/end differences are O(1e-8) mm
```

This means the old millimeter-level 10000-particle `By_ref` interval table was
at least partly a stale dense-output artifact.  Be careful: the comparison text
still says "one-step shift is about 2.99792458 mm for DT = 1e-11 s", which is
not correct for this `DT = 5e-12 s`, `STATDUMPFREQ = 10` sparse output.

The larger horizontal/R56 behavior remains after rerunning:

```text
OPALX R56 = -4.451959196508e-02 m, Delta vs reference = -1.186258631743e-03 m
OPAL  R56 = -4.447862731188e-02 m, Delta vs reference = -1.145293978543e-03 m
OPALX-OPAL R56 difference ~= -4.0965e-05 m
OPALX manufactured centered z residual RMS = 2.554285459839e-06 m
OPAL  manufactured centered z residual RMS = 2.935860475369e-06 m
OPALX y residual RMS = 1.455387964380e-08 m
OPAL  y residual RMS = 1.653817787715e-08 m
```

So the 10000 case is not failing because of noisy longitudinal fits.  Both
codes have small centered residuals but both fitted `R56` values remain about
`1.1e-3 m` more negative than the current small-angle reference.

Distribution diagnostics:

```text
small manufactured distribution:
  n = 125
  unique x,y,z,delta counts = 5,5,5,5
  sigma_x = 7.071067811865e-05 m
  sigma_y = 5.656854249492e-05 m
  sigma_z = 7.071067811865e-04 m
  sigma_delta = 1.414213562373e-03
  corr(x,y) = 0.6

large distribution:
  n = 10000
  unique x,y,z,delta counts = 10,10,10,10
  sigma_x = 6.382847385042e-05 m
  sigma_y = 5.106277908034e-05 m
  sigma_z = 6.382847385042e-04 m
  sigma_delta = 1.276569477008e-03
  corr(x,y) ~= 0
  correlations of x/y/z with delta are numerically zero
```

The 10000 distribution is therefore a different transverse/longitudinal grid,
not just a higher-statistics copy of the 125-particle distribution.

Per-initial-`z` R56 diagnostics show strong slice dependence in the large case:

```text
large OPALX ELEMEDGE per-z R56 range:
  min -4.505465445424e-02 m
  max -4.370799097248e-02 m

large OPAL 2022 per-z R56 range:
  min -4.504864151791e-02 m
  max -4.381170372386e-02 m
```

The small body-pose manufactured run also has some per-`z` variation, but it is
much closer to the reference:

```text
small OPALX body-pose per-z R56 range:
  min -4.334635442043e-02 m
  max -4.310741860762e-02 m
```

Particle-level OPALX-minus-OPAL comparison at the exit monitor for the fresh
10000 run:

```text
OPALX selected monitor: Step#0, SPOS = 9.250753364361e+00 m
OPAL  selected monitor: Step#0, SPOS = 9.255537067310e+00 m

rms(d_x)  = 3.293336007387e-05 m
rms(d_px) = 1.286926119514e-01
rms(d_y)  = 4.875653728863e-09 m
rms(d_py) = 2.703159808557e-06
rms(d_z)  = 7.482067864036e-04 m
```

A linear fit of particle-level `d_px = px_OPALX - px_OPAL` against standardized
initial `(x, y, z, delta)` gives `R^2 ~= 0.9706`, dominated by initial `z`:

```text
stdcoef z0    = +1.185505e-01
stdcoef delta = +1.196629e-02
stdcoef x0    = -1.937755e-04
stdcoef y0    ~= 0
```

Current interpretation:

- stale dense `.stat` output explains the old 10000 `By_ref` interval-shift
  table;
- it does not explain the larger horizontal/R56 differences;
- the 10000 grid exposes longitudinal-coordinate/monitor-plane dependence that
  the 125-particle diagnostic undersamples;
- next useful physics check is to separate true transport from temporal-monitor
  crossing conventions by comparing at matched monitor plane definitions and/or
  using a narrower initial `z` grid while keeping the 10000 transverse sampling.

### Beam-only comparison direction

User proposed not using temporal monitors for the OPALX-vs-OPAL comparison at
all.  The comparison should be based solely on beam quantities in the `.stat`
files:

```text
energy
dE
rms_x, rms_y, rms_s
rms_px, rms_py, rms_ps
emit_x, emit_y, emit_s
```

`compare_opal_opalx.py` now has a `--beam-only` mode for this.  In that mode it:

- writes `stat_summary.csv` and `summary.md` for only the beam quantities above;
- does not compute/write `by_field_intervals.csv`;
- does not add `By_ref`, reference-particle, temporal-monitor, or `R56`
  sections;
- writes the usual RMS position, RMS momentum, emittance, and energy plots plus
  a new `stat_dE_opalx_opal.png` plot.

Command run for the fresh 10000 `ELEMEDGE` stats:

```sh
/Users/adelmann/.venv-h6/bin/python \
  /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/compare_opal_opalx.py \
  --beam-only \
  --opalx-stat /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge/test-chicane-distribution-2.stat \
  --opal-stat /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opal_2022/test-chicane-distribution-1_opal.stat \
  --out-dir /Users/adelmann/git/opalx/sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/compare_elemedge_beam_only \
  --plots
```

Key fresh beam-only result:

```text
samples = 630
energy rms_diff = 2.176398868634e-10 MeV
dE rms_diff     = 1.747047e-11 MeV
rms_x rms_diff  = 4.647197478335e-06 m
rms_y rms_diff  = 2.095704012947e-09 m
rms_s rms_diff  = 1.640814146860e-07 m
rms_px rms_diff = 3.133895752589e-03
rms_py rms_diff = 4.837914516374e-07
rms_ps rms_diff = 1.249646840120e-06
emit_x rms_diff = 1.378381644701e-06
emit_y rms_diff = 5.201054430983e-12
emit_s rms_diff = 4.182755660300e-07
```

Interpretation with the beam-only target:

- energy and `dE` agree at roundoff/noise level;
- vertical RMS and vertical emittance agree closely;
- remaining differences are mainly horizontal RMS/momentum/emittance and some
  longitudinal RMS/emittance;
- monitor/R56 behavior should be treated as a separate diagnostic, not as the
  primary OPALX-vs-OPAL agreement criterion.

### Placement/interface follow-up

The phrase "driven by initial z" was based on a particle-level diagnostic for
the fresh 10000 `ELEMEDGE` run.  For each particle, compare
`px_OPALX - px_OPAL` at the exit monitor and regress that difference against
standardized initial `(x, y, z, delta)`.  The fit gives:

```text
R^2 = 0.970623
intercept = -4.398443559e-02
stdcoef x0    = -1.937755335e-04
stdcoef y0    ~= 0
stdcoef z0    = +1.185505439e-01
stdcoef delta = +1.196629152e-02
```

Mean `px_OPALX - px_OPAL` changes strongly by initial longitudinal slice:

```text
z0=-1.000000000e-03: mean_dpx=-2.286316101e-01
z0=-7.777777778e-04: mean_dpx=-2.148770443e-01
z0=-5.555555556e-04: mean_dpx=-1.485790116e-01
z0=-3.333333333e-04: mean_dpx=-6.715866092e-02
z0=-1.111111111e-04: mean_dpx=-5.085778106e-02
z0=+1.111111111e-04: mean_dpx=-3.447704695e-02
z0=+3.333333333e-04: mean_dpx=+1.975013997e-02
z0=+5.555555556e-04: mean_dpx=+4.554964777e-02
z0=+7.777777778e-04: mean_dpx=+7.714853708e-02
z0=+1.000000000e-03: mean_dpx=+1.622884742e-01
```

This does not prove a physical cause by itself.  It says the horizontal
momentum disagreement is phase-dependent through the bunch length: different
longitudinal slices see different horizontal kicks/edge timing.

The user's hint about `By` spikes at element interfaces is consistent with a
remaining element-placement or active-field-interval issue.  In the fresh sparse
10000 `.stat` files the sampled `By_ref` transitions line up almost exactly,
but the generated element position files still show millimeter-scale
differences between OPALX and OPAL `BEGIN`/`END` body/field markers.  Examples
from the generated `*_ElementPositions.txt` files:

```text
marker      OPALX first/horiz coord    OPAL first/horiz coord    first-coordinate diff
BEGIN B1    0.250000000 / 0.000000000  0.245000000 / 0.000000000  +5.000 mm
END B1      1.248750260 / 0.049979169  1.255101141 / 0.050613575  -6.351 mm
BEGIN B2    2.740841802 / 0.199687685  2.737242363 / 0.199326536  +3.599 mm
END B2      3.739592063 / 0.249666854  3.747350132 / 0.249807667  -7.758 mm
BEGIN B3    5.739175274 / 0.249666854  5.736933345 / 0.249804873  +2.242 mm
END B3      6.737925535 / 0.199687685  6.747034486 / 0.199191298  -9.109 mm
BEGIN B4    8.230017077 / 0.049979169  8.229175708 / 0.050478336  +0.841 mm
END B4      9.228767337 /-0.000000000  9.239283477 /-0.000002794 -10.516 mm
```

Some nominal edge markers agree better than these body markers, so the next
investigation should distinguish:

- nominal edge placement (`ENTRY EDGE`/`EXIT EDGE`);
- active field interval (`FIELD BEGIN`/`FIELD END`);
- generated body `BEGIN`/`END` markers;
- what the tracker actually uses for applying hard-edge `By`.

This is now the leading hypothesis for the remaining horizontal RMS/momentum
and emittance differences.

### Reproducible marker comparison and B1 focus

Added:

```text
sandbox/test-chicane-distribution-1/compare_element_markers.py
```

Purpose:

- parse generated OPALX and OPAL `*_ElementPositions.txt` files;
- compare common non-`MID` markers;
- write CSV and Markdown tables;
- support `--element B1` to focus on one bend.

Commands run:

```sh
/Users/adelmann/.venv-h6/bin/python compare_element_markers.py
/Users/adelmann/.venv-h6/bin/python compare_element_markers.py --element B1
```

Outputs:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/marker_compare/element_marker_comparison.csv
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/marker_compare/element_marker_comparison.md
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/marker_compare/element_marker_comparison_B1.csv
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/marker_compare/element_marker_comparison_B1.md
```

B1 marker comparison:

```text
ENTRY EDGE B1: OPALX - OPAL =  0.000000 mm
BEGIN B1:      OPALX - OPAL = +5.000000 mm
END B1:        OPALX - OPAL = -6.350881 mm in first coordinate,
                                 -0.634405 mm in second coordinate
EXIT EDGE B1:  OPALX - OPAL =  0.000000 mm
```

Important code-path finding:

- In `src/Elements/OpalBeamline.cpp`, bend `ENTRY EDGE`/`EXIT EDGE` markers are
  written from `getNominalEntryTransform()` / `getNominalExitTransform()`.
- Bend `BEGIN`/`END` markers are written from
  `placedElement.getNominalBodyTransform().transformFrom(designPath.front/back())`.
- `designPath` comes from `BendBase::getDesignPath()`, which samples
  `getTransform(s)` from `getEntrance()` to `getExit()`.
- For `SBEND`, `PlanarArcGeometry` uses the OPAL chord-to-arc normalization
  length from `OpalSBend::update()`: `L_norm = L * (theta/2)/sin(theta/2)`.
- `BendBase::updateFieldSupportExtent()` sets
  `fieldBegin_m = -entryFringe`, `fieldEnd_m = bodyLength + exitFringe`, with
  the default OPAL 1D profile fringe half-width.

For the B1 case with `GAP = 1e-3`, the generated files show:

```text
OPALX ENTRY EDGE B1 = 0.250000000 m
OPALX BEGIN B1      = 0.250000000 m
OPALX FIELD BEGIN B1= 0.245000000 m

OPAL ENTRY EDGE B1  = 0.250000000 m
OPAL BEGIN B1       = 0.245000000 m
```

So the apparent B1 `BEGIN` mismatch is at least partly a marker-definition
difference:

- OPAL's `BEGIN` appears to mean active field/support begin;
- OPALX's `BEGIN` means design-path begin in the bend body frame;
- OPALX writes the active support begin separately as `FIELD BEGIN`.

This explains the +5 mm B1 `BEGIN` discrepancy.  It does not by itself explain
all horizontal beam differences.  The next diagnostic should compare OPALX
`FIELD BEGIN/FIELD END` against OPAL `BEGIN/END` for all bends and check
whether the tracker uses those intervals consistently with old OPAL's hard-edge
and fringe profile semantics.

### Active field interval comparison

Extended `compare_element_markers.py` with:

```sh
--field-vs-opal-body
```

This compares OPALX `FIELD BEGIN/FIELD END` markers against OPAL `BEGIN/END`
markers for bends.  Commands run:

```sh
/Users/adelmann/.venv-h6/bin/python compare_element_markers.py --field-vs-opal-body
/Users/adelmann/.venv-h6/bin/python compare_element_markers.py --field-vs-opal-body --element B1
```

Outputs:

```text
marker_compare/element_marker_comparison_field_vs_opal_body.csv
marker_compare/element_marker_comparison_field_vs_opal_body.md
marker_compare/element_marker_comparison_field_vs_opal_body_B1.csv
marker_compare/element_marker_comparison_field_vs_opal_body_B1.md
```

Key result:

```text
OPALX FIELD BEGIN B1 vs OPAL BEGIN B1: 0.000000 mm
```

So the B1 entrance support agrees exactly if OPALX `FIELD BEGIN` is compared
to OPAL `BEGIN`.

However, raw `FIELD END` vs OPAL `END` is not directly interpretable for bends:

```text
OPALX FIELD END B1 = (1.255416788, 0.000000000)
OPAL  END B1       = (1.255101141, 0.050613575)
```

OPALX writes `FIELD END` by transforming local `(0, 0, fieldEnd)` with the
nominal entry transform.  For a curved bend this is an entry-axis display point,
not the curved reference-path endpoint.  The tracker path is different:

- `BendBase::getFieldCSTrafoLab2Local()` returns the nominal entry transform;
- `OpalBeamline::transformToLocal()` then calls
  `BendBase::transformFieldFrameToLocal()`;
- `BendBase::transformFieldFrameToLocal()` calls
  `convertEntryCartesianToFieldLocal()`;
- `BendBase::isInside()` checks the resulting local field coordinate against
  `fieldBegin_m <= z < fieldEnd_m`.

Thus the text marker for OPALX `FIELD END` is a poor geometry visualizer for a
curved bend.  It does not prove the tracker uses a straight field interval.

For later bends, OPALX field starts compared to OPAL bend starts still show
accumulated offsets:

```text
B1 FIELD BEGIN vs OPAL BEGIN:  0.000000 mm
B2 FIELD BEGIN vs OPAL BEGIN: -1.382489 mm
B3 FIELD BEGIN vs OPAL BEGIN: -2.761522 mm
B4 FIELD BEGIN vs OPAL BEGIN: -4.133652 mm
```

These offsets are now the most concrete placement discrepancy.  They are close
to accumulated chord/arc/reference-line convention offsets and should be
checked before changing tracking.

### GAP diagnostics

Added `run_gap0_elemedge_compare.py` as a reusable configurable-GAP diagnostic.
It now supports:

```sh
--gap <value>
--statdumpfreq <integer>
```

The default output directory is derived from the gap and stat dump frequency,
for example:

```text
comparison/dt_scan/finite_10000/dt_5p0em12_gap1p0em09_sdf20/
```

For future long beam-only scans, prefer `STATDUMPFREQ=20` unless the task is
specifically to inspect samples at element interfaces.  `STATDUMPFREQ=10` is
also acceptable, but `20` keeps the generated stat files and plots smaller
without changing tracking.

`GAP=0` is not a valid no-support diagnostic for old OPAL 2022.1.  It triggers
the unscaled default `1DPROFILE1` support behavior instead of suppressing the
support.  The beam-only comparison becomes much worse than the original
`GAP=1e-3` comparison:

```text
quantity  rms_diff GAP=1e-3  rms_diff GAP=0   ratio
energy    2.176399e-10       9.967616e-11     0.458
dE        1.747047e-11       1.000344e-11     0.573
rms_x     4.647197e-06       2.227305e-04    47.928
rms_y     2.095704e-09       6.572536e-08    31.362
rms_s     1.640814e-07       5.507249e-07     3.356
rms_px    3.133896e-03       8.748455e-02    27.916
rms_py    4.837915e-07       3.326040e-05    68.749
rms_ps    1.249647e-06       1.028082e-05     8.227
emit_x    1.378382e-06       6.915888e-05    50.174
emit_y    5.201054e-12       1.505977e-10    28.955
emit_s    4.182756e-07       1.412462e-06     3.377
```

`GAP=1e-12` was also attempted with `STATDUMPFREQ=20`.  OPALX wrote a stat
file, but old OPAL did not produce a usable stat file before the diagnostic
wrapper had to be stopped.  Treat this as a pathological near-zero-gap case,
not as a physics result.

`GAP=1e-9` with `STATDUMPFREQ=20` completed in old OPAL, but old OPAL deleted
the whole beam:

```text
ParallelTTracker > * Deleted 4000 particles, remaining 6000 particles
ParallelTTracker > * Deleted 6000 particles, remaining 0 particles
```

The OPAL stat file ends with `numParticles=0` and NaNs, so this case is also
invalid for beam RMS comparisons.  It is still useful for element-map geometry.

For `GAP=1e-9`, OPALX field starts compared to OPAL bend starts are:

```text
B1 FIELD BEGIN vs OPAL BEGIN:  0.000000 mm
B2 FIELD BEGIN vs OPAL BEGIN: -0.885434 mm, -0.088840 mm
B3 FIELD BEGIN vs OPAL BEGIN: -1.775313 mm, -0.088840 mm
B4 FIELD BEGIN vs OPAL BEGIN: -2.660746 mm,  0.000000 mm
```

The corresponding `ENTRY EDGE`/`BEGIN`/`EXIT EDGE` marker offsets are almost
the same, so the small positive-gap element map shows an accumulated reference
placement offset of about `0.89 mm` per bend after B1.

### Curved bend field-extent marker fix

Implemented the visualization fix discussed after the GAP diagnostics:

- Added `BendBase::getDesignPathPoint(double z)`.
- The helper accepts the bend field-local longitudinal coordinate, where
  `z=0` is body entry and `z=L` is body exit.
- Points inside the body are evaluated on the curved design path.
- Points outside the body use tangent extrapolation from the entry or exit
  frame, which represents fringe/support regions without changing tracking.
- `OpalBeamline::save3DLattice()` now writes bend `FIELD BEGIN`/`FIELD END`
  using that curved/tangent path visualization.
- Non-bend field extent markers still use the existing straight field-frame
  transform.

This is a marker-output/diagnostic change only.  It does not change
`fieldBegin_m`, `fieldEnd_m`, `BendBase::isInside()`, field evaluation, or
particle tracking.

Focused validation:

```sh
cmake --build build --target TestOpalBeamlinePlacement -j
ctest --test-dir build --output-on-failure -R TestOpalBeamlinePlacement
```

Result:

```text
1/1 Test #71: TestOpalBeamlinePlacement ........ Passed
```

Then rebuilt the actual executable and regenerated the 10000-particle
`ELEMEDGE` OPALX marker/stat output:

```sh
cmake --build build -j
/Users/adelmann/git/opalx/build/src/opalx test-chicane-distribution-2.in
```

The first attempt accidentally used the stale executable because `cmake
--build build --target opalx` only rebuilt the library target.  The executable
target is `opalx_exe`, and the full build updated `build/src/opalx`.

After regenerating with the updated executable, the bend field ends are on the
curved/tangent path.  Example OPALX marker rows:

```text
"FIELD END: B1"  1.253725281  0.05047833635  0
"FIELD END: B2"  3.744592063  0.2496668541   0
"FIELD END: B3"  6.742900556  0.1991885178   0
"FIELD END: B4"  9.233767337 -7.042977312e-16 0
```

The corrected field-vs-OPAL-body comparison is now:

```text
B1 FIELD BEGIN vs OPAL BEGIN:  0.000000 mm,  0.000000 mm
B1 FIELD END   vs OPAL END:   -1.375860 mm, -0.135238 mm
B2 FIELD BEGIN vs OPAL BEGIN: -1.375582 mm, -0.138019 mm
B2 FIELD END   vs OPAL END:   -2.758069 mm, -0.140813 mm
B3 FIELD BEGIN vs OPAL BEGIN: -2.758071 mm, -0.138019 mm
B3 FIELD END   vs OPAL END:   -4.133930 mm, -0.002780 mm
B4 FIELD BEGIN vs OPAL BEGIN: -4.133652 mm,  0.000000 mm
B4 FIELD END   vs OPAL END:   -5.516140 mm,  0.002794 mm
```

This removes the previous false `~50 mm` transverse discrepancy in the
`FIELD END` marker comparison.  The remaining geometry discrepancy is now a
clear accumulated longitudinal/reference placement offset of about `1.38 mm`
per bend support interval.

The OPALX stat file was rerun to normal process completion after an early stop
temporarily left it partial.  The beam-only comparison is back to 630 samples:

```text
rms_x  rms_diff = 4.647197e-06
rms_px rms_diff = 3.133896e-03
emit_x rms_diff = 1.378382e-06
energy rms_diff = 2.178407e-10
dE     rms_diff = 1.747829e-11
```

Recommended next discriminating checks:

1. Do not continue with `GAP=0`, `GAP=1e-12`, or `GAP=1e-9` as beam RMS
   comparisons in old OPAL.  They are diagnostics of old-OPAL support/profile
   behavior, not valid beam-code comparisons.
2. Compare scalar active intervals from the OPALX tracker
   (`fieldBegin_m`, `fieldEnd_m`, `getReferencePathLength()`) against old
   OPAL's effective Bend2D intervals, not against raw displayed marker
   coordinates.
3. Inspect the reference placement accumulation from B1 to B4.  The current
   evidence points to an accumulated chord/arc/reference-line convention
   mismatch of order millimeters, while the original `GAP=1e-3` beam
   disagreement remains dominated by horizontal RMS/momentum/emittance.
4. Next inspect the remaining `~1.38 mm` accumulated offset: compare OPALX
   body/reference length convention against old OPAL's effective Bend2D
   support/body convention bend-by-bend.

## 2026-06-10: Remaining sigma_x discrepancy / placement attack

The final `sigma_x` discrepancy in the completed 10000-particle beam-only
comparison remains about 20--25%:

```text
rms_x final OPALX = 9.166109e-05
rms_x final OPAL  = 7.418512e-05
diff              = 1.747598e-05
```

The marker evidence still points to placement, not monitor comparison:

- OPALX `FIELD BEGIN` now coincides with old OPAL `BEGIN` for B1.
- OPALX `FIELD END` is now written on the curved/tangent bend path.
- The remaining difference accumulates by about `1.38 mm` per bend support
  interval.

Old OPAL source inspection:

- `/Users/adelmann/git/old-opal/src/Classic/AbsBeamline/Bend2D.cpp`
  `Bend2D::initialise()` calls `setElementLength(endField_m - startField_m)`
  after `setupBendGeometry()`.
- `setupBendGeometry()` rescales the default `1DPROFILE1-DEFAULT` Enge support
  when `GAP` differs from the map default. For the chicane's `GAP=1e-3`, the
  default `+/-0.1 m` support is scaled to `+/-0.005 m`.
- Old OPAL then shifts/searches the Enge origins in
  `findBendEffectiveLength()` to match the requested bend angle.  This affects
  the displayed/effective support convention used by the old `BEGIN/END`
  markers and by the old placement cursor.

Rejected experiment:

- Tried changing OPALX idealized compatibility placement to advance the
  reference cursor by the bend chord length instead of `getArcLength()`.
- Unit-level direct object construction was not a valid regression path for
  this because it bypasses parser/compatibility semantics.
- A real chicane run in
  `comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge_chordcursor`
  showed the change was not a small correction: the trajectory progressed past
  the previous complete-run extent (`s = 18.9 m`, 686 stat rows) and had to be
  killed.  The code change was reverted.

Current verified code state after reverting that experiment:

```sh
cmake --build build --target TestOpalBeamlinePlacement -j
ctest --test-dir build --output-on-failure -R TestOpalBeamlinePlacement
```

Result:

```text
1/1 Test #71: TestOpalBeamlinePlacement ........ Passed
```

Next useful direction:

1. Do not make a blind cursor-length change.  The OPALX compatibility placement
   and tracking-time termination are coupled through the generated element
   transforms.
2. Instrument old OPAL effective Bend2D quantities from its output/logs or a
   minimal old-OPAL debug build: `startField_m`, `endField_m`,
   `elementEdge_m`, final shifted Enge origin delta, `getElementLength()`, and
   the `arcLength` used in `OpalBeamline::compute3DLattice()`.
3. Add the same one-line diagnostic to OPALX around
   `OpalBeamline::compileCompatibilityPlacement()` for B1--B4:
   `ELEMEDGE`, `thisLength`, `getArcLength()`, field support, computed
   `arcLength`, and the resulting `beginThis3D`.
4. Once the exact scalar convention is reproduced, adjust only the compatibility
   placement bridge or the comparison deck.  Do not alter the bend field chart
   or curved marker writer without a separate physics reason.

## 2026-06-10: DT-matched scalar diagnostics and old-gap test

Ran short 125-particle diagnostics at the same timestep as the 10k comparison
(`DT = 5e-12`) because old OPAL `Bend2D::findBendEffectiveLength()` uses the
reference pusher timestep.

Diagnostic directories:

```text
sandbox/test-chicane-distribution-1/diagnostics/old_opal_info_dt5
sandbox/test-chicane-distribution-1/diagnostics/opalx_elemedge_info_dt5
```

Generated reusable diagnostic:

```sh
/Users/adelmann/.venv-h6/bin/python \
  sandbox/test-chicane-distribution-1/diagnose_bend_intervals.py \
  --old-log sandbox/test-chicane-distribution-1/diagnostics/old_opal_info_dt5/old_opal_info.log \
  --old-markers sandbox/test-chicane-distribution-1/diagnostics/old_opal_info_dt5/data/test-chicane-distribution-1_opal_ElementPositions.txt \
  --opalx-markers sandbox/test-chicane-distribution-1/diagnostics/opalx_elemedge_info_dt5/data/test-chicane-distribution-2_ElementPositions.txt \
  --out-dir sandbox/test-chicane-distribution-1/diagnostics/bend_interval_compare_dt5
```

Key result from `bend_interfaces.md`:

```text
B1->D12   old projected gap = 0.964497 mm, OPALX = 0
B2->D23   old projected gap = 0.964497 mm, OPALX = 0
B3->D34   old projected gap = 0.964497 mm, OPALX = 0
B4->MEXIT old projected gap = 0.964496 mm, OPALX = 0
```

So old OPAL starts the following element about `0.965 mm` downstream of the bend
exit edge, while OPALX's idealized ELEMEDGE compatibility deck starts the next
element exactly at the bend exit edge.  This explains part of the interface
marker mismatch, but it is not the dominant final `sigma_x` discrepancy.

Direct old-gap hypothesis test:

- Created
  `comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge_oldgap/test-chicane-distribution-2.in`.
- Added cumulative `ELEMEDGE` offsets using
  `OLD_OPAL_GAP = 9.644967e-4 m` after each bend.
- Ran OPALX to normal completion.
- Compared beam-only stats in
  `comparison/dt_scan/finite_10000/dt_5p0em12/compare_elemedge_oldgap_beam_only`.

Result:

```text
rms_x final diff:  1.976340e-05  (worse than baseline 1.747598e-05)
rms_x rms diff:    8.724355e-06  (worse than baseline 4.647197e-06)
rms_px final diff: 2.203006e-03  (better than baseline ~1.35e-02)
rms_px rms diff:   4.866209e-03  (worse than baseline 3.133896e-03)
```

Marker check for old-gap deck:

```text
B2/B3/B4 FIELD BEGIN vs old OPAL BEGIN: within ~1--4 microns
B1/B2/B3/B4 FIELD END vs old OPAL END: still about 1.38 mm short
```

Conclusion:

- The old OPAL downstream interface gap is real and timestep-dependent.
- Adding that gap to OPALX ELEMEDGE placement does not explain or fix final
  `sigma_x`; it worsens beam-size agreement while improving only the final
  momentum-spread-like horizontal momentum endpoint.
- The remaining `sigma_x` issue is more likely from bend field/edge/fringe
  integrated optics, not simply from the post-bend element placement gap.

Next useful checks:

1. Compare bend-local field sampling/kicks along the reference particle, not
   just element placement.  Focus on B1 and B4 where final horizontal beam-size
   sensitivity is largest.
2. Generate side-by-side `B_y(s)` and integrated bend kick curves for OPALX's
   analytic fringe model versus old OPAL's effective default map at `DT=5e-12`.
3. Keep the old-gap deck as a diagnostic only; do not turn it into an OPALX code
   change.

## 2026-06-12 reference-field interface plots

User asked to step back and inspect the field spikes before further code
changes.

Added reusable diagnostic:

```sh
/Users/adelmann/.venv-h6/bin/python \
  sandbox/test-chicane-distribution-1/plot_reference_field_interfaces.py \
  --field Bz_ref

/Users/adelmann/.venv-h6/bin/python \
  sandbox/test-chicane-distribution-1/plot_reference_field_interfaces.py \
  --field By_ref
```

Outputs:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/reference_field_interface_plots/
```

Observation from stat-file reference fields:

```text
Bz_ref max abs OPALX-OPAL difference near bend interfaces: ~1e-21 T
By_ref max abs OPALX-OPAL difference near bend interfaces:
  B2 entry 9.956386e-11 T
  B2 exit  1.251156e-10 T
  B3 entry 3.239108e-08 T
  B3 exit  4.540363e-11 T
  B4 entry 3.700683e-10 T
  B4 exit  1.189152e-10 T
```

Interpretation:

- The `Bz_ref` spike visible when plotted on an auto-scaled axis is roundoff at
  essentially zero field.
- The relevant bend component in these stats is `By_ref`; it overlays between
  OPALX and OPAL, and the absolute edge differences are tiny compared with the
  bend field (~0.33 T).
- The previous visual spike is therefore most likely a relative-error artifact
  at a sharp field edge, not direct evidence of a millimeter-scale field
  placement mismatch in the stat-file field intervals.
- Still unresolved: the marker geometry and the actual integrated optics can
  differ even when the stat-sampled reference field interval overlays.  The next
  diagnostic should compare bend-local field/kick integrals along a prescribed
  common trajectory rather than only comparing each code's own stat reference
  particle sample.

## 2026-06-12 common-grid reference-field integrals

Added reusable diagnostic:

```sh
/Users/adelmann/.venv-h6/bin/python \
  sandbox/test-chicane-distribution-1/compare_reference_field_integrals.py \
  --field By_ref
```

Output directory:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/reference_field_integrals/
```

The script detects the four nonzero `By_ref` bend intervals from the stat files,
interpolates OPALX and old OPAL onto the same dense `s` grid for each bend, and
compares cumulative trapezoidal `int B_y ds`.

Result:

```text
B1 integral diff:  3.315125952e-13 T m  (1.00e-12 relative)
B2 integral diff:  9.388601008e-13 T m  (-2.84e-12 relative)
B3 integral diff:  2.444530134e-10 T m  (-7.41e-10 relative)
B4 integral diff: -2.164934898e-15 T m  (-6.56e-15 relative)
```

The worst case is B3, but even there the relative integrated-field difference is
`~7e-10`.  Therefore the stat-file reference `By_ref` spikes/edge differences do
not explain the remaining final `sigma_x` discrepancy.

Important caveat:

- This still compares each code's exported stat reference fields, then places
  those samples on a common `s` grid.
- It does not yet force OPAL and OPALX to evaluate the field on an externally
  prescribed identical 3D trajectory.
- If placement remains suspect, the next stronger diagnostic is a direct
  prescribed-trajectory field sampler or trajectory/kick comparison through the
  element-local coordinates.

## 2026-06-12 longitudinal slice comparison

Added reusable diagnostic:

```sh
/Users/adelmann/.venv-h6/bin/python \
  sandbox/test-chicane-distribution-1/compare_longitudinal_slices.py
```

Output directory:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/longitudinal_slice_compare/
```

The script groups the fresh 10000-particle OPALX/OPAL exit monitor H5 data by
initial `z0`.  Each longitudinal slice contains the full transverse and
momentum grid.  The selected monitor samples match the existing R56 diagnostic
convention:

```text
OPALX selected step: Step#0, SPOS = 9.250753364361 m
OPAL  selected step: Step#0, SPOS = 9.255537067310 m
```

Main result:

```text
largest |mean(OPALX-OPAL px)|:
  z0 = -1.0 mm, mean Delta px = -2.286316100580e-01

linear fit mean(Delta px) vs z0:
  slope     = 1.857330071435e+02 1/m
  intercept = -4.398443558713e-02
  R^2       = 0.974542
```

The horizontal kick difference is therefore strongly and almost linearly
dependent on initial longitudinal slice.

Horizontal `sigma_x` variance decomposition at the exit monitor:

```text
OPALX:
  total sigma_x          = 8.954720529e-05 m
  within-slice sigma     = 6.069065463e-05 m
  slice-centroid sigma   = 6.584334755e-05 m
  centroid variance frac = 0.540654389

OPAL 2022.1:
  total sigma_x          = 7.392105877e-05 m
  within-slice sigma     = 6.031931605e-05 m
  slice-centroid sigma   = 4.273058671e-05 m
  centroid variance frac = 0.334149915
```

Interpretation:

- The average within-slice horizontal size is essentially the same in OPALX and
  OPAL (`~6.05e-5 m`).
- The final `sigma_x` discrepancy comes mainly from larger OPALX
  longitudinal-slice centroid spread, not from broader individual slices.
- This is consistent with the earlier particle-level regression showing that
  `px_OPALX - px_OPAL` is dominated by initial `z`.
- The next discriminating check should track where the slice-centroid shear is
  created along the line: compare the same `z0`-grouped centroids after B1, B2,
  B3, and B4, or add equivalent per-section monitors in both decks.  Do not
  return to reference-field integral checks for this issue.

## 2026-06-12 section longitudinal-slice localization

Added section-localization diagnostics under:

```text
sandbox/test-chicane-distribution-1/
```

Scripts:

```text
run_section_stop_compare.py
compare_section_stop_slices.py
run_section_terminal_monitor_compare.py
compare_section_monitor_slices.py
```

Also added `run_section_monitor_compare.py`, but that first attempt inserted
multiple monitors into one line and perturbed/extended the tracking.  Treat the
`opalx_elemedge_section_monitors/` and `opal_2022_section_monitors/` outputs as
invalid for physics conclusions.

### Raw ZSTOP section dumps

`run_section_stop_compare.py` truncates by `ZSTOP` without changing the lattice.
Those H5 dumps are useful as raw tracking-coordinate checkpoints, but they do
not use the same temporal-monitor local coordinate convention as the final
exit-monitor comparison.  They therefore did not show the same strong final
`z0`-dependent `px` shear and should not be used as the primary explanation for
the final `sigma_x` discrepancy.

Output:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/section_stop_slice_compare/
```

### Terminal temporal section monitors

`run_section_terminal_monitor_compare.py` makes separate truncated decks with
the original upstream lattice ending at one temporal monitor:

```text
after_b1, after_b2, after_b3, after_b4
```

This is the better diagnostic because the exported coordinates use the same
kind of temporal-monitor local frame as the final exit monitor, while avoiding
upstream inserted monitors.

Output:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/section_monitor_slice_compare/
```

Important summary:

```text
section   dpx slope [1/m]   R^2       sigma_x OPALX   sigma_x OPAL    centroid sigma OPALX   centroid sigma OPAL
after_b1   9.626202694e+01  0.848348  9.051549502e-05 9.061533863e-05 9.427077541e-06       9.512189147e-06
after_b2  -1.226961835e+02  0.816073  3.237417463e-04 3.243036311e-04 2.187291203e-05       2.733236025e-05
after_b3  -5.074452371e+01  0.933698  2.511173597e-04 2.535000978e-04 3.996537796e-05       4.462316082e-05
after_b4   7.952424722e+01  0.824036  7.525249562e-05 7.392105877e-05 4.535636421e-05       4.273058671e-05
```

Additional plots now generated by `compare_section_monitor_slices.py`:

```text
monitor_centroid_sigma_x.png
monitor_dpx_slope.png
slice_mean_x_by_section.png
slice_diff_mean_x_by_section.png
slice_mean_px_by_section.png
slice_diff_mean_px_by_section.png
slice_rms_x_by_section.png
slice_diff_rms_x_by_section.png
```

Interpretation so far:

- The OPALX-vs-OPAL horizontal momentum-centroid difference is strongly
  longitudinal-slice dependent at each terminal section monitor, and its sign
  flips through the chicane.
- The final 20% OPALX `sigma_x` excess is not from wider individual
  longitudinal slices.  The within-slice widths remain close; the dominant term
  is the spread of slice centroids.
- OPAL has the larger slice-centroid contribution after B2 and B3.  OPALX only
  becomes larger by the final bend/exit region.
- The terminal `after_b4` monitor is close to, but not numerically identical
  to, the original baseline final `MEXIT` monitor in OPALX:

```text
baseline MEXIT OPALX: Step#0 SPOS = 9.250753364361 m and Step#1 SPOS = 9.25000388 m
terminal MSEC OPALX:  Step#0 SPOS = 9.251669235 m
OPAL MEXIT/MSEC:      Step#0 SPOS = 9.255537067 m in both cases
```

The baseline final exit monitor remains the authoritative final discrepancy
for now:

```text
baseline MEXIT OPALX total sigma_x        = 8.954720529e-05 m
baseline MEXIT OPAL total sigma_x         = 7.392105877e-05 m
baseline MEXIT OPALX centroid sigma_x     = 6.584334755e-05 m
baseline MEXIT OPAL centroid sigma_x      = 4.273058671e-05 m
baseline MEXIT mean(Delta px) vs z0 slope = 1.857330071435e+02 1/m
```

A clean OPALX baseline rerun was attempted in:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000/dt_5p0em12/opalx_elemedge_exact_rerun/
```

It was killed with `pkill` for this diagnostic only.  Later inspection of the
SDDS header showed the first stat column is `t`, not `s`; the killed run had
not passed `ZSTOP`, it was just slow.  Do not use that partial rerun for physics
conclusions.

Recommended next step:

1. Treat the original final `MEXIT` slice comparison as authoritative for the
   20% `sigma_x` discrepancy.
2. Inspect why OPALX final `MEXIT` writes two full-particle steps and why the
   nominally equivalent terminal `MSEC` writes one step at a different `SPOS`.
3. Then focus on the last-bend/exit-monitor placement and temporal-monitor
   crossing logic, not on the B1-B3 body fields.

## 2026-06-12 current-build DT and final-monitor check

Added scripts:

```text
sandbox/test-chicane-distribution-1/compare_dt_monitor_slices.py
sandbox/test-chicane-distribution-1/run_reduced_dt_monitor_compare.py
```

`compare_dt_monitor_slices.py` summarizes existing or newly produced
exit-monitor H5 files by timestep.  It records monitor step metadata, `SPOS`,
the fitted slope of slice-mean `Delta px` versus initial `z0`, and the
horizontal variance decomposition.

`run_reduced_dt_monitor_compare.py` can generate reduced-grid or full-grid
OPALX/OPAL reruns from the 10000-particle baseline decks.  It now supports
`--transverse-mode all` and `--skip-opal`.

### Existing 125-particle finite DT scan

Output:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite/dt_monitor_slice_compare/
```

The complete existing 125-particle ELEMEDGE/OPAL runs show that the
monitor-local slice-kick diagnostic is timestep-sensitive:

```text
DT = 2e-12:  dpx slope = 8.495468905e+01 1/m
DT = 5e-12:  dpx slope = 1.898318964e+02 1/m
DT = 1e-11:  dpx slope = 2.062170763e+02 1/m
```

The 125-particle sample is too sparse for the final `sigma_x` question, but it
does show that temporal-monitor crossing/drift is coupled to `DT`.

### Reduced 2500-particle current-build DT scan

Output:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_2500/dt_monitor_slice_compare/
```

This used a 2500-particle subset preserving all 10 `z0` and all 10 `delta`
slices.  It is useful for timestep sensitivity, but the transverse subset is
not centered and should not be used as a final physics result.

Result:

```text
DT = 2e-12: dpx slope = 3.642769548e+01 1/m, sigma_x(OPALX-OPAL)/OPAL = 1.805761%
DT = 5e-12: dpx slope = 7.955888371e+01 1/m, sigma_x(OPALX-OPAL)/OPAL = 1.765055%
```

The slope decreases with smaller timestep, while total `sigma_x` is already
close for the current build.

### Fresh full 10000-particle current-build OPALX run

Ran current OPALX only, reusing the existing OPAL 2022.1 `DT=5e-12` monitor as
reference:

```sh
/Users/adelmann/.venv-h6/bin/python \
  sandbox/test-chicane-distribution-1/run_reduced_dt_monitor_compare.py \
  --out-dir sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000_current \
  --dt 5.0e-12 \
  --transverse-mode all \
  --skip-opal \
  --statdumpfreq 20
```

Then generated fresh longitudinal-slice plots:

```text
sandbox/test-chicane-distribution-1/comparison/dt_scan/finite_10000_current/dt_5p0em12/longitudinal_slice_compare/
```

Fresh current-build full-grid result:

```text
OPALX selected step: Step#0, SPOS = 9.251669234960 m
OPAL  selected step: Step#0, SPOS = 9.255537067310 m

dpx slope vs z0 = 7.952424721976e+01 1/m, R^2 = 0.824036

OPALX total sigma_x        = 7.525249562012e-05 m
OPAL  total sigma_x        = 7.392105876699e-05 m
relative difference        = 1.801160421 %

OPALX within-slice sigma_x = 6.004780031643e-05 m
OPAL  within-slice sigma_x = 6.031931604890e-05 m

OPALX centroid sigma_x     = 4.535636420850e-05 m
OPAL  centroid sigma_x     = 4.273058671053e-05 m
centroid relative diff     = 6.144960086 %
```

Important conclusion:

- The earlier 10000-particle `20%` `sigma_x` discrepancy came from stale OPALX
  H5 output under `finite_10000/dt_5p0em12/opalx_elemedge/`.
- With the current executable and the same full 10000 distribution at
  `DT=5e-12`, the discrepancy is about `1.8%`, matching the terminal
  `after_b4` section-monitor result.
- The remaining issue is no longer the large 20% discrepancy.  What remains is
  a smaller monitor-local slice-centroid difference and clear timestep
  sensitivity of the `Delta px` slice slope.

## Pending workflow discussion: stalled executables

Before continuing with more long OPALX/OPAL comparison runs, explicitly discuss
the stalled-executable problem and change the modus operandi.  Stale or stalled
executables have caused misleading results in this investigation, including the
apparent 10000-particle `20%` `sigma_x` discrepancy from old OPALX H5 output.

Proposed points to settle before the next run:

- Always record the exact executable path, build timestamp or commit, and run
  start time next to generated H5/stat outputs.
- Treat unexpectedly slow or silent runs as a workflow event, not just as long
  jobs: check stat progress, executable identity, and output freshness before
  interpreting results.
- Prefer scripts that create fresh output directories with explicit run labels
  such as `current`, `rerun`, or the commit hash; avoid reusing ambiguous
  directories.
- Add a timeout/progress policy for sandbox comparisons so stale processes are
  killed deliberately and documented.
- Do not draw physics conclusions from H5/stat files until their producer run
  has been verified as fresh and complete.
