# Drift Experiment Task State

Last updated: 2026-07-04

## 2026-07-04 TestParticleOrbit witnesses and full c0 charge

Goal:

- Move the current BeamBeam witness deck from the transverse pair-4 smoke
  particles toward the `TestParticleOrbit.dat` model.

Changes:

- Updated
  `sandbox/Drift-Experiment/withness_t0_4ps_active_1600fs_release_shifted_c0/spacecharge_drift_withness.in`:
  - `primary_charge_scale = 1.0`;
  - `n_witness_per_species = 6`;
  - `witness_kinetic_energy = 0.00051099917 GeV`;
  - c1/c2 now read separated electron/positron files from
    `../../track12particles/track12_electrons.fromfile` and
    `../../track12particles/track12_positrons.fromfile`.
- Generated separated witness files under `sandbox/track12particles` from the
  workspace copy `sandbox/TestParticleOrbit.dat`:
  - `track12_electrons.fromfile`;
  - `track12_positrons.fromfile`;
  - `track12_witness_metadata.csv`.

Notes:

- Direct reads from `/Users/adelmann/Desktop/TestParticleOrbit.dat` are blocked
  by macOS sandbox permissions in this Codex process.  The generated files use
  the existing repository copy `sandbox/TestParticleOrbit.dat`, which matches
  the already parsed track12 initial-condition data.
- All six e-/e+ pairs are injected at the same absolute witness T0 of `4 ps`.
  The c0 centroid is aligned with the IP for the central `ct = 0` pair; exact
  pair-wise slide timing for nonzero `ct_i` still requires one pair per deck or
  explicit source timing per pair.

2026-07-04 rerun/update:

- Running the full-charge deck exposed an incorrect c0 energy jump:
  c0 was initialized at `<EKIN> = 245 MeV` but the old tracker path applied a
  self-space-charge kick to the BeamBeam source container, producing `E ~ 386
  MeV` in the step log.
- Fixed the tracker so a BeamBeam source with witness containers still deposits
  fields for c1/c2 but does not apply those source fields back to container 0.
  This is done both in the active BeamBeam window and in the pre-entry step
  where the BeamBeam geometry has been detected but the tail is not yet one
  mesh cell inside the window.
- Rebuilt `build_openmp` successfully.
- Reran:
  `env OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads ../../../build_openmp/src/opalx --info 2 spacecharge_drift_withness.in > spacecharge_drift_withness.out 2> spacecharge_drift_withness.err`
  in
  `sandbox/Drift-Experiment/withness_t0_4ps_active_1600fs_release_shifted_c0`.
- Runtime from `/usr/bin/time`: `real 365.96 s`, `user 1900.59 s`, `sys 29.63 s`.
- Parsed stdout c0 energy:
  - 1641 c0 log lines;
  - first: step 0, `E=245.000 MeV`;
  - last: step 1640, `E=245.000 MeV`;
  - min/max: `245.000 / 245.000 MeV`.
- `spacecharge_drift_withness.err` contains only timing output; no OPALX error.
- Relevant tests:
  `ctest --test-dir build_openmp -R BeamBeam --output-on-failure`
  passed `TestBeamBeamDiagnosticsWriter` and `TestBeamBeam`.

## 2026-07-04 near-IP active-field then ballistic model

Goal:

- Implement the agreed first version of the c0-active / c1-c2-ballistic model:
  - c0 is the only source that deposits charge;
  - c1/c2 are passive witnesses that gather/sample c0 fields;
  - c0 is retired after the near-IP active-field window;
  - c1/c2 continue ballistically after c0 retirement.

Implementation notes:

- No C++ change was needed for this first test after inspection:
  - the BeamBeam self-field path scatters `bunch.getParticleContainer()`, i.e.
    c0/container 0;
  - configured witness containers are gathered passively by
    `gatherBeamBeamFieldsToWitnessContainers()`;
  - the witness gather samples fields but does not deposit witness charge.
- The deck was changed to use the real-bunch near-IP cutoff estimate rather than
  the point-charge `51 fs` estimate:
  - `near_ip_active_time = 8.0e-13 s`;
  - `primary_retire_time = witness_t0 + near_ip_active_time`;
  - `fine_stop_s = bb_ip_s + 2.5e-4`;
  - `fine_steps = 1020`;
  - `c0_source_mesh_diameter = 5.0e-4`;
  - `APERTURE = "CIRCLE(0.0005)"`, about `250 um` radius.
- The analyzer now writes:
  - `near_ip_field_cutoff.csv`;
  - `near_ip_field_cutoff.png`;
  and marks the c0 retirement time on `witness_timing_overview.png`.

Verification:

- OPALX run directory:
  `sandbox/Drift-Experiment/withness_near_ip_active_1pct`
- Command:
  `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
- Analyzer:
  `env MPLBACKEND=Agg MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py --run-dir sandbox/Drift-Experiment/withness_near_ip_active_1pct --aperture-radius-m 2.5e-4 --near-ip-active-time-s 8.0e-13`
- `~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/analyze_withness_timing.py`

Results:

- OPALX completed cleanly.
- BeamBeam diagnostics:
  - active with c1/c2 empty;
  - active with c1/c2 each containing one particle;
  - completed with `source_retirement_pending=TRUE`;
  - completed with `retired_bunches=1`, `source_active=FALSE`, and c1/c2 still
    active.
- Pair-4 timing remains aligned:
  - `witness_t0 = 13.342592708656 ps`;
  - first c1/c2 H5 dump `t = 13.348 ps`;
  - first dump is `5.407 fs` after `T0`;
  - c0 is at `s = 4.001621 mm`, `+1.621 um` past IP.
- Configured retirement:
  - `retire_time = 14.142592708656 ps`;
  - `retire_time - witness_t0 = 800 fs`.
- Final written c1/c2 H5 sample:
  - `t = 14.220 ps`;
  - c1 `y = +206.3 um`;
  - c2 `y = -206.3 um`;
  - `max_abs_y_over_aperture = 0.825`, so both remain inside the 250 um mesh.
- `near_ip_field_cutoff.csv` shows that after the first post-retirement H5 dump
  at `14.148 ps`, c1/c2 have `E_abs = 0`, consistent with ballistic
  propagation after c0 is retired.

Output files:

- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/timing_mesh_summary.csv`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/witness_kinematics_summary.csv`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/witness_field_samples.csv`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/near_ip_field_cutoff.csv`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/witness_timing_overview.png`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/c0_injection_mesh_xy.png`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/witness_transverse_motion.png`
- `sandbox/Drift-Experiment/withness_near_ip_active_1pct/near_ip_field_cutoff.png`

Remaining caveat:

- This is a timing/containment/passive-gather smoke test only.  With `NXY = 16`
  and 250 um radius, the transverse cell size is `31.25 um`, much larger than
  the `1.944 um` c0 transverse rms size, so quantitative field accuracy still
  requires a better mesh/resolution study.

## 2026-07-04 transverse witness step 1/2

Goal:

- Start the `spacecharge_drift_withness.in` BeamBeam study by doing only the
  first two agreed steps:
  - correct c1/c2 witness kinematics to transverse 313 keV motion;
  - keep the BeamBeam mesh as a compact c0 source field mesh, not a 30 cm
    witness-domain mesh.

Snapshot before the change:

- Directory:
  `sandbox/Drift-Experiment/snapshots/2026-07-04_before_transverse_witness`
- Contains the staged pair-4 timing deck, old FROMFILE witnesses, and the
  previous timing summary.
- The old FROMFILE witnesses had longitudinal normalized momentum
  `pz = 1.7320511804180809` (`gamma = 2`, about `511 keV` kinetic energy).

Changed files:

- `sandbox/Drift-Experiment/spacecharge_drift_withness.in`
- `sandbox/Drift-Experiment/make_transverse_witness_fromfiles.py`
- `sandbox/Drift-Experiment/analyze_withness_timing.py`
- `sandbox/Drift-Experiment/track12_pair4_electron.fromfile`
- `sandbox/Drift-Experiment/track12_pair4_positron.fromfile`
- `HANDOFF.md`
- `sandbox/Drift-Experiment/TASK_STATE.md`

Design decisions:

- The new witness generator writes one electron and one positron witness with:
  - kinetic energy `313 keV`;
  - `gamma = 1.612525720454`;
  - `beta = 0.784487091427`;
  - normalized momentum magnitude `p/(mc) = 1.265005612290`.
- Current sign convention:
  - c1/e- has `py = +1.265005612290`, `px = pz = 0`;
  - c2/e+ has `py = -1.265005612290`, `px = pz = 0`.
- The deck keeps `APERTURE = "CIRCLE(0.0002)"`, `NXY = 16`, and `NZ = 32`.
  This is deliberately a compact c0-source mesh with 100 um radius and 12.5 um
  transverse cells.
- The deck documents that the 30 cm scale belongs to the witness observation
  domain. It should be handled by source-to-witness gathering, not by enlarging
  the c0 Poisson mesh.

Verification:

- Generated the transverse FROMFILEs:
  `~/.venv-h6/bin/python sandbox/Drift-Experiment/make_transverse_witness_fromfiles.py --output-dir sandbox/Drift-Experiment --kinetic-energy-kev 313 --direction y`
- Syntax checks:
  - `~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/make_transverse_witness_fromfiles.py`
  - `~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/analyze_withness_timing.py`
- OPALX smoke run:
  - run directory:
    `sandbox/Drift-Experiment/withness_transverse_pair4_step1_2`
  - command:
    `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
- Analyzer:
  `env MPLBACKEND=Agg MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py --run-dir sandbox/Drift-Experiment/withness_transverse_pair4_step1_2`

Results:

- OPALX completed and wrote c0/c1/c2 H5/stat output.
- BeamBeam diagnostics reached active BeamBeam with both witnesses emitted:
  `c1:n=1,c2:n=1`.
- Pair-4 timing remains aligned:
  - analytic `witness_t0 = 13.342592708656 ps`;
  - first c1/c2 H5 dump `t = 13.348 ps`;
  - first dump is `5.407 fs` after `T0`;
  - c0 is at `s = 4.001621 mm`, `+1.621 um` past the IP.
- H5 kinematics at first c1/c2 dump:
  - c1: `py = +1.264995`, `pz = -8e-6`;
  - c2: `py = -1.264987`, `pz = -1.5e-5`.
- The written smoke window ends with c1/c2 at about `y = +/-88.7 um`, so this
  first case remains inside the 100 um compact source aperture. A longer window
  or a true source-to-witness gather test is still needed for the 30 cm
  observation-domain problem.
- New output files in the run directory:
  - `timing_mesh_summary.csv`
  - `witness_field_samples.csv`
  - `witness_kinematics_summary.csv`
  - `witness_timing_overview.png`
  - `c0_injection_mesh_xy.png`
  - `witness_transverse_motion.png`

Next step:

- Exercise the field sampling after the witnesses move outside the compact
  source aperture, without enlarging the c0 mesh to the witness domain.
- The likely implementation direction is a source-field gather/evaluation at
  witness coordinates, using the compact c0 field solve as the source.

## Earlier c0-only drift validation

Goal: set up a c0-only, 30 cm drift OPALX input and quantitative
space-charge field validation against the analytic Gaussian solution using
`sandbox/compare-e-fields`.

Changed files:

- `src/Algorithms/ParallelTracker.cpp`
- `src/Algorithms/ParallelTracker.h`
- `sandbox/Drift-Experiment/spacecharge_drift_30cm.in`
- `sandbox/Drift-Experiment/run_spacecharge_drift_validation.py`
- `sandbox/Drift-Experiment/README.md`
- `sandbox/Drift-Experiment/TASK_STATE.md`

Design decisions:

- Added an opt-in normal space-charge diagnostic controlled by
  `OPALX_SC_FIELD_CSV` and `OPALX_SC_FIELD_STEPS`.
- The diagnostic dumps c0 fields after `computeSelfFields()` and after the
  beam-frame fields are transformed back to the reference frame.
- The deck contains only container 0, an isotropic Gaussian electron bunch, a
  30 cm drift, and one integration step.
- The runner compares OPALX fields against the analytic isotropic Gaussian
  with total charge `-1e-9 C` and `sigma = 1 mm`.

Tests run:

- `cmake --build build_openmp -j 8 --target opalx`
- `cmake --build build_openmp -j 8 --target opalx_exe`
- `PYTHONPYCACHEPREFIX=/tmp/opalx-beambeam-pycache ~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/run_spacecharge_drift_validation.py sandbox/compare-e-fields/compare_gaussian_pic_fields.py`
- `~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_drift_validation.py --force`
  - OPALX one-step run completed and wrote the c0 space-charge field dump.
  - Initial comparator invocation failed because `argparse` parsed negative
    `-1e-09` as an option token; fixed by passing `--charge=-1e-09`.
- `~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_drift_validation.py --skip-run`
  - Reused the raw OPALX dump.
  - Analytic Gaussian centered on the dumped c0 centroid.
  - `relative_vector_L2_vs_analytic = 1.0543e-2`.
  - Strict threshold `3/sqrt(400000) = 4.743e-3`.
  - Current result is `passed = False` under that strict total-particle
    statistical scale.
- Charge-sign sanity check:
  - Running the comparator with positive analytic charge gives
    `relative_vector_L2_vs_analytic = 2.003`, confirming that the electron
    field sign in the validation runner is correct.
- Convergence study:
  - Command:
    `~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_convergence.py --force --threads 8`
  - Initial run completed 24/27 cases and failed at `N_GRID=64`, `NPPG=14`,
    seed 42 because the filesystem was full while HDF5 wrote the initial H5
    particle dump.
  - Removed generated H5/raw artifacts from the current task, then resumed with:
    `~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_convergence.py --resume --threads 8`
  - All 27 requested runs completed.
  - Aggregate results are in
    `sandbox/Drift-Experiment/convergence_results/convergence_grouped_summary.csv`.
  - Mean relative vector L2 by `(N_GRID, NPPG)`:
    - `(16, 5)`: `1.06175e-1`
    - `(16, 8)`: `1.05434e-1`
    - `(16, 14)`: `1.04669e-1`
    - `(32, 5)`: `2.42086e-2`
    - `(32, 8)`: `2.44037e-2`
    - `(32, 14)`: `2.40839e-2`
    - `(64, 5)`: `6.68363e-3`
    - `(64, 8)`: `5.19025e-3`
    - `(64, 14)`: `4.39455e-3`
  - Aggregate plots:
    `relative_l2_vs_grid.png`, `relative_l2_vs_nppg.png`, and
    `relative_l2_vs_particles.png`.
- N_GRID=128 feasibility test:
  - Command:
    `~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_convergence.py --output-dir sandbox/Drift-Experiment/convergence_ngrid128_test --n-grid-values 128 --nppg-values 5 --seeds 42 --threads 8 --force`
  - Initial invocation exposed a relative-output-directory path bug in
    `run_spacecharge_convergence.py`; fixed the OPALX launch to pass
    `deck.name` from the case directory and an absolute
    `OPALX_SC_FIELD_CSV` path.
  - Case completed with `10,485,760` particles.
  - Result:
    `relative_vector_L2_vs_analytic = 6.94850e-3`.
  - Raw diagnostic CSV grew to about `3 GB` during the run and was removed by
    the driver after metrics were collected.
  - Output directory:
    `sandbox/Drift-Experiment/convergence_ngrid128_test/`.

Next step:

- Inspect whether the remaining `~1%` L2 residual is expected PIC/grid
  sampling noise for a 64^3 open solver with 400k macroparticles, or a
  branch-introduced error in the normal space-charge path.
- Run a multi-rank repeat if this validation becomes part of final acceptance.

## 2026-07-03 master merge reverted

Current repository state:

- The pending `origin/master` merge was aborted.
- The committed `origin/master` merge `902c9db0a` and its follow-up fix
  `c590470a4` were reverted together in commit `1177cfa4a`
  (`Revert origin/master merge`).
- The tracked tree is back on the BeamBeam branch lineage; the remaining
  Drift-Experiment files are untracked sandbox work.
- The local IPPL checkout in `build_openmp/_deps/ippl-src` was reset to
  `27d11d0b58bc5db6a582f8132c3d6c88285b26bd`, matching this branch's
  `solver_recv` API usage.
- A patch of the local Drift-Experiment/ParallelTracker edits that existed
  before the merge abort is saved at:
  `/tmp/opalx-beambeam-local-before-merge-abort.patch`.

Verification after revert:

- `cmake -S . -B build_openmp`
- `cmake --build build_openmp -j 8 --target opalx_exe`
  - Passed after resetting the local IPPL checkout to `27d11d0b`.

## 2026-07-03 redone 30 cm drift grid/NPPG study

Current goal:

- Redo `spacecharge_drift_30cm.in` with the full requested grid and particle
  number grid, and produce the radial/profile plots in the H5-based plotting
  style used in the later drift checks.

Changed files:

- `src/Algorithms/ParallelTracker.cpp`
- `src/Algorithms/ParallelTracker.h`
- `sandbox/Drift-Experiment/spacecharge_drift_30cm.in`
- `sandbox/Drift-Experiment/run_spacecharge_convergence.py`

Design decisions:

- Kept the reverted branch's solver setting `PARFFTZ = FALSE` in the 30 cm
  deck; only enabled `OPTION, EBDUMP = TRUE`.
- Restored an opt-in H5 space-charge diagnostic via
  `OPALX_SC_FIELD_H5_STEPS=0`; the diagnostic dumps the primary bunch after
  self-field computation and before `bunchUpdate()`.
- The convergence driver now reads H5, computes the analytic isotropic
  Gaussian field at the dumped particle coordinates, writes per-case radial
  profiles/plots, and deletes large H5 dumps unless `--keep-raw` is passed.
- Because only the grid/NPPG matrix was requested this turn and disk space was
  about 1.6 GiB, reran the 3x3 grid for seed 42 rather than the earlier
  three-seed/27-run ensemble.

Verification:

- `~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/run_spacecharge_convergence.py`
- `cmake --build build_openmp -j 8 --target opalx_exe`
- `~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_convergence.py --output-dir sandbox/Drift-Experiment/redone_full_grid_seed42 --n-grid-values 16,32,64 --nppg-values 5,8,14 --seeds 42 --threads 8 --force --sample-particles 50000`

Output:

- Output directory:
  `sandbox/Drift-Experiment/redone_full_grid_seed42`
- Raw H5 dumps were deleted after metric/plot extraction; no `.h5` files remain
  in the output directory.
- Relative vector L2 vs analytic:
  - `(16, 5)`: `1.07615e-1`
  - `(16, 8)`: `1.02523e-1`
  - `(16, 14)`: `1.06806e-1`
  - `(32, 5)`: `2.66366e-2`
  - `(32, 8)`: `2.26941e-2`
  - `(32, 14)`: `2.38546e-2`
  - `(64, 5)`: `6.74797e-3`
  - `(64, 8)`: `5.08562e-3`
  - `(64, 14)`: `4.18043e-3`
- Main plots:
  - `relative_l2_vs_grid.png`
  - `radial_field_vs_analytic_all_cases.png`
  - `signed_radial_field_vs_analytic_all_cases.png`
  - `radial_field_relative_error_all_cases.png`
- Per-case directories contain the rendered input, `metrics.csv`,
  `radial_profile.csv`, sampled particle comparison CSV, and the three
  per-case radial PNGs with particle scatter plus binned OPALX/analytic curves.

Remaining risks:

- This redo was single-rank/OpenMP only, matching the immediate plotting task.
  A multi-rank repeat is still needed if this diagnostic becomes an acceptance
  test.

## 2026-07-03 `spacecharge_drift_withness.in` timing smoke

Goal:

- Establish witness injection timing for
  `sandbox/Drift-Experiment/spacecharge_drift_withness.in`.
- Follow the `TestParticleOrbitSimulation.pptx`/track12 timing convention.
- Determine whether the initial BeamBeam mesh setup contains c0 and whether it
  can sample c0 fields at the pair-4 witness position.

Input/deck decisions:

- Use pair 4 first:
  - `witness_ct_m = 0`
  - `bb_ip_s = 4.0e-3 m`
  - `witness_t0 = bb_ip_s / (primary_beta * CLIGHT)`
  - analyzer value: `13.342592708656 ps`
- Use a 1 mm upstream drift and a 6 mm BeamBeam element:
  - BeamBeam element edge at `1 mm`
  - BeamBeam IP at `4 mm`
  - BeamBeam window spans `1 mm` to `7 mm`
- Keep the smoke run cheap:
  - `primary_macroparticles = 100`
  - `primary_charge_scale = 1.0e-5`
  - `NXY = 16`, `NZ = 32`
  - `APERTURE = "CIRCLE(0.0002)"`
  - `COPY_TIME = 0.0`

Important timing fix:

- The first staged run used `coarse_steps = 200`.
- OPALX reached `fine_start_s = 3.95 mm` around `13.2 ps`, but kept advancing
  the clock to `20 ps` before starting the fine 1 fs segment.
- Fixed by setting:
  - `coarse_steps = 132`
  - `fine_steps = 520`
- This aligns `MAXSTEPS` with the first `ZSTOP`:
  `fine_start_s / (primary_beta * CLIGHT * coarse_dt) = 131.76`.

C++ fix required for the no-copy case:

- `ParallelTracker::performBeamBeamWindowEntryTransition()` now skips the
  copied-source deposited-charge diagnostic if
  `BEAMBEAM::copyTimeReached(...)` is false.
- Without this, the run aborts at BeamBeam entry with:
  `Missing deposited-charge diagnostics for the pre-enlarge BeamBeam solve.`
- This is scoped to `COPY_TIME = 0.0`/no-copy cases; copied-source cases still
  run the transition diagnostic.

Verification:

- `cmake --build build_openmp -j 8 --target opalx_exe`
- OPALX run from:
  `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed`
- Command:
  `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
- Analyzer:
  `env MPLBACKEND=Agg MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py`
- `~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/analyze_withness_timing.py`

Output:

- `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/run.log`
- `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/timing_mesh_summary.csv`
- `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/witness_field_samples.csv`
- `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/witness_timing_overview.png`
- `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/c0_injection_mesh_xy.png`

Timing result:

- `BB-DIAG` transitions were:
  - inactive, c1/c2 empty
  - active, c1/c2 empty
  - active, c1/c2 each have one particle
  - active with source inactive after c0 exits the source-active interval
- First synchronized c0/witness H5 dump:
  - global step `279`
  - `t = 13.348 ps`
  - `t - witness_t0 = 5.407 fs`
  - c0 `s = 4.001621 mm`
  - c0 is `+1.621 um` past the IP

Mesh result:

- Longitudinal window coverage is sufficient:
  - c0 center at first witness dump: `4.001621 mm`
  - c0 longitudinal half extent: about `1.633 mm`
  - occupied span: approximately `2.37 mm` to `5.63 mm`
  - BeamBeam window: `1 mm` to `7 mm`
- Transverse containment is safe but under-resolved:
  - aperture half-width: `100 um`
  - cell size for `NXY=16`: `12.5 um`
  - c0 rms at injection: `2.125 um` by `1.841 um`
  - c0 max transverse extents: `4.87 um` by `5.55 um`
  - most of the primary is inside one central grid cell.

Remaining risks / next steps:

- The smoke run used only 100 primary macroparticles and `1e-5` primary charge
  scale, so it is valid for timing/mesh mechanics only.
- Prior track12 notes show that aperture/cell-size tuning alone can create
  large sign and magnitude pathologies in sampled witness fields.  The next
  step should be a controlled mesh/particle scan, not a blind aperture shrink.
- The next quantitative comparison should sample c0 fields at c1/c2 at the
  same H5 time and compare to the smooth analytic/manufactured pair-4
  reference after restoring adequate macroparticle statistics.

## 2026-07-04 Release T0=4 ps / 1600 fs witness run

Current goal:

- Rebuild OPALX in Release mode.
- Start c1/c2 witnesses at absolute `T0 = 4 ps`.
- Keep c0 active for `1600 fs` after witness T0, then retire c0 and verify
  post-retirement ballistic witness continuation.

Changed files:

- `sandbox/Drift-Experiment/spacecharge_drift_withness.in`
- `sandbox/Drift-Experiment/analyze_withness_timing.py`
- `HANDOFF.md`
- `sandbox/Drift-Experiment/TASK_STATE.md`

Build:

- `cmake -S . -B build_openmp -DBUILD_TYPE=Release -DCMAKE_BUILD_TYPE=Release`
- `cmake --build build_openmp -j 8 --target opalx_exe`
- Verified CMake cache says:
  - `BUILD_TYPE:STRING=Release`
  - `CMAKE_BUILD_TYPE:STRING=Release`

Deck decisions:

- `witness_t0 = 4.0e-12`
- `near_ip_active_time = 1.6e-12`
- `primary_retire_time = 5.6e-12`
- `APERTURE = "CIRCLE(0.0010)"`, giving a compact 500 um radius source field
  aperture for the 1600 fs active window.
- Added shifted source timing:
  `primary_source_r0z = bb_ip_s - primary_beta * CLIGHT * witness_t0 - witness_ct_m`
  and applied it as `R0Z = primary_source_r0z` on the primary emission source.
- `fine_start_s`/`fine_stop_s` now include `primary_source_r0z`, so staged
  `ZSTOP` values are absolute lattice positions.

Why the shift matters:

- A first Release T0=4 ps attempt in
  `withness_t0_4ps_active_1600fs_release` left c0 launched from `R0Z=0`.
- In that run, c0 was still upstream of the BeamBeam active window when c1/c2
  were injected, so no useful c0 field sampling occurred before c0 retirement.
- The corrected setup follows the earlier `track12particles` convention and
  makes c0 reach `bb_ip_s` at witness T0.

Corrected run:

- Directory:
  `sandbox/Drift-Experiment/withness_t0_4ps_active_1600fs_release_shifted_c0`
- OPALX command:
  `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
- Analyzer command:
  `env MPLBACKEND=Agg MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py --run-dir sandbox/Drift-Experiment/withness_t0_4ps_active_1600fs_release_shifted_c0 --witness-t0-s 4.0e-12 --near-ip-active-time-s 1.6e-12 --aperture-radius-m 5.0e-4`

Verification:

- `~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/analyze_withness_timing.py`
- OPALX run completed.
- `BB-DIAG` showed:
  - `BB-state=Active` with c0 active before witness population
  - `BB-state=Active` with `c1:active:n=1,c2:active:n=1`
  - `BB-state=Completed` with source retirement pending
  - final completed state with c0 source retired

Timing/mesh result:

- first witness H5 dump: `4.020 ps`, `20.0 fs` after T0
- configured retire time: `5.600 ps`, `1600 fs` after T0
- inferred `primary_source_r0z = 2.800833 mm`
- nearest c0/IP H5 dump: `t = 4.000 ps`, absolute `s = 4.000 mm`
- at first witness H5 dump, c0 absolute `s = 4.005996 mm`
  (`+5.996 um` downstream of IP)
- c0 rms at first witness dump:
  - `x = 2.125 um`
  - `y = 1.841 um`
  - `s = 0.573 mm`
- transverse cell size: `62.5 um`
- c1/c2 final H5 dump: `5.780 ps`
- sampled fields in `near_ip_field_cutoff.csv` are zero from `5.620 ps`
  onward, after c0 retirement.

Generated products:

- `witness_timing_overview.png`
- `c0_injection_mesh_xy.png`
- `witness_transverse_motion.png`
- `near_ip_field_cutoff.png`
- `timing_mesh_summary.csv`
- `witness_kinematics_summary.csv`
- `witness_field_samples.csv`
- `near_ip_field_cutoff.csv`

Remaining risks:

- The case is still only a timing/mesh smoke test:
  `NXY=16`, `NZ=32`, `primary_macroparticles=100`, and
  `primary_charge_scale=1e-5`.
- The sampled near-IP fields are visibly noisy and not yet a quantitative
  validation of the witness kicks.

## 2026-07-04 clean rerun with `--info 5`

Run directory:

- `/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/withness_t0_4ps_active_1600fs_release_shifted_c0`

Input:

- `spacecharge_drift_withness.in`
- User-edited high-statistics settings:
  - `NXY = 32`
  - `NZ = 64`
  - `primary_macroparticles = NXY*NXY*NZ*8 = 524288`
  - `witness_t0 = 4.0e-12`
  - `near_ip_active_time = 1.6e-12`

Command:

```sh
/usr/bin/time -p env OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads \
  ../../../build_openmp/src/opalx --info 5 spacecharge_drift_withness.in \
  > spacecharge_drift_withness.out \
  2> spacecharge_drift_withness.err
```

Result:

- OPALX exited with status `0`.
- `spacecharge_drift_withness.err` contains only:
  - `real 384.57`
  - `user 1953.02`
  - `sys 31.17`
- No `No space`, HDF5, MPI-IO, or failed-write messages were found in
  `spacecharge_drift_withness.out`/`.err`.
- `spacecharge_drift_withness.out` is about `12 MiB` because of `--info 5`.

Output coverage:

- Final stat rows:
  - c0: `t = 5.780 ps`, `numParticles = 0`, so c0 was retired.
  - c1: `t = 5.780 ps`, `numParticles = 1`, `mean_y = +418.531 um`.
  - c2: `t = 5.780 ps`, `numParticles = 1`, `mean_y = -418.487 um`.
- H5 coverage:
  - c0 H5: `82` dumps, from `2.000 ps` to `5.600 ps`; last dump is the
    retirement-time c0 source dump with `524288` particles.
  - c1 H5: `89` dumps, from `4.020 ps` to `5.780 ps`.
  - c2 H5: `89` dumps, from `4.020 ps` to `5.780 ps`.
- Key `--info 5` timers:
  - `BB window total`: `335.851 s`
  - `BB self fields`: `322.737 s`
  - `BB CIC scatter`: `144.841 s`
  - `Solve: Electric field`: `62.3396 s`
  - `FFT: Efield`: `58.9504 s`
  - `BB witness gather`: `0.717866 s`
