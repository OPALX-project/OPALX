# OPALX BeamBeam Handoff

Last updated: 2026-07-05

## Current Active Task: `spacecharge_drift_withness.in` Timing And Mesh

Goal:

- Continue with `sandbox/Drift-Experiment/spacecharge_drift_withness.in`.
- Establish the timing for witness-beam injection relative to the primary c0
  bunch, following `TestParticleOrbitSimulation.pptx`.
- Determine how the BeamBeam mesh must be configured to sample c0 fields at the
  witness position.

2026-07-05 edge-at-witness timing update:

- User observed that the previous `T0 = 4 ps` setup placed the witnesses in a
  strong c0 field immediately.  That was expected from the old deck formula:
  `primary_source_r0z = bb_ip_s - beta*c*T0`, which placed the c0 centroid at
  the IP at witness injection.
- Updated the active deck
  `sandbox/Drift-Experiment/spacecharge_drift_withness.in` so the c0
  `+3 sigma_z` leading edge reaches the first track12 witness at `T0`:
  - `primary_edge_sigma = 3.0`
  - `track12_first_ct_m = -1.5 * primary_sigma_z`
  - `track12_edge_reference_s = bb_ip_s + track12_first_ct_m`
  - `primary_centroid_at_witness_t0 =
    track12_edge_reference_s - primary_edge_sigma * primary_sigma_z`
  - `primary_source_r0z =
    primary_centroid_at_witness_t0 - primary_beta * CLIGHT * witness_t0`
- Important geometry point:
  - With the old `bb_ip_s = 4 mm` and BeamBeam range `1..7 mm`, this edge
    alignment put the c0 centroid at `1.3 mm` and its `-3 sigma_z` tail at
    about `-0.5 mm` at `T0`.
  - OPALX therefore kept BeamBeam `Inactive`, because the activation condition
    requires the source bunch tail to be inside the BeamBeam window:
    `bunchExtent.tail >= beginS + half longitudinal cell`.
  - Diagnostic directory for this failed activation test:
    `sandbox/Drift-Experiment/one_c0_track12_edge_t0_smoke`.
- Corrected active setup:
  - Kept the upstream drift at `1 mm`.
  - Increased `bb_length` to `10 mm`, so
    `bb_ip_s = drift_length + 0.5 * bb_length = 6 mm` and the BeamBeam range is
    `1..11 mm`.
  - At `T0 = 4 ps`, c0 centroid is `3.3 mm`, c0 `+3 sigma_z` edge is
    `5.1 mm`, and the first track12 witness is also at `5.1 mm`.
  - The c0 tail is then about `1.5 mm`, inside the `1..11 mm` BeamBeam window,
    so BeamBeam activates before witness injection.
- Corrected diagnostic run:
  - Directory:
    `sandbox/Drift-Experiment/one_c0_track12_edge_t0_ip6mm_smoke`
  - Run-copy-only output change: `C0PSDUMPFREQ = 20`, to save c0 near `T0`.
  - Command:
    `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads /usr/bin/time -p ../../../build_openmp/src/opalx --info 2 spacecharge_drift_withness.in > spacecharge_drift_withness.out 2> spacecharge_drift_withness.err`
  - Completed normally in `real 32.98 s`, `user 27.42 s`, `sys 22.32 s`.
  - BeamBeam diagnostics:
    - initially inactive before the source tail enters the window;
    - entered BeamBeam with
      `interaction_point_s=0.006 m, s_range=(0.001, 0.011) m,
      length=0.010 m, copy_time=NONE`;
    - c1/c2 each loaded 6 track12 particles and remained active;
    - c0 retired at `RETIRE_TIME` while witness containers stayed active.
- H5 timing summary from the corrected run:
  - first witness H5 dump: `4.02 ps`, `20 fs` after `T0`;
  - `primary_source_r0z = 2.100833 mm`;
  - `track12_edge_reference_s = 5.1 mm`;
  - c0 centroid at `T0`: `3.3 mm`;
  - c0 `+3 sigma_z` edge at `T0`: `5.1 mm`;
  - nearest saved c0 to the first witness dump is synchronized at global step
    `59`, because this diagnostic used `C0PSDUMPFREQ = 20`.
- Plotting updates:
  - `sandbox/Drift-Experiment/analyze_withness_timing.py` now uses Matplotlib
    `Agg`;
  - the timing plot includes c0 centroid, c0 `+3 sigma_z` edge, BeamBeam IP,
    and first track12 witness lines;
  - the default witness kinetic energy remains the raw track12
    `gamma = 2` value.
- Corrected-run plots:
  - `witness_timing_overview.png`
  - `near_ip_field_cutoff.png`
  - `c0_injection_mesh_xy.png`
  - `witness_transverse_motion.png`
  - `track12_figure1_opalx_h5.png`
  - `track12_figure1_opalx_h5_t.png`
- Current caveat:
  - This fixes the initial c0-edge timing for the current one-source `+z`
    plumbing model.
  - It still does not make the model equivalent to CAIN if CAIN used a
    truly oncoming c0 bunch moving in `-z`.

2026-07-05 one-c0 track12 smoke setup:

- User clarified the immediate model:
  - use only one physical c0 bunch;
  - use the raw `track12_electrons.fromfile` and
    `track12_positrons.fromfile` witnesses;
  - no mirrored/copy field source for now;
  - c0 may be modeled as coming from the left in the current OPALX `+z`
    convention.
- Updated the active deck
  `sandbox/Drift-Experiment/spacecharge_drift_withness.in` accordingly:
  - title/comment block now describes the one-c0/track12 model;
  - `OPTION, C0PSDUMPFREQ = 600`;
  - `primary_charge_scale = 1.0`;
  - `primary_macroparticles = 100` for the low-resolution smoke run;
  - `n_witness_per_species = 6`;
  - witness distributions use `track12_electrons.fromfile` and
    `track12_positrons.fromfile`;
  - `COPY_TIME = 0.0`, so `copy_time=NONE` and only container 0 deposits the
    BeamBeam source field;
  - `IP_S = bb_ip_s` pins the interaction point explicitly;
  - c0 source is shifted so the c0 centroid reaches `bb_ip_s = 4 mm` at the
    common witness `T0 = 4 ps`.
- Copied the full track12 FROMFILEs into `sandbox/Drift-Experiment/` for the
  base deck.
- Updated `sandbox/Drift-Experiment/analyze_withness_timing.py`:
  - forces Matplotlib `Agg` backend for batch/non-GUI runs;
  - default witness kinetic energy is now `EMASS_GEV`, matching the raw
    track12 `gamma = 2` (`Ps/(m_e c) = sqrt(3)`) setup.
- Smoke run directory:
  `sandbox/Drift-Experiment/one_c0_track12_smoke`
- Run command:
  `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads /usr/bin/time -p ../../../build_openmp/src/opalx --info 2 spacecharge_drift_withness.in > spacecharge_drift_withness.out 2> spacecharge_drift_withness.err`
- Run result:
  - completed normally in `real 31.12 s`, `user 25.43 s`, `sys 20.63 s`;
  - BeamBeam entered with
    `interaction_point_s=0.004 m, s_range=(0.001, 0.007) m,
    length=0.006 m, copy_time=NONE`;
  - `FromFile` loaded 6 electrons and 6 positrons;
  - c1/c2 became active at the common `T0 = 4 ps`;
  - c0 was retired at `RETIRE_TIME` with `marked 100, deleted 100,
    remaining 0`; witness containers remained active;
  - OPALX reports the witness kinetic energy around `511 keV`, as expected for
    the raw track12 `gamma = 2` file, not the older `313 keV` transverse-pair
    model.
- Analysis output in the smoke directory:
  - `timing_mesh_summary.csv`
  - `witness_kinematics_summary.csv`
  - `witness_field_samples.csv`
  - `witness_timing_overview.png`
  - `c0_injection_mesh_xy.png`
  - `witness_transverse_motion.png`
  - `near_ip_field_cutoff.png`
  - `track12_figure1_opalx_h5.png`
  - `track12_figure1_opalx_h5_t.png`
- Timing details from the H5/stat analysis:
  - first witness H5 dump is at `4.02 ps`, i.e. `20 fs` after witness `T0`;
  - because `C0PSDUMPFREQ = 600`, the nearest saved c0 H5 dump to the first
    witness dump is at `4.56 ps`, with c0 centroid about `167.9 um` beyond the
    IP in the shifted absolute convention;
  - to inspect c0 exactly at `T0`, reduce `C0PSDUMPFREQ` for the next short
    diagnostic run.
- Current physics caveat:
  - This smoke run injects all six track12 witnesses simultaneously and
    represents the PPTX timing values as initial `z` offsets.  It is not yet
    the exact pair-wise CAIN timing model.
  - With raw `Ps -> pz > 0` and c0 also moving in `+z`, this is a
    co-propagating one-source plumbing test.  The OPALX trajectories therefore
    should not yet be expected to reproduce the CAIN/head-on figure
    quantitatively.

2026-07-05 staged coarse-to-fine rerun:

- New run directory:
  `sandbox/Drift-Experiment/withness_t0_4ps_coarse_to_12p5ps_fine_to_30ps`
- Copied the completed centered-IP long-run input and changed only the staged
  timing:
  - `coarse_to_time = 12.5e-12`
  - `final_time = 30.0e-12`
  - `post_retire_observation_time = final_time - primary_retire_time`
  - `fine_start_s = primary_source_r0z + primary_beta * CLIGHT * coarse_to_time`
  - `fine_stop_s = primary_source_r0z + primary_beta * CLIGHT * final_time`
  - `coarse_steps = 125` with `coarse_dt = 100 fs`
  - `fine_steps = 17500` with `fine_dt = 1 fs`
- c0 source retirement is unchanged at `primary_retire_time = 20 ps`.
- BeamBeam geometry is unchanged:
  `IP_S = 4 mm`, active range `1 mm` to `17 mm`.
- Started OPALX with:
  `env OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads ../../../build_openmp/src/opalx --info 2 spacecharge_drift_withness.in > spacecharge_drift_withness.out 2> spacecharge_drift_withness.err`
- Run result:
  - Completed normally.
  - `/usr/bin/time` reported `real 2036.46 s`, `user 10587.46 s`,
    `sys 180.24 s`.
  - Final stderr contains only the timing line.
  - BeamBeam entered with `interaction_point_s=0.004 m,
    s_range=(0.001, 0.017) m, length=0.016 m`.
  - Witnesses reached active state `c1:active:n=6,c2:active:n=6`.
  - At `t=12.53 ps`, the run is in the intended `1 fs` fine segment.
  - Retirement diagnostic:
    `Retired BeamBeam source container[0] at RETIRE_TIME: marked 524288,
    deleted 524288, remaining 0. Witness containers remain active.`
  - Final BeamBeam state after retirement:
    `source_active=FALSE source_retirement_pending=FALSE`.
- H5 output:
  - c0 H5 has 12 saved dumps, from global step `599` at `t=12.975 ps`
    through step `7199` at `t=19.575 ps`.
  - The last saved c0 position is `s=5.8684 mm` in OPALX local SPOS, or
    `s=8.6693 mm` in the shifted absolute convention.
  - c1/c2 H5 have 880 saved dumps, from `t=6.0 ps` to `t=30.0 ps`.
  - c1/c2 post-retirement field samples are zero:
    `max_post_retire_E = 0.0`.
  - Because the first segment uses `100 fs` steps with `PSDUMPFREQ=20`, the
    first witness H5 dump is at `t=6.0 ps`, not immediately after the
    `4.0 ps` witness T0.
  - Because `C0PSDUMPFREQ=600`, the exact c0 state at witness injection and at
    the `12.5 ps` coarse/fine transition is not saved; the first c0 mesh
    snapshot is at `t=12.975 ps`, `s=6.6906 mm` in the shifted absolute
    convention.
- Witness extent by the final `30 ps` H5 dump:
  - c1 all-particle maxima: `max |x| = 110.0 um`,
    `max |y| = 676.6 um`, `max |z| = 9.3305 mm`.
  - c2 all-particle maxima: `max |x| = 1572.9 um`,
    `max |y| = 1025.5 um`, `max |z| = 9.2888 mm`.
  - Some witness particles are outside the compact `500 um` c0-source aperture
    during the post-retirement ballistic segment, as expected for this longer
    observation.
- Plots produced in the run directory:
  - `witness_timing_overview.png`
  - `c0_injection_mesh_xy.png`
  - `witness_transverse_motion.png`
  - `near_ip_field_cutoff.png`
  - `track12_figure1_opalx_h5.png`
  - `track12_figure1_opalx_h5_t.png`
- Plot observations:
  - The timing plot shows zero witness fields after c0 retirement at `20 ps`.
  - The sampled OPALX relative field before retirement still does not follow
    the simple line-Gaussian decay estimate.
  - The Figure 1 comparison remains qualitatively different from the
    `TestParticleOrbit.dat`/CAIN reference, especially for the blue trajectories
    after roughly `16 ps`.

2026-07-05 PPTX timing/retirement re-check:

- Extracted `sandbox/TestParticleOrbitSimulation.pptx` slide text and checked it
  against `sandbox/TestParticleOrbit.dat`.
- PPTX facts:
  - c0/source bunch: `sigma_z = 0.6 mm`, Gaussian longitudinal distribution cut
    at `3 sigma_z`.
  - Test particles are right-going and inserted at
    `ct/sigma_z = -1.5, -1.0, -0.5, 0, 0.5, 1.0`.
  - `ct/sigma_z = -1.5` is described as the moment the test particle touches
    the on-coming bunch edge.
- `TestParticleOrbit.dat` matches this:
  - initial `ct` values are approximately `-0.9, -0.5994, -0.3006, 0,
    +0.3006, +0.5994 mm`;
  - trajectories run to `ct = +1.8 mm = +3 sigma_z`.
- Current combined OPALX deck is not the same timing model:
  - it injects all 12 witnesses at one absolute `T0 = 4 ps`;
  - it encodes the six PPTX timings as initial witness `z` offsets;
  - the completed staged run used `track12_electrons.fromfile` and
    `track12_positrons.fromfile`, where `Ps/(m_e c) = sqrt(3)` was mapped into
    OPALX `pz`.  This is intentional if the goal is to reproduce the
    PPTX/CAIN track12 reference, because the PPTX says the particles are
    artificial right-going test particles and the file columns are
    `Px, Py, Ps`.
- Separate transverse BW-pair estimate:
  - If instead we switch back to the separate low-energy BW-pair model, the
    witness momentum should be transverse to the c0 `beta_z c` direction,
    i.e. `pz ~= 0` for the witness in the c0-z frame.
  - Under that transverse-witness interpretation, a co-propagating
    `beta_c0 - beta_w` estimate is not the right retirement criterion.
  - With transverse witnesses, the longitudinal separation from c0 grows as
    `beta_c0 c dt`.  At `dt = 16 ps` after pair-4 insertion, c0 has moved
    `4.797 mm` in z; relative to the initially leading `ct = +0.5994 mm`
    offset this is still `4.197 mm`, greater than `3 sigma_z = 1.8 mm`.
- Retirement implication:
  - If the intended physics is PPTX/oncoming, `RETIRE_TIME = 20 ps` is not too
    early; a `3 sigma_z` source cutoff passes a pair-4 witness in about
    `3.22 ps` after pair-4 insertion, and the CAIN output itself ends at
    `6.00 ps` after pair-4 insertion.
  - If the track12/PPTX right-going witness and c0 are both modeled in OPALX as
    moving in `+z`, `RETIRE_TIME = 20 ps` is early for the co-propagating
    longitudinal-support estimate.  With `gamma_w = 2`,
    `beta_w = 0.8660254`; with c0 at `245 MeV`, `beta_c0 = 0.9999978`.  The
    relative speed is only `(beta_c0 - beta_w)c`, so after `16 ps` c0 is only
    `0.643 mm` ahead of the `ct=0` witness and only `0.043 mm` ahead of the
    initially leading `ct=+0.5994 mm` pair-6 witness.
- Practical conclusion:
  - Do not rotate `TestParticleOrbit.dat` `Ps` into a transverse OPALX momentum
    if the goal is to reproduce the PPTX/CAIN track12 reference.  The PPTX says
    the test particles are artificial, right-going particles, and the file
    format names the momentum components `Px, Py, Ps`; for the OPALX local
    beamline frame this means `Ps -> pz`.
  - The earlier transverse-witness assumption belongs to the separate
    low-energy Breit-Wheeler pair model in `sandbox/note`, where the pair has
    about `313 keV` kinetic energy and is launched transversely in the `z=0`
    plane.  That is a different witness model from the raw
    `TestParticleOrbit.dat` reference, whose first row has
    `E ~= 1.022 MeV`, `Px = Py = 0`, and
    `Ps/(m_e c) = 1.73205118`.
  - Therefore the next decision is model selection:
    reproduce the track12/PPTX artificial right-going reference with
    `Ps -> pz`, or switch back to the separate transverse BW-pair model and use
    a different/generated FROMFILE.
  - The existing pair-wise generator
    `sandbox/track12particles/opalx/generate_timing_pair_inputs.py` already
    documents the timing-correct interpretation: the PPTX `ct` values are
    insertion times, not one simultaneous spatial snapshot.
- Agreed source-geometry decision:
  - Current simplification from user: use only one physical c0 bunch and the
    `track12_electrons.fromfile` / `track12_positrons.fromfile` witnesses.
    Do not add a mirrored/copy-only source mode for this stage.
  - Keep `COPY_TIME = 0` / copy disabled so the field source is only container
    0.  The earlier `COPY_ONLY` idea is deferred and should not be implemented
    for this immediate run.
  - Use OPALX's cleanly supported direction first: c0 moves in increasing `z`
    through the BeamBeam window.  This can be described as "from left to right"
    in the current local coordinate convention.
  - If the track12 witnesses also keep their raw `Ps -> pz > 0` momenta, this
    one-source OPALX run is co-propagating.  c0 catches them only with relative
    speed `(beta_c0 - beta_w)c`; this is acceptable as a plumbing/timing test,
    but it is not the same timescale as a head-on CAIN/oncoming-source model.
  - If we later need true CAIN head-on timing/force signs, revisit either a
    negative-direction source representation or a one-source mirrored-deposit
    implementation.  For now, prioritize the one physical c0 setup.

2026-07-04 long-run setup and result:

- User requested `C0PSDUMPFREQ = 600`, a simulation about `10x` longer, and
  all plots after the full run.
- C++ output-thinning change added:
  - New option `OPTION, C0PSDUMPFREQ = value;`
  - `-1` follows global `PSDUMPFREQ`, `0` disables c0 H5, positive values thin
    container-0 H5 on global phase-space dump steps.
  - Witness containers c1/c2 still use global `PSDUMPFREQ`.
- C++ BeamBeam geometry change added:
  - New optional `BEAMBEAM` attribute `IP_S`.
  - `IP_S = 0` preserves the old center-derived interaction point.
  - Positive `IP_S` pins the interaction point in path-length coordinates while
    the placed element range still defines the active BeamBeam window.
- Verification before the long run:
  - `cmake --build build_openmp -j 8` passed.
  - `ctest --test-dir build_openmp -R BeamBeam --output-on-failure` passed.
  - `~/.venv-h6/bin/python -m py_compile` passed for
    `sandbox/Drift-Experiment/analyze_withness_timing.py` and
    `sandbox/track12particles/note/plot_figure1_from_h5.py`.
- Current full-run directory:
  `sandbox/Drift-Experiment/withness_t0_4ps_active_16ps_centered_ip_release_shifted_c0`
- Input copied from the previous shifted-c0 1.6 ps run, then changed to:
  - `OPTION, PSDUMPFREQ = 20`
  - `OPTION, C0PSDUMPFREQ = 600`
  - `near_ip_active_time = 1.6e-11 s`
  - `post_retire_observation_time = 1.75e-12 s`
  - `fine_steps = 17800`
  - `bb_length = 1.6e-2 m`
  - `bb_ip_s = 4.0e-3 m`
  - `bb_elemedge = drift_length = 1.0e-3 m`
  - `IP_S = bb_ip_s`
- Geometry check:
  - `spacecharge_drift_withness_ips_smoke.in` completed after the C++ rebuild.
  - OPALX reported:
    `interaction_point_s=0.004 m, s_range=(0.001, 0.017) m, length=0.016 m`.
  - The full run therefore keeps the BeamBeam entry at `1 mm`, pins the IP at
    `4 mm`, and extends the active window to `17 mm`.
  - This avoids the failed negative-ELEMEDGE attempt, where the tracker entered
    BeamBeam immediately and the index-map split selected the wrong IP.
- Run command:
  `env OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads ../../../build_openmp/src/opalx --info 2 spacecharge_drift_withness.in > spacecharge_drift_withness.out 2> spacecharge_drift_withness.err`
- Run result:
  - Completed normally.
  - `/usr/bin/time` reported `real 4172.50 s`, `user 22434.58 s`,
    `sys 274.21 s`.
  - Final `spacecharge_drift_withness.err` contains only the timing line.
  - Early log confirmed `interaction_point_s=0.004 m,
    s_range=(0.001, 0.017) m, length=0.016 m`.
  - BeamBeam reached the intended sequence:
    inactive, active before witnesses, active with `c1:active:n=6,c2:active:n=6`,
    completed with `source_retirement_pending=TRUE`, then c0 retired.
  - Retirement diagnostic:
    `Retired BeamBeam source container[0] at RETIRE_TIME: marked 524288,
    deleted 524288, remaining 0. Witness containers remain active.`
  - c0 H5 thinning worked: final files were about `1.5 GB` for c0 and
    `9.2 MB` each for c1/c2.
  - c0 H5 has 26 saved dumps; the last saved c0 dump is at `t=19.56 ps`,
    before retirement, because `C0PSDUMPFREQ=600`.
  - Stdout shows c0 at `s=5.987 mm` at `t=19.969 ps`; it therefore reaches
    about `6 mm` at the configured `20 ps` retirement.
  - c1/c2 H5 continue to `t=21.80 ps`, global step `17839`.
  - H5 post-retirement c1/c2 field samples are zero:
    `max_post_retire_E = 0.0`.
- H5 analyzer summary:
  - First witness H5 dump: `t=4.02 ps`, `20 fs` after witness T0.
  - Because c0 is saved only every 600 global steps, the nearest saved c0 dump
    to witness injection is global step `599`, `540 fs` after the first
    witness dump; the exact c0-at-T0 state is not in c0 H5.
  - The sparse nearest c0 dump used for the mesh snapshot is at `t=4.56 ps`,
    `s=4.1679 mm` in the shifted absolute convention.
  - Final witness first-particle positions:
    c1 `(x,y,z) = (43.4, -225.1, -8154.7) um`,
    c2 `(x,y,z) = (-108.1, 191.6, -8052.3) um`.
  - Across all six particles, c1/c2 remain inside the compact transverse
    source mesh: max `|y|` is about `284 um` for c1 and `192 um` for c2
    versus the intended `500 um` radius.
- Plots after completion:
  - `witness_timing_overview.png`
  - `c0_injection_mesh_xy.png` (nearest saved c0 mesh snapshot; sparse c0 H5 may
    not contain the exact witness injection step)
  - `witness_transverse_motion.png`
  - `near_ip_field_cutoff.png`
  - `track12_figure1_opalx_h5.png`
  - `track12_figure1_opalx_h5_t.png`
- Plot observations:
  - The timing plot shows a large late sampled field spike before c0 retirement.
  - The near-IP cutoff plot confirms zero fields after retirement but is not
    consistent with the simple line-Gaussian decay estimate before retirement.
  - The H5 Figure 1 trajectory panel is qualitatively far from the
    `TestParticleOrbit.dat`/CAIN reference.  This is the next physics/debugging
    item; the current run verifies timing/output plumbing, not final model
    agreement.

2026-07-04 snapshot and step-1/2 status:

- Snapshot before changing witness kinematics:
  `sandbox/Drift-Experiment/snapshots/2026-07-04_before_transverse_witness`
  - Preserves the old staged pair-4 deck and longitudinal witness FROMFILEs.
  - Old witness files had `px = 0`, `py = 0`, `pz = 1.7320511804180809`
    (`gamma = 2`, about `511 keV` kinetic energy).
- Step 1 is now implemented in the working deck:
  - `sandbox/Drift-Experiment/make_transverse_witness_fromfiles.py` generates
    the c1/c2 FROMFILEs.
  - Current case uses `313 keV` transverse kinetic energy:
    `gamma = 1.612525720454`, `beta = 0.784487091427`,
    `p/(mc) = 1.265005612290`.
  - c1/e- starts with `py > 0`, c2/e+ starts with `py < 0`, and both have
    `px = pz = 0` in the input files.
- Step 2 is documented in the deck:
  - `APERTURE = "CIRCLE(0.0002)"` remains the compact c0-source mesh
    (100 um radius, not the 30 cm witness observation domain).
  - `c0_source_mesh_diameter = 2.0e-4` and
    `witness_observation_radius = 0.30` make this convention explicit.
  - The intended implementation direction is source-to-witness field gathering
    from the compact c0 field, not enlarging the field solve mesh to 30 cm.
- Quick verification run:
  - run directory:
    `sandbox/Drift-Experiment/withness_transverse_pair4_step1_2`
  - OPALX command:
    `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
  - Analyzer command:
    `env MPLBACKEND=Agg MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py --run-dir sandbox/Drift-Experiment/withness_transverse_pair4_step1_2`
  - `~/.venv-h6/bin/python -m py_compile` passed for
    `make_transverse_witness_fromfiles.py` and `analyze_withness_timing.py`.
- Transverse smoke result:
  - OPALX completed and wrote c0/c1/c2 H5/stat files.
  - `BB-DIAG` reached active BeamBeam with `c1:n=1,c2:n=1`.
  - First synchronized c0/witness H5 dump remains global step `279`,
    `t = 13.348 ps`, i.e. `5.407 fs` after the analytic pair-4 `T0`.
  - c0 is at `s = 4.001621 mm`, only `+1.621 um` past the IP.
  - H5 witness momenta at first dump:
    - c1: `py = +1.264995`, `pz = -8e-6`
    - c2: `py = -1.264987`, `pz = -1.5e-5`
  - Final written positions in this smoke window are about
    `y = +88.7 um` for c1 and `y = -88.7 um` for c2, still inside the
    current 100 um compact source aperture.
- New analyzer outputs for this run:
  - `timing_mesh_summary.csv`
  - `witness_field_samples.csv`
  - `witness_kinematics_summary.csv`
  - `witness_timing_overview.png`
  - `c0_injection_mesh_xy.png`
  - `witness_transverse_motion.png`

### Order of Magnitude time and space estimates around IP

- Ballistic time for c1/c2 to reach a 30 cm transverse aperture:
  - witness kinetic energy: `313 keV`
  - witness `beta = 0.784487`
  - aperture radius: `0.30 m`
  - `t_hit = 0.30 / (beta_witness * c) = 1.276 ns`
  - with the current `1 fs` fine step this is about `1.28e6` steps
    (`1.28e4` steps with a `100 fs` step)
  - in that time c0 moves downstream by
    `beta_c0 * c * t_hit = 0.382 m = 38.2 cm`
  - total c0-to-witness separation at the 30 cm aperture is roughly
    `sqrt((30 cm)^2 + (38.2 cm)^2) = 48.6 cm`
- Point-charge 1% field scale relative to the injection field:
  - use the initial transverse offset
    `b = sigma_x = 1.944325075701 um`
  - for an inverse-square point-charge estimate,
    `E / E0 = (b / r)^2`
  - `E / E0 = 0.01` gives `r = 10 b = 19.44 um`
  - with ballistic c0/witness motion this occurs after about `50.8 fs`,
    or about `51` steps at `1 fs`
  - at that time c0 is about `15.2 um` downstream of the IP and c1/c2 are
    about `+/-11.9 um` transversely from their injection positions
  - by the time c1/c2 reach the 30 cm aperture, the point-charge field is only
    about `(1.944 um / 48.6 cm)^2 = 1.6e-11` of the injection field
- Consequence for the next implementation step:
  - the physically relevant near-IP interaction is on the scale of tens of
    femtoseconds and tens of microns for this point-charge estimate
  - a full 30 cm ballistic witness flight at `1 fs` is not a practical or
    useful first validation target for the c0 field influence

### Agreed near-IP active-field then ballistic model

- The `51` fine-step cutoff is only a point-charge centroid estimate.
  It is useful as a sanity check, but it is not the right cutoff for the
  present c0 bunch because the deck uses `sigma_z = 0.6 mm`.
- A more relevant first estimate treats c0 near the IP as a long Gaussian
  bunch/line source.  With initial witness offset
  `b = sigma_x = 1.944325075701 um`, a simple estimate is
  `E(t) / E(0) ~= b / sqrt(b^2 + (beta_w c t)^2)
  * exp(-(beta_c0 c t)^2 / (2 sigma_z^2))`.
- For the current `313 keV` transverse witnesses this gives the following
  approximate field-retention times:
  - `10%`: `82 fs`, c0 `24.6 um` downstream, c1/c2 `19.3 um` transverse
  - `5%`: `165 fs`, c0 `49.3 um` downstream, c1/c2 `38.7 um` transverse
  - `2%`: `405 fs`, c0 `121 um` downstream, c1/c2 `95.2 um` transverse
  - `1%`: `768 fs`, c0 `230 um` downstream, c1/c2 `181 um` transverse
- Agreed model for the next implementation:
  - After c1/c2 are injected, c0 remains the only bunch that deposits charge.
  - c1/c2 are passive witness particles: they sample/gather the c0 space-charge
    field but do not contribute to the space-charge deposition.
  - During this near-IP active-field phase, the c0 field mesh must cover the
    c1/c2 witness positions so the field can be sampled at their coordinates.
  - For the first real-bunch cutoff test, use a cutoff around the `1%` estimate,
    i.e. about `0.8 ps` or `800` fine steps at `1 fs`, not the `51` step
    point-charge cutoff.
  - The transverse source-field aperture should therefore cover at least about
    `181 um` witness displacement plus margin; a first low-resolution test can
    use roughly a `250 um` radius compact source/witness-sampling mesh.
  - Once the cutoff is reached, c0 can be ignored for the witness dynamics and
    c1/c2 should be propagated ballistically unless another element acts.
- Validation expectations:
  - Plot OPALX sampled fields at c1/c2 versus the analytic/estimated c0 field
    through the near-IP active window.
  - Plot `E/E0` and mark the cutoff time.
  - Confirm c1/c2 stay inside the field mesh until the cutoff.
  - Confirm c1/c2 do not affect c0 or each other through space-charge
    deposition.

### 2026-07-04 near-IP active-field/ballistic smoke result

- Existing code inspection before the run:
  - BeamBeam self-field deposition uses the source container,
    `bunch.getParticleContainer()`, i.e. c0/container 0.
  - `WITNESS_CONTAINERS = "1,2"` are gathered passively by
    `gatherBeamBeamFieldsToWitnessContainers()`.
  - The passive gather temporarily maps c1/c2 into the c0 source frame, samples
    the current source field, and restores the witness frame.
  - The witnesses do not deposit charge onto the BeamBeam mesh in this path.
- Current deck changes:
  - `near_ip_active_time = 8.0e-13 s`
  - `primary_retire_time = witness_t0 + near_ip_active_time`
  - `fine_stop_s = bb_ip_s + 2.5e-4`
  - `fine_steps = 1020`
  - `c0_source_mesh_diameter = 5.0e-4`
  - `APERTURE = "CIRCLE(0.0005)"`, i.e. about `250 um` radius
- Verification run:
  - run directory:
    `sandbox/Drift-Experiment/withness_near_ip_active_1pct`
  - OPALX command:
    `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
  - Analyzer command:
    `env MPLBACKEND=Agg MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py --run-dir sandbox/Drift-Experiment/withness_near_ip_active_1pct --aperture-radius-m 2.5e-4 --near-ip-active-time-s 8.0e-13`
  - `~/.venv-h6/bin/python -m py_compile` passed for
    `sandbox/Drift-Experiment/analyze_withness_timing.py`.
- Run result:
  - OPALX completed cleanly.
  - `BB-DIAG` sequence:
    - inactive, witnesses empty
    - active, witnesses empty
    - active, c1/c2 each have one particle
    - completed with `source_retirement_pending=TRUE`
    - completed with `retired_bunches=1`, `source_active=FALSE`, c1/c2 active
  - First synchronized c0/witness H5 dump remains:
    - global step `279`
    - `t = 13.348 ps`
    - `t - witness_t0 = 5.407 fs`
    - c0 `s = 4.001621 mm`, `+1.621 um` past IP
  - Configured c0 retirement:
    - `retire_time = 14.142592708656 ps`
    - `retire_time - witness_t0 = 800 fs`
  - Final written witness H5 sample:
    - `t = 14.220 ps`
    - c1 `y = +206.3 um`
    - c2 `y = -206.3 um`
    - both remain within the `250 um` source/witness-sampling aperture
  - `near_ip_field_cutoff.csv` shows the first post-retirement H5 field samples
    at `14.148 ps` have `E_abs = 0`, consistent with ballistic continuation
    after c0 retirement.
- Generated plots:
  - `witness_timing_overview.png`
  - `c0_injection_mesh_xy.png`
  - `witness_transverse_motion.png`
  - `near_ip_field_cutoff.png`
- Current caveat:
  - With `NXY = 16` and a `250 um` radius, the transverse cell size is
    `31.25 um`, much larger than `sigma_xy = 1.944 um`.
  - This run verifies timing, retirement, passive witness behavior, and mesh
    containment.  It is not yet a quantitative field-accuracy validation.

Current implementation:

- The deck now models a 1 mm upstream drift followed by a 6 mm BeamBeam element:
  - `drift_length = 1.0e-3 m`
  - `bb_length = 6.0e-3 m`
  - `bb_ip_s = 4.0e-3 m`
  - BeamBeam window spans `1.0 mm` to `7.0 mm`
- Uses the PPTX/track12 pair-4 timing first:
  - `witness_ct_m = 0`
  - `witness_t0 = bb_ip_s / (primary_beta * CLIGHT)`
  - numeric value from the analyzer: `13.342592708656 ps`
- Uses no mirrored source:
  - `COPY_TIME = 0.0`, interpreted by OPALX as `copy_time=NONE`
  - `WITNESS_CONTAINERS = "1,2"`
- The smoke deck is intentionally cheap:
  - `primary_macroparticles = 100`
  - `primary_charge_scale = 1.0e-5`
  - `NXY = 16`, `NZ = 32`
  - `APERTURE = "CIRCLE(0.0005)"` gives a 250 um transverse half-width
  - `RETIRE_TIME = witness_t0 + 0.8 ps` retires c0 at the near-IP cutoff
- The witness FROMFILEs now use the 313 keV transverse pair-4 setup described
  above, generated by `make_transverse_witness_fromfiles.py`.

Important fix in the input:

- The first staged attempt used `coarse_steps = 200` with
  `fine_start_s = 3.95 mm`.  OPALX reached the first `ZSTOP` near `13.2 ps`
  but continued advancing time to `20 ps` before starting the fine segment.
  That caused the witnesses to emit while c0 was effectively parked upstream of
  the IP.
- The deck now uses:
  - `coarse_steps = 132`, because
    `fine_start_s / (primary_beta * CLIGHT * coarse_dt) = 131.76`
  - `fine_steps = 520`, enough for c0 to cross from `3.95 mm` to `4.10 mm`
    with `fine_dt = 1 fs`

C++ diagnostic/transition fix used for this run:

- `ParallelTracker::performBeamBeamWindowEntryTransition()` now skips the
  copied-source charge-conservation pre-solve when
  `BEAMBEAM::copyTimeReached(...)` is false.
- Rationale: with `COPY_TIME = 0.0`, OPALX intentionally has no copied source,
  so requiring copied-source deposited-charge diagnostics aborts the no-copy
  BeamBeam timing case before witness injection.

Completed smoke run:

```sh
cmake --build build_openmp -j 8 --target opalx_exe
mkdir -p sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed
cp sandbox/Drift-Experiment/spacecharge_drift_withness.in \
  sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/
cp sandbox/Drift-Experiment/track12_pair4_*.fromfile \
  sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed/
env OMP_NUM_THREADS=8 \
  ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1
env MPLBACKEND=Agg \
  MPLCONFIGDIR=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.matplotlib \
  XDG_CACHE_HOME=/Users/adelmann/git/opalx-beambeam/sandbox/Drift-Experiment/.cache \
  ~/.venv-h6/bin/python sandbox/Drift-Experiment/analyze_withness_timing.py
```

Run output:

- `sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed`
- `timing_mesh_summary.csv`
- `witness_field_samples.csv`
- `witness_timing_overview.png`
- `c0_injection_mesh_xy.png`
- c0/c1/c2 H5 and stat files

Timing result:

- `BB-DIAG` sequence:
  - inactive, witnesses empty
  - active, witnesses empty
  - active, c1/c2 each have one particle
  - active with source inactive after c0 leaves the source-active interval
- First synchronized c0/witness H5 dump:
  - global step `279`
  - `t = 13.348 ps`
  - `t - witness_t0 = 5.407 fs`
  - c0 `s = 4.001621 mm`
  - c0 is only `+1.621 um` past the IP
- This establishes the timing convention for pair 4:
  `witness_t0 = (bb_ip_s - witness_ct_m) / (primary_beta * CLIGHT)`.

Mesh result:

- Longitudinal containment is fine for this setup:
  - at first witness dump, c0 center is `4.001621 mm`
  - c0 longitudinal half-extent from H5 is about `1.633 mm`
  - occupied c0 span is approximately `2.37 mm` to `5.63 mm`
  - BeamBeam window is `1.0 mm` to `7.0 mm`
- Transverse containment is safe but field sampling is very coarse:
  - aperture half-width = `100 um`
  - `NXY = 16` gives transverse cell size `12.5 um`
  - c0 rms sizes at injection are `2.125 um` and `1.841 um`
  - c0 max transverse particle extents are about `4.87 um` and `5.55 um`
  - the bunch is well inside the aperture but mostly inside one central cell
- Prior notes in `sandbox/track12particles/opalx/beambeam_window_scan/README.md`
  show aperture/cell-size tuning alone can produce large sign and magnitude
  pathologies.  The next quantitative field comparison should therefore scan
  mesh resolution deliberately rather than simply shrinking the aperture.

Changed/added files for this task:

- `sandbox/Drift-Experiment/spacecharge_drift_withness.in`
- `sandbox/Drift-Experiment/analyze_withness_timing.py`
- `sandbox/Drift-Experiment/make_transverse_witness_fromfiles.py`
- `sandbox/Drift-Experiment/track12_pair4_electron.fromfile`
- `sandbox/Drift-Experiment/track12_pair4_positron.fromfile`
- `src/Algorithms/ParallelTracker.cpp`
- `src/Algorithms/ParallelTracker.h`
- `HANDOFF.md`
- `sandbox/Drift-Experiment/TASK_STATE.md`

Next steps:

- Use the fixed staged deck as the baseline for the next BeamBeam field checks.
- Add a small mesh scan around the witness location.  The safe starting point is
  to keep the 6 mm longitudinal BeamBeam window, but increase transverse
  resolution and/or reduce aperture with care.  The current low-resolution
  `100 um / 16` case has `12.5 um` cells, much larger than `sigma_x`.
- Once the timing/mesh mechanics are stable, increase particles/charge toward
  the track12 settings and compare the c0 field sampled at c1/c2 against the
  smooth analytic/manufactured reference.

## Recent Task: Drift Space-Charge 30 cm Redo

Goal:

- Re-establish the `sandbox/Drift-Experiment/spacecharge_drift_30cm.in`
  validation after reverting the last `origin/master` merge.
- Use only H5 particle/field output for the OPALX diagnostic path.
- Redo the full grid and particles-per-cell matrix:
  - `N_GRID = 16, 32, 64`
  - `NPPG = 5, 8, 14`
  - `N_PARTICLES = N_GRID * N_GRID * N_GRID * NPPG`
- Produce the same radial field comparison plots used in the later H5 drift
  checks.

Current repository state:

- Branch: `271-implement-interation-point-element`
- The branch is ahead of `origin/271-implement-interation-point-element`.
- The previous `origin/master` merge was reverted in commit `1177cfa4a`
  (`Revert origin/master merge`).
- The local IPPL checkout in `build_openmp/_deps/ippl-src` was reset to
  `27d11d0b58bc5db6a582f8132c3d6c88285b26bd` so this branch builds against the
  expected `solver_recv` API.
- A saved local patch from before aborting/reverting the merge remains at
  `/tmp/opalx-beambeam-local-before-merge-abort.patch`.

Changed files for the current task:

- `src/Algorithms/ParallelTracker.cpp`
- `src/Algorithms/ParallelTracker.h`
- `sandbox/Drift-Experiment/spacecharge_drift_30cm.in`
- `sandbox/Drift-Experiment/run_spacecharge_convergence.py`
- `sandbox/Drift-Experiment/TASK_STATE.md`
- `HANDOFF.md`

Implementation decisions:

- Kept the reverted branch's 30 cm drift field-solver settings, including
  `PARFFTZ = FALSE`.
- Enabled `OPTION, EBDUMP = TRUE` in
  `sandbox/Drift-Experiment/spacecharge_drift_30cm.in`.
- Restored an opt-in H5 diagnostic in `ParallelTracker`:
  - Environment variable: `OPALX_SC_FIELD_H5_STEPS`
  - Current runner uses `OPALX_SC_FIELD_H5_STEPS=0`.
  - The diagnostic writes the primary bunch after self-field computation and
    transformation back to the reference frame, before `bunchUpdate()`.
- `run_spacecharge_convergence.py` now:
  - reads OPALX H5 particle dumps directly,
  - computes the analytic isotropic Gaussian field at the same dumped particle
    coordinates,
  - writes metrics, sampled particle comparisons, per-case radial profiles, and
    PNG plots,
  - deletes large raw H5 dumps after each successful case unless `--keep-raw`
    is passed.
- The comparator helper `sandbox/compare-e-fields/compare_gaussian_pic_fields.py`
  was not present after the merge revert and is not tracked by Git. The
  convergence runner now contains a local fallback implementation of the
  Gaussian field and metric routines, preserving the previous metrics schema.

Commands run for the current task:

```sh
~/.venv-h6/bin/python -m py_compile sandbox/Drift-Experiment/run_spacecharge_convergence.py
cmake --build build_openmp -j 8 --target opalx_exe
~/.venv-h6/bin/python sandbox/Drift-Experiment/run_spacecharge_convergence.py \
  --output-dir sandbox/Drift-Experiment/redone_full_grid_seed42 \
  --n-grid-values 16,32,64 \
  --nppg-values 5,8,14 \
  --seeds 42 \
  --threads 8 \
  --force \
  --sample-particles 50000
```

Run status:

- Completed the 3x3 grid/NPPG matrix for seed 42.
- This turn did not rerun the earlier three-seed/27-run ensemble; disk space was
  about 1.6 GiB free, and the user requested the grid and particle-number grid.
- Raw H5 dumps were deleted after extraction. No `.h5` files remain in the new
  output directory.
- Output directory:
  `sandbox/Drift-Experiment/redone_full_grid_seed42`

Relative vector L2 vs analytic for the redone seed-42 matrix:

| N_GRID | NPPG | N_PARTICLES | relative vector L2 |
| ---: | ---: | ---: | ---: |
| 16 | 5 | 20480 | 0.107615 |
| 16 | 8 | 32768 | 0.102523 |
| 16 | 14 | 57344 | 0.106806 |
| 32 | 5 | 163840 | 0.026637 |
| 32 | 8 | 262144 | 0.022694 |
| 32 | 14 | 458752 | 0.023855 |
| 64 | 5 | 1310720 | 0.006748 |
| 64 | 8 | 2097152 | 0.005086 |
| 64 | 14 | 3670016 | 0.004180 |

Generated plots:

- `sandbox/Drift-Experiment/redone_full_grid_seed42/relative_l2_vs_grid.png`
- `sandbox/Drift-Experiment/redone_full_grid_seed42/relative_l2_vs_nppg.png`
- `sandbox/Drift-Experiment/redone_full_grid_seed42/relative_l2_vs_particles.png`
- `sandbox/Drift-Experiment/redone_full_grid_seed42/radial_field_vs_analytic_all_cases.png`
- `sandbox/Drift-Experiment/redone_full_grid_seed42/signed_radial_field_vs_analytic_all_cases.png`
- `sandbox/Drift-Experiment/redone_full_grid_seed42/radial_field_relative_error_all_cases.png`
- Each per-case directory also contains:
  - `metrics.csv`
  - `radial_profile.csv`
  - `particle_field_comparison_sample.csv`
  - `radial_field_vs_analytic.png`
  - `signed_radial_field_vs_analytic.png`
  - `radial_field_relative_error.png`

Remaining risks / next steps:

- The redo was single-rank/OpenMP only. Repeat at MPI multi-rank if this
  diagnostic becomes an acceptance test.
- If the user wants the earlier statistical ensemble again, rerun with
  `--seeds 42,43,44`, but watch disk space closely.
- The H5 diagnostic hook is opt-in and intended for validation; it should not
  affect normal runs unless `OPALX_SC_FIELD_H5_STEPS` is set.

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

## 2026-07-04 Release witness T0=4 ps / 1600 fs active window

Goal:

- Rebuild OPALX in Release mode.
- Start the pair-4 witness beams at absolute `T0 = 4 ps`.
- Keep the c0 BeamBeam source active for `1600 fs` after witness `T0`.

Build:

- Reconfigured with both build-type knobs because
  `cmake/OPALXOptions.cmake` derives `CMAKE_BUILD_TYPE` from `BUILD_TYPE`:
  `cmake -S . -B build_openmp -DBUILD_TYPE=Release -DCMAKE_BUILD_TYPE=Release`
- Verified cache:
  - `BUILD_TYPE:STRING=Release`
  - `CMAKE_BUILD_TYPE:STRING=Release`
  - `CMAKE_CXX_FLAGS_RELEASE:STRING=-O3 -DNDEBUG`
- Rebuilt target:
  `cmake --build build_openmp -j 8 --target opalx_exe`

Input/deck changes:

- Updated `sandbox/Drift-Experiment/spacecharge_drift_withness.in`:
  - `witness_t0 = 4.0e-12`
  - `near_ip_active_time = 1.6e-12`
  - `primary_retire_time = witness_t0 + near_ip_active_time`
  - `APERTURE = "CIRCLE(0.0010)"`, i.e. 500 um radius
  - `coarse_steps = 40`
  - `fine_steps = 1780`
- Added shifted primary source timing:
  - `primary_source_r0z = bb_ip_s - primary_beta * CLIGHT * witness_t0 - witness_ct_m`
  - `SOURCE_PRIMARY_ELECTRONS` now uses `R0Z = primary_source_r0z`
  - `fine_start_s`/`fine_stop_s` are absolute lattice positions computed from
    the shifted c0 source.
- Reason: the first T0=4 ps attempt left c0 launched at `R0Z=0`; c0 had not
  entered the BeamBeam source-active window when witnesses were injected, so
  `BB-DIAG` never reached the intended active-with-witness state.

Post-processing changes:

- Updated `sandbox/Drift-Experiment/analyze_withness_timing.py` to infer the
  shifted `primary_source_r0z` and plot/report absolute c0 lattice position as
  `absolute_spos_m`.
- The H5 `SPOS` remains the traveled path length; for shifted-source decks the
  absolute c0 position is `SPOS + primary_source_r0z`.

Corrected run:

- Directory:
  `sandbox/Drift-Experiment/withness_t0_4ps_active_1600fs_release_shifted_c0`
- Command:
  `env OMP_NUM_THREADS=8 ../../../build_openmp/src/opalx spacecharge_drift_withness.in > run.log 2>&1`
- Wall-clock by file timestamps: about 1 minute.
- BeamBeam diagnostics:
  - active before witnesses are populated
  - active with `c1:active:n=1,c2:active:n=1`
  - completed with source retirement pending
  - completed with c0 source retired

Timing summary:

- `witness_t0 = 4.0 ps`
- first witness H5 dump: `4.020 ps`, `20.0 fs` after T0
- configured c0 retire time: `5.600 ps`, `1600 fs` after T0
- inferred `primary_source_r0z = 2.800833 mm`
- nearest c0/IP H5 dump: `t = 4.000 ps`, `s = 4.000 mm`
- at first witness H5 dump, c0 absolute `s = 4.005996 mm`, i.e.
  `+5.996 um` downstream of the IP
- c1/c2 final H5 dump: `5.780 ps`
- sampled witness fields are zero from `5.620 ps` onward, after c0 retirement.

Plots/results:

- `witness_timing_overview.png`
- `c0_injection_mesh_xy.png`
- `witness_transverse_motion.png`
- `near_ip_field_cutoff.png`
- CSV summaries:
  - `timing_mesh_summary.csv`
  - `witness_kinematics_summary.csv`
  - `witness_field_samples.csv`
  - `near_ip_field_cutoff.csv`

Remaining caveats:

- This is still a low-resolution timing smoke case:
  `NXY=16`, `NZ=32`, `primary_macroparticles=100`,
  `primary_charge_scale=1e-5`.
- The near-IP field samples are noisy and should not be interpreted as a
  quantitative physics validation until the source charge/macroparticle count
  and mesh scan are restored.

## 2026-07-05 one-c0 track12 400k long active-window run

Goal:

- Run a 400k c0-source case for a physically sensible active time, keeping c0
  alive long enough for the track12 witnesses to sample the late source field.
- Use one physical c0 source with `COPY_TIME = 0`, raw
  `track12_electrons.fromfile` / `track12_positrons.fromfile`, and H5-only
  post-processing.

Input/deck changes:

- Updated `sandbox/Drift-Experiment/spacecharge_drift_withness.in` so BeamBeam
  geometry no longer forces premature completion:
  - `bb_start_s = 1.0e-3`
  - `bb_ip_s = 6.0e-3`
  - `bb_length = 1.9e-2`, i.e. BeamBeam range is about `1..20 mm`
  - `IP_S` remains fixed at `6 mm`; it is no longer tied to the element
    midpoint.
- Set `OPTION, BOUNDPDESTROY = 1.0e9` for this diagnostic.  With the default
  value `10`, OPALX started deleting c0 source macroparticles around `40 ps`,
  which changes the source charge before the requested `RETIRE_TIME`.

Aborted diagnostics:

- `one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs`
  - Had `near_ip_active_time = 50 ps`, but BeamBeam range was still `1..11 mm`.
  - BeamBeam completed before the requested retire time, so the run was stopped.
- `one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs_extended_bb`
  - Fixed the range to `1..20 mm`.
  - Still had `BOUNDPDESTROY = 10`; OPALX started deleting c0 particles around
    `40 ps`, so the run was stopped.

Completed clean run:

- Directory:
  `sandbox/Drift-Experiment/one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs_extended_bb_nobound`
- Command:
  `OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads /usr/bin/time -p ../../../build_openmp/src/opalx --info 2 spacecharge_drift_withness.in > spacecharge_drift_withness.out 2> spacecharge_drift_withness.err`
- Key settings:
  - `primary_macroparticles = 400000`
  - `near_ip_active_time = 5.0e-11`
  - `post_retire_observation_time = 1.0e-12`
  - `fine_dt = 5.0e-15`
  - `fine_steps = 10210`
  - `C0PSDUMPFREQ = 600`
- Runtime from `/usr/bin/time`: `real 1589.40 s` (`26.49 min`),
  `user 10382.11 s`, `sys 159.80 s`.
- Output sizes:
  - `spacecharge_drift_withness_c0.h5`: `733 MB`
  - `spacecharge_drift_withness_c1.h5`: `5.3 MB`
  - `spacecharge_drift_withness_c2.h5`: `5.3 MB`
  - `spacecharge_drift_withness.out`: `2.7 MB`
- `spacecharge_drift_withness.err` contains only the `/usr/bin/time` output.

Run diagnostics:

- BeamBeam range printed at startup:
  `s_range=(0.001, 0.020) m`, `interaction_point_s=0.006 m`.
- No `Marked` / `Deleted` boundary-clipping messages appeared before retirement.
- At `t = 54.000 ps`, c0 was still present at `s = 16.189 mm`.
- At `t = 54.005 ps`, BeamBeam completed with
  `source_retirement_pending=TRUE`, then retired all 400000 c0 source
  particles.  c1/c2 remained active through the post-retire observation.
- Final H5 dump time: `55.05 ps`.

Plots/results:

- `witness_timing_overview.png`
- `c0_injection_mesh_xy.png`
- `witness_transverse_motion.png`
- `near_ip_field_cutoff.png`
- `track12_figure1_opalx_h5.png`
- `track12_figure1_opalx_h5_t.png`
- CSV summaries:
  - `timing_mesh_summary.csv`
  - `witness_kinematics_summary.csv`
  - `witness_field_samples.csv`
  - `near_ip_field_cutoff.csv`

Remaining caveats:

- The run now verifies timing, source retention, and time-based retirement, but
  it does not reproduce the CAIN track12 Figure 1 qualitatively.
- The sampled near-IP field stays large until c0 retirement instead of following
  the simple line-Gaussian falloff estimate.
- The witness trajectories reach the compact `500 um` source-sampling aperture
  (`max_abs_y_over_aperture` is about `0.996` for c1 and `1.0002` for c2), so
  aperture/mesh extent is now an active modeling parameter, not just a
  numerical detail.

## 2026-07-05 analytic Gaussian witness model sandbox

Goal:

- Switch back to the analytic boosted-Gaussian witness model, but keep the same
  timing and witness setup as the completed OPALX run above.
- Use the same six raw `track12_electrons.fromfile` and
  `track12_positrons.fromfile` witnesses, same c0 charge/size, same witness
  bunch charge magnitude, same IP, same c0 edge-at-first-witness timing, and
  same c0 retirement time.

New sandbox:

- Directory: `sandbox/analytic-model`
- Copied provenance/input files:
  - `spacecharge_drift_withness.opalx-reference.in`
  - `track12_electrons.fromfile`
  - `track12_positrons.fromfile`
- Added reproducible driver:
  - `sandbox/analytic-model/run_analytic_witness.py`
  - Uses the triaxial rigid boosted Gaussian evaluator and Boris pusher from
    `sandbox/track12particles/track12particles.py`.
  - Advances the witnesses with `dt = 5 fs`, samples exactly at the OPALX H5
    witness dump times, and compares against the clean OPALX c1/c2 H5 files.
- Added documentation:
  - `sandbox/analytic-model/README.md`

Default analytic setup:

- c0 is a physical electron bunch with charge `-1.25e10 e =
  -2.0027207925e-09 C`.
- c0 kinetic energy is `245 MeV`, `beta = 0.999997833949`.
- c0 Gaussian sigmas:
  - `sigma_x = sigma_y = 1.944325075701 um`
  - `sigma_z = 0.6 mm`
- Witnesses:
  - 6 electrons and 6 positrons from the OPALX FROMFILEs.
  - `pz = sqrt(3)`, `px = py = 0`, matching the raw track12/PPTX convention.
  - bunch charge magnitude per witness species is `6 e =
    9.613059804e-19 C`.
- Timing/geometry:
  - witness `T0 = 4 ps`
  - BeamBeam/IP plane `S_IP = 6 mm`
  - first witness is at `S = 5.1 mm`
  - c0 centroid at `T0` is `3.3 mm`
  - c0 `+3 sigma_z` edge at `T0` is `5.1 mm`
  - c0 field is active until `54 ps`, then zeroed to mimic c0 retirement.

Generated default outputs:

- `setup_summary.json`
- `witness_initial_conditions.csv`
- `analytic_witness_trajectory.csv`
- `opalx_witness_trajectory_sampled.csv`
- `analytic_vs_opalx_samples.csv`
- `analytic_vs_opalx_summary.csv`
- `analytic_vs_opalx_side_by_side_x_vs_s.png`
- `analytic_vs_opalx_side_by_side_x_vs_t.png`
- `analytic_vs_opalx_overlay_x_vs_s.png`
- `analytic_vs_opalx_overlay_x_vs_t.png`
- `analytic_field_timing.png`
- `analytic_minus_opalx_differences.png`

First diagnostic observations:

- The default physical negative-c0 analytic field has the expected sign at the
  first saved witness sample.  For electron pair 1 at `t = 4.1 ps`:
  - analytic `Ex = -5.4949567e7 V/m`
  - OPALX H5 `Ex = +8.8042964e7 V/m`
- OPALX c1 then kicks in the direction implied by the positive H5 field, not by
  the physical negative-electron source field.
- A sign-flipped diagnostic was also run in
  `sandbox/analytic-model/effective-positive-c0` using
  `--source-charge-scale -1`.
  - This changes the first-sample analytic `Ex` sign, but does not remove the
    order-one trajectory/field differences.
  - Therefore the current OPALX-vs-analytic discrepancy is not explained by a
    single global source-charge sign flip alone.
- Default comparison summary highlights:
  - `relative_l2_difference(x_um) = 1.4705`
  - `relative_l2_difference(s_minus_ip_mm) = 1.6526`
  - `relative_l2_difference(E_abs_V_per_m) = 0.9934`

Run commands:

- Physical c0:
  `/Users/adelmann/.venv-h6/bin/python sandbox/analytic-model/run_analytic_witness.py`
- Sign diagnostic:
  `/Users/adelmann/.venv-h6/bin/python sandbox/analytic-model/run_analytic_witness.py --output-dir sandbox/analytic-model/effective-positive-c0 --source-charge-scale -1`
