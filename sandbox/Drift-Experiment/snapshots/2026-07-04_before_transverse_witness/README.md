# Snapshot: before transverse-witness kinematics

Date: 2026-07-04

Purpose:

- Preserve the staged pair-4 timing baseline before changing c1/c2 from the
  old track12 longitudinal test-particle kinematics to transverse witness
  kinematics.

Contents:

- `spacecharge_drift_withness.in`: staged 1 mm drift + 6 mm BeamBeam deck after
  the timing fix (`coarse_steps = 132`, `fine_steps = 520`).
- `track12_pair4_electron.fromfile`
- `track12_pair4_positron.fromfile`
- `timing_mesh_summary.csv`

Baseline state:

- Witness files still had `px = 0`, `py = 0`, `pz = 1.7320511804180809`.
- That corresponds to a longitudinal normalized momentum with `gamma = 2` and
  kinetic energy about `511 keV`.
- This is not the intended next model.  The next model should put c1/c2 in the
  transverse plane with about `313 keV` kinetic energy.

Timing result preserved by this snapshot:

- `witness_t0 = 13.342592708656 ps`
- first synchronized c0/c1/c2 H5 sample at `13.348 ps`
- c0 at `s = 4.001621 mm`, about `1.621 um` past the IP
