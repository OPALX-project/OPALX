# OPALX 12-Particle BeamBeam Experiment

This directory contains the OPALX version of the CAIN 12-particle comparison.

## Files

- `generate_track12_opalx_inputs.py` converts `sandbox/TestParticleOrbit.dat`
  into OPALX `FROMFILE` witness distributions.
- `track12_electrons.fromfile` and `track12_positrons.fromfile` are generated
  inputs with columns `x y z px py pz`.
- `track12_witness_metadata.csv` records the particle names and source-file
  timestamp, since OPALX `FROMFILE` does not carry particle labels.
- `track12_beambeam.in` is the intended comparison input with 400000 primary
  source macroparticles.
- `track12_beambeam_smoke.in` is a cheaper parser/tracking smoke run with 4000
  primary source macroparticles and a coarser time step.
- `track12_one_bunch_electrons.in` is the first validation run: one read-in
  electron bunch only, no BeamBeam interaction.
- `track12_primary_one_beam_smoke.in` tracks one 4000-macroparticle primary
  Gaussian electron beam in a drift, with no BeamBeam interaction.
- `track12_primary_one_beam.in` is the same source-only run at 400000
  macroparticles.
- `track12_primary_with_pairs_smoke.in` tracks the primary beam plus the read-in
  six-electron and six-positron CAIN test particles in a drift, still without a
  BeamBeam element.
- `track12_beambeam_one_source_quick.in` is a fast one-source BeamBeam
  diagnostic. It enables `WITNESS_CONTAINERS` but leaves `COPY_TIME = 0`, so no
  copied counter-propagating source is included.
- `track12_beambeam_copied_source_1fs_probe.in` is a staged-timestep diagnostic
  with a 1 fs fine segment through the copied-source BeamBeam window.
- `generate_timing_pair_inputs.py` creates the pair-wise timing-correct OPALX
  decks under `timing_pairs/`. These decks keep the numeric TestParticleOrbit
  witness coordinates but shift the primary source centroid per pair so the
  first field evaluation uses the insertion timing implied by
  `TestParticleOrbitSimulation.pptx`.

## Coordinate Convention

The CAIN columns are mapped as

```text
x_OPALX = x_CAIN
y_OPALX = y_CAIN
z_OPALX = s_CAIN
px_OPALX = Px_CAIN / (m_e c)
py_OPALX = Py_CAIN / (m_e c)
pz_OPALX = Ps_CAIN / (m_e c)
```

The witness emission sources apply `R0Z = bb_ip_s`, so the CAIN `s` coordinate
is interpreted as a longitudinal offset around the BeamBeam interaction point.

The PowerPoint source slides define the witness `ct` values as insertion times.
For an oncoming primary bunch centered at the IP at `ct = 0`, the source center
at insertion is `s_c = -ct_i`. Because the TestParticleOrbit initial rows use
`s_i = ct_i`, each pair needs a different primary source offset. The generated
pair decks use

```text
primary_source_r0z = bb_ip_s - ct_i
witness_r0z = bb_ip_s
```

so `s_i - s_c = 2 ct_i`; pair 1 starts at the `-3 sigma_z` source edge and
pair 4 starts at the bunch center. This fixes timing only; it intentionally does
not tune source direction, species mapping, field normalization, or pusher
details.

## Commands

Regenerate witness files:

```sh
~/.venv-h6/bin/python generate_track12_opalx_inputs.py
```

Generate pair-wise slide-timed BeamBeam inputs:

```sh
~/.venv-h6/bin/python generate_timing_pair_inputs.py
```

Run the smoke input from this directory:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_one_bunch_electrons.in
```

Run the source-only primary-beam smoke input:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_primary_one_beam_smoke.in
```

Run the primary-plus-e-/e+ no-BeamBeam smoke input:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_primary_with_pairs_smoke.in
```

Run the BeamBeam smoke input from this directory:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_beambeam_smoke.in
```

Run the 400000-particle comparison:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_beambeam.in
```
