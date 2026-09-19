# BeamBeam input inventory

This directory contains OPALX inputs for the gamma-gamma pair tracking studies.
The current production target is:

- BeamBeam window length: `bb_length = 0.32 m`
- BeamBeam cylinder aperture: `APERTURE="CIRCLE(0.30)"`, i.e. radius `15 cm`
- Interaction point: the BeamBeam element midpoint, 16 cm from either edge
- Witness population: 1,297 CAIN electrons and 1,297 CAIN positrons
- Primary source: anisotropic lab-frame Gaussian with
  `sigma_x = sigma_y = 1.944325075701 um` and `sigma_z = 0.6 mm`
- Witness injection: `tinj = primary_ip_time - 16.747e-12`
- Primary source retirement: `primary_retire_time = 1000e-12`
- Coarse debug timestep unless noted: `DT = 1.0e-12`

The `1000 ps` retirement is the present explicit deck choice, not a claimed
three-sigma physical cutoff. Copy and retirement times require convergence
checks for production conclusions. See
[`../BEAMBEAM_PHYSICS_AND_VALIDATION.md`](../BEAMBEAM_PHYSICS_AND_VALIDATION.md)
for the model contract.

## Current large-cylinder inputs

| file | purpose |
| --- | --- |
| `gamma_gamma_pairs-large-cylinder.in` | Nominal primary charge, 32 cm BeamBeam window, 15 cm cylinder radius, coarse `DT = 1 ps`, source retirement at `1000 ps`. Use this for visual field debugging and nominal crossing/histogram studies when staged-DT refinement is not needed. |
| `gamma_gamma_pairs-large-cylinder-q1em5.in` | Same as the nominal large-cylinder case, but with `primary_charge_scale = 1.0e-5`. Use for reduced-primary-charge comparisons against the nominal case. |
| `gamma_gamma_pairs-large-cylinder-retire1000ps.in` | Explicit-name alias of the nominal `1000 ps` retirement coarse case. Kept because existing EF/H5 debug artifacts and scripts use this stem. |
| `gamma_gamma_pairs-large-cylinder-retire1000ps-q1em5.in` | Explicit-name alias of the reduced-charge `1000 ps` retirement coarse case. Kept for existing comparison artifacts and scripts. |
| `gamma_gamma_pairs-large-cylinder-retire1000ps-omp2.in` | OpenMP benchmark copy of the nominal `1000 ps` coarse case. Run with `OMP_NUM_THREADS=2` so output stem does not overwrite the nominal plot run. |
| `gamma_gamma_pairs-large-cylinder-retire1000ps-omp4.in` | OpenMP benchmark copy of the nominal `1000 ps` coarse case. Run with `OMP_NUM_THREADS=4` so output stem does not overwrite the nominal plot run. |
| `gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt.in` | Production-style large-cylinder run with `1000 ps` source retirement and staged timesteps around the overlap region. Use this when the witness beam dynamics need the smaller overlap timestep. |
| `gamma_gamma_pairs-staged-dt-smoke.in` | Small/cheap staged-DT parser and tracking smoke test. It uses the large-cylinder geometry and `1000 ps` retirement but is intended only to verify timestep-vector mechanics. |

## Legacy/proof-of-principle inputs

| file | purpose |
| --- | --- |
| `gamma_gamma_pairs-2.in` | Historical 5 cm proof-of-principle cylinder with `CIRCLE(0.02)`, `121 ps` retirement, and a coarse 1 ps step. Do not use for current validation or production conclusions. |
| `gamma_gamma_pairs-3.in` | Historical 5 cm rectangular-aperture setup. Do not use for current validation or production conclusions. |

## Generated data note

The `data/` subdirectory and local `*.stat`/`*.h5` outputs are run artifacts.
The most recent post-overlap field-debug products were generated from
`gamma_gamma_pairs-large-cylinder-retire1000ps.in` and converted with
`sandbox/python/convert_efield_dumps_to_h5.py`.
