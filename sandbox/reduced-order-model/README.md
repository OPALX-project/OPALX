# BeamBeam reduced-order model

This directory is the reproducible workspace for comparing a rigid two-Gaussian
source model with OPALX witness-particle trajectories.

- `python/` contains the independent reference-field and trajectory scripts.
- `opalx/` contains OPALX input decks and their local input data.

The OPALX `BEAMBEAM` element provides the model switch:

```opal
IP: BEAMBEAM, L=bb_length, BBRIGID=TRUE;
```

`BBRIGID=FALSE` is the default and preserves the self-consistent source response.
With `BBRIGID=TRUE`, source container 0 still deposits charge and its solved mesh
field is gathered to configured witness containers, but that BeamBeam collective
field is removed from the source particles before their momentum update. The
source still drifts normally and still responds to external fields.

Consequently, an unchanged moving Gaussian is obtained only when the OPALX deck
also gives the primary source zero divergence and momentum spread and does not
apply external focusing inside the BeamBeam element. `BBRIGID` does not silently
overwrite those input choices.

The first comparison should use identical source geometry, charge, clock, pair
emission times, and initial pair phase space in both implementations. Compare the
electric and magnetic fields first, then the electron/positron trajectories.

## Stage 1: rigid-source fields

The authoritative parameters are in `parameters.json`. The initial field study
uses two identical 245 MeV electron bunches, each containing 1.25e10 electrons,
with lab-frame rms sizes

```text
sigma_x = sigma_y = 1.944325075701 um
sigma_z = 0.6 mm.
```

At every snapshot the IP is at the origin. The centroids are at z=-d/2 and
z=+d/2, moving in +z and -z respectively, for d/sigma_z = 3, 2, 1, and 0. The
Python evaluator uses the exact uniform-motion Lorentz transform of each rigid
Gaussian. Pair fields are not included.
