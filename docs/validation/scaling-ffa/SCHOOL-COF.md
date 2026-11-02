# Scaling FFA closed-orbit finder and the school ring

`COF` now accepts `SCALINGFFAMAGNET` in static magnetic rings. It evaluates the
existing FFA fields through the shared scalar ray tracker, including overlapping
Enge fringes. It also checks radial and vertical bounds independently of field
selection within each magnet's spiral angular window. A ray outside the aperture
must fail explicitly, rather than silently receiving zero field and continuing.
FFAs do not enter the conventional nearest-body aperture selection, which would
apply an inappropriate straight or opposite-arm window.

The existing four-dimensional damped Newton solver, return-plane localization,
integrators and nominal geometry closure checks are unchanged. Existing solver
and test tolerances are unchanged.

## Result at fixed 12 MeV

For the legacy-compatible 16-cell school field, in the section through the
positive radial axis with forward direction +Z:

| Quantity | Value |
|---|---:|
| Kinetic energy | 12 MeV |
| Radius at launch | **4.392925555470 m** |
| OPALX section x (relative to R0=4.42 m) | -0.0270744445304 m |
| Radial mechanical momentum px/(mc) | **0.0162803888787** |
| Tangential mechanical momentum pz/(mc) | 0.1596165907281 |
| Vertical y, py/(mc) | 0, 0 |
| Return time | approximately 5.901895189e-7 s |
| Closed-orbit path length | approximately 28.02969936 m |

The DOP853 result has a fresh full-turn position residual of 3.24e-14 m,
momentum residual below 2.51e-15 in p/(mc), and relative kinetic-energy drift
below 1.74e-14. These are numerical closure diagnostics, not an uncertainty
estimate; agreement with independent equations and timestep checks is reported
below.

`school-closed-orbit.json` contains this full result. `school-16-cell.in` is a
self-contained input from the displaced-start regression: starting from the
notebook's stored radius/momentum plus a 100 µm vertical displacement, the
native four-dimensional COF converged in two Newton corrections (28 map
evaluations). That RK4 solution agrees with the result above within 0.053 nm in
position and 4e-12 in p/(mc), meeting its configured 1e-10 m / 1e-11 mc tolerances.
The input includes no external field or distribution files. Run it in a fresh
directory because the solver rejects an existing `orbit.json`.

## Preserve the actual legacy field placement

Two details are essential when converting this notebook:

1. Old `MAGNET_END` advances the start of the next element. In each cell the
   original F element starts at lattice phase -6°, with `MAGNET_START=14.23°`
   and `MAGNET_END=18.73°`. The D start is therefore measured after 18.73°:
   `-6 + 18.73 + 1.52 = 14.25°`, not `-6 + 1.52°`.
2. The recovered old Enge implementation returns `getCentreLength()=x0/2`.
   With x0=L/2 and the old FFA centre offset of half that value, its field centre
   was `MAGNET_START + L/8`. Current OPALX uses the corrected centre length L
   and an entrance-to-centre offset L/2. To preserve the old field, translate
   each real magnet's new entrance by **-3L/8**. Do not change the modern Enge
   implementation or its profile width.

Consequently the first-cell new F entrance/centre are 6.5425° / 8.7925°, and
D entrance/centre are 13.40625° / 14.53125°. Repeat these by 22.5° for 16 cells.
These offsets are supported by local historical source at `265ec2879^`:
`src/Elements/OpalScalingFFAMagnet.cpp`,
`src/Classic/AbsBeamline/ScalingFFAMagnet.cpp`,
`src/Classic/AbsBeamline/EndFieldModel/{Enge,AsymmetricEnge}.h`, and
`src/Classic/AbsBeamline/Ring.cpp`. The same Enge length convention is still
visible immediately before `c27447726` (September 2026).

The generator retains the notebook's exact spiral argument **3.142/4**, field
index 7.717022064259446, fourth-order expansion, 0.12 m fringe parameter,
coefficients {-1.3,3.42,-1.2,0.23}, radial bounds R0-0.8 / R0+1.2 m and full
height 1 m. Converting the old horizontal x-y / vertical z system into OPALX's
horizontal x-z / vertical y system reverses magnetic-field signs: F B0 becomes
+0.6532097983351404 T; D B0 becomes -0.502765236428169 T.

There are 32 real magnets plus 33 explicitly named `B0=0` gap arcs. The gap arcs
carry only the circular nominal design geometry and the same broad aperture;
they add no field or bend force. This retains strict nominal end-to-start
closure and a 2πR0 design circumference without relaxing `GEOMTOL/ANGLETOL` or
truncating neighbouring magnetic tails. Actual motion through zero-field space
still follows the Lorentz equation.

The original executable/version used for the stored notebook outputs is absent.
This is a legacy-compatible reconstruction from available source history, not a
claim of a controlled binary-to-binary regression.

## The notebook's “12 MeV” orbit has higher particle energy

The old input file adds radial momentum to the full reference tangential
momentum, rather than reducing the tangential component to keep total energy
fixed. For its saved pr=0.016268421156887713 mc and reference p=0.16044471664008755
mc, the actual launch has approximately **12.1225864268454 MeV** kinetic energy.
COF instead holds the specified total energy fixed.

The notebook stores r=4.395245598494839 m. At its actual higher energy, the
native COF and independent polar equations for the reconstructed field give
r=4.395241339211451 m, pr=0.016269912780107 mc: a 4.26 µm radius difference.
Legacy stepping, output precision and the external optimiser are not a controlled
accuracy baseline. The fixed-12 MeV result above should not be forced to match
the notebook's higher-energy launch. `school-notebook-energy.json` records the
native high-order COF result for that energy.

## Validation and reproduction

Use `~/.venv-h6` for the Python scripts. Example commands from the OPALX root:

```sh
cmake --build build-omp --target opalx_exe TestScalingFFAMagnet TestClosedOrbitSolver TestClosedOrbitInitialState TestOneTurnMap -j 6
~/.venv-h6/bin/python docs/validation/scaling-ffa/school_cof.py --output /tmp/school-cof --dt 1e-10 --integrator DOP853
~/.venv-h6/bin/python docs/validation/scaling-ffa/school_cof.py --output /tmp/school-cof-mpi --ranks 2 --dt 1e-10 --track-turns 3
~/.venv-h6/bin/python docs/validation/scaling-ffa/school_reference.py --output /tmp/school-reference
~/.venv-h6/bin/python docs/validation/scaling-ffa/check_school_cof.py --cases /tmp/school-cof /tmp/school-cof-mpi --output /tmp/school-checks.json
```

A new directory is required for each run. `--initial x px y py` selects the
section seed; the default uses the stored notebook radius/momentum. `--energy-mev`
selects fixed total kinetic energy; it defaults to 12. To reproduce the
higher-energy comparison, pass `--energy-mev 12.122586426845434` to the independent
reference script. The scripts supply the local MPI slot-mapping workaround and
one OpenMP thread per rank.

`school_reference.py` independently integrates median-plane polar equations
using an explicit polynomial/logistic Enge sum. It uses neither OPALX fields nor
its tracker. With total normalized momentum p0 and radial momentum pr, it solves

```
dr/dphi  = r*pr/sqrt(p0^2-pr^2)
dpr/dphi = sqrt(p0^2-pr^2) - c*r*B/(mass_eV)
```

for a one-cell fixed point and verifies a full turn. The independent full-turn
residual is below 4e-14 in the respective coordinate units. The comparison
thresholds are 1 nm and 1e-10 mc, covering the configured Newton tolerances and
integration differences without weakening existing tests.

Recorded checks (`school-results.json`):

- A native four-dimensional solve converges from the displaced notebook seed.
- RK4 at 100, 50 and 25 ps, plus DOP853 at 100 ps, agree with the independent
  orbit to better than 0.06 nm / 4e-12 mc. Reducing RK4 from 50 to 25 ps lowers
  relative energy drift from 1.59e-10 to 1.01e-11 (15.8×), consistent with fourth
  order. The 25 ps case also halves FDSTEP to check derivative sensitivity.
- DOP853 one- and two-rank COF results agree, including the returned matrix.
- Eight particles launched through `TRACK, INITIALORBIT=ORBIT` complete three
  directed turns. Comparing the first two complete returns, the maximum radial
  return error is 3.09 µm at 100 ps on two ranks and 0.651 µm at 50 ps on one
  rank. The checks require <5 µm position and <5e-7 mc radial-momentum error;
  these allow the measured second-order Boris error and endpoint interpolation.
  They are transport checks, not tighter COF residual tolerances.
- An independent weak-focusing radial FFA ring (k=-0.5) converges from a 4D
  displacement to its known circular orbit within 3e-14 m. The runnable fixture
  is in `cof-analytic-ring/`.
- Initial radial and vertical aperture violations reject with an explicit loss
  diagnostic and produce no successful JSON result.
- All 67 tests in six relevant executables pass on one and two ranks
  (`school-unit-results.json`).

Raw successful runs are under `output/ffa-school-cof-legacy-*`; the 25 ps result
is in `output/ffa-school-cof-legacy-fine-v2`. An initial 25 ps attempt exhausted
its 20000-step budget before a complete turn. Increasing the generated safety
budget to 100000 fixed that test; no tolerance or physics change was made.
Earlier uncompensated-centre diagnostic runs are not school results.

The orbit figure is reproducible with:

```sh
~/.venv-h6/bin/python docs/validation/scaling-ffa/school_plot.py --track-dir /tmp/school-cof-mpi --output /tmp/school-figures
```

The figure overlays the independent polar orbit and the actual production Boris
reference traversal, in full-ring and one-cell views. The subsequent
[production tune analysis](TUNES.md) follows small oscillations about this orbit
and checks them against independent cell optics. Native spectral tune
support, a public probe interface, the complete notebook conversion and
space-charge validation remain separate work.
