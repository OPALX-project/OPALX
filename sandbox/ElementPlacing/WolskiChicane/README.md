# Wolski TESLA chicane benchmark

This sandbox case packages a simple four-bend electron chicane based on the
TESLA bunch-compressor numbers in Andy Wolski's bunch-compressor lecture notes.

Reference values:

- beam kinetic energy: `5.0 GeV`
- injected RMS bunch length: `sigma_z = 6.0 mm`
- injected uncorrelated RMS energy spread: `sigma_delta = 1.3e-3`
- target extracted RMS bunch length: `sigma_z,f = 0.3 mm`
- RF frequency: `1.3 GHz`

Using Wolski's linear compression model,

```math
\sigma_{z,f} \approx |R_{56}| \sigma_{\delta,u},
```

the target compressor strength is

```math
R_{56} = -0.230480589\ \mathrm{m},
\qquad
h = -\frac{1}{R_{56}} = 4.327913277\ \mathrm{m}^{-1}.
```

The equivalent ideal RF voltage for a `1.3 GHz` chirper at `5 GeV` is

```math
V = 794.229541\ \mathrm{MV}.
```

## Files

- `wolski_tesla_chicane.bmad`
  - Bmad lattice with an ideal `lcavity` chirper and a symmetric four-bend
    chicane. Includes zero-length `chicane_in` and `chicane_out` markers.
- `wolski_tesla_chicane_e0.bmad`
  - zero-pole-face-angle RF variant of the Bmad lattice with the same entrance
    and exit markers.
- `wolski_tesla_chicane_no_rf_e0.bmad`
  - zero-pole-face-angle no-RF Bmad transport lattice with the same entrance
    and exit markers.
- `wolski_tesla_chicane_beam_track.in`
  - `bmad_beam_rms_track` input for the Bmad lattice using an explicit
    unchirped particle file. Also names the per-step RMS output file.
- `wolski_tesla_chicane_e0_beam_track.in`
  - `bmad_beam_rms_track` input for the zero-face RF Bmad lattice using the same
    explicit unchirped particle file.
- `wolski_tesla_chicane_bmad_unchirped.beam`
  - BMAD canonical particle input used for the RF benchmarks. The cavity then
    generates the chirp internally.
- `wolski_tesla_chicane_no_rf_e0_negchirp_beam_track.in`
  - `bmad_beam_rms_track` input for the zero-face no-RF Bmad lattice with the
    embedded negative chirp.
- `bmad_beam_rms_track.f90`
  - custom Bmad driver that writes projected RMS phase-space values
    `rms_x, rms_px, rms_y, rms_py, rms_z, rms_pz` from the bunch covariance
    matrix at every saved beamline point.
- `build_bmad_beam_rms_track.sh`
  - local build helper for the custom Bmad RMS driver.
- `opalx/wolski_tesla_chicane_opalx_negchirp_dt_1em11_monitor_both.in`
  - OPALX no-RF deck with temporal monitors at the chicane entrance and exit.
- `opalx/c-chicane-configuration.pdf`
  - reference sketch for the C-chicane geometry and bend sign convention used
    in the OPALX-side studies.
- `opalx/`
  - active OPALX run directory. The `.in` decks, `.stat` files, monitor files,
    `.h5` dumps, and log-like runtime artifacts live here.
- `opalx/data/`
  - OPALX geometry/export artifacts written by the element-position tooling, for
    example `*_DesignPath.dat`, `*_ElementPositions.txt`, `*_ElementPositions.sdds`,
    `*_ElementPositions.py`, and `*_3D.opal`.
- `opalx/old/`
  - archived OPALX input decks from earlier chicane studies, including
    `wolski_tesla_chicane_opalx.in` and
    `wolski_tesla_chicane_opalx_signed_e0.in`.
- `generate_wolski_chicane_distribution.py`
  - writes `wolski_tesla_chicane_distribution.txt` for the OPALX run.
- `summarize_wolski_chicane.py`
  - compares the Bmad and OPALX output bunch lengths with Wolski's target.

## Model split

The Bmad case uses an ideal `lcavity` to generate the chirp internally while
tracking an explicit particle file, `wolski_tesla_chicane_bmad_unchirped.beam`.
The lattice now also contains explicit `chicane_in` and `chicane_out` markers
at the same planes sampled by the OPALX entrance and exit monitors.

The OPALX case uses the same four-bend chicane, but the longitudinal chirp is
embedded directly in the input distribution:

```math
\delta(z) = h z + \delta_u,
\qquad
\delta_u \sim \mathcal{N}(0, \sigma_{\delta,u}^2).
```

This keeps the chicane comparison reproducible without introducing a separate
RF-field map into the sandbox case.

With the custom `bmad_beam_rms_track` program, BMAD now writes both the final
particle file and a plain-text RMS table for the bunch. The RMS table is built
from Bmad's bunch covariance matrix \(\Sigma\) using
\(\mathrm{rms}(u_i) = \sqrt{\Sigma_{ii}}\) for
\(u = (x, p_x, y, p_y, z, p_z)\).

Important limitation: in the standard Bmad bunch-tracking path, these saved RMS
points are at the start/end of each tracked element and marker. Even when an
individual bend uses internal Runge-Kutta substeps, the bunch-level save path
does not emit a dense intra-element history.

The physical BMAD bend-body sign convention for the chicane is

```text
(b1, b2, b3, b4) = (+, -, -, +),
```

with the edge angles following the same entrance/exit face orientation.

## Run

Generate the OPALX input bunch:

```bash
python3 generate_wolski_chicane_distribution.py
```

Build the Bmad RMS driver:

```bash
/Users/adelmann/git/opalx/sandbox/ElementPlacing/WolskiChicane/build_bmad_beam_rms_track.sh
```

Run Bmad:

```bash
/Users/adelmann/git/opalx/sandbox/ElementPlacing/WolskiChicane/bmad_beam_rms_track \
  wolski_tesla_chicane_beam_track.in
```

Run OPALX:

```bash
cd opalx
/Users/adelmann/git/opalx/build/src/opalx \
  wolski_tesla_chicane_opalx_negchirp_dt_1em11_monitor_both.in --info 2
```

The OPALX run products are then written relative to `opalx/`:

- main runtime outputs such as `*.stat`, `*.h5`, `*_Monitors.stat`, and
  monitor `*.loss` files stay in `opalx/`
- geometry/export side products such as `*_DesignPath.dat`,
  `*_ElementPositions.txt`, `*_ElementPositions.sdds`, `*_ElementPositions.py`,
  and `*_3D.opal` are organized under `opalx/data/`

Summarize the results:

```bash
/Users/adelmann/git/opalx/.venv-h6/bin/python summarize_wolski_chicane.py
```

For the no-RF negative-chirp case, the BMAD per-point RMS file is:

```text
wolski_no_rf_neg_rms.dat
```

## Current status

The sandbox now contains three relevant reference levels:

- Wolski linear target: `0.300000 mm`
- Bmad with ideal RF chirper: `0.344638 mm`
- Bmad with the same chirp embedded directly in the input beam (`wolski_no_rf_neg.dat`):
  `0.418681 mm`

This means the chirped-input benchmark is a valid approximation to the RF case,
and OPALX should be compared first against the Bmad no-RF benchmark before
blaming the missing RF element.

Current OPALX status:

- the baseline `wolski_tesla_chicane_opalx*.stat` variants remain far from the
  Bmad no-RF benchmark
- the signed explicit-bend variants improve when `DT` is reduced, but they do
  not yet converge to the Bmad result
- the main static blocker is the signed `SBEND` exit semantics in the return
  bends: for the parser-created `B2` case, the exact downstream drift entry
  point still maps to field-local `s \approx 0.1974 m` instead of `0.2 m`

So the remaining gap is not just a chirp-sign issue. The current evidence
points to an explicit signed-bend placement / field-chart consistency problem
in OPALX.

## Zero pole-face-angle comparison

To separate pole-face-angle effects from the signed-bend transport problem, the
benchmark now also includes zero-face-angle variants on both the Bmad and OPALX
sides.

Key observation:

- in Bmad, setting all pole-face angles to zero changes nothing in this setup
  because the lattice already uses `fringe_type = none`
- in OPALX, the signed zero-face case still exhibits the same short overlap
  windows around `B2 -> D1B` and `B3 -> NDB`

The resulting RMS beam sizes are:

```text
Case                                  sigma_x    sigma_y  sigma_z/s
Wolski target                                -         -   0.300000
Bmad RF                             18.812797   0.148883   0.344638
Bmad RF zero-face                   18.812797   0.148883   0.344638
Bmad no-RF (-h)                     19.130069   0.000000   0.418681
Bmad no-RF (-h) zero-face           19.130069   0.000000   0.418681
OPALX no-RF (-h), dt=1e-11          19.363506   0.000000   2.010764
OPALX signed zero-face              61.701890   0.000000   6.259259
```

So the pole-face angle itself is not the dominant effect here. The much larger
`sigma_x` and `sigma_s` in the signed zero-face OPALX run are consistent with
the same signed-bend handoff problem that already showed up in the overlap
diagnostics.
