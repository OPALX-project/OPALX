# BeamBeam Window Scan: Pair-4 First Kick

Purpose: isolate how the BeamBeam frozen transverse window affects the pair-4
one-step witness field.  The scan keeps the source/witness timing fixed and
varies only the BeamBeam aperture and transverse grid size.

Command:

```sh
~/.venv-h6/bin/python sandbox/track12particles/opalx/scan_beambeam_window_convergence.py --force
```

Reference manufactured field, copropagating branch:

| quantity | value |
|---|---:|
| `Ex` | `-4.84687271e9 V/m` |
| `By` | `-16.16739265 T` |
| `|E|` | `4.84687288e9 V/m` |
| `|B|` | `16.16739265 T` |

OPALX BeamBeam first-kick results for c1/electron:

| aperture [um] | NXY | cell [um] | `Ex` [V/m] | `By` [T] | `|E|` ratio | c1/e- `dPx` |
|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 32 | 62.500 | `+6.268184e9` | `+20.9084` | 1.293 | `-4.92688e-4` |
| 1000 | 64 | 31.250 | `+5.02526e10` | `+167.624` | 10.368 | `-3.95039e-3` |
| 200 | 32 | 12.500 | `+7.10000e11` | `+2368.30` | 146.486 | `-5.71341e-2` |
| 200 | 64 | 6.250 | `-6.32980e11` | `-2111.39` | 130.757 | `+5.06991e-2` |
| 100 | 32 | 6.250 | `-6.35952e11` | `-2121.30` | 131.369 | `+5.09459e-2` |
| 100 | 64 | 3.125 | `-1.81548e13` | `-60557.8` | 3745.69 | `+3.94806` |
| 50 | 32 | 3.125 | `-1.81548e13` | `-60557.9` | 3745.69 | `+3.94807` |
| 50 | 64 | 1.562 | `-3.86323e13` | `-128863` | 7971.19 | `+9.77445` |
| 20 | 32 | 1.250 | `-4.23750e13` | `-141348` | 8743.91 | `+10.8574` |
| 20 | 64 | 0.625 | `-5.13254e13` | `-171203` | 10591.7 | `+13.4558` |

Findings:

- The original 1 mm aperture / NXY=32 case has the wrong sign relative to the
  manufactured electron-source field.
- Refining the transverse cell size eventually flips the sign to the
  manufactured sign, around cell size `6.25 um`, but the magnitude is then
  already about `130x` too large.
- Cases with the same transverse cell size are nearly identical even when
  aperture and NXY differ, e.g. `200 um / 64` and `100 um / 32`, or
  `100 um / 64` and `50 um / 32`.  The pathology is controlled primarily by
  transverse cell size.
- This is not a viable aperture tuning knob: smaller cells fix the sign but
  amplify the sampled field far above the smooth manufactured Gaussian model.

Artifacts:

- `pair4_beambeam_window_scan.csv`
- `pair4_ex_vs_cell_size.png`
- `pair4_by_vs_cell_size.png`
- `pair4_e_ratio_vs_cell_size.png`

Next checks:

- Repeat with larger primary macroparticle counts or a deterministic smooth
  source deposition to separate grid resolution from macroparticle shot noise.
- Scan particle count at fixed `(aperture, NXY)` near the sign transition,
  especially cell sizes `12.5 um` and `6.25 um`.
- Add a minimal BeamBeam regression test with a two-electron source/witness
  geometry so this sign behavior is pinned before changing physics code.
