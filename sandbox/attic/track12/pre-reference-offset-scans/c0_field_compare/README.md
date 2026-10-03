# Container 0 Primary Field Comparison

Purpose: compare OPALX fields on the primary electron bunch (`c0`) against the
smooth manufactured anisotropic Gaussian field at the same particle coordinates.

Setup:

- Pair 4, one OPALX step.
- Primary source: `400000` macroparticles.
- Field grid in the input deck: `NX=32`, `NY=32`, `NZ=64`.
- Two OPALX cases:
  - fixed BeamBeam interaction-window mesh,
  - BeamBeam active but default bunch-following mesh using
    `OPALX_BB_DISABLE_FROZEN_WINDOW_MESH=1`.
- The full c0 raw source kick CSVs were written in `/tmp/opalx-c0-analytic/`
  and are about `380 MB` each.  This directory stores a reproducible script and
  a deterministic `20000`-particle sample.

Commands:

```sh
# fixed BeamBeam window
env OMP_NUM_THREADS=8 \
  OPALX_BB_SOURCE_KICK_CSV=source_kicks.csv \
  OPALX_BB_SOURCE_KICK_STEPS=1 \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx signcheck.in

# BeamBeam active, but default bunch-following mesh
env OMP_NUM_THREADS=8 \
  OPALX_BB_DISABLE_FROZEN_WINDOW_MESH=1 \
  OPALX_BB_SOURCE_KICK_CSV=source_kicks.csv \
  OPALX_BB_SOURCE_KICK_STEPS=1 \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx signcheck.in

# analytic comparison
~/.venv-h6/bin/python sandbox/track12particles/opalx/compare_c0_fields_to_analytic.py \
  --fixed-csv /tmp/opalx-c0-analytic/fixed/source_kicks.csv \
  --free-csv /tmp/opalx-c0-analytic/free/source_kicks.csv
```

Summary over the sampled `20000` primary particles:

| case | quantity | median | p05 | p95 | mean |
|---|---|---:|---:|---:|---:|
| fixed window | `|E|_OPALX / |E|_analytic` | `1.41716` | `1.05512` | `3.12508` | `1.65897` |
| fixed window | `E` direction cosine | `-0.999993` | `-1.00000` | `-0.999826` | `-0.999895` |
| fixed window | `|B|_OPALX / |B|_analytic` | `1.41716` | `1.05512` | `3.12508` | `1.65897` |
| fixed window | `B` direction cosine | `-0.999993` | `-1.00000` | `-0.999826` | `-0.999895` |
| default-following mesh | `|E|_OPALX / |E|_analytic` | `6218.11` | `4001.34` | `7318.62` | `6012.46` |
| default-following mesh | `E` direction cosine | `0.999844` | `0.995535` | `0.999999` | `0.997940` |
| default-following mesh | `|B|_OPALX / |B|_analytic` | `6218.11` | `4001.34` | `7318.62` | `6012.46` |
| default-following mesh | `B` direction cosine | `0.999844` | `0.995535` | `0.999999` | `0.997940` |

Findings:

- On `c0`, the fixed BeamBeam window produces fields that are almost exactly
  anti-parallel to the analytic smooth Gaussian field.  The magnitude is only
  moderately high, with median `|E|` ratio `1.42`.
- With the default bunch-following mesh, OPALX fields are almost exactly
  parallel to the analytic field, but the magnitude is enormous: median `|E|`
  ratio `6218`.
- This mirrors the passive witness behavior: disabling the fixed BeamBeam mesh
  repairs the sign but leaves a large field-strength problem.

Artifacts:

- `c0_field_summary.csv`
- `fixed_window_c0_field_sample.csv`
- `default_following_mesh_c0_field_sample.csv`
- `c0_field_sample_combined.csv`
- `c0_e_abs_ratio_hist.png`
- `c0_e_direction_cosine_hist.png`
- `fixed_window_c0_e_xy.png`
- `default_following_mesh_c0_e_xy.png`

The `*_c0_e_xy.png` plots show one primary bunch only: OPALX `|E|`,
manufactured-solution `|E|`, and their magnitude ratio at the sampled c0
particle coordinates in the transverse `x,y` plane.
