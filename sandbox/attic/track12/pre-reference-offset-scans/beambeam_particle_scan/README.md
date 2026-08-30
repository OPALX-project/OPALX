# BeamBeam Particle-Count Scan: Pair-4 First Kick

Purpose: test whether the BeamBeam field blow-up seen after transverse-cell
refinement is mainly source macroparticle shot noise.

Command:

```sh
for np in 50000 100000 200000 400000 800000; do
  ~/.venv-h6/bin/python sandbox/track12particles/opalx/scan_beambeam_window_convergence.py \
    --force \
    --output-dir sandbox/track12particles/opalx/beambeam_particle_scan \
    --apertures-m 0.0002 \
    --nxy-values 32,64 \
    --particles "$np"
done

~/.venv-h6/bin/python sandbox/track12particles/opalx/scan_beambeam_window_convergence.py \
  --output-dir sandbox/track12particles/opalx/beambeam_particle_scan \
  --apertures-m 0.0002 \
  --nxy-values 32,64 \
  --particle-values 50000,100000,200000,400000,800000
```

Reference manufactured field, copropagating branch:

| quantity | value |
|---|---:|
| `Ex` | `-4.84687271e9 V/m` |
| `By` | `-16.16739265 T` |
| `|E|` | `4.84687288e9 V/m` |
| `|B|` | `16.16739265 T` |

Results:

| particles | cell [um] | `Ex` [V/m] | `By` [T] | `|E|` ratio | c1/e- `dPx` |
|---:|---:|---:|---:|---:|---:|
| 50k | 6.25 | `-5.32965e11` | `-1777.78` | 110.058 | `+0.0424584` |
| 100k | 6.25 | `-6.54512e11` | `-2183.21` | 135.070 | `+0.0524884` |
| 200k | 6.25 | `-5.96507e11` | `-1989.73` | 123.082 | `+0.0476774` |
| 400k | 6.25 | `-6.32980e11` | `-2111.39` | 130.757 | `+0.0506991` |
| 800k | 6.25 | `-6.33591e11` | `-2113.43` | 130.809 | `+0.0507488` |
| 50k | 12.5 | `+7.07653e11` | `+2360.47` | 146.038 | `-0.0569373` |
| 100k | 12.5 | `+7.15716e11` | `+2387.37` | 147.668 | `-0.0576153` |
| 200k | 12.5 | `+7.16090e11` | `+2388.61` | 147.743 | `-0.0576468` |
| 400k | 12.5 | `+7.10000e11` | `+2368.30` | 146.486 | `-0.0571341` |
| 800k | 12.5 | `+7.07808e11` | `+2360.99` | 146.034 | `-0.0569497` |

Findings:

- At fixed cell size, increasing the number of primary macroparticles from
  `50k` to `800k` does not drive the OPALX field toward the manufactured field.
- The `12.5 um` cell remains wrong-sign and about `146x` high.
- The `6.25 um` cell remains correct-sign but about `110x` to `135x` high.
- The field is therefore not dominated by ordinary finite-particle shot noise
  in this range.  The dominant effect appears deterministic in the BeamBeam
  window mesh/deposition/evaluation path.

Artifacts:

- `pair4_beambeam_window_scan.csv`
- `pair4_ex_vs_particles.png`
- `pair4_by_vs_particles.png`
- `pair4_e_ratio_vs_particles.png`

Next checks:

- Replace random macroparticle deposition by a deterministic smooth Gaussian
  source on the BeamBeam window mesh, or add a direct analytic grid deposition
  diagnostic, to isolate CIC/grid Green-function behavior.
- Add a minimal regression test around the two-electron BeamBeam sign probe so
  future changes cannot hide this behavior.
