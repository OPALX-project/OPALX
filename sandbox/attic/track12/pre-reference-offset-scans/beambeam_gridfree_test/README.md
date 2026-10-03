# BeamBeam Default-Mesh Quick Test

Purpose: test pair-4 one-step witness fields with the BeamBeam active path but
without replacing the normal bunch-following field domain by the frozen
BeamBeam interaction-window mesh.

Implementation:

- Added diagnostic runtime switch:
  `OPALX_BB_DISABLE_FROZEN_WINDOW_MESH=1`.
- With the switch set, the BeamBeam state remains active and passive witness
  gather still runs, but `computeSelfFields()` uses the usual bunch-following
  field domain.

Commands:

```sh
# fixed BeamBeam window baseline
env OMP_NUM_THREADS=8 \
  OPALX_BB_WITNESS_KICK_CSV=fixed_kicks.csv \
  OPALX_BB_WITNESS_KICK_STEPS=1 \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx signcheck.in

# BeamBeam active, but default bunch-following mesh
env OMP_NUM_THREADS=8 \
  OPALX_BB_DISABLE_FROZEN_WINDOW_MESH=1 \
  OPALX_BB_WITNESS_KICK_CSV=free_kicks.csv \
  OPALX_BB_WITNESS_KICK_STEPS=1 \
  /Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx signcheck.in
```

Pair-4 c1/electron first kick:

| case | `Ex` [V/m] | `By` [T] | `|E|` [V/m] | `|B|` [T] | c1/e- `dPx` |
|---|---:|---:|---:|---:|---:|
| fixed BeamBeam window | `+6.268184e9` | `+20.9084` | `6.268187e9` | `20.9084` | `-4.926876e-4` |
| default-following mesh | `-3.091357e13` | `-103116` | `3.091673e13` | `103127` | `+7.55208` |
| manufactured copropagating | `-4.846873e9` | `-16.1674` | `4.846873e9` | `16.1674` | expected positive |

Findings:

- Disabling the frozen BeamBeam window mesh fixes the pair-4 field sign.
- The magnitude remains far too large: `|E|` is about `6379x` the manufactured
  smooth-Gaussian value.
- Therefore the fixed window is responsible for the original wrong sign, but
  the larger problem is not only the fixed aperture.  The default mesh still
  resolves/evaluates the random 400k macroparticle source in a way that is much
  stronger than the manufactured smooth Gaussian at the witness point.

Artifacts:

- `fixed_window_pair4_kicks.csv`
- `default_following_mesh_pair4_kicks.csv`
