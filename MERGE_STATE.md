# BeamBeam master integration state

Updated: 2026-08-07

## Goal

Merge `origin/master` into `271-implement-interation-point-element` with the
master tree and architecture authoritative, then port the minimum BeamBeam
element, field-solver, tracker, diagnostics, and test integration needed to
retain the feature.

## Safety and merge state

- Pre-merge tracked work is in Git stash
  `pre-origin-master integration 2026-08-07 tracked work`.
- Binary patch backup:
  `/tmp/opalx-beambeam-pre-master-2026-08-07.patch`.
- Analytic/integration workspace backup:
  `/tmp/opalx-beambeam-pre-master-untracked-2026-08-07.tar.gz`.
- Merge target: `origin/master` at
  `0dea93e9a61bcee2081d3d9419239872d45b17b3`.
- Current IPPL dependency: `36a4ca62a52e36a8a2c945c9048d0f010971d309`.
- All master-owned paths were reset to the exact master version and master
  deletions were accepted. Obsolete `Component` and interaction-animation
  abstractions were not restored.
- The merge remains active until the reviewed integration is staged and
  committed. No push is authorized.

## Integrated design

- `BeamBeam` now uses the current `ElementBase` and `Geometry` interfaces.
- The tracker retains master's per-container orbit threaders, emission-time
  reselection, source-loss handling, P3M path, and load balancing. BeamBeam is
  added through focused lifecycle, coordinate-transform, fixed-mesh, field,
  witness, retirement, and diagnostics hooks.
- Container 0 is the self-consistent source. Configured witness containers
  gather the source field without depositing charge.
- The fixed BeamBeam mesh follows master's cell-centred convention:
  `h = span / (N - 1)` with the origin shifted down by half a cell.
- Source scatter/copy behavior is retained for the legacy non-binned solver.
  Field arrays are mirrored and restored in a CUDA-safe way.
- BeamBeam rho, phi, and E diagnostics support tagged ASCII/HDF5 output.
  `C0PSDUMPFREQ` independently filters container-0 particle dumps.
- Mesh visualization retains the element type and colors without restoring the
  removed legacy movie-script implementation.

## Validation completed

- Configured `build_openmp` with current IPPL and Kokkos 5.2.0.
- Built `opalx_exe`, `TestBeamBeam`, and
  `TestBeamBeamDiagnosticsWriter` successfully.
- Built the complete default target, including every configured unit-test
  executable, successfully. The full build exposed that the initial
  `BeamlineVisitor::visitBeamBeam` hook was pure virtual and therefore made
  existing specialized test visitors abstract. The hook now has a default
  no-op; BeamBeam-aware production visitors continue to override it.
- `TestMultipoleT`: 5/5 tests passed after the visitor compatibility fix.
- `TestVariableRFCavity`: 15/15 tests passed after the visitor compatibility
  fix.
- `TestBeamBeamDiagnosticsWriter`: 4/4 tests passed, one MPI rank.
- `TestBeamBeam`: 11/11 passed with one rank and 11/11 with two ranks.
- Manufactured static BeamBeam input completed normally with one rank and
  produced rho/phi/E diagnostics.
- Static copied-source BeamBeam smoke test completed normally with two ranks,
  including inactive, overlap, non-overlap, and completed lifecycle states.
- One-source plus electron/positron witness smoke test completed 700 steps with
  one rank. Both witness containers emitted six particles and remained active
  through BeamBeam completion.
- The maintained `fields` sandbox analysis matches its accepted baseline.

## Validation limitations and remaining risks

- CTest's local MPI wrapper cannot map even one rank on this machine
  (`all allocated nodes are already filled`). Direct `mpiexec` runs with an
  explicit localhost mapping were used instead.
- The full sandbox regression driver cannot run from the current checkout
  because `sandbox/OPALX-IMPACT/extracted_graph_data.csv` and
  `sandbox/track-e-p/gamma_gamma_pairs-2_c0.stat` are absent.
- Several accepted sandbox decks use pre-master syntax: negative
  `PSDUMPFREQ`, partial FFT decomposition, `PC` on `FROMFILE` beams, or the old
  `COPY` attribute. Temporary copies were adapted for smoke testing; repository
  decks and accepted baselines were not overwritten.
- Master's placement and cell-centred mesh conventions change the manufactured
  field snapshot indices and numerical baseline. Refreshing that accepted
  physics baseline requires an explicit acceptance decision.
- BeamBeam copied-source deposition is currently limited to the legacy
  non-binned field-solver path. Binned BeamBeam copy behavior remains a follow-up.

## Next step

Keep the sandbox baseline migration separate from the code integration. Do not
push without explicit permission.
