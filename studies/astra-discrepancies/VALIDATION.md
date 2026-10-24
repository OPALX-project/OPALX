# Current-branch validation

2026-10-24. This check applies to the issue-branch integration on base
`8c1a862`, including the final CM correction, robust `hypot` gamma evaluation
and new stationary/counterstreaming regressions. It is separate from the
historical frozen-executable validation described in [PROVENANCE.md](PROVENANCE.md).

**All eight selected CTest suites passed, 8/8, in 5.49 s.**
`TestSingleBinCMParallel` uses two MPI ranks, including the empty-rank case;
the other seven selected suites use one rank. OpenMP uses four threads.

A before-fix control reproduced the full-pipeline failure for the new
stationary ±15 normalized-momentum case. The robust gamma evaluation passes
that regression and the nearly stationary case with a 10⁻¹² perturbation.
These tests address floating-point robustness at the stationary limit;
they do not constitute new EuXFEL tracking or regenerate the archived tables.

From the repository root, with the existing configured OpenMP build:

```sh
cmake --build omp-build --target \
  TestSpaceChargeSolver TestParticleBinTraversal TestSingleBinCMParallel \
  TestCartesianPIC3DComponents TestCartesianPIC3DAlgorithms \
  TestCartesianPIC3DFields TestPoissonSolver TestP3MRestFrameParallel -j4

PRTE_MCA_rmaps_default_mapping_policy=slot:oversubscribe \
OMP_NUM_THREADS=4 OMP_PROC_BIND=false \
ctest --test-dir omp-build --output-on-failure \
  -R 'TestSpaceChargeSolver|TestParticleBinTraversal|TestSingleBinCMParallel|TestCartesianPIC3DComponents|TestCartesianPIC3DAlgorithms|TestCartesianPIC3DFields|TestPoissonSolver|TestP3MRestFrameParallel'
```

No tracking was performed for this integration/commit. The experimental spatial
operators remain evidence from earlier frozen study executables; their results
must not be read as new trajectories from the current branch.

Associated manual pages reside in the separate `opalx-manual` repository.
Those manual changes are not included in this repository commit or evidence
bundle. The current branch's implementation documentation and tests should be
reviewed together with that separate manual work.

The portable archive itself was checked for byte-for-byte equality of all
copied figures/CSVs, exact numeric identity of seven CSV copies, resolving local
Markdown links, absence of private machine paths and complete SHA-256 coverage.
Run `shasum -a 256 -c SHA256SUMS` from this directory to verify the bundle.
