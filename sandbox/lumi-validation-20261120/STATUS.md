# Lumi validation — 2026-11-20

Complete: build, relevant single-/two-GPU tests, all four scientific smokes,
and the full Track-12 run passed their completion/integrity checks. All outputs
were fetched and analyzed locally. No ID allocator change was necessary.

| Stage | Slurm job | Result |
|---|---:|---|
| Release HIP build + focused n1 tests | 23901108 | completed, 10m51s |
| n2 tests + experiments 1–4 | 23901410 | completed, 58s |
| Track-12, one MI250 | 23901513 | completed, 469.86s launch elapsed |

Track-12 retained all 12 witnesses through global step 1500: 6504 actual H5
records per species, finite phase space, no duplicate or disappearing IDs.
The CAIN comparison contains 13012 samples including four reconstructed birth
states. These are integrity/workflow checks, not a new convergence acceptance.

| Reference | relative L2 x | relative L2 s | y RMS difference |
|---|---:|---:|---:|
| Three-sigma-truncated manufactured model | 0.853738% | 0.117302% | 1.73728 micrometres |
| CAIN | 9.482848% | 2.967314% | 14.76368 micrometres |

The manufactured y reference is zero, so its relative L2 y is undefined.
CAIN's positron x relative error is 107.906% but its absolute RMS is only
5.76585 micrometres; do not describe every trajectory as uniformly few-percent
agreement. Later first kicks are under-resolved at this mesh (pair-6 OPALX/MSOL
px ratio about 0.797). No numerical tolerances were changed.

Numerical reports: `track12_summary.json`, `track12_validation.csv`, and
`track12/results/track12_comparison.json`. Comparison plots and pointwise CSV
are under `track12/results/`.

Two-rank smoke gate required two launcher-only adjustments: respect Slurm's
ROCR visibility mask rather than applying a second rank-based HIP mask, and
disable Cray GPU IPC (`MPICH_GPU_IPC_ENABLED=0`) after an HSA IPC-import failure.
GPU-aware MPI stays enabled. Distinct devices confirmed by PCI bus addresses
d1:00.0 and d6:00.0. Failed attempts 23901330/23901382 are preserved in the logs;
the successful full smoke campaign is 23901410 (58 seconds).

Local `analyze_smokes.py` validated every experiment; see `smoke_status.csv`.
Actual HIP executable SHA256 is
`31e97e63d5208ef210c4188b717b48ccd739fba37874f38735e08754372d1307`.
Experiment 3 first-electron momentum magnitude is 1.322336 mc, with all exact
births, finite phase space and persistent IDs. No invalid first-kick spike.

Source commit: `e24dc92185198c5f65718e057400e6c6398dce56`, pushed to
`271-implement-interation-point-element`. No source changes after that commit.

Isolated source: `/flash/project_465002453/opalx-beambeam-recovery-20261120`.
Run root: `/flash/project_465002453/bb-validation-e24dc9218-20261120`.
The existing dirty master checkout `~/opalx` was not modified.

Build: Slurm `23901108`, dev-g, one GPU, 20 CPUs, Release HIP/GFX90A,
LUMI25.03/AMD6.3.4/Cray-MPICH8.1.32. All build and scheduler logs are run-local.
The build script ran three focused n1 tests after compilation before creating
`build-and-unit.completed`. The initial dependent smoke submission hit the
four-submitted-job limit; it was retried after the build released its slot.
No unrelated jobs or the existing dirty master checkout were changed.

Prepared smoke decks (locally generated and transferred, no remote Python setup):

| Case | primary particles | mesh | steps |
|---|---:|---|---:|
| experiment 1 | 1000 | 8x8x16 | four field-separation cases |
| experiment 2 | 2000 deterministic | 16x16x32 | 3 |
| experiment 3 | 10000 | 16x16x16 | 160; all 1297+1297 births |
| experiment 4 | 32 | 4x4x4 | 1607; primary-only lifecycle |

Smoke script additionally runs emitted-file, BeamBeam, field-service and
diagnostic tests on two GPU ranks. Inputs, scripts and output belong to their
run directory. Scientific analysis/plotting is local with `~/.venv-h6`.

Only after local smoke validation was `smokes-validated-locally` created and
`track12/submit.sh` submitted from `track12`, with scheduler logs in that same
directory. Track12: one MI250 GPU, 400000 deterministic primary macroparticles,
256x256x128, 1501 steps, 12 witnesses, dt = 6.0041537135667364 fs, `--info 4`.
Transferred primary SHA256:
`625621deae71b5e66522dfc8b9c17d02a75bf9db3e7cd6a0d70f5e94aa59ccd5`.

ID concern resolved without allocator changes: birth allocation is MPI-rank
strided. Every prior witness ID persists in all 38 dumps at both 1 and 2 ranks;
zero duplicates. Enhanced per-row birth-ID/migration regression passes on
1/2/4 ranks. This does not imply same file row has the same ID across MPI sizes.

Local broad GCC tests: 147/148 pass. Unchanged master bend-field exact-zero
assertion differs by 1.39e-17; no code/tolerance change. Relevant BeamBeam,
emission and space-charge tests pass. Raw historical CAIN inputs retained with
original CRLF/whitespace. See root HANDOFF.md for live continuation details.

To reproduce local analysis after fetching the output:

```sh
~/.venv-h6/bin/python sandbox/lumi-validation-20261120/analyze_smokes.py
~/.venv-h6/bin/python sandbox/lumi-validation-20261120/analyze_track12.py
```

Lumi campaign scripts and reports are retained locally; generated data/binaries
are ignored by this directory's `.gitignore` and were not committed.
