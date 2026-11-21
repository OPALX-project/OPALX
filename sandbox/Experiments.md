# Scientific experiments

These are manual physics studies, **not** the short gtest/regression-tests-x
merge-gate tests. The frozen/rigid-primary approximation is unchanged. Local
coarse runs check the workflow; they do not replace accepted convergence studies
or establish agreement with CAIN's different space-charge model.

| Directory | Experiment | Local default and analysis |
|---|---|---|
| `experiment-1` | Independent anisotropic rigid two-Gaussian fields versus OPALX PIC fields | Four centroid separations; 20,000 primary particles, 32×32×64 mesh. Analytic field plots, OPALX field comparisons and error tables. |
| `experiment-2` | Six electron + six positron exact-timed CAIN witnesses | 2,000 deterministic primary particles, 16×16×32 mesh, all 1,501 steps. Reconstruct trajectories and compare all CAIN samples and the rigid-source reference. This intentionally coarse local profile is not an acceptance baseline. |
| `experiment-3` | Exact emission of all 1,297 CAIN electrons and 1,297 positrons | Convert and validate the raw table; run the 10,000-primary, 16³, 160-step input; analyze finite phase space, species and exact cumulative birth counts. |
| `experiment-4` | 32 cm production geometry and dynamic transverse field window | Default: 32-primary, 4³, 1,607-step lifecycle smoke, plus analysis of the full witness input. **No witnesses are tracked in this default smoke.** `--full` selects the full-population staged deck; it can be much slower. |
| `OPALX-IMPACT` | Independent historical drift/space-charge solver comparison | Restored byte-for-byte, including its original inputs, plotting script and outputs. Left unchanged: it is the explicit exception to the numbered layout and new launcher interface. |

## Running locally

From the repository root (or invoke a launcher by absolute path):

```bash
sandbox/experiment-1/run.sh
sandbox/experiment-2/run.sh
sandbox/experiment-3/run.sh
sandbox/experiment-4/run.sh
```

All launchers select `build_openmp/src/opalx` and `$HOME/.venv-h6/bin/python`.
Override these with `OPALX_EXE` and `OPALX_PYTHON`, or pass `--opalx PATH`.
They use at most four CPU threads; `--threads 1` through `--threads 4` is allowed.
No build is performed and no batch job is submitted implicitly.

Useful options:

```bash
sandbox/experiment-1/run.sh --smoke
sandbox/experiment-2/run.sh --smoke
sandbox/experiment-2/run.sh --particles 100000 --mesh 64 64 128
sandbox/experiment-3/run.sh --prepare-only
sandbox/experiment-4/run.sh --full
sandbox/experiment-2/run.sh --analyze-only /absolute/path/to/experiment-2/results/local-RUN
```

`--smoke` uses a smaller field solve in experiment 1 and only three integration
steps in experiment 2. A Track-12 smoke checks finite phase space and generates
diagnostic plots; it **does not** perform a complete CAIN trajectory comparison.
Experiments 3 and 4 already have small/default workflow profiles; `--smoke` does
not silently truncate the full-population birth schedule of experiment 3.

`--prepare-only` generates inputs (and input/analytic diagnostics) without running
OPALX. Each normal invocation creates a unique `results/local-…` directory.
`--output-dir` may name a new directory below that experiment's `results/`; an
existing directory is never overwritten. `--analyze-only` reanalyzes an existing
launcher-generated run and does not rerun OPALX.

`--merlin6` is reserved and currently exits with a clear "not implemented"
error. A100 submission will be added later. Older Slurm helpers are retained
under `scripts/merlin/` for reference, but are not the new remote interface.

## Layout and provenance

Each numbered experiment has `run.sh`, `inputs/`, `scripts/` and `results/`.
Shared launchers (`run_local.sh`, `run_experiment.py`) and the optional legacy
spherical model (`boosted_gaussian_witness.py`) are under `sandbox/scripts/`.
The spherical helper is needed only for experiment 2's `spherical` model option;
the current default uses the anisotropic model. All new logs,
generated input copies, H5/STAT output, plots, CSV/JSON reports and `REPORT.md`
belong under the run's `results/` tree. A manifest records the executable hash,
source/input hashes, Git state, commands, profile and completion/failure status.
No new HDF5 dependency was added to the regression suite; scientific H5 analysis
uses the existing study tools.

Previously retained results were moved without changing their bytes to
`results/legacy/`, preserving their old subdirectory structure. Their old path
strings and READMEs are historical provenance, not current invocation instructions.
Generated ElementPositions Python files inside those archived run directories
are output artifacts; executable analysis tools are in `scripts/`.

Migration initially preserved 2,056 historical files byte-for-byte. OPALX–IMPACT's 23
restored files were independently checked against the cleanup archive hashes.

On 2026-10-03, a subsequent user-approved disk cleanup removed 34 raw H5 files:
12 old primary (`c0`) snapshots from experiment 2's `timing_pairs` and
`timing_pairs_100`, 20 field/density/potential dumps from experiment 4's legacy
`track-e-p/data`, and two converted field-debug H5 files from `shared-data`.
This reclaimed about 1.50 GiB. No new backup was made; regenerate these raw
diagnostics before reanalyzing them. Witness trajectories, inputs, scripts,
plots, comparison tables and recent local runs were retained. Exact deleted
paths are recorded in `tmp/sandbox-cleanup-20261003/deleted-raw-output-files.txt`
relative to the repository root.

| Previous workspace | New location |
|---|---|
| `reduced-order-model` | `experiment-1` |
| `track12particles` | `experiment-2` |
| `TestParticleOrbit.dat`, reference presentation | `experiment-2/inputs` |
| `cain-opalx-reduced-order-model` | `experiment-3` |
| `track-e-p/fort98.txt` | `experiment-3/inputs/fort98.txt` |
| `track-e-p` and `track-e-p-new` | `experiment-4` |
| Shared `python/` utilities | `scripts/` |
| Shared generated `data/` | `experiment-4/results/legacy/shared-data` |
| `note/` | `note/` (kept outside experiment results) |

Optional cylinder diagnostics now live in `experiment-4/scripts/`:
`plot_cylinder_crossings.py`, `compare_cylinder_crossing_histograms.py`, and
`convert_efield_dumps_to_h5.py`. They are standalone tools, not called by `run.sh`;
their defaults read historical `results/legacy` data, and plots go to
`experiment-4/results/figures`. Their historical retirement-based defaults do not
describe the current geometric lifecycle. Use `--help` to inspect their options.

The six unused shared legacy tools (the GUI, its three helpers, the ASCII/H5
comparison, and the spherical-field movie) were archived before removal at
`tmp/sandbox-cleanup-20261003/shared-legacy-scripts.tar.gz`, relative to the
repository root. Archive contents were verified byte-for-byte before deletion.
This is a local recovery archive, not part of the runnable experiment workflow.

The earlier `RETIRE_TIME` production inputs are preserved under
`experiment-4/inputs/legacy/`; they are not silently rewritten or used by the
current launcher. The current dynamic-window deck supplies `--full`.

The full-population timed input now explicitly sets `ORBITTHREADER=FALSE` on
both passive witness beams, matching the working Track-12/production setup.
This reuses the primary design-orbit map; source geometry, charge, field model,
birth times and time steps are unchanged. The original input is preserved at
`experiment-3/results/legacy/opalx/gamma_gamma_pairs_timed.original.in`.

The scientific interpretation remains in
[the physics/validation report](BEAMBEAM_PHYSICS_AND_VALIDATION.md).
The [LaTeX note](note/bb-note.tex) references its relocated figure assets.
The accepted October CPU/A100 comparison remains separately under
`tmp/beambeam-cross-architecture-20261002/`; nothing there was moved or rerun.

## Scientific helper checks

These optional Python checks validate the relocated analysis/input tools, not
the short CI regression suite:

```bash
PYTHONDONTWRITEBYTECODE=1 ~/.venv-h6/bin/python -m unittest discover \
  -s sandbox/experiment-1/scripts -p 'test_*.py'
PYTHONDONTWRITEBYTECODE=1 ~/.venv-h6/bin/python -m unittest discover \
  -s sandbox/experiment-2/scripts -p 'test_*.py'
```
