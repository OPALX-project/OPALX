# Sandbox Regression Checks

This directory keeps the BeamBeam sandbox checks runnable from one command.
The checks use pandas dataframes for tabular data and compare compact metrics
against `sandbox_regression_baseline.json`.

## What Is Covered

- `fields`: compares the `BeamBeam-static-1V` ASCII `rho`, `phi`, and `E`
  dumps in `data/sandbox/` against the analytic symmetric Gaussian BeamBeam
  manufactured solution.
- `pairs`: checks the gamma-gamma `e-/e+` FROMFILE momentum summaries and the
  final OPALX `.stat` outputs for `track-e-p/gamma_gamma_pairs-2.in`.
- `impact`: compares `sandbox/OPALX-IMPACT/drift-1.stat` with the digitized
  IMPACT reference in `extracted_graph_data.csv`.

## Run

From the repository root:

```bash
source .venv-h6/bin/activate
python sandbox/regression/run_sandbox_regressions.py
```

Run one check:

```bash
python sandbox/regression/run_sandbox_regressions.py --check fields
python sandbox/regression/run_sandbox_regressions.py --check pairs
python sandbox/regression/run_sandbox_regressions.py --check impact
```

Each run compares against the accepted baseline in
`sandbox/regression/sandbox_regression_baseline.json` and writes the latest
flattened metric table to `sandbox/regression/current_metrics.csv`.  The CSV is
the latest run output; the JSON baseline is the comparison reference.  Use
`--csv` only to choose a different latest-run output path:

```bash
python sandbox/regression/run_sandbox_regressions.py \
  --csv sandbox/regression/current_metrics.csv
```

## Regenerate Outputs First

If an OPALX executable is available, the harness can rerun the sandbox inputs
before checking metrics:

```bash
python sandbox/regression/run_sandbox_regressions.py \
  --run-opalx \
  --opalx-exe /path/to/opalx
```

The default executable path is `build_openmp/src/opalx`.  If only the unit tests
have been built, that executable may not exist yet.

## Overview PDF

Build a compact human-readable overview from the current metrics:

```bash
python sandbox/regression/make_overview_pdf.py
```

This writes `sandbox/regression/sandbox_regression_overview.tex`,
`sandbox/regression/sandbox_regression_overview.pdf`, and
`sandbox/regression/opalx_impact_drift_comparison.png`.  The overview
summarizes what each check compares, includes selected latest-run metrics,
shows the accepted-baseline git provenance, reuses the figures from
`sandbox/note` as visual context, and overlays OPALX with the digitized IMPACT
reference for the drift check.

## Accept New Results

Only update the baseline after the physics or numerics change is intentional:

```bash
python sandbox/regression/run_sandbox_regressions.py --update-baseline
```

The baseline JSON records the git branch, commit, dirty flag, and
`git status --short` output present when the baseline was accepted.

For algorithmic changes, record the expected effect on conservation, stability,
floating-point reproducibility, and OpenMP/MPI execution order in the change
description.

## Pair Input Utilities

The pair regression uses the default six-column FROMFILE distributions:

```bash
cd sandbox/track-e-p
source ../../.venv-h6/bin/activate

python convert_fort98_to_fromfile.py fort98.txt
python analyze_pair_momenta.py
```

For creation-time diagnostics, regenerate the optional time-preserving
FROMFILE variants and run the analyzer on those files:

```bash
python convert_fort98_to_fromfile.py fort98.txt --add-time
python analyze_pair_momenta.py \
  --electron-file gamma_gamma_electrons-t.fromfile \
  --positron-file gamma_gamma_positrons-t.fromfile \
  --plot
```

The `--add-time` files use the header `x y z px py pz t bin_number`, with
`t=ct/c` in seconds and `bin_number=1` for every particle.  The analyzer reports
time summaries in picoseconds only and writes `pair_time_species_summary.csv`,
`pair_time_pair_summary.csv`, and `pair_momentum_time_histogram.png` when both
species inputs contain `t`.
