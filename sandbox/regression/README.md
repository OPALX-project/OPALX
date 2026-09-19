# Active sandbox regression outputs

This directory holds compact baseline data and locally generated figures for
the active broad sandbox checks. It is deliberately separate from the archived
LaTeX publication assets in `../attic/notes/`.

Run the lightweight checks from the repository root with:

```bash
source ~/.venv-h6/bin/activate
python sandbox/python/run_sandbox_regressions.py
```

The combined wrapper consumes historical comparison outputs rather than
regenerating them by default. In the current checkout the `fields` component is
available and passes, while the complete command is blocked until
`sandbox/OPALX-IMPACT/extracted_graph_data.csv` and
`sandbox/track-e-p/gamma_gamma_pairs-2_c0.stat` are regenerated or restored.
This does not replace the maintained focused tests listed below.

The default files are:

- `sandbox_regression_baseline.json`: accepted compact metrics;
- `current_metrics.csv`: latest flattened comparison, generated locally;
- `figures/`: locally generated diagnostic plots.

Updating a baseline accepts new numerical behavior and therefore requires a
review of the physics change and its tolerances. The focused manufactured,
timed-emission, and MPI tests documented in
[`../BEAMBEAM_PHYSICS_AND_VALIDATION.md`](../BEAMBEAM_PHYSICS_AND_VALIDATION.md)
remain the primary acceptance evidence.
