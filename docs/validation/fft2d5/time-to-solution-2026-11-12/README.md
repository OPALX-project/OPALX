# FFT2D5 kernel-time distribution on LUMI

The figures compare the means-only version with the subsequent compact symmetric
moment reduction on one MI250X GCD. Both use 100,000 matched DBA particles,
a 64×64×512 mesh, 3 ps steps, OPEN longitudinal fields and FFT batch size 64.
The 99-step window starts at the second reference-path segment build and ends
at the last traced kernel. CSV files contain summed execution times by function;
JSON files record the window, executable hashes and temporal activity.

These historical measurements include transform fusion, which was subsequently
removed from PR #682 after its whole-tracking gain was below 1%. The unfused
branch was correctness-tested locally; these are not new measurements of it.
They also predate master's rest-frame temperature diagnostics. The merged
implementation retains those three sums: its compact/dense accumulators now
contain 46/76 doubles rather than the measured 43/73.

| Figure | LUMI profile job | Trace PID |
|---|---:|---:|
| `before-compact.png` | 23687487 | 93994 |
| `after-compact.png` | 23688478 | 83806 |

Raw traces remain under the corresponding `means-only-final-validation-20261112`
and `compact-symmetric-validation-20261112` directories in `~/opalx/` on LUMI.
Only kernels were traced. The remaining time interval is unclassified and does
not measure GPU occupancy or establish idle time.
