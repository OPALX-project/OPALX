# Lumi Figure 6(b) cutoff rerun — 22 November 2026

Goal: rerun the exact previous one-GPU, 400000-primary, 256x256x128,
1501-step Track-12 case with the uncommitted longitudinal witness cutoff.
Do not commit or push. Preserve all previous outputs and executable.

Source checkout: /flash/project_465002453/opalx-beambeam-recovery-20261120
Verified base: e24dc92185198c5f65718e057400e6c6398dce56, branch
271-implement-interation-point-element. Tracked tree initially clean;
untracked build-hip retained. Dirty master ~/opalx is out of scope.
New campaign: /flash/project_465002453/bb-cutoff-track12-20261122
Baseline: /flash/project_465002453/bb-validation-e24dc9218-20261120/track12
Local baseline: sandbox/lumi-validation-20261120/track12

Only cutoff production header/algorithm and field-service test are transferred.
No local ID/H5-test changes, LaTeX, data, or unrelated files are applied remotely.
Old HIP executable will be copied to pre-cutoff/ before the focused patch and
incremental Release HIP rebuild. Build/new executable hashes and the exact patch
are retained. Slurm scripts and scheduler output live in each run directory.

GPU build/n1 tests -> n2 field-service tests -> one-GPU Track-12 (--info 4).
Same LUMI25.03/AMD6.3.4/Cray-MPICH8.1.32, existing site toolchain/link settings.
Keep GPU-aware MPI enabled and GPU IPC disabled, respecting Slurm ROCr masks.
Copy exact baseline input and deterministic primary files, verifying SHA256
before execution. Fetch results and postprocess only on the Mac with ~/.venv-h6.

Expected: all Track-12 witnesses remain inside the longitudinal force window,
so the ordinary gather path remains in use. Compare immutable H5 IDs, clocks,
all six phase-space components and frame metadata directly; do not reconnect
trajectories by proximity or assume bitwise GPU repeatability. If differences
exceed ordinary GPU run-to-run variability, use the saved old binary as control.

Build job 23927919 submitted from runroot/build (dev-g, 1 GPU, 20 CPUs).
Patch applied after git apply --check; all three remote source hashes match
source.sha256. Remote tracked changes are exactly the three intended files.
Five copied input hashes verified against baseline before submission.
Old executable saved with baseline SHA256
31e97e63d5208ef210c4188b717b48ccd739fba37874f38735e08754372d1307.
Build/n1 gate completed successfully. All three applicable field-service tests
pass on one GPU, including the new cutoff/re-entry/nonuniform-CIC case; the
four-rank-only mirror test is expected to skip on one/two ranks.
Two-GPU test job 23927942 submitted from runroot/unit-mpi2.
MPI2 job 23927942 completed exit0: all three applicable tests pass on both GPU
ranks. Track-12 job 23927947 submitted from runroot/track12 (small-g, 1 GPU,
4 CPUs, OMP1, --info4). All full-run inputs remain byte-identical to baseline.
Track-12 confirmed running on nid005026, at least step300 at elapsed2m22s.
GPU2 binding confirms distinct PCI buses c9:00.0 and d9:00.0, one visible GCD
per task. New executable SHA256:
f44bf28796a8134cc96d8a3afbb5e8a62a90adf1be72c152dd28123266da73c3.
Patch SHA256: ae3484a4948b94604485dc967b05afd657ae28e1235979ceb50b620d7a4c1595.
Track-12 job 23927947 completed exit0; all1501 steps. Slurm elapsed8m01s;
measured srun execution472.25s versus previous469.86s. Build23927919 completed
in2m20s and MPI2 gate23927942 in15s. Current state: fetching completed particle
data, then local ID-by-ID/CAIN/manufactured comparison. No source/kernel changes needed
for the HIP build. The saved pre-cutoff executable and old results are intact.

Build and MPI2 logs fetched locally. Archive-mode rsync returned directory
chmod warnings on the Mac; retry with -rz --no-perms --no-owner --no-group
succeeded. Particle data are not fetched until completed. Local analyze.py
checks exact input hashes, IDs/clocks/phase space, canonical bytewise hashes,
CAIN/manufactured errors, runtime, and executed-binary fingerprint. Baseline
integrity passes: max |s-IP| <1.65mm, zero outside-window records.

New Track-12 fetched and analyzed locally. Input hashes, IDs, clocks and
13,008 actual H5 witness records match in identity; all 12 witnesses remain
finite and inside the longitudinal window. Outputs are not byte-identical:
max lab-position difference3.024185e-10m, max momentum difference1.271378e-6mc,
RMS position3.759128e-11m. Curves visually overlap; CAIN/MSOL errors unchanged
to six displayed significant digits. First plotting attempt needed braces in
Matplotlib \mathbf labels; fixed, analysis and overlay PNG now succeed.

Additional old-binary repeat control job23928050 submitted from runroot/old-repeat
with the exact saved31e97e63... executable and identical5 input hashes, oneGPU,
small-g/OMP1/--info4. Purpose: quantify normal repeatability before attributing
the small nonzero difference to the cutoff revision. No tolerances changed.
Local analyze.py --control RUN will compare baseline->old-repeat and
old-repeat->cutoff directly by immutable IDs. Current state: control submitted;
main cutoff run and GPU tests completed successfully. No source changes after
the validated patch. Original Figure6b PNG and baseline data are untouched.

## Final result

All four jobs completed exit0: build/n1 23927919, GPU n2 23927942,
cutoff Track-12 23927947, saved-old-binary control 23928050. Final remote
source hashes still match the validated patch; only the three intended tracked
files are modified. No commit or push was made.

The control completed all1501 steps in475.05s. All five control input hashes
match, its executable hash is the original31e97e63..., all witness IDs/clocks
match, and phase space stays finite. Direct comparison of13,008 actual H5
records (positions transformed to lab coordinates) gives:

| Comparison | Max position difference [nm] | Max momentum difference [m_e c] | RMS position difference [nm] |
| --- | ---: | ---: | ---: |
| Original vs cutoff | 0.302419 | 1.271378e-6 | 0.0375913 |
| Original vs unchanged old executable repeat | 0.982497 | 2.734897e-6 | 0.110604 |
| Old executable repeat vs cutoff | 0.681065 | 1.934992e-6 | 0.0780638 |

The cutoff/original difference is smaller than the observed old-executable
repeat difference, consistent with GPU run-to-run/environment variability.
One control is not a statistical characterization; no tolerance was changed
and no Monte Carlo resampling occurred (the primary file is identical).
Runtime472.25s vs original469.86s (+0.51%); old control475.05s. These single
measurements do not establish a performance change.

Manufactured relative L2 errors, original -> cutoff:
x 0.00853737948661 -> 0.00853738358529;
s 0.00117302207984 -> 0.00117302302728.
CAIN relative L2 x 0.0948284803983 -> 0.0948284763466;
s 0.0296731387292 -> 0.0296731381476.
All reference errors are unchanged to six significant digits. All witnesses
remain within1.65mm ofIP, so this validates inside-window preservation, not
the physical accuracy of hard truncation at10mm. GPU n1/n2 regressions cover
outside-window zero fields, endpoints, re-entry and mixed nonuniform CIC.

Reproduce local analysis after fetching both completed runs:
`~/.venv-h6/bin/python sandbox/lumi-cutoff-20261122/analyze.py --control sandbox/lumi-cutoff-20261122/old-repeat`
(Use campaign-local MPLCONFIGDIR/XDG_CACHE_HOME if needed.)
Comparison plot and machine-readable metrics:
`track12/cutoff-comparison/previous_vs_cutoff.png`, `summary.json`,
`h5_differences.csv`, `baseline_vs_old_repeat.csv`, `old_repeat_vs_cutoff.csv`.
Updated CAIN/manufactured plots are in `track12/results/`; the previous
Figure6(b) assets and LaTeX figure reference have not been overwritten.
