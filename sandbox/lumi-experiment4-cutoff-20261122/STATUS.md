# Lumi loss-producing Experiment 4 — 22 November 2026

Goal: rerun the full loss case with the new longitudinal witness cutoff and
compare immutable-ID trajectories, particularly positron oscillations.
No commit or push. No physics input changes or source changes planned.

Plan: verify the already-tested Lumi cutoff binary and saved old binary;
copy the exact local cutoff-20261122-n1 deck and input files; submit independent
one-GPU cutoff and old-binary control jobs with --info 4; fetch data and audit
locally; plot all witnesses and selected fixed-ID positron histories, with
old GPU and previous CPU comparisons. Preserve all previous results.

Physics: evolving 10k Gaussian primaries, 1297+1297 simultaneous witnesses,
16^3 mesh, IP0.17m, 32cm element, 15cm radius, fixed IP±10mm longitudinal field
window, expanding transverse domain. DT1ps/1fs/1ps with550/25001/1025 steps.
The old local pre-cutoff run aborted at step25566 before losses. The local
cutoff run finished26576 steps at1.600001ns with643e-/731e+ survivors and no
repeated positron oscillations. This is NOT the rigid exact-timed Track12 model.

Remote source: /flash/project_465002453/opalx-beambeam-recovery-20261120
Existing validated binaries: /flash/project_465002453/bb-cutoff-track12-20261122
New runroot: /flash/project_465002453/bb-experiment4-cutoff-20261122
Each case retains its Slurm script/logs, input hashes and execution fingerprint.
Only postprocessing runs locally, using ~/.venv-h6.

Submitted cutoff job23928245 and matched old-binary job23928246, each oneGPU
small-g/4CPUs/OMP1/--info4. Exact5 input hashes verified remotely before both
submissions. Deck SHA9dd216372a61f377bbd32ceaf29221442350489a09de456d061fb64e932a191f
matches the previous local cutoff run; versus previous failed CPU run only an
output-policy comment differs. Both GPU cases use identical deck bytes.
Cutoff executable SHAf44bf28796a8134cc96d8a3afbb5e8a62a90adf1be72c152dd28123266da73c3;
old SHA31e97e63d5208ef210c4188b717b48ccd739fba37874f38735e08754372d1307.
SSH initially delayed several minutes, then recovered without configuration or
host-key changes. No duplicate submissions. Remote source remains intentionally
uncommitted; no pull/reset/rebuild. Local analyze.py prepared to audit and plot
old/new GPU and compare common steps with previous CPU old/cutoff results.

Both jobs initially pending Priority with no estimated start. Requested move
of existing cutoff job23928245 to dev-g (same oneGPU/4CPU/2h limit); short
validation fits dev-g's3h limit. Old control remains small-g. No duplicate job
or input/executable change; execution fingerprint records actual partition.
First partition-only update rejected time limit; explicitdev_g QOS rejected
as invalid for job. Read prior successfulbuild23927919: dev-g/QOSnormal/30min.
Retry uses those verified settings for existingcutoff23928245. Both remained
pending small-g after rejected updates. No script executed prematurely.
Final update found job alreadyrunning: BOTH started07:31:52 in small-g,
cutoff nid005026 andold nid005033. No partition/time change applied; bothretain
2h limit. Local compare_prefix smoke against previousCPUold/cutoff: allcommon
saved IDs/clocks andphasevalues exactlyequal forbothspecies. Five synthetic
loss-audit tests pass; submit.sh syntax andanalyze.py AST pass. Matplotlib
defaultcache permissionwarnings during smoke; finalanalysis willusecampaigncache.

Both jobs finished: cutoff23928245 COMPLETED0:0 (Slurm8m01s, timed477.98s),
old23928246 FAILED28:0 (Slurm7m54s, timed471.36s, run exit-status156).
Old reproduces witnessCIC outside-storage abort justafter completedstep25565,
t590.001ps. Cutoff continues through aperture losses and completes. No job
updates applied; both executed small-g with original2h limit. Outputs fetched
locally (excluding cores) with rsync-rz no-perms/owner/group; original data
remain untouched. Local analyze.py running with campaign-local caches and
OMP/BLAS1; it checks input/executable hashes, audits each case and compares
common immutable-ID saved histories with GPUold and CPUold/cutoff baselines.

## Completed analysis

Both GPU input manifests and executed-binary hashes verified. All analyses,
CSV comparisons and plots generated locally with analyze.py; no tolerances
changed. Final remote source hashes still match the validated cutoff patch,
and tracked remote changes remain exactlythe3 intended cutoff files. Queue
empty for bothjobs. No source edits, commits or pushes during this campaign.

| Measure | GPU pre-cutoff control | GPU cutoff | Previous CPU cutoff n1 |
| --- | ---: | ---: | ---: |
| Completed integration steps | 25566 attempted / last completed25565 | 26576 | 26576 |
| Last saved witness time [ps] | 584.001 | 1600.001 | 1600.001 |
| Final saved electrons / positrons | 1297 / 1297 | 643 / 731 | 643 / 731 |
| Witness losses e- / e+ | 0 / 0 | 654 / 566 | 654 / 566 |
| Loss snapshot intervals per species | 0 | 49 | 49 |
| Repeated positron oscillations near overlap | absent | absent | absent |
| Largest positron px reversal count near overlap | 1 | 1 | 1 |

GPU cutoff audit: zero duplicate/reappearing/replacement IDs, nonfinite
components or causality violations. At all49 loss intervals per species the
survivors are ballistic: maxdeltaP=0, maxdriftresidual1.28e-15m. Same-side
outside-window intervals (99932electron/102297positron) also have exactlyzero
momentum change. No post-loss trajectory jump is reproduced. PositronID706
turns once near overlap; no positron has repeated px reversals within20ps.
This does not establish convergence of the coarse16^3 model or bound omitted
impulse from hard cutoff. It is not the rigid/timed Track12 setup.

Matched GPU old->cutoff over4940 common saved container snapshots:
maxlabpositiondifference1.2490009027e-15m, maxmomentumdifference2.1670688906e-15mc,
identical IDs andclocks. Old abort leaves no later loss data for comparison.
The cutoff therefore preserves the pre-window-exit histories to roundoff.

CPU cutoff->GPU cutoff over5144common snapshots: maxpositiondifference
0.0004541508971m (0.454151mm), maxmomentumdifference0.001508411731mc.
Two snapshot population differences: electronstep26339 hasoneCPU-onlysurvivor;
positronstep26379 hasoneGPU-onlysurvivor. Finalcounts agree but phase space
and intermediate loss timing are NOT identical. CPUold->GPUold already has
maxpositiondifference7.888603379e-6m andmaxmomentumdifference0.001508411731mc
over4940common snapshots beforeabort. CPU/GPU difference therefore predates
the cutoff; its detailed source/solver/RNG attribution is not established
by this study. The Gaussian primaries are generated at runtime (unlike the
fixed FROMFILE Track12 control). No statistical error/tolerance is inferred.

Plots (relative toIP0.17m and nominalIP/c):
- comparison/all_trajectories.png: all sampled witnesses, old/new columns;
  top panels use different ranges because old stopped early. Bottom panels
  show all positrons near overlap, where old/new trajectories overlay.
- comparison/positron_oscillations.png: full saved near-overlap histories for
  six fixed IDs selected by highest momentum reversal/axis-crossing counts;
  dashedold/solidcutoff. No nearest-neighbour joining or old-run extrapolation.
- cutoff/loss-audit/loss_trajectory_audit.png: populations, causality atlosses,
  selectedpositron histories through fulltail.

Reproduce after fetching both outputs:
`MPLCONFIGDIR=/tmp/exp4-mpl XDG_CACHE_HOME=/tmp/exp4-cache OMP_NUM_THREADS=1 ~/.venv-h6/bin/python sandbox/lumi-experiment4-cutoff-20261122/analyze.py`
Detailed metrics: comparison/summary.json and threeper-snapshotCSV files;
cutoff/loss-audit/{summary.json,*_loss_intervals.csv,*_identity_summary.csv}.
Scripts andscheduler logs remain ineachremote run directory.

## Why Experiment 4 differs from Track-12 (diagnostic follow-up)

New reproducible read-only analysis: compare_models.py -> model-comparison/
{summary.json,primary_evolution_and_births.png,*_envelope.csv,*_births.csv}.
No tracking/source/deck changes. Experiment4 primary H5 is the previous CPU2
run (newGPU primary output disabled); GPU witness envelopes are measured from
the newcutoff run. Track12 primary evolution is a free-drift prediction from
its exactfixedfile, not a measured GPU history. Scope explicit in plot/report.

Track12 BBRIGID=TRUE suppresses source collective kicks; actual fixedfile
sigma_x=1.91903um, sigma_y=1.92068um (nominal1.94433um), rms_px~4.8e-10mc.
Predicted width change over9.012ps is<1e-18m. The flag does NOT freeze positions:
ordinary drift and external tracking remain active (OpalBeamBeam Doxygen).
Experiment4 defaultBBRIGID=FALSE, nominalentrance sigma1.94433um and
sigma_xprime=sqrt(epsilon/beta)=2.67445mrad. No focusing elements or specified
x-px/y-py correlations: defaults0. Unforced drift over17cm alone predicts
sigma(IP)=sqrt(epsilon*beta+L^2*epsilon/beta)=454.661um; savedCPU at567.050ps
shows435.529/444.206um, sigma_z588.927um. Growth therefore does not require
a self-field instability. The small nominal spot is specified at the source,
not transported as a converging beam to a waist atIP. Width~227x Track12.

Experiment4 injects all1297/species simultaneously atIP+inputz at
t=IP/c-16.747ps. Each primary centroid is nominally5.0206mm fromIP then.
Track12 artificialwitnesses have px=py=0, pz=1.732051mc and sixexplicitbirths
at -3.00208,-1.99938,-1.00269,0,+1.00269,+1.99938ps relative nominaloverlap,
at correlatedz offsets -0.9,-0.5994,-0.3006,0,+0.3006,+0.5994mm.
Thus they are NOT allborn exactly atIP. Sigma_z/c~2.00138ps.
Experiment4 actualpositrons: medianv_perp~0.498533c, medianradius2.127um at
firstsaved sample ->2.47117mm atnominaloverlap; free-flight input estimate
2.50295mm over16.747ps. They mostly leave the central core before overlap.
Track12 begins with no transverse witness momentum near a narrow rigid source.

These initial-condition/model differences strongly explain weak/absent repeated
focusing inExperiment4; they do not isolate relativecausal contributions of
sourcewidth,timing,witnessphase space and16^3vs256x256x128mesh. MatchedGPUold/new
commonhistories agree~1e-15m, so cutoff is not origin ofpre-overlap difference.
Next proposed (NOT run): separately control source preparation/waist and
witnessbirths/initialphase space, then mesh convergence. BBRIGIDalone cannot
remove Experiment4's initial divergence. No production model changed; no
commit/push and no LaTeX edit in this follow-up.

## Pre-commit verification

User subsequently authorized a local commit, not a push. The BeamBeam Lumi
checkout has exactlythree modified source/test files, each byte-identical to
the local file (SHA256 verified). Remote Experiment4 submit.sh and input hash
manifest also match locally. No remote-only BeamBeam source fix exists.
The separate unrelated ~/opalx master checkout hasfour pre-existing tracked
edits and unrelated experiment directories; it was not changed or merged.
Its tracked binary diff and base/status were backed up locally as ignored
unrelated-lumi-master.patch and unrelated-lumi-master-status.txt. That backup
and unrelated files are excluded from the BeamBeam commit.

clang-format dry-run and git diff --check pass. Final CTest four-entry subset
passes with PRTE_MCA_rmaps_default_mapping_policy=slot:OVERSUBSCRIBE and
PRTE_MCA_hwloc_default_binding_policy=none. Default launch failed in mapping
before tests started; no assertions/tolerances changed. Explicit mpiexec
--oversubscribe --bind-to none passes full19-test MultiContainer suite on n1,
deletion/H5 regression on n2, field-service suite on n1/n2/n4 (expected4-rank
mirror skip on n1/n2). Five synthetic loss-audit tests pass. Generated data,
PDF, local AGENTS.md/HANDOFF.md and historical PYTHON_SCRIPTS.md are excluded.
