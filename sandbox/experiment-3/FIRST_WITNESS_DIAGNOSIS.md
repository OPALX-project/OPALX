# First-witness gather failure — 2026-11-20

Status: coverage repair implemented and smoke-tested; no commit or push.

## Repair validation

The AfterEmission phase checks the witness transverse envelope in source axes.
If it exceeds the solved domain, it rebuilds the combined non-shrinking envelope
and repeats the primary/copy solve with the saved pre-emission step state.
This includes the newborn fractional half-drift position but not its charge.
Fixed longitudinal bounds, source model, units, and numerical tolerances are
unchanged. A device reduction and MPI collective reject any invalid local CIC
stencil before gather halo exchange. No particle-array host copies are added.
An extra solve may change floating-point results on the newly enlarged mesh;
this is intentional and replaces an invalid memory read, not a rounding fix.

GCC build and field-service tests pass on 1, 2 and 4 ranks, including collective
rejection when only rank 0 contains an uncovered particle and successful recovery
after restoring positions. The four-rank non-slab mirror test also passes.

Experiment 2 smoke results: `local-20261120-221513-ly7dacao` (one rank),
`local-20261120-221600-mz9s5nhy` (two ranks). Both preserve the October smoke
manifest and have maximum coordinate delta 3.34e-16 m and momentum delta
3.71e-10 mc. See each run's `former_smoke_*` CSV/JSON files.

Experiment 3 results: `local-20261120-221523-66s1bzsm` (one rank),
`local-20261120-221540-ptq9gur7` (two ranks). Both complete all 160 steps,
exact 1297+1297 birth schedules and finite phase space. One rank performs two
post-emission re-solves. First-electron momenta after the first kick:

| ranks | px [mc] | py [mc] | pz [mc] |
|---|---:|---:|---:|
| 1 | 1.07882238 | -0.75428244 | 0.12565862 |
| 2 | 1.07881516 | -0.75428139 | 0.12565843 |

The invalid ~5865 mc first kick is gone. October experiment-3 data are not a
valid trajectory acceptance baseline: birth-time gathers were uncovered there
too. One-rank historical deltas are saved in the new run's comparison CSVs.
Cross-rank later emitted ID sets differ (first detected at dump Step#10), so
the strict ID-based comparison deliberately refuses a full trajectory match;
only birth inventories, clocks, finiteness, and the first particle are compared
across ranks here. Primary Gaussian samples also vary with rank/backend.
No tolerance has been relaxed. Fine-grid accuracy and accelerator validation
remain separate work.

Reproduce local smokes with `OPALX_EXE=.../build_gcc/src/opalx` and
`sandbox/experiment-2/run.sh --smoke`, `sandbox/experiment-3/run.sh`.
For two ranks, use an executable wrapper invoking
`mpiexec --oversubscribe --bind-to none -n 2 .../build_gcc/src/opalx "$@"`.
Use `sandbox/scripts/compare_smoke_phase_space.py OLD NEW` for strict
ID/clock-matched comparisons. It rejects unmatched IDs rather than guessing.

## Original diagnosis

### ID audit follow-up

The differing one-/two-rank ID sets do not indicate renumbering. IPPL allocates
birth IDs as `rank + nranks * local_birth_counter`; emission capacity and rank
distribution determine which file row receives each ID. All 38 dumps of both
witness species were checked: zero duplicate IDs and zero previously emitted
IDs disappearing in either run. Final IDs are 0..1296 on one rank; maxima are
1736 (electrons) and 1731 (positrons) on two ranks, with 1297 distinct IDs each.
Do not overwrite birth IDs to force equality across different MPI allocations.

The timed-file regression now records each file row's initial ID, migrates the
particles after every emission batch, and checks that subsequent batches and
migration preserve that row/ID mapping. All ten emitted-file tests pass on
1, 2 and 4 ranks. Migration required fixing the test fixture's mesh/layout
lifetime: these borrowed objects now outlive the particle container.

Evidence run: `results/local-20261120-220751-7r4c2c0s/opalx.log`.
An opt-in temporary trace was added around BeamBeam gathering and Boris kicks,
then removed from source after collecting the evidence. GCC serial backend,
one MPI rank, existing 160-step experiment-3 smoke profile, `--info 4`.
The trace introduces host copies only for diagnostics; no numerical kernels,
units, tolerances, dependency sources or algorithms were changed.

At global step 122, start time 112.462423087 ps:

| Stage | Electron result |
|---|---|
| Emission | P = (1.0734722114, -0.7535157123, 0.1256345426) in mc; E=B=0 |
| Fractional first timestep | 0.7740923992 ps |
| Position at gather, source frame | (1.108170824, -0.154965226, 1.305625527) mm |
| Field after gather | E = (4.2157e10, -1.0151e11, -3.8605e13) V/m; B ~ (-1.0385e16, -1.0329e16, -1.0341e16) T |
| Momentum after first kick | (5865.0113785, 5835.1457756, 5840.7755116) in mc |

The fixed field domain has x=[-32.106284,31.954121] micrometres,
y=[-32.043189,32.017216] micrometres, and 16 cells/nodes per direction.
The electron is outside the mesh before gathering. With the actual mesh origin
(-3.4241630846e-5,-3.4178536164e-5,-9.3819798874e-3) m and spacing
(4.2706936886e-6,4.2706936886e-6,1.3333333333e-3) m, IPPL's expression
`l=(R-origin)/spacing+0.5` gives (268.00044,-27.78269,8.51570).
Integer indices (268,-27,8) cannot address a 16^3 field even with halos.
`ParticleAttrib::gather` converts these to unsigned stencil indices and directly
indexes the view without a coverage check. The resulting values are invalid
memory reads, not physical field values. The positron also lies outside the
mesh and gathers invalid values (including Bz about 480 T).

Cause: `computeWindowSelfFields()` builds the combined envelope and solves
BEFORE emission, when witnesses are empty. `AfterEmission` calls only
`gatherFieldsToWitnessContainers()`; newly born particles were not included in
the solved envelope. Migrating particles onto the field layout does not expand
the solved domain. The next step expands the mesh, but the first kick is already
corrupted. Historical smoke runs also had uncovered birth positions; a benign
old memory value is not evidence that this was safe before the merge.

Recommended repair (not implemented): establish coverage of all newly emitted
witnesses before gathering. If emission extends the current domain, enlarge it
and re-solve the primary/copy fields at the same field-evaluation time, without
depositing witness charge. Add a collective OPALX-side CIC-stencil coverage
guard and single-/multi-rank regression with an initially empty witness container
whose first birth lies outside the primary envelope. Do not clamp, wrap, or
silently zero uncovered fields. Simply enlarging the geometry without re-solving
would reinterpret cached fields on the wrong mesh and is invalid.

Scope: this establishes the gather defect. It does not establish numerical
agreement for the enlarged mesh, or explain every historical/current trajectory
difference; primary RNG/backend differences remain a separate comparison issue.
