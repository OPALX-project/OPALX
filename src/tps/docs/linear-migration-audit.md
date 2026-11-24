# First linear-map migration to TPS {#tps_linear_migration_audit}

Historical design audit of OPALX commit `34c835e96`, 2026-11-24, before the
linear migrations. The original proposal and code-path descriptions below record
that earlier state; they are not a current implementation inventory. For the
implemented storage, numerical backends, accepted Newton QR policy and validation
status, see the [Part 1 completion audit](part1-completion-audit.md).

## Implementation status

The first composition-only migration is now implemented in
LinearTransferMapBuilder, using opalx::maps::toTpsLinearTransferMap. Existing
finite-difference and Richardson construction is unchanged. The new unit tests
cover matrix/evaluation parity, noncommuting composition, chart guards, traversal
IDs, ambient truncation restoration, nonzero ownerless field tails, shared overlap
and nonunit-momentum drift. The 4D COF/eigenanalysis adapter is also implemented: it preserves the final
reference/residual, extracts a Jacobian and delegates to the original LAPACK
implementation. The Newton solver, full-tune winding and spectral path remain
unchanged.
The audit below records the historical design and proposed sequence. The current
completion audit is authoritative for remaining acceptance work.

## Recommendation

Start with the **ordered 6D segment-map composition in LinearTransferMapBuilder**,
called by OrbitThreader. Keep finite-difference map construction and existing
matrix-facing outputs. Use first-order context-bearing TPS maps to compute the
combined derivative. This replaces one real algebra operation and establishes the
coordinate adapter before changing tracking, storage, or eigensolvers.

The next small change should adapt the independent **4D COF map to eigenanalysis**.
Do not combine both coordinate migrations into one first change. Wrapping a COF
matrix alone is a valid smaller plumbing experiment, but the builder composition
pilot exercises TPS composition in the existing element-map workflow.

## Code paths at the original audit

| Path | What it does | Important boundary |
| --- | --- | --- |
| OrbitThreader::execute -> LinearTransferMapBuilder::build | Design-beam reference sampling, unique owner segmentation, perturbed-ray derivatives, ordered 6D product | OrbitThreader.cpp:288-305; LinearTransferMapBuilder.cpp:643-650 |
| OneTurnMap | Independent fixed-energy 4D directed-section return, central-difference Jacobian and optional phase winding | OneTurnMap.cpp:80-148, 150 onward |
| ClosedOrbitSolver -> CofCmd | Damped Newton fixed point; produces its own final finite-difference matrix; rank-zero analysis/reporting and broadcast | ClosedOrbitSolver.cpp:133-173; CofCmd.cpp:359-463 |
| LinearMapEigenAnalysis | Full coupled 4x4 spectrum and eigenvectors with LAPACKE_dgeev; packed-basis SVD and scaled eigenpair residual checks | LinearMapEigenAnalysis.cpp:33-137 |
| SpectralTunes | Separate two-ray, uniformly sampled, mean-subtracted Lomb--Scargle spectrum | SpectralTunes.cpp:18-64, 67-163 |

OrbitThreader is the orchestrator. LinearTransferMapBuilder constructs derivatives,
including Richardson refinement, and combines them. OrbitThreader attaches copies
to owners and stores the combined result. On RINGs it reuses the starting section
axes at return. A periodic reference return is not by itself a verified closed
orbit (OrbitThreader.h:108-118).

Builder composition uses the unique `result.segments` list. An overlap segment can
be attached to multiple elements, and ownerless intervals can contain field tails.
Rebuilding the product from element-owned copies risks duplicating overlaps and
omitting ownerless intervals. OpalBeamline's ordered attachment accessor deduplicates
owner copies but cannot recover ownerless intervals that were never attached.

## Coordinate contracts

### OrbitThreader/builder: six-dimensional differential transport

The implementation in LinearTransferMapBuilder::coordinates defines

\f[
 X=(x,x',y,y',\zeta,\delta),\qquad
 x'=\frac{\mathbf p\cdot e_x}{\mathbf p\cdot e_s},\quad
 y'=\frac{\mathbf p\cdot e_y}{\mathbf p\cdot e_s},
\f]
\f[
 \zeta=-\beta_0c(t-t_0),\qquad
 \delta=|\mathbf p|/p_0-1.
\f]

- x, y and zeta use metres; slopes and delta are dimensionless.
- x/y are transverse displacement from the boundary reference position.
- Momentum is mechanical p/(mc); p0 is its reference magnitude and
  beta0=p0/sqrt(1+p0^2). No vector-potential canonical momentum is supplied.
- Each boundary owns its frame, reference position/momentum and arrival time.
  Frames use minimum-rotation/Bishop transport, with the original section normal
  restored at a ring return. Returned momentum need not align with that normal.
- Matrices are differential maps, `dX_out=M*dX_in`; no affine offset is stored.
  Richardson levels improve the finite-difference derivative, not polynomial degree.

For TPS, use degree one, named coordinates in this exact order, numerical chart
scales all one in the declared units, `CoordinateConvention::mechanical`, and no
canonicalPairs. The six local variables represent perturbations in these reported
coordinates. An output row is `sum_j M(i,j)*variable(j)`, with zero constant by the
existing differential-map contract, not by asserting exact orbit closure.

Construct one chart per boundary and reuse it unchanged for both adjacent maps;
MapContext joins require exact metadata equality. Define frame.axes explicitly as
rows `(e_x,e_y,e_s)` expressed in lab coordinates, i.e. local=R*(lab-origin).
Keep referenceMomentum in p/(mc) and referenceTime in seconds. Retain the original
LinearTransferMapReference companion metadata, including compensated low parts,
path, integration method, finite-difference settings and diagnostics.

Do not assume every boundary's encoded reference vector is all zero: at a ring
return, its stored frame can disagree with the returned momentum. Preserve its
reference slope values when encoding the chart. Preserve the distinction between
perturbations and absolute encoded coordinates; this first stage is differential
composition, not a complete absolute-coordinate particle tracker.

### COF/OneTurnMap: four-dimensional fixed-energy return

\f[
 u=(x,u_x,y,u_y),\qquad u_x=p_x/(mc),\quad u_y=p_y/(mc),
 \qquad u_z=+\sqrt{u_0^2-u_x^2-u_y^2}.
\f]

Positions are metres, the momentum entries and fixed magnitude u0 are dimensionless.
The initial state is on section z=0 at T0; each perturbed trajectory reaches its own
first armed forward crossing, with independently determined arrival time. The
shared ray tracker supplies fields and losses. Return momentum is not projected
back onto the launch energy shell; integration drift is reported separately.

This is not the transverse block of the 6D slope map. Even at fixed energy,
`x'=u_x/u_z` is a nonlinear coordinate change. Different section/frame conventions
and longitudinal constraints must also be handled before comparing these maps.
A p0=1, zero-transverse-momentum example alone can hide an incorrect conversion.

For the later TPS adapter, use a mechanical 4D chart centered at the solved u*,
with equal input/output reference u*, the same fixed section, and unit chart scales.
Preserve the measured residual:

\f[
 r=T(u_*)-u_*,\qquad
 \widehat T_1(\xi)=r+M\xi,\qquad \xi=u-u_*.
\f]

Never silently set r to zero. With nonunit chart scaling D, the coefficients would
instead be D^-1*r and D^-1*M*D. Eigenanalysis already uses its independent
conditioning matrix S through S^-1*M*S; do not apply normalization twice.

Keep T0 and return duration as companion diagnostics. For an autonomous 4D section
chart, leave optional referenceTime unset rather than claim all trajectories share
a fixed output time. Do not discard time from the separate 6D boundary chart.
A nonzero residual also means default preserveOrigin composition rejects this
4D map; the first adapter is for local analysis, not repeated-turn composition.

Both existing coordinate systems are noncanonical. A nonzero canonical-J residual
does not alone establish unphysical tracking, and a unit determinant does not prove
symplecticity. Canonical normal forms require a separately defined transformation.

## Tune-analysis findings

- COF eigenanalysis consumes ClosedOrbitSolver::Result::matrix directly. The
  solver's OneTurnMap overload delegates to a general callback and computes its
  own eight probes; it does not call OneTurnMap::jacobian. A successful result
  carries the final Jacobian and freshly verified residual. Reusing those avoids
  another eight returns in a TPS adapter.
- Fractional modes from the positive-imaginary eigenvalues are in (0,1/2), with
  complementary branch 1-q. They do not identify planes or integer winding.
- OneTurnMap::fullTunes supplements the accurate matrix with continuous shadow-ray
  phase winding. It supports stable uncoupled planar optics and rejects coupling
  along the trajectory even when the endpoint map appears decoupled. Replacing it
  with eigenvalue arguments would lose information.
- SpectralTunes has a separate sampling/turn convention: nominal step counts,
  finite launch amplitudes and grid spacing 1/(4*legacyTurns). Its output is not
  generally identical to zero-amplitude directed-section eigenphases. Retain it as
  an independent diagnostic; do not replace the spectrum with matrix eigenvalues.
- TPS currently supplies no eigensolver or canonical normal-form service. Preserve
  LAPACKE_dgeev/dgesvd and the closed-orbit solver's pivoted QR initially. TPS's LU
  solver is not a substitute for the QR rank/conditioning contract.

## Smallest first implementation

1. Add a narrowly scoped builder/TPS adapter for boundary charts, matrix-to-first-
   order-polynomial conversion and Jacobian extraction. Validate finite entries,
   coordinate metadata and dimensions before entering polynomial operations.
2. In LinearTransferMapBuilder::build, replace only the final ordered matrix product
   with TPS composition of unique segments: first segment first, then outer=next,
   inner=accumulated. Reuse adjacent boundary charts; attach unique occurrence IDs
   scoped to the reference pass, not owner names or pointer allocation order.
3. Extract the combined Jacobian and fill existing result.combined and diagnostics.
   Keep the no-segment result absent, preserve element attachments and finite-
   difference/Richardson settings. No change to ray integration, support selection,
   event location, aperture handling or the existing input-difference solve.
4. Initially keep TPS temporaries inside the initialized host builder. This avoids
   adding persistent TPS ownership/lifetime requirements to ElementBase. Migrating
   stored segment maps to TPS can follow after composition equivalence is proven.

This step neither improves derivative accuracy nor constructs nonlinear maps. Its
acceptance criterion is preserved physical conventions and equivalent first-order
composition with checked chart joins/provenance. Floating-point operation ordering
can differ from the former matrix product, so use justified existing tolerances
rather than promise bitwise equality. No performance improvement is claimed.

Then add a 4D COF adapter and a TPS-facing LinearMapEigenAnalysis overload that
extracts the local Jacobian and delegates to the existing LAPACK implementation.
Route CofCmd through it while retaining matrix/JSON compatibility. Do not recompute
Newton derivatives, change the QR solve, or alter fullTunes/spectral sampling in
that change. New TPS use must occur within IPPL/Kokkos lifetime; existing pure
matrix-only eigenanalysis tests currently do not need TPS initialization.

## Tests required

| Test group | Required assertions |
| --- | --- |
| Adapter coefficients | Every matrix entry survives conversion/extraction; linear evaluation agrees with an independent dense product; degree is one; finite/dimension guards reject bad input. |
| Composition order | Use noncommuting drift/lens and coupled 6x6 examples; compare with independently formed M_N...M_1; preserve empty-result behavior. |
| Chart joins | Translated/rotated frames, boundary momentum/time changes, adjacency identity and mismatched units/scales; preserve nonzero reference slopes at a mismatched ring return. |
| Provenance/ownership | Overlap attached to two owners contributes once; ownerless fringe/tail interval contributes once; repeated physical names remain distinct occurrences. |
| Analytic drift | Retain M01=M23=L and M45=L/gamma0^2; use p0 different from one to distinguish slope and mechanical-momentum conventions. |
| Builder regression | Extend existing drift, overlap, neighbour fringe, ring-return, compensated-clock and Richardson tests; preserve element attachments, enabled/disabled behavior, integration choice and diagnostics. |
| Later COF bridge | Nonzero expansion reference/residual, same-section mechanical metadata, rotated section, fixed-energy pz constraint and independently located arrivals; no extra probe calls. |
| Later eigenanalysis bridge | Direct-matrix/TPS parity for coupled stable rotations, complementary branches, scaled matrices, real/complex instabilities, near resonances and defective bases; preserve report semantics and residual checks. |
| Later tune invariants | Keep integer/conjugate winding, intermediate-coupling rejection, coarse-grid/alias fallback and spectral offset/grid regressions. Do not assert arbitrary spectral/full-tune equality. |
| Runtime/MPI | Construction on host inside active IPPL/Kokkos lifetime; single-rank and dedicated two-rank tests. For future COF wiring, validate rank-zero exception-to-broadcast propagation so metadata errors cannot strand other ranks. |

For matrix conversion/composition use the existing small-case double tolerance
2e-13 where conditioning permits; preserve each existing physical/FD test's own
tolerance. Spectrum comparisons need existing eigenpair tolerances and set/phase-
invariant eigenvector checks, not arbitrary equality of LAPACK eigenvector signs.
No tolerance widening is proposed.

## Audit validation and limits

Reran nine current CTest executables on the existing OpenMP build with
OMP_NUM_THREADS=2: TestOrbitThreader, TestOneTurnMap, TestClosedOrbitSolver,
TestClosedOrbitInitialState, TestLinearMapEigenAnalysis, TestSpectralTunes,
TestCoordinateSystemTrafo, TestTPSIntegration and TestTPSIntegrationMPI. All pass.
Log: /tmp/opalx-tps-map-audit-tests.log. These exercise current behavior, not the
proposed migration, whose new tests are not yet implemented.

The opt-in PSI field-map closed-orbit benchmark is not a distributed default unit
fixture; it needs OPALX_COF_CYCLOTRON_MAP and runs on one rank. No new machine study,
DT/FDSTEP convergence study, GPU validation or performance benchmark was performed.
Host algebra/cache synchronization and cache recovery after allocation failure
remain the separately documented TPS limitations.
