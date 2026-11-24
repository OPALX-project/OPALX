# OPALX map workflow audit {#workflow_audit}

Scope: construct element/field-map transfer maps, compose them, normalize and
analyze them, and use either element maps or a composed map for tracking.
This is a source audit of the local TPS and OPALX checkouts on 2026-11-24.
Proposed function names below are design suggestions, not implemented APIs.

## What is present

| Capability | Current support | Qualification |
|---|---|---|
| Multivariate polynomial algebra | Tps, elementary functions, derivatives, integration | Dense, total-degree truncation; finite precision policy is documented |
| Nonlinear map composition | VpsMap::substitute(inner, order) | Outer acts after inner; not ordinary matrix multiplication |
| Single-point tracking | constantTerm(point), despite its misleading name | Generic scalar evaluation, with per-call scratch allocation |
| Jacobian | linearTerms(point), linearTerms() | Analytic derivative of stored polynomial; constant-component case now covered by regressions |
| Linear coordinate transforms | matrix substitution and matrix-on-map multiplication | Rectangular output transform repaired and regression-tested |
| Linear solve/inverse | Kokkos-backed LUFactorization | Pivoted Getrf/Getrs; no nonlinear normal-form machinery |
| Complex/nested coefficient arithmetic | Header-instantiated templates | Selected regression coverage, not every elementary function supports complex arguments |
| Eigenanalysis in OPALX | Algorithms/LinearMapEigenAnalysis.cpp | Already uses LAPACKE_dgeev and dgesvd; not part of TPS |

## Cleanup applied

1. Removed working TpsDouble.cpp, TpsComplex.cpp and TpsTps.cpp. They only forced
   explicit instantiation; there are no extern-template consumers. Templates
   remain in headers and all coefficient types remain available. The historical
   build retains its original translation units. Complex/nested regressions pass.
2. Removed four pre-standard iterator_category/value_type free helpers. The
   iterators supply standard traits and satisfy C++20 iterator concepts.
3. Removed unused LUMatrix includes from Vps.hpp and VpsMap.h, plus VpsMap's
   self-include. LU is retained as an explicit standalone include, because linear
   solves are relevant to closed-orbit and normalization work.
4. Removed the declaration-only cot function (no implementation) and the broken
   VpsMap(array) constructor, which called a nonexistent Vps constructor. They
   provided no working capability. A safe span-based component constructor can
   replace the latter when the public map API is consolidated.

Not deleted merely because current tests do not call it: complex coefficients,
scalar integration, matrix solves, filtering, polynomial I/O and nested TPS.
These support field integration, normal forms, persistence or parameter maps.
There are no external-consumer usage inventories, so this is not a claim that
all public features have known callers. Reachability alone cannot prove a
public template API is unnecessary.

## Missing interfaces, in implementation order

| Priority | Needed capability / suggested API | Ownership |
|---|---|---|
| P0 (core implemented) | MapContext now records coordinates, order, references, frames, units, scaling and canonical pairs; OPALX adapters remain | TPS core context + OPALX adapter |
| P0 (implemented) | Context-bearing Map, identityMap, evaluate, jacobian, compose, composeSequence with explicit order and checked boundaries | TPS public map layer |
| P0 | buildElementMap / buildSegmentMap along an OPALX reference trajectory, including field overlaps | OPALX adapter/integrator |
| P0 | TPS-valued field evaluation or local derivative/jet provider; cell/support boundary policy | OPALX field-map adapter |
| P0 | Map-vs-direct-ray validation, order/step convergence and a certified/tested validity region | Joint validation harness |
| P1 | Closed-orbit/fixed-point solve; recenter and reexpand about an orbit; input/output chart conversion | Analysis/adapter layer |
| P1 | Linear spectrum API returning eigenvectors, paired modes, tunes, residuals and conditioning | Reuse OPALX LAPACK backend |
| P1 | Canonical symplectic residuals and mode normalization; coordinate-scale/linear normalization | Analysis layer |
| P1 | Formal inverseMap (nonsingular linear part), coordinate conjugation A^-1 F A | TPS map/analysis layer |
| P1 | evaluateBatch with reusable workspace, no allocation per particle, loss/aperture/validity status | Tracking adapter + TPS evaluator |
| P1 | Versioned map serialization with coefficients AND context/reference metadata | Map layer |
| P2 | Poisson bracket, Lie action/exponential, homogeneous-order factorization, nonlinear normal forms | Nonlinear analysis layer |
| P2 | Resonance handling, tune shifts/detuning and chromatic coefficients | Nonlinear analysis + physics conventions |

A thread-safe context/workspace is needed before concurrent map construction or
parallel host tracking. The existing global truncation and mutable fixed cache
are not suitable as the final shared mutable execution model.

## Construction from elements and field maps

Current OPALX interfaces inspected:

- src/AbsBeamline/ElementBase.h: apply() uses concrete particle/vector types.
- src/Fields/Fieldmap.h: getFieldstrength and getFieldDerivative accept double
  positions and expose values/first spatial derivatives, not an arbitrary-order
  TPS field jet. A virtual double-only interface cannot propagate DA coefficients
  merely by replacing a particle scalar type.
- src/Structure/LinearTransferMap.h: current linear maps include reference frames,
  segment IDs/ownership and overlap information. Their reported slopes/mechanical
  momenta are not globally canonical. Existing ray construction excludes RF and
  collective fields; those assumptions must not silently become nonlinear-map
  guarantees.

Two construction routes need an explicit choice: propagate TPS through a
scalar-generic orbit/field integrator, or build local field jets and integrate
those. Sampling/fitting direct rays is another possible route but requires
conditioning and truncation-error studies; it is not automatically differential
algebra. Finite differences of field values alone do not supply reliable
arbitrary-order derivatives.

For tabulated fields, identify the interpolation cell from the reference point,
represent its interpolation polynomial consistently, and detect/support domain
crossings. A single TPS does not represent a branch across a nonsmooth cell
boundary. High-order maps need adequate interpolation smoothness, controlled
integration steps, support/aperture events and reference-plane termination.
Time-dependent/RF elements require a consistent time/energy coordinate and phase
convention. Collective effects require a stated frozen-field or other model;
a fixed single-particle map does not reproduce general self-consistent dynamics.

Maps must join matching output/input charts. Compose disjoint path-ordered
segments once; overlapping fields must be summed before integration, rather
than multiplying independent overlapping element maps or counting shared
segments multiple times.

## Normalization, decomposition and EV analysis

Separate three operations:

1. Recenter and convert physical coordinates to a defined canonical chart, with
   explicit reference momentum/time and scaling. Slopes/mechanical momenta alone
   do not justify a canonical J matrix.
2. Linear normal-mode analysis of the Jacobian at the fixed point. Reuse
   OPALX's LAPACK-backed LinearMapEigenAnalysis; verify eigen-residuals, mode
   pairing, stability and conditioning before symplectic normalization.
3. Nonlinear factorization/normal forms, which additionally need Poisson/Lie
   algebra, map inversion/conjugation and explicit resonance tolerances.

Do not use the older GSLEigen.h compatibility heuristic as a replacement for the
existing LAPACK analysis. Its own comments describe approximate eigenvectors.
The TPS library should expose matrix adapters/results, not import OPALX as a
mandatory dependency or fetch a duplicate LAPACK build.

## Correctness findings resolved after this audit

The rectangular output-transform loop now iterates over matrix output rows.
Tests cover reductions, square maps, expansions and incompatible dimensions.
Constant components now give zero Jacobian rows in both origin and point-based
calls. For an all-scalar legacy VpsMap the input dimension is not retained: use
the point-taking Jacobian or the new context-bearing Map to supply it.

The context API is now available in MapContext.hpp and Map.hpp: validated chart
metadata, immutable coefficient snapshots, explicit order, checked composition,
local evaluation and Jacobians. Remaining OPALX adapter and analysis gaps in the
table above are unchanged. Global scalar precision is scoped/restored, not made
thread-safe; the shared cache still requires external synchronization.

The existing regression suite does not establish correctness of arbitrary
six-dimensional maps, high orders, nonzero-origin composition, map inversion,
normal forms or field-map construction. Those need targeted analytic and OPALX
reference-ray tests as each capability is introduced.

## Tracking and benchmarking

Scalar evaluation already exists, so basic tracking is not missing mathematically.
What is missing is a stable tracking interface: batch/SoA input/output, reusable
monomial workspace, validity/loss handling, and preservation of element order
for aperture checks and diagnostics. Combining an entire lattice into one map
can lose intermediate aperture information. A truncated Taylor map is not
necessarily symplectic under indefinite iteration.

Benchmark both element-by-element and composed tracking; compare map construction
cost, composition cost, memory, throughput and accuracy against direct OPALX
tracking. Start with six variables and multiple orders, and separate field-map
interpolation error from map truncation and integration error. No benchmark or
speedup claim is made by this audit.

The former gap diagnostic was replaced by `map_defects` regressions. Context
and composition behavior is checked by `map_context` and `context_example`.


Subsequent Kokkos migration: dynamic storage now uses HostSpace Views and fixed
frame data uses Kokkos::Array. Hand-written Matrix/Vector kernels and LUMatrix
were replaced by checked Kokkos Kernels adapters and LUFactorization. Existing
OPALX/IPPL targets are reused; standalone discovery/fetch mirrors IPPL's 5.2
configuration. A complete device batch evaluator remains a separate milestone.
