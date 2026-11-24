# Critical code review {#code_review}

Reviewed in the standalone test bed on 2026-11-24 for the OPALX construction, composition, analysis and tracking
workflow. Scope: the standalone working library, build configuration, tests and
examples. The historical baseline is excluded from cleanup. This is a source and
regression review, not a proof of numerical correctness or GPU performance.

## Findings addressed

- **Readability:** compressed control flow obscured validation and ownership in
  the new APIs and tests. Apply a checked-in formatting policy throughout working
  C++ sources. PreparedMap now separates `buildPlan`, `ensureNode` and `upload`;
  comments explain parent-first ordering and detached host mirrors. Preserve the
  numerical operations and evaluation order.
- **Documentation:** legacy `//` continuations were omitted from Doxygen member
  documentation, and bracket markup was displayed literally. Convert declaration
  continuations to Doxygen comments and Markdown; correct stale unchecked-index
  claims in Vps and old RangeError names. Document borrowed array View lifetimes,
  the prepared recurrence, error reporting and scalar identity shifts in Matrix.
  Retain copyright notices; remove RCS metadata duplicated by the baseline.
- **Header dependencies:** TpsMath included only the scalar declarations while
  defining algorithms using array/series implementations. Include Tps.hpp directly;
  HostStorage now directly includes its standard exception/string dependencies.
  Compile every installed header in its own translation unit and instantiate
  elementary functions without the umbrella header.
- **Redundant code:** nine legacy exception header/implementation pairs repeated
  trivial copy/destructor definitions. Replace them with Error.hpp; retain the
  seven meaningful derived error types and their messages. Remove the unused
  ArithmeticError intermediate base; rename ClassicException to Error. No working
  throw/catch site needed the intermediate category. Baseline exception files and
  OPALX's separate global exception family remain unchanged.
- **Redundant special members:** TpsSubstitution is a documented aggregate with
  zero-initialized fields; remove hand-written constructor/copy/destructor/assignment.
  Remove unused scalar-function using declarations and a redundant empty statement.
  Default forwarding special members in TpsMonomial and empty array destructors.

## Retained deliberately

Array1D and Array2D supply checked value semantics and manage nontrivial coefficient
objects; they are not unused wrappers. SliceIterator/ConstSliceIterator implement
live column traversal and have C++20 iterator/conformance regressions. Matrix and
Vector delegate numerical operations to Kokkos Kernels and remain dependencies of
substitution, Jacobians and linear solves. Vps provides VpsMap component storage.
TpsMonomial, TpsData and TpsSubstitution implement coefficient/cache bookkeeping.
Complex/nested scalar support, serialization, differentiation and integration are
retained: internal call counts alone do not establish that public algebra is dead.
The hand-written LU implementation had already been removed before this review.

## Remaining concerns, in priority order

1. **High — concurrent host algebra is unsafe.** ScopedMapOrder still temporarily
   changes per-scalar global truncation, and TpsData::getTpsData/build mutates a
   shared cache. TpsRep reference counts are not atomic. OPALX must serialize host
   construction/composition/preparation, or replace this state with explicit owned
   algebra contexts before allowing concurrent host map algebra. Immutable prepared
   evaluation does not access these globals; each concurrent evaluation still
   needs disjoint workspace. An evaluation-only parallel test cannot validate host
   construction concurrency.
2. **Medium — cache rebuild failure is not transactional.** TpsData::build sets
   topOrder before allocating/filling all replacement tables. An allocation or
   coefficient-operation exception can leave partial state, and a later request
   for that order can skip rebuilding. Do not assume catching allocation failure
   makes the host cache reusable. Follow-up: build complete replacement tables
   before commit and test injected allocation failures. This pre-existing numerical
   ownership change is not hidden inside a readability cleanup.
3. **Medium — GPU compiler/runtime coverage is outstanding.** Device annotations
   and CPU Kokkos tests do not establish CUDA/HIP compilation, cross-stream lifetime
   correctness or device throughput. Run the same prepared-map tests on OPALX's
   selected device backend before enabling GPU tracking; benchmark realistic orders
   and sparsities before changing the dense output-by-node coefficient layout.
4. **Medium — integration contracts remain adapter responsibilities.** DeviceSpan
   cannot verify real allocation bounds, lifetime, memory accessibility or overlap.
   Chart matching is exact metadata equality; no physical conversion is implicit.
   OPALX must supply coordinate/reference conversion, scheduling/workspace,
   aperture/loss handling and space-charge splitting. Map accuracy outside a local
   expansion domain and symplecticity are not certified by metadata.
5. **Low — documentation coverage is not exhaustive.** Warning-free Doxygen with
   EXTRACT_ALL is not proof every legacy algebra operation has complete contracts.
   The new boundary APIs and math conventions are the preferred entry point;
   coefficient-type restrictions and compound-operation exception guarantees need
   further targeted review as OPALX analysis use cases are added.

6. **Low — GCC iterator diagnostic remains visible.** GCC 15 at Release optimization
   warns about a hypothetical ptrdiff_t position overflow in the inlined
   `std::ranges::sort` column-iterator regression. Iterator increment/equality
   expressions are unchanged by this review; the tested range has four elements.
   Bounded column traversal/sorting and sanitizer tests pass. Arbitrary positions
   beyond representable offsets remain outside the iterator contract. No warning
   suppression or numerical tolerance change was added.

## Verification

See HANDOFF-TPS.md at the OPALX repository root for integration results. The
standalone review and its original results remain in the external `opal-tps` test bed. No numerical
tolerance changes are part of this cleanup. MPI behavior is not applicable to this
standalone library, which has no MPI dependency; OPALX rank integration is separate.

## OPALX import

The reviewed working sources are now under src/tps and regression executables
under unit_tests/tps. Only the standalone baseline alternatives in test sources
were removed; the historical baseline remains in the external test bed. OPALX
provides the selected Kokkos runtime; TPS reuses it and resolves Kokkos Kernels.
Single- and two-rank integration tests exercise OPALX linkage and IPPL lifecycle.
The host concurrency and cache-recovery restrictions above remain unchanged.
