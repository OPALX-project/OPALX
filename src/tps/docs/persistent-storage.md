# Persistent linear TPS maps {#tps_persistent_storage}

The runtime `OpalBeamline` owns one immutable `TransferMapPass`: all ordered,
unique segment handles, including ownerless field-tail intervals, a combined
first-order TPS polynomial and the existing diagnostic values. `OrbitThreader`
also retains its completed pass. Neither the pass nor its segments owns elements.
Each element attachment has a shared segment handle and an independent `pass`
ordinal within that owner's attachment vector. An overlap has one segment object,
even when its ordinal differs between owners. Clones start without attachments.

`LinearTransferMapDraft` holds mutable finite-difference construction scratch.
Completed `LinearTransferMap` records retain TPS as their sole mathematical source
of truth. `matrix()` derives the Jacobian at zero in the existing 6D reporting
type; no mutable matrix cache competes with the polynomial. Original compensated
position, time and path residuals, full reference momentum, transported axes,
integration/refinement settings and contributor metadata remain in
`LinearTransferMapMetadata`. The TPS chart rounds compensated reference values to
double and does not replace that companion data. The differential coordinates
remain `(x,xprime,y,yprime,zeta,delta)` with mechanical slopes, zero affine offset
and order one. The 4D COF momentum chart remains a separate convention.

A segment ID is `reference-pass:N/segment:k`. The generation counter is checked
for exhaustion and serialized per process. Identical construction sequences give
identical IDs on MPI ranks. These IDs distinguish local recalculations; they are
not globally unique across processes/jobs with different construction sequences.
No element name or allocation address is used as traversal identity.

## Publication and invalidation

The builder returns a complete snapshot and transient owner references. It does
not change lattice state. Publication validates owner membership, pass membership
and consecutive owner ordinals, and stages all attachment vectors before no-throw
swaps. Invalid publication leaves the current pass and attachments unchanged.
An empty build has no pass; publishing a null pass with no attachments clears.

Map-enabled design recalculation clears old lattice results before starting.
A failed calculation therefore leaves no current pass; there is no rollback to
old optics. Previously copied immutable snapshots remain usable. Disabled map
calculation and secondary-species passes retain their previous behavior.
`clearTransferMaps()` clears the complete pass and every attachment. Individual
`ElementBase::clearLinearTransferMaps()` clears only that owner's view.
Swapping lattices swaps passes together with elements. Merging invalidates both
lattices' passes and attachments; adding an occurrence invalidates the mapped
geometry. Direct mutation of element field/geometry settings still requires the
caller to clear/recalculate maps; there is no automatic field-change observer.

`getLinearTransferMapsInReferenceOrder()` remains an owner-view query, with shared
handles instead of borrowed vector pointers. It excludes ownerless intervals.
Use `getTransferMapPass()->segments()` for complete future map traversal. Copies
of owner views pin both their representative element and map; complete pass
snapshots pin maps alone and do not keep elements alive.

## Runtime lifetime

TPS HostSpace coefficients require Kokkos. All construction, access, publication,
release and finalization here are serialized host operations. Before IPPL/Kokkos
finalization call `releaseTransferMapStorage()`; `Main.cpp` does so explicitly.
A Kokkos finalize hook provides the same release for tests/embedded callers.
The weak registry does not retain otherwise dead snapshots. Release is idempotent
and frees polynomial payloads in every live persistent handle while the runtime
is active, including external copies. New maps may be constructed after explicit
release while the runtime remains active; the finalize hook still covers them.

After release/finalization, metadata and IDs remain accessible, `available()` is
false and polynomial/Jacobian access throws `OpalException`. Late handle destruction
is safe. Borrowed polynomial references may not cross release/finalization. Copies
of raw `tps::Map` made by a caller carry the usual caller-owned Kokkos lifetime
obligation; the persistent store cannot revoke unrelated copies.

## Scope and verification

Finite-difference rays and Richardson extrapolation retain their existing arithmetic
and tolerances. The input solve and determinant/canonical-form diagnostics now use
TPS/Kokkos; see [6D numerical backend](linear-map-algebra.md) for policies and rounding.
LAPACK eigenanalysis, Newton and full tunes are unchanged. Tests cover complete
nonzero ownerless tails, shared overlap identity and distinct owner ordinals,
invalid publication, snapshot lifetime, cloning, clearing, swapping, merging,
geometry additions, failed recalculation and one/two-rank generation consistency.
A separate process tests explicit/idempotent release, fresh construction after
release, retained segments/passes across finalization and late destruction.

The 6D solve/diagnostic backend migration is also implemented. COF Newton QR
remains an accepted Part 1 backend; see the
[completion audit](part1-completion-audit.md) for the numerical-policy boundaries.
Nonlinear construction and device tracking remain later integration work.
