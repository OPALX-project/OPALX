// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_DEVICE_EXTERNAL_FIELD_H
#define OPALX_DEVICE_EXTERNAL_FIELD_H

#include <Kokkos_NumericTraits.hpp>
#include <limits>
#include <memory>
#include <set>
#include <vector>
#include "AbsBeamline/BendFieldModel.h"
#include "AbsBeamline/CyclotronSector.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/RFCavity.h"
#include "Algorithms/CompensatedSum.h"
#include "Algorithms/TrackReferenceStep.h"
#include "Steppers/BorisPusher.h"
#include "Utilities/ElementList.h"

class OpalBeamline;

/** Internal field-candidate selection and device geometry for ring elements.
 * Analytic geometry and coefficients are uploaded once. Mapped cyclotron
 * candidates are selected separately on the host using their native geometry;
 * this does not extend Builder eligibility. Production ParallelTracker uses
 * analytic descriptors only to select spatial field candidates and test support
 * changes; fields and particle updates use its ordinary element kernels and Boris/PIC
 * sequence. Positions [m], momenta p/(mc), time [s], mass energy [eV], B [T].
 * Analytic external-only transport helpers remain a unit-test oracle. The separate
 * CyclotronLattice capability transports production mapped-sector/thin-gap particles. Aperture-loss
 * ownership remains a separate tracking operation.
 */
namespace device_external {
    using Vector = Vector_t<double, 3>;

    /** Return selected occurrences in the beamline's prepared order.
     * Host-only ordering for spatial analytic-ring and mapped cyclotron callbacks:
     * repeated frame transformations must not depend on allocated element
     * addresses. Selection remains set-valued; distinct same-name occurrences
     * are preserved, while a
     * repeated pointer is applied once. Unknown candidates throw rather than
     * silently omitting a field. Does not change geometry, field formulas,
     * descriptor indices, or the legacy ordering used by other tracking modes.
     */
    std::vector<std::shared_ptr<ElementBase>> orderedCandidates(
            OpalBeamline&, const std::set<std::shared_ptr<ElementBase>>&);

    /** @brief Could a tracking-frame AABB intersect a cyclotron angular wedge?
     *
     * Host-side conservative selection using the current bunch bounds [m]. For
     * sector angle 0<a<=pi and local X=x+Rm, the wedge lies in z>=0 and
     * sin(a)*X-cos(a)*z>=0. Reject only when every transformed corner violates a
     * half-plane. The 64 epsilon scale margin broadens candidates only; native
     * device field kernels retain radial/vertical bounds and half-open seams.
     * Invalid bounds/geometry conservatively retain the candidate.
     */
    bool cyclotronSectorMayIntersect(
            const Geometry&, const CoordinateSystemTrafo& toLocal, const Vector& minimum,
            const Vector& maximum);

    /** @brief Replace mapped cyclotron membership using current spatial bounds.
     *
     * A periodic path index describes the threader's first orbit, not every
     * TRACK midpoint or particle orbit. Preserve all non-cyclotron candidates;
     * select supplied sectors using the same misaligned transform as ordinary
     * element callbacks. Bounds are already on the host, so selection requires
     * no particle copies or additional MPI collectives. This host capability is
     * independent of analytic Builder::supports and introduces no descriptors.
     */
    std::set<std::shared_ptr<ElementBase>> selectCyclotronCandidates(
            const ElementList& sectors, const CoordinateSystemTrafo& toLab, const Vector& minimum,
            const Vector& maximum, std::set<std::shared_ptr<ElementBase>> mapped);

    /// Device value representation of CoordinateSystemTrafo (parent -> local).
    struct Rigid {
        Vector origin = Vector(0);
        matrix3x3_t rotation;
        KOKKOS_INLINE_FUNCTION Vector pointTo(const Vector& r) const {
            return prod_vector(rotation, Vector(r - origin));
        }
        KOKKOS_INLINE_FUNCTION Vector pointFrom(const Vector& r) const {
            return Vector(prod_vector_transpose(rotation, r) + origin);
        }
        KOKKOS_INLINE_FUNCTION Vector vectorTo(const Vector& p) const {
            return prod_vector(rotation, p);
        }
        KOKKOS_INLINE_FUNCTION Vector vectorFrom(const Vector& p) const {
            return prod_vector_transpose(rotation, p);
        }
    };

    /** Explicit mapped-cyclotron capability, separate from analytic Builder.
     * Descriptors borrow immutable device map storage and perform no host calls.
     * Positions [m], momenta p/(mc), mass rest energy [eV]; laboratory SI fields.
     * Native sector kernels retain exact radial/vertical/angular support.
     * Magnetic/self-field kicks assume charge +1. Production tracking restricts
     * this capability to unpolarized positive protons, without restart.
     */
    struct MappedCyclotronSector {
        Rigid frame;
        CyclotronSector::DeviceField field;
    };
    struct CyclotronGap {
        Rigid frame;
        RFCavity::DeviceGap field;
    };
    struct CyclotronLattice {
        Kokkos::View<MappedCyclotronSector*> sectors;
        Kokkos::View<CyclotronGap*> gaps;
        /// Drift-kick-drift in the spatial magnetic map; no self fields or RF.
        KOKKOS_INLINE_FUNCTION bool magnetic(Vector& r, Vector& p, double dt, double mass) const {
            Vector x = r, v = p;
            x += (0.5 * dt * Physics::c / Kokkos::sqrt(1 + dot(v, v))) * v;
            Vector b(0);
            bool supported = false;
            for (size_t j = 0; j < sectors.extent(0); ++j) {
                const auto& sector = sectors(j);
                Vector localB(0);
                supported = sector.field.add(sector.frame.pointTo(x), localB) || supported;
                b += sector.frame.vectorFrom(localB);
            }
            if (!supported) return false;
            BorisPusher().kick(x, v, Vector(0), b, dt, mass, 1);
            x += (0.5 * dt * Physics::c / Kokkos::sqrt(1 + dot(v, v))) * v;
            r = x;
            p = v;
            return true;
        }

        /// Device policy: one visit per gap, no per-particle energy-target stop.
        struct EventPolicy {
            const CyclotronLattice& lattice;
            unsigned long long visitedMask = 0;
            KOKKOS_INLINE_FUNCTION size_t count() const { return lattice.gaps.extent(0); }
            KOKKOS_INLINE_FUNCTION bool magnetic(
                    Vector& r, Vector& p, double dt, double, double mass) const {
                return lattice.magnetic(r, p, dt, mass);
            }
            KOKKOS_INLINE_FUNCTION double distance(size_t j, const Vector& r) const {
                return lattice.gaps(j).frame.pointTo(r)(2);
            }
            KOKKOS_INLINE_FUNCTION bool supports(size_t j, const Vector& r) const {
                const auto& gap = lattice.gaps(j);
                return gap.field.supports(gap.frame.pointTo(r)(0));
            }
            KOKKOS_INLINE_FUNCTION bool kick(
                    size_t j, const Vector& r, Vector& p, double t, double mass) const {
                const auto& gap = lattice.gaps(j);
                auto localP     = gap.frame.vectorTo(p);
                if (!gap.field.apply(gap.frame.pointTo(r)(0), t, mass, localP)) return false;
                p = gap.frame.vectorFrom(localP);
                return true;
            }
            KOKKOS_INLINE_FUNCTION bool visited(size_t j) const {
                return visitedMask & (1ULL << j);
            }
            KOKKOS_INLINE_FUNCTION void markVisited(size_t j) { visitedMask |= 1ULL << j; }
            KOKKOS_INLINE_FUNCTION bool accepted(
                    size_t, const Vector&, const Vector&, double) const {
                return false;
            }
        };

        /** Magnetic transport with independently timed thin RF kicks.
         * The shared event search operates on copies, so a failed particle leaves
         * its live state intact. The caller owns bunch-level failure handling.
         * No allocations/collectives; at most 64 gap descriptors on this path.
         * terminalEndpoint includes a roundoff-short final gap crossing only
         * after the reference has completed its EKINSTOP kick. frameSteps is the
         * accepted global step count including this final step; rf_event::advance
         * uses it only for accumulated frame-arithmetic endpoint resolution.
         */
        KOKKOS_INLINE_FUNCTION bool advance(
                Vector& position, Vector& momentum, double time, double dt, double mass,
                rf_event::Diagnostic* diagnostic = nullptr, bool terminalEndpoint = false,
                unsigned long long frameSteps = 1) const {
            if (gaps.extent(0) > 64 || !(dt >= 0)) {
                if (diagnostic)
                    diagnostic->set(rf_event::Failure::InvalidStep, position, momentum, time);
                return false;
            }
            Vector r = position, p = momentum;
            EventPolicy policy{*this};
            if (!rf_event::advance(
                         r, p, time, dt, mass, policy, diagnostic, terminalEndpoint, frameSteps)
                         .success)
                return false;
            position = r;
            momentum = p;
            return true;
        }

        /** Frozen gathered self fields [V/m,T], split into Boris half kicks.
         * PIC prediction, solve/gather, and accepted duration belong to the tracker.
         * This thin-gap approximation does not solve self fields inside RF events.
         */
        KOKKOS_INLINE_FUNCTION bool advanceWithSelfFields(
                Vector& r, Vector& p, const Vector& e, const Vector& b, double t, double dt,
                double mass, rf_event::Diagnostic* diagnostic = nullptr,
                bool terminalEndpoint = false, unsigned long long frameSteps = 1) const {
            Vector x = r, v = p;
            BorisPusher().kick(x, v, e, b, .5 * dt, mass, 1);
            if (!advance(x, v, t, dt, mass, diagnostic, terminalEndpoint, frameSteps)) return false;
            BorisPusher().kick(x, v, e, b, .5 * dt, mass, 1);
            r = x;
            p = v;
            return true;
        }
    };

    /// No pointers to host elements or virtual calls in the device descriptor.
    struct Element {
        enum Kind { Passive, Multipole, SectorBend, RectangularBend, UniformRF };
        Kind kind = Passive;
        Rigid frame;
        double begin = 0, end = 0;
        ApertureType aperture = ApertureType::RECTANGULAR;
        double apertureX = 1e6, apertureY = 1e6;
        BendFieldModel::FieldInputs bend{};

        KOKKOS_INLINE_FUNCTION Vector chart(const Vector& lab) const {
            const auto local = frame.pointTo(lab);
            return kind == SectorBend
                           ? GeometryHelper::toBendArcCoords(local, bend.curvature, bend.bodyLength)
                           : local;
        }
        KOKKOS_INLINE_FUNCTION bool contains(const Vector& lab) const {
            const auto r = chart(lab);
            // VariableRFCavity's WIDTH/HEIGHT include their transverse faces.
            // This is field support; generic aperture-loss ownership is separate.
            if (kind == UniformRF)
                return r(2) >= begin && r(2) < end && Kokkos::abs(r(0)) <= apertureX
                       && Kokkos::abs(r(1)) <= apertureY;
            return r(2) >= begin && r(2) < end
                   && ApertureHelper::isInsideAperture(r(0), r(1), aperture, apertureX, apertureY);
        }
        /// Caller has checked support. Reuses the normal bend device field model.
        KOKKOS_INLINE_FUNCTION Vector magnetic(const Vector& lab) const {
            const auto r = chart(lab);
            Vector b(0);
            if (kind == Multipole) {
                b(0) = -bend.dipoleSkew + bend.quadNormal * r(1) - bend.quadSkew * r(0);
                b(1) = bend.dipoleNormal + bend.quadNormal * r(0) + bend.quadSkew * r(1);
            } else if (kind == SectorBend || kind == RectangularBend) {
                b = BendFieldModel::bendField(r, bend);
                if (kind == SectorBend)
                    b = GeometryHelper::rotateArcFieldToEntry(
                            b, r(2), bend.curvature, bend.bodyLength);
            }
            return frame.vectorFrom(b);
        }
    };

    struct State {
        Vector position = Vector(0), momentum = Vector(0), correction = Vector(0);
    };

    struct Lattice {
        Kokkos::View<const Element*> elements;
        double maximumStep = std::numeric_limits<double>::max();
        bool magneticOnly  = true;  ///< The external-only test oracle cannot transport RF.

        KOKKOS_INLINE_FUNCTION Vector magnetic(const Vector& r) const {
            Vector b(0);
            for (size_t i = 0; i < elements.extent(0); ++i)
                if (elements(i).contains(r)) b += elements(i).magnetic(r);
            return b;
        }

        KOKKOS_INLINE_FUNCTION bool differentSupport(const Vector& a, const Vector& b) const {
            for (size_t i = 0; i < elements.extent(0); ++i)
                if (elements(i).contains(a) != elements(i).contains(b)) return true;
            return false;
        }

        KOKKOS_INLINE_FUNCTION static void halfDrift(State& ray, double h) {
            const double factor =
                    0.5 * Physics::c * h / Kokkos::sqrt(1 + dot(ray.momentum, ray.momentum));
            for (unsigned d = 0; d < 3; ++d)
                compensated::add(factor * ray.momentum(d), ray.position(d), ray.correction(d));
        }

        /** External-only test oracle: depth-first subdivision with a bounded local stack.
         * No recursion, allocations, particle copies to host, or synchronization
         * between particles. False signals invalid input or exhausted refinement;
         * the caller must discard the partial state and report failure.
         */
        KOKKOS_INLINE_FUNCTION bool advance(
                State& state, double dt, double mass, double charge) const {
            if (!magneticOnly || !Kokkos::isfinite(dt) || dt == 0 || !Kokkos::isfinite(mass)
                || mass <= 0 || !Kokkos::isfinite(charge))
                return false;
            unsigned char pending[65];
            unsigned count = 1;
            pending[0]     = 0;
            while (count) {
                const unsigned depth = pending[--count];
                const double h       = Kokkos::ldexp(dt, -static_cast<int>(depth));
                const double scale =
                        Kokkos::fmax(1., Kokkos::sqrt(dot(state.position, state.position)));
                const double floor = 64 * Kokkos::epsilon_v<double> * scale / Physics::c;
                const bool resolved =
                        Kokkos::abs(h) <= Kokkos::fmax(1e-12 * Kokkos::abs(dt), floor);
                bool split  = !resolved && Kokkos::abs(h) > maximumStep;
                State trial = state;
                if (!split) {
                    halfDrift(trial, h);
                    const Vector midpoint = trial.position;
                    const Vector b        = magnetic(midpoint);
                    BorisPusher().kick(midpoint, trial.momentum, Vector(0), b, h, mass, charge);
                    halfDrift(trial, h);
                    split = !resolved
                            && (differentSupport(state.position, midpoint)
                                || differentSupport(state.position, trial.position));
                }
                if (split) {
                    if (depth >= 64 || count + 2 > 65) return false;
                    pending[count++] = depth + 1;
                    pending[count++] = depth + 1;
                } else {
                    for (unsigned d = 0; d < 3; ++d)
                        if (!Kokkos::isfinite(trial.position(d))
                            || !Kokkos::isfinite(trial.momentum(d)))
                            return false;
                    state = trial;
                }
            }
            return true;
        }

        /** Launch in the particle view's execution space. Only the scalar failure
         * count is reduced to the host. E/B diagnostics are endpoint fields in the
         * current bunch frame (E=0 in this static magnetic transport).
         */
        template <class View>
        int transport(
                View r, View p, View e, View b, size_t n, Rigid toLab, double dt, double mass,
                double charge) const {
            const Lattice lattice = *this;
            int failures          = 0;
            using Execution       = typename View::execution_space;
            Kokkos::parallel_reduce(
                    "Ring::deviceBoris", Kokkos::RangePolicy<Execution>(0, n),
                    KOKKOS_LAMBDA(size_t i, int& errors) {
                        State ray;
                        ray.position = toLab.pointTo(r(i));
                        ray.momentum = toLab.vectorTo(p(i));
                        if (!lattice.advance(ray, dt, mass, charge)) {
                            ++errors;
                            return;
                        }
                        r(i) = toLab.pointFrom(ray.position);
                        p(i) = toLab.vectorFrom(ray.momentum);
                        e(i) = Vector(0);
                        b(i) = toLab.vectorFrom(lattice.magnetic(ray.position));
                    },
                    failures);
            return failures;
        }
    };

    /// Host-only snapshot builder; never receives or reads particle arrays.
    struct Builder {
        /// Inspect eligibility without building descriptors or changing the lattice.
        /// Uses the existing non-const OpalBeamline::getElements() read accessor.
        static bool supports(OpalBeamline&);
        static Lattice build(OpalBeamline&);

    private:
        static const char* unsupportedReason(const ElementBase&);
    };
}  // namespace device_external
#endif
