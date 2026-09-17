// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPALX_DEVICE_EXTERNAL_FIELD_H
#define OPALX_DEVICE_EXTERNAL_FIELD_H

#include <limits>
#include <memory>
#include <set>
#include <vector>
#include "AbsBeamline/BendFieldModel.h"
#include "AbsBeamline/ElementBase.h"
#include "Algorithms/CompensatedSum.h"
#include "Steppers/BorisPusher.h"

class OpalBeamline;

/** Internal device geometry for analytic ring elements. Geometry and
 * coefficients are uploaded once. Production ParallelTracker uses these
 * descriptors only to select spatial field candidates and test support changes;
 * fields and particle updates use its ordinary element kernels and Boris/PIC
 * sequence. Positions [m], momenta p/(mc), time [s], mass energy [eV], B [T].
 * The independent external-only transport helpers remain a unit-test oracle;
 * they are not called by production particle tracking. Aperture-loss ownership
 * remains a separate tracking operation.
 */
namespace device_external {
    using Vector = Vector_t<double, 3>;

    /** Return selected occurrences in the beamline's prepared order.
     * Host-only ordering for the bare analytic-ring callbacks: repeated frame
     * transformations must not depend on allocated element addresses. Selection
     * remains set-valued; distinct same-name occurrences are preserved, while a
     * repeated pointer is applied once. Unknown candidates throw rather than
     * silently omitting a field. Does not change geometry, field formulas,
     * descriptor indices, or the legacy ordering used by other tracking modes.
     */
    std::vector<std::shared_ptr<ElementBase>> orderedCandidates(
            OpalBeamline&, const std::set<std::shared_ptr<ElementBase>>&);

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
                const double floor =
                        64 * std::numeric_limits<double>::epsilon() * scale / Physics::c;
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
