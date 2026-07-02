#ifndef OPALX_BendFieldModel_HH
#define OPALX_BendFieldModel_HH

//
// BendFieldModel
//   Stateless field-shape math shared by the analytic SBEND and RBEND fringe
//   fields: the OPAL default Enge longitudinal profile (amplitude and its first
//   two derivatives), the combined amplitude from the two pole-face distances,
//   the pole-face vertical edge-focusing coefficient, and the field evaluation
//   itself. All functions are device-callable and hold no state; the bends supply
//   their own geometry scalars (gap, curvature, face angles) and coefficients.
//   Coordinate conversions live in GeometryHelper (BeamlineGeometry/Geometry.h).
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "VectorMath.h"

#include <Kokkos_Core.hpp>

/**
 * @namespace BendFieldModel
 * @brief Stateless field-shape math for the analytic bend fringe model.
 *
 * The longitudinal field envelope is the OPAL default `FM1DPROFILE1` Enge
 * function
 * @f[
 *   F(u) = \frac{1}{1 + \exp\!\left(\sum_{i=0}^{5} c_i u^i\right)}, \qquad
 *   u = \frac{d}{g},
 * @f]
 * where @f$d@f$ is the (perpendicular) distance from a pole face and @f$g@f$ is
 * the full gap. @f$F \to 1@f$ deep inside the magnet, @f$F \to 0@f$ outside, and
 * @f$F \approx 0.3825@f$ at the nominal face.
 */
namespace BendFieldModel {

    /// Enge amplitude together with its first and second derivatives with respect
    /// to the distance coordinate passed to engeProfile().
    struct EngeValue {
        double value;             ///< F
        double firstDerivative;   ///< dF/dd
        double secondDerivative;  ///< d^2F/dd^2
    };

    /**
     * @brief OPAL default Enge longitudinal profile and its derivatives.
     *
     * @param coordinate Distance from the pole face in metres (negative inside the
     *        body, zero at the face, positive outside).
     * @param profileGap Full gap @f$g@f$ scaling the profile. A non-positive gap
     *        gives a hard edge (F = 1 everywhere).
     * @return F and dF/dd, d^2F/dd^2 with respect to @p coordinate.
     */
    KOKKOS_INLINE_FUNCTION EngeValue engeProfile(const double coordinate, const double profileGap) {
        if (profileGap <= 0.0) {
            return EngeValue{1.0, 0.0, 0.0};
        }

        const double c0 = 0.478959;
        const double c1 = 1.911289;
        const double c2 = -1.185953;
        const double c3 = 1.630554;
        const double c4 = -1.082657;
        const double c5 = 0.318111;

        const double u  = coordinate / profileGap;
        const double u2 = u * u;
        const double u3 = u2 * u;
        const double u4 = u2 * u2;
        const double u5 = u4 * u;

        const double exponent = c0 + c1 * u + c2 * u2 + c3 * u3 + c4 * u4 + c5 * u5;
        if (exponent > 80.0) {
            return EngeValue{0.0, 0.0, 0.0};
        }
        if (exponent < -80.0) {
            return EngeValue{1.0, 0.0, 0.0};
        }

        const double invGap    = 1.0 / profileGap;
        const double polyPrime  = c1 + 2.0 * c2 * u + 3.0 * c3 * u2 + 4.0 * c4 * u3 + 5.0 * c5 * u4;
        const double polySecond = 2.0 * c2 + 6.0 * c3 * u + 12.0 * c4 * u2 + 20.0 * c5 * u3;
        const double expPrime   = polyPrime * invGap;
        const double expSecond  = polySecond * invGap * invGap;

        const double engeExp = Kokkos::exp(exponent);
        const double f       = 1.0 / (1.0 + engeExp);
        const double f2      = f * f;
        const double dExp    = expPrime * engeExp;
        const double d2Exp   = (expSecond + expPrime * expPrime) * engeExp;

        const double firstDerivative  = -dExp * f2;
        const double secondDerivative = -d2Exp * f2 + 2.0 * dExp * dExp * f2 * f;
        return EngeValue{f, firstDerivative, secondDerivative};
    }

    /// Full gap used to scale the Enge profile: GAP if set, otherwise 2*HGAP.
    KOKKOS_INLINE_FUNCTION double profileGap(const double fullGap, const double halfGap) {
        return (fullGap > 0.0) ? fullGap : 2.0 * halfGap;
    }

    /// Half width of one fringe. OPAL's default map places its profile points at
    /// 0.1 m for a 0.02 m gap, i.e. five gaps.
    KOKKOS_INLINE_FUNCTION double fringeHalfWidth(const double profileGap) {
        return 5.0 * profileGap;
    }

    /// Combined fringe amplitude, with the derivatives of the limiting face.
    struct FringeAmplitude {
        double value;             ///< F in [0, 1]
        double firstDerivative;   ///< dF/dd of the limiting face
        double secondDerivative;  ///< d^2F/dd^2 of the limiting face
        int activeFace;           ///< 0 = entrance limits, 1 = exit limits
    };

    /**
     * @brief Fringe amplitude from the two pole-face distances.
     *
     * @param distPastEntrance Distance past the entrance face (negative inside).
     * @param distPastExit Distance past the exit face (negative inside the body).
     * @param profileGap Full gap.
     * @return The smaller of the entrance and exit Enge envelopes (1 deep inside
     *         the body, decaying to 0 across either fringe) plus which face limits
     *         it, so the caller can apply that face's derivative correction.
     */
    KOKKOS_INLINE_FUNCTION FringeAmplitude fringeAmplitude(
            const double distPastEntrance, const double distPastExit, const double profileGap) {
        const EngeValue entrance = engeProfile(distPastEntrance, profileGap);
        const EngeValue exit     = engeProfile(distPastExit, profileGap);
        if (exit.value < entrance.value) {
            return FringeAmplitude{exit.value, exit.firstDerivative, exit.secondDerivative, 1};
        }
        return FringeAmplitude{entrance.value, entrance.firstDerivative, entrance.secondDerivative, 0};
    }

    /**
     * @brief Vertical edge-focusing kick coefficient with the FINT correction.
     *
     * @f[
     *   k_y = -h \tan(E - \psi), \qquad
     *   \psi = h\,\mathrm{HGAP}\,\mathrm{FINT}\,\frac{1 + \sin^2 E}{\cos E},
     * @f]
     * with signed curvature @f$h@f$. The caller distributes it over the fringe as
     * @f$B_x \mathrel{+}= (B_0/h)\,k_y\,F'\,y@f$.
     */
    KOKKOS_INLINE_FUNCTION double edgeVerticalKickCoefficient(
            const double curvature, const double halfGap, const double fringeIntegral,
            const double edgeAngle) {
        const double cosE    = Kokkos::cos(edgeAngle);
        const double safeCos = (Kokkos::abs(cosE) > 1.0e-6) ? cosE : Kokkos::copysign(1.0e-6, cosE);
        const double sinE    = Kokkos::sin(edgeAngle);
        const double psi =
                curvature * halfGap * fringeIntegral * (1.0 + sinE * sinE) / safeCos;
        return -curvature * Kokkos::tan(edgeAngle - psi);
    }

    /**
     * @brief Pure-value inputs for one bend field evaluation.
     *
     * This is not element state: the bend fills it once (coefficients fixed at
     * parse time, geometry from its Geometry, edge coefficients precomputed on the
     * host) and passes it to bendField() on both the device and host field paths.
     * Absent multipole coefficients are simply zero.
     */
    struct FieldInputs {
        double dipoleNormal;          ///< normal dipole B0 (normal[0])
        double dipoleSkew;            ///< skew dipole (skew[0])
        double quadNormal;            ///< normal quadrupole (normal[1])
        double quadSkew;              ///< skew quadrupole (skew[1])
        double bodyLength;            ///< magnet body length (arc length)
        double curvature;             ///< reference-path curvature h (0 for a straight body)
        double faceAngle;             ///< entrance pole-face angle E1 (frame tilt vs design orbit)
        double profileGap;            ///< Enge gap (GAP or 2*HGAP); 0 => hard edge
        double cosEntrance;           ///< |cos E1|, pole-face projection
        double cosExit;               ///< |cos E2|
        double entryEdgeCoefficient;  ///< distributed vertical edge-focusing (entry)
        double exitEdgeCoefficient;   ///< distributed vertical edge-focusing (exit)
    };

    /**
     * @brief Analytic bend magnetic field at a point in bend coordinates.
     *
     * @p arc is (radial offset x, vertical y, arc-length s) from toBendArcCoords();
     * s is the longitudinal coordinate (entrance at 0, exit at @p bodyLength). The
     * dipole and multipoles are scaled by the Enge amplitude F, the source-free
     * fringe corrections @f$-\tfrac12 B_0 F'' y^2@f$ and @f$B_0 F' y@f$ are added,
     * and the vertical edge focusing is distributed over the fringe as
     * @f$B_x \mathrel{+}= c\,F'\,y@f$. With @c profileGap = 0 this reduces to the
     * hard-edge field (F = 1 inside, edge terms 0).
     */
    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> bendField(
            const Vector_t<double, 3>& arc, const FieldInputs& in) {
        const double x = arc(0);
        const double y = arc(1);
        const double z = arc(2);

        const FringeAmplitude fringe = fringeAmplitude(
                -z * in.cosEntrance, (z - in.bodyLength) * in.cosExit, in.profileGap);
        const double scale = fringe.value;

        // Chain rule from the active face's distance coordinate to z.
        const double projection = (fringe.activeFace == 0) ? -in.cosEntrance : in.cosExit;
        const double dScale     = projection * fringe.firstDerivative;
        const double d2Scale    = projection * projection * fringe.secondDerivative;

        Vector_t<double, 3> B(0.0);
        // Dipole with the source-free fringe correction.
        B(1) += in.dipoleNormal * (scale - 0.5 * d2Scale * y * y);
        B(2) += in.dipoleNormal * dScale * y;
        B(0) -= scale * in.dipoleSkew;
        // Upright / skew quadrupole scaled by the fringe.
        B(0) += scale * in.quadNormal * y;
        B(1) += scale * in.quadNormal * x;
        B(0) -= scale * in.quadSkew * x;
        B(1) += scale * in.quadSkew * y;

        // Vertical edge focusing distributed over the fringe ramp.
        const double edgeCoefficient =
                (fringe.activeFace == 0) ? in.entryEdgeCoefficient : -in.exitEdgeCoefficient;
        B(0) += edgeCoefficient * dScale * y;

        return B;
    }

}  // namespace BendFieldModel

#endif  // OPALX_BendFieldModel_HH
