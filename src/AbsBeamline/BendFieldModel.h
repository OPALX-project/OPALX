#ifndef OPALX_BendFieldModel_HH
#define OPALX_BendFieldModel_HH

//
// BendFieldModel
//   Stateless field-shape math shared by the analytic SBEND and RBEND fringe
//   fields: the OPAL default Enge longitudinal profile (amplitude and its first
//   two derivatives), the combined amplitude from the two pole-face distances,
//   and the pole-face vertical edge-focusing coefficient. All functions are
//   device-callable and hold no state; the bends supply their own geometry
//   (pole-face distances, gap, curvature) and coefficients.
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

    /// |cos(angle)| floored to a small positive value, for pole-face projection.
    KOKKOS_INLINE_FUNCTION double safeAbsCos(const double angle) {
        const double c = Kokkos::abs(Kokkos::cos(angle));
        return (c > 1.0e-6) ? c : 1.0e-6;
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

}  // namespace BendFieldModel

#endif  // OPALX_BendFieldModel_HH
