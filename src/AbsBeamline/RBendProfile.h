// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPALX_RBEND_PROFILE_H
#define OPALX_RBEND_PROFILE_H

#include <Kokkos_Core.hpp>
#include "VectorMath.h"

/// Value-only representation of a 1DProfile1 rectangular dipole field.
/// Parsing and strength calibration happen on the host; tracking uses this same
/// analytic evaluator on the host and device, without GSL or host pointers.
namespace RBendProfile {
    constexpr int maximumOrder = 15;
    struct Profile {
        int order                             = 5;
        double coefficients[maximumOrder + 1] = {0.478959, 1.911289,  -1.185953,
                                                 1.630554, -1.082657, 0.318111};
        double begin = -0.1, origin = 0.0, end = 0.1;
    };
    struct Inputs {
        bool enabled = false;
        Profile entrance, exit;
        double length = 0.0, gap = 0.02, amplitude = 0.0;
    };
    struct Value {
        double f, first, second;
    };

    /// Derivatives are with respect to the outward distance in metres.
    KOKKOS_INLINE_FUNCTION Value profile(double distance, double gap, const Profile& p) {
        const double u = distance / gap;
        double poly = p.coefficients[p.order], first = 0.0, second = 0.0;
        for (int i = p.order - 1; i >= 0; --i) {
            second = second * u + 2.0 * first;
            first  = first * u + poly;
            poly   = poly * u + p.coefficients[i];
        }
        if (poly > 80.0) return {0.0, 0.0, 0.0};
        if (poly < -80.0) return {1.0, 0.0, 0.0};
        // Stable logistic evaluation, including far into either saturated tail.
        const double e  = Kokkos::exp(-Kokkos::abs(poly));
        const double f  = poly >= 0.0 ? e / (1.0 + e) : 1.0 / (1.0 + e);
        const double df = -f * (1.0 - f);
        return {f, df * first / gap, df * (second + (2.0 * f - 1.0) * first * first) / (gap * gap)};
    }

    KOKKOS_INLINE_FUNCTION Value envelope(double z, const Inputs& in) {
        const bool entryRegion = z < in.entrance.end;
        const double ze        = z - in.length;
        const bool exitRegion  = ze >= in.exit.begin;
        if (entryRegion
            && (!exitRegion
                || Kokkos::abs(in.entrance.origin - z) < Kokkos::abs(ze - in.exit.origin))) {
            if (z < in.entrance.begin) return {0.0, 0.0, 0.0};
            Value v = profile(in.entrance.origin - z, in.gap, in.entrance);
            v.first = -v.first;
            return v;
        }
        if (exitRegion) {
            if (ze >= in.exit.end) return {0.0, 0.0, 0.0};
            return profile(ze - in.exit.origin, in.gap, in.exit);
        }
        return {1.0, 0.0, 0.0};
    }

    KOKKOS_INLINE_FUNCTION Vector_t<double, 3> field(
            const Vector_t<double, 3>& r, const Inputs& in) {
        const auto v = envelope(r(2), in);
        // Parallel rectangular pole faces: Maxwell fringe terms provide the
        // focusing. Adding the native FINT kick here would count it twice.
        return Vector_t<double, 3>(
                0.0, in.amplitude * (v.f - 0.5 * v.second * r(1) * r(1)),
                in.amplitude * v.first * r(1));
    }
}  // namespace RBendProfile
#endif
