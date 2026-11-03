// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
#ifndef OPALX_SCALING_FFA_SUPPORT_H
#define OPALX_SCALING_FFA_SUPPORT_H

#include "OPALTypes.h"

/** Exact spatial support of a scaling spiral FFA, in its entrance frame.
 * Shared by the host element and device ring candidate selector. Angular
 * extents follow the logarithmic spiral, not an entrance-frame z slab.
 */
struct ScalingFFASupport {
    double tanDelta_m        = 0.;
    double r0_m              = 1.;
    double rMin_m            = 0.;
    double rMax_m            = 0.;
    double phiStart_m        = 0.;
    double azimuthalExtent_m = 0.;
    double verticalExtent_m  = 0.;

    /** Shared cylindrical chart: [radius, vertical, azimuth, r/|R0|, spiral].
     * Field evaluation, support selection and aperture checks use this one
     * conversion; callers with an existing chart can use the overloads below.
     */
    KOKKOS_INLINE_FUNCTION void getCylindricalCoordinates(
            const Vector_t<double, 3>& point, Vector_t<double, 5>& cylindrical) const {
        const double x          = point(0) + r0_m;
        const double r          = Kokkos::sqrt(x * x + point(2) * point(2));
        const double normalized = Kokkos::abs(r / r0_m);
        const double phi        = Kokkos::atan2(point(2), x * Kokkos::copysign(1., r0_m));
        cylindrical             = {
                r, point(1), phi, normalized,
                phi - tanDelta_m * Kokkos::log(normalized) - phiStart_m};
    }

    KOKKOS_INLINE_FUNCTION bool contains(const Vector_t<double, 5>& cylindrical) const {
        const double r = cylindrical(0), spiral = cylindrical(4);
        return r > 0 && r >= rMin_m && r <= rMax_m
               && Kokkos::abs(cylindrical(1)) <= verticalExtent_m && spiral >= -azimuthalExtent_m
               && spiral <= azimuthalExtent_m;
    }

    KOKKOS_INLINE_FUNCTION bool contains(const Vector_t<double, 3>& point) const {
        Vector_t<double, 5> cylindrical;
        getCylindricalCoordinates(point, cylindrical);
        return contains(cylindrical);
    }

    /** Aperture loss independent of field selection, for static orbit optics.
     * Only the spiral angular window owns these radial/vertical limits. A ray
     * outside this window may legitimately be in another magnet or a gap.
     */
    KOKKOS_INLINE_FUNCTION bool outsideAperture(const Vector_t<double, 5>& cylindrical) const {
        const double r = cylindrical(0), spiral = cylindrical(4);
        if (!(r > 0)) return true;  // No cylindrical chart exists at the ring centre.
        if (spiral < -azimuthalExtent_m || spiral > azimuthalExtent_m) return false;
        return r < rMin_m || r > rMax_m || Kokkos::abs(cylindrical(1)) > verticalExtent_m;
    }

    KOKKOS_INLINE_FUNCTION bool outsideAperture(const Vector_t<double, 3>& point) const {
        Vector_t<double, 5> cylindrical;
        getCylindricalCoordinates(point, cylindrical);
        return outsideAperture(cylindrical);
    }
};
#endif
