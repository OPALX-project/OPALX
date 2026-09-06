#ifndef OPALX_CYCLOTRON_RF_PROFILE_H
#define OPALX_CYCLOTRON_RF_PROFILE_H

#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <string>
#include <utility>
#include "OPALTypes.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

/** @brief Legacy cyclotron voltage profile, not a continuous electromagnetic map.
 * Input is a point count followed by (u, V, dV/du), all dimensionless, with
 * u=(gap coordinate-RMIN)/(RMAX-RMIN). Cubic Hermite interpolation uses the
 * supplied derivatives without recomputing them. Both endpoints belong to the
 * support. Views are shared and must be treated as read-only after loading.
 * Host loading requires initialized Kokkos; evaluation transfers no data.
 */
class CyclotronRFProfile {
public:
    using View = Kokkos::View<double**>;
    View data;
    decltype(Kokkos::create_mirror_view(std::declval<View>())) host;

    /// Read a finite, strictly increasing grid spanning [0,1]; malformed input throws.
    explicit CyclotronRFProfile(const std::string& filename) {
        std::ifstream input(filename);
        int count = 0;
        auto fail = [&]() {
            throw OpalException("CyclotronRFProfile", "Invalid radial RF profile: " + filename);
        };
        if (!(input >> count) || count < 2 || count > 1000000) fail();
        data = View("cyclotron_rf_profile", count, 3);
        host = Kokkos::create_mirror_view(data);
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < 3; ++j)
                if (!(input >> host(i, j)) || !std::isfinite(host(i, j))) fail();
            if (i && host(i, 0) <= host(i - 1, 0)) fail();
        }
        if (host(0, 0) != 0 || host(count - 1, 0) != 1) fail();
        std::string extra;
        if (input >> extra) fail();
        Kokkos::deep_copy(data, host);
    }

    /** @brief Evaluate the profile and its derivative with respect to u.
     * @param grid Read-only profile view accessible in the calling execution space.
     * @param u Normalized gap coordinate; outside [0,1] contributes zero.
     * @param[out] value Dimensionless voltage multiplier.
     * @param[out] derivative Dimensionless dV/du.
     * @return Whether u belongs to the finite profile support.
     */
    template<class V>
    KOKKOS_INLINE_FUNCTION static bool evaluate(
            const V& grid, double u, double& value, double& derivative) {
        value = derivative = 0;
        if (!(u >= 0 && u <= 1)) return false;
        int low = 0, high = static_cast<int>(grid.extent(0)) - 1;
        while (high - low > 1) {
            const int mid = (low + high) / 2;
            if (u < grid(mid, 0)) high = mid;
            else low = mid;
        }
        const double dx = grid(high, 0) - grid(low, 0);
        const double dy = grid(high, 1) - grid(low, 1);
        const double a = grid(low, 2), b = grid(high, 2);
        const double c2 = 3 * dy - dx * (b + 2 * a);
        const double c3 = -2 * dy + dx * (b + a);
        const double t = (u - grid(low, 0)) / dx;
        value = grid(low, 1) + t * dx * a + t * t * c2 + t * t * t * c3;
        derivative = a + 2 * t / dx * c2 + 3 * t * t / dx * c3;
        return true;
    }
};

/** @brief Device-copyable, static-frequency SINGLEGAP kick for median-plane protons.
 * Local x runs along the gap; local z is the positive crossing direction and y
 * is vertical. Parameters use SI: voltage [V], omega [rad/s], phase [rad], width
 * and profileLength [m]. Momentum is mechanical beta*gamma, mass is rest energy
 * [eV]. The phase convention is omega*time-phase (time in seconds).
 *
 * With interpolated f(u), gamma gains V*f*cos(phase)*sinc(omega*width/(2*c*beta))/m.
 * The legacy focusing rotation is -f'(u)*V*sin(phase)*c /
 * (profileLength*omega*2*pi*|p|*m). The additional 2*pi and omission of the transit
 * factor from this rotation deliberately reproduce OPAL 2022.1, not a rederived
 * focusing model. Radial momentum is preserved before the rotation.
 *
 * This first contract is restricted to positive, forward, median-plane protons;
 * invalid states return false without modifying momentum. It neither detects
 * crossings nor advances time, and must be called once per accepted gap event.
 */
struct CyclotronRFKick {
    double voltage = 0, omega = 0, phase = 0;
    double width = 0, profileLength = 0;

    /// Validate host configuration before launching kernels. Zero width is allowed.
    void validate() const {
        if (!std::isfinite(voltage) || !std::isfinite(omega) || !std::isfinite(phase)
            || !std::isfinite(width) || !std::isfinite(profileLength)
            || omega <= 0 || width < 0 || profileLength <= 0)
            throw OpalException("CyclotronRFKick", "Invalid voltage, frequency or gap dimensions.");
    }

    /** @brief Apply one whole impulse; no fractional kick to meet an energy target.
     * @param grid Validated voltage profile in the calling memory space.
     * @param u Normalized gap coordinate [0,1].
     * @param time Absolute RF time [s].
     * @param mass Proton rest energy [eV].
     * @param[in,out] p Local mechanical momentum [beta*gamma].
     * @return False for unsupported/out-of-support states or an unphysical energy kick.
     */
    template<class V>
    KOKKOS_INLINE_FUNCTION bool apply(
            const V& grid, double u, double time, double mass, Vector_t<double, 3>& p) const {
        double f, derivative;
        if (!(mass > 0) || !Kokkos::isfinite(mass) || !Kokkos::isfinite(time)
            || !Kokkos::isfinite(p[0]) || !Kokkos::isfinite(p[2])
            || p[1] != 0 || !(p[2] > 0)
            || !CyclotronRFProfile::evaluate(grid, u, f, derivative)) return false;
        if (voltage == 0) return true;
        const double p2 = p[0] * p[0] + p[2] * p[2];
        const double bg = Kokkos::sqrt(p2), gamma = Kokkos::sqrt(1 + p2);
        const double beta = bg / gamma;
        const double a = omega * width / (2 * Physics::c * beta);
        const double transit = Kokkos::abs(a) < 1e-8
                ? 1 - a * a / 6 : Kokkos::sin(a) / a;
        const double phi = omega * time - phase;
        const double finalGamma = gamma + voltage * f * transit * Kokkos::cos(phi) / mass;
        const double longitudinal2 = finalGamma * finalGamma - 1 - p[0] * p[0];
        if (!(finalGamma >= 1) || !(longitudinal2 > 0)
            || !Kokkos::isfinite(longitudinal2)) return false;
        const double longitudinal = Kokkos::sqrt(longitudinal2);
        const double rotation = -derivative * voltage / profileLength * Kokkos::sin(phi)
                / (omega * Physics::two_pi) / (bg * mass / Physics::c);
        const double radial = p[0];
        p[0] = Kokkos::cos(rotation) * radial + Kokkos::sin(rotation) * longitudinal;
        p[2] = -Kokkos::sin(rotation) * radial + Kokkos::cos(rotation) * longitudinal;
        return true;
    }
};
#endif
