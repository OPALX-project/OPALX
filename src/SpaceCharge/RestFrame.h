/**
 * @file RestFrame.h
 * @brief Rest-frame Lorentz boost shared by the binned CartesianPIC3D and Barnes-Hut solves.
 */

#ifndef OPALX_SPACE_CHARGE_REST_FRAME_H
#define OPALX_SPACE_CHARGE_REST_FRAME_H

#include "Ippl.h"
#include "Physics/Physics.h"

#include <Kokkos_Core.hpp>

#include <array>
#include <functional>

namespace opalx::spacecharge {

    /**
     * @brief Field transformation from the rest frame of a bunch (or bin) to the lab frame.
     *
     * The frame moves with v = c <P> / gamma, where <P> is the mean normalized momentum
     * (beta gamma) in Cartesian solve axes and gamma the mean of the particles' gamma. A
     * rest-frame electrostatic field E' maps to E = gamma E' - (gamma - 1)(E'.v^)v^ and
     * B = (gamma / c^2) v x E'; image charges of a grounded plane move the other way and use
     * the opposite B sign.
     */
    struct RestFrameBoost {
        using Vector = ippl::Vector<double, 3>;

        Vector velocity{0.0, 0.0, 0.0};      ///< m/s, solve axes.
        Vector velocityUnit{0.0, 0.0, 0.0};  ///< Zero for a bunch at rest.
        double gamma        = 1.0;
        double gammaMinusOne = 0.0;
        double gammaOverCSq  = 0.0;  ///< gamma / c^2.

        static RestFrameBoost fromMoments(const std::array<double, 3>& meanMomentum, double gamma) {
            RestFrameBoost boost;
            const double invGamma = (gamma > 0.0) ? (1.0 / gamma) : 0.0;
            Vector momentum(0.0);
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                momentum[dimension] = meanMomentum[dimension];
            }
            boost.velocity            = Physics::c * momentum * invGamma;
            const double velocityNorm = Kokkos::sqrt(boost.velocity.dot(boost.velocity));
            boost.velocityUnit =
                    (velocityNorm > 0.0) ? Vector(boost.velocity / velocityNorm) : Vector(0.0);
            boost.gamma         = gamma;
            boost.gammaMinusOne = gamma - 1.0;
            boost.gammaOverCSq  = gamma / (Physics::c * Physics::c);
            return boost;
        }

        [[nodiscard]] bool moving() const {
            return velocity[0] != 0.0 || velocity[1] != 0.0 || velocity[2] != 0.0;
        }

        /** @brief Lab-frame E of the rest-frame field @p electricPrime. */
        KOKKOS_INLINE_FUNCTION Vector electric(const Vector& electricPrime) const {
            const double electricDotVelocityUnit = electricPrime.dot(velocityUnit);
            return gamma * electricPrime - gammaMinusOne * electricDotVelocityUnit * velocityUnit;
        }

        /** @brief Lab-frame B of the rest-frame field @p electricPrime; @p magneticSign is -1 for
         * image charges. */
        KOKKOS_INLINE_FUNCTION Vector magnetic(
                const Vector& electricPrime, double magneticSign) const {
            return magneticSign * gammaOverCSq * cross(velocity, electricPrime);
        }
    };

    /**
     * @brief Rest frame of the first @p count particles of every rank: mean momentum and mean
     * gamma over all ranks (the moments of a single CartesianPIC3D bin holding every particle).
     */
    template <class MomentumView>
    RestFrameBoost restFrameOf(const MomentumView& momentum, std::size_t count) {
        ippl::Vector<double, 4> sums(0.0);
        Kokkos::parallel_reduce(
                "restFrameOf", count,
                KOKKOS_LAMBDA(const std::size_t i, ippl::Vector<double, 4>& sum) {
                    const ippl::Vector<double, 3> p = momentum(i);
                    for (unsigned dimension = 0; dimension < 3; ++dimension) {
                        sum[dimension] += p[dimension];
                    }
                    sum[3] += Kokkos::sqrt(1.0 + p.dot(p));
                },
                sums);
        std::array<double, 5> global{sums[0], sums[1], sums[2], sums[3], static_cast<double>(count)};
        ippl::Comm->allreduce(global.data(), 5, std::plus<double>());
        if (!(global[4] > 0.0)) {
            return RestFrameBoost{};
        }
        const std::array<double, 3> meanMomentum{
                global[0] / global[4], global[1] / global[4], global[2] / global[4]};
        return RestFrameBoost::fromMoments(meanMomentum, global[3] / global[4]);
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_REST_FRAME_H
