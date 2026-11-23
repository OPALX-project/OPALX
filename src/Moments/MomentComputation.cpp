/**
 * @file MomentComputation.cpp
 * @brief Kokkos kernels and MPI reductions behind compute().
 *
 * Three passes cover all groups. The first-order pass reduces sums of R, P, gamma and kinetic
 * energy together with the position bounds; it reads only the arrays the request needs. The
 * second-order pass centres on the global means of the first pass (two-pass algorithm) and reduces
 * the 21 independent covariance entries, the kinetic-energy variance and the rest-frame
 * temperature sums. The polarization pass reduces sums and squares of Pol. Each pass runs one
 * kernel per container, joins the local results on the host and ends with one allreduce per MPI
 * operator; the bounds use a single MPI_MIN over (min, -max).
 */

#include "Moments/BeamMoments.h"

#include "Ippl.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>

namespace opalx::moments {

    namespace {
        using VectorView = ParticleContainer_t::particle_position_type::view_type;
        using SpinView   = ippl::ParticleAttrib<ParticleContainer_t::spin_vector_type>::view_type;
        using MassView   = ParticleContainer_t::qm_view_type;

        /// Mass in GeV for the kinetic energy: the reference mass if the container has one.
        struct KineticMass {
            MassView mass;
            bool scalarMass;
            bool useReference;
            double referenceGeV;

            KOKKOS_INLINE_FUNCTION double operator()(std::size_t i) const {
                return useReference ? referenceGeV : (scalarMass ? mass(0) : mass(i));
            }
        };

        KineticMass kineticMass(const ParticleContainer_t& particles) {
            const MassView mass       = particles.getMView();
            const PartData* reference = particles.getReference();
            return {mass, mass.extent(0) == 1, reference != nullptr,
                    reference ? reference->getM() * Units::eV2GeV : 0.0};
        }

        /// Sums [R (3), P (3), gamma, kinetic energy in MeV] and the position bounds.
        struct FirstOrderValue {
            double sum[8];
            double lowest[3];
            double highest[3];
        };

        struct FirstOrderPass {
            using value_type = FirstOrderValue;

            VectorView r;
            VectorView p;
            KineticMass mass;
            bool positions;
            bool momenta;
            bool bounds;

            KOKKOS_INLINE_FUNCTION void init(value_type& value) const {
                for (int k = 0; k < 8; ++k) {
                    value.sum[k] = 0.0;
                }
                for (int d = 0; d < 3; ++d) {
                    value.lowest[d]  = Kokkos::reduction_identity<double>::min();
                    value.highest[d] = Kokkos::reduction_identity<double>::max();
                }
            }

            KOKKOS_INLINE_FUNCTION void join(value_type& into, const value_type& from) const {
                for (int k = 0; k < 8; ++k) {
                    into.sum[k] += from.sum[k];
                }
                for (int d = 0; d < 3; ++d) {
                    into.lowest[d]  = Kokkos::min(into.lowest[d], from.lowest[d]);
                    into.highest[d] = Kokkos::max(into.highest[d], from.highest[d]);
                }
            }

            KOKKOS_INLINE_FUNCTION void accumulate(std::size_t i, value_type& value) const {
                if (positions || bounds) {
                    const Vector_t<double, 3> ri = r(i);
                    for (int d = 0; d < 3; ++d) {
                        if (positions) {
                            value.sum[d] += ri[d];
                        }
                        if (bounds) {
                            value.lowest[d]  = Kokkos::min(value.lowest[d], ri[d]);
                            value.highest[d] = Kokkos::max(value.highest[d], ri[d]);
                        }
                    }
                }
                if (momenta) {
                    const Vector_t<double, 3> pi = p(i);
                    for (int d = 0; d < 3; ++d) {
                        value.sum[3 + d] += pi[d];
                    }
                    const double gamma = Kokkos::sqrt(dot(pi) + 1.0);
                    value.sum[6] += gamma;
                    value.sum[7] += (gamma - 1.0) * mass(i) * Units::GeV2MeV;
                }
            }
        };

        /// [packed lower triangle of the 6x6 covariance (21), kinetic-energy variance,
        /// temperature sums (3)], contiguous so that one allreduce covers all of it.
        struct SecondOrderValue {
            static constexpr int energy      = 21;
            static constexpr int temperature = 22;
            static constexpr int size        = 25;
            double sum[size];
        };

        KOKKOS_INLINE_FUNCTION constexpr int packedIndex(int i, int j) {
            return i * (i + 1) / 2 + j;
        }

        struct SecondOrderPass {
            using value_type = SecondOrderValue;

            VectorView r;
            VectorView p;
            KineticMass mass;
            double meanR[3];
            double meanP[3];
            double meanGamma;
            double meanKineticEnergy;
            double betaC[3];
            double gammaC;
            double boostParallel;

            KOKKOS_INLINE_FUNCTION void init(value_type& value) const {
                for (int k = 0; k < value_type::size; ++k) {
                    value.sum[k] = 0.0;
                }
            }

            KOKKOS_INLINE_FUNCTION void join(value_type& into, const value_type& from) const {
                for (int k = 0; k < value_type::size; ++k) {
                    into.sum[k] += from.sum[k];
                }
            }

            KOKKOS_INLINE_FUNCTION void accumulate(std::size_t i, value_type& value) const {
                const Vector_t<double, 3> ri = r(i);
                const Vector_t<double, 3> pi = p(i);
                const double centred[6] = {ri[0] - meanR[0], pi[0] - meanP[0], ri[1] - meanR[1],
                                           pi[1] - meanP[1], ri[2] - meanR[2], pi[2] - meanP[2]};
                for (int a = 0; a < 6; ++a) {
                    for (int b = 0; b <= a; ++b) {
                        value.sum[packedIndex(a, b)] += centred[a] * centred[b];
                    }
                }

                const double gamma        = Kokkos::sqrt(dot(pi) + 1.0);
                const double massGeV      = mass(i);
                const double kinetic      = (gamma - 1.0) * massGeV * Units::GeV2MeV;
                const double energyOffset = kinetic - meanKineticEnergy;
                value.sum[value_type::energy] += energyOffset * energyOffset;

                // Boost into the centre-of-momentum frame, beta_c = <u> / <gamma>. Boosting the
                // mean 4-vector gives zero, so the boost acts on du = u - <u> and
                // dgamma = gamma - <gamma>:
                //   u' = du + [gamma_c^2 / (gamma_c + 1) (beta_c . du) - gamma_c dgamma] beta_c.
                const double betaDotDu =
                        betaC[0] * centred[1] + betaC[1] * centred[3] + betaC[2] * centred[5];
                const double shift    = boostParallel * betaDotDu - gammaC * (gamma - meanGamma);
                const double uRest[3] = {
                        centred[1] + shift * betaC[0], centred[3] + shift * betaC[1],
                        centred[5] + shift * betaC[2]};
                const double gammaRest = Kokkos::sqrt(
                        1.0 + uRest[0] * uRest[0] + uRest[1] * uRest[1] + uRest[2] * uRest[2]);
                for (int d = 0; d < 3; ++d) {
                    value.sum[value_type::temperature + d] +=
                            massGeV * uRest[d] * uRest[d] / gammaRest;
                }
            }
        };

        /// Sums [Pol (3), Pol^2 (3)].
        struct PolarizationValue {
            double sum[6];
        };

        struct PolarizationPass {
            using value_type = PolarizationValue;

            SpinView polarization;

            KOKKOS_INLINE_FUNCTION void init(value_type& value) const {
                for (int k = 0; k < 6; ++k) {
                    value.sum[k] = 0.0;
                }
            }

            KOKKOS_INLINE_FUNCTION void join(value_type& into, const value_type& from) const {
                for (int k = 0; k < 6; ++k) {
                    into.sum[k] += from.sum[k];
                }
            }

            KOKKOS_INLINE_FUNCTION void accumulate(std::size_t i, value_type& value) const {
                const auto pol = polarization(i);
                for (int d = 0; d < 3; ++d) {
                    value.sum[d] += pol[d];
                    value.sum[3 + d] += pol[d] * pol[d];
                }
            }
        };

        /// Adds the groups a requested group depends on.
        Request withDependencies(Request request) {
            if (request.contains(Moment::Spread)) {
                request |= Moment::MeanPosition | Moment::MeanMomentum;
            }
            return request;
        }

        /// Particles one work item reduces before it updates the Kokkos accumulator.
        constexpr std::size_t particlesPerItem = 32;

        /// Host backends (Serial, OpenMP) give each work item a contiguous run of particles.
        constexpr bool hostExecution =
                std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>;

        /**
         * @brief Reduces several particles per work item in registers before one join.
         *
         * Kokkos keeps the per-thread accumulator in memory (shared memory on GPUs) and the
         * compiler cannot rule out aliasing with the particle arrays, so a per-particle update
         * stores the whole value every time. A work item accumulates into a local value and joins
         * it once. On host backends item w takes the contiguous run [32 w, 32 w + 32); on devices
         * it takes w, w + W, ... so that neighbouring threads read neighbouring particles. The
         * pattern depends only on the local particle count and the backend, so every sum comes
         * out the same whichever other groups the request contains.
         */
        template <class Pass>
        struct Chunked {
            using value_type = typename Pass::value_type;

            Pass pass;
            std::size_t particles;
            std::size_t items;

            KOKKOS_INLINE_FUNCTION void init(value_type& value) const { pass.init(value); }
            KOKKOS_INLINE_FUNCTION void join(value_type& into, const value_type& from) const {
                pass.join(into, from);
            }
            KOKKOS_INLINE_FUNCTION void operator()(std::size_t item, value_type& value) const {
                value_type local;
                pass.init(local);
                if constexpr (hostExecution) {
                    const std::size_t end = Kokkos::min(particles, (item + 1) * particlesPerItem);
                    for (std::size_t i = item * particlesPerItem; i < end; ++i) {
                        pass.accumulate(i, local);
                    }
                } else {
                    for (std::size_t i = item; i < particles; i += items) {
                        pass.accumulate(i, local);
                    }
                }
                pass.join(value, local);
            }
        };

        /// Work items for @p particles; devices get at least as many as they run at once.
        std::size_t workItems(std::size_t particles) {
            const std::size_t chunks = (particles + particlesPerItem - 1) / particlesPerItem;
            if constexpr (hostExecution) {
                return chunks;
            } else {
                const auto resident =
                        static_cast<std::size_t>(Kokkos::DefaultExecutionSpace().concurrency());
                return std::max(chunks, std::min(particles, resident));
            }
        }

        /// Runs @p pass over every container and joins the local results into one value.
        template <class Pass, class MakePass>
        typename Pass::value_type reduceContainers(
                std::span<const ParticleContainer_t* const> containers, const char* label,
                MakePass makePass) {
            typename Pass::value_type total;
            Pass{}.init(total);
            for (const ParticleContainer_t* particles : containers) {
                const std::size_t localNum = particles->getLocalNum();
                const std::size_t items    = workItems(localNum);
                const Chunked<Pass> chunked{makePass(*particles), localNum, items};
                typename Pass::value_type local;
                Kokkos::parallel_reduce(label, Kokkos::RangePolicy<>(0, items), chunked, local);
                chunked.join(total, local);
            }
            return total;
        }
    }  // namespace

    BeamMoments compute(std::span<const ParticleContainer_t* const> containers, Request request) {
        static IpplTimings::TimerRef timer = IpplTimings::getTimer("computeMoments");
        IpplTimings::startTimer(timer);

        BeamMoments moments;
        for (const ParticleContainer_t* particles : containers) {
            if (particles == nullptr) {
                throw OpalException("opalx::moments::compute", "Null particle container.");
            }
            moments.count_m += particles->getTotalNum();
        }
        // An empty set divides by one, which keeps every mean and spread at zero.
        const double count   = static_cast<double>(std::max<std::size_t>(moments.count_m, 1));
        const Request groups = withDependencies(request);

        const bool positions = groups.contains(Moment::MeanPosition);
        const bool momenta   = groups.contains(Moment::MeanMomentum);
        const bool bounds    = groups.contains(Moment::Bounds);
        if (positions || momenta || bounds) {
            const FirstOrderValue local = reduceContainers<FirstOrderPass>(
                    containers, "opalx::moments::firstOrder", [&](const ParticleContainer_t& pc) {
                        return FirstOrderPass{pc.R.getView(), pc.P.getView(), kineticMass(pc),
                                              positions,      momenta,        bounds};
                    });
            if (positions || momenta) {
                double sum[8];
                ippl::Comm->allreduce(&local.sum[0], &sum[0], 8, std::plus<double>());
                for (unsigned d = 0; d < 3; ++d) {
                    moments.meanR_m[d] = sum[d] / count;
                    moments.meanP_m[d] = sum[3 + d] / count;
                }
                moments.meanGamma_m         = sum[6] / count;
                moments.meanKineticEnergy_m = sum[7] / count;
            }
            if (bounds) {
                double extremes[6];
                for (unsigned d = 0; d < 3; ++d) {
                    extremes[d]     = local.lowest[d];
                    extremes[3 + d] = -local.highest[d];
                }
                ippl::Comm->allreduce(&extremes[0], 6, std::less<double>());
                for (unsigned d = 0; d < 3; ++d) {
                    // No particle on any rank leaves min > max; report a zero-size box at the
                    // origin so that meshes built from it stay valid.
                    const bool empty  = extremes[d] > -extremes[3 + d];
                    moments.minR_m[d] = empty ? 0.0 : extremes[d];
                    moments.maxR_m[d] = empty ? 0.0 : -extremes[3 + d];
                }
            }
        }

        if (groups.contains(Moment::Spread)) {
            SecondOrderPass centre{};
            for (unsigned d = 0; d < 3; ++d) {
                centre.meanR[d] = moments.meanR_m[d];
                centre.meanP[d] = moments.meanP_m[d];
                centre.betaC[d] =
                        moments.meanGamma_m > 0.0 ? moments.meanP_m[d] / moments.meanGamma_m : 0.0;
            }
            centre.meanGamma         = moments.meanGamma_m;
            centre.meanKineticEnergy = moments.meanKineticEnergy_m;
            const double betaC2      = centre.betaC[0] * centre.betaC[0]
                                  + centre.betaC[1] * centre.betaC[1]
                                  + centre.betaC[2] * centre.betaC[2];
            centre.gammaC        = 1.0 / std::sqrt(1.0 - betaC2);
            centre.boostParallel = centre.gammaC * centre.gammaC / (centre.gammaC + 1.0);

            const SecondOrderValue local = reduceContainers<SecondOrderPass>(
                    containers, "opalx::moments::secondOrder", [&](const ParticleContainer_t& pc) {
                        SecondOrderPass pass = centre;
                        pass.r               = pc.R.getView();
                        pass.p               = pc.P.getView();
                        pass.mass            = kineticMass(pc);
                        return pass;
                    });
            SecondOrderValue global;
            ippl::Comm->allreduce(
                    &local.sum[0], &global.sum[0], SecondOrderValue::size, std::plus<double>());
            for (int a = 0; a < 6; ++a) {
                for (int b = 0; b <= a; ++b) {
                    moments.covariance_m[a][b] = global.sum[packedIndex(a, b)] / count;
                    moments.covariance_m[b][a] = moments.covariance_m[a][b];
                }
            }
            moments.stdKineticEnergy_m = std::sqrt(global.sum[SecondOrderValue::energy] / count);
            for (unsigned d = 0; d < 3; ++d) {
                moments.restFrameTemperature_m[d] = global.sum[SecondOrderValue::temperature + d]
                                                    * Units::GeV2eV / count / Physics::kB;
            }
        }

        const bool allSpin = std::all_of(
                containers.begin(), containers.end(), [](const ParticleContainer_t* particles) {
                    return particles->hasSpin();
                });
        const Request polarization =
                groups.contains(Moment::Polarization) && allSpin ? Moment::Polarization : Request{};
        if (!polarization.empty()) {
            const PolarizationValue local = reduceContainers<PolarizationPass>(
                    containers, "opalx::moments::polarization", [](const ParticleContainer_t& pc) {
                        return PolarizationPass{pc.Pol.getView()};
                    });
            double sum[6];
            ippl::Comm->allreduce(&local.sum[0], &sum[0], 6, std::plus<double>());
            for (unsigned d = 0; d < 3; ++d) {
                const double mean             = sum[d] / count;
                const double variance         = sum[3 + d] / count - mean * mean;
                moments.meanPolarization_m[d] = mean;
                moments.rmsPolarization_m[d]  = variance > 0.0 ? std::sqrt(variance) : 0.0;
            }
        }

        moments.available_m = groups.without(Moment::Polarization) | polarization;
        IpplTimings::stopTimer(timer);
        return moments;
    }

    BeamMoments compute(const ParticleContainer_t& particles, Request request) {
        const ParticleContainer_t* const single[] = {&particles};
        return compute(std::span<const ParticleContainer_t* const>(single), request);
    }

}  // namespace opalx::moments
