#include "Processes/GlobalProcesses/Decay.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <cstdint>
#include <string>

namespace {

uint64_t decayRngSeed(std::size_t containerIndex) {
    const uint64_t seedBase =
        static_cast<uint64_t>((Options::seed == -1) ? 1234567 : Options::seed);
    return seedBase + 104729ULL * static_cast<uint64_t>(containerIndex + 1)
           + 8191ULL * static_cast<uint64_t>(ippl::Comm->rank() + 1);
}

/// Compact views of the kinematics of parents marked for decay.
struct DecayedParentViews {
    Kokkos::View<ippl::Vector<double, 3>*> R;
    Kokkos::View<ippl::Vector<double, 3>*> P;
    Kokkos::View<double*> dt;
};

/**
 * @brief Gather R/P/dt of parents marked for decay into compact views.
 * @param nLocal          Number of local particles in the parent container.
 * @param localDestroyNum Number of local particles flagged in @p invalid.
 * @param invalid         Per-particle decay flag (length @p nLocal).
 * @param Rview           View of parent positions.
 * @param Pview           View of parent momenta.
 * @param dtView          View of parent per-particle time steps.
 * @return Three length-@p localDestroyNum views packed contiguously.
 *
 * @note Uses a parallel_scan to assign each marked particle a unique compact
 *   index, then a parallel_for to gather R/P/dt into the returned views.
 */
template <typename RView, typename PView, typename DtView>
DecayedParentViews collectDecayedParents(
    ippl::detail::size_type nLocal,
    ippl::detail::size_type localDestroyNum,
    Kokkos::View<bool*>     invalid,
    RView                   Rview,
    PView                   Pview,
    DtView                  dtView) {
    using pc_size_type = ippl::detail::size_type;

    Kokkos::View<pc_size_type*> compactIdx("Decay::compactIdx", nLocal);
    Kokkos::parallel_scan(
        "Decay::compact",
        nLocal,
        KOKKOS_LAMBDA(const pc_size_type i, pc_size_type& runningTotal, const bool isFinal) {
            if (invalid(i)) {
                if (isFinal) {
                    compactIdx(i) = runningTotal;
                }
                ++runningTotal;
            }
        });
    Kokkos::fence();

    DecayedParentViews out{
        Kokkos::View<ippl::Vector<double, 3>*>("Decay::parentR", localDestroyNum),
        Kokkos::View<ippl::Vector<double, 3>*>("Decay::parentP", localDestroyNum),
        Kokkos::View<double*>("Decay::parentDt", localDestroyNum)};

    Kokkos::parallel_for(
        "Decay::collectParents",
        nLocal,
        KOKKOS_LAMBDA(const pc_size_type i) {
            if (invalid(i)) {
                const pc_size_type j = compactIdx(i);
                out.R(j)  = Rview(i);
                out.P(j)  = Pview(i);
                out.dt(j) = dtView(i);
            }
        });
    Kokkos::fence();
    return out;
}

/**
 * @brief Mark particles for decay using the relativistic exponential law.
 * @param nLocal   Number of local particles in the parent container.
 * @param dt       Lab-frame time step [s].
 * @param tau0     Rest-frame mean lifetime [s].
 * @param pool     RNG pool (captured by value into the device lambda).
 * @param Pview    View of parent momenta [beta*gamma].
 * @param invalid  Output view, one bool per local particle — true if the particle decayed.
 * @return Local count of particles marked for decay.
 */
template <typename PView>
ippl::detail::size_type markDecayedParticles(
    ippl::detail::size_type nLocal,
    double dt,
    double tau0,
    Kokkos::Random_XorShift64_Pool<> pool,
    PView Pview,
    Kokkos::View<bool*> invalid) {
    using pc_size_type = ippl::detail::size_type;

    pc_size_type localDestroyNum = 0;
    Kokkos::parallel_reduce(
        "Decay::mark",
        nLocal,
        KOKKOS_LAMBDA(const pc_size_type i, pc_size_type& count) {
            auto generator = pool.get_state();
            const double p2 =
                Pview(i)[0] * Pview(i)[0] + Pview(i)[1] * Pview(i)[1] + Pview(i)[2] * Pview(i)[2];
            const double gamma  = Kokkos::sqrt(1.0 + p2);
            const double tauLab = gamma * tau0;
            const double pDecay = 1.0 - Kokkos::exp(-dt / tauLab);
            const bool decayed  = generator.drand(0.0, 1.0) < pDecay;
            invalid(i)          = decayed;
            count += decayed;
            pool.free_state(generator);
        },
        localDestroyNum);
    Kokkos::fence();
    return localDestroyNum;
}

}  // namespace

Decay::Decay(double restLifetimeSeconds, std::size_t containerIndex, double parentMassGeV)
    : tau0_m(restLifetimeSeconds),
      randPool_m(decayRngSeed(containerIndex)),
      parentMassGeV_m(parentMassGeV) {
}

void Decay::setDaughterContainer(std::shared_ptr<ParticleContainer<double, 3>> daughterPC,
                                 double daughterMassGeV) {
    if (daughterPC && daughterPC->Sp != allowedDaughterSpecies_m) {
        const auto expected = static_cast<ParticleType>(allowedDaughterSpecies_m);
        const auto actual   = static_cast<ParticleType>(daughterPC->Sp);
        throw OpalException(
            "Decay::setDaughterContainer",
            "Daughter container species mismatch: expected \""
                + ParticleProperties::getParticleTypeString(expected)
                + "\" but got \""
                + ParticleProperties::getParticleTypeString(actual) + "\".");
    }
    daughterPC_m = std::move(daughterPC);
    daughterMassGeV_m = daughterMassGeV;
}

size_t Decay::apply(ParticleContainer<double, 3>& pc,
                    double dt,
                    [[maybe_unused]] long long globalTrackStep,
                    [[maybe_unused]] size_t containerIdx) {
    using pc_size_type = ippl::detail::size_type;

    if (dt <= 0.0 || !(tau0_m > 0.0) || !std::isfinite(tau0_m)) {
        return 0;
    }

    const pc_size_type nLocal = pc.getLocalNum();
    if (nLocal == 0 && !daughterPC_m) {
        return 0;
    }

    /* Phase 1: Mark decayed particles using relativistic decay probability. */
    Kokkos::View<bool*> invalid("Decay::invalid", nLocal);
    const pc_size_type localDestroyNum =
        markDecayedParticles(nLocal, dt, tau0_m, randPool_m, pc.P.getView(), invalid);

    pc_size_type globalDestroyNum = 0;
    ippl::Comm->allreduce(localDestroyNum, globalDestroyNum, 1, std::plus<pc_size_type>());
    if (globalDestroyNum == 0) {
        return 0;
    }

    if (daughterPC_m) {
        /* Phase 2: Gather kinematics of decayed parents into compact views. */
        const DecayedParentViews parents = collectDecayedParents(
            nLocal, localDestroyNum, invalid,
            pc.R.getView(), pc.P.getView(), pc.dt.getView());

        /* Phase 3: Create daughters — subclass-specific momentum sampling. */
        const pc_size_type oldDaughterLocal = daughterPC_m->getLocalNum();
        daughterPC_m->create(localDestroyNum);

        if (localDestroyNum > 0) {
            createDaughterParticles(localDestroyNum, oldDaughterLocal,
                                    parents.R, parents.P, parents.dt);
        }
    }

    /* Phase 4: Destroy decayed parent particles. */
    pc.destroy(invalid, localDestroyNum);
    return static_cast<size_t>(globalDestroyNum);
}
