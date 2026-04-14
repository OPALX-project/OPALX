#include "Processes/GlobalProcesses/Decay.h"

#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "Utilities/Options.h"

#include <cmath>
#include <cstdint>

namespace {

uint64_t decayRngSeed(std::size_t containerIndex) {
    const uint64_t seedBase =
        static_cast<uint64_t>((Options::seed == -1) ? 1234567 : Options::seed);
    return seedBase + 104729ULL * static_cast<uint64_t>(containerIndex + 1)
           + 8191ULL * static_cast<uint64_t>(ippl::Comm->rank() + 1);
}

}  // namespace

Decay::Decay(double restLifetimeSeconds, std::size_t containerIndex, double parentMassGeV)
    : tau0_m(restLifetimeSeconds),
      randPool_m(decayRngSeed(containerIndex)),
      parentMassGeV_m(parentMassGeV) {
}

void Decay::setDaughterContainer(std::shared_ptr<ParticleContainer<double, 3>> daughterPC,
                                 double daughterMassGeV) {
    daughterPC_m = std::move(daughterPC);
    daughterMassGeV_m = daughterMassGeV;
}

size_t Decay::apply(ParticleContainer<double, 3>& pc,
                    double dt,
                    long long globalTrackStep,
                    size_t containerIdx) {
    using pc_size_type = ippl::detail::size_type;

    if (dt <= 0.0 || !(tau0_m > 0.0) || !std::isfinite(tau0_m)) {
        return 0;
    }

    const pc_size_type nLocal = pc.getLocalNum();
    if (nLocal == 0 && !daughterPC_m) {
        return 0;
    }

    (void)globalTrackStep;
    (void)containerIdx;

    const auto pool = randPool_m;

    /* Phase 1: Mark decayed particles using relativistic decay probability. */
    Kokkos::View<bool*> invalid("Decay::invalid", nLocal);
    auto Pview = pc.P.getView();
    const double tau0 = tau0_m;

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
            if (decayed) {
                ++count;
            }
            pool.free_state(generator);
        },
        localDestroyNum);
    Kokkos::fence();

    pc_size_type globalDestroyNum = 0;
    ippl::Comm->allreduce(localDestroyNum, globalDestroyNum, 1, std::plus<pc_size_type>());
    if (globalDestroyNum == 0) {
        return 0;
    }

    /* Phase 2 & 3: Create daughter particles (if a daughter container is set). */
    if (daughterPC_m) {
        /* Phase 2: Prefix scan to collect compact indices of decayed particles. */
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

        auto Rview = pc.R.getView();
        auto dtView = pc.dt.getView();

        /* Collect R, P, dt of decayed parents into compact temporary views. */
        using vector_view_type = Kokkos::View<ippl::Vector<double, 3>*>;
        vector_view_type parentR("Decay::parentR", localDestroyNum);
        vector_view_type parentP("Decay::parentP", localDestroyNum);
        Kokkos::View<double*> parentDt("Decay::parentDt", localDestroyNum);

        Kokkos::parallel_for(
            "Decay::collectParents",
            nLocal,
            KOKKOS_LAMBDA(const pc_size_type i) {
                if (invalid(i)) {
                    const pc_size_type j = compactIdx(i);
                    parentR(j) = Rview(i);
                    parentP(j) = Pview(i);
                    parentDt(j) = dtView(i);
                }
            });
        Kokkos::fence();

        /* Phase 3: Create daughters — subclass-specific momentum sampling. */
        const pc_size_type oldDaughterLocal = daughterPC_m->getLocalNum();
        daughterPC_m->create(localDestroyNum);

        if (localDestroyNum > 0) {
            createDaughterParticles(localDestroyNum, oldDaughterLocal,
                                    parentR, parentP, parentDt);
        }
    }

    /* Phase 4: Destroy decayed parent particles. */
    pc.destroy(invalid, localDestroyNum);
    return static_cast<size_t>(globalDestroyNum);
}
