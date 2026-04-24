#ifndef OPAL_GLOBALPROCESS_PION_DECAY_H
#define OPAL_GLOBALPROCESS_PION_DECAY_H

#include "Physics/ParticleProperties.h"
#include "Processes/GlobalProcesses/Decay.h"

/**
 * @brief Charged pion decay: pi -> mu + nu_mu (two-body).
 *
 * @note In the parent rest frame the daughter momentum is fixed:
 *   p = (M_pi^2 - m_mu^2) / (2 M_pi), with isotropic direction.
 *   Daughters are Lorentz-boosted to the lab frame.
 */
class PionDecay : public Decay {
public:
    PionDecay(double restLifetimeSeconds, std::size_t containerIndex, double parentMassGeV)
        : Decay(restLifetimeSeconds, containerIndex, parentMassGeV) {
        allowedDaughterSpecies_m = static_cast<short>(ParticleProperties::getParticleType("MUON"));
    }

    void createDaughterParticles(
            std::size_t localDestroyNum, std::size_t oldDaughterLocal,
            const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
            const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
            const Kokkos::View<double*>& parentDt) override;
};

#endif
