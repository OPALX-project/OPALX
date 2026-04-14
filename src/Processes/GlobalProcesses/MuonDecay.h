#ifndef OPAL_GLOBALPROCESS_MUON_DECAY_H
#define OPAL_GLOBALPROCESS_MUON_DECAY_H

#include "Processes/GlobalProcesses/Decay.h"

/**
 * @brief Muon decay: mu -> e + nu_e + nu_mu (three-body).
 *
 * @note Daughter electron energy is sampled from the Michel spectrum
 *   f(x) = 2 x^2 (3 - 2x), where x = E_e / E_max and E_max = m_mu / 2.
 *   The direction in the parent rest frame is isotropic; daughters are
 *   Lorentz-boosted to the lab frame.
 */
class MuonDecay : public Decay {
public:
    using Decay::Decay;

protected:
    void createDaughterParticles(
        std::size_t localDestroyNum,
        std::size_t oldDaughterLocal,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
        const Kokkos::View<double*>& parentDt) override;
};

#endif
