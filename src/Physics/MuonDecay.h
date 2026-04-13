#ifndef OPAL_PHYSICS_MUON_DECAY_H
#define OPAL_PHYSICS_MUON_DECAY_H

#include "Physics/Physics.h"

#include <Kokkos_Core.hpp>

namespace Physics {
namespace MuonDecay {

    /// Maximum electron energy in the muon rest frame [GeV]: \f$ E_\text{max} = m_\mu / 2 \f$.
    KOKKOS_INLINE_FUNCTION constexpr double maxElectronEnergy() {
        return Physics::m_mu / 2.0;
    }

    /// Minimum reduced energy \f$ x_\text{min} = m_e / E_\text{max} \f$.
    KOKKOS_INLINE_FUNCTION constexpr double minElectronX() {
        return Physics::m_e / maxElectronEnergy();
    }

    /// Unnormalized Michel spectrum: \f$ f(x) = 2 x^2 (3 - 2x) \f$ for \f$ x \in [x_\text{min}, 1] \f$.
    KOKKOS_INLINE_FUNCTION double michelSpectrum(double x) {
        return 2.0 * x * x * (3.0 - 2.0 * x);
    }

    /// Upper bound of the Michel spectrum, attained at \f$ x = 1 \f$: \f$ f(1) = 2 \f$.
    KOKKOS_INLINE_FUNCTION constexpr double michelUpperBound() {
        return 2.0;
    }

}  // namespace MuonDecay
}  // namespace Physics

#endif  // OPAL_PHYSICS_MUON_DECAY_H
