#ifndef OPAL_GLOBALPROCESS_DECAY_H
#define OPAL_GLOBALPROCESS_DECAY_H

#include "Processes/GlobalProcesses/GlobalProcess.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <cstddef>
#include <memory>

namespace ippl {
    template <typename T, unsigned Dim>
    class Vector;
}

class Decay : public GlobalProcess {
public:
    Decay(double restLifetimeSeconds, std::size_t containerIndex, double parentMassGeV);

    ~Decay() override = default;

    /// Set the daughter particle container and its rest mass.
    /// When set, decayed particles spawn a daughter in @p daughterPC.
    /// When not set (default), decayed particles are simply destroyed.
    void setDaughterContainer(std::shared_ptr<ParticleContainer<double, 3>> daughterPC,
                              double daughterMassGeV);

    size_t apply(ParticleContainer<double, 3>& pc,
                 double dt,
                 long long globalTrackStep,
                 size_t containerIdx) override;

protected:

    /// Create daughter particles from collected parent data.
    /// Called only when daughterPC_m is set and localDestroyNum > 0.
    /// Subclasses implement the physics-specific momentum sampling here.
    virtual void createDaughterParticles(
        std::size_t localDestroyNum,
        std::size_t oldDaughterLocal,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentR,
        const Kokkos::View<ippl::Vector<double, 3>*>& parentP,
        const Kokkos::View<double*>& parentDt) = 0;

    /// Mean life time at rest [s]
    double tau0_m;

    /// Random pool for decay sampling
    Kokkos::Random_XorShift64_Pool<> randPool_m;

    /// Daughter container for decay products (nullptr = destroy-only mode).
    std::shared_ptr<ParticleContainer<double, 3>> daughterPC_m;

    /// Rest mass of the daughter particle [GeV].
    double daughterMassGeV_m = 0.0;

    /// Rest mass of the parent particle [GeV].
    double parentMassGeV_m;
};

#endif
