#ifndef OPAL_GLOBALPROCESS_DECAY_H
#define OPAL_GLOBALPROCESS_DECAY_H

#include "Processes/GlobalProcesses/GlobalProcess.h"

#include <Kokkos_Random.hpp>

#include <cstddef>
#include <memory>

class Decay : public GlobalProcess {
public:
    Decay(double restLifetimeSeconds, std::size_t containerIndex);

    /// Set the daughter particle container and its rest mass.
    /// When set, decayed particles spawn a daughter in @p daughterPC.
    /// When not set (default), decayed particles are simply destroyed.
    void setDaughterContainer(std::shared_ptr<ParticleContainer<double, 3>> daughterPC,
                              double daughterMassGeV);

    size_t apply(ParticleContainer<double, 3>& pc,
                 double dt,
                 long long globalTrackStep,
                 size_t containerIdx) override;

private:
    double tau0_m;
    Kokkos::Random_XorShift64_Pool<> randPool_m;

    /// Daughter container for decay products (nullptr = destroy-only mode).
    std::shared_ptr<ParticleContainer<double, 3>> daughterPC_m;

    /// Rest mass of the daughter particle [GeV].
    double daughterMassGeV_m = 0.0;
};

#endif
