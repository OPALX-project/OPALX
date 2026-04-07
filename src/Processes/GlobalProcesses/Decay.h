#ifndef OPAL_GLOBALPROCESS_DECAY_H
#define OPAL_GLOBALPROCESS_DECAY_H

#include "Processes/GlobalProcesses/GlobalProcess.h"

#include <Kokkos_Random.hpp>

#include <cstddef>

class Decay : public GlobalProcess {
public:
    Decay(double restLifetimeSeconds, std::size_t containerIndex);

    size_t apply(ParticleContainer<double, 3>& pc,
                 double dt,
                 long long globalTrackStep,
                 size_t containerIdx) override;

private:
    double tau0_m;
    Kokkos::Random_XorShift64_Pool<> randPool_m;
};

#endif
