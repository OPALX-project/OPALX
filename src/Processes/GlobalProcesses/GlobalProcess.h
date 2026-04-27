#ifndef OPAL_GLOBALPROCESS_H
#define OPAL_GLOBALPROCESS_H

#include <cstddef>

template <typename T, unsigned Dim>
class ParticleContainer;

class GlobalProcess {
public:
    virtual ~GlobalProcess() = default;

    virtual size_t apply(
            ParticleContainer<double, 3>& pc, double dt, long long globalTrackStep,
            size_t containerIdx) = 0;
};

#endif
