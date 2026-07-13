#ifndef OPALX_PARTBUNCH_FWD_H
#define OPALX_PARTBUNCH_FWD_H

template <typename T, unsigned Dim>
class PartBunch;

template <typename T, unsigned Dim>
class ParticleContainer;

using PartBunch_t         = PartBunch<double, 3>;
using ParticleContainer_t = ParticleContainer<double, 3>;

#endif
