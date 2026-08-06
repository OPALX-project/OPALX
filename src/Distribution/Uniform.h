#ifndef OPALX_UNIFORM_H
#define OPALX_UNIFORM_H

#include <Kokkos_Random.hpp>

#include "Distribution/SamplingBase.hpp"

/**
 * One-shot cold distribution that is spatially uniform inside an ellipsoid.
 * SIGMAX, SIGMAY and SIGMAZ are interpreted as its semi-axes.
 */
class Uniform : public SamplingBase {
public:
    Uniform(std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            Distribution_t* opalDist);
    Uniform(std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& semiAxes,
            double avrgpz = 0.0);

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

private:
    void initRandomPool();

    Kokkos::Random_XorShift64_Pool<> randPool_m;
    Vector_t<double, 3> semiAxes_m;
    double avrgpz_m;
};

#endif  // OPALX_UNIFORM_H
