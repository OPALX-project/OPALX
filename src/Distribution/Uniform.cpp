#include "Distribution/Uniform.h"

#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>

Uniform::Uniform(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        Distribution_t* opalDist)
    : SamplingBase(pc, fc, opalDist),
      semiAxes_m(opalDist->getSigmaR()),
      avrgpz_m(opalDist->getAvrgpz()) {
    initRandomPool();
}

Uniform::Uniform(
        std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& semiAxes, double avrgpz)
    : SamplingBase(pc), semiAxes_m(semiAxes), avrgpz_m(avrgpz) {
    initRandomPool();
}

void Uniform::initRandomPool() {
    const size_t seed = Options::seed == -1
                                ? 1234567
                                : static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
    randPool_m        = Kokkos::Random_XorShift64_Pool<>(seed);
}

void Uniform::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> /*nr*/) {
    if (emissionModel_m != "NONE") {
        throw OpalException(
                "Uniform::generateParticles",
                "EMISSIONMODEL '" + emissionModel_m + "' is not supported for UNIFORM");
    }
    if (t0_m > 0.0) {
        throw OpalException(
                "Uniform::generateParticles", "Delayed emission is not supported for UNIFORM");
    }

    const size_t nlocal = computeLocalEmitCount(numberOfParticles);
    const size_t offset = pc_m->getLocalNum();
    pc_m->createParticles(nlocal);

    auto R                 = pc_m->R.getView();
    auto P                 = pc_m->P.getView();
    auto pool              = randPool_m;
    const auto a           = semiAxes_m;
    const auto R0          = R0_m;
    const auto P0          = P0_m;
    const double pz        = avrgpz_m;
    constexpr double twoPi = 6.283185307179586476925286766559;

    Kokkos::parallel_for(
            "Uniform::generateParticles", nlocal, KOKKOS_LAMBDA(const size_t j) {
                auto generator      = pool.get_state();
                const double radius = Kokkos::pow(generator.drand(), 1.0 / 3.0);
                const double z      = 2.0 * generator.drand() - 1.0;
                const double phi    = twoPi * generator.drand();
                const double rxy    = radius * Kokkos::sqrt(1.0 - z * z);
                const size_t i      = offset + j;

                R(i) = Vector_t<double, 3>(
                               a[0] * rxy * Kokkos::cos(phi), a[1] * rxy * Kokkos::sin(phi),
                               a[2] * radius * z)
                       + R0;
                P(i) = Vector_t<double, 3>(0.0, 0.0, pz) + P0;
                pool.free_state(generator);
            });
    Kokkos::fence();

    fillPolarization(offset, nlocal);
    pc_m->markMomentsDirty();
    hasEmittedOnce_m = true;
}
