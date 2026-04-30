#include "Gaussian.h"
#include <algorithm>
#include <memory>
#include "Distribution.h"
#include "SamplingBase.hpp"
#include "Utilities/OpalException.h"

/**
 * @brief Constructs a Gaussian sampler.
 *
 * @param pc Shared pointer to the particle container.
 * @param fc Shared pointer to the field container.
 * @param opalDist Borrowed distribution object.
 */
Gaussian::Gaussian(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        Distribution_t* opalDist)
    : SamplingBase(pc, fc, opalDist) {
    samperTimer_m = IpplTimings::getTimer("Sampling");
    initRandomPool();
    setSigmaP(opalDist->getSigmaP());
    setSigmaR(opalDist->getSigmaR());
    setAvrgpz(opalDist->getAvrgpz());
    setCutoffR(opalDist->getCutoffR());
}

Gaussian::Gaussian(
        std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& sigmaR,
        const Vector_t<double, 3>& sigmaP, double avrgpz, const Vector_t<double, 3>& cutoffR,
        bool fixMeanR)
    : SamplingBase(pc) {
    setSigmaR(sigmaR);
    setSigmaP(sigmaP);
    setAvrgpz(avrgpz);
    setCutoffR(cutoffR);
    setFixMeanR(fixMeanR);

    samperTimer_m = IpplTimings::getTimer("Sampling");
    initRandomPool();
}

void Gaussian::initRandomPool() {
    extern Inform* gmsg;
    size_t randInit;

    if (Options::seed == -1) {
        randInit = 1234567;
        *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
    } else {
        randInit = static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
    }

    GeneratorPool rand_pool64(randInit);
    randPool_m = rand_pool64;
    return;
}

/**
 * @brief Generates particles following a Gaussian distribution.
 *
 * @param numberOfParticles The total number of particles to generate.
 * @param nr The number of grid points in each dimension (not used here).
 */
void Gaussian::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> /*nr*/) {
    if (emissionModel_m != "NONE")
        throw OpalException(
                "Gaussian::generateParticles",
                "EMISSIONMODEL '" + emissionModel_m + "' is not supported for GAUSS distributions");

    // Only generate during initial sampling (t0 <= 0). For t0 > 0, this
    // distribution is time-independent and should not contribute here unless
    // explicitly triggered via emitParticles (which sets hasEmittedOnce_m).
    if (t0_m > 0.0 && !hasEmittedOnce_m) {  // YES this !hasEmittedOnce_m is correct!
        return;
    }
    auto rand_pool64 = randPool_m;

    IpplTimings::startTimer(samperTimer_m);

    double mu[3], sd[3];

    Vector_t<double, 3> rmin = -cutoffR_m;
    Vector_t<double, 3> rmax = cutoffR_m;

    for (int i = 0; i < 3; i++) {
        mu[i]   = 0.0;
        sd[i]   = sigmaR_m[i];
        rmin(i) = (rmin(i) + mu[i]) * sigmaR_m[i];
        rmax(i) = (rmax(i) + mu[i]) * sigmaR_m[i];
    }

    const double par[6] = {mu[0], sd[0], mu[1], sd[1], mu[2], sd[2]};

    using Dist_t     = ippl::random::NormalDistribution<double, 3>;
    using sampling_t = ippl::random::InverseTransformSampling<
            double, 3, Kokkos::DefaultExecutionSpace, Dist_t>;
    Dist_t dist(par);

    const int nranks = std::max(1, ippl::Comm->size());
    // Use computeLocalEmitCount to distribute particles across ranks or uniform fallback.
    size_t nlocal = pc_m ? computeLocalEmitCount(static_cast<size_t>(numberOfParticles))
                         : static_cast<size_t>(
                                   floor(numberOfParticles / nranks)
                                   + (ippl::Comm->rank() < static_cast<int>(
                                              numberOfParticles % static_cast<size_t>(nranks))
                                              ? 1
                                              : 0));

    sampling_t sampling(dist, rmax, rmin, rmax, rmin, nlocal);
    nlocal = sampling.getLocalSamplesNum();

    const size_t nlocalCurrent = pc_m->getLocalNum();
    pc_m->createParticles(nlocal);

    view_type RviewFull = pc_m->R.getView();
    auto Rview = Kokkos::subview(RviewFull, std::make_pair(nlocalCurrent, nlocalCurrent + nlocal));

    sampling.generate(Rview, rand_pool64);

    if (fixMeanR_m) {
        double meanR[3], loc_meanR[3];
        for (int i = 0; i < 3; i++) {
            meanR[i]     = 0.0;
            loc_meanR[i] = 0.0;
        }

        Kokkos::parallel_reduce(
                "calc moments of particle distr.", nlocal,
                KOKKOS_LAMBDA(const int k, double& cent0, double& cent1, double& cent2) {
                    cent0 += Rview(k)[0];
                    cent1 += Rview(k)[1];
                    cent2 += Rview(k)[2];
                },
                Kokkos::Sum<double>(loc_meanR[0]), Kokkos::Sum<double>(loc_meanR[1]),
                Kokkos::Sum<double>(loc_meanR[2]));
        Kokkos::fence();

        MPI_Allreduce(loc_meanR, meanR, 3, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        ippl::Comm->barrier();

        for (int i = 0; i < 3; i++) {
            meanR[i] = meanR[i] / (1.0 * numberOfParticles);
        }

        Kokkos::parallel_for(
                nlocal, KOKKOS_LAMBDA(const int k) {
                    Rview(k)[0] -= meanR[0];
                    Rview(k)[1] -= meanR[1];
                    Rview(k)[2] -= meanR[2];
                });
        Kokkos::fence();
    }

    for (int i = 0; i < 3; i++) {
        mu[i] = 0.0;
        sd[i] = sigmaP_m[i];
    }

    view_type PviewFull = pc_m->P.getView();
    auto Pview = Kokkos::subview(PviewFull, std::make_pair(nlocalCurrent, nlocalCurrent + nlocal));

    Kokkos::parallel_for(nlocal, ippl::random::randn<double, 3>(Pview, rand_pool64, mu, sd));
    Kokkos::fence();

    double avrgpz = avrgpz_m;
    Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(const size_t k) { Pview(k)[2] += avrgpz; });
    Kokkos::fence();

    // Apply per-emission-source offsets after all mean-fixing/corrections.
    // EMISSIONSOURCE offsets are expected to translate the generated bunch
    // without being affected by the internal "fix mean" logic.
    const Vector_t<double, 3> R0 = R0_m;
    const Vector_t<double, 3> P0 = P0_m;
    Kokkos::parallel_for(
            nlocal, KOKKOS_LAMBDA(const size_t k) {
                Rview(k) += R0;
                Pview(k) += P0;
            });

    pc_m->markMomentsDirty();

    IpplTimings::stopTimer(samperTimer_m);
}

void Gaussian::emitParticles(double t, double dt) {
    // One-shot delayed emission for GAUSS: emit once when [t, t+dt] crosses t0_m.
    const double tStart = t;
    const double tEnd   = t + dt;

    if (hasEmittedOnce_m) {
        return;
    }

    // Only meaningful for t0 > 0.
    if (t0_m <= 0.0) {
        return;
    }

    // Fire when the time interval [tStart, tEnd] crosses t0_m.
    if (!(tStart <= t0_m && t0_m < tEnd)) {
        return;
    }

    if (!opalDist_m) {
        return;
    }
    size_t Ndist = opalDist_m->getNumParticles();
    if (Ndist == 0) {
        return;
    }

    const size_t nlocalBefore = pc_m->getLocalNum();

    // Mark as emitted so generateParticles will not early-return.
    hasEmittedOnce_m = true;
    Vector_t<double, 3> dummyNr(0.0);
    generateParticles(Ndist, dummyNr);

    // Set per-particle dt for newly created particles. Note that t0 can be in
    // between particles, so we set the remaining fraction of the current time
    // step for the new particles to ensure correct sub-stepping in the
    // subsequent push.
    //
    // switchToUnitlessPositions(use_dt_per_particle=true) in pushParticles scales
    // R by 1/(c * dtview(i)).  New particles come with dtview = 0 (Kokkos
    // zero-init), so that division produces inf, and the subsequent multiply by
    // c*0 in switchOffUnitlessPositions gives inf*0 = NaN positions.
    // Assign the remaining fraction of the current timestep so the push in
    // timeIntegration2 moves them by the correct sub-step distance.
    //
    // Note that this is not necessary inside generateParticles, since for injection at start,
    // setTime is called before any calculations inside ParallelTracker::execute.
    const size_t nlocalAfter = pc_m->getLocalNum();
    const size_t nNew        = nlocalAfter - nlocalBefore;
    if (nNew > 0) {
        const double fracDt = tEnd - t0_m;
        auto dtview         = pc_m->dt.getView();
        const size_t offset = nlocalBefore;
        Kokkos::parallel_for(
                "Gaussian_setDt", nNew,
                KOKKOS_LAMBDA(const size_t j) { dtview(offset + j) = fracDt; });
        Kokkos::fence();
    }
}
