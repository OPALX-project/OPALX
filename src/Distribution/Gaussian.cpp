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
 * @param opalDist Borrowed distribution object.
 */
Gaussian::Gaussian(std::shared_ptr<ParticleContainer_t> pc, Distribution_t* opalDist)
    : SamplingBase(pc, opalDist) {
    samperTimer_m = IpplTimings::getTimer("Sampling");
    initRandomPool();
    setSigmaP(opalDist->getSigmaP());
    setSigmaR(opalDist->getSigmaR());
    setAvrgpz(opalDist->getAvrgpz());
    setCutoffR(opalDist->getCutoffR());
    setCutoffP(opalDist->getCutoffP());
}

Gaussian::Gaussian(
        std::shared_ptr<ParticleContainer_t> pc, const Vector_t<double, 3>& sigmaR,
        const Vector_t<double, 3>& sigmaP, double avrgpz, const Vector_t<double, 3>& cutoffR,
        bool fixMeanR, const Vector_t<double, 3>& cutoffP)
    : SamplingBase(pc) {
    setSigmaR(sigmaR);
    setSigmaP(sigmaP);
    setAvrgpz(avrgpz);
    setCutoffR(cutoffR);
    setCutoffP(cutoffP);
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
 * The algorithm follows OPAL's Distribution::generateGaussZ (rejection sampling
 * with elliptical transverse cuts) and Distribution::adjustPhaseSpace (exact
 * centering and average-energy matching).
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

    // The sampler expects strictly positive cutoffs (in units of sigma). "No cut" is
    // expressed as a very large value, e.g. std::numeric_limits<double>::max(), which
    // is what Distribution::setDistParametersGauss substitutes for a zero input cutoff.
    for (unsigned d = 0; d < 3; ++d) {
        if (!(cutoffR_m[d] > 0.0) || !(cutoffP_m[d] > 0.0)) {
            throw OpalException(
                    "Gaussian::generateParticles",
                    "All position and momentum cutoffs (CUTOFFX, CUTOFFY, CUTOFFLONG, CUTOFFPX, "
                    "CUTOFFPY, CUTOFFPZ) must be greater than zero.");
        }
    }

    auto rand_pool64 = randPool_m;

    IpplTimings::startTimer(samperTimer_m);

    const int nranks = std::max(1, ippl::Comm->size());
    // Use computeLocalEmitCount to distribute particles across ranks or uniform fallback.
    size_t nlocal = pc_m ? computeLocalEmitCount(static_cast<size_t>(numberOfParticles))
                         : static_cast<size_t>(
                                   floor(numberOfParticles / nranks)
                                   + (ippl::Comm->rank() < static_cast<int>(
                                              numberOfParticles % static_cast<size_t>(nranks))
                                              ? 1
                                              : 0));

    const size_t nlocalCurrent = pc_m->getLocalNum();
    pc_m->createParticles(nlocal);

    view_type RviewFull = pc_m->R.getView();
    view_type PviewFull = pc_m->P.getView();
    auto Rview = Kokkos::subview(RviewFull, std::make_pair(nlocalCurrent, nlocalCurrent + nlocal));
    auto Pview = Kokkos::subview(PviewFull, std::make_pair(nlocalCurrent, nlocalCurrent + nlocal));

    // Rejection sampling in normalised (unit-sigma) coordinates, as in OPAL:
    //   (x/cx)^2 + (y/cy)^2 <= 1,  (px/cpx)^2 + (py/cpy)^2 <= 1,  |z| <= cz,  |pz| <= cpz.
    // A zero sigma collapses the corresponding coordinate to exactly zero (point mass).
    const Vector_t<double, 3> sigmaR  = sigmaR_m;
    const Vector_t<double, 3> sigmaP  = sigmaP_m;
    const Vector_t<double, 3> cutoffR = cutoffR_m;
    const Vector_t<double, 3> cutoffP = cutoffP_m;
    const double avrgpz               = avrgpz_m;

    Kokkos::parallel_for(
            "Gaussian_rejection_sampling", nlocal, KOKKOS_LAMBDA(const size_t k) {
                typename GeneratorPool::generator_type gen = rand_pool64.get_state();

                double x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0;
                bool allow = false;
                while (!allow) {
                    x  = gen.normal(0.0, 1.0);
                    px = gen.normal(0.0, 1.0);
                    y  = gen.normal(0.0, 1.0);
                    py = gen.normal(0.0, 1.0);
                    z  = gen.normal(0.0, 1.0);
                    pz = gen.normal(0.0, 1.0);

                    const double rx  = x / cutoffR[0];
                    const double ry  = y / cutoffR[1];
                    const double rpx = px / cutoffP[0];
                    const double rpy = py / cutoffP[1];

                    const bool xAndYOk   = (rx * rx + ry * ry <= 1.0);
                    const bool pxAndPyOk = (rpx * rpx + rpy * rpy <= 1.0);
                    const bool zOk       = (Kokkos::fabs(z) <= cutoffR[2]);
                    const bool pzOk      = (Kokkos::fabs(pz) <= cutoffP[2]);

                    allow = (xAndYOk && pxAndPyOk && zOk && pzOk);
                }

                rand_pool64.free_state(gen);

                Rview(k)[0] = x * sigmaR[0];
                Rview(k)[1] = y * sigmaR[1];
                Rview(k)[2] = z * sigmaR[2];
                Pview(k)[0] = px * sigmaP[0];
                Pview(k)[1] = py * sigmaP[1];
                Pview(k)[2] = avrgpz + pz * sigmaP[2];
            });
    Kokkos::fence();

    if (fixMeanR_m) {
        // Equivalent of OPAL's Distribution::adjustPhaseSpace:
        //  - shift x, y, z, px, py such that their means are exactly zero,
        //  - shift pz by eps such that  mean( pz + (px^2 + py^2) / (2 pz) ) == avrgpz,
        //    i.e. the mean of |p| (Taylor expanded for small transverse momenta)
        //    matches the reference beam momentum exactly.
        double loc_avrg[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double avrg[6]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        Kokkos::parallel_reduce(
                "Gaussian_phase_space_moments", nlocal,
                KOKKOS_LAMBDA(
                        const size_t k, double& sx, double& sy, double& sz, double& spx,
                        double& spy, double& sp) {
                    sx += Rview(k)[0];
                    sy += Rview(k)[1];
                    sz += Rview(k)[2];
                    spx += Pview(k)[0];
                    spy += Pview(k)[1];
                    const double px = Pview(k)[0];
                    const double py = Pview(k)[1];
                    const double pz = Pview(k)[2];
                    // Taylor series of sqrt(px^2 + py^2 + pz^2) = pz * sqrt(1 + eps^2), eps << 1.
                    // Guard pz == 0 (e.g. avrgpz == 0 with zero momentum spread).
                    sp += (pz != 0.0) ? pz + (px * px + py * py) / (2.0 * pz) : pz;
                },
                Kokkos::Sum<double>(loc_avrg[0]), Kokkos::Sum<double>(loc_avrg[1]),
                Kokkos::Sum<double>(loc_avrg[2]), Kokkos::Sum<double>(loc_avrg[3]),
                Kokkos::Sum<double>(loc_avrg[4]), Kokkos::Sum<double>(loc_avrg[5]));
        Kokkos::fence();

        MPI_Allreduce(loc_avrg, avrg, 6, MPI_DOUBLE, MPI_SUM, ippl::Comm->getCommunicator());
        ippl::Comm->barrier();

        for (int i = 0; i < 6; i++) {
            avrg[i] = avrg[i] / (1.0 * numberOfParticles);
        }

        const double eps = avrgpz - avrg[5];

        Kokkos::parallel_for(
                "Gaussian_adjust_phase_space", nlocal, KOKKOS_LAMBDA(const size_t k) {
                    Rview(k)[0] -= avrg[0];
                    Rview(k)[1] -= avrg[1];
                    Rview(k)[2] -= avrg[2];
                    Pview(k)[0] -= avrg[3];
                    Pview(k)[1] -= avrg[4];
                    Pview(k)[2] += eps;
                });
        Kokkos::fence();
    }

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
    Kokkos::fence();

    fillPolarization(nlocalCurrent, nlocal);

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
