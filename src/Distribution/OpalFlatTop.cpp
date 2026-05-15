#include "OpalFlatTop.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>

#include "Distribution.h"
#include "Physics/Physics.h"
#include "SamplingBase.hpp"
#include "Utilities/OpalException.h"

OpalFlatTop::OpalFlatTop(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        Distribution_t* opalDist)
    : SamplingBase(pc, fc, opalDist),
      rand_pool_m(determineRandInit()),
      host_rng_m(determineHostSeed()) {
    setParameters(opalDist);
}

OpalFlatTop::OpalFlatTop(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
        double tPulseLengthFWHM, Vector_t<double, 3> sigmaR, double ftOscAmplitude,
        double ftOscPeriods)
    : SamplingBase(pc, fc), rand_pool_m(determineRandInit()), host_rng_m(determineHostSeed()) {
    setInternalVariables(
            emitting, sigmaTFall, sigmaTRise, cutoff, tPulseLengthFWHM, sigmaR, ftOscAmplitude,
            ftOscPeriods);
}

size_t OpalFlatTop::determineRandInit() {
    extern Inform* gmsg;
    size_t randInit;
    if (Options::seed == -1) {
        randInit = 1234567;
        *gmsg << "* Seed = " << randInit << " on all ranks" << endl;
    } else {
        randInit = static_cast<size_t>(Options::seed + 100 * ippl::Comm->rank());
    }
    return randInit;
}

size_t OpalFlatTop::determineHostSeed() {
    if (Options::seed == -1) {
        return 1234567;
    }
    return static_cast<size_t>(Options::seed);
}

void OpalFlatTop::setParameters(Distribution_t* opalDist) {
    setInternalVariables(
            opalDist->emitting_m, opalDist->getSigmaTFall(), opalDist->getSigmaTRise(),
            opalDist->getCutoffR(), opalDist->getTPulseLengthFWHM(), opalDist->getSigmaR(),
            opalDist->getFTOSCAmplitude(), opalDist->getFTOSCPeriods());

    emissionSteps_m = opalDist->getEmissionSteps();
    opalDist->setTEmission(emissionTime_m);

    if (fc_m) {
        fc_m->setDecomp({false, false, true});
    }
}

void OpalFlatTop::setInternalVariables(
        bool emitting, double sigmaTFall, double sigmaTRise, Vector_t<double, 3> cutoff,
        double tPulseLengthFWHM, Vector_t<double, 3> sigmaR, double ftOscAmplitude,
        double ftOscPeriods) {
    emitting_m         = emitting;
    sigmaTFall_m       = sigmaTFall;
    sigmaTRise_m       = sigmaTRise;
    cutoffR_m          = cutoff;
    sigmaR_m           = sigmaR;
    ftOscAmplitude_m   = ftOscAmplitude;
    ftOscPeriods_m     = ftOscPeriods;
    withDomainDecomp_m = false;

    fallTime_m = sigmaTFall_m * cutoffR_m[2];
    flattopTime_m =
            tPulseLengthFWHM - std::sqrt(2.0 * std::log(2.0)) * (sigmaTRise_m + sigmaTFall_m);
    if (flattopTime_m < 0.0) {
        flattopTime_m = 0.0;
    }
    riseTime_m     = sigmaTRise_m * cutoffR_m[2];
    emissionTime_m = fallTime_m + flattopTime_m + riseTime_m;

    normalizedFlankArea_m =
            0.5 * std::sqrt(Physics::two_pi) * std::erf(cutoffR_m[2] / std::sqrt(2.0));
    distArea_m = flattopTime_m + (sigmaTRise_m + sigmaTFall_m) * normalizedFlankArea_m;
}

void OpalFlatTop::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) {
    nr_m = nr;

    if (!emitting_m) {
        throw OpalException(
                "OpalFlatTop::generateParticles",
                "OPALFLATTOP only supports emitted beams. Set EMITTED = true.");
    }
    if (numberOfParticles == 0) {
        throw OpalException(
                "OpalFlatTop::generateParticles",
                "OPALFLATTOP requires a positive NPARTDIST value.");
    }

    totalN_m = numberOfParticles;
    buildBirthTimeInventory(numberOfParticles);
}

void OpalFlatTop::buildBirthTimeInventory(size_t numberOfParticles) {
    if (distArea_m <= 0.0) {
        throw OpalException(
                "OpalFlatTop::buildBirthTimeInventory",
                "The flat-top distribution area is zero. Check TPULSEFWHM, TRISE, TFALL, and "
                "CUTOFFLONG.");
    }

    birthTimes_m.clear();
    birthTimes_m.reserve(numberOfParticles);
    nextGlobalIndex_m = 0;
    inventoryBuilt_m  = false;

    const size_t numRise = static_cast<size_t>(
            static_cast<double>(numberOfParticles) * sigmaTRise_m * normalizedFlankArea_m
            / distArea_m);
    const size_t numFall = static_cast<size_t>(
            static_cast<double>(numberOfParticles) * sigmaTFall_m * normalizedFlankArea_m
            / distArea_m);
    const size_t numFlat =
            (numRise + numFall <= numberOfParticles) ? numberOfParticles - numRise - numFall : 0;

    const double fallLimit = sigmaTFall_m * cutoffR_m[2];
    for (size_t i = 0; i < numFall; ++i) {
        const double tailTime      = sampleTruncatedHalfGaussian(sigmaTFall_m, fallLimit);
        const double opalPulseTime = -tailTime + fallTime_m;
        birthTimes_m.push_back(toBirthTime(opalPulseTime));
    }

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double modulationAmp = std::min(1.0, std::abs(ftOscAmplitude_m) / 100.0);
    const double numPeriods    = std::abs(ftOscPeriods_m);
    const double modulationPeriod =
            (numPeriods > 0.0 && flattopTime_m > 0.0) ? flattopTime_m / numPeriods : 0.0;

    for (size_t i = 0; i < numFlat; ++i) {
        double tFlat = 0.0;
        if (flattopTime_m > 0.0) {
            if (modulationAmp == 0.0 || modulationPeriod == 0.0) {
                tFlat = flattopTime_m * unit(host_rng_m);
            } else {
                bool allow = false;
                while (!allow) {
                    tFlat                      = flattopTime_m * unit(host_rng_m);
                    const double randFuncValue = unit(host_rng_m);
                    const double funcValue =
                            (1.0
                             + modulationAmp * std::sin(Physics::two_pi * tFlat / modulationPeriod))
                            / (1.0 + modulationAmp);
                    allow = randFuncValue <= funcValue;
                }
            }
        }
        birthTimes_m.push_back(toBirthTime(fallTime_m + tFlat));
    }

    const double riseLimit = sigmaTRise_m * cutoffR_m[2];
    for (size_t i = 0; i < numRise; ++i) {
        const double tailTime      = sampleTruncatedHalfGaussian(sigmaTRise_m, riseLimit);
        const double opalPulseTime = tailTime + fallTime_m + flattopTime_m;
        birthTimes_m.push_back(toBirthTime(opalPulseTime));
    }

    while (birthTimes_m.size() < numberOfParticles) {
        birthTimes_m.push_back(toBirthTime(fallTime_m));
    }

    std::sort(birthTimes_m.begin(), birthTimes_m.end());
    inventoryBuilt_m = true;
}

double OpalFlatTop::sampleTruncatedHalfGaussian(double sigma, double limit) {
    if (sigma <= 0.0 || limit <= 0.0) {
        return 0.0;
    }

    std::normal_distribution<double> normal(0.0, sigma);
    for (;;) {
        const double value = std::abs(normal(host_rng_m));
        if (value <= limit) {
            return value;
        }
    }
}

double OpalFlatTop::toBirthTime(double opalPulseTime) const {
    return t0_m + 0.5 * emissionTime_m - opalPulseTime;
}

Vector_t<double, 3> OpalFlatTop::getInitialReferenceMomentum() const {
    if (emissionModel_m == "ASTRA") {
        return Vector_t<double, 3>({0.0, 0.0, 0.5 * euclidean_norm(P0_m)});
    }
    return P0_m;
}

std::pair<size_t, size_t> OpalFlatTop::computeLocalEmitRange(size_t totalToEmit) const {
    if (!pc_m || totalToEmit == 0) {
        return {0, 0};
    }

    const int nranks = ippl::Comm->size();
    const int rank   = ippl::Comm->rank();
    if (nranks <= 0) {
        return {0, totalToEmit};
    }

    const size_t capacity  = pc_m->R.size();
    const size_t localNum  = pc_m->getLocalNum();
    const size_t spaceLeft = (localNum <= capacity) ? (capacity - localNum) : size_t(0);

    std::vector<unsigned long> spaceLeftAll(static_cast<size_t>(nranks), 0);
    unsigned long mySpace = static_cast<unsigned long>(spaceLeft);
    MPI_Allgather(
            &mySpace, 1, MPI_UNSIGNED_LONG, spaceLeftAll.data(), 1, MPI_UNSIGNED_LONG,
            ippl::Comm->getCommunicator());

    const size_t nranksU = static_cast<size_t>(nranks);
    const size_t base    = totalToEmit / nranksU;
    const size_t rem     = totalToEmit % nranksU;

    std::vector<size_t> nlocal(nranksU, 0);
    size_t assigned = 0;
    for (int r = 0; r < nranks; ++r) {
        const size_t ideal             = base + (static_cast<size_t>(r) < rem ? 1 : 0);
        const size_t cap               = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
        nlocal[static_cast<size_t>(r)] = std::min(ideal, cap);
        assigned += nlocal[static_cast<size_t>(r)];
    }

    size_t deficit = totalToEmit - assigned;
    while (deficit > 0) {
        int chosen = -1;
        for (int r = 0; r < nranks; ++r) {
            const size_t cap = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
            if (nlocal[static_cast<size_t>(r)] < cap) {
                chosen = r;
                break;
            }
        }
        if (chosen < 0) {
            throw OpalException(
                    "OpalFlatTop::computeLocalEmitRange",
                    "Particle container capacity is insufficient for OPALFLATTOP emission. "
                    "Increase BEAM::NALLOC or reduce NPARTDIST.");
        }
        nlocal[static_cast<size_t>(chosen)] += 1;
        --deficit;
    }

    size_t offset = 0;
    for (int r = 0; r < rank; ++r) {
        offset += nlocal[static_cast<size_t>(r)];
    }
    return {offset, nlocal[static_cast<size_t>(rank)]};
}

void OpalFlatTop::emitParticles(double t, double dt) {
    if (!inventoryBuilt_m || birthTimes_m.empty()) {
        return;
    }
    if (nextGlobalIndex_m >= birthTimes_m.size()) {
        return;
    }

    const double tEnd = t + dt;
    auto emitEndIt =
            std::upper_bound(birthTimes_m.begin() + nextGlobalIndex_m, birthTimes_m.end(), tEnd);
    const size_t globalEnd = static_cast<size_t>(emitEndIt - birthTimes_m.begin());
    const size_t totalNew  = globalEnd - nextGlobalIndex_m;

    if (totalNew == 0) {
        return;
    }

    const double overdueTolerance = std::max(1.0e-18, std::abs(dt) * 1.0e-12);
    if (birthTimes_m[nextGlobalIndex_m] < t - overdueTolerance) {
        throw OpalException(
                "OpalFlatTop::emitParticles",
                "OPALFLATTOP found an overdue birth time. This means the tracker "
                "time shift or emission dt skipped particles that should have been "
                "emitted in an earlier step.");
    }

    const auto [localOffset, nNew] = computeLocalEmitRange(totalNew);
    const size_type nlocalBefore   = pc_m->getLocalNum();
    pc_m->createParticles(static_cast<size_type>(nNew));

    generateLocalParticles(nlocalBefore, nextGlobalIndex_m + localOffset, nNew, t, dt);
    nextGlobalIndex_m = globalEnd;
}

void OpalFlatTop::generateLocalParticles(
        size_type nlocalBefore, size_t globalBegin, size_t nNew, double tStart, double dt) {
    if (nNew == 0) {
        return;
    }

    Kokkos::View<double*> stepDt("OpalFlatTop_stepDt", nNew);
    auto hStepDt                  = Kokkos::create_mirror_view(stepDt);
    const double tEnd             = tStart + dt;
    const double overdueTolerance = std::max(1.0e-18, std::abs(dt) * 1.0e-12);
    for (size_t i = 0; i < nNew; ++i) {
        const double birthTime = birthTimes_m[globalBegin + i];
        if (birthTime < tStart - overdueTolerance) {
            throw OpalException(
                    "OpalFlatTop::generateLocalParticles",
                    "OPALFLATTOP would need to pre-age a particle born before the "
                    "current step. Old OPAL emits these particles in earlier tracker "
                    "steps, so this indicates an inconsistent tracker time shift.");
        }
        const double effectiveBirthTime = birthTime < tStart ? tStart : birthTime;
        hStepDt(i)                      = tEnd - effectiveBirthTime;
    }
    Kokkos::deep_copy(stepDt, hStepDt);

    auto randPool = rand_pool_m;
    auto Rview    = pc_m->R.getView();
    auto Pview    = pc_m->P.getView();
    auto dtView   = pc_m->dt.getView();

    const Vector_t<double, 3> sigmaR = sigmaR_m;
    const Vector_t<double, 3> R0     = R0_m;
    const Vector_t<double, 3> P0     = P0_m;
    const double pTot                = euclidean_norm(P0);
    const bool useAstra              = emissionModel_m == "ASTRA";
    const bool useNone               = emissionModel_m == "NONE";
    if (!useAstra && !useNone) {
        throw OpalException(
                "OpalFlatTop::generateLocalParticles",
                "EMISSIONMODEL '" + emissionModel_m
                        + "' is not supported for OPALFLATTOP distributions");
    }

    const size_type offset = nlocalBefore;
    const double twoPi     = Physics::two_pi;
    const double c         = Physics::c;
    Kokkos::parallel_for(
            "OpalFlatTop_generateLocalParticles", nNew, KOKKOS_LAMBDA(const size_t i) {
                auto generator     = randPool.get_state();
                const double r     = Kokkos::sqrt(generator.drand(0.0, 1.0));
                const double theta = twoPi * generator.drand(0.0, 1.0);

                Vector_t<double, 3> p = P0;
                if (useAstra) {
                    const double rand1 = generator.drand(0.0, 1.0);
                    const double rand2 = generator.drand(0.0, 1.0);
                    const double phi   = 2.0 * Kokkos::acos(Kokkos::sqrt(rand1));
                    const double angle = twoPi * rand2;
                    p[0]               = pTot * Kokkos::sin(phi) * Kokkos::cos(angle);
                    p[1]               = pTot * Kokkos::sin(phi) * Kokkos::sin(angle);
                    p[2]               = pTot * Kokkos::fabs(Kokkos::cos(phi));
                }
                randPool.free_state(generator);

                const size_t j = offset + i;
                Rview(j)[0]    = r * Kokkos::cos(theta) * sigmaR[0] + R0[0];
                Rview(j)[1]    = r * Kokkos::sin(theta) * sigmaR[1] + R0[1];
                Rview(j)[2]    = R0[2];
                Pview(j)       = p;
                dtView(j)      = stepDt(i);

                const double gamma = Kokkos::sqrt(1.0 + p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                const double drift = 0.5 * c * stepDt(i) / gamma;
                Rview(j)[0] += p[0] * drift;
                Rview(j)[1] += p[1] * drift;
                Rview(j)[2] += p[2] * drift;
            });
    Kokkos::fence();
    pc_m->markMomentsDirty();
}

void OpalFlatTop::initDomainDecomp(double BoxIncr) {
    if (!fc_m) {
        return;
    }

    auto* mesh = &fc_m->getMesh();
    auto* FL   = &fc_m->getFL();
    ippl::Vector<double, 3> o;
    ippl::Vector<double, 3> e;
    const double tol = 1e-15;
    o[0]             = -sigmaR_m[0] - tol;
    e[0]             = sigmaR_m[0] + tol;
    o[1]             = -sigmaR_m[1] - tol;
    e[1]             = sigmaR_m[1] + tol;
    o[2]             = 0.0 - tol;
    e[2]             = Physics::c * emissionTime_m + tol;

    const ippl::Vector<double, 3> l = e - o;
    hr_m                            = (1.0 + BoxIncr / 100.0) * (l / nr_m);
    mesh->setMeshSpacing(hr_m);
    mesh->setOrigin(o - 0.5 * hr_m * BoxIncr / 100.0);
    pc_m->getLayout().updateLayout(*FL, *mesh);
}

void OpalFlatTop::setWithDomainDecomp(bool withDomainDecomp) {
    withDomainDecomp_m = withDomainDecomp;
}

void OpalFlatTop::setBirthTimesForTest(std::vector<double> birthTimes) {
    birthTimes_m = std::move(birthTimes);
    std::sort(birthTimes_m.begin(), birthTimes_m.end());
    totalN_m          = birthTimes_m.size();
    nextGlobalIndex_m = 0;
    inventoryBuilt_m  = true;
}
