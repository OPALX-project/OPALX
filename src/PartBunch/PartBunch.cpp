/**
 * @file PartBunch.cpp
 * @brief Template method definitions for PartBunch.
 */

#include "PartBunch/PartBunch.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include "Particle/ParticleAttrib.h"
#include "Physics/ParticleProperties.h"
#include "Structure/Beam.h"
#include "Utilities/Util.h"

#undef doDEBUG

/**
 * @copybrief PartBunch::PartBunch
 */
template <typename T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(
        std::vector<double> qi, std::vector<double> mi, const std::vector<Beam*>& beams,
        std::vector<size_t> totalParticlesPerBeam, double lbt, std::string integration_method,
        opalx::spacecharge::CartesianDomainConfig<T, Dim> domainConfig)
    : dt_m(0),
      it_m(0),
      integration_method_m(integration_method),
      lbt_m(lbt),
      boundingBoxIncreasePercent_m(0.0),
      t_m(0.0),
      globalTrackStep_m(0),
      rmsDensity_m(0.0) {
    qi_m         = qi;
    mi_m         = mi;
    bunchState_m = std::make_shared<BunchStateHandler>();

    Inform m("PartBunch::PartBunch");
    m << level4 << "PartBunch Constructor" << endl;

    const size_t num_containers = beams.size();
    if (num_containers == 0) {
        throw OpalException("PartBunch::PartBunch", "num_containers must be > 0.");
    }
    if (qi.size() != num_containers) {
        throw OpalException("PartBunch::PartBunch", "qi size must match num_containers.");
    }
    if (mi.size() != num_containers) {
        throw OpalException("PartBunch::PartBunch", "mi size must match num_containers.");
    }
    if (totalParticlesPerBeam.size() != num_containers) {
        throw OpalException(
                "PartBunch::PartBunch", "totalParticlesPerBeam size must match num_containers.");
    }
    for (size_t i = 0; i < num_containers; ++i) {
        if (beams[i] == nullptr) {
            throw OpalException("PartBunch::PartBunch", "beams must not contain null pointers.");
        }
    }

    // The initial geometry only provides a valid mesh and layout while particles are sampled or
    // restored. CartesianPIC3D replaces the physical bounds before its first runtime solve, while
    // the mesh and FieldLayout object addresses remain stable for every borrowing container.
    cartesianDomain_m            = std::make_unique<CartesianDomain_t>(domainConfig);
    boundingBoxIncreasePercent_m = domainConfig.boundingBoxIncreasePercent;
    for (unsigned dimension = 0; dimension < Dim; ++dimension) {
        nr_m[dimension]     = static_cast<int>(domainConfig.meshSize[dimension]);
        domain_m[dimension] = ippl::Index(nr_m[dimension]);
        decomp_m[dimension] = domainConfig.decomposition[dimension];
    }

    const bool useOverlap =
            domainConfig.layoutType == opalx::spacecharge::ParticleLayoutType::SpatialOverlap;
    T overlapCutoff = opalx::spacecharge::resolveOverlapCutoff<Dim>(
            domainConfig.overlapCutoff, domainConfig.overlapCutoffCells,
            cartesianDomain_m->spacing());
    // A mesh-tied cutoff follows the placeholder spacing, which can exceed half a rank region
    // when an axis has few cells per rank. The runtime update replaces it before the first solve.
    if (useOverlap && domainConfig.overlapCutoffCells > T(0)) {
        overlapCutoff = std::min(overlapCutoff, cartesianDomain_m->largestFittingOverlapCutoff());
    }
    if (useOverlap && !(overlapCutoff > T(0))) {
        throw OpalException(
                "PartBunch::PartBunch",
                "A spatial-overlap particle layout requires a positive cutoff.");
    }
    const auto layoutType           = useOverlap ? ParticleContainer_t::LayoutType::SpatialOverlap
                                                 : ParticleContainer_t::LayoutType::Spatial;
    const bool isAllPeriodic        = domainConfig.periodicParticleBoundary;
    const ippl::BC particleBoundary = isAllPeriodic ? ippl::BC::PERIODIC : ippl::BC::NO;
    m << level5 << "* Initial Cartesian domain set to isAllPeriodic = " << isAllPeriodic << endl;

    this->hr_m     = cartesianDomain_m->spacing();
    this->origin_m = cartesianDomain_m->origin();
    this->dt_m     = 0.5 / this->nr_m[2];

    std::ostringstream fieldLayout;
    fieldLayout << cartesianDomain_m->layout();
    initialFieldLayout_m = fieldLayout.str();

    pcontainer_m = std::make_shared<ParticleContainer_t>(
            cartesianDomain_m->mesh(), cartesianDomain_m->layout(), beams[0]->hasPolarization(),
            layoutType, overlapCutoff, particleBoundary);
    pcontainers_m.push_back(pcontainer_m);
    pcontainer_m->setBunchStateHandler(bunchState_m);
    /// \todo if we want, we could also have a separate BunchStateHandler for each container later?
    /// But I think it could also make sense to only have one global handler.
    for (size_t i = 1; i < num_containers; ++i) {
        auto pc = std::make_shared<ParticleContainer_t>(
                cartesianDomain_m->mesh(), cartesianDomain_m->layout(), beams[i]->hasPolarization(),
                layoutType, overlapCutoff, particleBoundary);
        pc->setBunchStateHandler(bunchState_m);
        pcontainers_m.push_back(std::move(pc));
    }
    const auto& containers = this->getParticleContainers();
    particleNames_m.resize(containers.size());
    pcRoot_m.resize(containers.size());
    for (size_t i = 0; i < containers.size(); ++i) {
        pcRoot_m[i] = i;  // independent until TrackRun links decay chains
    }
    for (size_t i = 0; i < containers.size(); ++i) {
        containers[i]->setQ(qi[i]);
        containers[i]->setM(mi[i]);
        containers[i]->setReference(&beams[i]->getReference());
        particleNames_m[i] = beams[i]->getParticleName();
        containers[i]->Sp  = static_cast<short>(
                ParticleProperties::getParticleType(beams[i]->getParticleName()));
    }

    // Pre-allocate per-rank capacity without bumping localNum_m. Subsequent emission /
    // distribution loaders call createParticles() to fill the buffer non-destructively.
    const double nRanks = static_cast<double>(ippl::Comm->size());
    for (size_t i = 0; i < num_containers; ++i) {
        const size_t maxLocalNum =
                static_cast<size_t>(totalParticlesPerBeam[i] / nRanks + 2 * nRanks + 1);
        containers[i]->allocateParticles(maxLocalNum);
        *gmsg << level3 << "* Container " << i << ": capacity for " << maxLocalNum
              << " particles allocated." << endl;
    }

    globalPartPerNode_m = std::make_unique<size_t[]>(ippl::Comm->size());

    resetPcActive();

    m << level5 << "* PartBunch constructor done." << endl;
}

/**
 * @copybrief PartBunch::resetPcActive
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::resetPcActive() {
    const auto& containers = this->getParticleContainers();
    const size_t n         = containers.size();
    pcActive_m.resize(n);
    pcAtSStop_m.resize(n);
    for (size_t i = 0; i < n; ++i) {
        pcAtSStop_m[i] = false;
        const auto& pc = containers[i];
        if (!pc || pc->getTotalNum() == 0) {
            pcActive_m[i] = false;
        } else {
            pcActive_m[i] = true;
        }
    }
}

/**
 * @copybrief PartBunch::setPcAtSStop
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setPcAtSStop(size_t i) {
    if (i >= pcActive_m.size()) {
        return;
    }
    pcActive_m[i]  = false;
    pcAtSStop_m[i] = true;
}

/**
 * @copybrief PartBunch::refreshPcActiveAfterEmit
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::refreshPcActiveAfterEmit() {
    const auto& containers = this->getParticleContainers();
    const size_t n         = containers.size();
    if (pcActive_m.size() != n) {
        return;
    }
    for (size_t i = 0; i < n; ++i) {
        if (pcAtSStop_m[i]) {
            continue;
        }
        const auto& pc = containers[i];
        if (pc && pc->getTotalNum() > 0) {
            pcActive_m[i] = true;
        }
    }
}

/**
 * @copybrief PartBunch::gatherLoadBalanceStatistics
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::gatherLoadBalanceStatistics() {
    std::fill_n(globalPartPerNode_m.get(), ippl::Comm->size(), 0);  // Fill the array with zeros
    globalPartPerNode_m[ippl::Comm->rank()] = this->getParticleContainer()->getLocalNum();
    ippl::Comm->allreduce(globalPartPerNode_m.get(), ippl::Comm->size(), std::plus<size_t>());
}

/**
 * @copybrief PartBunch::print
 */
template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os, std::span<const opalx::moments::BeamMoments> moments) {
    Inform::FmtFlags_t ff = os.flags();

    const auto& containers = this->getParticleContainers();
    if (moments.size() != containers.size()) {
        throw OpalException("PartBunch::print", "Expected one moment snapshot per container.");
    }
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        const auto& m  = moments[ci];

        // ParticleContainer tracks charge/mass storage mode for QM attributes.
        std::string qmStorageModeStr = "SINGLE";
        const auto qmMode            = pc->getQMStorageMode();
        if (qmMode == ParticleContainer_t::QMStorageMode::Attributes) {
            qmStorageModeStr = "ATTRIBUTES";
        }

        os << level1 << std::scientific << "\n"
           << "* ************** B U N C H "
              "********************************************************* \n"
           << "* CONTAINER       = " << ci << "\n"
           << "* PARTICLES       = " << pc->getTotalNum() << "\n"
           << "* CHARGE          = " << pc->getTotalCharge() << " (Cb) \n"
           << "* QM STORAGE MODE = " << qmStorageModeStr << "\n"
           << "* <EKIN>          = " << Util::getEnergyString(m.meanKineticEnergy()) << "\n"
           << "* <dEKIN>         = " << Util::getEnergyString(m.stdKineticEnergy()) << "\n"
           << "* INTEGRATOR      = " << integration_method_m << "\n"
           << "* LB Threshold    = " << lbt_m << "\n"
           << "* MIN R (origin)  = " << Util::getLengthString(m.minR(), 5) << "\n"
           << "* MAX R (max ext) = " << Util::getLengthString(m.maxR(), 5) << "\n"
           << "* RMS R           = " << Util::getLengthString(m.rmsR(), 5) << "\n"
           << "* RMS P           = " << m.rmsP() << " [beta gamma]\n"
           << "* Mean R          = " << m.meanR() << " [m]\n"
           << "* Mean P          = " << m.meanP() << " [beta gamma]\n"
           << "* MESH SPACING    = " << Util::getLengthString(hr_m, 5) << "\n"
           << "* COMPDOM INCR    = " << boundingBoxIncreasePercent_m << " (%) \n"
           << "* FIELD LAYOUT    = " << initialFieldLayout_m << "\n"
           << "* Cov Matrix : \n* ";
        for (unsigned int i = 0; i < 2 * Dim; i++) {
            for (unsigned int j = 0; j < 2 * Dim; j++) {
                os << level1 << m.covariance(i, j) << " ";
            }
            os << level1 << "\n* ";
        }
        os << level1
           << "* "
              "********************************************************************************"
              "** \n"
           << endl;
    }

    os.flags(ff);
    return os;
}

/** Explicit instantiation for 3D double (OPAL-T). */
template class PartBunch<double, 3>;
