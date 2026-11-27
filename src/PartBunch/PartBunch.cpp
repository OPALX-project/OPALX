/**
 * @file PartBunch.cpp
 * @brief Template method definitions for PartBunch.
 */

#include "PartBunch/PartBunch.h"
#include <algorithm>
#include "Physics/ParticleProperties.h"
#include "Processes/GlobalProcesses/GlobalProcess.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

/**
 * @copybrief PartBunch::PartBunch
 */
template <typename T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(
        const std::vector<BunchBeam>& beams, double lbt,
        opalx::spacecharge::CartesianDomainConfig<T, Dim> domainConfig, bool qmAttributes,
        bool aggressiveStateSync)
    : dt_m(0.0), lbt_m(lbt), boundingBoxIncreasePercent_m(0.0), t_m(0.0), globalTrackStep_m(0) {
    bunchState_m = std::make_shared<BunchStateHandler>(aggressiveStateSync);

    Inform m("PartBunch::PartBunch");

    const size_t num_containers = beams.size();
    if (num_containers == 0) {
        throw OpalException("PartBunch::PartBunch", "num_containers must be > 0.");
    }

    // The initial geometry only provides a valid mesh and layout while particles are sampled or
    // restored. CartesianPIC3D replaces the physical bounds before its first runtime solve, while
    // the mesh and FieldLayout object addresses remain stable for every borrowing container.
    cartesianDomain_m            = std::make_unique<CartesianDomain_t>(domainConfig);
    boundingBoxIncreasePercent_m = domainConfig.boundingBoxIncreasePercent;

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

    // All containers share one BunchStateHandler.
    for (size_t i = 0; i < num_containers; ++i) {
        auto pc = std::make_shared<ParticleContainer_t>(
                cartesianDomain_m->mesh(), cartesianDomain_m->layout(), beams[i].hasSpin,
                layoutType, overlapCutoff, particleBoundary, qmAttributes);
        pc->setBunchStateHandler(bunchState_m);
        pcontainers_m.push_back(std::move(pc));
    }
    particleNames_m.resize(num_containers);
    pcRoot_m.resize(num_containers);
    pcActive_m.resize(num_containers);
    pcAtSStop_m.resize(num_containers);
    for (size_t i = 0; i < num_containers; ++i) {
        pcRoot_m[i] = i;  // buildTracking() links decay chains with setPcRoot()
    }
    for (size_t i = 0; i < num_containers; ++i) {
        pcontainers_m[i]->setQ(beams[i].chargePerParticle);
        pcontainers_m[i]->setM(beams[i].massPerParticle);
        pcontainers_m[i]->setReference(beams[i].reference);
        particleNames_m[i] = beams[i].particleName;
        pcontainers_m[i]->Sp =
                static_cast<short>(ParticleProperties::getParticleType(beams[i].particleName));
    }

    // Reserve room for each beam's particles on this rank; the containers stay empty. Samplers
    // and emission add particles later with createParticles(), which keeps existing particles.
    const double nRanks = static_cast<double>(ippl::Comm->size());
    for (size_t i = 0; i < num_containers; ++i) {
        const size_t maxLocalNum =
                static_cast<size_t>(beams[i].numParticles / nRanks + 2 * nRanks + 1);
        pcontainers_m[i]->allocateParticles(maxLocalNum);
        *gmsg << level3 << "* Container " << i << ": capacity for " << maxLocalNum
              << " particles allocated." << endl;
    }

    resetPcActive();
}

/**
 * @copybrief PartBunch::resetPcActive
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::resetPcActive() {
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        pcAtSStop_m[i] = false;
        pcActive_m[i]  = pcontainers_m[i]->getTotalNum() > 0;
    }
}

/**
 * @copybrief PartBunch::setPcAtSStop
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setPcAtSStop(size_t i) {
    pcActive_m.at(i)  = false;
    pcAtSStop_m.at(i) = true;
}

/**
 * @copybrief PartBunch::refreshPcActiveAfterEmit
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::refreshPcActiveAfterEmit() {
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        if (!pcAtSStop_m[i] && pcontainers_m[i]->getTotalNum() > 0) {
            pcActive_m[i] = true;
        }
    }
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::setTimeStep(double dt) {
    this->setdT(dt);
    for (const auto& pc : pcontainers_m) {
        pc->dt = dt;
    }
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::pushParticles() {
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        if (pcActive_m[i]) {
            pcontainers_m[i]->pushParticles();
        }
    }
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::kickParticles() {
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        if (pcActive_m[i]) {
            pcontainers_m[i]->kickParticles();
        }
    }
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::evolveSpinHalfStep() {
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        if (pcActive_m[i] && pcontainers_m[i]->hasSpin()) {
            pcontainers_m[i]->evolveSpinHalfStep();
        }
    }
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::resetFields() {
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        if (pcActive_m[i]) {
            pcontainers_m[i]->zeroFields();
        }
    }
}

template <typename T, unsigned Dim>
size_t PartBunch<T, Dim>::applyGlobalProcesses() {
    // Decay writes daughters in the shared reference frame of parent and daughter, so every
    // container must be in its reference frame.
    requireReferenceFrame("PartBunch::applyGlobalProcesses");
    size_t nMarked = 0;
    for (size_t ci = 0; ci < pcontainers_m.size(); ++ci) {
        // A beam that has passed its ZSTOP is frozen and does not decay any more. An empty
        // container is not skipped: a daughter container that an earlier container of this loop
        // fills applies its own processes in the same call.
        if (pcAtSStop_m[ci]) {
            continue;
        }
        const auto& pc = pcontainers_m[ci];
        for (const auto& process : pc->getGlobalProcesses()) {
            if (process) {
                nMarked += process->apply(*pc, dt_m, globalTrackStep_m, ci);
            }
        }
    }
    // A decay may have put the first particles into an empty daughter container. Switch it on
    // now, so its particles take part in the next step's first push. nMarked is the global
    // count, so every rank does the same.
    if (nMarked > 0) {
        this->refreshPcActiveAfterEmit();
    }
    return nMarked;
}

template <typename T, unsigned Dim>
std::vector<opalx::moments::BeamMoments> PartBunch<T, Dim>::markOutliers(
        double sigmas,
        const std::function<opalx::moments::Request(const ParticleContainer_t&)>& request,
        size_t& nMarked) {
    using opalx::moments::Moment;
    std::vector<opalx::moments::BeamMoments> moments(pcontainers_m.size());
    if (sigmas <= 0.0) {
        return moments;
    }
    Inform m("PartBunch::markOutliers");
    for (size_t i = 0; i < pcontainers_m.size(); ++i) {
        const auto& pc = pcontainers_m[i];
        if (!pcActive_m[i] || pc->getTotalNum() == 0) {
            continue;
        }
        moments[i] =
                opalx::moments::compute(*pc, Moment::MeanPosition | Moment::Spread | request(*pc));
        const Vector_t<double, 3> halfWidth = sigmas * moments[i].rmsR();
        const size_t nContainerMarked       = pc->markParticlesOutside(
                moments[i].meanR() - halfWidth, moments[i].meanR() + halfWidth);
        nMarked += nContainerMarked;
        if (nContainerMarked > 0) {
            m << level2 << "Marked " << nContainerMarked << " particles outside " << sigmas
              << "-sigma boundary for deletion (container " << i << ")." << endl;
        }
    }
    return moments;
}

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::deleteInvalidParticles(
        bool activeOnly, Inform& m, const std::string& reason) {
    for (size_t ci = 0; ci < pcontainers_m.size(); ++ci) {
        if (activeOnly && !pcActive_m[ci]) {
            continue;
        }
        const auto& pc                 = pcontainers_m[ci];
        const size_t nContainerDeleted = pc->deleteInvalidParticles();
        if (nContainerDeleted > 0) {
            m << level2 << "Deleted " << nContainerDeleted << " particles marked by " << reason
              << ", " << pc->getTotalNum() << " remaining (container " << ci << ")." << endl;
        }
    }
}

/**
 * @copybrief PartBunch::print
 */
template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(
        Inform& os, std::span<const opalx::moments::BeamMoments> moments) const {
    Inform::FmtFlags_t ff = os.flags();

    if (moments.size() != pcontainers_m.size()) {
        throw OpalException("PartBunch::print", "Expected one moment snapshot per container.");
    }
    for (size_t ci = 0; ci < pcontainers_m.size(); ++ci) {
        const auto& pc = pcontainers_m[ci];
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
           << "* LB Threshold    = " << lbt_m << "\n"
           << "* MIN R (origin)  = " << Util::getLengthString(m.minR(), 5) << "\n"
           << "* MAX R (max ext) = " << Util::getLengthString(m.maxR(), 5) << "\n"
           << "* RMS R           = " << Util::getLengthString(m.rmsR(), 5) << "\n"
           << "* RMS P           = " << m.rmsP() << " [beta gamma]\n"
           << "* Mean R          = " << m.meanR() << " [m]\n"
           << "* Mean P          = " << m.meanP() << " [beta gamma]\n"
           << "* MESH SPACING    = " << Util::getLengthString(cartesianDomain_m->spacing(), 5)
           << "\n"
           << "* COMPDOM INCR    = " << boundingBoxIncreasePercent_m << " (%) \n"
           << "* FIELD LAYOUT    = " << cartesianDomain_m->layout() << "\n"
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

/** Explicit instantiation for 3D double. */
template class PartBunch<double, 3>;
