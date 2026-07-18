/**
 * @file PartBunch.cpp
 * @brief Template method definitions for PartBunch.
 */

#include "PartBunch/PartBunch.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include "Algorithms/Matrix.h"
#include "PartBunch/BCHandler.hpp"
#include "Particle/ParticleAttrib.h"
#include "Physics/ParticleProperties.h"
#include "SpaceCharge/Pic3D/PicWorkspace.h"
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
        FieldSolverCmd* OPALFieldSolver)
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
    if (OPALFieldSolver == nullptr) {
        throw OpalException("PartBunch::PartBunch", "OPALFieldSolver must not be null.");
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

    //  get the needed information from OPAL FieldSolver command

    boundingBoxIncreasePercent_m = OPALFieldSolver->getBoxIncr();
    nr_m                         = Vector_t<int, Dim>(
            OPALFieldSolver->getNX(), OPALFieldSolver->getNY(), OPALFieldSolver->getNZ());

    const Vector_t<bool, 3> domainDecomposition = OPALFieldSolver->getDomainDecomposition();

    for (unsigned i = 0; i < Dim; i++) {
        this->domain_m[i] = ippl::Index(nr_m[i]);
        this->decomp_m[i] = domainDecomposition[i];
    }

    const BCHandler<Dim> bcHandler = OPALFieldSolver->constructBCHandler();

    // TODO: support mixed periodic/open per axis; currently all periodic or all open.
    const bool isAllPeriodic = bcHandler.isAll(BCHandler<Dim>::PERIODIC);
    m << level5 << "* Initial PIC workspace set to isAllPeriodic = " << isAllPeriodic << endl;

    // Set up the initial particle layout. Pic3DSolver replaces this geometry with the
    // physical computational domain before the first runtime solve.

    Vector_t<double, Dim> length(6.0);
    this->hr_m     = length / this->nr_m;
    this->origin_m = -3.0;
    this->dt_m     = 0.5 / this->nr_m[2];

    rmin_m = origin_m;
    rmax_m = origin_m + length;

    picWorkspace_m = std::make_shared<PicWorkspace_t>(
            hr_m, rmin_m, rmax_m, decomp_m, domain_m, origin_m, isAllPeriodic);
    std::ostringstream fieldLayout;
    fieldLayout << picWorkspace_m->layout();
    initialFieldLayout_m = fieldLayout.str();

    pcontainer_m = std::make_shared<ParticleContainer_t>(
            picWorkspace_m->mesh(), picWorkspace_m->layout(), beams[0]->hasPolarization());
    pcontainers_m.push_back(pcontainer_m);
    pcontainer_m->setBunchStateHandler(bunchState_m);
    /// \todo if we want, we could also have a separate BunchStateHandler for each container later?
    /// But I think it could also make sense to only have one global handler.
    for (size_t i = 1; i < num_containers; ++i) {
        auto pc = std::make_shared<ParticleContainer_t>(
                picWorkspace_m->mesh(), picWorkspace_m->layout(), beams[i]->hasPolarization());
        pc->setBunchStateHandler(bunchState_m);
        pcontainers_m.push_back(std::move(pc));
    }
    const auto& containers = this->getParticleContainers();
    particleNames_m.resize(containers.size());
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
 * @copybrief PartBunch::calcBeamParameters
 */
template <typename T, unsigned Dim>
void PartBunch<T, Dim>::calcBeamParameters() {
    Inform m("PartBunch::calcBeamParameters");
    std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;

    using view_type = ippl::ParticleAttrib<Vector_t<double, 3>>::view_type;
    view_type Rview = pc->R.getView();
    view_type Pview = pc->P.getView();
    this->getParticleContainer()->updateMoments();
    m << level5 << "Moments updated." << endl;

    ////////////////////////////////////
    //// Calculate Moments of R and P //
    ////////////////////////////////////

    using MomentsVec = ippl::Vector<double, 2 * Dim>;
    using MomentsMat = ippl::Vector<MomentsVec, 2 * Dim>;

    MomentsVec loc_centroid(0.0);
    MomentsMat loc_moment(MomentsVec(0.0));

    MomentsVec centroid(0.0);
    MomentsMat moment(MomentsVec(0.0));

    for (unsigned i = 0; i < 2 * Dim; ++i) {
        const size_t nLocal = pc->getLocalNum();
        Kokkos::parallel_reduce(
                "calc moments of particle distr.", nLocal,
                KOKKOS_LAMBDA(
                        const size_t k, double& cent, double& mom0, double& mom1, double& mom2,
                        double& mom3, double& mom4, double& mom5) {
                    double part[2 * Dim];
                    part[0] = Rview(k)[0];
                    part[1] = Pview(k)[0];
                    part[2] = Rview(k)[1];
                    part[3] = Pview(k)[1];
                    part[4] = Rview(k)[2];
                    part[5] = Pview(k)[2];

                    cent += part[i];
                    mom0 += part[i] * part[0];
                    mom1 += part[i] * part[1];
                    mom2 += part[i] * part[2];
                    mom3 += part[i] * part[3];
                    mom4 += part[i] * part[4];
                    mom5 += part[i] * part[5];
                },
                Kokkos::Sum<T>(loc_centroid[i]), Kokkos::Sum<T>(loc_moment[i][0]),
                Kokkos::Sum<T>(loc_moment[i][1]), Kokkos::Sum<T>(loc_moment[i][2]),
                Kokkos::Sum<T>(loc_moment[i][3]), Kokkos::Sum<T>(loc_moment[i][4]),
                Kokkos::Sum<T>(loc_moment[i][5]));
        Kokkos::fence();
    }
    m << level5 << "Local moments calculated." << endl;

    moment   = loc_moment;
    centroid = loc_centroid;
    ippl::Comm->allreduce(moment, 1, std::plus<MomentsMat>());
    ippl::Comm->allreduce(centroid, 1, std::plus<MomentsVec>());
    ippl::Comm->barrier();
    m << level5 << "Global moments calculated." << endl;

    ippl::Vector<double, Dim> rmax_loc(0.0);
    ippl::Vector<double, Dim> rmin_loc(0.0);
    ippl::Vector<double, Dim> rmax(0.0);
    ippl::Vector<double, Dim> rmin(0.0);

    // TODO: fuse min/max reductions with ippl::Vector reductions.
    for (unsigned d = 0; d < Dim; ++d) {
        Kokkos::parallel_reduce(
                "rel max", pc->getLocalNum(),
                KOKKOS_LAMBDA(const int i, double& mm) {
                    double tmp_vel = Rview(i)[d];
                    mm             = tmp_vel > mm ? tmp_vel : mm;
                },
                Kokkos::Max<T>(rmax_loc[d]));

        Kokkos::parallel_reduce(
                "rel min", pc->getLocalNum(),
                KOKKOS_LAMBDA(const int i, double& mm) {
                    double tmp_vel = Rview(i)[d];
                    mm             = tmp_vel < mm ? tmp_vel : mm;
                },
                Kokkos::Min<T>(rmin_loc[d]));
    }
    m << level5 << "Local min/max calculated." << endl;
    Kokkos::fence();
    rmax = rmax_loc;
    rmin = rmin_loc;
    ippl::Comm->allreduce(rmax, 1, std::greater<ippl::Vector<double, Dim>>());
    ippl::Comm->allreduce(rmin, 1, std::less<ippl::Vector<double, Dim>>());
    ippl::Comm->barrier();
    m << level5 << "Global min/max calculated." << endl;

    rmax_m = rmax;
    rmin_m = rmin;
}

/**
 * @copybrief PartBunch::print
 */
template <typename T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os) {
    // if (pc->getLocalNum() != 0) {  // to suppress Nans
    Inform::FmtFlags_t ff = os.flags();

    const auto& containers = this->getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc) {
            os << level1 << "Skipping null particle container: " << ci << endl;
            continue;
        }

        const double ek  = pc->getMeanKineticEnergy();
        const double dek = pc->getStdKineticEnergy();

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
           << "* <EKIN>          = " << Util::getEnergyString(ek) << "\n"
           << "* <dEKIN>         = " << Util::getEnergyString(dek) << "\n"
           << "* INTEGRATOR      = " << integration_method_m << "\n"
           << "* LB Threshold    = " << lbt_m << "\n"
           << "* MIN R (origin)  = " << Util::getLengthString(pc->getMinR(), 5) << "\n"
           << "* MAX R (max ext) = " << Util::getLengthString(pc->getMaxR(), 5) << "\n"
           << "* RMS R           = " << Util::getLengthString(pc->getRmsR(), 5) << "\n"
           << "* RMS P           = " << pc->getRmsP() << " [beta gamma]\n"
           << "* Mean R          = " << pc->getMeanR() << " [m]\n"
           << "* Mean P          = " << pc->getMeanP() << " [beta gamma]\n"
           << "* MESH SPACING    = " << Util::getLengthString(hr_m, 5) << "\n"
           << "* COMPDOM INCR    = " << boundingBoxIncreasePercent_m << " (%) \n"
           << "* FIELD LAYOUT    = " << initialFieldLayout_m << "\n"
           << "* Centroid : \n* ";
        for (unsigned int i = 0; i < 2 * Dim; i++) {
            os << level1 << pc->getCentroid()[i] << " ";
        }
        os << level1 << endl << "* Cov Matrix : \n* ";
        for (unsigned int i = 0; i < 2 * Dim; i++) {
            for (unsigned int j = 0; j < 2 * Dim; j++) {
                os << level1 << pc->getCovMatrix()(i, j) << " ";
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

template <typename T, unsigned Dim>
void PartBunch<T, Dim>::updateAllParticleMoments() {
    const auto& containers = this->getParticleContainers();
    for (size_t i = 0; i < containers.size(); ++i) {
        if (!containers[i]) {
            continue;
        }
        this->getParticleContainer(i)->updateMoments();
    }
}

template <typename T, unsigned Dim>
std::shared_ptr<typename PartBunch<T, Dim>::PicWorkspace_t> PartBunch<T, Dim>::takePicWorkspace() {
    if (!picWorkspace_m) {
        throw OpalException(
                "PartBunch::takePicWorkspace",
                "The setup-time PIC workspace was already transferred.");
    }
    return std::move(picWorkspace_m);
}

/** Explicit instantiation for 3D double (OPAL-T). */
template class PartBunch<double, 3>;
