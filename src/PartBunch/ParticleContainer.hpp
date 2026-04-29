
#ifndef OPAL_PARTICLE_CONTAINER_H
#define OPAL_PARTICLE_CONTAINER_H

// #include <functional>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "Manager/BaseManager.h"

#include "PartBunch/FieldContainer.hpp"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/DistributionMoments.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"

#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

// #include <Kokkos_Core.hpp>

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

using size_type = ippl::detail::size_type;

#include "Processes/GlobalProcesses/GlobalProcess.h"

/**
 * @class ParticleContainer
 * @brief Container for all per-particle (and per-simulation) fields tracked during OPALX tracking.
 *
 * The values tracked in Kokkos::Views during the simulation are:
 * R  - Position (from base class)
 * P  - Momentum [beta*gamma]
 * dt - Time step
 * Phi- Scalar potential
 * Bin- Energy bin
 * E  - Electric field
 * B  - Magnetic field
 *
 * Charge (Q) and mass (M):
 * - Default: QM_MODE="SINGLE" -> QMStorageMode=SingleValue
 *   * Q and M are stored as a single shared value per container (memory-efficient).
 *   * Access via `getQView()` / `getMView()` returns the shared views.
 * - Alternative: QM_MODE="ATTRIBUTES" -> QMStorageMode=Attributes
 *   * Q and M are stored as per-particle attributes.
 *   * Access via `getQView()` / `getMView()` returns the per-particle attribute views.
 *
 * Access to Q/M should be done with `getQView()` / `getMView()`.
 * They automatically select the correct underlying storage mode.
 */
template <typename T, unsigned Dim = 3>
class ParticleContainer : public ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>> {
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

private:
    /**
     * Forbid access of the following functions outside of the ParticleContainer wrappers! This is
     * for safety, since it might lead to undefined behaviour. The idea is to handle particle count
     * changes completely through the ParticleContainer!
     */
    using Base::alloc;
    using Base::create;
    using Base::destroy;

public:
    enum class QMStorageMode { SingleValue, Attributes };

    /// Defines which type to use as a particle bin.
    using bin_index_type = short int;  // Needed in AdaptBins class

    /// View types of Q and M values
    using qm_view_type = typename ippl::ParticleAttrib<double>::view_type;

public:
    /// Charge view in [Cb].
    /// In `SingleValue` mode this is a rank-1 view of length 1.
    /// In `Attributes` mode this is the per-particle attribute view.
    qm_view_type getQView() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return QAttr.getView();
        }
        return QView_m;
    }

    /// Mass view in [GeV].
    /// In `SingleValue` mode this is a rank-1 view of length 1.
    /// In `Attributes` mode this is the per-particle attribute view.
    qm_view_type getMView() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return MAttr.getView();
        }
        return MView_m;
    }

    /// timestep in [s]
    ippl::ParticleAttrib<double> dt;

    /// the scalar potential in [Cb/s]
    ippl::ParticleAttrib<double> Phi;

    /// the energy bin the particle is in
    ippl::ParticleAttrib<bin_index_type> Bin;

    /// the particle specis
    short Sp = 0;

    /// particle momenta [\beta\gamma]
    typename Base::particle_position_type P;

    /// electric field at particle position
    typename Base::particle_position_type E;

    /// electric field for gun simulation with bins
    // typename Base::particle_position_type Etmp; // TODO: might not need this...

    /// magnetic field at particle position
    typename Base::particle_position_type B;

    ParticleContainer(Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL)
        : pl_m(FL, mesh),
          qmStorageMode_m(
                  Options::useQMAttributes ? QMStorageMode::Attributes
                                           : QMStorageMode::SingleValue),
          distMoments_m(),
          QView_m("ParticleContainer::QView_m", 1),
          MView_m("ParticleContainer::MView_m", 1) {
        this->initialize(pl_m);
        registerAttributes();
        setupBCs();
        Kokkos::deep_copy(QView_m, 0.0);
        Kokkos::deep_copy(MView_m, 0.0);
    }

    ~ParticleContainer() {}

    void registerAttributes() {
        // register the particle attributes
        this->addAttribute(dt);
        this->addAttribute(Phi);
        this->addAttribute(Bin);
        this->addAttribute(P);
        this->addAttribute(E);
        this->addAttribute(B);

        if (qmStorageMode_m == QMStorageMode::Attributes) {
            this->addAttribute(QAttr);
            this->addAttribute(MAttr);
        }
    }

    void setupBCs() { setBCAllPeriodic(); }

    /// Apply coordinate transform to local particles: translate R, rotate P, E, B.
    void transformBunch(const CoordinateSystemTrafo& trafo) {
        const size_t nLoc = this->getLocalNum();
        trafo.transformBunchTo(this->R.getView(), nLoc);
        trafo.rotateBunchTo(this->P.getView(), nLoc);
        trafo.rotateBunchTo(this->E.getView(), nLoc);
        trafo.rotateBunchTo(this->B.getView(), nLoc);
    }
    PLayout_t<T, Dim>& getPL() { return pl_m; }

    void setBunchStateHandler(std::shared_ptr<BunchStateHandler> handler) {
        // We only keep the slot: per-container flags own their own sync, so
        // no back-reference to the handler is needed. The handler itself
        // manages bunch-wide state, which ParticleContainer never touches.
        containerState_m = handler->registerContainer();
        distMoments_m.setContainerState(containerState_m);
    }

    // -- per-container state pass-throughs --------------------------------
    // Thin wrappers around the slot's own API so callers don't need to reach
    // into `containerState_m` directly.

    bool isUnitlessPositions() const { return containerState_m->unitlessPositions; }

    bool isMomentsDirty() const { return containerState_m->momentsDirty; }
    void markMomentsDirty() { containerState_m->markMomentsDirty(); }
    void markMomentsClean() { containerState_m->markMomentsClean(); }

    void updateMoments() {
        /*
        Quick check that this container was registered with a BunchStateHandler. If not, throw an
        error. This seems to be an easy to make error when interacting (e.g. through unit tests)
        with the ParticleContainer.
        */
        if (!containerState_m) {
            throw OpalException(
                    "ParticleContainer::updateMoments",
                    "BunchStateHandler not set in ParticleContainer (containerState is null).");
        }

        size_t Np = this->getTotalNum();
        Np = (Np == 0) ? 1 : Np;  // only used for normalization in the moments class --> avoid
                                  // division by zero

        size_t Nlocal = this->getLocalNum();

        distMoments_m.computeMoments(this->R.getView(), this->P.getView(), getMView(), Np, Nlocal);
    }

    void setEnergyReferenceMass(double referenceMassGeV, bool rescaleToReference = true) {
        distMoments_m.setEnergyReferenceMass(referenceMassGeV, rescaleToReference);
    }

    Vector_t<double, 3> getMeanP() const { return distMoments_m.getMeanMomentum(); }

    Vector_t<double, 3> getRmsP() const { return distMoments_m.getStandardDeviationMomentum(); }

    Vector_t<double, 3> getMeanR() const { return distMoments_m.getMeanPosition(); }

    Vector_t<double, 3> getRmsR() const { return distMoments_m.getStandardDeviationPosition(); }

    Vector_t<double, 3> getRmsRP() const { return distMoments_m.getStandardDeviationRP(); }

    void computeMinMaxR() {
        size_t Nlocal = this->getLocalNum();
        distMoments_m.computeMinMaxPosition(this->R.getView(), Nlocal);
    }

    Vector_t<double, 3> getMinR() const { return distMoments_m.getMinPosition(); }

    Vector_t<double, 3> getMaxR() const { return distMoments_m.getMaxPosition(); }

    matrix6x6_t getCovMatrix() const { return distMoments_m.getMoments6x6(); }

    double getMeanKineticEnergy() const { return distMoments_m.getMeanKineticEnergy(); }

    double getStdKineticEnergy() const { return distMoments_m.getStdKineticEnergy(); }

    Vector_t<double, 6> getMeans() const { return distMoments_m.getMeans(); }

    Vector_t<double, 6> getCentroid() const { return distMoments_m.getCentroid(); }

    Vector_t<double, 3> getNormEmit() const { return distMoments_m.getNormalizedEmittance(); }

    Vector_t<double, 3> getGeometricEmit() const { return distMoments_m.getGeometricEmittance(); }

    double getDx() const { return distMoments_m.getDx(); }

    double getDDx() const { return distMoments_m.getDDx(); }

    double getDy() const { return distMoments_m.getDy(); }

    double getDDy() const { return distMoments_m.getDDy(); }

    double getDebyeLength() const { return distMoments_m.getDebyeLength(); }

    double getMeanGammaZ() const { return distMoments_m.getMeanGammaZ(); }

    double getTemperature() const { return distMoments_m.getTemperature(); }

    double getPlasmaParameter() const { return distMoments_m.getPlasmaParameter(); }

    double computeDebyeLength(double density) {
        size_t Np = this->getTotalNum();
        Np = (Np == 0) ? 1 : Np;  // only used for normalization in the moments class --> avoid
                                  // division by zero

        size_t Nlocal = this->getLocalNum();
        distMoments_m.computeDebyeLength(this->P.getView(), Np, Nlocal, density);
        return distMoments_m.getDebyeLength();
    }

    /**
     * @brief Set particle charge for the active Q storage mode.
     * @param q Charge value in [Cb].
     *
     * In `QMStorageMode::Attributes`, this assigns `q` to every local particle.
     * In `QMStorageMode::SingleValue`, this updates the shared scalar charge view.
     */
    void setQ(double q) {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view         = QAttr.getView();
            const size_type n = this->getLocalNum();
            if (n == 0) {
                return;
            }
            Kokkos::parallel_for(
                    "ParticleContainer::setQ", n,
                    KOKKOS_LAMBDA(const size_type i) { view(i) = q; });
            Kokkos::fence();
        } else {
            Kokkos::deep_copy(QView_m, q);
        }
    }

    /// @brief Get charge per particle [Cb].
    double getChargePerParticle() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = QAttr.getView();
            if (view.extent(0) == 0) {
                return 0.0;
            }
            double q = 0.0;
            Kokkos::deep_copy(q, Kokkos::subview(view, 0));
            return q;
        }
        double q = 0.0;
        Kokkos::deep_copy(q, Kokkos::subview(QView_m, 0));
        return q;
    }

    /// @brief Get total charge [Cb] in this container.
    double getTotalCharge() const { return getChargePerParticle() * this->getTotalNum(); }

    /**
     * @brief Set particle mass for the active M storage mode.
     * @param m Mass value in [GeV].
     *
     * In `QMStorageMode::Attributes`, this assigns `m` to every local particle.
     * In `QMStorageMode::SingleValue`, this updates the shared scalar mass view.
     */
    void setM(double m) {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view         = MAttr.getView();
            const size_type n = view.extent(0);

            Kokkos::parallel_for(
                    "ParticleContainer::setM", n,
                    KOKKOS_LAMBDA(const size_type i) { view(i) = m; });
            Kokkos::fence();
        } else {
            Kokkos::deep_copy(MView_m, m);
        }
    }

    /// @brief Get mass per particle [GeV].
    double getMassPerParticle() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = MAttr.getView();
            if (view.extent(0) == 0) {
                return 0.0;
            }
            double m = 0.0;
            Kokkos::deep_copy(m, Kokkos::subview(view, 0));
            return m;
        }
        double m = 0.0;
        Kokkos::deep_copy(m, Kokkos::subview(MView_m, 0));
        return m;
    }

    /// @brief Get total mass [GeV] in this container.
    double getTotalMass() const { return getMassPerParticle() * this->getTotalNum(); }

    /// @brief Get the reference particle position (const).
    const Vector_t<double, Dim>& getRefPartR() const { return refPartR_m; }

    /// @brief Get the reference particle position.
    Vector_t<double, Dim>& getRefPartR() { return refPartR_m; }

    /// @brief Set the reference particle position.
    void setRefPartR(const Vector_t<double, Dim>& refPartR) { refPartR_m = refPartR; }

    /// @brief Get the reference particle momentum (const).
    const Vector_t<double, Dim>& getRefPartP() const { return refPartP_m; }

    /// @brief Get the reference particle momentum.
    Vector_t<double, Dim>& getRefPartP() { return refPartP_m; }

    /// @brief Set the reference particle momentum.
    void setRefPartP(const Vector_t<double, Dim>& refPartP) { refPartP_m = refPartP; }

    /// @brief Set reference particle data.
    void setReference(const PartData* ref) {
        reference_m = ref;
        if (reference_m) {
            // PartData mass is stored in eV; DistributionMoments expects GeV.
            setEnergyReferenceMass(reference_m->getM() * Units::eV2GeV, true);
        }
    }

    /// @brief Get reference particle data.
    const PartData* getReference() const { return reference_m; }

    /// @brief Set longitudinal position along design trajectory.
    void set_sPos(double sPos) { sPos_m = sPos; }

    /// @brief Get longitudinal position along design trajectory.
    double get_sPos() const { return sPos_m; }

    /// @brief Set global-to-local rotation quaternion.
    void setGlobalToLocalQuaternion(const Quaternion_t& globalToLocalQuaternion) {
        globalToLocalQuaternion_m = globalToLocalQuaternion;
    }

    /// @brief Get global-to-local rotation quaternion.
    Quaternion_t getGlobalToLocalQuaternion() const { return globalToLocalQuaternion_m; }

    /// @brief Get local-to-lab coordinate transformation (const).
    const CoordinateSystemTrafo& getToLabTrafo() const { return toLabTrafo_m; }

    /// @brief Get local-to-lab coordinate transformation.
    CoordinateSystemTrafo& getToLabTrafo() { return toLabTrafo_m; }

    /// @brief Set local-to-lab coordinate transformation.
    void setToLabTrafo(const CoordinateSystemTrafo& toLabTrafo) { toLabTrafo_m = toLabTrafo; }

    /// @brief Advance reference/lab transform state and map bunch accordingly.
    void updateRefToLabCSTrafo(double bunchDT) {
        Vector_t<double, 3> R = toLabTrafo_m.transformFrom(refPartR_m);
        Vector_t<double, 3> P = toLabTrafo_m.rotateFrom(refPartP_m);

        const double ds =
                std::copysign(1.0, bunchDT) * std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
        sPos_m += ds;

        CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t<double, 3>(0, 0, 1)));
        transformBunch(update);
        toLabTrafo_m = toLabTrafo_m * update.inverted();
    }

    /// @brief Apply a fractional Boris step and update reference/lab transform state.
    template <typename Pusher>
    void applyFractionalStep(const Pusher& pusher, double tau, double pathLengthTarget) {
        refPartR_m /= (Physics::c * 2 * tau);
        pusher.push(refPartR_m, refPartP_m, tau);
        refPartR_m *= (Physics::c * 2 * tau);

        sPos_m = pathLengthTarget;
        toLabTrafo_m.transformFrom(refPartR_m);
        Vector_t<double, 3> R = refPartR_m;
        toLabTrafo_m.rotateFrom(refPartP_m);
        Vector_t<double, 3> P = refPartP_m;
        CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t<double, 3>(0, 0, 1)));
        toLabTrafo_m = toLabTrafo_m * update.inverted();
    }

    /**
     * @brief Scale particle time-step weights by charge before scatter.
     *
     * Multiplies each local `dt[i]` by the corresponding charge used for deposition.
     * In `QMStorageMode::Attributes`, the per-particle `QAttr(i)` is used.
     * In `QMStorageMode::SingleValue`, the shared scalar `QView_m(0)` is used.
     *
     */
    void scaleDtByCharge() {
        auto dtView       = dt.getView();
        const size_type n = this->getLocalNum();
        if (n == 0) {
            return;
        }

        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto qView = QAttr.getView();
            Kokkos::parallel_for(
                    "ParticleContainer::scaleDtByCharge(attrs)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) *= qView(i); });
        } else {
            auto QView = QView_m;
            Kokkos::parallel_for(
                    "ParticleContainer::scaleDtByCharge(single)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) *= QView(0); });
        }
        Kokkos::fence();
    }

    /**
     * @brief Restore original `dt` values after `scaleDtByCharge()`.
     *
     * Divides each local `dt[i]` by the same charge factor applied in
     * `scaleDtByCharge()`.
     *
     */
    void unscaleDtByCharge() {
        auto dtView       = dt.getView();
        const size_type n = this->getLocalNum();
        if (n == 0) {
            return;
        }

        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto qView = QAttr.getView();
            Kokkos::parallel_for(
                    "ParticleContainer::unscaleDtByCharge(attrs)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) /= qView(i); });
        } else {
            auto QView = QView_m;
            Kokkos::parallel_for(
                    "ParticleContainer::unscaleDtByCharge(single)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) /= QView(0); });
        }
        Kokkos::fence();
    }

    /**
     * @brief Transform positions to unitless coordinates using each particle's dt[i].
     *
     * Applies \f$ R'_i = R_i / (c \, dt_i) \f$. Requires valid non-zero dt values per particle.
     *
     * @throws OpalException if this container is already in unitless positions.
     */
    void switchToUnitlessPositions() {
        if (containerState_m->unitlessPositions) {
            throw OpalException(
                    "ParticleContainer::switchToUnitlessPositions",
                    "ParticleContainer is already in unitless positions!");
        }
        auto Rview             = this->R.getView();
        auto dtview            = this->dt.getView();
        const size_type nLocal = this->getLocalNum();
        Kokkos::parallel_for(
                "ParticleContainer::switchToUnitlessPositions", nLocal,
                KOKKOS_LAMBDA(const size_type i) { Rview(i) *= 1.0 / (Physics::c * dtview(i)); });
        Kokkos::fence();
        containerState_m->setUnitlessPositions(true);
    }

    /**
     * @brief Restore physical positions from unitless form using each particle's dt[i].
     *
     * Applies \f$ R_i = R'_i \, c \, dt_i \f$.
     *
     * @throws OpalException if this container is not currently in unitless positions.
     */
    void switchOffUnitlessPositions() {
        if (!containerState_m->unitlessPositions) {
            throw OpalException(
                    "ParticleContainer::switchOffUnitlessPositions",
                    "ParticleContainer is already in physical positions!");
        }
        auto Rview             = this->R.getView();
        auto dtview            = this->dt.getView();
        const size_type nLocal = this->getLocalNum();
        Kokkos::parallel_for(
                "ParticleContainer::switchOffUnitlessPositions", nLocal,
                KOKKOS_LAMBDA(const size_type i) { Rview(i) *= Physics::c * dtview(i); });
        Kokkos::fence();
        containerState_m->setUnitlessPositions(false);
    }
    QMStorageMode getQMStorageMode() const { return qmStorageMode_m; }

    void setGlobalProcesses(std::vector<std::unique_ptr<GlobalProcess>> processes) {
        globalProcesses_m = std::move(processes);
    }

    const std::vector<std::unique_ptr<GlobalProcess>>& getGlobalProcesses() const {
        return globalProcesses_m;
    }

    /**
     * @brief Delete particles whose position is more than sigmasAway standard deviations
     *        from the bunch mean in any spatial dimension.
     *
     * Recomputes distribution moments, marks all particles outside the
     * [mean - sigmasAway*sigma, mean + sigmasAway*sigma] hyper-rectangle as
     * invalid, then calls the IPPL ParticleBase::destroy to compact the
     * particle arrays.
     *
     * @param sigmasAway Number of standard deviations defining the boundary.
     * @return Global number of particles destroyed (across all MPI ranks).
     */
    size_type deleteParticlesOutside(double sigmasAway) {
        size_type nLocal = this->getLocalNum();

        if (nLocal == 0 && this->getTotalNum() == 0) return 0;
        if (sigmasAway <= 0.0) return 0;

        // Force fresh moments: cached values from bunchUpdate may not reflect the exact
        // particle state at this point (emission/migration ordering can shift R between
        // the last bunchUpdate and this call). Safety-critical deletion always recomputes.
        markMomentsDirty();
        updateMoments();

        Vector_t<double, Dim> meanR = getMeanR();
        Vector_t<double, Dim> rmsR  = getRmsR();

        double lb0 = meanR[0] - sigmasAway * rmsR[0];
        double lb1 = meanR[1] - sigmasAway * rmsR[1];
        double lb2 = meanR[2] - sigmasAway * rmsR[2];
        double ub0 = meanR[0] + sigmasAway * rmsR[0];
        double ub1 = meanR[1] + sigmasAway * rmsR[1];
        double ub2 = meanR[2] + sigmasAway * rmsR[2];

        Kokkos::View<bool*> invalid("deleteParticlesOutside::invalid", nLocal);
        auto Rview = this->R.getView();

        size_type localDestroyNum = 0;
        Kokkos::parallel_reduce(
                "ParticleContainer::deleteParticlesOutside::mark", nLocal,
                KOKKOS_LAMBDA(const size_type i, size_type& count) {
                    bool outside = (Rview(i)[0] < lb0 || Rview(i)[0] > ub0)
                                   || (Rview(i)[1] < lb1 || Rview(i)[1] > ub1)
                                   || (Rview(i)[2] < lb2 || Rview(i)[2] > ub2);
                    invalid(i) = outside;
                    count += outside ? 1 : 0;
                },
                localDestroyNum);
        Kokkos::fence();

        size_type globalDestroyNum = 0;
        ippl::Comm->allreduce(localDestroyNum, globalDestroyNum, 1, std::plus<size_type>());

        if (globalDestroyNum == 0) return 0;

        destroyParticles(invalid, localDestroyNum);

        // Only called if globalDestroyNum > 0, i.e. if any particles were destroyed --> statistics
        // changed --> moments are dirty
        markMomentsDirty();

        return globalDestroyNum;
    }

    /**
     * @brief Create/allocate a specified number of particles.
     *
     * This function creates a given number of particles in the container. It's a wrapper around the
     * non destructive IPPL particle create function, but will print out a warning if the create
     * call led to unnecessary reallocation (i.e. if the new total number of particles exceeds the
     * previous capacity).
     *
     * @note The underlying `create` is a collective call, so all MPI ranks must call this function.
     * IPPL automatically handles the short-circuit if internal capacity is sufficient.
     *
     * @param numParticles The number of particles to create.
     */
    void createParticles(size_type numParticles) {
        Inform m("ParticleContainer::createParticles");

        // Total allocated capacity of the underlying view
        size_type oldCapacity = this->R.size();
        this->create(numParticles, true);  // non_destructive = true
        size_type newCapacity = this->R.size();

        // Pretty print numParticles, newCapacity and new totalNum + localNum after creation
        constexpr int labelWidth = 32;
        m << level4 << std::left << std::setw(labelWidth) << "Requested creation:" << numParticles
          << " particles" << endl
          << std::setw(labelWidth) << "New total number:" << this->getTotalNum()
          << " (local: " << this->getLocalNum() << ")" << endl
          << std::setw(labelWidth) << "Underlying view capacity:" << newCapacity << endl;

        if (newCapacity != oldCapacity) {
            m << level1
              << "WARNING: createParticles triggered a reallocation of the underlying particle "
                 "views! This can be a costly operation. To avoid this, consider increasing "
                 "preallocation (BEAM::NALLOC) or the overallocation factor."
              << endl;
        }
    }

    void allocateParticles(size_type numParticles) {
        Inform m("ParticleContainer::allocateParticles");

        // Total allocated capacity of the underlying view
        size_type oldCapacity = this->R.size();
        if (oldCapacity != 0) {
            throw OpalException(
                    "ParticleContainer::allocateParticles",
                    "Underlying views already allocated. This function is meant to be called on an "
                    "empty container, since it is destructive on existing particles. If you want "
                    "to create particles without deallocating existing ones, use createParticles() "
                    "instead.");
        }

        this->alloc(numParticles);  // alloc is always destructive

        m << level4 << std::left << std::setw(32) << "Requested allocation:" << numParticles
          << " particles" << endl
          << std::setw(32) << "Size of underlying view:" << this->R.size() << endl;
    }

    /**
     * @brief Destroy the particles marked invalid in this container.
     *
     * Wraps `ippl::ParticleBase::destroy` with input validation: throws if
     * `localDestroyNum` exceeds the local particle count, or if the `invalid`
     * mask is smaller than the local particle count. The underlying call is
     * collective (allreduce of `localNum_m`) so all MPI ranks must call this
     * function, even with `localDestroyNum == 0`.
     *
     * @note Does NOT mark moments dirty automatically. Callers that depend on
     * moment freshness must call `markMomentsDirty()` themselves.
     *
     * @tparam Properties Kokkos view properties of the invalid mask.
     * @param invalid Boolean mask of length >= getLocalNum(); true entries are removed.
     * @param localDestroyNum Number of true entries in `invalid` on this rank.
     */
    template <typename... Properties>
    void destroyParticles(
            const Kokkos::View<bool*, Properties...>& invalid, size_type localDestroyNum) {
        Inform m("ParticleContainer::destroyParticles");

        const size_type nLocal = this->getLocalNum();
        if (localDestroyNum > nLocal) {
            throw OpalException(
                    "ParticleContainer::destroyParticles",
                    "localDestroyNum (" + std::to_string(localDestroyNum)
                            + ") exceeds local particle count (" + std::to_string(nLocal) + ").");
        }
        if (invalid.extent(0) < nLocal) {
            throw OpalException(
                    "ParticleContainer::destroyParticles",
                    "invalid mask extent (" + std::to_string(invalid.extent(0))
                            + ") is smaller than local particle count ("
                            + std::to_string(nLocal) + ").");
        }

        this->destroy(invalid, localDestroyNum);

        constexpr int labelWidth = 32;
        m << level4 << std::left << std::setw(labelWidth)
          << "Requested destruction:" << localDestroyNum << " particles" << endl
          << std::setw(labelWidth) << "New total number:" << this->getTotalNum()
          << " (local: " << this->getLocalNum() << ")" << endl;
    }

private:
    void setBCAllPeriodic() { this->setParticleBC(ippl::BC::PERIODIC); }

    PLayout_t<T, Dim> pl_m;

    QMStorageMode qmStorageMode_m = QMStorageMode::SingleValue;

    DistributionMoments distMoments_m;

    /// Per-container state slot allocated by the handler at `setBunchStateHandler`.
    /// Owned here as the only strong reference; the handler keeps a weak_ptr, so
    /// destroying this container automatically releases the slot. The slot's own
    /// methods handle MPI consistency, so no direct handler reference is needed.
    std::shared_ptr<BunchStateHandler::ContainerState> containerState_m;

    // Single shared scalar mode stored as a length-1 Kokkos view.
    qm_view_type QView_m;
    qm_view_type MView_m;

    // Per-particle attributes mode
    ippl::ParticleAttrib<double> QAttr;
    ippl::ParticleAttrib<double> MAttr;

    // Reference particle information
    Vector_t<double, Dim> refPartR_m;
    Vector_t<double, Dim> refPartP_m;

    // Global to local quaternion
    Quaternion_t globalToLocalQuaternion_m;

    // Reference particle to lab transformation
    CoordinateSystemTrafo toLabTrafo_m;

    // Particle reference data (!= reference particle)
    const PartData* reference_m = nullptr;

    // Distance along the beamline
    double sPos_m = 0.0;

    /// Global physics processes attached to this container.
    std::vector<std::unique_ptr<GlobalProcess>> globalProcesses_m;
};

#endif
