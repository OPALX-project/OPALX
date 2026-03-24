
#ifndef OPAL_PARTICLE_CONTAINER_H
#define OPAL_PARTICLE_CONTAINER_H

// #include <functional>
#include <memory>
// #include <vector>

#include "Manager/BaseManager.h"

#include "Algorithms/DistributionMoments.h"

#include "Utilities/Options.h"

// #include <Kokkos_Core.hpp>

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

using size_type = ippl::detail::size_type;

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
    //typename Base::particle_position_type Etmp; // TODO: might not need this...

    /// magnetic field at particle position
    typename Base::particle_position_type B;

    ParticleContainer(Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL)
        : pl_m(FL, mesh),
          qmStorageMode_m(
              Options::useQMAttributes ? QMStorageMode::Attributes : QMStorageMode::SingleValue),
          distMoments_m(),
          QView_m("ParticleContainer::QView_m", 1),
          MView_m("ParticleContainer::MView_m", 1) {
        this->initialize(pl_m);
        registerAttributes();
        setupBCs();
        Kokkos::deep_copy(QView_m, 0.0);
        Kokkos::deep_copy(MView_m, 0.0);
    }

    ~ParticleContainer() {
    }

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

    void setupBCs() {
        setBCAllPeriodic();
    }

    /*!
     * Create nNew new particle slots. If current capacity is sufficient, calls base create.
     * Otherwise copies existing attribute data to temporaries, calls base create (which
     * reallocs), then copies data back so existing particles are preserved.
     */
    /*void createKeepData(size_type nNew) {
        const size_type required = this->getLocalNum() + nNew;
        if (this->R.getView().extent(0) >= required) {
            this->create(nNew);
            return;
        }
        const size_type nOld = this->getLocalNum();
        std::vector<std::function<void()>> copyBack;
        auto saveAttr = [&](auto& attr, const char* name) {
            auto& v = attr.getView();
            using VType = std::remove_reference_t<decltype(v)>;
            VType temp(std::string(name), nOld);
            Kokkos::deep_copy(
                temp,
                Kokkos::subview(v, std::make_pair(size_type(0), nOld)));
            copyBack.push_back([&attr, nOld, temp]() {
                Kokkos::deep_copy(
                    Kokkos::subview(attr.getView(),
                                   std::make_pair(size_type(0), nOld)),
                    temp);
            });
        };
        saveAttr(this->R, "createKeepData_R");
        saveAttr(this->ID, "createKeepData_ID");
        saveAttr(Q, "createKeepData_Q");
        saveAttr(M, "createKeepData_M");
        saveAttr(dt, "createKeepData_dt");
        saveAttr(Phi, "createKeepData_Phi");
        saveAttr(Bin, "createKeepData_Bin");
        saveAttr(P, "createKeepData_P");
        saveAttr(E, "createKeepData_E");
        saveAttr(B, "createKeepData_B");
        this->create(nNew);
        for (auto& f : copyBack) {
            f();
        }
    }*/

    PLayout_t<T, Dim>& getPL() {
        return pl_m;
    }

    void updateMoments(){
        size_t Np = this->getTotalNum();
        Np = (Np == 0) ? 1 : Np; // only used for normalization in the moments class --> avoid division by zero

        size_t Nlocal = this->getLocalNum();

        distMoments_m.computeMoments(this->R.getView(), this->P.getView(), getMView(), Np, Nlocal);
    }

    void setEnergyReferenceMass(double referenceMassGeV, bool rescaleToReference = true) {
        distMoments_m.setEnergyReferenceMass(referenceMassGeV, rescaleToReference);
    }

    Vector_t<double, 3> getMeanP() const{
         return distMoments_m.getMeanMomentum();
    }

    Vector_t<double, 3> getRmsP() const{
         return distMoments_m.getStandardDeviationMomentum();
    }

    Vector_t<double, 3> getMeanR() const{
         return distMoments_m.getMeanPosition();
    }

    Vector_t<double, 3> getRmsR() const{
         return distMoments_m.getStandardDeviationPosition();
    }

    Vector_t<double, 3> getRmsRP() const{
         return distMoments_m.getStandardDeviationRP();
    }

    void computeMinMaxR(){
        size_t Nlocal = this->getLocalNum();
        distMoments_m.computeMinMaxPosition(this->R.getView(), Nlocal);
    }

    Vector_t<double, 3> getMinR() const {
         return distMoments_m.getMinPosition();
    }

    Vector_t<double, 3> getMaxR() const {
         return distMoments_m.getMaxPosition();
    }

    matrix6x6_t getCovMatrix() const {
         return distMoments_m.getMoments6x6();
    }

    double getMeanKineticEnergy() const {
          return distMoments_m.getMeanKineticEnergy();
    }

    double getStdKineticEnergy() const {
          return distMoments_m.getStdKineticEnergy();
    }

    Vector_t<double, 6> getMeans() const {
        return distMoments_m.getMeans();
    }

    Vector_t<double, 6> getCentroid() const {
        return distMoments_m.getCentroid();
    }

    Vector_t<double, 3> getNormEmit() const {
        return distMoments_m.getNormalizedEmittance();
    }

    Vector_t<double, 3> getGeometricEmit() const {
        return distMoments_m.getGeometricEmittance();
    }

   double getDx() const {
       return distMoments_m.getDx();
   }

   double getDDx() const {
       return distMoments_m.getDDx();
   }

   double getDy() const {
       return distMoments_m.getDy();
   }

   double getDDy() const {
       return distMoments_m.getDDy();
   }

   double getDebyeLength() const {
       return distMoments_m.getDebyeLength();
   }

   double getMeanGammaZ() const {
       return distMoments_m.getMeanGammaZ();
   }

    double getTemperature() const {
        return distMoments_m.getTemperature();
    }

    double getPlasmaParameter() const {
        return distMoments_m.getPlasmaParameter();
    }

    double computeDebyeLength(double density){
        size_t Np = this->getTotalNum();
        Np = (Np == 0) ? 1 : Np; // only used for normalization in the moments class --> avoid division by zero

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
            auto view = QAttr.getView();
            const size_type n = view.extent(0);
            using exec_space = typename ippl::ParticleAttrib<double>::execution_space;
            Kokkos::parallel_for(
                "ParticleContainer::setQ", Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { view(i) = q; });
            Kokkos::fence();
        } else {
            Kokkos::deep_copy(QView_m, q);
        }
    }

    /**
     * @brief Set particle mass for the active M storage mode.
     * @param m Mass value in [GeV].
     *
     * In `QMStorageMode::Attributes`, this assigns `m` to every local particle.
     * In `QMStorageMode::SingleValue`, this updates the shared scalar mass view.
     */
    void setM(double m) {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = MAttr.getView();
            const size_type n = view.extent(0);
            using exec_space = typename ippl::ParticleAttrib<double>::execution_space;
            Kokkos::parallel_for(
                "ParticleContainer::setM", Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { view(i) = m; });
            Kokkos::fence();
        } else {
            Kokkos::deep_copy(MView_m, m);
        }
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
        auto dtView = dt.getView();
        const size_type n = dtView.extent(0);
        if (n == 0) {
            return;
        }

        using exec_space = typename ippl::ParticleAttrib<double>::execution_space;
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto qView = QAttr.getView();
            Kokkos::parallel_for(
                "ParticleContainer::scaleDtByCharge(attrs)",
                Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { dtView(i) *= qView(i); });
        } else {
            auto QView = QView_m;
            Kokkos::parallel_for(
                "ParticleContainer::scaleDtByCharge(single)",
                Kokkos::RangePolicy<exec_space>(0, n),
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
        auto dtView      = dt.getView();
        const size_type n = dtView.extent(0);
        if (n == 0) {
            return;
        }

        using exec_space = typename ippl::ParticleAttrib<double>::execution_space;
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto qView = QAttr.getView();
            Kokkos::parallel_for(
                "ParticleContainer::unscaleDtByCharge(attrs)",
                Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { dtView(i) /= qView(i); });
        } else {
            auto QView = QView_m;
            Kokkos::parallel_for(
                "ParticleContainer::unscaleDtByCharge(single)",
                Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { dtView(i) /= QView(0); });
        }
        Kokkos::fence();
    }

    QMStorageMode getQMStorageMode() const {
        return qmStorageMode_m;
    }

private:
    void setBCAllPeriodic() {
        this->setParticleBC(ippl::BC::PERIODIC);
    }

    PLayout_t<T, Dim> pl_m;

    QMStorageMode qmStorageMode_m = QMStorageMode::SingleValue;

    DistributionMoments distMoments_m;

    // Single shared scalar mode stored as a length-1 Kokkos view.
    qm_view_type QView_m;
    qm_view_type MView_m;

    // Per-particle attributes mode
    ippl::ParticleAttrib<double> QAttr;
    ippl::ParticleAttrib<double> MAttr;
};

#endif
