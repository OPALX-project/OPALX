
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

// Define the ParticlesContainer class
template <typename T, unsigned Dim = 3>
class ParticleContainer : public ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>> {
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

public:
    enum class QMStorageMode { SingleValue, Attributes };

    /// Defines which type to use as a particle bin.
    using bin_index_type = short int;  // Needed in AdaptBins class

public:
    /// Charge (single value) in [Cb] when using `QMStorageMode::SingleValue`.
    /// When using `QMStorageMode::Attributes`, `Q(i)` returns per-particle values.
    double getQScalar() const {
        return Q_m;
    }

    /// Mass (single value) in [GeV] when using `QMStorageMode::SingleValue`.
    /// When using `QMStorageMode::Attributes`, `M(i)` returns per-particle values.
    double getMScalar() const {
        return M_m;
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
          distMoments_m() {
        this->initialize(pl_m);
        registerAttributes();
        setupBCs();
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
            std::cout << "ATTRIBUTE MODE" << std::endl;
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

        if (qmStorageMode_m == QMStorageMode::Attributes) {
            distMoments_m.computeMoments(this->R.getView(), this->P.getView(), this->MAttr.getView(),
                                         Np, Nlocal);
        } else {
            distMoments_m.computeMoments(this->R.getView(), this->P.getView(), M_m, Np, Nlocal);
        }
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

    /// Per-particle charge accessor. In `SingleValue` mode this returns `getQScalar()`.
    double Q(size_type i) const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return QAttr(i);
        }
        return Q_m;
    }

    /// Per-particle mass accessor. In `SingleValue` mode this returns `getMScalar()`.
    double M(size_type i) const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return MAttr(i);
        }
        return M_m;
    }

    void setQ(double q) {
        Q_m = q;
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = QAttr.getView();
            const size_type n = view.extent(0);
            using exec_space = typename ippl::ParticleAttrib<double>::execution_space;
            Kokkos::parallel_for(
                "ParticleContainer::setQ", Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { view(i) = q; });
            Kokkos::fence();
        }
    }

    void setM(double m) {
        M_m = m;
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = MAttr.getView();
            const size_type n = view.extent(0);
            using exec_space = typename ippl::ParticleAttrib<double>::execution_space;
            Kokkos::parallel_for(
                "ParticleContainer::setM", Kokkos::RangePolicy<exec_space>(0, n),
                KOKKOS_LAMBDA(const size_type i) { view(i) = m; });
            Kokkos::fence();
        }
    }

    /// Scale `dt[i]` by `Q(i)` so that `scatter(*dt, ...)` uses the correct per-particle charge.
    void scaleDtByCharge() {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            dt = dt * QAttr;
        } else {
            dt = dt * Q_m;
        }
    }

    /// Undo `scaleDtByCharge()` after scattering.
    void unscaleDtByCharge() {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            dt = dt / QAttr;
        } else {
            dt = dt / Q_m;
        }
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

    // Single shared scalar mode
    double Q_m = 0.0;
    double M_m = 0.0;

    // Per-particle attributes mode
    ippl::ParticleAttrib<double> QAttr;
    ippl::ParticleAttrib<double> MAttr;

    DistributionMoments distMoments_m;
};

#endif
