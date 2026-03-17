
#ifndef OPAL_PARTICLE_CONTAINER_H
#define OPAL_PARTICLE_CONTAINER_H

// #include <functional>
#include <memory>
// #include <vector>

#include "Manager/BaseManager.h"

#include "Algorithms/DistributionMoments.h"

// #include <Kokkos_Core.hpp>

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

using size_type = ippl::detail::size_type;

// Define the ParticlesContainer class
template <typename T, unsigned Dim = 3>
class ParticleContainer : public ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>> {
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

public:
    /// Defines which type to use as a particle bin.
    using bin_index_type = short int;  // Needed in AdaptBins class

public:
    /// charge in [Cb]
    ippl::ParticleAttrib<double> Q;

    /// mass
    ippl::ParticleAttrib<double> M;

    /// timestep in [s]
    ippl::ParticleAttrib<double> dt;

    /// the scalar potential in [Cb/s]
    ippl::ParticleAttrib<double> Phi;

    /// the energy bin the particle is in
    ippl::ParticleAttrib<bin_index_type> Bin;

    /// the particle specis
    ippl::ParticleAttrib<short> Sp;

    /// particle momenta [\beta\gamma]
    typename Base::particle_position_type P;

    /// electric field at particle position
    typename Base::particle_position_type E;

    /// electric field for gun simulation with bins
    //typename Base::particle_position_type Etmp; // TODO: might not need this...

    /// magnetic field at particle position
    typename Base::particle_position_type B;

    ParticleContainer(Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL) : pl_m(FL, mesh), distMoments_m() {
        this->initialize(pl_m);
        registerAttributes();
        setupBCs();
    }

    ~ParticleContainer() {
    }

    void registerAttributes() {
        // register the particle attributes
        this->addAttribute(Q);
        this->addAttribute(M);
        this->addAttribute(dt);
        this->addAttribute(Phi);
        this->addAttribute(Bin);
        this->addAttribute(Sp);
        this->addAttribute(P);
        this->addAttribute(E);
        //this->addAttribute(Etmp);
        this->addAttribute(B);
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
        saveAttr(Sp, "createKeepData_Sp");
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
        distMoments_m.computeMoments(this->R.getView(), this->P.getView(), this->M.getView(), Np, Nlocal);
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

private:
    void setBCAllPeriodic() {
        this->setParticleBC(ippl::BC::PERIODIC);
    }

    PLayout_t<T, Dim> pl_m;

    DistributionMoments distMoments_m;
};

#endif
