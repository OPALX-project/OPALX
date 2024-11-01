
#ifndef OPAL_PARTICLE_CONTAINER_H
#define OPAL_PARTICLE_CONTAINER_H

#include <memory>

#include "Manager/BaseManager.h"

#include "Algorithms/DistributionMoments.h"

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

using size_type = ippl::detail::size_type;

// Define the ParticlesContainer class
template <typename T, unsigned Dim = 3>
class ParticleContainer : public ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>> {
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

public:
    /// charge in [Cb]
    ippl::ParticleAttrib<double> Q;

    /// mass
    ippl::ParticleAttrib<double> M;

    /// timestep in [s]
    ippl::ParticleAttrib<double> dt;

    // time [s]
    ippl::ParticleAttrib<ippl::Vector<double, 1>> t;

    /// the scalar potential in [Cb/s]
    ippl::ParticleAttrib<double> Phi;

    /// the energy bin the particle is in
    ippl::ParticleAttrib<short int> Bin;

    /// the particle specis
    ippl::ParticleAttrib<short> Sp;

    /// particle momenta [\beta\gamma]
    typename Base::particle_position_type P;

    /// electric field at particle position
    typename Base::particle_position_type E;

    /// electric field for gun simulation with bins
    typename Base::particle_position_type Etmp;

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
        this->addAttribute(t);
        this->addAttribute(Phi);
        this->addAttribute(Bin);
        this->addAttribute(Sp);
        this->addAttribute(P);
        this->addAttribute(E);
        this->addAttribute(Etmp);
        this->addAttribute(B);
    }
    void setupBCs() {
        setBCAllPeriodic();
    }

    PLayout_t<T, Dim>& getPL() {
        return pl_m;
    }

    void updateMoments(){
        size_t Np = this->getTotalNum();
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

    matrix_t getCovMatrix() const {
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
        size_t Nlocal = this->getLocalNum();
        distMoments_m.computeDebyeLength(this->R.getView(), this->P.getView(), Np, Nlocal, density);
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
