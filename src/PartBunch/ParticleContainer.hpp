
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
        distMoments_m.computeMoments(this->R.getView(), this->P.getView(), Np);
    }

    Vector_t<double, 3> get_pmean(){
         return distMoments_m.getMeanMomentum();
    }

    void computeMinMaxR(){
        distMoments_m.computeMinMaxPosition(this->R.getView());
    }

    Vector_t<double, 3> getMinR() const {
         return distMoments_m.getMinPosition();
    };

    Vector_t<double, 3> getMaxR() const {
         return distMoments_m.getMaxPosition();
    }

private:
    void setBCAllPeriodic() {
        this->setParticleBC(ippl::BC::PERIODIC);
    }

    PLayout_t<T, Dim> pl_m;

    DistributionMoments distMoments_m;
};

#endif
