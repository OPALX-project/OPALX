#ifndef PART_BUNCH_BASE_H
#define PART_BUNCH_BASE_H

#include "Particle/AbstractParticle.h"

#include "Algorithms/PBunchDefs.h"

template <class T, unsigned Dim>
class PartBunchBase
{
public:
    typedef typename AbstractParticle<T, Dim>::ParticlePos_t ParticlePos_t;
    typedef typename AbstractParticle<T, Dim>::ParticleIndex_t ParticleIndex_t;
    typedef typename AbstractParticle<T, Dim>::UpdateFlags UpdateFlags;
    typedef typename AbstractParticle<T, Dim>::Position_t Position_t;
    
    static const unsigned Dimension = Dim;
    
    PartBunchBase(AbstractParticle<T, Dim>* pb) : pbase(pb),
                                                  R(*(pb->R_p)),
                                                  ID(*(pb->ID_p))
    {
        pb->addAttribute(P);
        pb->addAttribute(Q);
        pb->addAttribute(M);
        pb->addAttribute(Phi);
        pb->addAttribute(Ef);
        pb->addAttribute(Eftmp);
    
        pb->addAttribute(Bf);
        pb->addAttribute(Bin);
        pb->addAttribute(dt);
        pb->addAttribute(PType);
        pb->addAttribute(TriID);
    }
    
    /*
     * Wrapped member functions of IpplParticleBase
     */
    
    size_t getTotalNum() const {
        return pbase->getTotalNum();
    }
        
    size_t getLocalNum() const {
        return pbase->getLocalNum();
    }
    
    
    size_t getDestroyNum() const {
        return pbase->getDestroyNum();
    }
    
    size_t getGhostNum() const {
        return pbase->getGhostNum();
    }
    
    void setTotalNum(size_t n) {
        pbase->setTotalNum(n);
    }
    
    void setLocalNum(size_t n) {
        pbase->setLocalNum(n);
    }
    
    bool getUpdateFlag(UpdateFlags f) const {
        return pbase->getUpdateFlag(f);
    }
    
    void setUpdateFlag(UpdateFlags f, bool val) {
        pbase->setUpdateFlag(f, val);
    }
    
    ParticleBConds<Position_t, Dimension>& getBConds() {
        return pbase->getBConds();
    }
    
    void setBConds(const ParticleBConds<Position_t, Dimension>& bc) {
        pbase->setBConds(bc);
    }
    
    bool singleInitNode() const {
        return pbase->singleInitNode();
    }
    
    void resetID() {
        pbase->resetID();
    }
    
    void update() {
        pbase->update();
    }
    
    void update(const ParticleAttrib<char>& canSwap) {
        pbase->update(canSwap);
    }
    
    void createWithID(unsigned id) {
        pbase->createWithID(id);
    }
    
    void create(size_t M) {
        pbase->create(M);
    }
    
    void globalCreate(size_t np) {
        pbase->globalCreate(np);
    }
    
    void destroy(size_t M, size_t I, bool doNow = false) {
        pbase->destroy(M, I, doNow);
    }
    
    void ghostDestroy(size_t M, size_t I) {
        pbase->ghostDestroy(M, I);
    }
    
    /*
     * (Pure) virtual member functions 
     */
    
    virtual Mesh_t &getMesh() = 0;

    virtual FieldLayout_t &getFieldLayout() = 0;
    
    
    
    /*
     * Bunch attributes
     */
    
    
    ParticlePos_t& R;
    ParticleIndex_t& ID;
    
    
    
    // Particle container attributes
    ParticleAttrib< Vector_t > P;      // particle momentum //  ParticleSpatialLayout<double, 3>::ParticlePos_t P;
    ParticleAttrib< double >   Q;      // charge per simulation particle, unit: C.
    ParticleAttrib< double >   M;      // mass per simulation particle, for multi-species particle tracking, unit:GeV/c^2.
    ParticleAttrib< double >   Phi;    // the electric potential
    ParticleAttrib< Vector_t > Ef;     // e field vector
    ParticleAttrib< Vector_t > Eftmp;  // e field vector for gun simulations

    ParticleAttrib< Vector_t > Bf;   // b field vector
    ParticleAttrib< int >      Bin;   // holds the bin in which the particle is in, if zero particle is marked for deletion
    ParticleAttrib< double >   dt;   // holds the dt timestep for particle

    ParticleAttrib< short >    PType; // we can distinguish dark current particles from primary particle
    ParticleAttrib< int >      TriID; // holds the ID of triangle that the particle hit. Only for BoundaryGeometry case.
    
// protected:
    /*std::unique_ptr<*/AbstractParticle<T, Dim>*/* >*/ pbase;
};

#endif
