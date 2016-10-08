#ifndef AMRPARTBUNCH_H
#define AMRPARTBUNCH_H


#include <map>
#include <tuple>

#include "PartBunchBase.h"

#include <Particles.H>


/// Particle bunch class for BoxLib
class AmrPartBunch : public PartBunchBase,
                    public ParticleContainer<7, 0>
{
public:
    typedef std::map<int, std::tuple<int, int, int> > map_t;

public:
    
    /// Just calls constructor of ParticleContainer
    AmrPartBunch(const Geometry & geom, const DistributionMapping & dmap,
                 const BoxArray & ba);
    
    /// Does nothing.
    virtual ~AmrPartBunch();
    
    // ------------------------------------------------------------------------
    // INHERITED MEMBER FUNCTIONS
    // ------------------------------------------------------------------------
    
    void myUpdate();
    
    
    void create(size_t m);
    
    ///@todo Implement
    void gatherStatistics();
    
    size_t getLocalNum() const;
    
    ///@todo Implement
    Vector_t getRMin();
    
    ///@todo Implement
    Vector_t getRMax();
    
    ///@todo Implement
    Vector_t getHr();
    
    ///@todo Implement
    double scatter();
    
    ///@todo Implement
    void initFields();
    
    ///@todo Implement
    void gatherCIC();
    
    Vector_t& getR(int i);
    
    ///@todo
    double& getQM(int i);
    
    ///@todo
    Vector_t& getP(int i);
    
    ///@todo
    Vector_t& getE(int i);
    
    ///@todo
    Vector_t& getB(int i);
    
    inline void setR(Vector_t pos, int i);
    
    
private:
    /// Create the index mapping in order to have random access
    void buildIndexMapping_m();
    
    int nLocalParticles_m;
    
private:
    /* Mapping for
     * index --> (level, grid, deque length)
     * and vice-versa
     */
    map_t idxMap_m;
    
    // remove
    Vector_t R, P, E, B;
    double qm;
};



// ----------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// ----------------------------------------------------------------------------

AmrPartBunch::AmrPartBunch(const Geometry& geom,
                           const DistributionMapping & dmap,
                           const BoxArray & ba)
    : ParticleContainer(geom, dmap, ba)
{ }


AmrPartBunch::~AmrPartBunch()
{ }

// ----------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// ----------------------------------------------------------------------------

void AmrPartBunch::buildIndexMapping_m() {
    
//     idxMap_m.clear();
    
    int i = 0;
    nLocalParticles_m = 0;
    for (int l = 0; l < m_particles.size(); ++l) {
        for (unsigned int g = 0; g <  m_particles[l].size(); ++g) {
            
            nLocalParticles_m += m_particles[l][g].size();
            
            for (unsigned int d = 0; d < m_particles[l][g].size(); ++d) {
                idxMap_m[i] = std::tuple<int, int, int>(l,g,d);
                ++i;
            }
        }
    }
}


// ----------------------------------------------------------------------------
// INHERITED MEMBER FUNCTIONS
// ----------------------------------------------------------------------------


void AmrPartBunch::myUpdate() {
    buildIndexMapping_m();
    
//     for (int i = 0; i < nLocalParticles_m; i++) {
//         int l,g,dq;
//         std::tie(l,g,dq) = idxMap_m[i];
//         std::cout << "i= " << i << " ---> ";
//         std::cout << "level= " << l << " grid= " << g << " dqIndex= " << dq;
//         std::cout << std::endl;
//     }
}


void AmrPartBunch::create(size_t m) {
    // dummy values, just to reserve space
    // InitRandom(total count, seed, particle mass, serialize)
    ParticleContainer<7, 0>::InitRandom (m * ParallelDescriptor::NProcs(), 42, 1.0, false);
    
    buildIndexMapping_m();
}


void AmrPartBunch::gatherStatistics() { }


size_t AmrPartBunch::getLocalNum() const {
    return nLocalParticles_m;
}


Vector_t AmrPartBunch::getRMin() {
    return Vector_t(0.0, 0.0, 0.0);
}


Vector_t AmrPartBunch::getRMax() {
    return Vector_t(0.0, 0.0, 0.0);
}


Vector_t AmrPartBunch::getHr() {
    return Vector_t(0.0, 0.0, 0.0);
}

double AmrPartBunch::scatter() {
    return 0;
}


void AmrPartBunch::initFields() {}


void AmrPartBunch::gatherCIC() {}


Vector_t& AmrPartBunch::getR(int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    R = Vector_t(m_particles[l][g][dq].m_pos[0],
                 m_particles[l][g][dq].m_pos[1],
                 m_particles[l][g][dq].m_pos[2]);
    
    return R;
}


double& AmrPartBunch::getQM(int i) {
    return qm;
}


Vector_t& AmrPartBunch::getP(int i) {
    return P;
}


Vector_t& AmrPartBunch::getE(int i) {
    return E;
}


Vector_t& AmrPartBunch::getB(int i) {
    return B;
}


void AmrPartBunch::setR(Vector_t pos, int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    
    for (int d = 0; d < 3; ++d)
        m_particles[l][g][dq].m_pos[d] = pos(d);
}

#endif