#ifndef AMRPARTBUNCH_H
#define AMRPARTBUNCH_H


#include <map>
#include <tuple>

#include "PartBunchBase.h"

#include <Particles.H>


class AmrPartBunch : public PartBunchBase,
                    public ParticleContainer<7, 0>
//                      public AmrParticleContainer
{
private:
    /* Mapping for
     * index --> (level, grid, deque length)
     * and vice-versa
     */
    std::map<int, std::tuple<int, int, int> > idxMap_m;
    
    
    Vector_t R, P, E, B;
    double qm;

public:
    
    AmrPartBunch(const Geometry            & geom, 
                 const DistributionMapping & dmap,
                 const BoxArray            & ba) : ParticleContainer(geom, dmap, ba)
                 {}
    
    
    virtual ~AmrPartBunch() {}
    
    void InitRandom (long icount, unsigned long iseed, Real particleMass, bool serialize = false, RealBox bx = RealBox()) {
        ParticleContainer<7, 0>::InitRandom (icount, iseed, particleMass, serialize, bx);
    }
    
    
    //BEGIN PartBunchBase functions
    void myUpdate() {
        buildIndexMapping_m();
        
        
        for (int i=0; i< nLocalParticles; i++) {   
            int l,g,dq;
            std::tie(l,g,dq) = idxMap_m[i];
            std::cout << "i= " << i << " ---> ";
            std::cout << "level= " << l << " grid= " << g << " dqIndex= " << dq;
            std::cout << std::endl;
        }
        
    }
    
    
    void create(int nloc) {
        
    }
    
    void gatherStatistics() {
        
    }
    
    size_t getLocalNum() const {
        return nLocalParticles;
    }
    
    Vector_t getRMin() {
        return Vector_t(0.0, 0.0, 0.0);
    }
    
    Vector_t getRMax() {
        return Vector_t(0.0, 0.0, 0.0);
    }
    
    Vector_t getHr() {
        return Vector_t(0.0, 0.0, 0.0);
    }
    
    
    double scatter() {
        return 0;
    }
    
    void initFields() {}
    
    void gatherCIC() {}
    
    
    Vector_t& getR(int i) {
        int l, g, dq;
        std::tie(l,g,dq) = idxMap_m[i];
        R = Vector_t(m_particles[l][g][dq].m_pos[0],
                     m_particles[l][g][dq].m_pos[1],
                     m_particles[l][g][dq].m_pos[2]);
        
        return R;
    }
    
    double& getQM(int i) {
        return qm;
    }
    
    Vector_t& getP(int i) {
        return P;
    }
    
    Vector_t& getE(int i) {
        return E;
    }
    
    Vector_t& getB(int i) {
        return B;
    }
    //END PartBunchBase functions
    
private:
    void buildIndexMapping_m();
    int nLocalParticles;
    
};



void AmrPartBunch::buildIndexMapping_m() {
    int i = 0;
    nLocalParticles = 0;
    for (int l = 0; l < m_particles.size(); ++l) {
        for (unsigned int g = 0; g <  m_particles[l].size(); ++g) {
            
            nLocalParticles += m_particles[l][g].size();
            
            for (unsigned int d = 0; d < m_particles[l][g].size(); ++d) {
                idxMap_m[i] = std::tuple<int, int, int>(l,g,d);
                ++i;
            }
        }
    }
}



#endif