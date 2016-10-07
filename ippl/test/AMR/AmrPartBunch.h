#ifndef AMRPARTBUNCH_H
#define AMRPARTBUNCH_H


#include <map>
#include <tuple>

#include "PartBunchBase.h"

#include <Particles.H>

class AmrParticleContainer : public ParticleContainer<7, 0>
{};


class AmrPartBunch : public PartBunchBase,
                     public AmrParticleContainer
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
    
    virtual ~AmrPartBunch() {}
    
    //BEGIN PartBunchBase functions
    void myUpdate() {
        
    }
    
    
    void create(int nloc) {
        
    }
    
    void gatherStatistics() {
        
    }
    
    size_t getLocalNum() const {
        return 0;
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
    
// private:
//     void buildIndexMapping_m()
    
};



// void AmrPartBunch::buildIndexMapping_m() {
//     int i = 0;
//     for (int l = 0; l < nLevels; ++l) {
//         for (int g = 0; g < nGrids; ++g) {  
//             for (int d = 0; d < dqs; ++d) {
//                 idxMap_m[i] = std::tuple<int, int, int>(l,g,d);
//                 ++i;
//             }
//         }
//     }
// }



#endif