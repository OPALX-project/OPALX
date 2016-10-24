#ifndef AMRPARTBUNCH_H
#define AMRPARTBUNCH_H

/*!
 * @file AmrPartBunch.h
 * @details Particle bunch class
 * for BoxLib
 * @authors Matthias Frey \n
 *          Andreas Adelmann \n
 *          Ann Almgren \n
 *          Weiqun Zhang
 * @date LBNL, October 2016
 */


#include <map>
#include <tuple>

#include "PartBunchBase.h"

#include <AmrParGDB.H>

#include <Particles.H>


/// Particle bunch class for BoxLib
class AmrPartBunch : public PartBunchBase,
                     public ParticleContainer<10 /*real attributes*/, 0>
{
public:
    typedef std::map<int, std::tuple<int, int, int> > map_t;
    
    static size_t nAttributes;

public:
    
    /// Just calls constructor of ParticleContainer
    AmrPartBunch(const Geometry & geom, const DistributionMapping & dmap,
                 const BoxArray & ba);
    
    /// Just calls constructor of ParticleContainer
    AmrPartBunch(const Array<Geometry>& geom,
                 const Array<DistributionMapping>& dmap,
                 const Array<BoxArray>& ba,
                 const Array<int>& rr);
                 
    
    /// Does nothing.
    virtual ~AmrPartBunch();
    
    // ------------------------------------------------------------------------
    // INHERITED MEMBER FUNCTIONS
    // ------------------------------------------------------------------------
    
    void myUpdate();
    
    void create(size_t m);
    
    void gatherStatistics();
    
    size_t getLocalNum() const;
    
    size_t getTotalNum() const;
    
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
    
    inline Vector_t getR(int i);
    
    inline double getQM(int i);
    
    inline Vector_t getP(int i);
    
    inline Vector_t getE(int i);
    
    inline Vector_t getB(int i);
    
    inline void setR(Vector_t pos, int i);
    
    inline void setQM(double q, int i);
    
    inline void setP(Vector_t v, int i);
    
    inline void setE(Vector_t Ef, int i);
    
    inline void setB(Vector_t Bf, int i);
    
//     void setParGDB(AmrParGDB* gdb) {
//         Define(gdb);
//     }
    
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
};

#endif