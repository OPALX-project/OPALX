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
 * @brief Amr particle bunch class
 */


#include <map>
#include <tuple>

// #include "Ippl.h"

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
    
    
    void destroyAll() {
        for (unsigned int i = 0; i < m_particles.size(); ++i)
            this->RemoveParticlesAtLevel(i);
        nLocalParticles_m = 0;
    }
    
//     void setParGDB(AmrParGDB* gdb) {
//         Define(gdb);
//     }

    void python_format(int step) {
        std::string st = std::to_string(step);
        Inform::WriteMode wm = Inform::OVERWRITE;
        for (int i = 0; i < Ippl::getNodes(); ++i) {
            if ( i == Ippl::myNode() ) {
                wm = (i == 0) ? wm : Inform::APPEND;
                
                const ParGDBBase* gdb = GetParGDB();
                
                std::string grid_file = "pyplot_grids_" + st + ".dat";
                Inform msg("", grid_file.c_str(), wm, i);
                for (int l = 0; l < gdb->finestLevel() + 1; ++l) {
                    Geometry geom = gdb->Geom(l);
                    for (int g = 0; g < gdb->ParticleBoxArray(l).size(); ++g) {
                        msg << l << ' ';
                        RealBox loc = RealBox(gdb->boxArray(l)[g],geom.CellSize(),geom.ProbLo());
                        for (int n = 0; n < BL_SPACEDIM; n++)
                            msg << loc.lo(n) << ' ' << loc.hi(n) << ' ';
                        msg << endl;
                    }
                }
                
                std::string particle_file = "pyplot_particles_" + st + ".dat";
                Inform msg2all("", particle_file.c_str(), wm, i);
                int l, g, dq;
                for (size_t i = 0; i < this->getLocalNum(); ++i) {
                    
                    std::tie(l,g,dq) = idxMap_m[i];
                    
                    msg2all << m_particles[l][g][dq].m_pos[0] << " "
                            << m_particles[l][g][dq].m_pos[1] << " "
                            << m_particles[l][g][dq].m_pos[2]
                            << endl;
                }
            }
            Ippl::Comm->barrier();
        }
    }
    
    Vector_t interpolate(int i, MultiFab& quantity) {
        
        int lev, grid, dq;
        std::tie(lev, grid, dq) = idxMap_m[i];
        
        const Real strttime  = ParallelDescriptor::second();
        
        MultiFab* ac_pointer;
        if (OnSameGrids(lev,quantity)) {
            ac_pointer = &quantity;
        } else {
            ac_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),quantity.nComp(),quantity.nGrow(),
                                      m_gdb->ParticleDistributionMap(lev),Fab_allocate);
            for (MFIter mfi(*ac_pointer); mfi.isValid(); ++mfi)
                ac_pointer->setVal(0.);
            ac_pointer->copy(quantity,0,0,quantity.nComp());
            ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
        }
        
        const FArrayBox& gfab = (*ac_pointer)[grid];
        Real grav[1] = { 0.0 };
        
        int idx[1] = { 0 };

    ParticleBase::Interp(m_particles[lev][grid][dq],
                         m_gdb->Geom(lev),
                         gfab,
                         idx,
                         grav,
                         1);
        
        return Vector_t(grav[0], 0, 0);
    }
    
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