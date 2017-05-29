#ifndef BOXLIB_LAYOUT_H
#define BOXLIB_LAYOUT_H

#include "AmrParticle/ParticleAmrLayout.h"
#include "AmrParticle/AmrParticleBase.h"

#include "Amr/AmrDefs.h"

#include <AMReX_ParGDB.H>

template<class T, unsigned Dim>
class BoxLibLayout : public ParticleAmrLayout<T, Dim>,
                     public amrex::ParGDB
{
    
public:
    typedef typename ParticleAmrLayout<T, Dim>::pair_t pair_t;
    typedef typename ParticleAmrLayout<T, Dim>::pair_iterator pair_iterator;
    typedef typename ParticleAmrLayout<T, Dim>::SingleParticlePos_t SingleParticlePos_t;
    typedef typename ParticleAmrLayout<T, Dim>::Index_t Index_t;
    
    typedef amr::AmrField_t AmrField_t;
    typedef amr::AmrFieldContainer_t AmrFieldContainer_t;
    typedef typename ParticleAmrLayout<T, Dim>::ParticlePos_t ParticlePos_t;
    typedef ParticleAttrib<Index_t> ParticleIndex_t;
    
    typedef amr::AmrProcMap_t           AmrProcMap_t;
    typedef amr::AmrGrid_t              AmrGrid_t;
    typedef amr::AmrGeometry_t          AmrGeometry_t;
    typedef amr::AmrGeomContainer_t     AmrGeomContainer_t;
    typedef amr::AmrGridContainer_t     AmrGridContainer_t;
    typedef amr::AmrProcMapContainer_t  AmrProcMapContainer_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    typedef amr::AmrIntVectContainer_t  AmrIntVectContainer_t;
    typedef amr::AmrIntArray_t          AmrIntArray_t;
    typedef amr::AmrDomain_t            AmrDomain_t;
    typedef amr::AmrBox_t               AmrBox_t;
    typedef amr::AmrReal_t              AmrReal_t;
    
    static const Vector_t lowerBound;
    static const Vector_t upperBound;
    
//     static bool do_tiling;
//     static AmrIntVect_t tile_size;

public:
    
    /*!
     * Initializes default Geometry, DistributionMapping and BoxArray.
     */
    BoxLibLayout();
    
    BoxLibLayout(int nGridPoints, int maxGridSize);
    
    BoxLibLayout(const AmrGeometry_t &geom,
                 const AmrProcMap_t &dmap,
                 const AmrGrid_t &ba);
    
    BoxLibLayout(const AmrGeomContainer_t &geom,
                 const AmrProcMapContainer_t &dmap,
                 const AmrGridContainer_t &ba,
                 const AmrIntArray_t &rr);
    
    void update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData, const ParticleAttrib<char>* canSwap = 0);
    
    void update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData, int lev_min = 0,
                const ParticleAttrib<char>* canSwap = 0);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    AmrIntVect_t Index (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                        const unsigned int ip, int level) const;

    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    AmrIntVect_t Index (SingleParticlePos_t &R, int lev) const;
    
    
    void initDefaultBox(int nGridPoints, int maxGridSize);
    
    void resize(int maxLevel) {
        int length = maxLevel + 1;
        this->m_geom.resize(length);
        this->m_dmap.resize(length);
        this->m_ba.resize(length);
        this->m_nlevels = length;
        this->refRatio_m.resize(maxLevel);
//         this->m_rr.resize(maxLevel);
        maxLevel_m = maxLevel;
        
        for (int i = 0; i < length; ++i) {
            std::cout << "resize BA: " << this->m_ba[i] << std::endl;
            
        }
    }
    
//     void define(const AmrGeomContainer_t& geom,
//                 const AmrGridContainer_t& ba,
//                 const AmrProcMapContainer_t& dmap,
//                 const AmrIntArray_t & rr)
//     {
//         maxLevel_m = ba.size() - 1;
//         this->m_nlevels = ba.size();
//         this->m_geom.resize( this->m_nlevels );
//         this->m_ba.resize( this->m_nlevels );
//         this->m_dmap.resize( this->m_nlevels );
//         this->m_rr.resize(maxLevel_m);
//         
//         for (int i = 0; i < this->m_nlevels; ++i) {
//             this->m_geom[i] = geom[i];
//             this->m_ba[i]   = ba[i];
//             this->m_dmap[i] = dmap[i];
//             this->m_rr[i]   = rr[i];
//         }
//     }
    
    void define(const AmrGeomContainer_t& geom) {
//         this->m_geom.resize( geom.size() );
        std::cout << "this->m_geom.size() = " << this->m_geom.size() << std::endl;
        std::cout << "geom.size() = " << geom.size() << std::endl;
        for (unsigned int i = 0; i < geom.size(); ++i)
            this->m_geom[i] = geom[i];
    }
    
    void define(const AmrIntVectContainer_t& refRatio) {
        for (unsigned int i = 0; i < refRatio.size(); ++i) {
            refRatio_m[i] = refRatio[i];
        }
    }
    
    
    void define(const AmrGeomContainer_t& geom,
                const AmrGridContainer_t& ba,
                const AmrProcMapContainer_t& dmap)
    {
        std::cout << "define: " << std::endl;
        std::cout << geom.size() << " " << ba.size() << " " << dmap.size() << std::endl;
        
        for (int i = 0; i < this->m_nlevels; ++i) {
            std::cout << "AMR BA: " << ba[i] << " " << ba[i].ok() << std::endl;
        }
        
        std::cout << this->m_geom.size() << " " << this->m_ba.size() << " " << this->m_dmap.size() << std::endl;
        
        maxLevel_m = ba.size() - 1;
        this->m_nlevels = ba.size();
        this->m_geom.resize( this->m_nlevels );
        this->m_ba.resize( this->m_nlevels );
        this->m_dmap.resize( this->m_nlevels );
        
        for (int i = 0; i < this->m_nlevels; ++i) {
            this->m_geom[i] = geom[i];
            this->m_ba[i]   = ba[i];
            this->m_dmap[i] = dmap[i];
            
            std::cout << "BA: " << this->m_ba[i] << std::endl;
            
        }
        
        std::cout << this->m_geom.size() << " " << this->m_ba.size() << " " << this->m_dmap.size() << std::endl;
        
        std::cout << "maxLevel_m = " << maxLevel_m << std::endl;
        std::cout << "this->m_nlevels = " << this->m_nlevels << std::endl;
        std::cout << "BoxLibLayout::define" << std::endl;
    }
                   
    
    inline void setFinestLevel(int finestLevel);
    
    inline void setMaxLevel(int maxLevel);
    
    // overwritten functions
    inline bool LevelDefined (int level) const;
    
    inline int finestLevel () const;
    
    inline int maxLevel () const;
    
    inline AmrIntVect_t refRatio (int level) const;
    
    inline int MaxRefRatio (int level) const;
    
    
    /*!
     * Linear mapping to AMReX computation domain [a,b]^3. Each dimension
     * is mapped independently. The potential and electric field need to be scaled
     * afterwards appropriately.
     * 
     * [a, b] --> [c, d]
     * 
     * y = (c - d) / (a - b) * x + (a * d - b * c) / (a - b)
     * 
     * @param PData is the particle data
     * @param lold old lower bound
     * @param uold old upper bound
     * @param lnew new lower bound
     * @param unew new upper bound
     * @returns scaling factor
     */
    Vector_t domainMapping(AmrParticleBase< BoxLibLayout<T,Dim> >& PData,
                           const Vector_t& lold,
                           const Vector_t& uold,
                           const Vector_t& lnew,
                           const Vector_t& unew);
    
private:
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    bool Where (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                const unsigned int ip, 
                int lev_min = 0, int lev_max = -1, int nGrow = 0) const;

    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    bool EnforcePeriodicWhere (AmrParticleBase< BoxLibLayout<T,Dim> >& prt,
                               const unsigned int ip,
                               int lev_min = 0, int lev_max = -1) const;

    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Returns true if the particle was shifted.
    bool PeriodicShift (SingleParticlePos_t R) const;
    
    
//     int getTileIndex(const IntVect& iv, const Box& box, Box& tbx);
    
    void locateParticle(AmrParticleBase< BoxLibLayout<T,Dim> >& p, 
                        const unsigned int ip,
                        int lev_min, int lev_max, int nGrow,
                        bool &particleLeftDomain) const
    {
        bool outside = D_TERM( p.R[ip](0) <  AmrGeometry_t::ProbLo(0)
                            || p.R[ip](0) >= AmrGeometry_t::ProbHi(0),
                            || p.R[ip](1) <  AmrGeometry_t::ProbLo(1)
                            || p.R[ip](1) >= AmrGeometry_t::ProbHi(1),
                            || p.R[ip](2) <  AmrGeometry_t::ProbLo(2)
                            || p.R[ip](2) >= AmrGeometry_t::ProbHi(2));
        
        if ( outside ) {
            std::cout << "outside: " << outside << std::endl; std::cin.get();
        }
        
        bool success;
        if (outside)
        {
            // Note that EnforcePeriodicWhere may shift the particle if it is successful.
            success = EnforcePeriodicWhere(p, ip, lev_min, lev_max);
            if (!success && lev_min == 0)
            {
                // The particle has left the domain; invalidate it.
                particleLeftDomain = true;
                p.destroy(1, ip);
                success = true;
            }
        }
        else
        {
            success = Where(p, ip, lev_min, lev_max);
        }
    
        if (!success)
        {
            success = (nGrow > 0) && Where(p, ip, lev_min, lev_min, nGrow);
        }
    
        if (!success)
        {
            amrex::Abort("ParticleContainer<NR, NI, NA>::locateParticle(): invalid particle.");
        }
    }
    
private:
    int finestLevel_m;
    int maxLevel_m;
    // don't use m_rr from ParGDB since it is the same refinement in all directions
    AmrIntVectContainer_t refRatio_m;    // Refinement ratios [0:finest_level-1]
};

#include "BoxLibLayout.hpp"

#endif
