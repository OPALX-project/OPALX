#ifndef BOXLIB_LAYOUT_H
#define BOXLIB_LAYOUT_H

#include "AmrParticle/ParticleAmrLayout.h"
#include "AmrParticle/AmrParticleBase.h"

#include "Amr/AmrDefs.h"

#include <ParGDB.H>

template<class T, unsigned Dim>
class BoxLibLayout : public ParticleAmrLayout<T, Dim>,
                     public ParGDB
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
    
    static bool do_tiling;
    static IntVect tile_size;

public:
    
    /*!
     * Initializes default Geometry, DistributionMapping and BoxArray.
     */
    BoxLibLayout();
    
    BoxLibLayout(int nGridPoints, int maxGridSize, double lower, double upper);
    
    BoxLibLayout(const Geometry &geom,
                 const DistributionMapping &dmap,
                 const BoxArray &ba);
    
    BoxLibLayout(const Array<Geometry> &geom,
                 const Array<DistributionMapping> &dmap,
                 const Array<BoxArray> &ba,
                 const Array<int> &rr);
    
    void update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData, const ParticleAttrib<char>* canSwap = 0);
    
    void update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData, int lev_min = 0,
                const ParticleAttrib<char>* canSwap = 0);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    IntVect Index (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                   const unsigned int ip, int leve) const;

    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    IntVect Index (SingleParticlePos_t &R, int lev) const;
    
    
    void initDefaultBox(int nGridPoints, int maxGridSize,
                        double lower, double upper);
    
    void resize(int maxLevel) {
        int length = maxLevel + 1;
        this->m_geom.resize(length);
        this->m_dmap.resize(length);
        this->m_ba.resize(length);
        this->m_nlevels = length;
        maxLevel_m = maxLevel;
        
        for (int i = 0; i < length; ++i) {
            std::cout << "resize BA: " << this->m_ba[i] << std::endl;
            
        }
    }
    
    void define(const Array<Geometry>& geom,
                const Array<BoxArray>& ba,
                const Array<DistributionMapping>& dmap,
                const Array<int> & rr)
    {
        maxLevel_m = ba.size() - 1;
        this->m_nlevels = ba.size();
        this->m_geom.resize( this->m_nlevels );
        this->m_ba.resize( this->m_nlevels );
        this->m_dmap.resize( this->m_nlevels );
        this->m_rr.resize(maxLevel_m);
        
        for (int i = 0; i < this->m_nlevels; ++i) {
            this->m_geom[i] = geom[i];
            this->m_ba[i]   = ba[i];
            this->m_dmap[i] = dmap[i];
            this->m_rr[i]   = rr[i];
        }
    }
    
    void define(const Array<Geometry>& geom) {
//         this->m_geom.resize( geom.size() );
        std::cout << "this->m_geom.size() = " << this->m_geom.size() << std::endl;
        std::cout << "geom.size() = " << geom.size() << std::endl;
        for (unsigned int i = 0; i < geom.size(); ++i)
            this->m_geom[i] = geom[i];
    }
    
    
    void define(const Array<Geometry>& geom,
                const Array<BoxArray>& ba,
                const Array<DistributionMapping>& dmap)
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
    
//     inline IntVect refRatio (int level) const;
    
//     inline int MaxRefRatio (int level) const;
    
private:
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Checks/sets a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    bool Where (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                const unsigned int ip, 
                int lev_min = 0, int lev_max = -1, int nGrow=0) const;

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
        bool outside = D_TERM( p.R[ip](0) <  Geometry::ProbLo(0)
                            || p.R[ip](0) >= Geometry::ProbHi(0),
                            || p.R[ip](1) <  Geometry::ProbLo(1)
                            || p.R[ip](1) >= Geometry::ProbHi(1),
                            || p.R[ip](2) <  Geometry::ProbLo(2)
                            || p.R[ip](2) >= Geometry::ProbHi(2));
        
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
            BoxLib::Abort("ParticleContainer<NR, NI, NA>::locateParticle(): invalid particle.");
        }
    }
    
private:
    int finestLevel_m;
    int maxLevel_m;
};

#include "BoxLibLayout.hpp"

#endif
