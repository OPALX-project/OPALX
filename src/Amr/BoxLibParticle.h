#ifndef BOXLIB_PARTICLE_H
#define BOXLIB_PARTICLE_H

#include "AmrParticle/AmrParticleBase.h"

#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>

#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Particles.H>
#include <AMReX_RealBox.H>

#include <AMReX_BLFort.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_F.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>
    
template<class PLayout>
class BoxLibParticle : public virtual AmrParticleBase<PLayout>
{
public:
    typedef typename AmrParticleBase<PLayout>::ParticlePos_t        ParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::ParticleIndex_t      ParticleIndex_t;
    typedef typename AmrParticleBase<PLayout>::SingleParticlePos_t  SingleParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::AmrField_t           AmrField_t;
    typedef typename AmrParticleBase<PLayout>::AmrFieldContainer_t  AmrFieldContainer_t; // Array<std::unique_ptr<MultiFab> >
    
    typedef typename PLayout::AmrProcMap_t  AmrProcMap_t;
    typedef typename PLayout::AmrGrid_t     AmrGrid_t;
    typedef typename PLayout::AmrGeometry_t AmrGeometry_t;
    typedef typename PLayout::AmrIntVect_t  AmrIntVect_t;
    typedef typename PLayout::AmrBox_t      AmrBox_t;
    typedef typename PLayout::AmrReal_t     AmrReal_t;
    
    typedef amrex::FArrayBox                FArrayBox_t;
    
public:
    BoxLibParticle();
    
    BoxLibParticle(PLayout *layout);
    
    // scatter the data from the given attribute onto the given Field, using
    // the given Position attribute
    
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 int lbase, int lfine);
    
    
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrField_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 int level = 0);
    
    // gather the data from the given Field into the given attribute, using
    // the given Position attribute
    template <class FT, unsigned Dim, class PT>
    void gather(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                ParticleAttrib<Vektor<PT, Dim> >& pp,
                int lbase, int lfine);
    
//     // BoxLib specific functions
//     void resizeContainerGDB(int length) {
//         PLayout *layout_p = &this->getLayout();
//         layout_p->resizeGDB(length);
//     }
    
//     void define(const Array<Geometry>& geom,
//                 const Array<BoxArray>& ba,
//                 const Array<DistributionMapping>& dmap,
//                 const Array<int> & rr)
//     {
//         PLayout *layout_p = &this->getLayout();
//         layout_p->setGDB(geom, ba, dmap, rr);
//     }
    
//     void updateGDB(const Array<Geometry>& geom,
//                    const Array<BoxArray>& ba,
//                    const Array<DistributionMapping>& dmap) {
//         PLayout *layout_p = &this->getLayout();
//         layout_p->updateGDB(geom, ba, dmap);
//     }
    
private:
    
    
    template <class AType>
    void AssignDensityFort(ParticleAttrib<AType> &pa,
                           AmrFieldContainer_t& mf_to_be_filled, 
                           int lev_min, int ncomp, int finest_level) const;
    
    template <class AType>
    void InterpolateFort(ParticleAttrib<AType> &pa,
                         AmrFieldContainer_t& mesh_data, 
                         int lev_min, int lev_max);
    
    template <class AType>
    void InterpolateSingleLevelFort(ParticleAttrib<AType> &pa, AmrField_t& mesh_data, int lev);
    
    template <class AType>
    void AssignCellDensitySingleLevelFort(ParticleAttrib<AType> &pa, AmrField_t& mf, int level,
                                          int ncomp=1, int particle_lvl_offset = 0) const;
    
private:
    bool allow_particles_near_boundary_m;
    
    IpplTimings::TimerRef AssignDensityTimer_m;
};


#include "BoxLibParticle.hpp"

#endif
