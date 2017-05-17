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

#include <ParmParse.H>

#include <ParGDB.H>
#include <REAL.H>
#include <IntVect.H>
#include <Array.H>
#include <Utility.H>
#include <Geometry.H>
#include <VisMF.H>
#include <Particles.H>
#include <RealBox.H>

#include <BLFort.H>
#include <MultiFabUtil.H>
#include <MultiFabUtil_F.H>
#include <Interpolater.H>
#include <FillPatchUtil.H>

    
template<class PLayout>
class BoxLibParticle : public virtual AmrParticleBase<PLayout>
{
public:
    typedef typename AmrParticleBase<PLayout>::ParticlePos_t        ParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::ParticleIndex_t      ParticleIndex_t;
    typedef typename AmrParticleBase<PLayout>::SingleParticlePos_t  SingleParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::AmrField_t           AmrField_t;
    typedef typename AmrParticleBase<PLayout>::AmrFieldContainer_t  AmrFieldContainer_t; // Array<std::unique_ptr<MultiFab> >
    
public:
    BoxLibParticle();
    
    BoxLibParticle(PLayout *layout);
    
    // scatter the data from the given attribute onto the given Field, using
    // the given Position attribute
    
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 int lbase = 0, int lfine = -1);
    
    
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrField_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 int level = 0);
    
    // gather the data from the given Field into the given attribute, using
    // the given Position attribute
    template <class FT, unsigned Dim, class PT>
    void gather(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                ParticleAttrib<Vektor<PT, Dim> >& pp,
                int lbase = 0, int lfine = -1);
    
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
    
    // amrex repository AMReX_MultiFabUtil.H (missing in BoxLib repository)
    void sum_fine_to_coarse(/*const */MultiFab& S_fine, MultiFab& S_crse,
                            int scomp, int ncomp, const IntVect& ratio,
                            const Geometry& cgeom, const Geometry& fgeom) const
    {
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());
        BL_ASSERT(ratio == ratio[0]);
        BL_ASSERT(S_fine.nGrow() % ratio[0] == 0);

        const int nGrow = S_fine.nGrow() / ratio[0];

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA, ncomp, nGrow, S_fine.DistributionMap());

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_S_fine, true); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& tbx = mfi.growntilebox(nGrow);
            
            BL_FORT_PROC_CALL(BL_AVGDOWN, bl_avgdown)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi] , scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi], 0),
                 ratio.getVect(),&ncomp);
        }
        
        S_crse.copy(crse_S_fine, 0, scomp, ncomp, nGrow, 0,
                    cgeom.periodicity(), FabArrayBase::ADD);
    }
    
private:
    bool allow_particles_near_boundary_m;
    
    IpplTimings::TimerRef AssignDensityTimer_m;
};


#include "BoxLibParticle.hpp"

#endif
