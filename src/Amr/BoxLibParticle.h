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

    
template<class PLayout>
class BoxLibParticle : public virtual AmrParticleBase<PLayout>
{
public:
    typedef typename AmrParticleBase<PLayout>::ParticlePos_t        ParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::ParticleIndex_t      ParticleIndex_t;
    typedef typename AmrParticleBase<PLayout>::SingleParticlePos_t  SingleParticlePos_t;
    typedef typename AmrParticleBase<PLayout>::AmrField_t           AmrField_t;
    typedef typename AmrParticleBase<PLayout>::AmrFieldContainer_t  AmrFieldContainer_t;
    
    typedef double RealType;
    typedef std::deque<Particle<1,0> > C;
    typedef C PBox;
    typedef Particle<1,0> ParticleType;
    typedef typename std::map<int,PBox> PMap;
    
public:
    BoxLibParticle();
    
    BoxLibParticle(PLayout *layout);
    
    // scatter the data from the given attribute onto the given Field, using
    // the given Position attribute
    
    template <class FT, unsigned Dim, class PT>
    void scatter(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                 ParticleAttrib<Vektor<PT, Dim> >& pp,
                 int lbase = 0, int lfine = -1);
    
    // gather the data from the given Field into the given attribute, using
    // the given Position attribute
    template <class FT, unsigned Dim, class PT>
    void gather(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                ParticleAttrib<Vektor<PT, Dim> >& pp,
                int lbase = 0, int lfine = -1);
    
private:
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    static void CIC_Cells_Fracs_Basic (const SingleParticlePos_t &R,
                                       const Real* plo, const Real* dx,
                                       Real* fracs,  IntVect* cells);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    static int CIC_Cells_Fracs (const SingleParticlePos_t &R,
                                const Real*         plo,
                                const Real*         dx_geom,
                                const Real*         dx_part,
                                Array<Real>&        fracs,
                                Array<IntVect>&     cells);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    bool FineToCrse (const int ip,
                     int                                flev,
                     PLayout*                           layout_p,
                     const Array<IntVect>&              fcells,
                     const BoxArray&                    fvalid,
                     const BoxArray&                    compfvalid_grown,
                     Array<IntVect>&                    ccells,
                     Array<Real>&                       cfracs,
                     Array<int>&                        which,
                     Array<int>&                        cgrid,
                     Array<IntVect>&                    pshifts,
                     std::vector< std::pair<int,Box> >& isects);
    
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    void FineCellsToUpdateFromCrse (const int ip,
                                    int lev,
                                    PLayout* layout_p,
                                    const IntVect& ccell,
                                    const IntVect& cshift,
                                    Array<int>& fgrid,
                                    Array<Real>& ffrac,
                                    Array<IntVect>& fcells,
                                    std::vector< std::pair<int,Box> >& isects);
    
    
    //Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    //sends/receivs the particles that are needed by other processes to during AssignDensity
    void AssignDensityDoit(int level, AmrFieldContainer_t* mf, PMap& data,
                           int ncomp, int lev_min = 0);
    
    // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
    // Assign values from grid back to particles
    void Interp(const SingleParticlePos_t &R, const Geometry &geom, const FArrayBox& fab, 
                const int* idx, Real* val, int cnt);
    
    template<class AType>
    void AssignDensity(ParticleAttrib<AType> &pa,
                       bool sub_cycle,
                       AmrFieldContainer_t& mf_to_be_filled,
                       int lev_min, int finest_level);
    
    template<class AType>
    void AssignDensitySingleLevel (ParticleAttrib<AType> &pa,
                                   AmrField_t& mf_to_be_filled,
                                   int lev,
                                   int particle_lvl_offset = 0);
    
    template<class AType>
    void AssignCellDensitySingleLevel(ParticleAttrib<AType> &pa,
                                      AmrField_t& mf_to_be_filled,
                                      int lev, int particle_lvl_offset = 0);
    
    template<class AType>
    void NodalDepositionSingleLevel(ParticleAttrib<AType> &pa,
                                    AmrField_t& mf_to_be_filled,
                                    int lev, int particle_lvl_offset = 0);
    
    template<class AType>
    void GetGravity(ParticleAttrib<AType> &pa, AmrFieldContainer_t &mf);
    
private:
    bool allow_particles_near_boundary_m;
    
    IpplTimings::TimerRef AssignDensityTimer_m;
};


#include "BoxLibParticle.hpp"

#endif
