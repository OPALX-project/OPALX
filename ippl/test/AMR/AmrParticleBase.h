#ifndef AMRPARTBASE_H
#define AMRPARTBASE_H

#include "Ippl.h"

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



//AMRPArticleBase class definition. Template parameter is the specific AmrParticleLayout-derived
//class which determines how the particles are distribute amoung processors.
//AmrParticleBase is derived from IpplParticle base and 

typedef double RealType;

typedef long                    SortListIndex_t;
typedef std::vector<SortListIndex_t> SortList_t;
typedef std::vector<ParticleAttribBase *>      attrib_container_t;

template<class PLayout>
class AmrParticleBase : public IpplParticleBase<PLayout> {

public:
  typedef typename PLayout::ParticlePos_t   ParticlePos_t;
  typedef typename PLayout::ParticleIndex_t ParticleIndex_t;
  typedef typename PLayout::SingleParticlePos_t SingleParticlePos_t;

  ParticleIndex_t m_lev;
  ParticleIndex_t m_grid;
  ParticleAttrib<double> qm;

private:


  bool allow_particles_near_boundary;

  static void CIC_Cells_Fracs_Basic (const SingleParticlePos_t &R, const Real* plo, 
				     const Real* dx, Real* fracs,  IntVect* cells);

  static int CIC_Cells_Fracs (const SingleParticlePos_t &R,
			      const Real*         plo,
			      const Real*         dx_geom,
			      const Real*         dx_part,
			      Array<Real>&        fracs,
			      Array<IntVect>&     cells);
  
public: 

  AmrParticleBase() : allow_particles_near_boundary(false) { }

  // destructor - delete the layout if necessary
  ~AmrParticleBase() { }

  void initializeAmr() {
    this->addAttribute(m_lev);
    this->addAttribute(m_grid);
    this->addAttribute(qm);
  }

  void setAllowParticlesNearBoundary(bool value) {
    allow_particles_near_boundary = value;
  }

  // Update the particle object after a timestep.  This routine will change
  // our local, total, create particle counts properly.
  void update() {
    
    std::cout << "update Amr" << std::endl;
    // make sure we've been initialized
    PLayout *Layout = &this->getLayout();

    PAssert(Layout != 0);
    
    // ask the layout manager to update our atoms, etc.
    Layout->update(*this);
    INCIPPLSTAT(incParticleUpdates);

  }

  void update(const ParticleAttrib<char>& canSwap) {
    // make sure we've been initialized
    
    std::cout << "update Amr swap" << std::endl;
    PLayout *Layout = &this->getLayout();
    PAssert(Layout != 0);
    
    // ask the layout manager to update our atoms, etc.
    Layout->update(*this, &canSwap);
    INCIPPLSTAT(incParticleUpdates);
    
  }

  //sort particles based on the grid and level that they belong to
  void sort() {
    size_t LocalNum = this->getLocalNum();

    SortList_t slist;
    m_grid.calcSortList(slist);
    this->sort(slist);
    m_lev.calcSortList(slist);
    this->sort(slist);
  }

  void sort(SortList_t &sortlist) {
    attrib_container_t::iterator abeg = this->begin();
    attrib_container_t::iterator aend = this->end();
    for ( ; abeg != aend; ++abeg )
      (*abeg)->sort(sortlist);
  }

  template <class AType> 
  void AssignDensitySingleLevel (ParticleAttrib<AType> &pa, 
				 MultiFab& mf_to_be_filled,
				 int lev,
				 int particle_lvl_offset = 0)
  {
    std::cout << "Assign density single level" << std::endl;

    if (mf_to_be_filled.is_nodal()) {
      std::cout << "Nodal" << std::endl;
      NodalDepositionSingleLevel(pa, mf_to_be_filled, lev, particle_lvl_offset);
    } else if (mf_to_be_filled.boxArray().ixType().cellCentered()) {
      std::cout << "Cell" << std::endl;
      AssignCellDensitySingleLevel(pa, mf_to_be_filled, lev, particle_lvl_offset);
    } else {
      BoxLib::Abort("AssignCellDensitySingleLevel: mixed type not supported");
    }
  }

  template <class AType> 
  void AssignDensitySingleLevel (ParticleAttrib<AType> &pa, 
				 int lev,
				 int particle_lvl_offset = 0)
  {
    std::cout << "Assign density single level" << std::endl;

  }

  template <class AType>
  void AssignCellDensitySingleLevel(ParticleAttrib<AType> &pa,
				    MultiFab& mf_to_be_filled,
				    int lev,
				    int particle_lvl_offset = 0)
  {
    
    MultiFab* mf_pointer;
    PLayout *Layout = &this->getLayout();
    const ParGDBBase* m_gdb = Layout->GetParGDB();

    if ( m_gdb->OnSameGrids(lev, mf_to_be_filled) ) {
      // If we are already working with the internal mf defined on the 
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    } else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need 
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), 0, mf_to_be_filled.nGrow(),
				m_gdb->ParticleDistributionMap(lev), Fab_allocate);
    }
  
    // We must have ghost cells for each FAB so that a particle in one grid can spread its effect to an
    //    adjacent grid by first putting the value into ghost cells of its own grid.  The mf->sumBoundary call then
    //    adds the value from one grid's ghost cell to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
      BoxLib::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
      BoxLib::Error("AssignDensity: problem must be periodic in no or all directions");
    }

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
      (*mf_pointer)[mfi].setVal(0);
    }
    
    //loop trough particles and distribute values on the grid
    size_t LocalNum = this->getLocalNum();
    for (size_t ip = 0; ip < LocalNum; ++ip) 
    {
      
      //if particle doesn't belong on this level exit loop
      if (m_lev[ip] != (unsigned)lev)
	break;

      FArrayBox& fab = (*mf_pointer)[m_grid[ip]];

      Array<Real> fracs;
      Array<IntVect> cells;

      const int M = CIC_Cells_Fracs(this->R[ip], plo, dx, dx_particle, fracs, cells);
      //
      // If this is not fully periodic then we have to be careful that the
      // particle's support leaves the domain unless we specifically want to ignore
      // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
      // We test this by checking the low and high corners respectively.
      //
      if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
	if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
	  BoxLib::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
	}
      }

      for (int i = 0; i < M; ++i)
      {
	if ( !fab.box().contains(cells[i])) {
	  continue;
	}

	// If the domain is not periodic and we want to let particles
	//    live near the boundary but "throw away" the contribution that 
	//    does not fall into the domain ...
	if ( ! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i])) {
	  continue;
	}

	//sum up mass in the particle attribute
	fab(cells[i], 0) += pa[ip] * fracs[i];

      }
      
    }

    mf_pointer->SumBoundary(gm.periodicity());

    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
      mf_to_be_filled.copy(*mf_pointer,0,0,0);
      delete mf_pointer;
    }
    
  }

  template<class AType>
  void NodalDepositionSingleLevel(ParticleAttrib<AType> &pa,
				  MultiFab& mf_to_be_filled,
				  int lev,
				  int particle_lvl_offset = 0)
  {
    MultiFab* mf_pointer;
    PLayout *Layout = &this->getLayout();
    const ParGDBBase* m_gdb = Layout->GetParGDB();

    if ( m_gdb->OnSameGrids(lev, mf_to_be_filled) )
    {
      // If we are already working with the internal mf defined on the 
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else
    {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need 
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(BoxLib::convert(m_gdb->ParticleBoxArray(lev),
						mf_to_be_filled.boxArray().ixType()),
				0, mf_to_be_filled.nGrow(),
				m_gdb->ParticleDistributionMap(lev), Fab_allocate);
    }

    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     dx          = gm.CellSize();
    
    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) 
      BoxLib::Error("AssignDensity: problem must be periodic in no or all directions");

    mf_pointer->setVal(0.0);

    Array<IntVect> cells;
    cells.resize(8);
    
    Array<Real> fracs;
    fracs.resize(8);

    Array<Real> sx;
    sx.resize(2);
    Array<Real> sy;
    sy.resize(2);
    Array<Real> sz;
    sz.resize(2);

    //loop trough particles and distribute values on the grid
    size_t LocalNum = this->getLocalNum();
    for (size_t ip = 0; ip < LocalNum; ++ip) 
    {
      FArrayBox& fab = (*mf_pointer)[m_grid[ip]];
      
      IntVect m_cell = Layout->Index(this->R[ip], gm);
      cells[0] = m_cell;
      cells[1] = m_cell+IntVect(1,0,0);
      cells[2] = m_cell+IntVect(0,1,0);
      cells[3] = m_cell+IntVect(1,1,0);
      cells[4] = m_cell+IntVect(0,0,1);
      cells[5] = m_cell+IntVect(1,0,1);
      cells[6] = m_cell+IntVect(0,1,1);
      cells[7] = m_cell+IntVect(1,1,1);

      Real x = this->R[ip][0] / dx[0];
      Real y = this->R[ip][1] / dx[1];
      Real z = this->R[ip][2] / dx[2];

      int i = m_cell[0];
      int j = m_cell[1];
      int k = m_cell[2];

      Real xint = x - i;
      Real yint = y - j;
      Real zint = z - k;

      sx[0] = 1.0-xint;
      sx[1] = xint;
      sy[0] = 1.0-yint;
      sy[1] = yint;
      sz[0] = 1.0-zint;
      sz[1] = zint;

      fracs[0] = sx[0] * sy[0] * sz[0];
      fracs[1] = sx[1] * sy[0] * sz[0];
      fracs[2] = sx[0] * sy[1] * sz[0];
      fracs[3] = sx[1] * sy[1] * sz[0];
      fracs[4] = sx[0] * sy[0] * sz[1];
      fracs[5] = sx[1] * sy[0] * sz[1];
      fracs[6] = sx[0] * sy[1] * sz[1];
      fracs[7] = sx[1] * sy[1] * sz[1];

      for (int i = 0; i < 8; i++)
	fab(cells[i],0) += pa[ip] * fracs[i];

    }

    mf_pointer->SumBoundary(gm.periodicity());
    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,0);
	delete mf_pointer;
    }

  }

  
};

#include "AmrParticleBase.hpp"

#endif
