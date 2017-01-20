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
#include <Particles_F.H>
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

  ParticleIndex_t m_lev;
  ParticleIndex_t m_grid;
  
public: 

  AmrParticleBase()  { }

  // destructor - delete the layout if necessary
  ~AmrParticleBase() { }

  void initializeAmr() {
    this->addAttribute(m_lev);
    this->addAttribute(m_grid);
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


  
  //
  // Checks/sets a particle's location on a specific level.
  // (Yes this is distict from the functionality provided above)
  //
  static bool SingleLevelWhere (AmrParticleBase& p, const ParGDBBase* gdb, int level);
  //
  // Updates a particle's location (Where), tries to periodic shift any particles
  // that have left the domain. May need work (see inline comments)
  //
  static void Reset (AmrParticleBase& prt, const ParGDBBase* gdb, bool update, bool verbose=true); 
  
  static Real InterpDoit (const FArrayBox& fab, const Real* fracs, const IntVect* cells, int comp);

  static Real InterpDoit (const FArrayBox& fab, const IntVect& hi, const Real* frac, int comp);
  
  static void Interp (const AmrParticleBase& prt, const Geometry& geom, 
		      const FArrayBox& fab, const int* idx, Real* val, int cnt);

  static const std::string& Version ();

  static const std::string& DataPrefix ();

  static void GetGravity (const FArrayBox& gfab, const Geometry& geom, 
			  const AmrParticleBase& p, Real* grav);

  static int MaxReaders ();

  static long MaxParticlesPerRead ();
  //
  // Returns the next particle ID for this processor.
  // Particle IDs start at 1 and are never reused.
  // The pair, consisting of the ID and the CPU on which the particle is "born",
  // is a globally unique identifier for a particle.  The maximum of this value
  // across all processors must be checkpointed and then restored on restart
  // so that we don't reuse particle IDs.
  //
  static int NextID ();
  // This version can only be used inside omp critical.
  static int UnprotectedNextID ();
  //
  // Reset on restart.
  //
  static void NextID (int nextid);
  //
  // Used by AssignDensity.
  //
  static bool CrseToFine (const BoxArray&       cfba, 
			  const Array<IntVect>& cells, 
			  Array<IntVect>&       cfshifts, 
			  const Geometry&       gm, 
			  Array<int>&           which, 
			  Array<IntVect>&       pshifts);

  static bool FineToCrse (const AmrParticleBase&                p, 
			  int                                flev, 
			  const ParGDBBase*                  gdb, 
			  const Array<IntVect>&              fcells, 
			  const BoxArray&                    fvalid, 
			  const BoxArray&                    compfvalid_grown, 
			  Array<IntVect>&                    ccells, 
			  Array<Real>&                       cfracs, 
			  Array<int>&                        which, 
			  Array<int>&                        cgrid, 
			  Array<IntVect>&                    pshifts, 
			  std::vector< std::pair<int,Box> >& isects);

  static void FineCellsToUpdateFromCrse (const AmrParticleBase&                p, 
					 int lev, const ParGDBBase*         gdb, 
					 const IntVect&                     ccell,
					 const IntVect&                     cshift, 
					 Array<int>&                        fgrid, 
					 Array<Real>&                       ffrac, 
					 Array<IntVect>&                    fcells, 
					 std::vector< std::pair<int,Box> >& isects);

  static void CIC_Fracs (const Real* frac, Real* fracs);
  
  static void CIC_Cells (const IntVect& hicell, IntVect* cells);
  //
  // Old, *-based CIC for use in Interp.
  //
  static void CIC_Cells_Fracs_Basic (const AmrParticleBase& p, const Real* plo, 
				     const Real* dx, Real* fracs,  IntVect* cells);
  //
  // Wraps the arbitrary dx function.
  //
    static int CIC_Cells_Fracs (const AmrParticleBase& p, 
                                const Real*         plo, 
                                const Real*         dx, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
    //
    // Does CIC computations for arbitrary particle/grid dx's.
    //
    static int CIC_Cells_Fracs (const AmrParticleBase& p, 
                                const Real*         plo, 
                                const Real*         dx_geom, 
                                const Real*         dx_part, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
    //
    // Useful for sorting particles into lexicographic order of their cell position.
    //
    class Compare
    {
    public:
        bool operator () (const AmrParticleBase& lhs,
                          const AmrParticleBase& rhs) const
        {
            return lhs.m_cell.lexLT(rhs.m_cell);
        }
    };
};


#endif
