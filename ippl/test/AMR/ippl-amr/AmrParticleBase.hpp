#include "Ippl.h"

#include <Particles.H>
#include <ParmParse.H>
#include <limits>

// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
template<class PLayout>
void AmrParticleBase<PLayout>::Interp(const SingleParticlePos_t &R,
				      const Geometry &geom,
				      const FArrayBox& fab,
				      const int* idx,
				      Real* val,
				      int cnt)
{
  const int       M   = D_TERM(2,+2,+4);
  const Real*     plo = geom.ProbLo();
  const Real*     dx  = geom.CellSize();
  
  Real    fracs[M];
  IntVect cells[M];
  //
  // Get "fracs" and "cells".
  //
  CIC_Cells_Fracs_Basic(R, plo, dx, fracs, cells);
  
  for (int i = 0; i < cnt; i++)
  {
    BL_ASSERT(idx[i] >= 0 && idx[i] < fab.nComp());
    val[i] = ParticleBase::InterpDoit(fab,fracs,cells,idx[i]);
  }
}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
template<class PLayout>
void AmrParticleBase<PLayout>::CIC_Cells_Fracs_Basic(const SingleParticlePos_t &R, 
						     const Real* plo, 
						     const Real* dx, 
						     Real* fracs,  
						     IntVect* cells)
{

  const Real len[BL_SPACEDIM] = { D_DECL((R[0]-plo[0])/dx[0] + Real(0.5),
					 (R[1]-plo[1])/dx[1] + Real(0.5),
					 (R[2]-plo[2])/dx[2] + Real(0.5)) };

  const IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));
  
  const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

  ParticleBase::CIC_Fracs(frac, fracs);
  ParticleBase::CIC_Cells(cell, cells);

}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
template<class PLayout>
int AmrParticleBase<PLayout>::CIC_Cells_Fracs (const SingleParticlePos_t &R,
					       const Real*         plo,
					       const Real*         dx_geom,
					       const Real*         dx_part,
					       Array<Real>&        fracs,
					       Array<IntVect>&     cells)
{
    BL_PROFILE("AmrParticleBase::CIC_Cells_Fracs()");
    if (dx_geom == dx_part)
    {
        const int M = D_TERM(2,+2,+4);
        fracs.resize(M);
        cells.resize(M);
        AmrParticleBase::CIC_Cells_Fracs_Basic(R,plo,dx_geom,fracs.dataPtr(),cells.dataPtr());
        return M;
    }
    //
    // The first element in fracs and cells is the lowest corner, the last is the highest.
    //
    const Real hilen[BL_SPACEDIM] = { D_DECL((R[0]-plo[0]+dx_part[0]/2)/dx_geom[0],
                                             (R[1]-plo[1]+dx_part[1]/2)/dx_geom[1],
                                             (R[2]-plo[2]+dx_part[2]/2)/dx_geom[2]) };

    const Real lolen[BL_SPACEDIM] = { D_DECL((R[0]-plo[0]-dx_part[0]/2)/dx_geom[0],
                                             (R[1]-plo[1]-dx_part[1]/2)/dx_geom[1],
                                             (R[2]-plo[2]-dx_part[2]/2)/dx_geom[2]) };

    const IntVect hicell(D_DECL(floor(hilen[0]), floor(hilen[1]), floor(hilen[2])));
    
    const IntVect locell(D_DECL(floor(lolen[0]), floor(lolen[1]), floor(lolen[2])));
    
    const Real cell_density = D_TERM(dx_geom[0]/dx_part[0],*dx_geom[1]/dx_part[1],*dx_geom[2]/dx_part[2]);
    
    const int M = D_TERM((hicell[0]-locell[0]+1),*(hicell[1]-locell[1]+1),*(hicell[2]-locell[2]+1));

    fracs.resize(M);
    cells.resize(M);
    //
    // This portion might be slightly inefficient. Feel free to redo it if need be.
    //
    int i = 0;
    for (int zi = locell[2]; zi <= hicell[2]; zi++)
    {
        const Real zf = std::min(hilen[2]-zi,Real(1))-std::max(lolen[2]-zi,Real(0));
        for (int yi = locell[1]; yi <= hicell[1]; yi++)
        {
            const Real yf = std::min(hilen[1]-yi,Real(1))-std::max(lolen[1]-yi,Real(0));
            for (int xi = locell[0]; xi <= hicell[0]; xi++)
            {
                cells[i][0] = xi;
                cells[i][1] = yi;
                cells[i][2] = zi;
                fracs[i] = zf * yf * (std::min(hilen[0]-xi,Real(1))
				      -std::max(lolen[0]-xi,Real(0))) * cell_density;
                i++;
            }
        }
    }

    return M;
}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
template<class PLayout>
bool AmrParticleBase<PLayout>::FineToCrse (const int ip,
					   int                               flev,
					   const ParGDBBase*                  gdb,
					   const Array<IntVect>&              fcells,
					   const BoxArray&                    fvalid,
					   const BoxArray&                    compfvalid_grown,
					   Array<IntVect>&                    ccells,
					   Array<Real>&                       cfracs,
					   Array<int>&                        which,
					   Array<int>&                        cgrid,
					   Array<IntVect>&                    pshifts,
					   std::vector< std::pair<int,Box> >& isects)
{

  PLayout *Layout = &this->getLayout();

  BL_PROFILE("ParticleBase::FineToCrse()");
  BL_ASSERT(gdb != 0);
  BL_ASSERT(flev > 0);
  //
  // We're in AssignDensity(). We want to know whether or not updating
  // with a particle we'll cross a fine->crse boundary.  Note that crossing
  // a periodic boundary, where the periodic shift lies in our valid region,
  // is not considered a Fine->Crse crossing.
  //
  const int M = fcells.size();

  which.resize(M);
  cgrid.resize(M);
  ccells.resize(M);
  cfracs.resize(M);

  for (int i = 0; i < M; i++)
  {
    cgrid[i] = -1;
    which[i] =  0;
  }

  const Box& ibx = BoxLib::grow(gdb->ParticleBoxArray(flev)[m_grid[ip]],-1);
  IntVect m_cell = Layout->Index(this->R[ip], flev);

  BL_ASSERT(ibx.ok());

  if (ibx.contains(m_cell))
    //
    // We're strictly contained in our valid box.
    // We can't cross a fine->crse boundary.
    //
    return false;

  if (!compfvalid_grown.contains(m_cell))
    //
    // We're strictly contained in our "valid" region. Note that the valid
    // region contains any periodically shifted ghost cells that intersect
    // valid region.
    //
    return false;
  //
  // Otherwise ...
  //
  const Geometry& cgm = gdb->Geom(flev-1);
  const IntVect&  rr  = gdb->refRatio(flev-1);
  const BoxArray& cba = gdb->ParticleBoxArray(flev-1);
  
  CIC_Cells_Fracs(this->R[ip], cgm.ProbLo(), cgm.CellSize(), cgm.CellSize(), cfracs, ccells);

  bool result = false;

  for (int i = 0; i < M; i++)
  {
    IntVect ccell_refined = ccells[i]*rr;
    //
    // We've got to protect against the case when we're at the low
    // end of the domain because coarsening & refining don't work right
    // when indices go negative.
    //
    for (int dm = 0; dm < BL_SPACEDIM; dm++)
      ccell_refined[dm] = std::max(ccell_refined[dm], -1);

    if (!fvalid.contains(ccell_refined))
    {
      result   = true;
      which[i] = 1;

      Box cbx(ccells[i],ccells[i]);
    
      if (!cgm.Domain().contains(ccells[i]))
      {
	//
	// We must be at a periodic boundary.
	// Find valid box into which we can be periodically shifted.
	//
	BL_ASSERT(cgm.isAnyPeriodic());

	cgm.periodicShift(cbx, cgm.Domain(), pshifts);

	BL_ASSERT(pshifts.size() == 1);

	cbx -= pshifts[0];

	ccells[i] -= pshifts[0];
	BL_ASSERT(cbx.ok());
	BL_ASSERT(cgm.Domain().contains(cbx));
      }
      //
      // Which grid at the crse level do we need to update?
      //
      cba.intersections(cbx,isects,true,0);

      BL_ASSERT(!isects.empty());

      cgrid[i] = isects[0].first;  // The grid ID at crse level that we hit.
    }
  }

  return result;

}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
template<class PLayout>
void AmrParticleBase<PLayout>::FineCellsToUpdateFromCrse (
  const int ip,
  int lev,
  const ParGDBBase* gdb,
  const IntVect& ccell,
  const IntVect& cshift,
  Array<int>& fgrid,
  Array<Real>& ffrac,
  Array<IntVect>& fcells,
  std::vector< std::pair<int,Box> >& isects)
{
  BL_PROFILE("ParticleBase::FineCellsToUpdateFromCrse()");
  BL_ASSERT(lev >= 0);
  BL_ASSERT(lev < gdb->finestLevel());

  const Box&      fbx = BoxLib::refine(Box(ccell,ccell),gdb->refRatio(lev));
  const BoxArray& fba = gdb->ParticleBoxArray(lev+1);
  const Real*     plo = gdb->Geom(lev).ProbLo();
  const Real*     dx  = gdb->Geom(lev).CellSize();
  const Real*     fdx = gdb->Geom(lev+1).CellSize();

  if (cshift == IntVect::TheZeroVector())
  {
    BL_ASSERT(fba.contains(fbx));
  }
  //
  // Instead of clear()ing these we'll do a resize(0).
  // This'll preserve their capacity so that we'll only need
  // to do any memory allocation when their capacity needs to increase.
  //
  fgrid.resize(0);
  ffrac.resize(0);
  fcells.resize(0);
  //
  // Which fine cells does particle "p" (that wants to update "ccell") do we
  // touch at the finer level?
  //
  for (IntVect iv = fbx.smallEnd(); iv <= fbx.bigEnd(); fbx.next(iv))
  {
    bool touches = true;

    for (int k = 0; k < BL_SPACEDIM; k++)
    {
      const Real celllo = iv[k]  * fdx[k] + plo[k];
      const Real cellhi = celllo + fdx[k];

      if ((this->R[ip][k] < celllo) && (celllo > (this->R[ip][k] + dx[k]/2)))
	touches = false;

      if ((this->R[ip][k] > cellhi) && (cellhi < (this->R[ip][k] - dx[k]/2)))
	touches = false;
    }

    if (touches)
    {
      fcells.push_back(iv);
    }
  }

  Real sum_fine = 0;
  //
  // We need to figure out the fine fractions and the fine grid needed updating.
  //
  for (unsigned j = 0; j < fcells.size(); j++)
  {
    IntVect& iv = fcells[j];

    Real the_frac = 1;

    for (int k = 0; k < BL_SPACEDIM; k++)
    {
      const Real celllo = (iv[k] * fdx[k] + plo[k]);

      if (this->R[ip][k] <= celllo)
      {
	const Real isecthi = this->R[ip][k] + dx[k]/2;

	the_frac *= std::min((isecthi - celllo),fdx[k]);
      }
      else
      {
	const Real cellhi  = (iv[k]+1) * fdx[k] + plo[k];
	const Real isectlo = this->R[ip][k] - dx[k]/2;

	the_frac *= std::min((cellhi - isectlo),fdx[k]);
      }
    }

    ffrac.push_back(the_frac);

    sum_fine += the_frac;

    if (cshift != IntVect::TheZeroVector())
    {
      //
      // Update to the correct fine cell needing updating.
      // Note that "cshift" is from the coarse perspective.
      //
      const IntVect& fshift = cshift * gdb->refRatio(lev);
      //
      // Update fcells[j] to indicate a shifted fine cell needing updating.
      //
      iv -= fshift;
    }

    fba.intersections(Box(iv,iv),isects,true,0);

    BL_ASSERT(!isects.empty());

    fgrid.push_back(isects[0].first);
  }

  BL_ASSERT(ffrac.size() == fcells.size());
  BL_ASSERT(fgrid.size() == fcells.size());
  //
  // Now adjust the fine fractions so they sum to one.
  //
  for (unsigned j = 0; j < ffrac.size(); j++)
    ffrac[j] /= sum_fine;
}


template<class PLayout>
void AmrParticleBase<PLayout>::AssignDensityDoit(int rho_index,
						 PArray<MultiFab>* mf,
						 PMap&             data,
						 int               ncomp,
						 int               lev_min)
{
  if (rho_index != 0) BoxLib::Abort("AssignDensityDoit only works if rho_index = 0");

  BL_PROFILE("ParticleContainer<NR,NI,C>::AssignDensityDoit()");
  BL_ASSERT(NR >= ncomp);

  const int NProcs = ParallelDescriptor::NProcs();

  if (NProcs == 1)
  {
    BL_ASSERT(data.empty());
    return;
  }

  //
  // We may have data that needs to be sent to another CPU.
  //
  const int MyProc = ParallelDescriptor::MyProc();

  Array<int> Snds(NProcs,0), Rcvs(NProcs,0);

  int NumSnds = 0, NumRcvs = 0;

  for (const auto& kv : data)
  {
    NumSnds       += kv.second.size();
    Snds[kv.first] = kv.second.size();
  }

  ParallelDescriptor::ReduceIntMax(NumSnds);

  if (NumSnds == 0) {
    //
    // There's no parallel work to do.
    //
    return;
  }

  BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
		  ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

  BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
			       1,
			       ParallelDescriptor::Mpi_typemap<int>::type(),
			       Rcvs.dataPtr(),
			       1,
			       ParallelDescriptor::Mpi_typemap<int>::type(),
			       ParallelDescriptor::Communicator()) );
  BL_ASSERT(Rcvs[MyProc] == 0);

  BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
		  ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

  typedef std::map<int,int> IntIntMap;

  IntIntMap SndCnts, RcvCnts, rOffset;

  for (int i = 0; i < NProcs; i++) {
    if (Snds[i] > 0) {
      SndCnts[i] = Snds[i];
    }
  }

  for (int i = 0; i < NProcs; i++)
  {
    if (Rcvs[i] > 0)
    {
      RcvCnts[i] = Rcvs[i];
      rOffset[i] = NumRcvs;
      NumRcvs   += Rcvs[i];
    }
  }
  //
  // Don't need these anymore.
  //
  Array<int>().swap(Snds);
  Array<int>().swap(Rcvs);
  //
  // The data we want to receive.
  //
  // We only use: m_lev, m_grid, m_cell & m_data[0..ncomp-1] from the particles.
  //
  const int iChunkSize = 2 + BL_SPACEDIM;
  const int rChunkSize = ncomp;

  Array<int>                    irecvdata (NumRcvs*iChunkSize);
  Array<ParticleBase::RealType> rrecvdata (NumRcvs*rChunkSize);

  Array<int>         index(2*RcvCnts.size());
  Array<MPI_Status>  stats(2*RcvCnts.size());
  Array<MPI_Request> rreqs(2*RcvCnts.size());

  const int SeqNumI = ParallelDescriptor::SeqNum();
  const int SeqNumR = ParallelDescriptor::SeqNum();
  //
  // Post the receives.
  //
  int idx = 0;
  for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
  {
    const int Who  = it->first;
    const int iCnt = it->second   * iChunkSize;
    const int rCnt = it->second   * rChunkSize;
    const int iIdx = rOffset[Who] * iChunkSize;
    const int rIdx = rOffset[Who] * rChunkSize;

    BL_ASSERT(Who >= 0 && Who < NProcs);
    BL_ASSERT(iCnt > 0);
    BL_ASSERT(rCnt > 0);
    BL_ASSERT(iCnt < std::numeric_limits<int>::max());
    BL_ASSERT(rCnt < std::numeric_limits<int>::max());

    rreqs[2*idx+0] = ParallelDescriptor::Arecv(&irecvdata[iIdx],iCnt,Who,SeqNumI).req();
    rreqs[2*idx+1] = ParallelDescriptor::Arecv(&rrecvdata[rIdx],rCnt,Who,SeqNumR).req();
  }
  //
  // Send the data.
  //
  Array<int>                    isenddata;
  Array<ParticleBase::RealType> rsenddata;

  for (const auto& kv : SndCnts)
  {
    const int Who  = kv.first;
    const int iCnt = kv.second * iChunkSize;
    const int rCnt = kv.second * rChunkSize;

    BL_ASSERT(iCnt > 0);
    BL_ASSERT(rCnt > 0);
    BL_ASSERT(Who >= 0 && Who < NProcs);
    BL_ASSERT(iCnt < std::numeric_limits<int>::max());
    BL_ASSERT(rCnt < std::numeric_limits<int>::max());

    isenddata.resize(iCnt);
    rsenddata.resize(rCnt);

    PBox& pbox = data[Who];

    int ioff = 0, roff = 0;
    for (const auto& p : pbox)
    {
      isenddata[ioff+0] = p.m_lev  - lev_min;
      isenddata[ioff+1] = p.m_grid;

      D_TERM(isenddata[ioff+2] = p.m_cell[0];,
	     isenddata[ioff+3] = p.m_cell[1];,
	     isenddata[ioff+4] = p.m_cell[2];);

      ioff += iChunkSize;

      for (int n = 0; n < ncomp; n++) {
	rsenddata[roff+n] = p.m_data[n];
      }

      roff += ncomp;
    }

    PBox().swap(pbox);

    ParallelDescriptor::Send(isenddata.dataPtr(),iCnt,Who,SeqNumI);
    ParallelDescriptor::Send(rsenddata.dataPtr(),rCnt,Who,SeqNumR);
  }
  //
  // Receive the data.
  //
  for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
  {
    ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
  }
  //
  // Now update "mf".
  //
  if (NumRcvs > 0)
  {
    const int*                    idata = irecvdata.dataPtr();
    const ParticleBase::RealType* rdata = rrecvdata.dataPtr();

    for (int i = 0; i < NumRcvs; i++)
    {
      const int     lev  = idata[0];
      const int     grd  = idata[1];
      const IntVect cell = IntVect(D_DECL(idata[2],idata[3],idata[4]));

      BL_ASSERT((*mf)[lev].DistributionMap()[grd] == MyProc);
      BL_ASSERT((*mf)[lev][grd].box().contains(cell));

      for (int n = 0; n < ncomp; n++) {
	(*mf)[lev][grd](cell,n) += rdata[n];
      }

      idata += iChunkSize;
      rdata += rChunkSize;
    }
  }
}
  

