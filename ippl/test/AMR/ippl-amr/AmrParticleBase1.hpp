#include "Ippl.h"

#include <AMReX_Particles.H>
#include <AMReX_ParmParse.H>
#include <limits>

// Function from BoxLib adjusted to work with Ippl AmrParticleBase1 class
template<class PLayout>
void AmrParticleBase1<PLayout>::Interp(const SingleParticlePos_t &R,
				      const amrex::Geometry &geom,
				      const amrex::FArrayBox& fab,
				      const int* idx,
				      amrex::Real* val,
				      int cnt)
{
  const int       M   = D_TERM(2,+2,+4);
  const amrex::Real*     plo = geom.ProbLo();
  const amrex::Real*     dx  = geom.CellSize();
  
  amrex::Real    fracs[M];
  amrex::IntVect cells[M];
  //
  // Get "fracs" and "cells".
  //
  CIC_Cells_Fracs_Basic(R, plo, dx, fracs, cells);
  
  for (int i = 0; i < cnt; i++)
  {
    BL_ASSERT(idx[i] >= 0 && idx[i] < fab.nComp());
    // dummy template values
    val[i] = amrex::Particle<0, 0>::InterpDoit(fab,fracs,cells,idx[i]);
  }
}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase1 class
template<class PLayout>
void AmrParticleBase1<PLayout>::CIC_Cells_Fracs_Basic(const SingleParticlePos_t &R, 
						     const amrex::Real* plo, 
						     const amrex::Real* dx, 
						     amrex::Real* fracs,  
						     amrex::IntVect* cells)
{

  const amrex::Real len[AMREX_SPACEDIM] = { D_DECL((R[0]-plo[0])/dx[0] + amrex::Real(0.5),
					 (R[1]-plo[1])/dx[1] + amrex::Real(0.5),
					 (R[2]-plo[2])/dx[2] + amrex::Real(0.5)) };

  const amrex::IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));
  
  const amrex::Real frac[AMREX_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

  // dummy template values
  amrex::Particle<0, 0>::CIC_Fracs(frac, fracs);
  amrex::Particle<0, 0>::CIC_Cells(cell, cells);

}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase1 class
template<class PLayout>
int AmrParticleBase1<PLayout>::CIC_Cells_Fracs (const SingleParticlePos_t &R,
					       const amrex::Real*         plo,
					       const amrex::Real*         dx_geom,
					       const amrex::Real*         dx_part,
					       amrex::Vector<amrex::Real>&        fracs,
					       amrex::Vector<amrex::IntVect>&     cells)
{
    BL_PROFILE("AmrParticleBase1::CIC_Cells_Fracs()");
    if (dx_geom == dx_part)
    {
        const int M = D_TERM(2,+2,+4);
        fracs.resize(M);
        cells.resize(M);
        AmrParticleBase1::CIC_Cells_Fracs_Basic(R,plo,dx_geom,fracs.dataPtr(),cells.dataPtr());
        return M;
    }
    //
    // The first element in fracs and cells is the lowest corner, the last is the highest.
    //
    const amrex::Real hilen[AMREX_SPACEDIM] = { D_DECL((R[0]-plo[0]+dx_part[0]/2)/dx_geom[0],
                                             (R[1]-plo[1]+dx_part[1]/2)/dx_geom[1],
                                             (R[2]-plo[2]+dx_part[2]/2)/dx_geom[2]) };

    const amrex::Real lolen[AMREX_SPACEDIM] = { D_DECL((R[0]-plo[0]-dx_part[0]/2)/dx_geom[0],
                                             (R[1]-plo[1]-dx_part[1]/2)/dx_geom[1],
                                             (R[2]-plo[2]-dx_part[2]/2)/dx_geom[2]) };

    const amrex::IntVect hicell(D_DECL(floor(hilen[0]), floor(hilen[1]), floor(hilen[2])));
    
    const amrex::IntVect locell(D_DECL(floor(lolen[0]), floor(lolen[1]), floor(lolen[2])));
    
    const amrex::Real cell_density = D_TERM(dx_geom[0]/dx_part[0],*dx_geom[1]/dx_part[1],*dx_geom[2]/dx_part[2]);
    
    const int M = D_TERM((hicell[0]-locell[0]+1),*(hicell[1]-locell[1]+1),*(hicell[2]-locell[2]+1));

    fracs.resize(M);
    cells.resize(M);
    //
    // This portion might be slightly inefficient. Feel free to redo it if need be.
    //
    int i = 0;
    for (int zi = locell[2]; zi <= hicell[2]; zi++)
    {
        const amrex::Real zf = std::min(hilen[2]-zi,amrex::Real(1))-std::max(lolen[2]-zi,amrex::Real(0));
        for (int yi = locell[1]; yi <= hicell[1]; yi++)
        {
            const amrex::Real yf = std::min(hilen[1]-yi,amrex::Real(1))-std::max(lolen[1]-yi,amrex::Real(0));
            for (int xi = locell[0]; xi <= hicell[0]; xi++)
            {
                cells[i][0] = xi;
                cells[i][1] = yi;
                cells[i][2] = zi;
                fracs[i] = zf * yf * (std::min(hilen[0]-xi,amrex::Real(1))
				      -std::max(lolen[0]-xi,amrex::Real(0))) * cell_density;
                i++;
            }
        }
    }

    return M;
}

// Function from BoxLib adjusted to work with Ippl AmrParticleBase1 class
template<class PLayout>
bool AmrParticleBase1<PLayout>::FineToCrse (const int ip,
					   int                               flev,
					   const amrex::Vector<amrex::IntVect>&              fcells,
					   const amrex::BoxArray&                    fvalid,
					   const amrex::BoxArray&                    compfvalid_grown,
					   amrex::Vector<amrex::IntVect>&                    ccells,
					   amrex::Vector<amrex::Real>&                       cfracs,
					   amrex::Vector<int>&                        which,
					   amrex::Vector<int>&                        cgrid,
					   amrex::Vector<amrex::IntVect>&                    pshifts,
					   std::vector< std::pair<int,amrex::Box> >& isects)
{
    PLayout *Layout = &this->getLayout();
    const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
    
    BL_PROFILE("AmrParticleBase1::FineToCrse()");
    BL_ASSERT(m_gdb != 0);
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
    
//     Layout->Where(this->R, ip); //FIXME Is this really needed?
    
    const amrex::Box& ibx = amrex::grow(m_gdb->ParticleBoxArray(flev)[this->m_grid[ip]],-1);
    amrex::IntVect m_cell = Layout->Index(this->R[ip], flev);
    
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
    const amrex::Geometry& cgm = m_gdb->Geom(flev-1);
    const amrex::IntVect&  rr  = m_gdb->refRatio(flev-1);
    const amrex::BoxArray& cba = m_gdb->ParticleBoxArray(flev-1);
    
    CIC_Cells_Fracs(this->R[ip], cgm.ProbLo(), cgm.CellSize(), cgm.CellSize(), cfracs, ccells);
    
    bool result = false;
    
    for (int i = 0; i < M; i++)
    {
        amrex::IntVect ccell_refined = ccells[i]*rr;
        //
        // We've got to protect against the case when we're at the low
        // end of the domain because coarsening & refining don't work right
        // when indices go negative.
        //
        for (int dm = 0; dm < AMREX_SPACEDIM; dm++)
            ccell_refined[dm] = std::max(ccell_refined[dm], -1);
    
        if (!fvalid.contains(ccell_refined))
        {
            result   = true;
            which[i] = 1;
            
            amrex::Box cbx(ccells[i],ccells[i]);
    
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

// Function from BoxLib adjusted to work with Ippl AmrParticleBase1 class
template<class PLayout>
void AmrParticleBase1<PLayout>::FineCellsToUpdateFromCrse (
  const int ip,
  int lev,
  const amrex::IntVect& ccell,
  const amrex::IntVect& cshift,
  amrex::Vector<int>& fgrid,
  amrex::Vector<amrex::Real>& ffrac,
  amrex::Vector<amrex::IntVect>& fcells,
  std::vector< std::pair<int,amrex::Box> >& isects)
{
  BL_PROFILE("ParticleContainer<NR, NI, NA>::FineCellsToUpdateFromCrse()");
    
    PLayout *Layout = &this->getLayout();
    const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
    
    BL_ASSERT(lev >= 0);
    BL_ASSERT(lev < m_gdb->finestLevel());

    const amrex::Box&      fbx = amrex::refine(amrex::Box(ccell,ccell),m_gdb->refRatio(lev));
    const amrex::BoxArray& fba = m_gdb->ParticleBoxArray(lev+1);
    const amrex::Real*     plo = m_gdb->Geom(lev).ProbLo();
    const amrex::Real*     dx  = m_gdb->Geom(lev).CellSize();
    const amrex::Real*     fdx = m_gdb->Geom(lev+1).CellSize();

    if (cshift == amrex::IntVect::TheZeroVector())
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
    for (amrex::IntVect iv = fbx.smallEnd(); iv <= fbx.bigEnd(); fbx.next(iv))
    {
        bool touches = true;

        for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            const amrex::Real celllo = iv[k]  * fdx[k] + plo[k];
            const amrex::Real cellhi = celllo + fdx[k];

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

    amrex::Real sum_fine = 0;
    //
    // We need to figure out the fine fractions and the fine grid needed updating.
    //
    for (unsigned int j = 0; j < fcells.size(); j++)
    {
        amrex::IntVect& iv = fcells[j];

        amrex::Real the_frac = 1;

        for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            const amrex::Real celllo = (iv[k] * fdx[k] + plo[k]);

            if (this->R[ip][k] <= celllo)
            {
                const amrex::Real isecthi = this->R[ip][k] + dx[k]/2;

                the_frac *= std::min((isecthi - celllo),fdx[k]);
            }
            else
            {
                const amrex::Real cellhi  = (iv[k]+1) * fdx[k] + plo[k];
                const amrex::Real isectlo = this->R[ip][k] - dx[k]/2;

                the_frac *= std::min((cellhi - isectlo),fdx[k]);
            }
        }

        ffrac.push_back(the_frac);

        sum_fine += the_frac;

        if (cshift != amrex::IntVect::TheZeroVector())
        {
            //
            // Update to the correct fine cell needing updating.
            // Note that "cshift" is from the coarse perspective.
            //
            const amrex::IntVect& fshift = cshift * m_gdb->refRatio(lev);
            //
            // Update fcells[j] to indicate a shifted fine cell needing updating.
            //
            iv -= fshift;
        }

        fba.intersections(amrex::Box(iv,iv),isects,true,0);

        BL_ASSERT(!isects.empty());

        fgrid.push_back(isects[0].first);
    }

    BL_ASSERT(ffrac.size() == fcells.size());
    BL_ASSERT(fgrid.size() == fcells.size());
    //
    // Now adjust the fine fractions so they sum to one.
    //
    for (unsigned int j = 0; j < ffrac.size(); j++)
        ffrac[j] /= sum_fine;
}


template<class PLayout>
void AmrParticleBase<PLayout>::AssignDensityDoit(int rho_index,
						 amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf,
						 PMap&             data,
						 int               ncomp,
						 int               lev_min)
{
  if (rho_index != 0) amrex::Abort("AssignDensityDoit only works if rho_index = 0");

  BL_PROFILE("ParticleContainer<NR,NI,C>::AssignDensityDoit()");
  BL_ASSERT(1 >= ncomp);

  const int NProcs = amrex::ParallelDescriptor::NProcs();

  if (NProcs == 1)
  {
    BL_ASSERT(data.empty());
    return;
  }

  //
  // We may have data that needs to be sent to another CPU.
  //
  const int MyProc = amrex::ParallelDescriptor::MyProc();

  amrex::Vector<int> Snds(NProcs,0), Rcvs(NProcs,0);

  int NumSnds = 0, NumRcvs = 0;

  for (const auto& kv : data)
  {
    NumSnds       += kv.second.size();
    Snds[kv.first] = kv.second.size();
  }

  amrex::ParallelDescriptor::ReduceIntMax(NumSnds);

  if (NumSnds == 0) {
    //
    // There's no parallel work to do.
    //
    return;
  }

  BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
		  amrex::ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

  BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
			       1,
			       amrex::ParallelDescriptor::Mpi_typemap<int>::type(),
			       Rcvs.dataPtr(),
			       1,
			       amrex::ParallelDescriptor::Mpi_typemap<int>::type(),
			       amrex::ParallelDescriptor::Communicator()) );
  BL_ASSERT(Rcvs[MyProc] == 0);

  BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
		  amrex::ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

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
  amrex::Vector<int>().swap(Snds);
  amrex::Vector<int>().swap(Rcvs);
  //
  // The data we want to receive.
  //
  // We only use: m_lev, m_grid, m_cell & m_data[0..ncomp-1] from the particles.
  //
  const int iChunkSize = 2 + AMREX_SPACEDIM;
  const int rChunkSize = ncomp;

  amrex::Vector<int>                    irecvdata (NumRcvs*iChunkSize);
  amrex::Vector<amrex::ParticleCommData::RealType> rrecvdata (NumRcvs*rChunkSize);

  amrex::Vector<int>         index(2*RcvCnts.size());
  amrex::Vector<MPI_Status>  stats(2*RcvCnts.size());
  amrex::Vector<MPI_Request> rreqs(2*RcvCnts.size());

  const int SeqNumI = amrex::ParallelDescriptor::SeqNum();
  const int SeqNumR = amrex::ParallelDescriptor::SeqNum();
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

    rreqs[2*idx+0] = amrex::ParallelDescriptor::Arecv(&irecvdata[iIdx],iCnt,Who,SeqNumI).req();
    rreqs[2*idx+1] = amrex::ParallelDescriptor::Arecv(&rrecvdata[rIdx],rCnt,Who,SeqNumR).req();
  }
  //
  // Send the data.
  //
  amrex::Vector<int>                    isenddata;
  amrex::Vector<amrex::ParticleCommData::RealType> rsenddata;

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

    amrex::ParallelDescriptor::Send(isenddata.dataPtr(),iCnt,Who,SeqNumI);
    amrex::ParallelDescriptor::Send(rsenddata.dataPtr(),rCnt,Who,SeqNumR);
  }
  //
  // Receive the data.
  //
  for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
  {
    amrex::ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
  }
  //
  // Now update "mf".
  //
  if (NumRcvs > 0)
  {
    const int*                    idata = irecvdata.dataPtr();
    const amrex::ParticleCommData::RealType* rdata = rrecvdata.dataPtr();

    for (int i = 0; i < NumRcvs; i++)
    {
      const int     lev  = idata[0];
      const int     grd  = idata[1];
      const amrex::IntVect cell = amrex::IntVect(D_DECL(idata[2],idata[3],idata[4]));

      BL_ASSERT((*mf[lev]).DistributionMap()[grd] == MyProc);
      BL_ASSERT((*mf[lev])[grd].box().contains(cell));

      for (int n = 0; n < ncomp; n++) {
	(*mf[lev])[grd](cell,n) += rdata[n];
      }

      idata += iChunkSize;
      rdata += rChunkSize;
    }
  }
}


template<class PLayout>
template <class AType>
void AmrParticleBase<PLayout>::AssignDensityFort (ParticleAttrib<AType> &pa,
                                                  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_to_be_filled, 
                                                  int lev_min, int ncomp, int finest_level) const
{
//     BL_PROFILE("AssignDensityFort()");
    IpplTimings::startTimer(AssignDensityTimer_m);
    const PLayout *Layout = &this->getLayout();
    const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
    
    // not done in amrex
    int rho_index = 0;
    
    amrex::PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    amrex::Vector<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
    amrex::PCInterp mapper;
    
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const amrex::BoxArray& ba = mf_to_be_filled[lev]->boxArray();
        const amrex::DistributionMapping& dm = mf_to_be_filled[lev]->DistributionMap();
        tmp[lev].reset(new amrex::MultiFab(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }
    
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        AssignCellDensitySingleLevelFort(pa, *mf_to_be_filled[lev], lev, 1, 0);

        if (lev < finest_level) {
            amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *mf_to_be_filled[lev],
                                         rho_index, rho_index, ncomp, 
                                         m_gdb->Geom(lev), m_gdb->Geom(lev+1),
                                         cphysbc, fphysbc,
                                         m_gdb->refRatio(lev), &mapper, bcs);
        }

        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in 
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            sum_fine_to_coarse(*mf_to_be_filled[lev],
                               *mf_to_be_filled[lev-1],
                               rho_index, 1, m_gdb->refRatio(lev-1),
                               m_gdb->Geom(lev-1), m_gdb->Geom(lev));
        }
        
        mf_to_be_filled[lev]->plus(*tmp[lev], rho_index, ncomp, 0);
    }
    
    for (int lev = finest_level - 1; lev >= lev_min; --lev) {
        amrex::average_down(*mf_to_be_filled[lev+1], 
                            *mf_to_be_filled[lev], rho_index, ncomp, m_gdb->refRatio(lev));
    }
    
    IpplTimings::stopTimer(AssignDensityTimer_m);
}

// This is the single-level version for cell-centered density
template<class PLayout>
template <class AType>
void AmrParticleBase1<PLayout>::AssignCellDensitySingleLevelFort (ParticleAttrib<AType> &pa,
                                                                 amrex::MultiFab& mf_to_be_filled,
                                                                 int       lev,
                                                                 int       ncomp,
                                                                 int       particle_lvl_offset) const
{
//     BL_PROFILE("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignCellDensitySingleLevelFort()");
    
    const PLayout *Layout = &this->getLayout();
    const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
    
    int rho_index = 0;
    
//     if (rho_index != 0) amrex::Abort("AssignCellDensitySingleLevel only works if rho_index = 0");

    amrex::MultiFab* mf_pointer;

    if (m_gdb->OnSameGrids(lev, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the 
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need 
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new amrex::MultiFab(m_gdb->ParticleBoxArray(lev),
                                m_gdb->ParticleDistributionMap(lev),
                                ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread 
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

    const amrex::Real      strttime    = amrex::ParallelDescriptor::second();
    const amrex::Geometry& gm          = m_gdb->Geom(lev);
    const amrex::Real*     plo         = gm.ProbLo();
    const amrex::Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
    const amrex::Real*     dx          = gm.CellSize();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
      amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }
    
    for (amrex::MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }
    
    //loop trough particles and distribute values on the grid
    size_t LocalNum = this->getLocalNum();
    
    amrex::Real inv_dx[3] = { 1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = {0, 0, 0};
    
    size_t lBegin = LocalNumPerLevel_m.begin(lev);
    size_t lEnd   = LocalNumPerLevel_m.end(lev);
    
    for (size_t ip = lBegin; ip < lEnd; ++ip) {
        
        const int grid = m_grid[ip];
        amrex::FArrayBox& fab = (*mf_pointer)[grid];
        const amrex::Box& box = fab.box();
        
        // not callable:
        // begin amrex_deposit_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);
        for (int i = 0; i < 3; ++i) {
            lxyz[i] = ( this->R[ip](i) - plo[i] ) * inv_dx[i] + 0.5;
            ijk[i] = lxyz[i];
            wxyz_hi[i] = lxyz[i] - ijk[i];
            wxyz_lo[i] = 1.0 - wxyz_hi[i];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        amrex::IntVect i1(i-1, j-1, k-1);
        amrex::IntVect i2(i-1, j-1, k);
        amrex::IntVect i3(i-1, j,   k-1);
        amrex::IntVect i4(i-1, j,   k);
        amrex::IntVect i5(i,   j-1, k-1);
        amrex::IntVect i6(i,   j-1, k);
        amrex::IntVect i7(i,   j,   k-1);
        amrex::IntVect i8(i,   j,   k);
        
        fab(i1, 0) += wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*pa[ip];
        fab(i2, 0) += wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*pa[ip];
        fab(i3, 0) += wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*pa[ip];
        fab(i4, 0) += wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*pa[ip];
        fab(i5, 0) += wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*pa[ip];
        fab(i6, 0) += wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*pa[ip];
        fab(i7, 0) += wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*pa[ip];
        fab(i8, 0) += wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*pa[ip];
        // end of amrex_deposit_cic
    }
    
    mf_pointer->SumBoundary(gm.periodicity());

    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const amrex::Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1.0/vol, 0, 1, mf_pointer->nGrow());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }
}

template<class PLayout>
template <class AType>
void AmrParticleBase<PLayout>::InterpolateFort (ParticleAttrib<AType> &pa,
                                                amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mesh_data, 
                                                int lev_min, int lev_max)
{
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        InterpolateSingleLevelFort(pa, *mesh_data[lev], lev); 
    }
}

template<class PLayout>
template <class AType>
void AmrParticleBase1<PLayout>::InterpolateSingleLevelFort (ParticleAttrib<AType> &pa,
                                                           amrex::MultiFab& mesh_data, int lev)
{
    if (mesh_data.nGrow() < 1)
        amrex::Error("Must have at least one ghost cell when in InterpolateSingleLevelFort");
    
    PLayout *Layout = &this->getLayout();
    const amrex::ParGDBBase* m_gdb = Layout->GetParGDB();
    
    const amrex::Geometry& gm          = m_gdb->Geom(lev);
    const amrex::Real*     plo         = gm.ProbLo();
    const amrex::Real*     dx          = gm.CellSize();
    
    //loop trough particles and distribute values on the grid
    size_t LocalNum = this->getLocalNum();
    
    amrex::Real inv_dx[3] = { 1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = {0, 0, 0};
    
    size_t lBegin = this->LocalNumPerLevel_m.begin(lev);
    size_t lEnd   = this->LocalNumPerLevel_m.end(lev);
    
    for (size_t ip = lBegin; ip < lEnd; ++ip) {
        
        const int grid = m_grid[ip];
        amrex::FArrayBox& fab = mesh_data[grid];
        const amrex::Box& box = fab.box();
        int nComp = fab.nComp();
        
        // not callable
        // begin amrex_interpolate_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), nComp, plo, dx);
        for (int i = 0; i < 3; ++i) {
            lxyz[i] = ( this->R[ip](i) - plo[i] ) * inv_dx[i] + 0.5;
            ijk[i] = lxyz[i];
            wxyz_hi[i] = lxyz[i] - ijk[i];
            wxyz_lo[i] = 1.0 - wxyz_hi[i];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        amrex::IntVect i1(i-1, j-1, k-1);
        amrex::IntVect i2(i-1, j-1, k);
        amrex::IntVect i3(i-1, j,   k-1);
        amrex::IntVect i4(i-1, j,   k);
        amrex::IntVect i5(i,   j-1, k-1);
        amrex::IntVect i6(i,   j-1, k);
        amrex::IntVect i7(i,   j,   k-1);
        amrex::IntVect i8(i,   j,   k);
        
        for (int nc = 0; nc < nComp; ++nc) {
            pa[ip](nc) = wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*fab(i1, nc) +
                         wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*fab(i2, nc) +
                         wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*fab(i3, nc) +
                         wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*fab(i4, nc) +
                         wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*fab(i5, nc) +
                         wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*fab(i6, nc) +
                         wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*fab(i7, nc) +
                         wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*fab(i8, nc);
        }
        // end amrex_interpolate_cic
    }
}
